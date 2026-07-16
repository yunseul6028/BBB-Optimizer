"""도구 활용 기반 분자 최적화 에이전트 (LLM tool-use loop).

Claude(Anthropic API)가 **옵티마이저**로서, cargo에 붙일 링커·셔틀을 설계하는
분자 최적화 루프를 자율 수행한다. 우리 로컬 엔진(deepB3P=BBB, ToxinPred3=독성)이
`evaluate_candidates` 툴의 뒷단에서 **목적함수 평가**를 담당한다.

    후보 제안 → 평가(도구 호출) → 결과 관찰 → 개선 후보 재제안 → 수렴
    목적함수: BBB 투과율(deepB3P) 최대화, 독성(ToxinPred3) ≤ 임계값

best-so-far를 추적해 매 라운드의 추론·평가·최적값 갱신을 UI로 스트리밍한다.
⚠️ ANTHROPIC_API_KEY(또는 LLM_API_KEY)가 있을 때만 동작.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Generator

from .config import MODEL_MAX_LEN, Settings, bbb_scoring_seq
from .predictors import get_predictor
from .toxicity import get_toxicity_predictor

MAX_CANDIDATES_PER_STEP = 24


@dataclass
class AgentEvent:
    kind: str  # reasoning | text | evaluation | progress | optimum | final | error
    text: str = ""
    data: dict = field(default_factory=dict)


EVALUATE_TOOL = {
    "name": "evaluate_candidates",
    "description": (
        "목적함수 평가 도구. 후보 융합체(construct)들을 한 번에 조립·측정한다. "
        "각 construct = 고정 화물(cargo) + linker + shuttle. deepB3P로 BBB 투과 확률을, "
        "ToxinPred3로 독성 확률을 실측해 돌려준다. 이전 결과에서 좋았던 방향으로 후보를 "
        "개선해 반복 제안하라. 표준 20종 아미노산 1글자 코드만 사용."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "rationale": {"type": "string",
                          "description": "이 후보들을 제안한 근거(직전 결과 대비 무엇을 개선했는지)."},
            "constructs": {
                "type": "array",
                "description": "평가할 후보 융합체 목록.",
                "items": {
                    "type": "object",
                    "properties": {
                        "label": {"type": "string"},
                        "linker": {"type": "string", "description": "링커 서열(없으면 빈 문자열)"},
                        "shuttle": {"type": "string", "description": "셔틀 서열(없으면 빈 문자열)"},
                    },
                    "required": ["label", "linker", "shuttle"],
                },
            },
        },
        "required": ["rationale", "constructs"],
    },
}


def _system_prompt(cargo: str, tox_threshold: float, max_rounds: int) -> str:
    return (
        "당신은 **도구를 활용하는 분자 최적화 에이전트**입니다. 뇌혈관장벽(BBB)을 통과하는 "
        "항체-셔틀 융합 단백질을 설계합니다.\n\n"
        f"고정 화물(cargo): {cargo}\n"
        "목적함수(최대화): construct = cargo + linker + shuttle 의 BBB 투과율(deepB3P).\n"
        f"제약조건: 독성(ToxinPred3) ≤ {tox_threshold:.2f}. 이를 초과하면 후보 탈락.\n\n"
        "최적화 루프를 자율 수행하라:\n"
        "1) 제안: 지금까지의 평가 결과에서 무엇이 목적함수를 높였는지 근거를 들어 "
        "**개선된 후보군**을 만든다(탐험+활용 균형; 변수 하나씩 바꾼 대조도 유용).\n"
        "2) 평가: evaluate_candidates 도구로 BBB·독성을 실측한다.\n"
        "3) 관찰·개선: 좋았던 방향으로 다음 후보를 개선해 재제안한다. best-so-far를 계속 갱신.\n"
        f"4) 종료: 최대 {max_rounds}라운드. 소진하거나 수렴하면, **최적 융합체**(전체 서열, "
        "예측 BBB/독성)와 최적화 과정 요약(어떤 특성이 목적함수를 끌어올렸는지)으로 보고서를 쓴다.\n\n"
        f"참고: construct가 {MODEL_MAX_LEN}aa를 넘으면 **연결부위(링커+셔틀=C말단)** 기준으로 BBB를 "
        "계산한다(화물 벌크 제외). linker·shuttle은 항상 C말단이라 보존된다. deepB3P는 짧은 "
        "펩타이드로 학습돼 절대값보다 **후보 간 상대 비교**가 신뢰된다. "
        "모든 서술은 한국어로 간결하게."
    )


class OptimizationAgent:
    def __init__(self, settings: Settings, max_rounds: int = 6):
        self.settings = settings
        self.max_rounds = max_rounds
        self.predictor = get_predictor(settings)
        self.tox_predictor = get_toxicity_predictor(settings)

    def _evaluate(self, cargo: str, tool_input: dict) -> tuple[str, list[dict]]:
        constructs = (tool_input.get("constructs") or [])[:MAX_CANDIDATES_PER_STEP]
        clean = []
        for c in constructs:
            linker = "".join(ch for ch in (c.get("linker", "") or "").upper() if ch.isalpha())
            shuttle = "".join(ch for ch in (c.get("shuttle", "") or "").upper() if ch.isalpha())
            clean.append({"label": c.get("label", "?"), "linker": linker, "shuttle": shuttle,
                          "sequence": cargo + linker + shuttle})

        seqs = [c["sequence"] for c in clean]
        # BBB는 연결부위 윈도우, 독성은 전체 서열
        bbb_seqs = [bbb_scoring_seq(cargo, c["linker"], c["shuttle"]) for c in clean]
        preds = self.predictor.predict_many(bbb_seqs)
        toxes = self.tox_predictor.predict_many(seqs)

        rows, lines = [], ["label | BBB% | 독성% | 판정 | 길이 | 서열"]
        thr = self.settings.toxicity_threshold
        for c, p, t in zip(clean, preds, toxes):
            toxic = t.risk > thr
            trunc = len(c["sequence"]) > MODEL_MAX_LEN
            rows.append({"label": c["label"], "linker": c["linker"], "shuttle": c["shuttle"],
                         "sequence": c["sequence"], "bbb": round(p.bbb_permeability, 4),
                         "tox": round(t.risk, 4), "toxic": toxic, "truncated": trunc})
            lines.append(f"{c['label']} | {p.bbb_permeability*100:.1f} | {t.risk*100:.0f} | "
                         f"{'독성탈락' if toxic else '통과'} | {len(c['sequence'])}aa"
                         f"{'⚠' if trunc else ''} | {c['sequence']}")
        text = ("\n".join(lines)
                + f"\n(독성 임계값 {thr:.2f} 초과=탈락. BBB=deepB3P, 독성=ToxinPred3 실측.)")
        return text, rows

    def run(self, cargo: str) -> Generator[AgentEvent, None, None]:
        cargo = cargo.strip().upper()
        try:
            import anthropic
        except ImportError:
            yield AgentEvent("error", "anthropic SDK 미설치: pip install anthropic")
            return
        if not self.settings.llm_api_key:
            yield AgentEvent("error", "ANTHROPIC_API_KEY(또는 LLM_API_KEY) 미설정.")
            return

        client = anthropic.Anthropic(api_key=self.settings.llm_api_key)
        model = self.settings.llm_model
        system = _system_prompt(cargo, self.settings.toxicity_threshold, self.max_rounds)
        messages: list[dict] = [{
            "role": "user",
            "content": f"화물 `{cargo}`의 BBB 투과율을 최대화(독성 제약 하)하는 링커·셔틀을 "
                       "찾아라. 첫 후보군을 제안하고 평가하라.",
        }]
        best: dict | None = None

        def _call(tools_on: bool):
            kwargs = dict(model=model, max_tokens=4096, system=system, messages=messages,
                          thinking={"type": "adaptive", "display": "summarized"},
                          output_config={"effort": "high"})
            if tools_on:
                kwargs["tools"] = [EVALUATE_TOOL]
            return client.messages.create(**kwargs)

        try:
            for rnd in range(self.max_rounds):
                resp = _call(tools_on=True)
                for b in resp.content:
                    if b.type == "thinking" and getattr(b, "thinking", ""):
                        yield AgentEvent("reasoning", b.thinking)
                    elif b.type == "text" and b.text.strip():
                        yield AgentEvent("text", b.text)

                if resp.stop_reason == "refusal":
                    yield AgentEvent("error", "안전 분류기에 의해 요청이 거부되었습니다.")
                    return
                if resp.stop_reason == "end_turn":
                    yield AgentEvent("final", "\n\n".join(b.text for b in resp.content if b.type == "text"))
                    if best:
                        yield AgentEvent("optimum", data=best)
                    return

                tool_uses = [b for b in resp.content if b.type == "tool_use"]
                messages.append({"role": "assistant", "content": resp.content})
                results = []
                for tu in tool_uses:
                    text, rows = self._evaluate(cargo, tu.input)
                    yield AgentEvent("evaluation", text=tu.input.get("rationale", ""),
                                     data={"rows": rows})
                    for r in rows:
                        if not r["toxic"] and (best is None or r["bbb"] > best["bbb"]):
                            best = r
                    results.append({"type": "tool_result", "tool_use_id": tu.id, "content": text})
                messages.append({"role": "user", "content": results})
                yield AgentEvent("progress",
                                 data={"round": rnd + 1, "best_bbb": best["bbb"] if best else 0.0})

            # 라운드 소진 → 최종 보고 강제
            messages.append({"role": "user",
                             "content": "최적화 예산을 모두 사용했습니다. 지금까지의 최적 융합체"
                                        "(전체 서열, 예측 BBB/독성)와 무엇이 목적함수를 끌어올렸는지 "
                                        "요약해 최종 보고서를 쓰세요. 더 이상 도구를 호출하지 마세요."})
            resp = _call(tools_on=False)
            for b in resp.content:
                if b.type == "thinking" and getattr(b, "thinking", ""):
                    yield AgentEvent("reasoning", b.thinking)
            yield AgentEvent("final", "\n\n".join(b.text for b in resp.content if b.type == "text"))
            if best:
                yield AgentEvent("optimum", data=best)
        except Exception as exc:  # noqa: BLE001
            yield AgentEvent("error", f"LLM 호출 오류: {type(exc).__name__}: {exc}")


def get_optimization_agent(settings: Settings, max_rounds: int = 6) -> OptimizationAgent | None:
    return OptimizationAgent(settings, max_rounds) if settings.use_llm_agent else None
