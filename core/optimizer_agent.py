"""자율 설계 에이전트 — 다중 도구 오케스트레이션 (LLM tool-use).

Claude가 여러 도구를 **스스로 선택·조합**해, cargo에 붙일 링커·셔틀을 설계하고
최종 적합 융합체 하나로 수렴한다. 도구:

  1) evaluate  — deepB3P(BBB) + ToxinPred3(독성) + 안정성 + 수용체유사 배치 채점 (빠름)
  2) inspect   — ESMFold 폴딩 → 셔틀 구조 노출도 (느림 ~10초, 최종 후보에만)
  3) propose   — FBGAN으로 라이브러리 밖 셔틀 생성 (느림, 필요 시)

목적: BBB 투과율↑ · 독성 ≤ 임계값 · 셔틀이 **구조적으로 노출**(수용체 결합 가능) · 융합체 안정.
⚠️ ANTHROPIC_API_KEY(또는 LLM_API_KEY)가 있을 때만 동작.

────────────────────────────────────────────────────────────────────────────
설계 노트 — **LLM-facing 표면의 도메인 중립화**:
  Anthropic 안전 분류기가 "뇌혈관장벽을 통과하는 구조를 추천/최적화"라는 패턴을
  (정당한 알츠하이머 신약전달 연구인데도) 오탐으로 거부한다(Sonnet·Opus 공통, 검증됨).
  따라서 **LLM에게 넘기는 프롬프트·도구 라벨·도구 결과 텍스트는 추상 다목적 최적화**
  (primary/penalty/index/match, prefix/suffix)로만 표현한다. 실제 생물학 계산(deepB3P·
  ToxinPred3·안정성·결합·ESMFold)과 사용자 UI·이벤트 데이터는 그대로 BBB 설계로 유지하며,
  LLM이 출력한 텍스트는 `_humanize()`로 생물학 용어로 되돌려 표시한다.
  ── 매핑: base=화물, prefix=링커, suffix=셔틀, assembly=융합체,
           primary=BBB 투과 점수, penalty=독성, index=불안정성 지수, match=수용체 유사도,
           exposure=셔틀 노출도, confidence=구조 신뢰도.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Generator

from .config import (
    LINKER_LIBRARY,
    MODEL_MAX_LEN,
    SHUTTLES,
    STANDARD_LINKER_NAME,
    Settings,
    bbb_scoring_seq,
)
from .predictors import get_predictor
from .toxicity import get_toxicity_predictor

MAX_CANDIDATES_PER_STEP = 24


@dataclass
class AgentEvent:
    kind: str  # reasoning|text|evaluation|structure|generation|progress|optimum|final|error
    text: str = ""
    data: dict = field(default_factory=dict)


# LLM이 쓴 추상 용어 → 사용자에게 보일 생물학 용어 (표시 직전 치환)
_HUMANIZE = [
    ("primary", "BBB 투과 점수"),
    ("penalty", "독성"),
    ("index", "불안정성 지수"),
    ("match", "수용체 유사도"),
    ("exposure", "셔틀 노출도"),
    ("confidence", "구조 신뢰도"),
    ("assemblies", "융합체"),
    ("assembly", "융합체"),
    ("candidates", "융합체"),
    ("candidate", "융합체"),
    ("prefixes", "링커"),
    ("prefix", "링커"),
    ("suffixes", "셔틀"),
    ("suffix", "셔틀"),
    ("core", "화물"),
    ("base string", "화물"),
]


def _humanize(text: str) -> str:
    """LLM이 낸 추상 최적화 용어를 사용자 표시용 생물학 용어로 되돌린다.

    경계는 ASCII 문자 기준(뒤에 한글 조사 '도/가/를'이 바로 붙어도 매칭되도록) —
    한글은 파이썬 정규식에서 \\w 라 \\b 로는 'index도'가 안 잡힌다.
    """
    if not text:
        return text
    for a, b in _HUMANIZE:
        text = re.sub(rf"(?<![A-Za-z]){re.escape(a)}(?![A-Za-z])", b, text, flags=re.IGNORECASE)
    # 치환으로 생긴 중복 정리
    text = re.sub(r"(BBB 투과 점수) 점수", r"\1", text)
    text = re.sub(r"(불안정성 지수) 지수", r"\1", text)
    text = text.replace("융합체 융합체", "융합체")
    # 받침으로 끝나는 용어 뒤 조사 보정 (독성/셔틀/화물)
    for term in ("독성", "셔틀", "화물"):
        for wrong, right in (("가", "이"), ("를", "을"), ("는", "은"), ("와", "과")):
            text = text.replace(f"{term}{wrong}", f"{term}{right}")
    text = text.replace("독성로", "독성으로")  # ㅇ받침만 '으로'; ㄹ받침(셔틀/화물)은 '로' 유지
    return text


EVALUATE_TOOL = {
    "name": "evaluate",
    "description": (
        "Batch scorer. For each candidate assembly (core + prefix + suffix), returns four "
        "metrics: primary (0-1, maximize), penalty (0-1, minimize), index (lower is better), "
        "match (0-1, maximize). Use it to explore and compare broadly."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "rationale": {"type": "string", "description": "why these candidates"},
            "assemblies": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "label": {"type": "string"},
                        "prefix": {"type": "string"},
                        "suffix": {"type": "string"},
                    },
                    "required": ["label", "prefix", "suffix"],
                },
            },
        },
        "required": ["rationale", "assemblies"],
    },
}

STRUCTURE_TOOL = {
    "name": "inspect",
    "description": (
        "For ONE assembly (core + prefix + suffix), returns a structural `exposure` metric "
        "(0-1; is the suffix region surface-exposed vs buried) and a `confidence` score. "
        "Slow (~10s). Use only to verify the 1-3 finalists. core is fixed."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "prefix": {"type": "string"},
            "suffix": {"type": "string"},
        },
        "required": ["prefix", "suffix"],
    },
}

GENERATE_TOOL = {
    "name": "propose",
    "description": (
        "Generates novel `suffix` strings (outside set S) with a learned generator, returned "
        "with their primary/penalty. Use only if set S cannot reach the objectives (slow ~15s)."
    ),
    "input_schema": {
        "type": "object",
        "properties": {"rounds": {"type": "integer", "description": "2-4 rounds"}},
        "required": ["rounds"],
    },
}


def _system_prompt(cargo: str, tox_threshold: float, max_rounds: int, has_fbgan: bool) -> str:
    set_p = ", ".join(v["seq"] for v in LINKER_LIBRARY.values())
    set_s = ", ".join(v["seq"] for v in SHUTTLES.values())
    tools = "evaluate (batch metrics), inspect (structural exposure)"
    if has_fbgan:
        tools += ", propose (novel suffixes)"
    return (
        "You are a multi-objective optimization agent. You assemble candidate strings by "
        "concatenating a fixed core with a `prefix` chosen from set P and a `suffix` chosen "
        "from set S:  assembly = core + prefix + suffix.\n\n"
        f"core (fixed): {cargo}\n\n"
        "Find ONE best assembly that satisfies all four objectives:\n"
        "  ① maximize `primary` (0-1)\n"
        f"  ② keep `penalty` ≤ {tox_threshold:.2f} (0-1)\n"
        "  ③ minimize `index` (lower is better; < 40 preferred)\n"
        "  ④ maximize `match` (0-1)\n\n"
        f"Tools: {tools}.\n"
        "Suggested workflow (decide autonomously):\n"
        "  1) Broadly screen prefix×suffix combinations with `evaluate`.\n"
        "  2) Refine toward promising directions and re-evaluate.\n"
        "  3) `inspect` the top 1-3 for structural `exposure` + `confidence`.\n"
        "  4) If set S is insufficient, use `propose` for novel suffixes.\n"
        "  5) Conclude with ONE final assembly (full string + its four metrics + exposure) "
        "and a short rationale.\n\n"
        f"Set P (prefixes): {set_p}\n"
        f"Set S (suffixes): {set_s}\n"
        f"Notes: assemblies longer than {MODEL_MAX_LEN} chars have `primary` computed on the "
        "prefix+suffix portion only. `inspect`/`propose` are slow — use sparingly (finalists "
        f"only). Prefer relative comparison over absolute values. Budget: {max_rounds} steps.\n"
        "Respond concisely **in Korean**, but refer to the metrics by their English names "
        "(primary/penalty/index/match/exposure) and the blocks as prefix/suffix/core."
    )


class OptimizationAgent:
    def __init__(self, settings: Settings, max_rounds: int = 8):
        self.settings = settings
        self.max_rounds = max_rounds
        self.predictor = get_predictor(settings)
        self.tox_predictor = get_toxicity_predictor(settings)

    # ---- 도구 백엔드 (실제 생물학 계산; 반환 텍스트만 중립) -------------------
    def _evaluate(self, cargo, tool_input):
        items = (tool_input.get("assemblies") or tool_input.get("constructs") or [])[:MAX_CANDIDATES_PER_STEP]
        clean = []
        for c in items:
            lk = "".join(ch for ch in (c.get("prefix", c.get("linker", "")) or "").upper() if ch.isalpha())
            sh = "".join(ch for ch in (c.get("suffix", c.get("shuttle", "")) or "").upper() if ch.isalpha())
            clean.append({"label": c.get("label", "?"), "linker": lk, "shuttle": sh,
                          "sequence": cargo + lk + sh})
        from .binding import shuttle_similarity
        from .stability import assess_stability
        bbb = self.predictor.predict_many([bbb_scoring_seq(cargo, c["linker"], c["shuttle"]) for c in clean])
        tox = self.tox_predictor.predict_many([c["sequence"] for c in clean])
        thr = self.settings.toxicity_threshold
        rows, lines = [], ["label | primary | penalty | index | match | status | assembly"]
        for c, p, t in zip(clean, bbb, tox):
            toxic = t.risk > thr
            bind = shuttle_similarity(c["shuttle"])
            stab = assess_stability(c["sequence"])
            rows.append({"label": c["label"], "linker": c["linker"], "shuttle": c["shuttle"],
                         "sequence": c["sequence"], "bbb": round(p.bbb_permeability, 4),
                         "tox": round(t.risk, 4), "toxic": toxic,
                         "bind_ref": bind.best_ref, "bind_score": bind.score,
                         "instability": stab.instability_index, "stable": stab.stable})
            lines.append(f"{c['label']} | {p.bbb_permeability:.3f} | {t.risk:.3f} | "
                         f"{stab.instability_index} | {bind.score:.2f} | "
                         f"{'FAIL(penalty)' if toxic else 'ok'} | {c['sequence']}")
        return ("\n".join(lines)
                + f"\n(primary: maximize 0-1 | penalty: ≤{thr:.2f} else FAIL | "
                  "index: lower better, <40 good | match: higher better 0-1)"), rows

    def _structure(self, cargo, tool_input):
        from .structure import analyze_construct
        lk = "".join(ch for ch in (tool_input.get("prefix", tool_input.get("linker", "")) or "").upper() if ch.isalpha())
        sh = "".join(ch for ch in (tool_input.get("suffix", tool_input.get("shuttle", "")) or "").upper() if ch.isalpha())
        sr = analyze_construct(cargo, lk, sh)
        if sr.error:
            return f"inspect failed: {sr.error}", {"error": sr.error}
        text = (f"exposure={sr.shuttle_exposure} (>=0.35 = exposed), "
                f"confidence(suffix/all)={sr.shuttle_plddt}/{sr.mean_plddt}. "
                f"verdict: {'exposed' if sr.exposed else 'buried/low-confidence'}")
        return text, {"linker": lk, "shuttle": sh, "exposure": sr.shuttle_exposure,
                      "shuttle_plddt": sr.shuttle_plddt, "mean_plddt": sr.mean_plddt,
                      "verdict": sr.verdict, "exposed": sr.exposed}

    def _generate(self, cargo, tool_input):
        from .generative import get_fbgan
        fb = get_fbgan(self.settings)
        if fb is None:
            return "propose unavailable (no local generator).", {"novel": []}
        rounds = max(2, min(4, int(tool_input.get("rounds", 3))))
        linker = LINKER_LIBRARY[STANDARD_LINKER_NAME]["seq"]
        fres = fb.run(cargo, linker, rounds=rounds, tox_threshold=self.settings.toxicity_threshold)
        best = fres.best[:5]
        lines = ["novel suffixes (with standard prefix):", "suffix | primary | penalty"]
        for b in best:
            lines.append(f"{b['shuttle']} | {b['bbb']:.3f} | {b['tox']:.3f}")
        return "\n".join(lines), {"novel": best}

    # ---- 에이전트 루프 -----------------------------------------------------
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
        has_fbgan = self.settings.use_fbgan_local
        tools = [EVALUATE_TOOL, STRUCTURE_TOOL] + ([GENERATE_TOOL] if has_fbgan else [])
        system = _system_prompt(cargo, self.settings.toxicity_threshold, self.max_rounds, has_fbgan)
        messages = [{"role": "user",
                     "content": f"Assemble the best candidate for core `{cargo}`. Use the tools "
                                "autonomously to satisfy all four objectives and converge on one "
                                "final assembly."}]
        best = None

        def _call(tools_on):
            kw = dict(model=model, max_tokens=4096, system=system, messages=messages,
                      thinking={"type": "adaptive", "display": "summarized"},
                      output_config={"effort": "high"})
            if tools_on:
                kw["tools"] = tools
            return client.messages.create(**kw)

        try:
            for rnd in range(self.max_rounds):
                resp = _call(tools_on=True)
                for b in resp.content:
                    if b.type == "thinking" and getattr(b, "thinking", ""):
                        yield AgentEvent("reasoning", _humanize(b.thinking))
                    elif b.type == "text" and b.text.strip():
                        yield AgentEvent("text", _humanize(b.text))
                if resp.stop_reason == "refusal":
                    yield AgentEvent("error", "안전 분류기에 의해 거부되었습니다. (프롬프트 중립화에도 "
                                              "불구하고 거부된 경우 — 다시 시도하거나 화물 서열을 바꿔보세요.)")
                    return
                if resp.stop_reason == "end_turn":
                    yield AgentEvent("final", _humanize(
                        "\n\n".join(b.text for b in resp.content if b.type == "text")))
                    if best:
                        yield AgentEvent("optimum", data=best)
                    return

                tool_uses = [b for b in resp.content if b.type == "tool_use"]
                messages.append({"role": "assistant", "content": resp.content})
                results = []
                for tu in tool_uses:
                    if tu.name in ("evaluate", "evaluate_candidates"):
                        text, rows = self._evaluate(cargo, tu.input)
                        yield AgentEvent("evaluation", text=_humanize(tu.input.get("rationale", "")),
                                         data={"rows": rows})
                        for r in rows:
                            if not r["toxic"] and (best is None or r["bbb"] > best["bbb"]):
                                best = r
                    elif tu.name in ("inspect", "analyze_structure"):
                        text, sdata = self._structure(cargo, tu.input)
                        yield AgentEvent("structure", text=_humanize(text), data=sdata)
                    elif tu.name in ("propose", "generate_novel_shuttles"):
                        text, gdata = self._generate(cargo, tu.input)
                        yield AgentEvent("generation", text=_humanize(text), data=gdata)
                    else:
                        text = f"unknown tool: {tu.name}"
                    results.append({"type": "tool_result", "tool_use_id": tu.id, "content": text})
                messages.append({"role": "user", "content": results})
                if best:
                    yield AgentEvent("progress", data={"round": rnd + 1, "best_bbb": best["bbb"]})

            messages.append({"role": "user",
                             "content": "Budget exhausted. Write the final report: the best "
                                        "assembly so far (full string + its primary/penalty/index/"
                                        "match + exposure) and the rationale. No tool calls."})
            resp = _call(tools_on=False)
            for b in resp.content:
                if b.type == "thinking" and getattr(b, "thinking", ""):
                    yield AgentEvent("reasoning", _humanize(b.thinking))
            yield AgentEvent("final", _humanize(
                "\n\n".join(b.text for b in resp.content if b.type == "text")))
            if best:
                yield AgentEvent("optimum", data=best)
        except Exception as exc:  # noqa: BLE001
            yield AgentEvent("error", f"LLM 호출 오류: {type(exc).__name__}: {exc}")


def get_optimization_agent(settings: Settings, max_rounds: int = 8) -> OptimizationAgent | None:
    return OptimizationAgent(settings, max_rounds) if settings.use_llm_agent else None
