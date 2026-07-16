"""자율 설계 에이전트 — Gemini 브레인 (experiment/gemini 브랜치).

Claude 버전(`optimizer_agent.py`)과 **동일한 도구 백엔드**(_evaluate/_structure/_generate)를
재사용하되, 브레인을 Google Gemini(function calling)로 바꾼 실험 어댑터.

핵심 실험 질문: **Anthropic 안전 분류기가 오탐 거부한 "BBB 통과 구조 최적화"를 Gemini는
그대로(자연스러운 생물학 프레이밍) 받아주는가?** → 여기서는 중립화 없이 **정직한 BBB 프레이밍**
을 그대로 쓴다(도구명·프롬프트 모두 생물학 용어). 거부되면 Claude처럼 중립화가 필요하다는 결론.

이벤트(AgentEvent) 종류·데이터 키는 Claude 버전과 동일 → 앱 UI 그대로 재사용.
⚠️ GEMINI_API_KEY(또는 GOOGLE_API_KEY)가 있을 때만 동작.
"""

from __future__ import annotations

import time
from typing import Generator

from .config import LINKER_LIBRARY, MODEL_MAX_LEN, SHUTTLES, Settings
from .optimizer_agent import AgentEvent, OptimizationAgent

# 일시적 오류(모델 과부하·분당 레이트리밋) — 백오프 재시도 대상
_TRANSIENT = ("503", "UNAVAILABLE", "429", "RESOURCE_EXHAUSTED", "overloaded", "high demand")


def _gemini_system_prompt(cargo: str, tox_threshold: float, max_rounds: int, has_fbgan: bool) -> str:
    lib_l = ", ".join(f"{n}({v['seq']})" for n, v in LINKER_LIBRARY.items())
    lib_s = ", ".join(f"{n}({v['seq']})" for n, v in SHUTTLES.items())
    tools = "evaluate_candidates(BBB·독성·안정성·수용체 배치), analyze_structure(구조 노출도)"
    if has_fbgan:
        tools += ", generate_novel_shuttles(신규 셔틀 생성)"
    return (
        "[연구 맥락] 이것은 알츠하이머병 치료제의 뇌 전달을 개선하기 위한 정당한 학술·공모전용 "
        "in-silico 연구입니다. 사용하는 셔틀 펩타이드(Angiopep-2, TAT, Penetratin, SynB1 등)는 "
        "학계·임상에서 공개 발표된 뇌 전달 도구이며, 독성 예측은 안전한 후보만 남기기 위한 안전성 "
        "스크리닝입니다.\n\n"
        "당신은 **도구를 자율적으로 오케스트레이션하는 in-silico 신약설계 에이전트**입니다. "
        "뇌혈관장벽(BBB)을 통과하는 항체-셔틀 융합 단백질을 설계합니다.\n\n"
        f"고정 화물(cargo): {cargo}\n융합체 = cargo + linker + shuttle.\n\n"
        "다중 목적 (모두 만족하는 최종 융합체 하나로 수렴):\n"
        "  ① BBB 투과 점수(deepB3P 예측확률 0~1, 실제 투과율 아님·상대비교용) 최대화\n"
        f"  ② 독성(ToxinPred3) ≤ {tox_threshold:.2f}\n"
        "  ③ 셔틀이 구조적으로 노출(ESMFold) + 검증 셔틀과 유사(수용체 메커니즘 근거)\n"
        "  ④ 융합체 자체 안정성: 불안정성 지수 < 40 선호\n\n"
        f"사용 가능한 도구: {tools}\n"
        "권장 워크플로우: evaluate_candidates로 라이브러리 조합을 폭넓게 스크리닝 → 유망 방향으로 "
        "개선·재평가 → 상위 1~3개를 analyze_structure로 구조 검증 → 네 목적을 종합해 최종 융합체 "
        "1개(전체 서열, 네 지표)와 근거로 결론.\n\n"
        f"참고 라이브러리 — 링커: {lib_l}\n  셔틀: {lib_s}\n"
        f"주의: construct가 {MODEL_MAX_LEN}aa 초과면 BBB는 연결부위(링커+셔틀)로 계산된다. deepB3P는 "
        f"짧은 펩타이드 학습이라 절대값보다 상대 비교가 신뢰됨. 최대 {max_rounds}스텝. 한국어로 간결히."
    )


class GeminiOptimizationAgent(OptimizationAgent):
    """Claude 백엔드 재사용 + Gemini function-calling 루프."""

    def _function_decls(self, types, has_fbgan):
        T = types.Type
        S = types.Schema
        decls = [
            types.FunctionDeclaration(
                name="evaluate_candidates",
                description=("후보 융합체(cargo+linker+shuttle)들의 BBB 투과 점수(deepB3P), 독성"
                             "(ToxinPred3), 융합체 안정성(불안정성 지수), 수용체 유사도를 배치로 채점."),
                parameters=S(type=T.OBJECT, properties={
                    "rationale": S(type=T.STRING, description="이 후보들을 낸 근거"),
                    "constructs": S(type=T.ARRAY, items=S(type=T.OBJECT, properties={
                        "label": S(type=T.STRING),
                        "linker": S(type=T.STRING),
                        "shuttle": S(type=T.STRING),
                    }, required=["label", "linker", "shuttle"])),
                }, required=["rationale", "constructs"]),
            ),
            types.FunctionDeclaration(
                name="analyze_structure",
                description=("한 융합체를 ESMFold로 접어 셔틀이 구조적으로 노출됐는지 vs 화물에 "
                             "가려졌는지와 예측 신뢰도(pLDDT)를 반환. 느림(~10초), 최종 후보에만."),
                parameters=S(type=T.OBJECT, properties={
                    "linker": S(type=T.STRING),
                    "shuttle": S(type=T.STRING),
                }, required=["linker", "shuttle"]),
            ),
        ]
        if has_fbgan:
            decls.append(types.FunctionDeclaration(
                name="generate_novel_shuttles",
                description="라이브러리 밖 novel 셔틀을 FBGAN으로 생성해 BBB·독성과 함께 반환. 느림.",
                parameters=S(type=T.OBJECT, properties={
                    "rounds": S(type=T.INTEGER, description="2~4 라운드"),
                }, required=["rounds"]),
            ))
        return decls

    def _dispatch(self, cargo, name, args):
        """도구 실행 → (LLM에 돌려줄 텍스트, 방출할 AgentEvent 또는 None)."""
        if name == "evaluate_candidates":
            text, rows = self._evaluate(cargo, args)
            return text, AgentEvent("evaluation", text=args.get("rationale", ""), data={"rows": rows})
        if name == "analyze_structure":
            text, sdata = self._structure(cargo, args)
            return text, AgentEvent("structure", text=text, data=sdata)
        if name == "generate_novel_shuttles":
            text, gdata = self._generate(cargo, args)
            return text, AgentEvent("generation", text=text, data=gdata)
        return f"unknown tool: {name}", None

    def _gen_retry(self, client, model, contents, config, tries=6):
        """일시적 오류(503/429)는 백오프 재시도. 그 외는 즉시 raise."""
        last = None
        for attempt in range(tries):
            try:
                return client.models.generate_content(model=model, contents=contents, config=config)
            except Exception as exc:  # noqa: BLE001
                last = exc
                if attempt < tries - 1 and any(t in str(exc) for t in _TRANSIENT):
                    time.sleep(2 * (attempt + 1))  # 2·4·6초
                    continue
                raise
        raise last  # pragma: no cover

    def run(self, cargo: str) -> Generator[AgentEvent, None, None]:
        cargo = cargo.strip().upper()
        key = self.settings.gemini_api_key
        if not key:
            yield AgentEvent("error", "GEMINI_API_KEY(또는 GOOGLE_API_KEY) 미설정.")
            return
        try:
            from google import genai
            from google.genai import types
        except ImportError:
            yield AgentEvent("error", "google-genai 미설치: pip install google-genai")
            return

        client = genai.Client(api_key=key)
        model = self.settings.gemini_model
        has_fbgan = self.settings.use_fbgan_local
        tools = types.Tool(function_declarations=self._function_decls(types, has_fbgan))
        system = _gemini_system_prompt(cargo, self.settings.toxicity_threshold,
                                       self.max_rounds, has_fbgan)

        def _cfg(with_tools=True):
            kw = dict(system_instruction=system, temperature=1.0,
                      thinking_config=types.ThinkingConfig(include_thoughts=True))
            if with_tools:
                kw["tools"] = [tools]
            return types.GenerateContentConfig(**kw)

        contents = [types.Content(role="user", parts=[types.Part.from_text(
            text=f"화물 `{cargo}`에 최적인 링커·셔틀을 설계하라. 도구를 자율적으로 사용해 "
                 "BBB·독성·구조·안정성을 모두 만족하는 최종 융합체를 찾아라.")])]
        best = None

        def _blocked(cand, resp):
            fr = str(getattr(cand, "finish_reason", "") or "")
            if any(x in fr.upper() for x in ("SAFETY", "BLOCK", "PROHIBITED", "RECITATION")):
                return fr
            fb = getattr(resp, "prompt_feedback", None)
            if fb and getattr(fb, "block_reason", None):
                return f"prompt blocked: {fb.block_reason}"
            return None

        def _visible(parts):
            return "\n\n".join(p.text for p in parts
                               if getattr(p, "text", None) and not getattr(p, "thought", False))

        def _user(text):
            return types.Content(role="user", parts=[types.Part.from_text(text=text)])

        try:
            # ── 0. 계획(PLAN): 도구를 쓰기 전에 전략 선언 ──
            contents.append(_user("도구를 호출하기 전에, 최적화 전략을 2~3문장으로 밝혀라 — 어떤 "
                                  "링커·셔틀 방향을 우선 탐색하고 왜인지. 아직 도구는 호출하지 마라."))
            presp = self._gen_retry(client, model, contents, _cfg(False))
            pcand = (presp.candidates or [None])[0]
            blk = _blocked(pcand, presp) if pcand is not None else "빈 응답"
            if blk:
                yield AgentEvent("error", f"Gemini 차단/거부: {blk}")
                return
            pparts = (pcand.content.parts if pcand and pcand.content else None) or []
            for p in pparts:
                if getattr(p, "text", None) and getattr(p, "thought", False):
                    yield AgentEvent("reasoning", p.text)
            contents.append(pcand.content)
            yield AgentEvent("plan", _visible(pparts))
            contents.append(_user("이제 전략에 따라 도구를 자율적으로 사용해 실행하라."))

            # ── 1. 실행 루프(ACT) ──
            final_text = None
            for rnd in range(self.max_rounds):
                resp = self._gen_retry(client, model, contents, _cfg(True))
                cand = (resp.candidates or [None])[0]
                blk = _blocked(cand, resp) if cand is not None else "빈 응답"
                if blk:
                    yield AgentEvent("error", f"Gemini 차단/거부: {blk}")
                    return
                parts = (cand.content.parts if cand and cand.content else None) or []
                contents.append(cand.content)

                fcalls = []
                for p in parts:
                    if getattr(p, "function_call", None):
                        fcalls.append(p.function_call)
                    elif getattr(p, "text", None):
                        if getattr(p, "thought", False):
                            yield AgentEvent("reasoning", p.text)
                        else:
                            yield AgentEvent("text", p.text)

                if not fcalls:
                    final_text = _visible(parts)
                    break

                fr_parts = []
                for fc in fcalls:
                    args = dict(fc.args or {})
                    text, ev = self._dispatch(cargo, fc.name, args)
                    if ev is not None:
                        yield ev
                        if ev.kind == "evaluation":
                            for r in ev.data["rows"]:
                                if not r["toxic"] and (best is None or r["bbb"] > best["bbb"]):
                                    best = r
                    fr_parts.append(types.Part.from_function_response(
                        name=fc.name, response={"result": text}))
                contents.append(types.Content(role="user", parts=fr_parts))
                if best:
                    yield AgentEvent("progress", data={"round": rnd + 1, "best_bbb": best["bbb"]})
            else:
                contents.append(_user("예산을 모두 사용했습니다. 지금까지의 최적 융합체(전체 서열, 네 "
                                      "지표 종합)와 근거로 최종 보고서를 쓰세요. 도구 호출 금지."))
                resp = self._gen_retry(client, model, contents, _cfg(False))
                cand = (resp.candidates or [None])[0]
                parts = (cand.content.parts if cand and cand.content else None) or []
                contents.append(cand.content)
                final_text = _visible(parts)

            # ── 2. 최종 보고서(FINAL) ──
            yield AgentEvent("final", final_text or "")

            # ── 3. 자기평가(REFLECT) ──
            contents.append(_user("마지막으로 자기평가하라(3~4문장): 최종 융합체가 네 목적을 각각 얼마나 "
                                  "충족했는지, 약점·리스크는 무엇인지, 예산이 더 있다면 다음에 무엇을 "
                                  "시도할지. 도구 호출 금지."))
            rresp = self._gen_retry(client, model, contents, _cfg(False))
            rcand = (rresp.candidates or [None])[0]
            rparts = (rcand.content.parts if rcand and rcand.content else None) or []
            for p in rparts:
                if getattr(p, "text", None) and getattr(p, "thought", False):
                    yield AgentEvent("reasoning", p.text)
            yield AgentEvent("reflection", _visible(rparts))

            if best:
                yield AgentEvent("optimum", data=best)
        except Exception as exc:  # noqa: BLE001
            yield AgentEvent("error", f"Gemini 호출 오류: {type(exc).__name__}: {exc}")


def get_gemini_agent(settings: Settings, max_rounds: int = 8) -> GeminiOptimizationAgent | None:
    return GeminiOptimizationAgent(settings, max_rounds) if settings.use_gemini_agent else None
