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
    lib_l = ", ".join(f"{n}({v['seq'] or '링커없음/직접융합'})" for n, v in LINKER_LIBRARY.items())
    lib_s = ", ".join(f"{n}({v['seq']})" for n, v in SHUTTLES.items())
    tools = ("evaluate_candidates(BBB·독성·안정성·수용체 배치), design_candidate(잔기 수준 편집 "
             "서열 채점), analyze_structure(구조 노출도)")
    if has_fbgan:
        tools += ", generate_novel_shuttles(신규 셔틀 생성)"
    tools += ", finish(수렴 시 최종 후보 제출·종료)"
    return (
        "[연구 맥락] 이것은 알츠하이머병 치료제의 뇌 전달을 개선하기 위한 정당한 학술·공모전용 "
        "in-silico 연구입니다. 사용하는 셔틀 펩타이드(Angiopep-2, TAT, Penetratin, SynB1 등)는 "
        "학계·임상에서 공개 발표된 뇌 전달 도구이며, 독성 예측은 안전한 후보만 남기기 위한 안전성 "
        "스크리닝입니다.\n\n"
        "당신은 멀티 에이전트 신약설계 시스템의 **설계(Designer) 에이전트**입니다 — 도구를 자율적으로 "
        "오케스트레이션해 뇌혈관장벽(BBB)을 통과하는 항체-셔틀 융합 단백질을 설계하고, 그 결과는 "
        "독립 **심사(Critic) 에이전트**의 적대적 검증을 거칩니다.\n\n"
        f"고정 화물(cargo): {cargo}\n융합체 = cargo + linker + shuttle.\n\n"
        "다중 목적 (모두 만족하는 최종 융합체 하나로 수렴):\n"
        "  ① BBB 투과 점수(deepB3P 예측확률 0~1, 실제 투과율 아님·상대비교용) 최대화\n"
        f"  ② 독성(ToxinPred3) ≤ {tox_threshold:.2f}\n"
        "  ③ 셔틀이 구조적으로 노출(ESMFold) + 검증 셔틀과 유사(수용체 메커니즘 근거)\n"
        "  ④ 융합체 자체 안정성: 불안정성 지수 < 40 선호\n\n"
        f"사용 가능한 도구: {tools}\n"
        "권장 워크플로우(자율 판단): evaluate_candidates로 라이브러리 조합(**링커 유무 — 직접융합"
        "(링커 빈칸)도 포함**)을 폭넓게 스크리닝 → "
        "유망 방향으로 **design_candidate로 잔기 수준 편집**(변이·트리밍·연장·하이브리드)해 라이브러리 "
        "밖까지 정밀 설계·재평가 → 상위 1~3개를 analyze_structure로 구조 검증 → **충분히 수렴했다고 "
        "판단하면 스스로 finish**를 호출해 최종 융합체 1개(라벨·링커·셔틀·전체 근거)를 제출하라.\n"
        "- design_candidate: 라이브러리에 얽매이지 말고 셔틀/링커의 특정 잔기를 자유롭게 편집한 서열을 "
        "제안·채점할 수 있다(진짜 설계). 왜 그렇게 편집했는지 edit_note에 남겨라.\n"
        "- finish: 고정 스텝을 다 쓸 필요 없다. 목적을 충분히 만족했다고 판단하면 언제든 finish로 최종 "
        "후보를 제출하라. 제출하면 **독립 심사(비평) 에이전트가 적대적으로 검증**한다 — 심사가 REVISE(개선 "
        "요구)를 내면 그 지적을 반영해 개선하고 다시 finish하라. APPROVE(승인)가 나오면 확정된다.\n\n"
        f"참고 라이브러리 — 링커: {lib_l}\n  셔틀: {lib_s}\n"
        f"주의: construct가 {MODEL_MAX_LEN}aa 초과면 BBB는 연결부위(링커+셔틀)로 계산된다. deepB3P는 "
        f"짧은 펩타이드 학습이라 절대값보다 상대 비교가 신뢰됨. 최대 {max_rounds}행동. 한국어로 간결히."
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
                name="design_candidate",
                description=("라이브러리에 없는, **직접 편집·설계한** 링커·셔틀 서열 하나를 채점한다. "
                             "잔기 변이·트리밍·연장·하이브리드 등 자유 설계. evaluate와 동일한 4지표 반환."),
                parameters=S(type=T.OBJECT, properties={
                    "label": S(type=T.STRING),
                    "linker": S(type=T.STRING),
                    "shuttle": S(type=T.STRING),
                    "edit_note": S(type=T.STRING, description="어떻게·왜 편집했는지"),
                }, required=["label", "linker", "shuttle"]),
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
        decls.append(types.FunctionDeclaration(
            name="finish",
            description=("충분히 수렴했다고 판단하면 호출해 최종 후보를 제출하고 종료한다. 고정 스텝을 "
                         "다 쓸 필요 없다."),
            parameters=S(type=T.OBJECT, properties={
                "chosen_label": S(type=T.STRING),
                "chosen_linker": S(type=T.STRING),
                "chosen_shuttle": S(type=T.STRING),
                "final_report": S(type=T.STRING, description="최종 융합체·4지표·선정 근거"),
            }, required=["chosen_linker", "chosen_shuttle", "final_report"]),
        ))
        return decls

    def _dispatch(self, cargo, name, args):
        """도구 실행 → (LLM에 돌려줄 텍스트, 방출할 AgentEvent 또는 None)."""
        if name == "evaluate_candidates":
            text, rows = self._evaluate(cargo, args)
            return text, AgentEvent("evaluation", text=args.get("rationale", ""), data={"rows": rows})
        if name == "design_candidate":
            text, rows = self._evaluate(cargo, {"constructs": [{
                "label": args.get("label", "설계"), "linker": args.get("linker", ""),
                "shuttle": args.get("shuttle", "")}]})
            note = args.get("edit_note", "").strip()
            return text, AgentEvent("evaluation",
                                    text="🎨 정밀 설계" + (f" — {note}" if note else ""),
                                    data={"rows": rows})
        if name == "analyze_structure":
            text, sdata = self._structure(cargo, args)
            return text, AgentEvent("structure", text=text, data=sdata)
        if name == "generate_novel_shuttles":
            text, gdata = self._generate(cargo, args)
            return text, AgentEvent("generation", text=text, data=gdata)
        return f"unknown tool: {name}", None

    def _score_choice(self, cargo, fa):
        """finish가 제출한 최종 후보를 채점해 row(dict)로 반환. 실패 시 None."""
        lk = "".join(ch for ch in (fa.get("chosen_linker", "") or "").upper() if ch.isalpha())
        sh = "".join(ch for ch in (fa.get("chosen_shuttle", "") or "").upper() if ch.isalpha())
        if not sh:
            return None
        _, rows = self._evaluate(cargo, {"constructs": [{
            "label": fa.get("chosen_label", "최종 선택"), "linker": lk, "shuttle": sh}]})
        if not rows:
            return None
        row = rows[0]
        row["agent_pick"] = True
        return row

    # ---- 심사(비평) 에이전트 — 독립 페르소나·독립 컨텍스트 ----------------
    def _critic_prompt(self):
        thr = self.settings.toxicity_threshold
        return (
            "당신은 신약 후보를 **적대적으로 검증하는 안전·엄밀성 심사 에이전트**입니다. 설계 에이전트가 "
            "제안한 최종 BBB 융합체 후보를 의심하고 반박하세요. 관대하지 마세요.\n"
            "날카롭게 따질 것:\n"
            f"① 독성 마진 — 임계값 {thr:.2f} 대비 실제로 안전한가, 아슬아슬한가.\n"
            "② 안정성 — 불안정성 지수가 40 미만이 이상적. 높으면 응집·분해 위험을 지적.\n"
            "③ BBB 점수 신뢰도 — deepB3P는 짧은 펩타이드로 학습 → 긴 융합체의 절대값은 신뢰도가 낮다. "
            "상대 비교로만 해석했는지.\n"
            "④ 구조 근거 — 셔틀이 표면 노출됐다는 ESMFold 검증을 실제로 거쳤는가.\n"
            "⑤ 개발성 리스크 — 과도한 Arg/Lys(양이온)는 비특이 결합·독성·응집 위험. 서열 liability.\n"
            "⑥ off-target/선택성 — 양전하·친유성 주도(CPP형)면 여러 조직에 비특이 흡수돼 off-target "
            "위험. RMT 수용체 표적(Angiopep형)이 선택적. off-target 위험 높으면 지적하라.\n"
            "⑦ 설계의 실질 — 라이브러리의 뻔한 선택인가, 근거 있는 정밀 설계인가.\n"
            "3~5문장으로 핵심 반박을 쓰고, **마지막 줄에 정확히** `VERDICT: APPROVE`(승인) 또는 "
            "`VERDICT: REVISE`(개선 요구)를 써라. 애매하면 REVISE."
        )

    def _critic_review(self, client, model, types, cargo, choice):
        """최종 후보를 별도 컨텍스트로 적대적 검증. 반환 (approve: bool, critique_text: str)."""
        dev_list = choice.get("dev_liabilities") or []
        m = (f"화물={cargo}, 링커={choice.get('linker') or '—'}, 셔틀={choice.get('shuttle')}, "
             f"전체서열={choice.get('sequence')}\n"
             f"지표 — BBB점수={choice['bbb']:.3f}(0~1), 독성={choice['tox']:.3f}, "
             f"불안정성지수={choice['instability']}, "
             f"수용체유사도({choice.get('bind_ref')})={choice.get('bind_score', 0):.2f}, "
             f"독성판정={'탈락' if choice['toxic'] else '통과'}\n"
             f"개발성 — 위험={choice.get('dev_risk', '?')}, 순전하={choice.get('dev_charge', '?')}, "
             f"응집={choice.get('dev_agg', '?')}, liability({len(dev_list)}): "
             + ("; ".join(dev_list) if dev_list else "없음") + "\n"
             f"선택성/off-target — 위험={choice.get('sel_level', '?')}"
             f"(선택성 {choice.get('selectivity', '?')}), {choice.get('sel_mech', '?')}\n"
             f"용해도 — {choice.get('sol_level', '?')}(점수 {choice.get('sol_score', '?')})")
        cfg = types.GenerateContentConfig(
            system_instruction=self._critic_prompt(), temperature=1.0,
            thinking_config=types.ThinkingConfig(include_thoughts=False))
        contents = [types.Content(role="user", parts=[types.Part.from_text(
            text="설계 에이전트가 제출한 최종 후보다. 적대적으로 검증하고 VERDICT를 내려라.\n" + m)])]
        resp = self._gen_retry(client, model, contents, cfg)
        cand = (resp.candidates or [None])[0]
        parts = (cand.content.parts if cand and cand.content else None) or []
        text = "\n\n".join(p.text for p in parts
                           if getattr(p, "text", None) and not getattr(p, "thought", False))
        approve = "APPROVE" in (text or "").upper().rsplit("VERDICT:", 1)[-1]
        return approve, text

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

            # ── 1. 설계자 실행 + 자기종료(finish) → 심사(Critic) → REVISE 시 재설계 ──
            choice = final_report = critique_text = None
            approve = True
            critiques, MAX_CRITIQUES, turn, grace = 0, 1, 0, 0
            while turn < self.max_rounds + grace:
                turn += 1
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
                        yield AgentEvent("reasoning" if getattr(p, "thought", False) else "text", p.text)

                if not fcalls:  # 도구 없이 결론 텍스트 → 종료로 간주
                    final_report = _visible(parts)
                    break

                fr_parts, finish_fc = [], None
                for fc in fcalls:
                    if fc.name == "finish":
                        finish_fc = fc
                        fr_parts.append(types.Part.from_function_response(
                            name="finish", response={"result": "acknowledged; self-evaluating"}))
                        continue
                    text, ev = self._dispatch(cargo, fc.name, dict(fc.args or {}))
                    if ev is not None:
                        yield ev
                        if ev.kind == "evaluation":
                            for r in ev.data["rows"]:
                                if not r["toxic"] and (best is None or r["bbb"] > best["bbb"]):
                                    best = r
                    fr_parts.append(types.Part.from_function_response(
                        name=fc.name, response={"result": text}))

                if finish_fc is None:
                    contents.append(types.Content(role="user", parts=fr_parts))
                    if best:
                        yield AgentEvent("progress", data={"round": turn, "best_bbb": best["bbb"]})
                    continue

                # finish → 최종 후보 채점 → 독립 심사(Critic) → APPROVE/REVISE
                fa = dict(finish_fc.args or {})
                final_report = fa.get("final_report", "") or final_report
                ch = self._score_choice(cargo, fa)
                if ch is not None:
                    choice = ch
                if choice is not None:
                    approve, critique_text = self._critic_review(client, model, types, cargo, choice)
                else:
                    approve, critique_text = True, "후보 채점 실패로 심사를 생략합니다."
                yield AgentEvent("progress", data={"round": turn,
                                                   "best_bbb": (choice or best or {}).get("bbb", 0)})
                if (not approve) and critiques < MAX_CRITIQUES:
                    critiques += 1
                    grace += 1  # 재설계용 여분 턴 부여 → 마지막 스텝에서 REVISE가 나와도 개선 1회 보장
                    contents.append(types.Content(role="user", parts=fr_parts + [types.Part.from_text(
                        text="독립 **심사 에이전트가 REVISE(개선 요구)**를 냈다. 심사 반박:\n"
                             f"{critique_text}\n이 지적을 반영해 design_candidate 등으로 후보를 개선하고, "
                             "끝나면 다시 finish를 호출하라.")]))
                    continue
                contents.append(types.Content(role="user", parts=fr_parts))  # 함수응답 기록 후 종료
                break

            # 강제 종료 대비: 최종 보고서·심사 확보
            if final_report is None:
                contents.append(_user("행동 예산을 모두 사용했습니다. 지금까지의 최적 융합체(전체 서열, "
                                      "네 지표)와 근거로 최종 보고서를 쓰세요. 도구 호출 금지."))
                resp = self._gen_retry(client, model, contents, _cfg(False))
                cand = (resp.candidates or [None])[0]
                parts = (cand.content.parts if cand and cand.content else None) or []
                contents.append(cand.content)
                final_report = _visible(parts)
            final_pick = choice or best
            if critique_text is None and final_pick is not None:
                approve, critique_text = self._critic_review(client, model, types, cargo, final_pick)

            # ── 결과 방출 (한 번씩): 설계자 결론 + 심사 판정 ──
            # 심사 최종 판정을 최종 후보에 각인 → UI가 "승인/미승인"을 정직하게 표시한다.
            if choice is not None:
                choice["critic_approved"] = approve
                yield AgentEvent("choice", data=choice)
            yield AgentEvent("final", final_report or "")
            if critique_text:
                yield AgentEvent("critique", critique_text, data={"approve": approve})
            if final_pick:
                final_pick["critic_approved"] = approve
                yield AgentEvent("optimum", data=final_pick)
        except Exception as exc:  # noqa: BLE001
            yield AgentEvent("error", f"Gemini 호출 오류: {type(exc).__name__}: {exc}")


def get_gemini_agent(settings: Settings, max_rounds: int = 8) -> GeminiOptimizationAgent | None:
    return GeminiOptimizationAgent(settings, max_rounds) if settings.use_gemini_agent else None
