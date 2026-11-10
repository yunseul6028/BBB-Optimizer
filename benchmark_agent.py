"""에이전트 벤치마크 — LLM 필요 지표 (M4 효율/최적성 갭, M5 심사 효과).

benchmark.py(M1~M3, 로컬)와 짝. 이건 **LLM 호출이 필요**하므로 쿼터/크레딧이 있을 때만.
기본 브레인은 Gemini(설정값). 실행:  python benchmark_agent.py

M4) 에이전트 효율/최적성 갭:
    에이전트를 1회 돌려 (a) 평가한 후보 수(탐색 예산)와 (b) 최종 선택이 **도달 가능
    전선에서 어디인지**를 집계. ⚠️ 옛 `최종 BBB / oracle BBB` 백분율은 폐지(2026-11-10):
    BBB 로는 후보를 줄 세울 수 없고(M6·M8), 그 분모가 근거 등급 D(Tet1)였다.
M5) 심사(Critic) 효과:
    일부러 나쁜 후보(독성·근거 미등재·등급 D·CPP·퇴행)와 좋은 후보를 심사 에이전트에
    직접 넣어 REVISE/APPROVE 판정이 맞는지 + **의도한 게이트가 걸렸는지** 확인.
"""

from __future__ import annotations

from core.config import TOXICITY_THRESHOLD, get_settings
from core.predictors import get_predictor
from core.toxicity import get_toxicity_predictor
from core.optimizer_agent_gemini import get_gemini_agent
from benchmark import _library_scores, is_cpp_row


# max_rounds=6 은 의도적 선택 — "타이트한 예산에서도 에이전트가 효율적인가"를 보는 실측이다.
# (운영 기본값은 10 = 안전 상한, get_gemini_agent 참고. 여기 6은 수렴 바닥을 재는 프로브다.)
def m4_agent_efficiency(cargo="GSNKGAIIGLM", max_rounds=6):
    """에이전트 효율/최적성 갭.

    ⚠️ 2026-11-10 재정의: 옛 판정은 **`최종 BBB / 후보군 oracle BBB`** 였다. 폐지했다 —
    분모와 분자가 모두 BBB 이고, BBB 로는 후보를 줄 세울 수 없음이 두 번 확인됐다
    (M6 셔틀 간 · M8 링커 간 0/108). 게다가 그 분모(후보군 oracle = `직접융합+Tet1`,
    BBB 0.568)는 **근거 등급 D**(BBB 통과 단독 근거 없음)여서, '최적'이 문헌상 가장
    근거 없는 셔틀이었다(§15-4).

    새 판정은 **도달 가능 전선 대비 위치**다 — 백분율이 아니라 집합 소속이다:
      ① 최종 후보가 **달성 가능한 최상위 근거 그룹**인가
      ② 전선(`front_top`)의 어느 후보에게도 **지배당하지 않는가**
      ③ 전선 안에 실제로 들어갔나 / 전선 밖을 **설계해서** 비지배로 만들었나
    ③이 이 시스템이 기여를 주장할 수 있는 유일한 자리다(부모 근거를 유지한 liability 스크러빙).
    """
    settings = get_settings()
    agent = get_gemini_agent(settings, max_rounds=max_rounds, use_ledger=False)
    if agent is None:
        print("⚠️ Gemini 에이전트 비활성(GEMINI_API_KEY 필요) — M4 생략.")
        return
    from core.evidence import group_of, group_rank
    from core.optimizer_agent import dominates, is_candidate, rank_key

    front = agent.reachable_front(cargo)
    ft = front.get("front_top") or []
    fr = front.get("front") or []
    print("── M4. 에이전트 효율 (라이브 실행) ──")
    print(f"  도달 가능 집합: 후보군 {front.get('n_gate')}개 · 근거 그룹 "
          f"{front.get('by_group')} · 최상위 그룹 `{front.get('best_group')}`")
    print(f"  최상위 그룹 전선 {len(ft)}개 / 전체 전선 {len(fr)}개")

    n_eval, choice, best_seen = 0, None, None
    for ev in agent.run(cargo):
        if ev.kind == "evaluation":
            n_eval += len(ev.data.get("rows", []))
            for r in ev.data["rows"]:
                if is_candidate(r) and (best_seen is None or rank_key(r) < rank_key(best_seen)):
                    best_seen = r
        elif ev.kind in ("choice", "optimum"):
            choice = ev.data
        elif ev.kind == "error":
            print(f"  에이전트 오류: {ev.text}")
            return
    pick = choice if choice is not None else best_seen
    if pick is None:
        print("  ⚠️ 최종 후보 없음 — 모든 후보가 게이트에 걸렸다.")
        return

    g = group_of(pick.get("tier"))
    on_top = group_rank(pick.get("tier")) <= group_rank((ft[0].get("tier") if ft else None))
    dom_by = [r["label"] for r in ft if dominates(r, pick)]
    labels = {r["label"] for r in fr}
    in_front = f"{pick.get('linker') or '직접융합'}+{pick.get('shuttle', '?')}" in labels

    # 재현을 위해 **무엇을 골랐는지** 서열까지 남긴다 — 설계체면 라벨만으론 복원이 안 된다.
    print(f"  에이전트 최종 선택 서열: 링커 `{pick.get('linker') or '직접융합'}` + 셔틀 "
          f"`{pick.get('shuttle')}`  (근거 상속: 부모 {pick.get('tier_parent')}, "
          f"유사도 {pick.get('tier_sim')})")
    print(f"  에이전트 최종 선택: 근거 **{pick.get('tier')}**({g}) · "
          f"liab {pick.get('dev_liab')} · 불안정성 {pick.get('instability')} · "
          f"off {pick.get('sel_off')} · 응집핫스팟 {pick.get('dev_agg_hot')} "
          f"[BBB 참고 {pick.get('bbb', 0):.3f} — 순위에 쓰지 않음]")
    _top_g = front.get("best_group")
    print("  ① 최상위 근거 그룹인가 : "
          + ("예" if on_top else f"아니오 (최상위 = {_top_g})"))
    print(f"  ② 전선에 지배당하나   : {'아니오(비지배)' if not dom_by else '예 — ' + ', '.join(dom_by[:3])}")
    print(f"  ③ 전선 안인가         : {'라이브러리 전선 내' if in_front else '라이브러리 밖(설계체)'}")
    print(f"  에이전트가 평가한 후보 수: {n_eval}  (라이브러리 전수 {front.get('n_total')}개 중)")
    if not on_top:
        verdict = "근거 하향 — 더 단단한 근거의 대안이 있었다"
    elif dom_by:
        verdict = "지배당함 — 모든 축에서 더 나은 달성 가능 후보가 있었다"
    elif not in_front:
        verdict = "전선 초과(설계로 라이브러리 밖에서 비지배 달성)"
    else:
        verdict = "전선 도달(라이브러리 내 비지배)"
    print(f"  판정: {verdict}")


def m5_critic_effectiveness():
    settings = get_settings()
    agent = get_gemini_agent(settings, use_ledger=False)
    if agent is None:
        print("⚠️ Gemini 에이전트 비활성 — M5 생략.")
        return
    from google import genai
    from google.genai import types
    client = genai.Client(api_key=settings.gemini_api_key)
    model = settings.gemini_model

    # 나쁜 후보 5종(각기 **다른 게이트**를 때린다) + 좋은 후보 1종.
    # ⚠️ 2026-11-10 개편: 게이트 체인이 바뀌었다. 독성 → **근거 등급 보유** → 등급 D →
    #    CPP 비교군 → 근거 하향 → 지배 → 래칫. 옛 프로브가 겨냥했던 'BBB 하한'·'순서 특이성'·
    #    '전선 이탈'은 전부 폐지됐고, 그 역할은 **근거 등급 보유 게이트**가 승계했다.
    CARGO = "GSNKGAIIGLM"
    probes = [
        # (라벨, 링커, 셔틀, 기대, 겨냥한 게이트, best_seen 로 줄 더 좋은 후보)
        ("독성 셔틀(멜리틴)", "GGGGS", "GIGAVLKVLTTGLPALISWIKRKRQQ", "REVISE", "독성", None),
        # 옛 "저BBB"·"순서 측정불가" 프로브의 후계. poly-Ser 는 어떤 라이브러리 셔틀과도
        # 유사도가 상속 하한(0.70)에 못 미친다 → **근거 미등재**로 걸린다.
        ("근거 미등재(poly-Ser)", "GGGGS", "SSSSSSSSSSSS", "REVISE", "근거 미등재", None),
        # 옛 "조성쌓기 설계체" 프로브의 후계 — 에이전트가 실제로 올렸던 ApoE 변이(R 치환).
        # ⚠️ 실측(로컬): 부모 ApoE(159-167)2 와 유사도 **0.767 ≥ 0.70** 이라 등급 B 를
        #    **상속한다**. 상속 게이트는 모티프를 보존한 과편집을 막지 않는다 — 의도된 설계다
        #    (merit 의 출처는 부모 모티프이고, '얼마나 편집했나'가 아니라 **위험 축**이 심사한다).
        #    실제로 이 후보는 off 0.36(vs 전선 0.14)·liability 3(vs 1)로 **전선에 지배당한다.**
        #    그래서 이 프로브가 겨냥하는 게이트는 '근거 미등재'가 아니라 **'지배'** 다.
        ("과편집 설계체(ApoE변이)", "", "LRVYQAGRRLRVYQAGRR", "REVISE", "지배", None),
        # 등급 D 게이트 — Tet1 은 라이브러리 셔틀이라 상속은 되지만 물려받은 근거가 '없음'이다.
        # 우리 효능 축에서는 후보군 1위였던 조합(M2·M4 oracle)이라 **가장 중요한 프로브**다.
        ("근거 없음(등급 D·Tet1)", "", "HLNILSTLWKYR", "REVISE", "근거 등급 D", None),
        ("비교군(SynB1·CPP)", "EAAAKEAAAK", "RGGRLSYSRRRFSTSTGR", "REVISE", "CPP 비교군", None),
        # 근거 하향 게이트 — B6 은 등급 **C(주의)** 인데 최상위 그룹 '근거 있음'(A·B)이
        # 달성 가능하다. 등급 D 도 CPP 도 아니므로 앞선 게이트를 전부 지나 여기서 걸린다.
        ("근거 하향(등급 C·B6)", "AEAAAKEAAAKEAAAKA", "GHKAKGPRK", "REVISE", "근거 하향", None),
        # 래칫 프로브 — 앞선 게이트를 **전부** 통과해야 래칫을 때린다. 그래서 제출 후보를
        # **전선 멤버**(GGS+T7, 지배당하지 않음)로 두고, best_seen 만 그것을 지배하는
        # **설계체**(개발성 스크러빙 변이 HAVYPNH — T7 등급 B 상속, off 0.14→0.03·응집
        # 0.59→0.58, liability·불안정성 동일)로 준다. 로컬 실측으로 지배 성립을 확인했다.
        # ⚠️ 바로 아래 '양호' 프로브는 **같은 후보에 best_seen 만 없는 것**이다 —
        #    둘의 판정 차이가 곧 래칫 게이트의 순효과다(다른 변수가 없다).
        ("퇴행(더 좋은 걸 봤음)", "GGS", "HAIYPRH", "REVISE", "래칫", ("GGS", "HAVYPNH")),
        ("양호(GGS+T7, 전선 멤버)", "GGS", "HAIYPRH", "APPROVE?", "—", None),
    ]

    # 어느 게이트가 실제로 걸렸나 — 심사문의 지표 문구로 식별한다.
    # ⚠️ 이게 없으면 "REVISE 나왔다"만 알고 **의도한 게이트가 걸렸는지는 모른다**.
    #    실제로 예전 래칫 프로브는 BBB 하한에 먼저 걸려 있었는데 통과로 집계됐다.
    GATE_MARKERS = [
        ("독성", "안전성 하드 게이트"),
        ("근거 미등재", "문헌 근거 미등재"),
        ("근거 등급 D", "BBB 통과 단독 근거 없음)**이다"),
        ("CPP 비교군", "비특이 흡수 비교군"),
        ("근거 하향", "근거 하향"),
        # ⚠️ 두 게이트의 문구가 **"모든 축에서 지배당한다"를 공유**한다 → 그 구절을 마커로
        #    쓰면 래칫 메시지가 '지배'로도 잡혀 `지배+래칫` 으로 찍힌다(실측 2026-11-11).
        #    그래서 각 게이트의 **고유 머리말**을 마커로 쓴다. 어느 쪽도 다른 쪽에 없다.
        ("지배", "제출 후보가 전선의 **`"),
        ("래칫", "이미 채점해 본"),
    ]

    def fired_gate(txt: str) -> str:
        hit = [nm for nm, mk in GATE_MARKERS if mk in (txt or "")]
        return "+".join(hit) if hit else "—"
    print("\n── M5. 심사(Critic) 효과 (라이브) ──")
    hits = n_bad = 0
    mis_gate = []   # 의도와 다른 게이트에 걸린 프로브
    for label, lk, sh, expect, gate, better in probes:
        choice = agent._score_choice(CARGO, {"chosen_label": label,
                                             "chosen_linker": lk, "chosen_shuttle": sh})
        if choice is None:
            print(f"  {label}: 채점 실패")
            continue
        best_seen = None
        if better is not None:
            best_seen = agent._score_choice(CARGO, {"chosen_label": "본최고",
                                                    "chosen_linker": better[0],
                                                    "chosen_shuttle": better[1]})
        approve, critique, tier = agent._critic_review(client, model, types, CARGO, choice,
                                                       best_seen=best_seen)
        got = fired_gate(critique)
        ok_gate = (gate == "—") or (gate in got.split("+"))
        if expect.startswith("REVISE"):
            n_bad += 1
            # **의도한 게이트가 걸렸을 때만** 검출로 집계한다. REVISE 자체는 쉽다.
            hits += (not approve) and ok_gate
            if not ok_gate:
                mis_gate.append(f"{label}(기대 {gate} / 실제 {got})")
        print(f"  {label:26s} [겨냥 {gate:9s}] BBB={choice['bbb']:.2f} 독성={choice['tox']:.2f} "
              f"II={choice['instability']} → {tier}  실제게이트 {got}  "
              f"{'✅' if (tier.startswith('REVISE') == expect.startswith('REVISE')) and ok_gate else '❌'}")
    print(f"  나쁜 후보 검출(의도한 게이트까지 일치): {hits}/{n_bad}")
    if mis_gate:
        print("  ⚠️ 엉뚱한 게이트에 걸린 프로브: " + "; ".join(mis_gate))
        print("     → REVISE는 났지만 **겨냥한 게이트를 검증하지 못했다**. 프로브를 고쳐야 한다.")


def main():
    print("=" * 68)
    print("BBB-Optimizer 에이전트 벤치마크 — M4/M5 (LLM 필요)")
    print("=" * 68)
    m4_agent_efficiency()
    m5_critic_effectiveness()


if __name__ == "__main__":
    main()
