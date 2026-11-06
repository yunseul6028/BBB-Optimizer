"""에이전트 벤치마크 — LLM 필요 지표 (M4 효율/최적성 갭, M5 심사 효과).

benchmark.py(M1~M3, 로컬)와 짝. 이건 **LLM 호출이 필요**하므로 쿼터/크레딧이 있을 때만.
기본 브레인은 Gemini(설정값). 실행:  python benchmark_agent.py

M4) 에이전트 효율/최적성 갭:
    에이전트를 1회 돌려 (a) 평가한 후보 수(탐색 예산)와 (b) 최종 선택 BBB를 집계.
    라이브러리 oracle 최적(benchmark._library_scores)과 무작위 탐색 기준선에 대비.
    → 에이전트가 더 적은 평가로 최적에 근접/초과하면 "에이전트의 가치" 정량 증명.
M5) 심사(Critic) 효과:
    일부러 나쁜 후보(독성·불안정·저BBB)와 좋은 후보를 심사 에이전트에 직접 넣어
    REVISE/APPROVE 판정이 맞는지 확인. → 심사 에이전트가 실제로 걸러내는가.
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
    settings = get_settings()
    agent = get_gemini_agent(settings, max_rounds=max_rounds)
    if agent is None:
        print("⚠️ Gemini 에이전트 비활성(GEMINI_API_KEY 필요) — M4 생략.")
        return
    predictor = get_predictor(settings)
    tox_pred = get_toxicity_predictor(settings)
    rows = _library_scores(cargo, predictor, tox_pred, TOXICITY_THRESHOLD)
    valid = [r for r in rows if not r["toxic"]]
    # oracle 분모는 **후보군(RMT)** 기준이어야 한다. 비독성 전수의 최고는 CPP(SynB1)인데
    # CPP는 비교군이라 에이전트가 고를 수 없다 → 도달 불가능한 기준으로 나누면 항상 미달로 보인다.
    # (같은 오류를 애블레이션 집계에서도 고쳤다 — docs/experiments.md §11.)
    cand = [r for r in valid if not is_cpp_row(r)]
    oracle = max(cand, key=lambda r: r["bbb"])["bbb"] if cand else 0.0
    oracle_all = max(valid, key=lambda r: r["bbb"])["bbb"] if valid else 0.0

    n_eval, choice_bbb, best_seen = 0, None, 0.0
    print("── M4. 에이전트 효율 (라이브 실행) ──")
    for ev in agent.run(cargo):
        if ev.kind == "evaluation":
            n_eval += len(ev.data.get("rows", []))
            for r in ev.data["rows"]:
                if not r["toxic"]:
                    best_seen = max(best_seen, r["bbb"])
        elif ev.kind in ("choice", "optimum"):
            choice_bbb = ev.data.get("bbb", choice_bbb)
        elif ev.kind == "error":
            print(f"  에이전트 오류: {ev.text}")
            return
    final = choice_bbb if choice_bbb is not None else best_seen
    print(f"  후보군(RMT) oracle 최적 BBB: {oracle:.3f}   "
          f"[참고: CPP 포함 전수 {oracle_all:.3f} — 비교군이라 선택 불가]")
    print(f"  에이전트 최종 선택 BBB     : {final:.3f}  (최적 대비 {100*final/oracle:.1f}%)")
    print(f"  에이전트가 평가한 후보 수  : {n_eval}  (전수 {len(rows)}개 중)")
    verdict = "최적 초과(설계 도구로 라이브러리 밖 발견)" if final > oracle + 1e-6 else \
              "최적 근접" if final >= 0.9 * oracle else "최적 미달"
    print(f"  판정: {verdict}")


def m5_critic_effectiveness():
    settings = get_settings()
    agent = get_gemini_agent(settings)
    if agent is None:
        print("⚠️ Gemini 에이전트 비활성 — M5 생략.")
        return
    from google import genai
    from google.genai import types
    client = genai.Client(api_key=settings.gemini_api_key)
    model = settings.gemini_model

    # 나쁜 후보 4종(각기 **다른 게이트**를 때린다) + 좋은 후보 1종.
    # ⚠️ 예전 '좋은 후보'는 SynB1(CPP)이었는데, CPP는 이제 비교군 게이트로 자동 미승인된다
    #    → 항상 실패하는 프로브였다. 후보군(RMT) 최고 조합인 T7 직접융합으로 교체했다.
    CARGO = "GSNKGAIIGLM"
    probes = [
        # (라벨, 링커, 셔틀, 기대, 겨냥한 게이트, best_seen 로 줄 더 좋은 후보)
        ("독성 셔틀(멜리틴)", "GGGGS", "GIGAVLKVLTTGLPALISWIKRKRQQ", "REVISE", "독성", None),
        ("저BBB(poly-Ser)", "GGGGS", "SSSSSSSSSSSS", "REVISE", "BBB 하한", None),
        ("비교군(SynB1·CPP)", "EAAAKEAAAK", "RGGRLSYSRRRFSTSTGR", "REVISE", "CPP 비교군", None),
        # 래칫 프로브는 **앞선 게이트를 모두 통과**해야 래칫을 때린다. BBB 하한(0.15)과
        # 전선 하한(≈0.42)을 넘으면서 best_seen 보다만 낮아야 하므로, T7 직접융합(eff 0.491)을
        # 제출하고 더 좋은 `A+T7`(eff 0.635)을 이미 본 것으로 준다. 차이 0.144 > 여유폭 0.02.
        ("퇴행(더 좋은 걸 봤음)", "", "HAIYPRH", "REVISE", "래칫", ("A", "HAIYPRH")),
        ("양호(A+T7 설계체)", "A", "HAIYPRH", "APPROVE?", "—", None),
    ]
    print("\n── M5. 심사(Critic) 효과 (라이브) ──")
    hits = n_bad = 0
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
        approve, _, tier = agent._critic_review(client, model, types, CARGO, choice,
                                                best_seen=best_seen)
        if expect.startswith("REVISE"):
            n_bad += 1
            hits += (not approve)
        print(f"  {label:20s} [{gate:9s}] BBB={choice['bbb']:.2f} 독성={choice['tox']:.2f} "
              f"II={choice['instability']} → 심사 {tier}  (기대 {expect})")
    print(f"  나쁜 후보 REVISE 검출: {hits}/{n_bad}")


def main():
    print("=" * 68)
    print("BBB-Optimizer 에이전트 벤치마크 — M4/M5 (LLM 필요)")
    print("=" * 68)
    m4_agent_efficiency()
    m5_critic_effectiveness()


if __name__ == "__main__":
    main()
