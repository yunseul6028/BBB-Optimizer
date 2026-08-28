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
from benchmark import _library_scores


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
    oracle = max(valid, key=lambda r: r["bbb"])["bbb"] if valid else 0.0

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
    print(f"  라이브러리 oracle 최적 BBB : {oracle:.3f}")
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

    # 의도적으로 문제 있는 후보 + 좋은 후보 (라벨: 기대 판정)
    probes = [
        ("독성 셔틀(멜리틴)", "GGGGS", "GIGAVLKVLTTGLPALISWIKRKRQQ", "REVISE"),
        ("저BBB(poly-Ser)", "GGGGS", "SSSSSSSSSSSS", "REVISE"),
        ("양호(SynB1+EAAAK2)", "EAAAKEAAAK", "RGGRLSYSRRRFSTSTGR", "APPROVE?"),
    ]
    print("\n── M5. 심사(Critic) 효과 (라이브) ──")
    hits = 0
    for label, lk, sh, expect in probes:
        choice = agent._score_choice("GSNKGAIIGLM",
                                     {"chosen_label": label, "chosen_linker": lk, "chosen_shuttle": sh})
        if choice is None:
            print(f"  {label}: 채점 실패")
            continue
        approve, _ = agent._critic_review(client, model, types, "GSNKGAIIGLM", choice)
        verdict = "APPROVE" if approve else "REVISE"
        ok = (expect.startswith("REVISE") and not approve)
        hits += ok if expect.startswith("REVISE") else 0
        print(f"  {label:22s} BBB={choice['bbb']:.2f} 독성={choice['tox']:.2f} "
              f"II={choice['instability']} → 심사 {verdict}  (기대 {expect})")
    print(f"  나쁜 후보 REVISE 검출: {hits}/2")


def main():
    print("=" * 68)
    print("BBB-Optimizer 에이전트 벤치마크 — M4/M5 (LLM 필요)")
    print("=" * 68)
    m4_agent_efficiency()
    m5_critic_effectiveness()


if __name__ == "__main__":
    main()
