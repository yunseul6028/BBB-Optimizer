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

import time

from core.config import TOXICITY_THRESHOLD, get_settings
from core.predictors import get_predictor
from core.toxicity import get_toxicity_predictor
from core.optimizer_agent_gemini import GeminiOptimizationAgent, get_gemini_agent
from benchmark import _library_scores, is_cpp_row


# max_rounds=6 은 의도적 선택 — "타이트한 예산에서도 에이전트가 효율적인가"를 보는 실측이다.
# (운영 기본값은 10 = 안전 상한, get_gemini_agent 참고. 여기 6은 수렴 바닥을 재는 프로브다.)
def m4_agent_efficiency(cargo="GSNKGAIIGLM", max_rounds=6, use_ledger=False,
                        seed_ledger=False, tag=""):
    """에이전트 효율/최적성 갭.

    `use_ledger`: **콜드 스타트(False)가 M4 의 정의다.** True 는 누적 모드(원장 ON)의
    별도 실험용이며, 두 모드는 **비교 불가**이므로 반드시 나눠 보고한다(core/ledger.py).
    `seed_ledger`: 콜드 세션의 실측을 원장에 **기록만** 한다(읽지는 않는다). 누적 모드
    실험의 바닥을 실제 세션 결과로 깔기 위한 것이고, 콜드 측정 자체에는 영향이 없다
    (원장은 `use_ledger=True` 일 때만 읽힌다).

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
    agent = get_gemini_agent(settings, max_rounds=max_rounds, use_ledger=use_ledger)
    if agent is None:
        print("⚠️ Gemini 에이전트 비활성(GEMINI_API_KEY 필요) — M4 생략.")
        return None
    from core.evidence import group_of, group_rank
    from core.optimizer_agent import dominates, is_candidate, rank_key

    front = agent.reachable_front(cargo)
    ft = front.get("front_top") or []
    fr = front.get("front") or []
    mode = "누적 모드(원장 ON)" if use_ledger else "콜드 스타트"
    print(f"── M4. 에이전트 효율 (라이브 실행 · {mode}{tag}) ──")
    if use_ledger:
        from core import ledger as _ledger
        _d = _ledger.load(cargo)
        _fl = _ledger.floor_entry(cargo)
        print(f"  원장: {_d.get('sessions')}세션 기록 · 승자 {len(_d.get('wins') or [])}개 · "
              f"실패 {len(_d.get('fails') or [])}개 · 바닥 "
              + (f"`{(_fl.get('linker') or '직접융합')}+{_fl.get('shuttle')}`"
                 f"(근거 {_fl.get('tier')}·liab {_fl.get('dev_liab')}·off {_fl.get('sel_off')})"
                 if _fl else "없음(빈 원장 → 사실상 콜드)"))
    print(f"  도달 가능 집합: 후보군 {front.get('n_gate')}개 · 근거 그룹 "
          f"{front.get('by_group')} · 최상위 그룹 `{front.get('best_group')}`")
    print(f"  최상위 그룹 전선 {len(ft)}개 / 전체 전선 {len(fr)}개")

    n_eval, choice, best_seen = 0, None, None
    rows_seen = []
    for ev in agent.run(cargo):
        if ev.kind == "evaluation":
            rows = ev.data.get("rows", [])
            n_eval += len(rows)
            rows_seen.extend(rows)
            for r in rows:
                if is_candidate(r) and (best_seen is None or rank_key(r) < rank_key(best_seen)):
                    best_seen = r
        elif ev.kind in ("choice", "optimum"):
            choice = ev.data
        elif ev.kind == "error":
            print(f"  에이전트 오류: {ev.text}")
            return None
    pick = choice if choice is not None else best_seen
    if pick is None:
        print("  ⚠️ 최종 후보 없음 — 모든 후보가 게이트에 걸렸다.")
        return None

    # 콜드 세션의 실측을 원장에 **기록만** 한다(이 세션은 원장을 읽지 않았다).
    # 누적 모드 실험의 바닥이 실제 세션에서 나오게 하는 용도다.
    if seed_ledger and not use_ledger:
        from core.ledger import failure_notes, record
        d = record(cargo, rows_seen, pick, failure_notes(rows_seen))
        print(f"  📒 원장 시딩: 누적 {d.get('sessions')}세션 · 승자 {len(d.get('wins') or [])}개 · "
              f"실패 {len(d.get('fails') or [])}개")

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
    # 역할 분리 게이트의 거부 횟수. `scripts/ablation_tools.py` 는 전부터 이 값을 기록하는데
    # 메인 벤치마크는 계측만 하고 **보고하지 않았다** — 계측된 값이 아무 데도 안 나오면
    # 그 자체가 §37 계열의 "조용함"이다. 세는 것은 **독성 하드 게이트가 아니라**
    # `evaluate_candidates` 가 편집 서열을 되돌려 보낸 횟수다(= 역할 분리가 실제로 작동한 횟수).
    gate_rej = getattr(agent, "gate_rejections", 0)
    print(f"  역할 분리 게이트 거부: {gate_rej}회  "
          f"(편집 서열을 `evaluate_candidates` 에 넣으려 한 횟수 — 독성 게이트와 무관)")
    if not on_top:
        verdict = "근거 하향 — 더 단단한 근거의 대안이 있었다"
    elif dom_by:
        verdict = "지배당함 — 모든 축에서 더 나은 달성 가능 후보가 있었다"
    elif not in_front:
        verdict = "전선 초과(설계로 라이브러리 밖에서 비지배 달성)"
    else:
        verdict = "전선 도달(라이브러리 내 비지배)"
    print(f"  판정: {verdict}")
    return {"verdict": verdict, "q1_on_top": on_top, "q2_nondominated": not dom_by,
            "q3_outside_front": not in_front, "n_eval": n_eval,
            "linker": pick.get("linker") or "", "shuttle": pick.get("shuttle"),
            "tier": pick.get("tier"), "group": g, "dev_liab": pick.get("dev_liab"),
            "instability": pick.get("instability"), "sel_off": pick.get("sel_off"),
            "dev_agg_hot": pick.get("dev_agg_hot"), "bbb_report": pick.get("bbb"),
            "dominated_by": dom_by, "use_ledger": use_ledger,
            "gate_rejections": gate_rej}


# ── M5 프로브 표 · 게이트 마커 (모듈 수준) ────────────────────────────────
# 모듈 수준에 두는 이유: `--dry-run-m5` 검증기와 라이브 M5 가 **같은 표 하나**를
# 공유해야 한다. 표를 복제하면 검증한 표와 돌리는 표가 갈라져, 이 주석이 기록한
# 두 사건(프로브 3개 교정 · 간접성 하향 후 2개 재기초)이 또 조용히 재발한다.
# 나쁜 후보 5종(각기 **다른 게이트**를 때린다) + 좋은 후보 1종.
# ⚠️ 2026-11-10 개편: 게이트 체인이 바뀌었다. 독성 → **근거 등급 보유** → 등급 D →
#    CPP 비교군 → 근거 하향 → 지배 → 래칫. 옛 프로브가 겨냥했던 'BBB 하한'·'순서 특이성'·
#    '전선 이탈'은 전부 폐지됐고, 그 역할은 **근거 등급 보유 게이트**가 승계했다.
# 🔴 2026-11-11: **간접성 하향**(GRADE indirectness)으로 등급이 바뀌면서 프로브 2개를
#    재기초했다 — Angiopep-2 A→B · T7/Leptin30/SynB1 → C, front_top 11→5개.
#    교훈: **프로브는 등급 체계와 전선에 의존한다.** 목적함수·등급 규칙을 건드리면
#    유료 실행 전에 프로브를 로컬로 전수 재검증해야 한다(지난번 3개 교정에 이어 두 번째).
M5_CARGO = "GSNKGAIIGLM"
M5_PROBES = [
    # (라벨, 링커, 셔틀, 기대, 겨냥한 게이트, best_seen 로 줄 더 좋은 후보)
    ("독성 셔틀(멜리틴)", "GGGGS", "GIGAVLKVLTTGLPALISWIKRKRQQ", "REVISE", "독성", None),
    # 옛 "저BBB"·"순서 측정불가" 프로브의 후계. poly-Ser 는 어떤 라이브러리 셔틀과도
    # 유사도가 상속 하한(0.70)에 못 미친다 → **근거 미등재**로 걸린다.
    ("근거 미등재(poly-Ser)", "GGGGS", "SSSSSSSSSSSS", "REVISE", "근거 미등재", None),
    # 옛 "조성쌓기 설계체" 프로브의 후계 — 양이온·방향족 과편집(보상해킹의 전형).
    # ⚠️ 실측(로컬): 부모 Angiopep-2 와 유사도 **0.745 ≥ 0.70** 이라 등급 B 를
    #    **상속한다**. 상속 게이트는 모티프를 보존한 과편집을 막지 않는다 — 의도된 설계다
    #    (merit 의 출처는 부모 모티프이고, '얼마나 편집했나'가 아니라 **위험 축**이 심사한다).
    #    이 후보는 liab 5·불안정·off 0.29·응집 0.81 로 **front_top 멤버 5개 전부에게
    #    지배당한다**(로컬 확인). 그래서 겨냥하는 게이트는 '근거 미등재'가 아니라 **'지배'** 다.
    # 🔴 2026-11-11 교체: 옛 프로브(ApoE 변이 `LRVYQAGRR…`)는 **더 이상 지배당하지 않는다.**
    #    간접성 하향으로 T7 이 front_top 에서 빠지자(응집 0.59 멤버가 사라짐) 그 후보의
    #    응집 0.59 가 오히려 전선보다 좋아져 지배 관계가 끊겼다. 프로브는 전선에 의존한다.
    ("과편집 설계체(양이온 과편집)", "", "TFFYWGSRGKRNNKKTERW", "REVISE", "지배", None),
    # 등급 D 게이트 — Tet1 은 라이브러리 셔틀이라 상속은 되지만 물려받은 근거가 '없음'이다.
    # 우리 효능 축에서는 후보군 1위였던 조합(M2·M4 oracle)이라 **가장 중요한 프로브**다.
    ("근거 없음(등급 D·Tet1)", "", "HLNILSTLWKYR", "REVISE", "근거 등급 D", None),
    ("비교군(SynB1·CPP)", "EAAAKEAAAK", "RGGRLSYSRRRFSTSTGR", "REVISE", "CPP 비교군", None),
    # 근거 하향 게이트 — B6 은 등급 **C(주의)** 인데 최상위 그룹 '근거 있음'(A·B)이
    # 달성 가능하다. 등급 D 도 CPP 도 아니므로 앞선 게이트를 전부 지나 여기서 걸린다.
    ("근거 하향(등급 C·B6)", "AEAAAKEAAAKEAAAKA", "GHKAKGPRK", "REVISE", "근거 하향", None),
    # 래칫 프로브 — 앞선 게이트를 **전부** 통과해야 래칫을 때린다. 그래서 제출 후보를
    # **front_top 멤버**(`A(EAAAK)3A`+Angiopep-2, 축 3·0·0.09·0.75)로 두고, best_seen 만
    # 그것을 지배하는 **설계체**(개발성 Pareto 스크럽 변이 `TAFYGGSRGKGNNFKTEEY` —
    # Angiopep-2 등급 B 상속(유사도 0.858), liab 3→2·off 0.09→0.05·응집 0.75→0.65,
    # 불안정성 동일 = **3축에서 지배**)로 준다. 로컬 실측으로 지배 성립을 확인했다.
    # ⚠️ 바로 아래 '양호' 프로브는 **같은 후보에 best_seen 만 없는 것**이다 —
    #    둘의 판정 차이가 곧 래칫 게이트의 순효과다(다른 변수가 없다).
    # 🔴 2026-11-11 교체: 옛 프로브는 `GGS+T7` 이었는데 간접성 하향으로 **T7 이 등급 C**
    #    가 되어 '양호' 프로브가 **근거 하향 게이트에 먼저 걸린다**(= 오탈락). 래칫 프로브도
    #    래칫 이전에 근거 하향이 발동해 게이트 식별이 어긋난다. 그래서 최상위 그룹
    #    (Angiopep-2·ApoE) 기반으로 재기초했다. **프로브는 등급 체계에 의존한다.**
    ("퇴행(더 좋은 걸 봤음)", "AEAAAKEAAAKEAAAKA", "TFFYGGSRGKRNNFKTEEY", "REVISE",
     "래칫", ("AEAAAKEAAAKEAAAKA", "TAFYGGSRGKGNNFKTEEY")),
    ("양호(A(EAAAK)3A+Angiopep-2, 전선 멤버)", "AEAAAKEAAAKEAAAKA",
     "TFFYGGSRGKRNNFKTEEY", "APPROVE?", "—", None),
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


def _fired_gate(txt: str) -> str:
    hit = [nm for nm, mk in GATE_MARKERS if mk in (txt or "")]
    return "+".join(hit) if hit else "—"


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

    # 프로브 표·게이트 마커는 모듈 수준(`M5_PROBES`/`GATE_MARKERS`)이 유일한 출처다.
    # `--dry-run-m5` 가 **이 표 그대로** 검증한다.
    CARGO, probes, fired_gate = M5_CARGO, M5_PROBES, _fired_gate
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


# ── `--dry-run-m5` : 프로브 표 로컬 검증 (LLM 0회 · 네트워크 0회) ──────────
# 왜 있나: 위 주석이 기록한 두 사건 — 프로브 3개 교정, 간접성 하향 후 2개 재기초 —
# 둘 다 **누군가 손으로 로컬 점검을 기억했기 때문에** 잡혔다. 한 번은 래칫 프로브가
# 래칫이 아니라 근거 하향에 먼저 걸리고 있었는데 '통과'로 집계됐다(엉뚱한 이유의 통과).
# 이 모드는 그 손 점검을 코드로 고정해 **유료 실행 전에** 프로브 전수를 검증한다.
#
# 🔑 설계 원칙 — **게이트 조건을 다시 쓰지 않는다.**
#   어느 게이트가 먼저 걸리는지는 `GeminiOptimizationAgent._critic_review` 의 if/elif
#   체인이 **유일한 출처**다. 그래서 여기서는 조건을 재구현하지 않고, LLM 호출만 더미로
#   막은 뒤 그 함수를 **그대로 호출**한다(아래 `_DryTypes`/`_DryClient`).
#   체인이 elif 라 발동 게이트는 최대 1개이고, 더미 응답의 본문은 빈 문자열이므로
#   반환 텍스트에 남는 건 **실제로 발동한 게이트 문구뿐**이다.
#   → 복제한 게이트 조건은 **없다**. (복제했다면 반드시 여기 적어야 한다 —
#     폐지한 BBB 하한이 ledger.py 의 사본에서 살아남은 전례가 있다.)


class _DryParts:
    """`types.Part` 자리. 프롬프트는 어차피 버려지므로 내용을 만들지 않는다."""

    @staticmethod
    def from_text(text=""):
        return None


class _DryTypes:
    """`google.genai.types` 자리를 메우는 더미 — 인자를 받아 버리기만 한다."""

    GenerateContentConfig = staticmethod(lambda **kw: None)
    ThinkingConfig = staticmethod(lambda **kw: None)
    Content = staticmethod(lambda **kw: None)
    Part = _DryParts


class _DryModels:
    def generate_content(self, **kw):
        # 네트워크·유료 호출 없음. candidates 를 비워 심사문 **본문을 빈 문자열**로 만든다
        # → 결정론적 게이트 문구만 남아 "어느 게이트가 걸렸나"가 깨끗이 분리된다.
        class _R:
            candidates = []
        return _R()


class _DryClient:
    models = _DryModels()


def _pad(text: str, width: int) -> str:
    """표 정렬용 — 한글은 폭 2칸이므로 `f"{s:28s}"` 로는 열이 어긋난다."""
    from unicodedata import east_asian_width
    w = sum(2 if east_asian_width(c) in "WF" else 1 for c in text)
    return text + " " * max(0, width - w)


def _dry_score(agent, cargo, label, linker, shuttle):
    """드라이런용 안전 채점 — 실패하면 (None, 사유) 를 돌려준다.

    ⚠️ 왜 try 가 필요한가(실측 2026-11-12): 로컬 ToxinPred3 는 vendor 저장소를 **cwd 로
    공유**하며 중간 파일(seq.aac·seq.dpc)을 그 안에 쓴다. 같은 저장소에서 다른 스크립트가
    동시에 채점하면 그 파일이 섞여 추론이 터진다(배열 행 수 불일치). 검증기가 여기서
    죽으면 "프로브가 깨졌다"와 "환경이 겹쳤다"를 구분할 수 없으므로, 한 번 재시도하고
    그래도 실패하면 **사유를 적어** 깨진 프로브로 집계한다.
    """
    last = ""
    tries = 3            # 충돌은 짧게 기다리면 풀린다 → 2·4초 백오프로 세 번 시도.
    for attempt in range(tries):
        try:
            row = agent._score_choice(cargo, {"chosen_label": label, "chosen_linker": linker,
                                              "chosen_shuttle": shuttle})
            if row is not None:
                return row, ""
            last = "빈 결과"
        except Exception as exc:  # noqa: BLE001 — 사유를 보고하고 계속 간다.
            last = str(exc).split("\n")[0]
        if attempt < tries - 1:
            time.sleep(2 * (attempt + 1))
    return None, last


def dry_run_m5(cargo: str = None) -> int:
    """M5 프로브 표를 **로컬로** 전수 검증한다. 반환: 종료 코드(0=정상, 1=깨진 프로브 있음).

    프로브는 **등급 체계와 Pareto 전선에 의존**해서 조용히 썩는다. 유료 M5 를 돌리기 전에
    이걸 통과시켜라(종료 코드로 게이팅할 수 있다).
    """
    from core.evidence import group_of
    from core.optimizer_agent import dominates

    CARGO = cargo or M5_CARGO
    settings = get_settings()
    # 팩토리 `get_gemini_agent` 는 API 키를 요구한다 → 드라이런은 클래스를 직접 만든다.
    # 쓰는 건 로컬 예측기(deepB3P·ToxinPred3)와 결정론적 게이트뿐이라 **키가 없어도 돈다**.
    agent = GeminiOptimizationAgent(settings, max_rounds=6, use_ledger=False,
                                    ledger_write=False)

    print("=" * 68)
    print("M5 프로브 로컬 검증 (--dry-run-m5) — LLM 호출 0회")
    print("=" * 68)
    print(f"  화물: {CARGO}  (프로브 표는 이 화물로 교정되어 있다)")

    # 🔴 캐시 무효화는 선택이 아니다. 전선은 `.cache/reachable/` 에 캐시되고, 묵은 캐시로
    #    검증하면 **거짓 통과**가 난다 — 이 함정이 이 프로젝트에서 실제로 틀린 결론을 냈다.
    #    여기서 refresh 로 다시 깔아두면 `_critic_review` 내부의 `reachable_front` 호출도
    #    같은 신선한 캐시를 읽는다(두 경로가 같은 파일 하나를 본다).
    front = agent.reachable_front(CARGO, refresh=True)
    ft = front.get("front_top") or []
    print(f"  전선 재계산(refresh=True): 후보군 {front.get('n_gate')}개 · 최상위 그룹 "
          f"`{front.get('best_group')}` · front_top {len(ft)}개")
    for r in ft:
        print(f"    · {r['label']}: 근거 {r['tier']} · liab {r['dev_liab']} · 불안정성 "
              f"{r['instability']} · off {r['sel_off']} · 응집 {r['dev_agg_hot']}")

    rows = []
    for label, lk, sh, expect, gate, better in M5_PROBES:
        choice, why = _dry_score(agent, CARGO, label, lk, sh)
        if choice is None:
            print(f"\n  🔴 {label}: 채점 실패 — {why}")
            rows.append({"label": label, "gate": gate, "got": "채점실패",
                         "ok": False, "verdict": "🔴 채점 실패", "lk": lk, "sh": sh,
                         "hits": [], "better": better})
            continue
        best_seen = None
        if better is not None:
            best_seen, why = _dry_score(agent, CARGO, "본최고", better[0], better[1])
            if best_seen is None:
                print(f"\n  🔴 {label}: best_seen 채점 실패 — {why} "
                      f"(래칫 게이트를 검증할 수 없다)")
                rows.append({"label": label, "gate": gate, "got": "채점실패",
                             "ok": False, "verdict": "🔴 best_seen 채점 실패",
                             "lk": lk, "sh": sh, "hits": [], "better": better})
                continue

        # ① 결정론적 사실 — 게이트 판정의 입력이 그대로 보이게 찍는다.
        tier = choice.get("tier")
        sim = choice.get("tier_sim")
        if tier is None:
            ev = (f"미등재(최근접 {choice.get('tier_parent') or '없음'}·유사도 "
                  f"{(sim or 0.0):.3f})")
        elif (sim or 1.0) >= 1.0:
            ev = f"{tier}({group_of(tier)}) · 라이브러리 셔틀 그대로"
        else:
            ev = (f"{tier}({group_of(tier)}) · 상속 ← 부모 {choice.get('tier_parent')}"
                  f"·유사도 {sim:.3f}")
        dom_by = [r["label"] for r in ft if dominates(r, choice)]
        print(f"\n  {label}  [겨냥 게이트: {gate}]")
        print(f"    서열: 링커 `{lk or '직접융합'}` + 셔틀 `{sh}`")
        print(f"    독성: {choice['tox']:.3f} "
              f"({'탈락' if choice.get('toxic') else '통과'}, 임계 "
              f"{TOXICITY_THRESHOLD:.2f}) · 메커니즘 {choice.get('mechanism')}")
        print(f"    근거: {ev}")
        print(f"    순위 4축: liab {choice.get('dev_liab')} · 불안정성 "
              f"{choice.get('instability')} · off {choice.get('sel_off')} · 응집 "
              f"{choice.get('dev_agg_hot')}")
        print(f"    전선 지배: " + (f"지배당함 ← {', '.join(dom_by)}" if dom_by else "비지배"))
        if best_seen is not None:
            print(f"    best_seen: `{better[0] or '직접융합'}+{better[1]}` "
                  f"({agent._summary(best_seen)}) · 제출후보 지배? "
                  f"{'예' if dominates(best_seen, choice) else '아니오'}")

        # ② 어느 게이트가 먼저 걸리나 — **실제 체인**을 그대로 돌려서 알아낸다.
        _approve, critique, _tier = agent._critic_review(
            _DryClient(), settings.gemini_model, _DryTypes, CARGO, choice,
            best_seen=best_seen)
        hits = [nm for nm, mk in GATE_MARKERS if mk in (critique or "")]
        got = "+".join(hits) if hits else "—"

        # ③ 판정. (드라이런은 LLM 본문이 없으므로 REVISE/APPROVE 자체는 의미가 없다 —
        #    **게이트 발동 여부만** 판정한다. 결정론적 게이트가 LLM 판정을 덮어쓰므로
        #    게이트가 걸리면 라이브에서도 반드시 REVISE 다.)
        if gate == "—":
            ok = not hits
            verdict = ("✅ 게이트 없음(설계대로 — 승인은 LLM 판단에 맡겨진다)" if ok else
                       f"🔴 걸려선 안 되는 게이트 발동 — `{got}` (오탈락)")
        elif not hits:
            ok = False
            verdict = (f"⚠️ 아무 게이트도 안 걸린다 — 겨냥 `{gate}`. 나쁜 프로브인데 "
                       f"승인될 수 있다(거짓 통과)")
        elif hits == [gate]:
            ok = True
            verdict = f"✅ 겨냥한 게이트 `{gate}` 발동"
        elif gate in hits:
            ok = False
            verdict = (f"🔴 마커 중복 — 겨냥 `{gate}` 는 걸렸지만 `{got}` 로 복수 식별된다"
                       f"(게이트 마커를 고유 문구로 고쳐야 한다)")
        else:
            ok = False
            verdict = f"🔴 다른 게이트 발동 — 겨냥 `{gate}` / 실제 `{got}`"
        print(f"    → {verdict}")
        rows.append({"label": label, "gate": gate, "got": got, "ok": ok,
                     "verdict": verdict, "lk": lk, "sh": sh, "hits": hits,
                     "better": better})

    # ── 구조적 불변식: '양호' 와 '래칫' 은 **같은 구성체**이고 best_seen 만 다르다 ──
    # 그래서 두 판정의 차이는 **정확히 래칫 게이트 하나**여야 한다. 다른 게이트가 끼면
    # 래칫의 순효과를 분리해 측정하지 못한다(2026-11-11 사건이 바로 이것이었다).
    inv_ok = True
    rt = next((r for r in rows if r["gate"] == "래칫"), None)
    gd = next((r for r in rows if r["gate"] == "—"), None)
    print("\n── 구조적 불변식: '양호' vs '래칫' (같은 구성체 · best_seen 만 다름) ──")
    if rt is None or gd is None:
        inv_ok = False
        print("  🔴 래칫 프로브 또는 양호 프로브를 표에서 찾지 못했다.")
    else:
        same = (rt["lk"], rt["sh"]) == (gd["lk"], gd["sh"])
        print(f"  구성체 동일성   : {'✅ 같다' if same else '🔴 다르다'} "
              f"(래칫 `{rt['lk'] or '직접융합'}+{rt['sh']}` vs 양호 "
              f"`{gd['lk'] or '직접융합'}+{gd['sh']}`)")
        print(f"  best_seen 차이  : 래칫 {'있음' if rt.get('better') else '없음'} / "
              f"양호 {'있음' if gd.get('better') else '없음'}")
        print(f"  게이트 차이     : 래칫 `{rt['got']}` − 양호 `{gd['got']}`")
        diff_only_ratchet = (rt["hits"] == ["래칫"] and gd["hits"] == [])
        print("  → 차이가 래칫 하나인가: "
              + ("✅ 예 (래칫 게이트의 순효과가 분리된다)" if diff_only_ratchet else
                 "🔴 아니오 — 다른 게이트가 끼어 래칫의 순효과를 분리할 수 없다"))
        inv_ok = same and diff_only_ratchet

    print("\n── 프로브 검증 요약 ──")
    print(f"  {_pad('프로브', 38)} {_pad('겨냥', 14)} {_pad('실제', 14)} 판정")
    for r in rows:
        print(f"  {_pad(r['label'], 38)} {_pad(r['gate'], 14)} {_pad(r['got'], 14)} "
              + (r["verdict"][:2] if r["ok"] else r["verdict"]))
    bad = [r for r in rows if not r["ok"]]
    print(f"  정상 {len(rows) - len(bad)}/{len(rows)} 프로브"
          + (f" · 깨진 프로브 {len(bad)}개: " + "; ".join(r["label"] for r in bad)
             if bad else ""))
    code = 0 if (not bad and inv_ok) else 1
    if code:
        print("  🔴 유료 M5 실행 금지 — 프로브를 고치기 전에는 '검출 n/n' 이 "
              "**엉뚱한 이유의 통과**일 수 있다.")
    else:
        print("  ✅ 프로브 표가 겨냥한 게이트를 그대로 때린다 — 유료 M5 실행 가능.")
    return code


def _summarize(results, title):
    """여러 세션의 M4 결과를 **분산까지** 요약한다.

    1세션치로 효율을 주장하면 안 된다 — 옛 기준 4세션 재현에서 유효점수가 0.408~0.635 로
    흩어졌다(docs/experiments.md §11.2). 판정은 백분율이 아니라 예/아니오 세 개이므로,
    요약도 **세션별 판정의 분포**와 평가 횟수 범위로 낸다(평균 하나로 뭉개지 않는다).
    """
    ok = [r for r in results if r]
    if not ok:
        print(f"\n── {title}: 유효 세션 없음 ──")
        return
    print(f"\n── {title} — {len(ok)}세션 요약 ──")
    print(f"  {'#':>2}  {'①근거최상위':<12} {'②비지배':<9} {'③전선밖':<9} {'평가':>4}  "
          f"판정 / 최종 선택")
    for i, r in enumerate(ok, 1):
        print(f"  {i:>2}  {'예' if r['q1_on_top'] else '아니오':<12} "
              f"{'예' if r['q2_nondominated'] else '아니오':<9} "
              f"{'예' if r['q3_outside_front'] else '아니오':<9} {r['n_eval']:>4}  "
              f"{r['verdict']}  [{r['linker'] or '직접융합'}+{r['shuttle']} · 근거 {r['tier']} · "
              f"liab {r['dev_liab']} · II {r['instability']} · off {r['sel_off']}]")
    n = len(ok)
    evals = [r["n_eval"] for r in ok]
    print(f"  ① 최상위 근거 그룹 : {sum(r['q1_on_top'] for r in ok)}/{n}")
    print(f"  ② 비지배           : {sum(r['q2_nondominated'] for r in ok)}/{n}")
    print(f"  ③ 전선 밖(설계체)  : {sum(r['q3_outside_front'] for r in ok)}/{n}")
    import statistics
    print(f"  평가 횟수          : {min(evals)}~{max(evals)} (중앙값 "
          f"{statistics.median(evals):g})")
    uniq = {(r["linker"], r["shuttle"]) for r in ok}
    print(f"  서로 다른 최종 선택: {len(uniq)}종 / {n}세션"
          + ("  ← 세션 간 변동 있음" if len(uniq) > 1 else "  ← 전 세션 동일"))


def main(argv=None):
    import argparse

    # 🔴 **측정 전 전제조건 검사** (§32·§37). 참조 명부가 라이브러리 셔틀을 덮지 않으면
    #    `rmt_similarity` 가 조용히 0.0 이 되고, off 가 과대평가되며 `rmt_floor` 제약이
    #    시드 자신을 거부해 개발성 최적화가 **아무것도 하지 않는다** — 그런데 예외도
    #    경고도 없어서, 측정은 멀쩡히 끝나고 숫자만 틀린다. 유료 실행 앞에서 막는다.
    #    (이 검사를 "정의만 하고 부르지 않으면" `gate_rejections` 와 같은 꼴이 된다.)
    from core.binding import assert_roster_covers_library
    assert_roster_covers_library()

    p = argparse.ArgumentParser(description="M4/M5 에이전트 벤치마크 (LLM 필요)")
    p.add_argument("--sessions", type=int, default=1,
                   help="M4 콜드 스타트 세션 수 (분산 추정에 2~3 권장)")
    p.add_argument("--cumulative", type=int, default=0, metavar="N",
                   help="누적 모드(원장 ON) 세션 수. 콜드와 **따로** 보고한다")
    p.add_argument("--seed-ledger", action="store_true",
                   help="콜드 세션 결과를 원장에 기록(읽지는 않음) — 누적 실험의 바닥용")
    p.add_argument("--cargo", default="GSNKGAIIGLM")
    p.add_argument("--skip-m4", action="store_true")
    p.add_argument("--skip-m5", action="store_true")
    p.add_argument("--dry-run-m5", action="store_true",
                   help="M5 프로브 표를 로컬 검증만 한다(LLM 0회). 깨진 프로브가 있으면 "
                        "종료 코드 1 — 유료 실행을 게이팅하는 용도")
    args = p.parse_args(argv)

    # 드라이런은 **유료 경로를 전혀 건드리지 않는다** — M4 전에 바로 반환한다.
    if args.dry_run_m5:
        return dry_run_m5()

    print("=" * 68)
    print("BBB-Optimizer 에이전트 벤치마크 — M4/M5 (LLM 필요)")
    print("=" * 68)

    cold = []
    if not args.skip_m4:
        for i in range(max(1, args.sessions)):
            cold.append(m4_agent_efficiency(cargo=args.cargo, seed_ledger=args.seed_ledger,
                                            tag=f" · 세션 {i + 1}/{args.sessions}"))
            print()
        if args.sessions > 1:
            _summarize(cold, "M4 콜드 스타트")

    warm = []
    for i in range(args.cumulative):
        warm.append(m4_agent_efficiency(cargo=args.cargo, use_ledger=True,
                                        tag=f" · 세션 {i + 1}/{args.cumulative}"))
        print()
    if warm:
        _summarize(warm, "누적 모드(원장 ON) — 콜드와 비교 불가, 별도 보고")

    if not args.skip_m5:
        m5_critic_effectiveness()
    return 0


if __name__ == "__main__":
    # 종료 코드를 흘려보낸다 — `--dry-run-m5` 가 유료 실행의 게이트로 쓰이기 때문이다.
    raise SystemExit(main())
