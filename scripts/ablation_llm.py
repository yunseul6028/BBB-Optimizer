#!/usr/bin/env python3
"""LLM 애블레이션 — **유전알고리즘이 LLM 의 잔기 편집을 대체하는가?**

## 질문

이 시스템은 "LLM 설계 에이전트"로 소개되지만, **주장을 떠받치는 결정은 거의 전부
결정론적**이다: 채점(7 엔진) · 순위(사전식 = 근거 그룹 → 위험 4축 Pareto) · 하드
게이트(독성) · 래칫(퇴화한 최종 선택을 러닝 베스트로 **교체**). 그리고
`reachable_front` 는 라이브러리 설계 공간을 **전수 열거**한다(156 조합 → 비독성
132 → 후보 89 → 전선 5). 즉 라이브러리 안에서는 LLM 이 할 일이 없다.

남는 **비중복 역할은 `design_candidate` 하나**다 — 라이브러리 **밖**의 서열 편집.
아무도 측정하지 않은 열린 질문:

> **평범한 유전알고리즘(`pareto_scrub`)이 같은 축에서 LLM 의 편집과 맞먹는가?**

맞먹으면 LLM 은 불필요했다. 지면 LLM 의 비중복 기여와 그 **크기**를 말할 수 있다.
어느 쪽이든 결과다 — **어느 쪽으로도 튜닝하지 않는다.**

## 두 팔

- **A 팔 — LLM 0회.** 도달 가능 전선 전수 열거(`refresh=True`) → 후보군 셔틀마다
  `pareto_scrub`(시드 고정, 기본 rounds/pop) → 근거 등급 **상속 필터**(파이프라인과
  동일) → 살아남은 변이체를 **순위와 같은 코드 경로**(`OptimizationAgent._evaluate`
  → `risk_axes`/`dominates`/`rank_key`)로 채점. 축을 다시 계산하지 **않는다**
  (중복 축 로직은 이 프로젝트의 알려진 버그 계열이다).
- **B 팔 — 에이전트.** §22 의 콜드 3세션 기록 + `.cache/ledger` 의 후보별 기록.
  **재실행하지 않는다**(라이브 LLM 호출 0회 — 이 스크립트는 돈을 쓰지 않는다).

## ⚠️ 두 가지 한계 — 먼저 읽을 것

1. 🔴 **통화 단위가 다르다.** `pareto_scrub` 은 셔틀 하나당 목적함수 평가를
   rounds×pop = **~1200회** 쓰고, 에이전트는 세션당 **7~10회의 LLM 매개 평가 호출**을
   썼다. 이 둘은 **같은 단위가 아니다.** 이 비교의 질문은 **"같은 오라클을 쥐었을 때
   나오는 Pareto 집합의 품질"**이고, **"평가를 덜 썼다"가 아니다.** 효율 주장은
   이 데이터로 할 수 없다.
2. ⚠️ **범위 한계.** `pareto_scrub` 은 **셔틀만** 편집한다. 에이전트는 **링커도**
   설계할 수 있다(§23). 따라서 이 비교는 **셔틀 편집 축에서만 유효**하고, **링커
   설계는 LLM 전용 능력으로 남으며 그 축의 대응 널은 미측정**이다.

사용: /usr/bin/python3 scripts/ablation_llm.py [--seed 0] [--cargo GSNKGAIIGLM]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core import ledger                                      # noqa: E402
from core.binding import shuttle_similarity                  # noqa: E402
from core.config import (                                    # noqa: E402
    DEFAULT_CARGO,
    LINKER_LIBRARY,
    SHUTTLES,
    get_settings,
)
from core.developability_opt import (                        # noqa: E402
    OFF_BAND,
    off_target,
    pareto_scrub,
)
from core.evidence import (                                  # noqa: E402
    INHERIT_SIM_MIN,
    SHUTTLE_EVIDENCE,
    group_of,
    group_rank,
    inherit,
)
from core.optimizer_agent import (                           # noqa: E402
    MAX_CANDIDATES_PER_STEP,
    OptimizationAgent,
    dominates,
    in_library,
    is_candidate,
    rank_key,
    risk_axes,
)

# `pareto_scrub` 의 **기본값**을 그대로 쓴다(튜닝 금지 — 질문은 "평범한 GA 가
# 맞먹나"이므로 GA 쪽을 손보면 답이 오염된다). 평가 예산 보고용으로만 읽는다.
SCRUB_ROUNDS, SCRUB_POP = 40, 30


# ── B 팔: §22 콜드 3세션 기록 (docs/experiments.md §22-2 / §22-3) ─────────────
# 🔴 **재실행하지 않는다.** 여기 적힌 수는 그 세션들이 실제로 채점한 값이다.
#    §21-2 세션은 서열이 기록되기 전이라 서열이 없다(그 자체가 §21-6 의 남은 항목이었다).
AGENT_COLD = [
    {"session": "§21-2", "linker": None, "shuttle": None, "tier": "A",
     "dev_liab": 3, "instability": 25.8, "sel_off": 0.11, "dev_agg_hot": 0.51,
     "n_eval": 7, "note": "서열 미기록"},
    {"session": "콜드 1", "linker": "", "shuttle": "TYRYGGSRGKRNNFKTEQY", "tier": "A",
     "dev_liab": 3, "instability": 20.1, "sel_off": 0.20, "dev_agg_hot": 0.52,
     "n_eval": 9, "note": ""},
    {"session": "콜드 2", "linker": "EAAAK", "shuttle": "TYRYGGSRGKRNNFKTEQY", "tier": "A",
     "dev_liab": 3, "instability": 19.5, "sel_off": 0.20, "dev_agg_hot": 0.52,
     "n_eval": 10, "note": ""},
]
# 부모 전선 멤버 — §22-3 이 설계체와 직접 대조한 바로 그 멤버.
AGENT_PARENT_FRONT = {"label": "A(EAAAK)3A+Angiopep-2", "tier": "A", "dev_liab": 3,
                      "instability": 38.9, "sel_off": 0.09, "dev_agg_hot": 0.75}


# ── 예측기 재시도 래퍼 ──────────────────────────────────────────────────────
# ToxinPred3 는 `cwd=vendor/toxinpred3` 에서 중간파일(`seq.aac`)을 **저장소 디렉터리에
# 쓴다**. 그래서 같은 저장소에서 다른 프로세스가 동시에 돌면 서로 파일을 지워
# `FileNotFoundError: seq.aac not found` 로 **조용히가 아니라 시끄럽게** 깨진다.
# core/toxicity.py 를 고치는 것이 옳은 해법이지만 이 스크립트의 소관이 아니므로,
# **여기서만** 재시도로 감싼다(결과값은 바꾸지 않는다 — 실패한 호출만 다시 한다).
class _Retry:
    def __init__(self, inner, tries: int = 8, wait: float = 4.0):
        self._inner, self._tries, self._wait = inner, tries, wait

    def __getattr__(self, k):
        return getattr(self._inner, k)

    def predict_many(self, seqs):
        last = None
        for i in range(self._tries):
            try:
                return self._inner.predict_many(seqs)
            except Exception as e:                            # noqa: BLE001
                last = e
                msg = " ".join(str(e).split())[:70]
                print(f"    ⚠️ 예측기 충돌 재시도 {i + 1}/{self._tries}: {msg}…")
                time.sleep(self._wait)
        raise RuntimeError(f"예측기 재시도 소진: {last}")


def _score(agent: OptimizationAgent, cargo: str, combos: list[dict]) -> list[dict]:
    """구성체들을 **에이전트와 동일한 경로**(`_evaluate`)로 채점한다.

    🔴 축을 여기서 다시 계산하지 않는다. `_evaluate` 가 돌려주는 행을 그대로
    `risk_axes`/`dominates`/`rank_key` 에 넣는다 — 이 프로젝트는 축 로직을 복제해
    두 번 데였다(§19-4 상수 복제 · MIN_BBB_SIGNAL 의 ledger.py 사본).
    """
    rows: list[dict] = []
    n = len(combos)
    for i in range(0, n, MAX_CANDIDATES_PER_STEP):
        _, chunk = agent._evaluate(cargo, {"constructs": combos[i:i + MAX_CANDIDATES_PER_STEP]})
        rows.extend(chunk)
        print(f"    채점 {min(i + MAX_CANDIDATES_PER_STEP, n)}/{n}")
    return rows


def _axrow(r: dict) -> str:
    ax = risk_axes(r)
    return (f"{r.get('dev_liab')} | {float(r.get('instability') or 0):.1f}"
            f"{'(불안정)' if ax[1] else ''} | {float(r.get('sel_off') or 0):.2f} | "
            f"{float(r.get('dev_agg_hot') or 0):.2f}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="LLM 애블레이션 — GA(pareto_scrub) vs LLM 잔기 편집",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=0, help="GA 난수 시드(기본 0)")
    ap.add_argument("--cargo", default=DEFAULT_CARGO, help=f"화물 서열(기본 {DEFAULT_CARGO})")
    args = ap.parse_args()
    cargo, seed = args.cargo, args.seed

    print("=" * 78)
    print("LLM 애블레이션 — 평범한 GA 가 LLM 의 잔기 편집과 맞먹는가")
    print("=" * 78)
    print(f"화물 `{cargo}` · **GA 시드 {seed}** · pareto_scrub(rounds={SCRUB_ROUNDS}, "
          f"pop={SCRUB_POP}) = **함수 기본값** · 라이브 LLM 호출 **0회**")
    print("\n🔴 **통화 단위가 다르다**: GA 는 셔틀 하나당 목적함수 평가 "
          f"{SCRUB_ROUNDS}×{SCRUB_POP}=**{SCRUB_ROUNDS * SCRUB_POP}회**, 에이전트는 "
          "세션당 **LLM 매개 평가 7~10회**. 같은 단위가 아니다 — 이 비교는 "
          "**같은 오라클을 쥐었을 때 나오는 Pareto 집합의 품질**이고, "
          "**'평가를 덜 썼다'가 결코 아니다.**")
    print("⚠️ **범위 한계**: `pareto_scrub` 은 **셔틀만** 편집한다. 에이전트는 "
          "**링커도** 설계한다 → 이 비교는 **셔틀 편집 축에서만** 유효하고, "
          "**링커 설계는 LLM 전용 능력이며 그 대응 널은 미측정**이다.\n")

    settings = get_settings()
    agent = OptimizationAgent(settings)
    agent.predictor = _Retry(agent.predictor)
    agent.tox_predictor = _Retry(agent.tox_predictor)

    # ── 1. 도달 가능 전선 전수 열거 ────────────────────────────────────────
    # 🔴 `refresh=True` 필수 — 이 전선은 `.cache/reachable/` 에 캐시되고, **낡은 캐시가
    #    이미 한 번 틀린 결론을 만들었다**(축 정의가 바뀌면 캐시 키를 올려야 하는데
    #    올리지 않으면 옛 BBB 기반 전선이 섞인다 — reachable_front 주석).
    print("[1] 도달 가능 전선 전수 열거 (refresh=True — 낡은 캐시 금지)")
    t0 = time.time()
    res = agent.reachable_front(cargo, refresh=True)
    print(f"    전수 {res['n_total']} → 비독성 {res['n_nontoxic']} → 후보군 {res['n_gate']} "
          f"(RMT {res['n_gate_rmt']}) → **전선 {len(res['front'])}** "
          f"(최상위 그룹 전선 {len(res['front_top'])}) · {time.time() - t0:.0f}s")
    print(f"    최상위 근거 그룹 = **{res['best_group']}** · 그룹별 {res['by_group']}")
    print("\n    전선(비지배 집합) — liab | 불안정성 | off | 응집  "
          "(★ = `front_top` = 최상위 근거 그룹 전선):")
    front, front_top = res["front"], res["front_top"]
    top_labels = {f["label"] for f in front_top}
    for f in front:
        print(f"      {'★' if f['label'] in top_labels else ' '} {f['label']:<28} "
              f"[{f['tier']}·{f['group']}] {_axrow(f)}")
    print("    ⚠️ 두 집합을 섞지 말 것 — `front`(전체 후보군 비지배)와 "
          "`front_top`(최상위 근거 그룹만)은 다른 기준점이다. 아래 판정은 **둘 다** 적는다.")

    # ── 2. 후보군 셔틀마다 pareto_scrub ───────────────────────────────────
    # 후보군 셔틀 = **비교군(CPP) 아님 + 근거 그룹 '근거 있음/주의'**(= `is_candidate` 의
    # 셔틀 수준 조건). 독성은 구성체 수준 조건이라 셔틀을 지우지 못한다 —
    # 아래 [4] 에서 부모 구성체를 함께 채점해 각 셔틀이 실제로 후보군에 있는지 확인한다.
    cand_shuttles = [n for n, v in SHUTTLES.items()
                     if group_rank(SHUTTLE_EVIDENCE[n].tier) <= 1
                     and "CPP" not in shuttle_similarity(v["seq"]).mechanism.upper()]
    print(f"\n[2] 후보군 셔틀 {len(cand_shuttles)}종에 pareto_scrub (시드 {seed}, "
          f"rounds={SCRUB_ROUNDS}·pop={SCRUB_POP} = 기본값)")
    print("    제약: rmt_sim ≥ 0.85 · 문헌 보존 위치 동결 · off 양측 밴드(시드 ±0.05)")
    variants: list[dict] = []          # {name, parent, seq, objectives, off_target, n_edits}
    for n in cand_shuttles:
        pseq = SHUTTLES[n]["seq"]
        vs = pareto_scrub(pseq, seed=seed)
        edited = [v for v in vs if v["n_edits"] > 0]
        for j, v in enumerate(edited, 1):
            # `vid` — 같은 부모·같은 편집수 변이체가 여러 개라서 라벨이 충돌한다.
            # 충돌하면 출력에서 서로 다른 서열이 같은 이름으로 보인다(실제로 봤다).
            variants.append({"parent": n, "parent_seq": pseq, "vid": f"{n}#v{j}", **v})
        print(f"      · {n:<16} GA 전선 {len(vs):>2}개 (편집된 변이체 {len(edited):>2}개) "
              f"· 목적평가 ~{SCRUB_ROUNDS * SCRUB_POP}회")
    print(f"    편집 변이체 합계 **{len(variants)}개** · "
          f"목적평가 총 ~{len(cand_shuttles) * SCRUB_ROUNDS * SCRUB_POP}회")

    # ── 3. 근거 등급 상속 필터 (파이프라인과 동일) ────────────────────────
    print(f"\n[3] 근거 등급 **상속** 필터 — 설계체는 유사도 ≥ {INHERIT_SIM_MIN} 일 때만 "
          "부모 등급을 물려받는다 (못 받으면 '미등재' = 최종 후보 불가)")
    kept, dropped = [], []
    for v in variants:
        tier, parent, sim = inherit(v["seq"])
        v["tier"], v["tier_parent"], v["tier_sim"] = tier, parent, sim
        (kept if (tier is not None and group_rank(tier) <= 1) else dropped).append(v)
    print(f"    상속 통과 **{len(kept)}/{len(variants)}** · 탈락 {len(dropped)}개"
          + (f" (최고 유사도 {max(d['tier_sim'] for d in dropped):.3f})" if dropped else ""))
    if not kept:
        print("    🔴 통과한 변이체가 없다 — A 팔은 **양의 신호를 가진 설계체를 하나도 "
              "만들지 못했다**. 여기서 중단한다.")
        return

    # ── 4. 순위와 같은 코드 경로로 채점 ───────────────────────────────────
    # 축 4개 중 3개(liability·불안정성·응집)는 **링커에 의존**한다(off 만 셔틀 주도).
    # 그래서 에이전트가 쥐었던 것과 같은 운동장을 준다 — 라이브러리 링커 13종 전부.
    # 부모 셔틀(편집 0)도 같이 채점해 **편집 전/후**를 같은 자에서 읽는다.
    linkers = [v["seq"] or "" for v in LINKER_LIBRARY.values()]
    combos = []
    for v in kept:
        for lk in linkers:
            combos.append({"label": f"{v['vid']}({v['n_edits']}ed)/{lk or '직접융합'}",
                           "linker": lk, "shuttle": v["seq"]})
    for n in cand_shuttles:
        for lk in linkers:
            combos.append({"label": f"PARENT:{n}/{lk or '직접융합'}",
                           "linker": lk, "shuttle": SHUTTLES[n]["seq"]})
    print(f"\n[4] 변이체 {len(kept)}종 × 라이브러리 링커 {len(linkers)}종 + 부모 대조 "
          f"= 구성체 **{len(combos)}개**를 `_evaluate`(에이전트와 동일 경로)로 채점")
    rows = _score(agent, cargo, combos)

    designed = [r for r in rows if not str(r["label"]).startswith("PARENT:")]
    parents = [r for r in rows if str(r["label"]).startswith("PARENT:")]
    # 🔴 후보 자격 판정도 **그 함수**를 쓴다 (비독성 + 비교군 아님 + 근거 있음).
    ok = [r for r in designed if is_candidate(r)]
    par_ok = [r for r in parents if is_candidate(r)]
    par_shuttles = {r["shuttle"] for r in par_ok}
    print(f"    후보 자격 통과: 설계체 **{len(ok)}/{len(designed)}** · "
          f"부모 대조 {len(par_ok)}/{len(parents)}")
    print(f"    부모 셔틀 {len(par_shuttles)}/{len(cand_shuttles)}종이 실제로 후보군에 "
          "있다(독성 게이트는 구성체 수준이라 셔틀을 통째로 지우지 않는지 확인)")
    if not ok:
        print("    🔴 A 팔에 후보 자격을 가진 설계체가 없다. 여기서 중단한다.")
        return

    # ── 5. 전선 대조 — 비지배 / 지배 ──────────────────────────────────────
    print("\n[5] 전선 대조 (`dominates` = 근거 그룹 + 위험 4축, 같은 함수)")
    nondom = [r for r in ok if not any(dominates(f, r) for f in front)]
    nondom_top = [r for r in ok if not any(dominates(f, r) for f in front_top)]
    domin = [(r, [f for f in front if dominates(r, f)]) for r in ok]
    domin = [(r, fs) for r, fs in domin if fs]
    domin_top = [(r, [f for f in fs if f["label"] in top_labels]) for r, fs in domin]
    domin_top = [(r, fs) for r, fs in domin_top if fs]
    print(f"    **비지배**: front(22급) **{len(nondom)}/{len(ok)}** · "
          f"front_top **{len(nondom_top)}/{len(ok)}**")
    print(f"    전선 멤버를 **지배**: front **{len(domin)}/{len(ok)}** · "
          f"front_top **{len(domin_top)}/{len(ok)}**")
    print("    ← 에이전트는 이 값이 **0** 이었다(§22-3, 콜드 3세션 전부)")

    # 🔴 지배의 **출처**를 분해한다. off-target 축은 §26 의 대응 널 심사에서
    #    **악화 검출은 신뢰 / 개선 신호는 불신** 판정을 받았다. 그러므로 "off 가
    #    내려갔다"에만 의존하는 지배는 **주장으로 쓸 수 없다.** 축을 뺀 3축
    #    (liability·불안정성·응집)에서도 지배가 유지되는 쌍만 **견고한 지배**다.
    def _no_off(r: dict) -> dict:
        return {**r, "sel_off": 0.0}

    robust = [(r, [f for f in fs if dominates(_no_off(r), _no_off(f))])
              for r, fs in domin]
    robust = [(r, fs) for r, fs in robust if fs]
    robust_top = [(r, [f for f in fs if f["label"] in top_labels]) for r, fs in robust]
    robust_top = [(r, fs) for r, fs in robust_top if fs]
    print(f"\n    그중 **off 개선에 의존하지 않는 지배**(liability·불안정성·응집만으로도 "
          f"성립): front **{len(robust)}/{len(domin)}** · front_top "
          f"**{len(robust_top)}/{len(domin_top) if domin_top else 0}**")
    print(f"    = off 개선에 **의존하는** 지배 {len(domin) - len(robust)}개는 "
          "**주장으로 쓸 수 없다**(§26: off 는 악화 검출만 신뢰).")

    if domin:
        print("\n    견고한 지배 사례 (3축만으로도 성립 — 설계체 → 지배당한 전선 멤버):")
        seen_f: dict[str, int] = {}
        for r, fs in sorted(robust, key=lambda x: rank_key(x[0]))[:8]:
            print(f"      · {r['label']:<34} [{r['tier']}] {_axrow(r)}")
            for f in fs:
                seen_f[f["label"]] = seen_f.get(f["label"], 0) + 1
                print(f"          ▸ 지배{'★' if f['label'] in top_labels else ''}: "
                      f"{f['label']:<28} [{f['tier']}] {_axrow(f)}")
        if len(robust) > 8:
            print(f"      … 외 {len(robust) - 8}개")
        print(f"\n    (견고한 지배로) 지배당한 전선 멤버 빈도: {seen_f}")
        allf: dict[str, int] = {}
        for _r, fs in domin:
            for f in fs:
                allf[f["label"]] = allf.get(f["label"], 0) + 1
        print(f"    (off 의존 포함 전체) 지배당한 전선 멤버 빈도: {allf}")
    else:
        print("    (A 팔도 전선 멤버를 하나도 지배하지 못했다 — 에이전트와 같다.)")

    # A 팔의 '최종 선택' — LLM 이 없으므로 시스템의 **결정론 선택자**(`rank_key`)를 쓴다.
    # ⚠️ `rank_key` 의 축 순서는 **표시 관례이고 주장이 아니다**(optimizer_agent 주석).
    #    그래서 아래 세 질문의 답은 '최종 1개'와 '비지배 집합 전체'를 **둘 다** 적는다.
    best_g = min(group_rank(r.get("tier")) for r in ok)
    final_a = min([r for r in ok if group_rank(r.get("tier")) == best_g], key=rank_key)

    # ── 6. B 팔 — 기록된 콜드 세션 (재실행 없음) ──────────────────────────
    print("\n[6] B 팔 — 에이전트 콜드 3세션 **기록** (재실행 없음, LLM 0회)")
    led = ledger.load(cargo)
    wins = led.get("wins") or []
    led_designed = [w for w in wins if not in_library(w.get("linker") or "", w.get("shuttle") or "")]
    print(f"    원장 `.cache/ledger` : 세션 {led.get('sessions')} · 기록 {len(wins)}개 · "
          f"그중 **라이브러리 밖 설계체 {len(led_designed)}개** → "
          + ("**사용 가능**" if led_designed else "**없음 → 문서로 대체**"))
    for w in led_designed:
        print(f"      · {(w.get('linker') or '직접융합'):<12}+{w['shuttle']:<30} "
              f"[{w.get('tier')}←{w.get('tier_parent')}·{w.get('tier_sim')}] {_axrow(w)}")
    print("\n    §22-2 콜드 3세션 최종 선택:")
    for a in AGENT_COLD:
        sh = a["shuttle"] or "(미기록)"
        lk = "직접융합" if a["linker"] == "" else (a["linker"] or "?")
        print(f"      · {a['session']:<7} {lk:<10}+{sh:<30} [{a['tier']}] {_axrow(a)} "
              f"· 평가 {a['n_eval']}회 {a['note']}")
    # 원장 기록으로 문서 표를 **교차 검증**한다(두 출처가 어긋나면 그 자체가 결과다).
    xchk = []
    for a in AGENT_COLD:
        if not a["shuttle"]:
            continue
        m = [w for w in led_designed
             if w["shuttle"] == a["shuttle"] and (w.get("linker") or "") == (a["linker"] or "")]
        if m:
            w = m[0]
            same = (w["dev_liab"] == a["dev_liab"]
                    and abs(w["instability"] - a["instability"]) < 0.05
                    and abs(w["sel_off"] - a["sel_off"]) < 0.005
                    and abs(w["dev_agg_hot"] - a["dev_agg_hot"]) < 0.005)
            xchk.append((a["session"], "일치" if same else "**불일치**"))
        else:
            xchk.append((a["session"], "원장에 없음"))
    print(f"    문서 ↔ 원장 교차검증: {dict(xchk)}")
    print("    ⚠️ 등급 글자 비교 주의 — 기록 당시 Angiopep-2 는 **A**, §24(GRADE 간접성)"
          f" 이후 **{SHUTTLE_EVIDENCE['Angiopep-2'].tier}** 다. 다만 순위에 쓰는 "
          f"**그룹은 둘 다 '{group_of(SHUTTLE_EVIDENCE['Angiopep-2'].tier)}'(0)** 이므로 "
          "질문 ①의 답은 바뀌지 않는다.")

    # B 팔의 세 질문 — 기록값을 **같은 함수**에 넣어 판정한다.
    b_rows = [a for a in AGENT_COLD]
    b_nondom = sum(1 for a in b_rows if not any(dominates(f, a) for f in front))
    b_dom = sum(1 for a in b_rows if any(dominates(a, f) for f in front))
    b_top = sum(1 for a in b_rows if group_rank(a["tier"]) == 0)
    b_outside = sum(1 for a in b_rows
                    if a["shuttle"] is None or not in_library(a["linker"] or "", a["shuttle"]))

    # ── 7. 맞대결 ─────────────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("[7] 맞대결 — M4 의 세 질문 + 위험 4축")
    print("=" * 78)
    a_top = group_rank(final_a.get("tier")) == 0
    a_nondom_final = not any(dominates(f, final_a) for f in front)
    a_outside = not in_library(final_a.get("linker") or "", final_a.get("shuttle") or "")
    print("\n| 질문 | A 팔 (GA · LLM 0회) | B 팔 (에이전트 3세션) |")
    print("|---|---|---|")
    print(f"| ① 최상위 근거 그룹인가 | {'예' if a_top else '아니오'} "
          f"({final_a.get('tier')}·{group_of(final_a.get('tier'))}) | {b_top}/3 |")
    print(f"| ② 전선에 비지배인가 | {'예' if a_nondom_final else '아니오'} "
          f"(비지배 집합 front {len(nondom)} · front_top {len(nondom_top)}) | {b_nondom}/3 |")
    print(f"| ③ 라이브러리 밖인가 | {'예' if a_outside else '아니오'} | {b_outside}/3 |")
    print(f"| ④ 전선 멤버를 **지배**하나 | **{len(domin)}개**(front) · "
          f"**{len(domin_top)}개**(front_top) · 그중 off 개선에 **의존하지 않는** "
          f"견고한 지배 **{len(robust)}개** | **{b_dom}개** (0 — §22-3) |")

    print("\n위험 4축 — liab | 불안정성 | off | 응집  (전부 ↓ 좋음)")
    print("|  | 근거 | liab | 불안정성 | off | 응집 |")
    print("|---|---|---|---|---|---|")
    print(f"| 부모 전선 멤버 {AGENT_PARENT_FRONT['label']} | A/B | "
          f"{_axrow(AGENT_PARENT_FRONT)} |")
    for a in AGENT_COLD:
        sh = (a["shuttle"] or "미기록")[:18]
        print(f"| B 팔 {a['session']} ({sh}) | {a['tier']} | {_axrow(a)} |")
    print(f"| **A 팔 최종**(rank_key) {final_a['label']} | {final_a.get('tier')} | "
          f"{_axrow(final_a)} |")
    print("\nA 팔 비지배 집합 상위 10 (rank_key 순):")
    for r in sorted(nondom, key=rank_key)[:10]:
        v = next((x for x in kept if x["seq"] == r["shuttle"]), None)
        ed = f"{v['n_edits']}편집·rmt {v['rmt_sim']}" if v else "?"
        print(f"  · {r['label']:<34} [{r['tier']}←{r['tier_parent']}·{r['tier_sim']}] "
              f"{_axrow(r)}  ({ed})")

    # 축별 최선 — "GA 가 어디까지 갔나"를 축마다 하나씩.
    print("\nA 팔 축별 최선 (후보 자격 통과 설계체 중):")
    axes_spec = (("dev_liab", "liability", "{:.0f}"),
                 ("instability", "불안정성", "{:.1f}"),
                 ("sel_off", "off-target", "{:.2f}"),
                 ("dev_agg_hot", "응집 핫스팟", "{:.2f}"))
    for key, lab, fmt in axes_spec:
        b = min(ok, key=lambda r: float(r.get(key) or 0))
        print(f"  · {lab:<11} {fmt.format(float(b.get(key) or 0)):>6}  ← {b['label']}")

    # ── 7-b. **같은 부모**에서의 직접 대조 ────────────────────────────────
    # 에이전트의 콜드 3세션은 **전부 Angiopep-2 로 수렴**했다(§22-4). 다른 부모까지
    # 섞으면 "부모가 달라서 생긴 차이"를 "편집기가 달라서 생긴 차이"로 오독한다.
    print("\n" + "-" * 78)
    print("[7-b] **같은 부모(Angiopep-2)**에서의 직접 대조 — 부모 효과를 제거한 비교")
    print("-" * 78)
    ang = SHUTTLES["Angiopep-2"]["seq"]
    ang_ok = [r for r in ok if r["shuttle"] != ang
              and str(r["label"]).startswith("Angiopep-2#")]
    off_seed = off_target(ang)
    print(f"    시드 셔틀 Angiopep-2 `{ang}` · 시드 off **{off_seed:.2f}** · "
          f"현행 밴드 제약 |off − 시드| ≤ **{OFF_BAND}**")
    print(f"    A 팔 Angiopep-2 설계체(후보 자격 통과) {len(ang_ok)}개")
    for key, lab, fmt in axes_spec:
        if not ang_ok:
            break
        b = min(ang_ok, key=lambda r: float(r.get(key) or 0))
        print(f"      · {lab:<11} 최선 {fmt.format(float(b.get(key) or 0)):>6}  "
              f"← {b['label']} ({_axrow(b)})")
    # 🔴 운동장이 같지 않다 — off 밴드 제약은 **에이전트의 실패를 보고 나중에 추가**됐다.
    print("\n    🔴 **운동장이 같지 않다**: off 양측 밴드 제약(§26, OFF_BAND=0.05)은 "
          "바로 이 콜드 세션들의 off 악화를 보고 **나중에** 도입됐다. "
          "B 팔의 설계체를 현행 제약에 넣어 보면:")
    for a in AGENT_COLD:
        if not a["shuttle"]:
            print(f"      · {a['session']:<7} (서열 미기록 — 판정 불가, off {a['sel_off']:.2f})")
            continue
        d = abs(a["sel_off"] - off_seed)
        print(f"      · {a['session']:<7} off {a['sel_off']:.2f} · |Δ| {d:.2f} → "
              + ("**밴드 통과**" if d <= OFF_BAND + 1e-9
                 else f"🔴 **밴드 위반**(> {OFF_BAND}) = 현행 GA 는 이 서열을 "
                      "**애초에 만들 수 없다**"))
    print("    ⇒ A 팔은 **더 좁은 제약** 아래 탐색했고 B 팔은 그 제약 없이 탐색했다. "
          "A 팔에 유리한 비교가 아니라 **서로 다른 문제를 푼 것**이다 — "
          "이 항목을 빼고 지배 개수만 인용하면 과장이 된다.")

    # ── 7-c. A 팔 ↔ B 팔 **직접** 지배 — 전선을 매개하지 않는 유일한 맞대결 ──
    # 전선 대조는 둘 다 '전선 대비'이므로 간접 비교다. 두 팔을 **서로** 넣어 보는 것이
    # 이 질문("GA 가 LLM 의 편집과 맞먹나")에 대한 **가장 직접적인** 답이다.
    print("\n" + "-" * 78)
    print("[7-c] A 팔 ↔ B 팔 **직접** 대조 (전선을 매개하지 않는 맞대결)")
    print("-" * 78)
    print("    ⚠️ B 팔 기록의 등급 글자는 A(현행 B)지만 `dominates` 는 **그룹**을 쓰고 "
          "둘 다 0 이라 비교가 성립한다.")
    b_named = [a for a in AGENT_COLD if a["shuttle"]]
    for a in b_named:
        beats = [r for r in ok if dominates(r, a)]
        beats_rb = [r for r in beats if dominates(_no_off(r), _no_off(a))]
        loses = [r for r in ok if dominates(a, r)]
        print(f"\n    B 팔 {a['session']} ({'직접융합' if a['linker'] == '' else a['linker']}"
              f"+{a['shuttle']}) {_axrow(a)}")
        print(f"      · 이 설계체를 **지배하는 A 팔 설계체**: **{len(beats)}개** "
              f"(그중 off 개선에 의존하지 않는 견고한 지배 **{len(beats_rb)}개**)")
        print(f"      · 반대로 이 설계체가 **지배하는 A 팔 설계체**: **{len(loses)}개**")
        for r in sorted(beats_rb, key=rank_key)[:3]:
            print(f"          ▸ 예: {r['label']:<30} {_axrow(r)}")
    # 요약 — 어느 쪽이 상대를 지배하는가를 한 줄로.
    tot_beats = sum(1 for a in b_named for r in ok if dominates(r, a))
    tot_loses = sum(1 for a in b_named for r in ok if dominates(a, r))
    print(f"\n    합계: A→B 지배 **{tot_beats}쌍** vs B→A 지배 **{tot_loses}쌍**")

    # ── 8. 한계 재고지 ────────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("🔴 이 표를 읽는 규칙")
    print("=" * 78)
    print(f"1. **통화 단위가 다르다.** GA 는 셔틀당 목적함수 평가 "
          f"~{SCRUB_ROUNDS * SCRUB_POP}회(총 "
          f"~{len(cand_shuttles) * SCRUB_ROUNDS * SCRUB_POP}회), 에이전트는 세션당 "
          "LLM 매개 평가 7~10회. 이 비교는 **같은 오라클을 쥐었을 때의 Pareto 집합 "
          "품질**이고 **'평가를 덜 썼다'가 아니다.** 효율 주장은 불가.")
    print("2. **셔틀 편집 축에서만 유효하다.** `pareto_scrub` 은 셔틀만 편집한다 — "
          "**링커 설계는 LLM 전용 능력**이고 그 축의 대응 널은 **미측정**이다.")
    print("3. **비지배 ≠ 우위.** 전선과 비지배라는 것은 교환관계라는 뜻이다(§22-3). "
          "양쪽 팔 모두 '넘어섰다'로 쓰면 안 된다 — 지배 개수(질문 ④)만이 우위 주장이다.")
    print("4. B 팔은 **기록**이고 A 팔은 **이번 실행**이다. 그 사이 §24(GRADE 간접성)로 "
          "등급 글자가 바뀌었다(A→B). 순위 그룹은 불변이라 ①은 안전하지만, "
          "**글자끼리 비교하지 말 것.**")
    print("5. **운동장이 다르다**(7-b): off 양측 밴드 제약은 B 팔의 off 악화를 보고 "
          "**나중에** 추가됐다. B 팔 설계체 일부는 현행 제약에서 **생성 불가**다. "
          "두 팔은 엄밀히 말해 **같은 탐색 문제를 풀지 않았다.**")
    print("6. **지배 개수를 그대로 인용하지 말 것** — off 개선에 의존하는 지배는 "
          "§26 에서 신뢰할 수 없다고 판정된 신호다. 인용할 수 있는 수는 "
          f"**견고한 지배 {len(robust)}개**(liability·불안정성·응집만으로 성립)다.")
    print(f"7. 재현: /usr/bin/python3 scripts/ablation_llm.py --seed {seed} --cargo {cargo}")


if __name__ == "__main__":
    main()
