#!/usr/bin/env python3
"""대기열 ⑴ 의 **탐색 경로** A/B — `dev_liab` 응집 중복 제거가 `pareto_scrub` 을 바꾸는가.

    PYTHONPATH=. /usr/bin/python3 scripts/ab_liab_hotspot_search.py

왜 별도 측정인가
---------------------------------------------------------------------------
`scripts/ab_axis_variants.py` 는 **심사** 경로에서 이 변경의 효과가 **정확히 0**
(라이브러리 156/156 에서 상수 −1, 전선 23개 집합 비트 단위 동일)임을 보였다.
그런데 `core/developability_opt.py` 의 `objectives()` 가 **같은** `n_liabilities` 를
탐색 목적으로 쓴다. 그리고 탐색은 **셔틀 단독**을 채점하므로(화물 미전달, §43-6)
핫스팟 플래그가 **상수가 아니다** — 셔틀 12종 중 5종만 단독 핫스팟 ≥ 0.70 이다
(Angiopep-2 0.72 · Tet1 0.73 · Leptin30 0.76 · RVG29 0.76 · Penetratin 0.80).
즉 심사에서 죽은 항이 탐색에서는 **살아서 응집을 두 번 누른다**. 그 크기를 잰다.

두 팔
---------------------------------------------------------------------------
  prod : 프로덕션 그대로
  v1   : `n_liabilities` 에서 응집 핫스팟 플래그 1개를 뺀다(§43-3 중복 계수 제거)

🔴 측정 위생 — §39·§43 과 같은 규율
  ⑴ 패치는 `core.developability` 와 `core.developability_opt` **두 바인딩 모두**에
     꽂는다. `developability_opt` 가 모듈 수준 `from ... import` 를 쓰므로 한쪽만
     바꾸면 **조용히 무효화**된다.
  ⑵ 호출 계수기로 패치 도달을 **증명**한다. 0 이면 "차이 없음" 결론은 무효다.
  ⑶ 결과 객체를 제자리에서 고치지 않고 `replace` 로 새 객체를 만든다(두 번 거치면
     중복 차감된다).
  ⑷ 심사는 **프로덕션 정의로만** 한다 — 심사 쪽 효과가 상수임이 이미 증명됐으므로
     (구성체는 화물 때문에 항상 플래그가 켜진다) 두 팔의 심사 순위는 정의상 같다.

유료 LLM 0회. 시드 0~19 고정.
"""
from __future__ import annotations

import json
import os
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, os.getcwd())
sys.path.append(str(Path(__file__).resolve().parent.parent))

import Bio  # noqa: F401  🔴 없으면 불안정성·off 가 경고 없이 대체값으로 떨어진다

import core.developability as dv
import core.developability_opt as o
from core.config import DEFAULT_CARGO, SHUTTLES
from core.optimizer_agent import dominates, rank_key, risk_axes
from core.stability import assess_stability
from core.selectivity import assess_selectivity
from core.binding import shuttle_similarity
from core.evidence import inherit as evidence_inherit

HOT_FLAG = "응집 핫스팟"
SEEDS = list(range(20))
ROUNDS, POP, FLOOR = 25, 24, 0.85      # 프로덕션 호출값(optimizer_agent_gemini.py:430)
_ORIG = dv.assess_developability
_PROD = _ORIG                           # 심사는 항상 이것으로


class Counter:
    n = 0


def _v1(seq, *a, **k):
    """응집 핫스팟 플래그를 **순위·최소화 축에서만** 뺀다 (사람이 읽는 목록은 유지)."""
    r = _ORIG(seq, *a, **k)
    Counter.n += 1
    n_hot = sum(1 for x in r.liabilities if HOT_FLAG in x)
    return replace(r, n_liabilities=r.n_liabilities - n_hot)


def _set(fn) -> None:
    dv.assess_developability = fn
    o.assess_developability = fn


def judge_row(seq: str, cargo: str, linker: str) -> dict:
    """심사 행 — 키 이름을 `_evaluate` 와 같게 둔다(어댑터 금지). 항상 프로덕션 정의."""
    full, module = cargo + linker + seq, linker + seq
    df = _PROD(full)
    dm = _PROD(module) if module else df
    tier, parent, tsim = evidence_inherit(seq)
    return {"shuttle": seq, "tier": tier, "dev_liab": df.n_liabilities,
            "instability": assess_stability(full).instability_index,
            "sel_off": assess_selectivity(seq).off_target_risk,
            "dev_agg_hot": dm.agg_hotspot}


def run_arm(fn, shuttles: dict) -> dict:
    _set(fn)
    out = {}
    try:
        for name, seq in shuttles.items():
            for sd in SEEDS:
                front = o.pareto_scrub(seq, rounds=ROUNDS, pop=POP,
                                       rmt_floor=FLOOR, seed=sd)
                out[f"{name}|{sd}"] = sorted(m["seq"] for m in front)
    finally:
        _set(_ORIG)
    return out


def main() -> None:
    usable = {}
    for n, v in SHUTTLES.items():
        if shuttle_similarity(v["seq"]).score >= FLOOR:
            usable[n] = v["seq"]
        else:
            print(f"  제외 {n} — rmt_sim {shuttle_similarity(v['seq']).score:.3f} < {FLOOR}")
    print(f"셔틀 {len(usable)}종 × 시드 {len(SEEDS)} = {len(usable)*len(SEEDS)} 셀\n")

    prod = run_arm(_ORIG, usable)
    Counter.n = 0
    v1 = run_arm(_v1, usable)
    assert Counter.n > 0, "🔴 패치가 한 번도 안 불렸다 — 이 측정은 무효다"
    print(f"패치 호출 계수 {Counter.n:,} (0 이면 무효)\n")

    agg = {"n_cells": 0, "same_set": 0, "same_judge": 0, "n_shrink": 0, "n_grow": 0,
           "per_shuttle": {}, "direction": {}}
    detail = []
    for k in prod:
        name, sd = k.split("|")
        a, b = set(prod[k]), set(v1[k])
        same = a == b
        ra = min([judge_row(s, DEFAULT_CARGO, "") for s in a], key=rank_key)
        rb = min([judge_row(s, DEFAULT_CARGO, "") for s in b], key=rank_key)
        ja, jb = ra["shuttle"], rb["shuttle"]
        # 🔴 방향 분해는 **관례 무관** `dominates` 로만 한다 (§43-2). `rank_key` 의 축
        #    순서는 표시 관례이고 주장이 아니므로 "더 좋아졌다"의 근거가 못 된다.
        if ja == jb:
            direction = "동일"
        elif dominates(rb, ra):
            direction = "v1 지배"
        elif dominates(ra, rb):
            direction = "역행(prod 지배)"
        else:
            direction = "교환관계"
        agg["n_cells"] += 1
        agg["same_set"] += same
        agg["same_judge"] += (ja == jb)
        agg["n_shrink"] += len(b) < len(a)
        agg["n_grow"] += len(b) > len(a)
        p = agg["per_shuttle"].setdefault(name, {"n": 0, "same_set": 0, "same_judge": 0,
                                                 "size_prod": 0, "size_v1": 0})
        p["n"] += 1
        p["same_set"] += same
        p["same_judge"] += (ja == jb)
        p["size_prod"] += len(a)
        p["size_v1"] += len(b)
        agg["direction"][direction] = agg["direction"].get(direction, 0) + 1
        p.setdefault("dir", {})[direction] = p.get("dir", {}).get(direction, 0) + 1
        detail.append({"cell": k, "same_set": same, "same_judge": ja == jb,
                       "n_prod": len(a), "n_v1": len(b), "direction": direction,
                       "axes_prod": list(risk_axes(ra)), "axes_v1": list(risk_axes(rb))})

    print(f"{'셔틀':24s}{'출력집합 동일':>14s}{'심사1위 동일':>14s}{'전선크기 prod→v1':>20s}")
    for n, p in agg["per_shuttle"].items():
        print(f"{n:24s}{p['same_set']:>8d}/{p['n']:<5d}{p['same_judge']:>8d}/{p['n']:<5d}"
              f"{p['size_prod']/p['n']:>11.1f} → {p['size_v1']/p['n']:<6.1f}")
    print(f"\n{'합계':24s}{agg['same_set']:>8d}/{agg['n_cells']:<5d}"
          f"{agg['same_judge']:>8d}/{agg['n_cells']:<5d}  축소 {agg['n_shrink']} · 확대 {agg['n_grow']}")
    print(f"\n방향 분해(관례 무관 `dominates`): {agg['direction']}")
    print("  🔴 'v1 지배' 와 '역행' 이 비슷하면 이 변경은 **개선이 아니라 다른 선택**이다(§43-2 와 같은 꼴).")

    p = Path(os.environ.get("AB_OUT", "/tmp/ab_liab_search.json"))
    p.write_text(json.dumps({"agg": agg, "detail": detail,
                             "patch_calls": Counter.n}, ensure_ascii=False, indent=1))
    print(f"\n저장: {p}")


if __name__ == "__main__":
    main()
