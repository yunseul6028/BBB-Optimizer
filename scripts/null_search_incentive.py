#!/usr/bin/env python3
"""§39 — 탐색 유인 대응 널: §28-C 가 "양전하를 지워라"는 끌림을 없앴는가.

§38-4 는 **순위**를 재측정했고 전선이 불변이었다. 이 스크립트는 §28 이 실제로
겨냥했던 곳, 즉 **탐색**(`pareto_scrub`)을 잰다.

질문은 "GA 가 K/R 을 지우는가"가 아니다 — 길이 보존 치환이므로 무작위 표류만으로도
지워진다. 질문은 **"같은 편집 수의 무작위보다 더 지우는가"**이고, 그 초과분이 유인이다.

세 팔, 제약은 동일(문헌 보존 위치 동결 · `rmt_sim ≥ rmt_floor` · off 양측 밴드):
  OLD  : liability 축 = 화학적 취약성 + 전하 주도 항목 (§28-C 이전)
  NEW  : liability 축 = 화학적 취약성만 (현재 코드)
  NULL : 같은 편집 수 · 같은 제약의 **무작위 치환** (대응 널)

판정은 **시드 단위 쌍대 비교 + 양측 부호검정**으로 한다. 아카이브 멤버는 같은 시드 안에서
서로 상관되므로 멤버를 독립 표본으로 세면 p 가 과대평가된다(멤버 단위로는 n≈285 지만
독립 단위는 시드 60개다).

실행: /usr/bin/python3 scripts/null_search_incentive.py
"""
from __future__ import annotations

import random
import statistics as st
import sys
from dataclasses import replace
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import core.developability as dv
import core.developability_opt as o
from core.config import SHUTTLES
from core.evidence import CONSERVED_POSITIONS

_ORIG = dv.assess_developability
POS, NEG = "KR", "DE"
RMT_FLOOR = 0.85
N_SEEDS = 60
TARGETS = ("Angiopep-2", "B6", "RVG29")


def _old_def(seq, *a, **k):
    """§28-C 이전 정의 복원 — 전하 항목을 순위·최소화 축에 다시 합산.

    🔴 제자리 변형이 아니라 `replace` 로 **새 객체**를 만든다. 결과 객체를 제자리에서
    고치면 같은 객체를 두 번 거칠 때 전하 항목이 두 번 더해지고, 두 팔이 서로를
    오염시켜 "두 팔이 동일하다"는 **가짜 결론**이 나온다.
    """
    r = _ORIG(seq, *a, **k)
    return replace(r, n_liabilities=r.n_liabilities + r.n_charge_liabilities)


def _set_arm(fn) -> None:
    """탐색 경로가 닿는 **모든** 바인딩을 바꾼다.

    `developability_opt` 는 모듈 수준에서 `from .developability import
    assess_developability` 로 가져오므로(opt.py:22), `core.developability` 쪽만
    바꾸면 **조용히 무효화**된다. 두 곳을 모두 바꾼다.
    """
    dv.assess_developability = fn
    o.assess_developability = fn


def n_pos(s: str) -> int:
    return sum(c in POS for c in s)


def net_charge(s: str) -> int:
    return n_pos(s) - sum(c in NEG for c in s)


def matched_null(seed_seq: str, k: int, rng: random.Random, frozen,
                 off0: float, tries: int = 3000):
    """시드에서 **정확히 k 자리** 다르고 **같은 제약을 통과하는** 무작위 변이 1개."""
    idx = [i for i in range(len(seed_seq)) if i not in frozen]
    if not idx:
        return None
    for _ in range(tries):
        s = list(seed_seq)
        for i in rng.sample(idx, min(k, len(idx))):
            s[i] = rng.choice(o._STD)
        cand = "".join(s)
        if sum(a != b for a, b in zip(cand, seed_seq)) != k:
            continue                                    # 시드로 되돌아온 충돌
        _, sim = o._fitness(cand, RMT_FLOOR)
        if sim < RMT_FLOOR:
            continue
        if abs(o.off_target(cand) - off0) > o.OFF_BAND + 1e-9:
            continue
        return cand
    return None


def sign_test(diffs) -> tuple[int, int, float]:
    """양측 부호검정. 반환 (GA>널 시드 수, 동점 제외 시드 수, p)."""
    wins = sum(1 for d in diffs if d > 0)
    m = sum(1 for d in diffs if d != 0)
    if m == 0:
        return wins, 0, 1.0
    tail = sum(comb(m, i) for i in range(wins, m + 1)) / 2 ** m
    return wins, m, min(1.0, 2 * tail)


def main() -> None:
    print(f"대응 널: 같은 편집 수 · 동결 위치 · rmt_floor {RMT_FLOOR} · off 밴드 ±{o.OFF_BAND}")
    print(f"시드 {N_SEEDS}개 × 셔틀 {len(TARGETS)}종 × 팔 2개(+각각의 널)\n")

    agg: dict = {}
    for name in TARGETS:
        seq0 = SHUTTLES[name]["seq"]
        frozen = o._frozen_idx(seq0)
        off0 = o.off_target(seq0)
        pins = sorted(CONSERVED_POSITIONS.get(name) or ())
        sites = [i + 1 for i, c in enumerate(seq0) if c in POS and i not in frozen]
        print(f"── {name} ({len(seq0)}aa) · 시드 K/R {n_pos(seq0)} · 순전하 "
              f"{net_charge(seq0):+d} · 동결 {pins or '없음'} · 편집 가능한 양전하 자리 {sites}")

        for label, fn in (("OLD", _old_def), ("NEW", _ORIG)):
            _set_arm(fn)
            rng = random.Random(777)
            per_seed = {p: [] for p in sites}          # 시드별 (GA율 − 널율)
            d_kr, nd_kr, d_net, nd_net, sizes, edits = [], [], [], [], [], []
            for sd in range(N_SEEDS):
                try:
                    arc = o.pareto_scrub(seq0, seed=sd, rmt_floor=RMT_FLOOR)
                except ValueError as exc:              # 시드 밀반입 차단(§33-4)
                    print(f"   [{label}] seed {sd}: {exc}")
                    continue
                if not arc:
                    continue
                nulls = [matched_null(seq0, m["n_edits"], rng, frozen, off0) for m in arc]
                nulls = [x for x in nulls if x]
                if not nulls:
                    continue
                sizes.append(len(arc))
                edits += [m["n_edits"] for m in arc]
                d_kr += [n_pos(m["seq"]) - n_pos(seq0) for m in arc]
                nd_kr += [n_pos(x) - n_pos(seq0) for x in nulls]
                d_net += [net_charge(m["seq"]) - net_charge(seq0) for m in arc]
                nd_net += [net_charge(x) - net_charge(seq0) for x in nulls]
                for p in sites:
                    g = 100 * sum(m["seq"][p - 1] not in POS for m in arc) / len(arc)
                    n = 100 * sum(x[p - 1] not in POS for x in nulls) / len(nulls)
                    per_seed[p].append((g, n))
            agg[(name, label)] = dict(
                per_seed=per_seed, d_kr=d_kr, nd_kr=nd_kr, d_net=d_net,
                nd_net=nd_net, sizes=sizes, edits=edits)
        print()

    _set_arm(_ORIG)                                     # 원상복귀

    print("=" * 100)
    print("① 총량 — 양전하(K/R) 개수 변화. 초과 = GA − 널. **양수면 GA 가 무작위보다 전하를 지킨다**")
    print(f"{'셔틀':12} {'팔':4} {'후보':>5} {'편집':>5} | {'ΔK/R':>8} {'널':>8} {'초과':>8} | "
          f"{'Δ순전하':>8} {'널':>8}")
    for name in TARGETS:
        for label in ("OLD", "NEW"):
            d = agg[(name, label)]
            if not d["d_kr"]:
                continue
            m, nm = st.mean(d["d_kr"]), st.mean(d["nd_kr"])
            print(f"{name:12} {label:4} {len(d['d_kr']):>5} {st.mean(d['edits']):>5.1f} | "
                  f"{m:>8.3f} {nm:>8.3f} {m - nm:>+8.3f} | "
                  f"{st.mean(d['d_net']):>8.3f} {st.mean(d['nd_net']):>8.3f}")

    print()
    print("② 자리별 — 시드 단위 쌍대 비교 + 양측 부호검정")
    print(f"{'셔틀':12} {'자리':5} {'팔':4} {'GA%':>6} {'널%':>6} {'초과pp':>8} {'GA>널':>9} {'p':>8}")
    for name in TARGETS:
        seq0 = SHUTTLES[name]["seq"]
        for p in sorted(agg[(name, "NEW")]["per_seed"]):
            for label in ("OLD", "NEW"):
                rec = agg[(name, label)]["per_seed"][p]
                if not rec:
                    continue
                g = st.mean(x for x, _ in rec)
                n = st.mean(y for _, y in rec)
                wins, m, pv = sign_test([x - y for x, y in rec])
                flag = " 🔴" if pv < 0.05 and g > n else ""
                print(f"{name:12} {seq0[p-1]}{p:<4} {label:4} {g:>6.1f} {n:>6.1f} "
                      f"{g - n:>+8.1f} {wins:>5}/{m:<3} {pv:>8.4f}{flag}")

    print()
    print("③ 탐색력 — 아카이브(비지배 집합) 크기. 줄어들면 축 하나가 기울기를 잃은 것이다")
    for name in TARGETS:
        a, b = agg[(name, "OLD")]["sizes"], agg[(name, "NEW")]["sizes"]
        if a and b:
            print(f"   {name:12} OLD {st.mean(a):5.2f}개 → NEW {st.mean(b):5.2f}개 "
                  f"({100*(st.mean(b)/st.mean(a)-1):+.1f}%)")


if __name__ == "__main__":
    main()
