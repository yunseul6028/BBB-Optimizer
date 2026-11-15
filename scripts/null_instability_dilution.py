#!/usr/bin/env python3
"""§43-6 차단 사유 **#2** — "구성체 불안정성은 길이로 희석돼 탐색이 셔틀 쪽 기울기를
잃는다"는 가설을 잰다. 대기열 ⑵ 의 선결 조건이고, §52-1 이 사유 #1 을 해소한 뒤
**유일하게 남은 미측정 차단 사유**다.

가설 (§43-6)
-----------
`pareto_scrub` 의 `objectives()` 는 **셔틀 단독** 불안정성을 쓰고, 심사
(`optimizer_agent.risk_axes`)는 **구성체 전체**를 쓴다. 탐색을 구성체 범위로 맞추면
같은 잔기 치환의 효과가 **길이로 희석**되어 GA 가 내려갈 방향을 잃는다 — 이것이
§43-2 의 **역행 91/395** 를 설명한다는 가설.

🔴 먼저 논리를 정리한다 — **희석 자체는 탐색을 바꿀 수 없다**
-----------------------------------------------------------
`pareto_scrub` 은 목적 벡터를 **Pareto 지배**로만 비교한다. 지배는 각 축의 **순서**에만
의존하고 **크기에는 의존하지 않는다.** 따라서 구성체 범위가 축을 단순히 **배율로
축소**하기만 한다면 탐색 결과는 **한 치도 안 바뀐다.** 희석이 문제가 되려면 순서를
바꿔야 하고, 그 경로는 둘뿐이다:

  ⒜ **반전** — 셔틀 범위에서 A < B 인데 구성체 범위에서 A > B
  ⒝ **동점 생성** — `_instability()` 가 **소수 1자리로 반올림**한다
     (`developability_opt.py:168`). 배율로 줄인 뒤 반올림하면 서로 달랐던 값이
     **같은 칸에 떨어진다** → GA 가 그 쌍에서 기울기를 잃는다.

⒝ 가 이 가설의 실제 기전 후보다. **희석비가 아니라 ⒜⒝ 의 발생률을 재야 한다.**

재는 것 (전부 결정론 · LLM 0회)
------------------------------
1. **희석비** — 치환 하나의 ΔII(셔틀) 대 ΔII(구성체). 순수 길이 효과면
   `L_셔틀 / L_구성체` 에 맞아야 한다. 맞는지 본다(기전 확인).
2. **순서 보존** — Spearman ρ, **반전쌍 비율**, **새로 생긴 동점쌍 비율**.
   이것이 판정이다.
3. **이진 축 교차** — 심사는 `II < 40` 이진이다(`II_STABLE`). 범위를 바꾸면
   **판정이 뒤집히는** 변이체 수.
4. 🔴 **대응 널** — 화물을 **조성 보존 셔플**로 바꿔 전부 다시 잰다. 효과가 같으면
   원인은 **화물의 정체가 아니라 길이**다(이 프로젝트의 중심 규율: 같은 크기의
   무작위와 비교).
5. **링커 3종으로 길이를 바꿔** 1~3 을 반복 — 희석이 길이의 함수인지 직접 본다.

사용: /usr/bin/python3 scripts/null_instability_dilution.py [--cargo GSNKGAIIGLM]
"""
from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import DEFAULT_CARGO, SHUTTLES, VALID_AMINO_ACIDS  # noqa: E402
from core.developability_opt import _frozen_idx, _instability  # noqa: E402
from core.optimizer_agent import II_STABLE  # noqa: E402

LINKERS = {"직접융합": "", "EAAAK": "EAAAK", "A(EAAAK)3A": "AEAAAKEAAAKEAAAKA"}


def _spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def variants(shuttle: str):
    """비동결 위치의 **모든 단일 치환**. 결정론."""
    frozen = _frozen_idx(shuttle)
    out = []
    for i, c0 in enumerate(shuttle):
        if i in frozen:
            continue
        for aa in VALID_AMINO_ACIDS:
            if aa == c0:
                continue
            out.append(shuttle[:i] + aa + shuttle[i + 1:])
    return out


def measure(shuttle: str, cargo: str, linker: str):
    vs = variants(shuttle)
    pre = cargo + linker
    s0, c0 = _instability(shuttle), _instability(pre + shuttle)
    ii_s = [_instability(v) for v in vs]
    ii_c = [_instability(pre + v) for v in vs]
    ds = [a - s0 for a in ii_s]
    dc = [a - c0 for a in ii_c]

    ratios = [b / a for a, b in zip(ds, dc) if abs(a) > 1e-9]
    ratios.sort()
    med = ratios[len(ratios) // 2] if ratios else float("nan")

    n = len(vs)
    inv = tie_new = ordered = 0
    for i in range(n):
        for j in range(i + 1, n):
            if ii_s[i] == ii_s[j]:
                continue
            ordered += 1
            if ii_c[i] == ii_c[j]:
                tie_new += 1
            elif (ii_s[i] < ii_s[j]) != (ii_c[i] < ii_c[j]):
                inv += 1
    flip = sum(1 for a, b in zip(ii_s, ii_c) if (a < II_STABLE) != (b < II_STABLE))
    return {
        "n": n, "L_s": len(shuttle), "L_c": len(pre + shuttle),
        "len_ratio": len(shuttle) / len(pre + shuttle),
        "ratio_med": med, "rho": _spearman(ii_s, ii_c),
        "ordered": ordered, "inv": inv, "tie_new": tie_new, "flip": flip,
        "s0": s0, "c0": c0,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cargo", default=DEFAULT_CARGO)
    a = ap.parse_args()
    cargo = a.cargo
    rng = random.Random(0)
    sh_cargo = "".join(rng.sample(list(cargo), len(cargo)))   # 조성 보존 셔플(널)

    print("=" * 78)
    print("§43-6 차단 사유 #2 — 구성체 불안정성의 **길이 희석**이 탐색의 기울기를 죽이나")
    print("=" * 78)
    print(f"화물 `{cargo}` (L={len(cargo)}) · 대응 널 = 조성 보존 셔플 `{sh_cargo}`")
    print(f"심사 임계 II_STABLE = {II_STABLE} (이진)")
    print("🔴 Pareto 지배는 **순서**에만 의존한다 — 희석이 배율이기만 하면 탐색은 안 바뀐다.")
    print("   판정은 희석비가 아니라 **반전 · 새 동점 · 이진 교차**다.\n")

    names = sorted(SHUTTLES)
    for lname, lseq in LINKERS.items():
        print("-" * 78)
        print(f"[링커 `{lname}`]  구성체 = 화물({len(cargo)}) + 링커({len(lseq)}) + 셔틀")
        print("-" * 78)
        print(f"  {'셔틀':<16} {'n':>5} {'길이비':>6} {'희석비':>6} {'ρ':>6} "
              f"{'반전':>10} {'새동점':>12} {'이진교차':>8}")
        agg = {"ordered": 0, "inv": 0, "tie_new": 0, "flip": 0, "n": 0}
        for nm in names:
            seq = SHUTTLES[nm]["seq"]
            m = measure(seq, cargo, lseq)
            for k in agg:
                agg[k] += m[k]
            print(f"  {nm:<16} {m['n']:>5} {m['len_ratio']:>6.2f} {m['ratio_med']:>6.2f} "
                  f"{m['rho']:>6.3f} {m['inv']:>5}/{m['ordered']:<5} "
                  f"{m['tie_new']:>6}/{m['ordered']:<5} {m['flip']:>8}")
        print(f"  {'합계':<16} {agg['n']:>5} {'':>6} {'':>6} {'':>6} "
              f"{agg['inv']:>5}/{agg['ordered']:<5} "
              f"{agg['tie_new']:>6}/{agg['ordered']:<5} {agg['flip']:>8}")
        print(f"    → 반전 {100 * agg['inv'] / max(1, agg['ordered']):.3f}% · "
              f"새 동점 {100 * agg['tie_new'] / max(1, agg['ordered']):.3f}% · "
              f"이진 교차 {agg['flip']}/{agg['n']}")

    # ── 대응 널 ──
    print("\n" + "=" * 78)
    print("🔴 대응 널 — 화물을 **조성 보존 셔플**로 바꿔 같은 것을 잰다 (링커 `EAAAK`)")
    print("=" * 78)
    for tag, cg in (("실제 화물", cargo), ("셔플 널", sh_cargo)):
        agg = {"ordered": 0, "inv": 0, "tie_new": 0, "flip": 0, "n": 0}
        for nm in names:
            m = measure(SHUTTLES[nm]["seq"], cg, "EAAAK")
            for k in agg:
                agg[k] += m[k]
        print(f"  {tag:<10} 반전 {agg['inv']:>5}/{agg['ordered']} "
              f"({100 * agg['inv'] / max(1, agg['ordered']):.3f}%) · "
              f"새 동점 {agg['tie_new']:>5}/{agg['ordered']} "
              f"({100 * agg['tie_new'] / max(1, agg['ordered']):.3f}%) · "
              f"이진 교차 {agg['flip']}/{agg['n']}")
    print("\n  → 두 줄이 같으면 효과의 원인은 **화물의 정체가 아니라 길이**다.")
    print("\n재현: /usr/bin/python3 scripts/null_instability_dilution.py --cargo " + cargo)


if __name__ == "__main__":
    main()
