#!/usr/bin/env python3
"""§47 보조 — `null_toxicity_family4.py` 가 남긴 세 질문만 좁게 측정한다.

본 측정(`null_toxicity_family4.py`)이 낸 세 수치는 해석이 열려 있었다:

1. **"독소 조각 검출 8/24"** 에서 **어느 조각이 잡히고 어느 조각이 새는지** 모른다.
   조각별로 쪼개야 "게이트가 독소를 알아보나"와 "게이트가 소수성 덩어리를 알아보나"가
   갈린다.
2. **"무작위 29/120 vs 독소 8/24"** 의 초과가 우연인지 안 쟀다. 🔴 대응 널의 핵심은
   **초과분이지 절대값이 아니다**(프로토콜 §5: "널을 나란히 쓴다").
3. **Bombolitin-I 가 0.380 으로 임계에 정확히 앉아 통과했다.** §42-4b 가 `off_target`
   에서 발견한 **눈금 정착**과 같은 모양인지 — 즉 ML Score 가 양자화돼 있고 **0.38 이
   하필 눈금 위**인지 — 를 확인해야 한다. 그렇다면 `>` 와 `>=` 중 어느 비교를 쓰는지가
   설계 결정이 된다.

표본은 좁다(≈ 350 서열, 1분 미만). 🔴 재는 양은 여기서도 **ToxinPred3 ML Score** 이고
독성이 아니다.
"""

from __future__ import annotations

import random
import sys
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import LINKER_LIBRARY, SHUTTLES, get_settings  # noqa: E402
from core.toxicity import get_toxicity_predictor  # noqa: E402

SEED = 0
AA = "ACDEFGHIKLMNPQRSTVWY"
THRESH = 0.38
CARGOS = {"Aβ25-35": "GSNKGAIIGLM", "DSIP": "WAGGDASGE"}

TOXIC_FULL = {
    "Melittin":     "GIGAVLKVLTTGLPALISWIKRKRQQ",
    "Pardaxin":     "GFFALIPKIISSPLFKTLLSAVGSALSSSGGQE",
    "Apamin":       "CNCKAPETALCARRCQQH",
    "Bombolitin-I": "IKITTMLAKLGKVLAHV",
}
# 본 측정과 **같은** 6셀(링커 길이 17 통과 셀을 점수 오름차순으로 균등 추출)
PICK = [("Aβ25-35", "Leptin30"), ("DSIP", "ApoE(159-167)2"), ("DSIP", "T7"),
        ("DSIP", "TAT"), ("Aβ25-35", "RVG29"), ("Aβ25-35", "CDX")]


def fisher_one_sided(a: int, b: int, c: int, d: int) -> float:
    """2×2 단측 Fisher 정확검정. (a: 독소 탈락, b: 독소 통과, c: 널 탈락, d: 널 통과)
    반환 = P(독소 탈락 ≥ a). 정수 조합으로 계산해 부동소수 누적오차를 피한다."""
    n = a + b + c + d
    r1, c1 = a + b, a + c
    total = comb(n, c1)
    p = Fraction(0)
    for k in range(a, min(r1, c1) + 1):
        p += Fraction(comb(r1, k) * comb(n - r1, c1 - k), total)
    return float(p)


def main() -> None:
    pred = get_toxicity_predictor(get_settings())
    print("=" * 78)
    print("§47 보조 — 조각별 검출 · 초과분의 유의성 · 임계의 눈금 정착")
    print(f"예측기: {pred.source} · 임계 {THRESH} (비교는 `score > {THRESH}`)")
    print("=" * 78)

    # ── 1. 조각별 검출 ───────────────────────────────────────────────────
    print("\n## 1. 독소 조각 **조각별** 검출 — 어느 것이 새나\n")
    frags = {k: v[:17] for k, v in TOXIC_FULL.items()}
    for k, v in frags.items():
        print(f"- `{k}` 17aa 조각: `{v}`")
    seqs, keys = [], []
    for cn, sn in PICK:
        cg, ss = CARGOS[cn], SHUTTLES[sn]["seq"]
        for fk, fv in frags.items():
            seqs.append(cg + fv + ss)
            keys.append((cn, sn, fk))
    sc = [r.risk for r in pred.predict_many(seqs)]
    byfrag: dict[str, list[float]] = {k: [] for k in frags}
    for (cn, sn, fk), v in zip(keys, sc):
        byfrag[fk].append(v)
    print("\n| 조각 | 전장 ML Score(참조) | 삽입 후 점수 6셀 | **탈락(검출)** |")
    print("|---|---|---|---|")
    full_sc = [r.risk for r in pred.predict_many(list(TOXIC_FULL.values()))]
    for (fk, vals), fs in zip(byfrag.items(), full_sc):
        hit = sum(1 for v in vals if v > THRESH)
        mark = f"**{hit}/6**" + (" 🔴" if hit < 3 else "")
        print(f"| {fk} | {fs:.3f} | {' · '.join(f'{v:.3f}' for v in vals)} | {mark} |")

    # ── 2. 초과분의 유의성 ───────────────────────────────────────────────
    print("\n## 2. 초과분 — 독소 조각 vs 같은 길이 무작위 (계열 1 대조)\n")
    rng = random.Random(SEED)
    rseqs = []
    for cn, sn in PICK:
        cg, ss = CARGOS[cn], SHUTTLES[sn]["seq"]
        for _ in range(20):
            rseqs.append(cg + "".join(rng.choice(AA) for _ in range(17)) + ss)
    rsc = [r.risk for r in pred.predict_many(rseqs)]
    tox_fail = sum(1 for v in sc if v > THRESH)
    rnd_fail = sum(1 for v in rsc if v > THRESH)
    a, b = tox_fail, len(sc) - tox_fail
    c, d = rnd_fail, len(rsc) - rnd_fail
    p = fisher_one_sided(a, b, c, d)
    print(f"| | 탈락 | 전체 | 비율 |")
    print(f"|---|---|---|---|")
    print(f"| **독소 조각 삽입** | {a} | {a+b} | **{100*a/(a+b):.1f}%** |")
    print(f"| 같은 길이 무작위(널) | {c} | {c+d} | {100*c/(c+d):.1f}% |")
    print(f"\n**초과 {100*a/(a+b) - 100*c/(c+d):+.1f}pp · Fisher 단측 p = {p:.4f}**")
    print(f"평균 점수: 독소 {sum(sc)/len(sc):.3f} vs 무작위 {sum(rsc)/len(rsc):.3f}")

    # ── 3. 눈금 정착 ─────────────────────────────────────────────────────
    print("\n## 3. 🔴 ML Score 의 눈금과 임계 0.38 의 위치\n")
    cells, cseqs = [], []
    for cn, cg in CARGOS.items():
        for sn, sm in SHUTTLES.items():
            for ln, lm in LINKER_LIBRARY.items():
                cells.append((cn, sn, ln))
                cseqs.append(cg + (lm["seq"] or "") + sm["seq"])
    bsc = [r.risk for r in pred.predict_many(cseqs)]
    allv = sorted(set(bsc + sc + rsc + full_sc))
    diffs = sorted({round(b - a, 6) for a, b in zip(allv, allv[1:])})
    steps = Counter(round(b - a, 6) for a, b in zip(allv, allv[1:]))
    print(f"관측된 고유 점수 {len(allv)}개 · 최소 간격 **{diffs[0]}** · "
          f"가장 흔한 간격 {steps.most_common(3)}")
    at = sum(1 for v in bsc if abs(v - THRESH) < 1e-9)
    near = sum(1 for v in bsc if THRESH - 0.02 <= v <= THRESH + 0.02)
    print(f"\n라이브러리 전수 {len(bsc)}셀의 점수 분포:")
    print(f"  · **정확히 0.380 = {at}셀** → `score > 0.38` 이므로 **통과한다**")
    print(f"  · 0.360~0.400 구간 = {near}셀 ({100*near/len(bsc):.1f}%)")
    print(f"  · 탈락(>0.38) = {sum(1 for v in bsc if v > THRESH)}셀")
    print(f"\n임계 근방 점수값과 셀 수:")
    for v in [x for x in sorted(set(bsc)) if THRESH - 0.03 <= x <= THRESH + 0.03]:
        n = sum(1 for z in bsc if abs(z - v) < 1e-9)
        verdict = "탈락" if v > THRESH else "**통과**"
        print(f"  {v:.3f} : {n:3d}셀  → {verdict}")
    print("\n" + "=" * 78)
    print("판정은 docs/experiments.md §47 에 적는다.")
    print("=" * 78)


if __name__ == "__main__":
    main()
