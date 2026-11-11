"""off-target(선택성) 축의 **대응 널 심사** — 축 자체를 피고석에 세운다 (§26).

왜 필요한가
-----------
off-target 축의 유일한 검증은 `bench_selectivity.py` 의 **AUC 1.000(RMT>CPP)** 이다.
그런데 그 AUC 가 복원하는 라벨(RMT/CPP)은 **근거 등급이 쓰는 라벨과 같고**, 축의
계산식 자체가 RMT 참조셋과의 유사도(`rmt_sim`)를 읽는다 → **거의 순환**이다.
게다가 2026-11-11 목적함수 4축 전환 이후 이 축은 **탐색 방향까지** 정하고,
라이브러리 대표 셔틀(Angiopep-2)의 전선 생존을 **혼자 결정**한다(§24-6).
검증 부담이 커졌는데 검증은 순환 라벨 하나뿐이다.

대응 널의 원칙(`matched-null-principle`): **같은 크기의 무작위 편집과 비교**하고,
후보가 아니라 **축 자체**를 심사하고, **검출과 복구를 한 쌍**으로 본다.

네 가지 널
----------
A. **라벨 조회 분해** — AUC 1.000 이 물리화학(전하·GRAVY)에서 오나, `rmt_sim`(=라벨
   조회)에서 오나. 항을 하나씩 끄고 같은 AUC 를 다시 센다.
B. 🔴 **방향성 대응 널(결정적)** — RMT 셔틀에 **같은 크기 무작위 치환**(k=1,2,3)을
   주면 off 가 어디로 가나. 축이 "수용체 선택성 상실 위험"을 잰다면 무작위 편집은
   결합을 망치므로 off 가 **올라야** 한다. 반대로 off 가 **내려가면** 그 축은
   off-target 위험계가 아니라 **양이온성 측정기**다.
C. **순서 무관성** — 조성 보존 셔플(전하·GRAVY 가 **정확히** 불변, 순서만 섞음).
   off 가 안 바뀌면 축은 서열이 아니라 **조성**을 보고 있다.
D. **검출 + 복구 한 쌍** — ⑴ 검출: 알려진 기능 상실체 Angiopep-7(K10R/K15R)을
   Angiopep-2 와 구별하나. ⑵ 복구(해킹): CPP 의 양전하만 지운 "가짜 RMT"가 축을
   속여 선택적으로 승격되나.

실행: 프로젝트 루트에서 Biopython 이 설치된 파이썬으로 (서브프로세스 0회·결정론).
  $ /usr/bin/python3 scripts/null_selectivity.py
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import SHUTTLES
from core.selectivity import assess_selectivity
from scripts.bench_selectivity import HOLDOUT_CPP, auc

_STD = "ACDEFGHIKLMNPQRSTVWY"


# --------------------------------------------------------------------------
# 축의 세 항을 분해해 다시 조립한다 (core.selectivity.assess_selectivity 와 동일 식)
# --------------------------------------------------------------------------
def components(seq: str) -> tuple[float, float, float]:
    """(cationicity, lipophilicity, rmt_sim) — 축을 만드는 세 재료."""
    r = assess_selectivity(seq)
    cat = max(0.0, min(1.0, r.charge_density / 0.5)) if r.charge_density > 0 else 0.0
    lip = max(0.0, min(1.0, r.gravy / 0.8)) if r.gravy > 0 else 0.0
    return cat, lip, r.rmt_similarity


def off_variant(seq: str, *, use_physchem: bool = True, use_rmt: bool = True) -> float:
    """항을 켜고 끈 off-target. 전체 = `assess_selectivity().off_target_risk`."""
    cat, lip, rmt = components(seq)
    raw = 1 - (1 - cat) * (1 - lip) if use_physchem else 1.0
    if use_rmt:
        raw *= (1 - 0.6 * rmt)
    return round(max(0.0, min(1.0, raw)), 4)


def groups() -> tuple[dict, dict]:
    rmt = {n: SHUTTLES[n]["seq"] for n in SHUTTLES if "RMT" in SHUTTLES[n]["note"]}
    cpp = {n: SHUTTLES[n]["seq"] for n in SHUTTLES if "CPP" in SHUTTLES[n]["note"]}
    return rmt, {**cpp, **HOLDOUT_CPP}


# --------------------------------------------------------------------------
def null_a() -> None:
    print("\n=== 널 A — 라벨 조회 분해: AUC 1.000 은 어디서 오나 ===")
    rmt, cpp = groups()
    rs, cs = list(rmt.values()), list(cpp.values())
    arms = [
        ("전체 축 (물리화학 × rmt_sim)", dict(use_physchem=True, use_rmt=True)),
        ("물리화학만 (rmt_sim 제거)", dict(use_physchem=True, use_rmt=False)),
        ("rmt_sim 만 (물리화학 제거)", dict(use_physchem=False, use_rmt=True)),
    ]
    print(f"RMT {len(rs)} vs CPP {len(cs)} (홀드아웃 {len(HOLDOUT_CPP)} 포함)")
    print(f"{'조건':<32}{'RMT평균off':>11}{'CPP평균off':>11}{'AUC(선택성)':>13}")
    for label, kw in arms:
        a = [off_variant(s, **kw) for s in rs]
        b = [off_variant(s, **kw) for s in cs]
        # 선택성 = 1 - off 이므로 AUC(RMT 선택성 > CPP 선택성) = AUC(RMT off < CPP off)
        print(f"{label:<32}{sum(a)/len(a):>11.3f}{sum(b)/len(b):>11.3f}"
              f"{auc([-x for x in a], [-y for y in b]):>13.3f}")


def null_b(n: int = 120, seed: int = 0) -> None:
    print("\n=== 널 B 🔴 방향성 대응 널: 같은 크기 무작위 편집에서 off 가 어디로 가나 ===")
    print("축이 '수용체 선택성 상실 위험'을 재면 무작위 편집은 결합을 망치므로 off ↑ 가 기대값이다.")
    rmt, _ = groups()
    rng = random.Random(seed)
    print(f"{'셔틀':<16}{'k':>2}{'기준off':>8}{'평균Δoff':>10}{'상승':>6}{'하락':>6}{'동일':>6}")
    agg = {1: [0, 0, 0, 0.0], 2: [0, 0, 0, 0.0], 3: [0, 0, 0, 0.0]}
    nf = {1: [0, 0, 0, 0.0], 2: [0, 0, 0, 0.0], 3: [0, 0, 0, 0.0]}
    for name, seq in rmt.items():
        base = assess_selectivity(seq).off_target_risk
        for k in (1, 2, 3):
            up = dn = eq = 0
            tot = 0.0
            for _ in range(n):
                s = list(seq)
                for i in rng.sample(range(len(s)), min(k, len(s))):
                    s[i] = rng.choice(_STD)
                d = assess_selectivity("".join(s)).off_target_risk - base
                tot += d
                up += d > 1e-9
                dn += d < -1e-9
                eq += abs(d) <= 1e-9
            agg[k][0] += up; agg[k][1] += dn; agg[k][2] += eq; agg[k][3] += tot / n
            # ⚠️ 기준 off 가 바닥(≤0.06)인 셔틀은 **내려갈 데가 없다** → 상승 쪽으로 기운다.
            #    바닥 효과를 섞으면 "축이 올바른 방향으로 반응한다"는 **가짜 양성**이 된다.
            if base > 0.06:
                nf[k][0] += up; nf[k][1] += dn; nf[k][2] += eq; nf[k][3] += tot / n
            print(f"{name:<16}{k:>2}{base:>8.2f}{tot/n:>+10.3f}{up:>6}{dn:>6}{eq:>6}")
    n_nf = sum(1 for s in rmt.values() if assess_selectivity(s).off_target_risk > 0.06)
    print(f"\n{'합계':<24}{'k':>2}{'평균Δoff':>10}{'상승':>6}{'하락':>6}{'동일':>6}{'상승비':>8}")
    for k in (1, 2, 3):
        for lbl, (up, dn, eq, tot), m in ((f"RMT 전체 ({len(rmt)}종)", agg[k], len(rmt)),
                                          (f"바닥 제외 ({n_nf}종)", nf[k], n_nf)):
            share = up / (up + dn) if up + dn else 0.5
            print(f"{lbl:<24}{k:>2}{tot/m:>+10.3f}{up:>6}{dn:>6}{eq:>6}{share:>8.1%}")


def null_c(n: int = 200, seed: int = 0) -> None:
    print("\n=== 널 C — 순서 무관성: 조성 보존 셔플(전하·GRAVY 정확히 불변) ===")
    rmt, _ = groups()
    rng = random.Random(seed)
    print(f"{'셔틀':<16}{'기준off':>8}{'셔플평균':>9}{'변화한 셔플':>12}{'최대|Δ|':>8}")
    for name, seq in rmt.items():
        base = assess_selectivity(seq).off_target_risk
        vals = []
        for _ in range(n):
            s = list(seq)
            rng.shuffle(s)
            vals.append(assess_selectivity("".join(s)).off_target_risk)
        changed = sum(1 for v in vals if abs(v - base) > 1e-9)
        print(f"{name:<16}{base:>8.2f}{sum(vals)/len(vals):>9.3f}"
              f"{f'{changed}/{n}':>12}{max(abs(v - base) for v in vals):>8.2f}")


def null_d() -> None:
    print("\n=== 널 D — 검출 + 복구 한 쌍 ===")
    print("⑴ 검출: 알려진 기능 상실체를 구별하나")
    pairs = [
        ("Angiopep-2 (기능 O)", "TFFYGGSRGKRNNFKTEEY"),
        ("Angiopep-7 K10R/K15R (뇌 전달 ~9.5배↓)", "TFFYGGSRGRRNNFRTEEY"),
    ]
    for lbl, s in pairs:
        r = assess_selectivity(s)
        print(f"  {lbl:<42} off {r.off_target_risk:.2f} · 전하밀도 {r.charge_density:+.2f}"
              f" · rmt_sim {r.rmt_similarity:.2f}")

    print("⑵ 복구(해킹): CPP 의 양전하만 지우면 축이 속나")
    hacks = [
        ("TAT (원본·CPP)", "GRKKRRQRRRPPQ"),
        ("TAT 전하 제거 (K/R→Q)", "GQQQQQQQQQPPQ"),
        ("R9 (원본·CPP)", "RRRRRRRRR"),
        ("R9 전하 제거 (R→Q)", "QQQQQQQQQ"),
        ("MAP (원본·CPP)", "KLALKLALKALKAALKLA"),
        ("MAP 전하 제거 (K→Q)", "QLALQLALQALQAALQLA"),
    ]
    for lbl, s in hacks:
        r = assess_selectivity(s)
        print(f"  {lbl:<42} off {r.off_target_risk:.2f} · 위험 {r.risk_level} · {r.mechanism[:34]}")


def main() -> None:
    null_a()
    null_b()
    null_c()
    null_d()


if __name__ == "__main__":
    main()
