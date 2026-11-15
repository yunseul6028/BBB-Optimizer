#!/usr/bin/env /usr/bin/python3
"""근거 등급 **상속 게이트**의 기능 보존 대응 널 (§55).

왜 — `docs/protocol.md` §8 의 "**등급 상속의 기능 보존 미심사.** M9 는 음성 대조만 봤다"
를 메운다. §31 이 게이트를 감사해 `_pins_intact` 를 심었지만 그 독스트링 자신이
**"실효 범위는 1/12"**(문헌 핀이 있는 셔틀은 Angiopep-2 하나)라고 적고 있다.
그래서 묻는다 — **핀이 없는 11종에서, 등급을 상속받으면서 기능을 얼마나 깰 수 있나.**

기능의 조작적 정의: `binding.shuttle_similarity(seq).score`(= `rmt_sim`). 이 값은
탐색이 **기능 보존 하한**으로 쓰는 바로 그 양이다(`pareto_scrub(rmt_floor=0.85)`).
⚠️ 이 대리값의 한계는 §32 가 판정했다 — THR·Tet1 은 참조셋 미등재라 0.0 이고, 그 둘에서는
"기능 손실"을 정의할 수 없어 **제외**한다. 즉 이 측정은 축을 **그 축 자신의 기준으로**
심사한다(외부 라벨이 아니다 — §55-4 한계).

핵심 교차 수치: **상속 통과 ∧ rmt_sim < 0.85** = 게이트가 통과시킨 기능 상실체.
"""
from __future__ import annotations

import random
import sys
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from core.binding import shuttle_similarity
from core.config import SHUTTLES
from core.developability_opt import pareto_scrub
from core.evidence import INHERIT_SIM_MIN, inherit

AA = "ACDEFGHIKLMNPQRSTVWY"
FLOOR = 0.85          # pareto_scrub 의 기능 보존 하한과 같은 값
SEED = 20261115
N_RANDOM = 40         # 셔틀 × k 당 무작위 변이체 수
K_RANGE = (1, 2, 3, 4, 5, 6)


def rmt(seq: str) -> float:
    return shuttle_similarity(seq).score or 0.0


def inherited(seq: str) -> bool:
    tier, _parent, _sim = inherit(seq)
    return tier is not None


def main() -> None:
    rng = random.Random(SEED)
    usable = {n: v["seq"] for n, v in SHUTTLES.items() if rmt(v["seq"]) > 0.0}
    skipped = [n for n in SHUTTLES if n not in usable]

    print("=" * 86)
    print(f"§55 상속 게이트 기능 보존 널 — 시드 {SEED} · 기능 하한 {FLOOR} · 임계 {INHERIT_SIM_MIN}")
    print(f"대상 셔틀 {len(usable)}종 · 제외 {skipped} (rmt_sim = 0.0, §32 명부 미등재)")
    print("=" * 86)

    # ── ① 단일 치환 전수 ────────────────────────────────────────────────
    print("\n── ① 단일 치환 **전수** — 한 수로 등급을 지킨 채 기능을 얼마나 깨나 ──")
    print(f"{'셔틀':<18}{'원본 rmt':>9}{'전수 n':>8}{'상속':>10}{'상속∧붕괴':>12}{'최저 rmt(상속자)':>18}")
    tot_n = tot_inh = tot_bad = 0
    worst_rows = []
    for name, seq in usable.items():
        base = rmt(seq)
        n = inh = bad = 0
        lo, lo_seq = 1.0, None
        for i in range(len(seq)):
            for a in AA:
                if a == seq[i]:
                    continue
                mut = seq[:i] + a + seq[i + 1:]
                n += 1
                if not inherited(mut):
                    continue
                inh += 1
                r = rmt(mut)
                if r < lo:
                    lo, lo_seq = r, (i + 1, seq[i], a)
                if r < FLOOR:
                    bad += 1
        tot_n += n; tot_inh += inh; tot_bad += bad
        worst_rows.append((name, base, lo, lo_seq))
        print(f"{name:<18}{base:>9.3f}{n:>8}{inh:>7}/{n:<3}{bad:>9} ({100*bad/max(inh,1):>4.1f}%)"
              f"{lo:>12.3f} {str(lo_seq):>}")
    print(f"{'합계':<18}{'':>9}{tot_n:>8}{tot_inh:>7}/{tot_n:<3}"
          f"{tot_bad:>9} ({100*tot_bad/max(tot_inh,1):>4.1f}%)")

    # ── ② 대응 널 — 같은 편집 수 k 의 무작위 치환 ───────────────────────
    print(f"\n── ② 대응 널 — 무작위 k 치환 {N_RANDOM}개/셔틀 ──")
    print(f"{'k':>3}{'생성':>8}{'상속':>12}{'상속∧붕괴':>14}{'붕괴(전체)':>13}{'상속자 rmt 중위':>16}")
    for k in K_RANGE:
        n = inh = bad = broke = 0
        inh_r = []
        for name, seq in usable.items():
            for _ in range(N_RANDOM):
                idx = rng.sample(range(len(seq)), min(k, len(seq)))
                lst = list(seq)
                for i in idx:
                    lst[i] = rng.choice([a for a in AA if a != seq[i]])
                mut = "".join(lst)
                n += 1
                r = rmt(mut)
                if r < FLOOR:
                    broke += 1
                if inherited(mut):
                    inh += 1
                    inh_r.append(r)
                    if r < FLOOR:
                        bad += 1
        print(f"{k:>3}{n:>8}{inh:>7}/{n:<4}{bad:>9} ({100*bad/max(inh,1):>4.1f}%)"
              f"{broke:>8} ({100*broke/n:>4.1f}%){median(inh_r) if inh_r else float('nan'):>16.3f}")

    # ── ③ 양성 대조 — pareto_scrub 변이체는 전부 통과해야 한다 ──────────
    print("\n── ③ 양성 대조 — `pareto_scrub` 변이체(기능 하한 제약 통과) ──")
    pos = fail = 0
    errs = []
    for name, seq in usable.items():
        try:
            front = pareto_scrub(seq, seed=0)
        except Exception as exc:  # noqa: BLE001
            errs.append(f"{name}: {type(exc).__name__}")
            continue
        for row in front:
            if row["seq"] == seq:
                continue
            pos += 1
            if not inherited(row["seq"]):
                fail += 1
    print(f"  변이체 {pos}개 중 상속 거부 **{fail}개** (오탈락률 {100*fail/max(pos,1):.1f}%)")
    for e in errs:
        print(f"  ⚠️ {e}")

    print("\n해석: ①②의 '상속∧붕괴'가 **게이트가 통과시킨 기능 상실체**다. "
          "③이 0 이 아니면 게이트가 정당한 변이체를 막는 것이다.")


if __name__ == "__main__":
    main()
