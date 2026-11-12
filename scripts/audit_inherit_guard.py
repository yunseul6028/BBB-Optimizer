#!/usr/bin/env /usr/bin/python3
"""큐 #15 — 근거 등급 **상속**이 기능 보존을 보장하나?

`core.evidence.inherit()` 는 설계체에 **문헌 근거 등급**을 준다. 등급은 이 시스템에서
유일한 양(陽)의 신호이고 사전식 순위의 **1순위**다(`optimizer_agent.rank_key`). 따라서
거짓 상속은 이 시스템이 할 수 있는 가장 비싼 오류다.

이 스크립트는 네 가지를 **측정**한다(수정은 하지 않는다 — 문서 모드):
  A1. 결정적 1사례 — Angiopep-7 = K10R/K15R (문헌상 뇌 전달 ~9.5× 낮음).
      유사도는 얼마이고 `inherit` 은 등급을 주는가?
  A2. 문헌 고정 위치(K10·K15) 파괴 변이의 **전수** 상속률 (단일·이중 치환 전수).
  A3. 라이브러리 12종 전체 — 치환 k 를 0.70 임계 양쪽으로 쓸어 보며
      (상속 여부) × (고정 위치 손상 여부) 교차표.
  A4. 길이 변화 편집 — `conserved_positions()` 가 () 를 돌려주는 경로의 빈도.

실행:  /usr/bin/python3 scripts/audit_inherit_guard.py
⚠️ `python3 -I` 로 돌리지 말 것(Biopython 미탐지 → 엔진 조용한 퇴화).
시드: SEED 상수(기본 20261112). 서브프로세스 0회(전부 BLOSUM 로컬 계산).
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import SHUTTLES  # noqa: E402
from core.evidence import (  # noqa: E402
    CONSERVED_POSITIONS,
    INHERIT_SIM_MIN,
    conserved_positions,
    inherit,
    similarity,
)

SEED = 20261112
_STD = "ACDEFGHIKLMNPQRSTVWY"
ANG2 = SHUTTLES["Angiopep-2"]["seq"]          # TFFYGGSRGKRNNFKTEEY
ANG7 = ANG2[:9] + "R" + ANG2[10:14] + "R" + ANG2[15:]   # K10R / K15R
PINNED = CONSERVED_POSITIONS["Angiopep-2"]    # (10, 15)


def damaged_pins(parent_name: str, var: str) -> tuple[int, ...]:
    """부모의 문헌 고정 위치 중 **실제로 바뀐** 1-based 위치들.

    `conserved_positions` 와 달리 길이가 달라도 '손상'으로 본다(길이가 변하면 그 위치의
    잔기가 보존됐다고 말할 수 없다). 즉 **진실**(문헌이 집어 준 자리가 살아 있나)이고,
    `conserved_positions` 는 **시스템이 지키려는 집합**이다. 둘의 차이가 이 감사의 대상.
    """
    pins = CONSERVED_POSITIONS.get(parent_name, ())
    if not pins:
        return ()
    ref = SHUTTLES[parent_name]["seq"]
    if len(var) != len(ref):
        return pins                      # 길이 변화 = 위치 대응 자체가 깨짐
    return tuple(p for p in pins if var[p - 1] != ref[p - 1])


def hdr(t: str) -> None:
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ---------------------------------------------------------------------------
def a1_angiopep7() -> None:
    hdr("A1. 결정적 1사례 — Angiopep-7 (K10R/K15R, 문헌 음성 대조군)")
    print(f"Angiopep-2 : {ANG2}")
    print(f"Angiopep-7 : {ANG7}")
    sim = similarity(ANG7, ANG2)
    tier, parent, isim = inherit(ANG7)
    cp = conserved_positions(ANG7)
    print(f"\n유사도(BLOSUM62 전역·자기정규화) = {sim:.4f}   (임계 {INHERIT_SIM_MIN})")
    print(f"inherit()  → 등급 {tier!r} · 부모 {parent!r} · 유사도 {isim}")
    print(f"conserved_positions(Angiopep-7) = {cp}  ← 시스템이 '지켜야 한다'고 보는 자리")
    print(f"실제 손상된 문헌 고정 위치 = {damaged_pins('Angiopep-2', ANG7)}")
    print(f"\n판정: 상속 허용? {'예' if tier else '아니오'} / 문헌 고정 위치 파괴? "
          f"{'예' if damaged_pins('Angiopep-2', ANG7) else '아니오'}")
    # 참고: 보존적 치환(K→R)의 BLOSUM62 점수가 왜 유사도를 거의 깎지 않는지
    print(f"\n참고 — 같은 자리를 **비보존적**으로 바꾼 경우의 유사도:")
    for aa in ("D", "E", "P", "G", "W"):
        v = ANG2[:9] + aa + ANG2[10:14] + aa + ANG2[15:]
        t2, _, s2 = inherit(v)
        print(f"   K10{aa}/K15{aa}: 유사도 {s2:.3f} → 등급 {t2!r}")


# ---------------------------------------------------------------------------
def a2_exhaustive_pin_breakers() -> None:
    hdr("A2. 🔴 결정적 비율 — 문헌 고정 위치를 파괴하면서 등급을 받는가 (전수)")
    rows = []

    # (a) 단일 치환: 위치 10 또는 15, 19개 대체 잔기씩
    single = []
    for p in PINNED:
        for aa in _STD:
            if aa == ANG2[p - 1]:
                continue
            v = ANG2[:p - 1] + aa + ANG2[p:]
            tier, parent, sim = inherit(v)
            single.append((f"K{p}{aa}", v, tier, parent, sim))
    # (b) 이중 치환: 10 과 15 동시
    double = []
    for a in _STD:
        if a == ANG2[9]:
            continue
        for b in _STD:
            if b == ANG2[14]:
                continue
            v = ANG2[:9] + a + ANG2[10:14] + b + ANG2[15:]
            tier, parent, sim = inherit(v)
            double.append((f"K10{a}/K15{b}", v, tier, parent, sim))

    for name, data in (("단일 치환(K10* 또는 K15*)", single),
                       ("이중 치환(K10*/K15* 동시)", double)):
        n = len(data)
        granted = [d for d in data if d[2]]
        wrong_parent = [d for d in granted if d[3] != "Angiopep-2"]
        sims = sorted(d[4] for d in data)
        print(f"\n{name}:  n = {n}")
        print(f"   등급을 받은 것(상속 성공) = {len(granted)}/{n} = "
              f"{100 * len(granted) / n:.1f}%")
        print(f"   그중 문헌 고정 위치 파괴 = {len(granted)}/{len(granted)} = 100.0% "
              f"(정의상 전부 — 설계가 K10/K15 치환이다)")
        print(f"   부모 오지목 = {len(wrong_parent)}")
        print(f"   유사도 범위 = {sims[0]:.3f} ~ {sims[-1]:.3f} "
              f"(중위 {sims[n // 2]:.3f}) / 임계 {INHERIT_SIM_MIN}")
        below = [d for d in data if d[4] < INHERIT_SIM_MIN]
        print(f"   임계 미달(= 상속 차단된 것) = {len(below)}  "
              + (", ".join(d[0] for d in below[:10]) if below else "없음"))
        rows.append((name, n, len(granted)))

    tot_n = sum(r[1] for r in rows)
    tot_g = sum(r[2] for r in rows)
    print(f"\n🔴 합계: 문헌 고정 위치(K10·K15)를 파괴한 전수 변이 {tot_n}개 중 "
          f"**{tot_g}개({100 * tot_g / tot_n:.1f}%)가 등급을 상속받았다.**")
    print("   (`inherit()` 은 CONSERVED_POSITIONS 를 조회하지 않는다 — 구조적으로 100%)")


# ---------------------------------------------------------------------------
def a3_sweep() -> None:
    hdr("A3. 라이브러리 12종 — 0.70 임계를 가로지르는 치환 스윕 (k = 치환 개수)")
    rng = random.Random(SEED)
    reps = 40
    print(f"시드 {SEED} · 각 (셔틀, k) 당 무작위 치환 {reps}회 · 치환은 서로 다른 위치")
    print("\n셔틀 | L | k | 유사도 중위 | 상속률 | 고정위치손상/상속 | conserved_positions")
    print("---|---|---|---|---|---|---")
    total = {"n": 0, "granted": 0, "pin_damaged_and_granted": 0, "has_pins": 0}
    for name, meta in SHUTTLES.items():
        ref = meta["seq"]
        L = len(ref)
        pins = CONSERVED_POSITIONS.get(name, ())
        ks = sorted({1, 2, 3, max(1, L // 4), max(1, L // 3), max(1, L // 2),
                     max(1, int(L * 0.6)), max(1, int(L * 0.75))})
        for k in ks:
            if k > L:
                continue
            sims, granted, pin_hit_granted, cp_empty = [], 0, 0, 0
            for _ in range(reps):
                idx = rng.sample(range(L), k)
                s = list(ref)
                for i in idx:
                    alts = [a for a in _STD if a != ref[i]]
                    s[i] = rng.choice(alts)
                v = "".join(s)
                tier, parent, sim = inherit(v)
                sims.append(sim)
                dmg = damaged_pins(name, v)
                if tier:
                    granted += 1
                    if dmg:
                        pin_hit_granted += 1
                if not conserved_positions(v):
                    cp_empty += 1
                total["n"] += 1
                total["granted"] += 1 if tier else 0
                if tier and dmg:
                    total["pin_damaged_and_granted"] += 1
                if pins:
                    total["has_pins"] += 1
            sims.sort()
            print(f"{name} | {L} | {k} | {sims[len(sims) // 2]:.3f} | "
                  f"{granted}/{reps} | "
                  f"{(str(pin_hit_granted) + '/' + str(granted)) if pins else '—(고정위치 없음)'} | "
                  f"{reps - cp_empty}/{reps} 비어있지않음")
    print(f"\n합계: 변이체 {total['n']}개 · 상속 성공 {total['granted']} "
          f"({100 * total['granted'] / total['n']:.1f}%)")
    print(f"      그중 '문헌 고정 위치 손상 + 등급 상속' = "
          f"{total['pin_damaged_and_granted']}")
    print(f"      ⚠️ 문헌 고정 위치가 **선언된** 셔틀은 12종 중 "
          f"{len(CONSERVED_POSITIONS)}종({', '.join(CONSERVED_POSITIONS)})뿐 — "
          f"나머지 11종은 '보존할 필요 없음'이 아니라 '어디가 중요한지 모른다'다.")


# ---------------------------------------------------------------------------
def a4_length_change() -> None:
    hdr("A4. 길이 변화 편집 — 위치 가드가 **완전히 사라지는** 경로")
    print("`conserved_positions()` 는 길이가 다르면 () 를 돌려준다 → 동결 집합이 빈다.")
    print("그런데 `inherit()` 은 길이와 무관하게 유사도만 보므로 등급은 계속 나온다.\n")

    # (a) 전수 단일 결손 (19개)
    dels = []
    for i in range(len(ANG2)):
        v = ANG2[:i] + ANG2[i + 1:]
        tier, parent, sim = inherit(v)
        dels.append((i + 1, v, tier, parent, sim, conserved_positions(v)))
    gd = [d for d in dels if d[2]]
    print(f"(a) 단일 결손 전수: n = {len(dels)}")
    print(f"    등급 상속 = {len(gd)}/{len(dels)} = {100 * len(gd) / len(dels):.1f}%")
    print(f"    conserved_positions() == () 인 것 = "
          f"{sum(1 for d in dels if not d[5])}/{len(dels)} = 100.0%")
    print(f"    유사도 범위 = {min(d[4] for d in dels):.3f} ~ {max(d[4] for d in dels):.3f}")
    for p in PINNED:
        d = dels[p - 1]
        print(f"    · 문헌 고정 위치 K{p} 자체를 **삭제**: 유사도 {d[4]:.3f} → 등급 {d[2]!r} "
              f"· 동결집합 {d[5]}")

    # (b) 전수 단일 삽입 (20 잔기 × 20 위치)
    ins = []
    for i in range(len(ANG2) + 1):
        for aa in _STD:
            v = ANG2[:i] + aa + ANG2[i:]
            tier, parent, sim = inherit(v)
            ins.append((i, aa, tier, sim, conserved_positions(v)))
    gi = [d for d in ins if d[2]]
    print(f"\n(b) 단일 삽입 전수: n = {len(ins)}")
    print(f"    등급 상속 = {len(gi)}/{len(ins)} = {100 * len(gi) / len(ins):.1f}%")
    print(f"    conserved_positions() == () 인 것 = "
          f"{sum(1 for d in ins if not d[4])}/{len(ins)} = "
          f"{100 * sum(1 for d in ins if not d[4]) / len(ins):.1f}%")
    print(f"    유사도 범위 = {min(d[3] for d in ins):.3f} ~ {max(d[3] for d in ins):.3f}")

    # (c) 삽입 + 고정위치 치환 조합 — 가드가 없으니 둘을 같이 할 수 있나
    print("\n(c) 길이 변화 + 고정 위치 파괴 **동시**: 삽입 1개 + K10R/K15R")
    combo_granted = 0
    combo = []
    for i in range(len(ANG7) + 1):
        for aa in "GSAEK":
            v = ANG7[:i] + aa + ANG7[i:]
            tier, _p, sim = inherit(v)
            combo.append((i, aa, tier, sim, conserved_positions(v)))
            combo_granted += 1 if tier else 0
    print(f"    n = {len(combo)} · 등급 상속 = {combo_granted}/{len(combo)} = "
          f"{100 * combo_granted / len(combo):.1f}%")
    print(f"    동결집합이 빈 것 = {sum(1 for c in combo if not c[4])}/{len(combo)}")
    print(f"    유사도 범위 = {min(c[3] for c in combo):.3f} ~ "
          f"{max(c[3] for c in combo):.3f}")

    # (d) 실제 파이프라인이 이 경로를 쓰는가
    print("\n(d) 파이프라인에서 이 경로가 실제로 밟히는가:")
    print("    · `developability_opt._mutate` = **단일 치환만**(길이 보존) → pareto_scrub /")
    print("      scrub_liabilities 는 길이 변화 편집을 만들지 **않는다** → A4 경로 0%.")
    print("    · 그러나 Gemini 에이전트의 `design_candidate` 는 **LLM 이 서열을 직접 쓴다** —")
    print("      길이 제한 게이트가 없으므로(상속 게이트는 유사도만 본다) 이 경로가 열려 있다.")


def a5_separability() -> None:
    """임계값을 올려서 고칠 수 있는 문제인가 — 두 분포의 겹침."""
    hdr("A5. 임계값으로 고칠 수 있나 — '정당한 변이체' vs '핀 파괴 변이체' 유사도 분포")
    from core.developability_opt import pareto_scrub

    # 정당한 변이체 = pareto_scrub 산출물(K10/K15 동결되어 있으므로 **정의상 정당**)
    legit = set()
    for sd in (0, 1, 7, 42, 2022, 31337):
        for m in pareto_scrub(ANG2, rounds=25, pop=24, rmt_floor=0.85, seed=sd):
            if m["seq"] != ANG2:
                legit.add(m["seq"])
    ls = sorted(similarity(s, ANG2) for s in legit)
    broke = [s for s in legit if len(s) == len(ANG2)
             and (s[9] != ANG2[9] or s[14] != ANG2[14])]

    pin = []
    for p in PINNED:
        for aa in _STD:
            if aa != ANG2[p - 1]:
                pin.append(ANG2[:p - 1] + aa + ANG2[p:])
    for a in _STD:
        if a == ANG2[9]:
            continue
        for b in _STD:
            if b == ANG2[14]:
                continue
            pin.append(ANG2[:9] + a + ANG2[10:14] + b + ANG2[15:])
    ps = sorted(similarity(s, ANG2) for s in pin)

    print(f"정당한 개발성 변이체(pareto_scrub 6시드, K10/K15 동결): n={len(ls)}")
    print(f"   유사도 {ls[0]:.3f} ~ {ls[-1]:.3f} (중위 {ls[len(ls) // 2]:.3f}) · "
          f"핀 손상 {len(broke)}개(동결 작동 확인)")
    print(f"핀 파괴 변이체(전수): n={len(ps)}")
    print(f"   유사도 {ps[0]:.3f} ~ {ps[-1]:.3f} (중위 {ps[len(ps) // 2]:.3f})")
    ov = [x for x in ps if x >= ls[0]]
    print(f"\n🔴 겹침: 핀 파괴 변이체 중 '정당 변이체 최소 유사도({ls[0]:.3f})' 이상 = "
          f"{len(ov)}/{len(ps)} = {100 * len(ov) / len(ps):.1f}%")
    need = similarity(ANG7, ANG2)
    blocked = [x for x in ls if x < need]
    print(f"   Angiopep-7(={need:.3f})을 막으려면 임계 > {need:.3f} 가 필요하고, "
          f"그러면 **정당** 변이체 {len(blocked)}/{len(ls)} = "
          f"{100 * len(blocked) / len(ls):.1f}% 가 함께 탈락한다.")
    print("   → 두 분포가 **완전히 역전**돼 있어 어떤 임계값도 분리하지 못한다. "
          "유사도는 애초에 분리 가능한 축이 아니다(CONSERVED_POSITIONS 독스트링과 같은 결론).")


if __name__ == "__main__":
    print(f"[seed={SEED}] INHERIT_SIM_MIN={INHERIT_SIM_MIN} · "
          f"CONSERVED_POSITIONS={CONSERVED_POSITIONS}")
    a1_angiopep7()
    a2_exhaustive_pin_breakers()
    a3_sweep()
    a4_length_change()
    a5_separability()
