#!/usr/bin/env python3
"""§30-4a 비대칭 검사 — 기록된 에이전트 설계가 **현 GA 의 제약 안**에 있는가?

§30-4a: "off ±0.05 밴드는 그 세션들 이후에, 그 세션들 때문에 추가됐다. 에이전트
최종은 밴드를 위반한다(|Δ|=0.11) → 현재 GA 는 그것을 생성조차 할 수 없다."
→ 그렇다면 "B 가 A 를 지배"는 **규칙 밖 설계가 규칙 안 설계를 이긴 것**이 된다.

2판(수리 후)에도 같은 문제가 있는지 **직접 계산**한다. 유료 호출 0회
(off_target·shuttle_similarity 는 로컬, deepB3P·ToxinPred3 는 건드리지 않는다).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
sys.path.insert(0, "/Users/yunji/Claude/BBB-Optimizer")

from core.developability_opt import OFF_BAND, off_target
from core.binding import shuttle_similarity
from core.evidence import SHUTTLES, CONSERVED_POSITIONS, INHERIT_SIM_MIN

SEED_NAME = "Angiopep-2"
seed = SHUTTLES[SEED_NAME]["seq"] if isinstance(SHUTTLES[SEED_NAME], dict) else SHUTTLES[SEED_NAME]
off0 = off_target(seed)
pins = CONSERVED_POSITIONS[SEED_NAME]
FLOOR = 0.85   # pareto_scrub 의 rmt_floor 기본값

CASES = [
    ("1판 (11-11, 도구 1개 사망)", "TYRYGGSRGKRNNFKTEQY", 0.20),
    ("2판 (11-12, 수리 후)",       "TFEYGGSRGKRNHFKTEEY", 0.06),
]

print("=" * 96)
print(f"§30-4a 검사 — 시드 {SEED_NAME} `{seed}` · off0 = {off0:.4f} · 밴드 ±{OFF_BAND}")
print(f"문헌 보존 위치 {pins} · 기능 하한 rmt_sim ≥ {FLOOR} · 상속 임계 {INHERIT_SIM_MIN}")
print("=" * 96)

for label, seq, recorded in CASES:
    off = off_target(seq)
    d = abs(off - off0)
    rmt = shuttle_similarity(seq).score or 0.0
    edits = [(i + 1, a, b) for i, (a, b) in enumerate(zip(seed, seq)) if a != b]
    pin_ok = all(seq[p - 1] == seed[p - 1] for p in pins)
    print(f"\n── {label}")
    print(f"   서열      {seq}")
    print(f"   편집      {len(edits)}수  {['%d%s>%s' % e for e in edits]}")
    print(f"   off       계산 {off:.4f} · 기록 {recorded:.2f} · |Δ시드| = {d:.4f}")
    print(f"   🔴 밴드    {'통과 (GA 도 생성 가능)' if d <= OFF_BAND + 1e-9 else '🔴 위반 — 현 GA 는 생성 불가'}"
          f"   (허용 {OFF_BAND})")
    print(f"   보존 위치  {'통과' if pin_ok else '🔴 위반'}  (K10={seq[9]} · K15={seq[14]})")
    print(f"   rmt_sim   {rmt:.4f}  → 기능 하한 {'통과' if rmt >= FLOOR else '🔴 위반'}"
          f" · 상속 임계 {'통과' if rmt >= INHERIT_SIM_MIN else '🔴 위반'}")
    legal = d <= OFF_BAND + 1e-9 and pin_ok and rmt >= FLOOR
    print(f"   ⇒ **현 GA 제약 하에서 적법한 설계인가: {'예' if legal else '아니오'}**")
