#!/usr/bin/env /usr/bin/python3
"""축 leave-one-out 애블레이션 — 4축 중 하나를 빼면 전선·front_top·M4 지배가 바뀌나.

`docs/protocol.md` §8 의 "§24-6 축 제거 애블레이션 미동봉" 공백을 메운다.
🔴 `reachable_front(refresh=True)` 를 쓰지 않는다 — 그 함수는 `.cache/reachable/` 에
   결과를 **쓰므로** 패치된 축의 전선이 프로덕션 캐시를 오염시킨다. 게이트·Pareto 를
   같은 헬퍼(`is_candidate`/`is_comparison_only`/`group_rank`)로 여기서 다시 조립한다.
결정론 · LLM 0회.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core.config import DEFAULT_CARGO, LINKER_LIBRARY, SHUTTLES, get_settings
from core.evidence import group_rank
from core.optimizer_agent import (OptimizationAgent, MAX_CANDIDATES_PER_STEP,
                                  is_candidate, is_comparison_only, II_STABLE)

AXES = ["liab", "불안정성(이진)", "off", "응집핫스팟"]

def axis_tuple(r, drop=None):
    ii = float(r.get("instability") or 999.0)
    full = (int(r.get("dev_liab") or 0), 0 if ii < II_STABLE else 1,
            float(r.get("sel_off") or 0.0), float(r.get("dev_agg_hot") or 0.0))
    return tuple(v for i, v in enumerate(full) if i != drop)

def dom(q, r, drop=None):
    a = (group_rank(q.get("tier")),) + axis_tuple(q, drop)
    b = (group_rank(r.get("tier")),) + axis_tuple(r, drop)
    return all(x <= y for x, y in zip(a, b)) and a != b

ag = OptimizationAgent(get_settings())
combos = [{"label": f"{ln}+{sn}", "linker": LINKER_LIBRARY[ln]["seq"],
           "shuttle": SHUTTLES[sn]["seq"]} for ln in LINKER_LIBRARY for sn in SHUTTLES]
rows = []
for i in range(0, len(combos), MAX_CANDIDATES_PER_STEP):
    _, chunk = ag._evaluate(DEFAULT_CARGO, {"constructs": combos[i:i + MAX_CANDIDATES_PER_STEP]})
    rows.extend(chunk)
gate_all = [r for r in rows if not r.get("toxic")]
gate = [r for r in gate_all if is_candidate(r)]
print(f"전수 {len(rows)} · 게이트 통과 {len(gate_all)} · 후보군 {len(gate)}")
assert (len(rows), len(gate_all), len(gate)) == (156, 132, 89), "기준선이 공표치와 다르다 (§52-0)"

def front_of(drop=None):
    f = [r for r in gate if not any(dom(q, r, drop) for q in gate if q is not r)]
    bg = min((group_rank(r.get("tier")) for r in gate), default=None)
    return f, [r for r in f if group_rank(r.get("tier")) == bg]

base_f, base_t = front_of(None)
print(f"\n기준선(4축): 전선 {len(base_f)} · front_top {len(base_t)}")
assert (len(base_f), len(base_t)) == (23, 5), "기준선 전선이 공표치와 다르다"
base_tl = {r["label"] for r in base_t}

print(f"\n{'빼는 축':<16}{'전선':>6}{'front_top':>11}   front_top 변화")
print("-" * 72)
for d, name in enumerate(AXES):
    f, t = front_of(d)
    tl = {r["label"] for r in t}
    added, gone = tl - base_tl, base_tl - tl
    ch = "불변" if not added and not gone else \
         f"+{len(added)} / -{len(gone)}  들어옴 {sorted(added)}  빠짐 {sorted(gone)}"
    print(f"{name:<16}{len(f):>6}{len(t):>11}   {ch}")

# ── M4 설계체의 지배가 어느 축에 얹혀 있나 ──────────────────────────────
DESIGNS = [
    ("A1 §49 EAAAKEAAAK+TFEY…", dict(tier="B", dev_liab=1, instability=33.9, sel_off=0.06, dev_agg_hot=0.55)),
    ("A2 §49 EAAAK+TFEY…",      dict(tier="B", dev_liab=1, instability=37.3, sel_off=0.06, dev_agg_hot=0.55)),
    ("A3 §49 EAAAK+TFEY…",      dict(tier="B", dev_liab=1, instability=37.3, sel_off=0.06, dev_agg_hot=0.55)),
    ("B1 §53 EAAAK+TFGY…",      dict(tier="B", dev_liab=1, instability=34.9, sel_off=0.11, dev_agg_hot=0.59)),
    ("B2 §53 EAAAK+TFEY…",      dict(tier="B", dev_liab=1, instability=37.3, sel_off=0.06, dev_agg_hot=0.55)),
    ("B3 §53 EAAAK+TFEY…",      dict(tier="B", dev_liab=1, instability=37.3, sel_off=0.06, dev_agg_hot=0.55)),
    ("C1 EAAAK+TFEY…",          dict(tier="B", dev_liab=1, instability=37.3, sel_off=0.06, dev_agg_hot=0.55)),
    ("C2 AEAAAKEAAAK+TFEY…",    dict(tier="B", dev_liab=1, instability=36.3, sel_off=0.06, dev_agg_hot=0.55)),
    ("C3 직접융합+ApoE변이체(LAVYQAGARLAVYQYGAR)", dict(tier="B", dev_liab=1, instability=-4.6, sel_off=0.20, dev_agg_hot=0.66)),
]
print(f"\n{'설계체':<40}" + "".join(f"{'4축' if d is None else '−'+AXES[d]:>16}" for d in [None, 0, 1, 2, 3]))
print("-" * 120)
for label, dd in DESIGNS:
    cells = []
    for d in [None, 0, 1, 2, 3]:
        _, t = front_of(d)
        win = sum(1 for b in t if dom(dd, b, d))
        lose = sum(1 for b in t if dom(b, dd, d))
        cells.append(f"{win}/{len(t)}" + (f" (역{lose})" if lose else ""))
    print(f"{label:<40}" + "".join(f"{c:>16}" for c in cells))
