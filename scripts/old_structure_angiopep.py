#!/usr/bin/env python3
"""수정 전 구조에서 같은 입력을 돌린다 — 정답지가 붙은 before/after.

🔴 이 스크립트는 **옛 커밋의 작업 트리에서** 실행한다. 현재 코드에는 여기서 쓰는
   `MIN_BBB_SIGNAL` 도 BBB 순위 축도 없다(폐지됨).

실행 방법 (유료 LLM 0회):
    git worktree add --detach /tmp/bbb-old d671c60
    ln -s "$PWD/vendor" /tmp/bbb-old/vendor          # vendor/ 는 .gitignore 대상
    cd /tmp/bbb-old && PYTHONPATH=. /usr/bin/python3 \
        /path/to/scripts/old_structure_angiopep.py

기준 커밋 — `d671c60` (= `9b43eaf` 의 부모)
    9b43eaf "셔틀 순위를 점수에서 문헌 인용으로 이전 + BBB 하한 폐지" **직전** 상태.
    따라서 이 트리에는 ⑴ BBB 하한 하드 게이트 ⑵ BBB 를 1번 축으로 쓰는 Pareto
    ⑶ 문헌 근거 등급 없음 ⑷ 보존 위치 동결 없음 — 즉 5절이 "뺐다"고 말하는 것들이
    **전부 살아 있다**.

왜 이 비교가 성립하는가 (§30 의 단위 불일치를 피하는 근거)
--------------------------------------------------------------------------
비교하는 것은 **수치가 아니라 판정**이다. 그리고 판정을 가르는 입력인 deepB3P 는
`core/predictors.py` 가 **전 기간 커밋 2개**(최초 + 주석 수정 1건)뿐이라 채점이
바뀐 적이 없다. liability 처럼 정의가 바뀐 양은 이 스크립트에서 **쓰지 않는다**.

🔴 주장 범위: 이것은 "옛 시스템이 Angiopep-7 을 최종 출력했을 것"이 아니다.
   옛 시스템은 셔틀 12종 전수를 돌렸다. 여기서 보이는 것은 **정답이 붙은 한 계열
   안에서 옛 규칙이 매긴 순서**뿐이다.
"""
from core.config import (DEFAULT_CARGO, LINKER_LIBRARY, MIN_BBB_SIGNAL,
                         OFF_TARGET_PENALTY, SHUTTLES, get_settings)
from core.optimizer_agent import OptimizationAgent

A2 = SHUTTLES["Angiopep-2"]["seq"]
assert A2[9] == "K" and A2[14] == "K"
A5 = A2[:9] + "R" + A2[10:]
A7 = A5[:14] + "R" + A5[15:]
PANEL = [("Angiopep-2", A2), ("Angiopep-5", A5), ("Angiopep-7", A7)]
CARGO = DEFAULT_CARGO
II_STABLE = 40.0

ag = OptimizationAgent(get_settings())
cons = [{"label": f"{n}+{lk}", "linker": lv["seq"], "shuttle": s}
        for n, s in PANEL for lk, lv in LINKER_LIBRARY.items()]
from core.optimizer_agent import MAX_CANDIDATES_PER_STEP as _MAX  # 옛 코드와 같은 분할
rows = []
for i in range(0, len(cons), _MAX):                          # optimizer_agent 와 동일
    _, chunk = ag._evaluate(CARGO, {"constructs": cons[i:i + _MAX]})
    rows.extend(chunk)
assert len(rows) == len(cons), f"평가 누락: {len(rows)}/{len(cons)}"
for r in rows:
    r["_fam"] = r["label"].split("+")[0]

def axes(r):
    ii = float(r.get("instability") or 999.0)
    return (r.get("bbb") or 0.0, r.get("sel_off") or 0.0, 0 if ii < II_STABLE else 1)

def dominated(r, pool):
    b1, o1, i1 = axes(r)
    for q in pool:
        if q is r: continue
        b2, o2, i2 = axes(q)
        if b2 >= b1 and o2 <= o1 and i2 <= i1 and (b2 > b1 or o2 < o1 or i2 < i1):
            return True
    return False

print("=" * 78)
print("수정 전 구조 (커밋 d671c60 · BBB 가 순위 1번 축) — Angiopep 3종 × 링커 13종")
print("=" * 78)
print(f"하드 게이트: 독성 + BBB 하한 {MIN_BBB_SIGNAL}\n")

print("[1] 계열별 요약 — 셔틀 단독 BBB 와 구성체 BBB 최대")
print("-" * 78)
print(f"{'계열':<12}{'셔틀단독BBB':>12}{'구성체BBB최대':>14}{'하한통과':>10}{'최고 링커':>16}")
for n, s in PANEL:
    sub = [r for r in rows if r["_fam"] == n]
    best = max(sub, key=lambda r: r.get("bbb") or 0.0)
    sb = ag.predictor.predict(s).bbb_permeability
    npass = sum(1 for r in sub if (r.get("bbb") or 0) >= MIN_BBB_SIGNAL)
    print(f"{n:<12}{sb:>12.3f}{best['bbb']:>14.3f}{npass:>8}/13{best['label'].split('+')[1]:>16}")

nontoxic = [r for r in rows if not r.get("toxic")]
gate = [r for r in nontoxic if (r.get("bbb") or 0.0) >= MIN_BBB_SIGNAL]
print()
print("[2] 옛 하드 게이트의 판정")
print("-" * 78)
print(f"  독성 통과            : {len(nontoxic)}/39")
print(f"  BBB 하한({MIN_BBB_SIGNAL}) 통과 : **{len(gate)}/39**")
if not gate:
    print("  🔴 문헌으로 검증된 Angiopep-2 조합 13개가 **전부 탈락**한다.")
    print("     즉 옛 하드 게이트는 이 계열을 **전멸**시킨다.")

print()
print("[3] 게이트를 무시하고 옛 1번 축(BBB)으로만 줄을 세우면")
print("-" * 78)
ranked = sorted(rows, key=lambda r: -(r.get("bbb") or 0.0))
for i, r in enumerate(ranked[:6], 1):
    print(f"  {i}. {r['label']:<28} BBB {r['bbb']:.3f}  off {r.get('sel_off') or 0:.2f}")
print("   ...")
worst = [r for r in ranked if r["_fam"] == "Angiopep-2"]
print(f"  Angiopep-2 계열 최상위 순위: {ranked.index(worst[0]) + 1}위 / 39")
fam_best = {n: max((r.get('bbb') or 0) for r in rows if r['_fam'] == n) for n, _ in PANEL}
order = sorted(fam_best, key=lambda k: -fam_best[k])
print(f"  🔴 계열 순서: **{' > '.join(order)}**")
print(f"     문헌 순서: **Angiopep-2 > Angiopep-5 > Angiopep-7** (−기준 / −60% / −85%)")
print(f"     → **정확히 거꾸로다.**")

front = [r for r in rows if not dominated(r, rows)]
print()
print("[4] 게이트를 무시하고 옛 3축 Pareto 전선을 구하면")
print("-" * 78)
for r in sorted(front, key=lambda r: -(r.get("bbb") or 0)):
    print(f"  {r['label']:<28} BBB {r['bbb']:.3f}  off {r.get('sel_off') or 0:.2f}  "
          f"II {r.get('instability') or 0:.1f}")
fams = sorted({r['_fam'] for r in front})
print(f"  전선 {len(front)}개 · 포함된 계열: {fams}")
if "Angiopep-2" not in fams:
    print("  🔴 **문헌이 유일하게 작동한다고 말한 Angiopep-2 가 전선에 없다.**")
