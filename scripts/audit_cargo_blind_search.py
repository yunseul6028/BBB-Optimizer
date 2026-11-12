#!/usr/bin/env python3
"""§35-5 후속 — **탐색이 화물을 안 보는 것에 대가가 있나**.

질문
----
§35-5 는 맹점의 **존재**를 확인했다: `pareto_scrub` 내부의 `objectives()`·
`off_target()`·`shuttle_similarity()`·`_frozen_idx()` 는 전부 셔틀 서열만 받는데,
심사의 `dev_liab`·`instability` 는 **구성체 전체**(화물+링커+셔틀)로 계산된다.
이 스크립트는 존재를 다시 증명하지 않는다. **결과를 잰다**:

  탐색이 구성체 수준 목적을 최적화했다면 **다른 서열을 골랐을까?**

재는 것 (전부 결정론·재실행 가능)
-------------------------------
A. **전선 불일치(사후 재채점)** — 프로덕션 아카이브의 각 멤버에 두 벡터를 매긴다.
   비지배 집합이 달라지는가. **집합 차이 자체**를 적는다(크기만이 아니라 서열).
B. **대표가 뒤집히나** — ⑴ 탐색 자신의 사전식 1위, ⑵ **실제 심사**
   (`optimizer_agent.rank_key` = 근거그룹 + 4축, 그리고 `dominates`).
C. **반사실 탐색** — `o.objectives` 를 구성체 수준으로 몽키패치하고 **같은 시드로
   다시 돌린다**. 사후 재채점은 셔틀 전용 탐색이 방문한 후보만 보지만, 패치된 탐색은
   궤적 자체가 달라진다. 출력 서열 집합을 비교한다.
C-2. **대응 널** — 화물을 **조성 보존 셔플**로 바꿔 같은 반사실 탐색을 돌린다.
   목적을 바꾸면 궤적이 달라지는 것만으로도 출력 집합이 흔들리므로, "많이 다르다"가
   **화물 정보** 때문인지 **궤적 요동** 때문인지 가른다.
D. **축 분해** — 한 축씩만 구성체 범위로 바꿔 전선을 다시 만든다. 어느 축이
   불일치를 만드는가. D-2 는 liability Δ 를 **항목별**로 쪼갠다 — 지배 관계는
   상수 더하기에 불변이므로, 멤버마다 값이 **달라지는** 항목만 순위를 움직인다.
E. **화물 의존성** — cargo A(§35 의 Aβ(25-35)) · cargo B(§35 의 DSIP `WAGGDASGE`)
   둘 다. 불일치가 두 화물에서 **동일하면 그렇게 적는다** — 그건 불일치마저
   화물에 눈을 감았다는 뜻이다.

두 벡터의 정의 — 🔴 **범위만 다르고 양(quantity)은 같다**
-------------------------------------------------------
| 축 | V_search (탐색이 실제로 쓰는 것) | V_construct (심사 범위로 넓힌 것) |
|---|---|---|
| liability | `n_liabilities(셔틀)` | `n_liabilities(화물+링커+셔틀)` ← 심사 `dev_liab` 범위 |
| 응집 | `agg_hotspot_a3v(셔틀)` **원값** | `agg_hotspot_a3v(링커+셔틀)` **원값** ← 심사 `dev_agg_hot` 범위 |
| 불안정성 | `instability_index(셔틀)` **연속** | `instability_index(화물+링커+셔틀)` **연속** |

응집은 심사가 **정규화값**(`agg_hotspot`, 0~1·소수 2자리)을 쓰고 탐색은 **원값**
(`agg_hotspot_a3v`)을 쓴다(§36, 의도된 불일치). 여기서는 **범위 효과만 분리**하려고
두 벡터 모두 **원값**으로 통일했다 — 정규화 반올림이 만드는 동점은 별도 문제다.
심사 1위 계산(B-⑵)에서는 반대로 **심사의 양을 글자 그대로** 쓴다(정규화 응집·이진 불안정성).

off-target 은 어느 벡터에도 없다 — 목적이 아니라 **밴드 제약**이고(§26) 두 벡터
모두에서 **셔틀 단독**이다. 즉 구조적으로 화물 불변이며 이 측정으로 움직일 수 없다.

이 측정이 **볼 수 없는 것**
---------------------------
⒜ **독성 하드 게이트를 돌리지 않았다.** §35-3 이 "가장 화물 의존적인 장치"라고
   밝힌 것이 바로 독성 게이트(뒤집힘 24/156)인데, 여기서는 순위 축만 본다.
   따라서 "대표가 안 바뀐다"는 결론은 **게이트 통과를 가정한 조건부**다.
⒝ **BBB·ESMFold 를 돌리지 않았다** — 순위 축이 아니다(M6·M8).
⒞ **링커는 고정 2종**(직접융합 `""` · `GGGGS`)이다. 링커·셔틀 공진화는 범위 밖이다.
⒟ 반사실 탐색은 `objectives` 만 구성체 수준으로 바꾼다. `rmt_sim` 하한·동결 위치·
   off 밴드는 **그대로 셔틀 전용**이다(그게 실제 코드의 모습이다).
⒠ **실측이 아니다.** 네 축 전부 규칙/조성 기반 계산이며 실험 검증이 아니다.
⒡ 화물 n = 2. §35-6 ⓐ 와 같은 한계다.

실행: `/usr/bin/python3 scripts/audit_cargo_blind_search.py`   (⚠️ `python3 -I` 금지)
`core/` 는 건드리지 않는다 — 변주는 모듈 속성 몽키패치로만 하고 `finally` 에서 복원,
`main()` 끝에서 프로덕션 기본값과 같은지 **단언**한다.
"""

from __future__ import annotations

import sys
import unicodedata
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import core.developability_opt as o                      # noqa: E402
from core.binding import shuttle_similarity              # noqa: E402
from core.config import SHUTTLES                         # noqa: E402
from core.developability import (                        # noqa: E402
    _CHARGE_LIAB_LABELS,
    assess_developability,
    junction_liabilities,
)
from core.evidence import group_rank, inherit             # noqa: E402
from core.optimizer_agent import dominates, rank_key, risk_axes  # noqa: E402
from core.selectivity import assess_selectivity           # noqa: E402
from core.stability import assess_stability               # noqa: E402

# ---------------------------------------------------------------------------
# 측정 조건 — 모든 난수 출처를 여기 적는다
# ---------------------------------------------------------------------------
#: `pareto_scrub` 시드. 0~4 는 §25-2·§26-5·§29·§35 와 같은 집합이고 그 위로 확장했다.
SEEDS: tuple[int, ...] = tuple(range(20))
#: 셔틀 — `rmt_sim(자기) = 1.0` 이라 시드 밀반입 차단(§33-4)에 걸리지 않는 후보군 RMT.
#: Tet1 은 `REFERENCE_SHUTTLES` 에 없어 `sim0 = 0.0` → `pareto_scrub` 가 예외를 던진다(§32).
TARGETS: tuple[str, ...] = ("Angiopep-2", "ApoE(159-167)2", "T7", "B6",
                            "RVG29", "CDX", "THR", "Leptin30")
#: 화물 — §35-1 과 **같은 쌍**이라 수치를 나란히 둘 수 있다.
CARGO_A = "GSNKGAIIGLM"        # Aβ(25-35) · 11aa · 순전하 +1 · GRAVY +0.65
CARGO_B = "WAGGDASGE"          # DSIP · 9aa · 순전하 −2 · GRAVY −0.70
CARGOS = (("A", CARGO_A), ("B", CARGO_B))
#: 링커 — 직접융합(§35-2 대표의 링커)과 표준 GS 링커. `""` 에서는 설계모듈 = 셔틀이므로
#: 응집 축의 범위 차이가 **정의상 0** 이 된다(그 사실 자체가 결과다).
LINKERS = (("직접융합", ""), ("GGGGS", "GGGGS"))
#: 🔴 **대응 널(화물 정체성 대조군)** — cargo A 의 **조성·길이를 정확히 보존한 셔플**.
#: 왜 필요한가: 목적함수를 바꾸면 아카이브가 바뀌고, 아카이브가 바뀌면 부모 표본이
#: 바뀌어 GA 궤적 전체가 달라진다. 그래서 "출력 서열 집합이 많이 다르다"만으로는
#: **화물 정보 때문인지 궤적 요동(churn) 때문인지** 구분할 수 없다. 조성 보존 셔플은
#: 전하·GRAVY·길이·liability 조성을 그대로 두고 **순서만** 바꾼 화물이다(= `scramble_delta`
#: 계열의 널). `Jaccard(실제화물, 셔플화물)` 이 `Jaccard(프로덕션, 실제화물)` 만큼 낮으면
#: 차이의 대부분은 **화물 정체성이 아니라 궤적 요동**이다.
#: 셔플은 측정 전에 고정했다: `random.Random(20261112).sample(CARGO_A, len(CARGO_A))`.
CARGO_A_SCRAMBLE = "KGMIIASGLGN"   # = GSNKGAIIGLM 의 조성 보존 셔플(위 시드로 재생성 가능)
RMT_FLOOR = 0.85
ROUNDS, POP = 40, 30           # §26-5·§29·§35 와 동일

_PROD_OBJECTIVES = o.objectives     # 복원용 원본 참조


# ---------------------------------------------------------------------------
# 표시 보조 (한글·이모지는 2칸)
# ---------------------------------------------------------------------------
def _w(s: str, width: int, right: bool = False) -> str:
    disp = sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in str(s))
    pad = " " * max(0, width - disp)
    return (pad + str(s)) if right else (str(s) + pad)


# ---------------------------------------------------------------------------
# 두 목적 벡터 — 범위만 다르다
# ---------------------------------------------------------------------------
def v_search(seq: str) -> tuple:
    """탐색이 **실제로** 최소화하는 벡터 = `core.developability_opt.objectives` 그대로."""
    return _PROD_OBJECTIVES(seq)


def v_construct(seq: str, cargo: str, linker: str) -> tuple:
    """같은 세 양을 **심사의 범위**로 계산한 벡터 (liability·불안정성 ⟵ 화물 포함)."""
    full = cargo + linker + seq
    module = linker + seq
    return (assess_developability(full).n_liabilities,
            round(assess_developability(module).agg_hotspot_a3v, 3),
            o._instability(full))


def v_axis_swap(seq: str, cargo: str, linker: str, axis: int) -> tuple:
    """**한 축만** 구성체 범위로 바꾼 벡터 (axis 0=liability · 1=응집 · 2=불안정성)."""
    s = list(v_search(seq))
    c = v_construct(seq, cargo, linker)
    s[axis] = c[axis]
    return tuple(s)


def nondominated(members: list, vec) -> list:
    """`vec(seq)` 기준 비지배 집합 — 지배 판정은 `o._dominates`(탐색이 쓰는 그것)."""
    vals = {m: vec(m) for m in members}
    return [m for m in members
            if not any(o._dominates(vals[q], vals[m]) for q in members if q != m)]


# ---------------------------------------------------------------------------
# 심사 쪽 — `optimizer_agent` 의 실제 행(row)을 재현한다
# ---------------------------------------------------------------------------
def judge_row(seq: str, cargo: str, linker: str) -> dict:
    """`OptimizationAgent._evaluate` 가 만드는 행의 **순위 관련 키만** 같은 호출로 채운다.

    🔴 키 이름을 `_evaluate` 와 **같게** 둔다 — `risk_axes`·`dominates`·`rank_key` 를
    어댑터 없이 그대로 먹여야 "실제 비교"이기 때문이다.
    ⚠️ 독성(ToxinPred3)·BBB(deepB3P)는 **부르지 않는다** — 순위 축이 아니다. 독성은
    순위 축이 아니라 **하드 게이트**이고 §35-3 에서 가장 화물 의존적인 장치였다.
    그래서 아래 1위 판정은 **"독성 게이트를 통과했다고 가정할 때"** 의 조건부다.
    """
    full = cargo + linker + seq
    module = linker + seq
    dev_full = assess_developability(full)
    dev_mod = assess_developability(module) if module else dev_full
    tier, parent, tsim = inherit(seq)
    return {"shuttle": seq, "tier": tier, "tier_parent": parent, "tier_sim": tsim,
            "dev_liab": dev_full.n_liabilities,
            "instability": assess_stability(full).instability_index,
            "sel_off": assess_selectivity(seq).off_target_risk,
            "dev_agg_hot": dev_mod.agg_hotspot,
            "mechanism": shuttle_similarity(seq).mechanism}


def judge_best(members: list, cargo: str, linker: str) -> tuple:
    """심사의 결정론 1위 = `rank_key` 최소. 반환 (서열, row, rank_key)."""
    rows = [judge_row(m, cargo, linker) for m in members]
    best = min(rows, key=rank_key)
    return best["shuttle"], best, rank_key(best)


# ---------------------------------------------------------------------------
# 반사실 탐색 — `o.objectives` 만 구성체 수준으로 바꿔 같은 시드로 재실행
# ---------------------------------------------------------------------------
class PatchCounter:
    """🔴 패치가 **실제로 불렸는지** 세는 계수기. 0 이면 '차이 없음' 결론은 무효다.

    왜 필요한가: `developability_opt` 는 모듈 수준에서 심볼을 당겨 쓰므로
    (`from .developability import assess_developability`, opt.py:22) 한쪽만 바꾸면
    **조용히 무효화**된다. 이 프로젝트에서 그 함정이 가짜 '동일' 결론을 만든 적이 있다.
    여기서는 `objectives` 를 패치하므로 `pareto_scrub` 안의 전역 조회가 잡히지만,
    믿지 않고 센다.
    """

    def __init__(self, cargo: str, linker: str):
        self.cargo, self.linker, self.n = cargo, linker, 0

    def __call__(self, seq: str) -> tuple:
        self.n += 1
        return v_construct(seq, self.cargo, self.linker)


def scrub(shuttle: str, seed: int) -> list:
    """프로덕션 `pareto_scrub` 한 번 → 서열 리스트."""
    return [m["seq"] for m in o.pareto_scrub(shuttle, rounds=ROUNDS, pop=POP,
                                             rmt_floor=RMT_FLOOR, seed=seed)]


def scrub_cargo_aware(shuttle: str, seed: int, cargo: str, linker: str) -> tuple:
    """`objectives` 만 구성체 수준으로 패치한 `pareto_scrub`. 반환 (서열들, 호출 횟수)."""
    ctr = PatchCounter(cargo, linker)
    try:
        o.objectives = ctr
        seqs = scrub(shuttle, seed)
    finally:
        o.objectives = _PROD_OBJECTIVES
    return seqs, ctr.n


# ---------------------------------------------------------------------------
# 본 측정
# ---------------------------------------------------------------------------
def main() -> None:
    assert o.objectives is _PROD_OBJECTIVES
    print("§35-5 후속 — 탐색의 화물 맹점에 **대가가 있나**")
    print(f"시드 {SEEDS[0]}~{SEEDS[-1]} ({len(SEEDS)}개) · 셔틀 {len(TARGETS)}종 · "
          f"화물 2종 · 링커 {len(LINKERS)}종 · rounds={ROUNDS} pop={POP} "
          f"rmt_floor={RMT_FLOOR} off 밴드 ±{o.OFF_BAND}")
    print("⚠️ 독성 하드 게이트·BBB 는 돌리지 않는다(순위 축이 아니다) → 1위 판정은 조건부.\n")

    # ── 0) 프로덕션 아카이브 (화물·링커와 무관하게 **한 번만** 돌면 된다) ──
    prod: dict = {}
    for name in TARGETS:
        seq0 = SHUTTLES[name]["seq"]
        for sd in SEEDS:
            prod[(name, sd)] = scrub(seq0, sd)
    sizes = [len(v) for v in prod.values()]
    print(f"프로덕션 아카이브 {len(prod)}개 · 멤버 총 {sum(sizes)} · "
          f"평균 {sum(sizes)/len(sizes):.2f} · 범위 {min(sizes)}~{max(sizes)}")

    # 하네스 검증 ①: 아카이브가 V_search 기준으로 정말 비지배 집합인가
    bad = [(k, len(v), len(nondominated(v, v_search))) for k, v in prod.items()
           if len(nondominated(v, v_search)) != len(v)]
    print(f"하네스 검증 — V_search 비지배 아님인 아카이브: {len(bad)}개"
          + (f" (예: {bad[:3]})" if bad else " 🟢"))

    # 하네스 검증 ②: 두 벡터가 실제로 다른 숫자를 돌려주는가(패치 전 sanity)
    probe = SHUTTLES["Angiopep-2"]["seq"]
    print(f"하네스 검증 — Angiopep-2: V_search={v_search(probe)} · "
          f"V_construct(A,직접융합)={v_construct(probe, CARGO_A, '')} · "
          f"V_construct(B,직접융합)={v_construct(probe, CARGO_B, '')}\n")

    # ── A) 전선 불일치 (사후 재채점) + B-⑴ 탐색 자신의 사전식 1위 ──
    print("=" * 104)
    print("표 A — **사후 재채점**: 같은 아카이브를 두 범위로 재채점했을 때 비지배 집합의 차이")
    print("  '탈락' = V_search 전선에 있었지만 V_construct 에서 지배당한 멤버")
    print("  '추가' = 그 반대. 사전식1위 = 각 벡터의 튜플 최소(탐색 자신의 tie-break)")
    print(_w("셔틀", 16) + _w("화물", 6) + _w("링커", 10) + _w("멤버", 6, True)
          + _w("탈락", 7, True) + _w("추가", 7, True) + _w("불일치시드", 12, True)
          + _w("사전1위변경", 13, True))
    detail: dict = {}
    for name in TARGETS:
        for cname, cargo in CARGOS:
            for lname, linker in LINKERS:
                n_mem = n_drop = n_add = n_seed_diff = n_lexflip = 0
                drops, lexflips = [], []
                for sd in SEEDS:
                    mem = prod[(name, sd)]
                    n_mem += len(mem)
                    fs = set(nondominated(mem, v_search))
                    fc = set(nondominated(mem, lambda s: v_construct(s, cargo, linker)))
                    d, a = fs - fc, fc - fs
                    n_drop += len(d)
                    n_add += len(a)
                    if d or a:
                        n_seed_diff += 1
                        drops.append((sd, sorted(d), sorted(a)))
                    lx_s = min(mem, key=v_search)
                    lx_c = min(mem, key=lambda s: v_construct(s, cargo, linker))
                    if lx_s != lx_c:
                        n_lexflip += 1
                        lexflips.append((sd, lx_s, lx_c))
                detail[(name, cname, lname)] = {"drops": drops, "lexflips": lexflips}
                print(_w(name, 16) + _w(cname, 6) + _w(lname, 10) + _w(n_mem, 6, True)
                      + _w(n_drop, 7, True) + _w(n_add, 7, True)
                      + _w(f"{n_seed_diff}/{len(SEEDS)}", 12, True)
                      + _w(f"{n_lexflip}/{len(SEEDS)}", 13, True))

    print("\n집합 차이의 **내용**(앞 12건) — '→' 왼쪽이 탈락, 오른쪽이 추가")
    shown = 0
    for key, v in detail.items():
        for sd, d, a in v["drops"]:
            if shown >= 12:
                break
            print(f"  {key[0]}·cargo{key[1]}·{key[2]}·seed{sd}: "
                  f"{[x for x in d]} → {[x for x in a]}")
            shown += 1
    if shown == 0:
        print("  (없음 — 사후 재채점으로는 비지배 집합이 한 멤버도 달라지지 않았다)")

    print("\n사전식 1위가 바뀐 사례(앞 12건)")
    shown = 0
    for key, v in detail.items():
        for sd, a, b in v["lexflips"]:
            if shown >= 12:
                break
            print(f"  {key[0]}·cargo{key[1]}·{key[2]}·seed{sd}: {a} → {b}")
            shown += 1
    if shown == 0:
        print("  (없음)")

    # ── D) 축 분해 — 한 축씩만 구성체 범위로 ──
    print("\n" + "=" * 104)
    print("표 D — **축 분해**: 한 축만 구성체 범위로 바꿨을 때. "
          "Δ≠0 = 그 축의 값이 범위 때문에 달라진 멤버 수 / 전체 멤버 수")
    print(_w("화물", 6) + _w("링커", 10) + _w("축", 14) + _w("Δ≠0", 12, True)
          + _w("평균|Δ|", 10, True) + _w("전선탈락", 10, True) + _w("사전1위변경", 13, True))
    axis_names = ("liability", "응집(원값)", "불안정성")
    for cname, cargo in CARGOS:
        for lname, linker in LINKERS:
            for ax in range(3):
                n_mem = n_ne = n_drop = n_lex = 0
                sum_abs = 0.0
                for (name, sd), mem in prod.items():
                    for m in mem:
                        n_mem += 1
                        a, b = v_search(m)[ax], v_construct(m, cargo, linker)[ax]
                        if a != b:
                            n_ne += 1
                        sum_abs += abs(b - a)
                    fs = set(nondominated(mem, v_search))
                    fx = set(nondominated(mem, lambda s: v_axis_swap(s, cargo, linker, ax)))
                    n_drop += len(fs - fx)
                    if min(mem, key=v_search) != min(
                            mem, key=lambda s: v_axis_swap(s, cargo, linker, ax)):
                        n_lex += 1
                print(_w(cname, 6) + _w(lname, 10) + _w(axis_names[ax], 14)
                      + _w(f"{n_ne}/{n_mem}", 12, True)
                      + _w(f"{sum_abs/n_mem:.3f}", 10, True)
                      + _w(n_drop, 10, True)
                      + _w(f"{n_lex}/{len(prod)}", 13, True))
    print("  off-target 은 이 표에 **없다** — 목적이 아니라 밴드 제약이고 두 범위 모두 "
          "셔틀 단독이다(§35-3 의 `sel_off` 0/156 과 같은 이유: 정의상 불변).")

    # ── D-2) liability Δ 를 **항목별**로 쪼갠다 ──
    # 왜: 지배 관계는 **상수 더하기에 불변**이다. 모든 후보에 같은 수를 더하는 항은
    # 순위를 못 바꾼다. 그래서 "liability 축이 화물을 본다"가 실제로 순위를 움직이려면
    # Δ 가 **후보마다 달라야** 한다 — 어느 항목이 그런지 본다.
    # 범위(min~max)가 한 점인 항목은 **순수 상수**이고 지배 관계에 무력하다.
    print("\n표 D-2 — liability 범위 차이의 **항목별** 분해. "
          "Δ = (화물+링커+셔틀의 항목 수) − (셔틀의 항목 수)")
    print("  🔴 범위가 한 점 = 순수 상수 = **지배 관계를 바꿀 수 없다**. "
          "범위가 벌어진 항목만 순위를 움직인다")
    print(_w("화물", 6) + _w("링커", 10) + _w("항목", 22) + _w("평균Δ", 9, True)
          + _w("Δ≠0", 11, True) + _w("범위", 12, True))

    def _cat(x: str) -> str:
        for k in ("응집 핫스팟", "자유 시스테인", "산화 취약 잔기"):
            if x.startswith(k):
                return k
        return "모티프:" + x.split(" ")[0]

    for cname, cargo in CARGOS:
        for lname, linker in LINKERS:
            per: dict = {}
            jn_any = 0
            for mem in prod.values():
                for m in mem:
                    full = assess_developability(cargo + linker + m)
                    sh = assess_developability(m)
                    cf_, cs = {}, {}
                    for x in full.liabilities:
                        if x not in full.charge_liabilities:
                            cf_[_cat(x)] = cf_.get(_cat(x), 0) + 1
                    for x in sh.liabilities:
                        if x not in sh.charge_liabilities:
                            cs[_cat(x)] = cs.get(_cat(x), 0) + 1
                    for k in set(cf_) | set(cs):
                        per.setdefault(k, []).append(cf_.get(k, 0) - cs.get(k, 0))
                    jn_any += any(x["label"] not in _CHARGE_LIAB_LABELS
                                  for x in junction_liabilities(cargo, linker, m))
            n_tot = sum(len(v) for v in prod.values())
            for k in sorted(per, key=lambda z: -abs(sum(per[z]))):
                v = per[k] + [0] * (n_tot - len(per[k]))
                rng = f"{min(v)}~{max(v)}" + (" 🔴상수" if min(v) == max(v) else "")
                print(_w(cname, 6) + _w(lname, 10) + _w(k, 22)
                      + _w(f"{sum(v)/n_tot:+.3f}", 9, True)
                      + _w(f"{sum(1 for x in v if x)}/{n_tot}", 11, True)
                      + _w(rng, 12, True))
            print(_w("", 16) + _w("접합부 모티프 보유 멤버", 22)
                  + _w(f"{jn_any}/{n_tot}", 21, True))

    # ── C) 반사실 탐색 — objectives 를 구성체 수준으로 패치해 재실행 ──
    print("\n" + "=" * 104)
    print("표 C — **반사실 탐색**: `objectives` 를 구성체 수준으로 바꿔 **같은 시드로 재실행**.")
    print("  Jaccard = |교집합|/|합집합| (서열 집합) · '신규' = 프로덕션이 못 찾은 서열")
    print(_w("셔틀", 16) + _w("화물", 6) + _w("링커", 10) + _w("프로덕션", 10, True)
          + _w("반사실", 8, True) + _w("공통", 7, True) + _w("신규", 7, True)
          + _w("Jaccard", 9, True) + _w("패치호출", 10, True))
    cf: dict = {}
    total_calls = 0
    for name in TARGETS:
        seq0 = SHUTTLES[name]["seq"]
        for cname, cargo in CARGOS:
            for lname, linker in LINKERS:
                np_, nc, ni, nn, calls = 0, 0, 0, 0, 0
                for sd in SEEDS:
                    pset = set(prod[(name, sd)])
                    cseqs, k = scrub_cargo_aware(seq0, sd, cargo, linker)
                    cf[(name, cname, lname, sd)] = cseqs
                    cset = set(cseqs)
                    np_ += len(pset)
                    nc += len(cset)
                    ni += len(pset & cset)
                    nn += len(cset - pset)
                    calls += k
                total_calls += calls
                union = np_ + nc - ni
                print(_w(name, 16) + _w(cname, 6) + _w(lname, 10) + _w(np_, 10, True)
                      + _w(nc, 8, True) + _w(ni, 7, True) + _w(nn, 7, True)
                      + _w(f"{ni/union:.3f}" if union else "n/a", 9, True)
                      + _w(calls, 10, True))
    print(f"  🔴 패치 호출 총합 {total_calls} — 0 이면 위 결과는 전부 무효다.")
    assert total_calls > 0, "패치가 한 번도 불리지 않았다 → 측정 무효"
    assert o.objectives is _PROD_OBJECTIVES, "objectives 복원 실패"

    # ── C-2) 대응 널: 화물 **정체성**이 한 일 vs 궤적 요동이 한 일 ──
    # 직접융합만 돌린다(비용·§35-2 대표의 링커). 셔플 화물은 조성·길이·전하·GRAVY 보존.
    print("\n표 C-2 — **대응 널(조성 보존 화물 셔플)** · 링커 = 직접융합만")
    print(f"  cargo A = {CARGO_A} · 셔플 = {CARGO_A_SCRAMBLE} (같은 조성·같은 길이)")
    print("  J(프로덕션,A) 가 낮은 것은 **목적을 바꿨다는 사실** 때문일 수도 있다. "
          "J(A,셔플) 이 그만큼 낮으면 차이의 대부분은 **화물 정체성이 아니라 궤적 요동**이다.")
    print(_w("셔틀", 16) + _w("J(프로덕션,A)", 15, True) + _w("J(프로덕션,셔플)", 17, True)
          + _w("J(A,셔플)", 12, True) + _w("A·셔플 집합 동일", 17, True))
    for name in TARGETS:
        seq0 = SHUTTLES[name]["seq"]
        ip = isc = ia = up = usc = ua = same = 0
        for sd in SEEDS:
            pset = set(prod[(name, sd)])
            aset = set(cf[(name, "A", "직접융합", sd)])
            sset, _k = scrub_cargo_aware(seq0, sd, CARGO_A_SCRAMBLE, "")
            sset = set(sset)
            ip += len(pset & aset); up += len(pset | aset)
            isc += len(pset & sset); usc += len(pset | sset)
            ia += len(aset & sset); ua += len(aset | sset)
            same += (aset == sset)
        print(_w(name, 16) + _w(f"{ip/up:.3f}", 15, True) + _w(f"{isc/usc:.3f}", 17, True)
              + _w(f"{ia/ua:.3f}", 12, True) + _w(f"{same}/{len(SEEDS)}", 17, True))
    assert o.objectives is _PROD_OBJECTIVES, "objectives 복원 실패"

    # ── B-⑵ 실제 심사(rank_key·dominates)로 대표 비교 ──
    print("\n" + "=" * 104)
    print("표 B — **실제 심사**(`rank_key` = 근거그룹 + 4축)로 본 대표. "
          "비교 대상은 ⑴ 프로덕션 아카이브, ⑵ 반사실 아카이브.")
    print("  '심사1위 변경' = 두 아카이브의 `rank_key` 최소 서열이 다른 시드 수")
    print("  '지배' = 반사실 1위가 프로덕션 1위를 `dominates()` 하는 시드 수(그 역도)")
    print(_w("셔틀", 16) + _w("화물", 6) + _w("링커", 10) + _w("심사1위변경", 13, True)
          + _w("반사실▷프로덕션", 16, True) + _w("프로덕션▷반사실", 16, True)
          + _w("비교불가", 10, True))
    trades: list = []
    regress: list = []
    gated_top = 0            # 심사 1위가 실제 게이트(CPP·근거없음)에 걸릴 사례 수
    for name in TARGETS:
        for cname, cargo in CARGOS:
            for lname, linker in LINKERS:
                flip = dom_c = dom_p = incomp = 0
                for sd in SEEDS:
                    pm, cm = prod[(name, sd)], cf[(name, cname, lname, sd)]
                    ps, prow, pk = judge_best(pm, cargo, linker)
                    cs, crow, ck = judge_best(cm, cargo, linker)
                    gated_top += sum(1 for r in (prow, crow)
                                     if group_rank(r.get("tier")) > 1
                                     or str(r.get("mechanism") or "").upper()
                                     .startswith("CPP"))
                    if ps != cs:
                        flip += 1
                        if dominates(crow, prow):
                            dom_c += 1
                        elif dominates(prow, crow):
                            dom_p += 1
                            if len(regress) < 8:
                                regress.append((name, cname, lname, sd, prow, crow))
                        else:
                            incomp += 1
                            if len(trades) < 10:
                                trades.append((name, cname, lname, sd, prow, crow))
                print(_w(name, 16) + _w(cname, 6) + _w(lname, 10)
                      + _w(f"{flip}/{len(SEEDS)}", 13, True)
                      + _w(dom_c, 16, True) + _w(dom_p, 16, True)
                      + _w(incomp, 10, True))

    if trades:
        print("\n**교환관계**(어느 쪽도 지배하지 않음) 사례 — 개선 축과 악화 축을 같은 줄에")
        for name, cname, lname, sd, prow, crow in trades:
            a, b = risk_axes(prow), risk_axes(crow)
            names = ("dev_liab", "불안정성(이진)", "sel_off", "dev_agg_hot")
            moves = [f"{names[i]} {a[i]}→{b[i]}" for i in range(4) if a[i] != b[i]]  # noqa: E501
            print(f"  {name}·cargo{cname}·{lname}·seed{sd}: {prow['shuttle']} → "
                  f"{crow['shuttle']} | " + " · ".join(moves))
    else:
        print("\n(교환관계 사례 없음)")

    names4 = ("dev_liab", "불안정성(이진)", "sel_off", "dev_agg_hot")
    if regress:
        print("\n🔴 **역행** 사례 — 화물 인식 탐색의 1위가 심사 4축에서 **지배당한** 경우")
        print("  (화물 인식이 곧 심사 정렬이 아니라는 증거. 이유 가설: 구성체 불안정성은 "
              "길이로 **희석**되어 셔틀 쪽 기울기가 약해진다 — 가설이고 미측정)")
        for name, cname, lname, sd, prow, crow in regress:
            a, b = risk_axes(prow), risk_axes(crow)
            moves = [f"{names4[i]} {a[i]}→{b[i]}" for i in range(4) if a[i] != b[i]]
            print(f"  {name}·cargo{cname}·{lname}·seed{sd}: {prow['shuttle']} → "
                  f"{crow['shuttle']} | " + " · ".join(moves))
    print(f"\n무결성 — 심사 1위가 실제 게이트(CPP·근거없음)에 걸릴 사례: {gated_top} "
          f"(0 이어야 정상: `rank_key` 가 그런 멤버를 맨 뒤로 보내므로)")

    # ── E) 화물 의존성 — 불일치 자체가 화물에 반응하나 ──
    print("\n" + "=" * 104)
    print("표 E — **불일치의 화물 의존성**. 같은 (셔틀·시드·링커)에서 cargo A 와 B 의 결과가 같은가")
    print(_w("링커", 10) + _w("항목", 30) + _w("A와 B 동일", 12, True) + _w("전체", 8, True))
    for lname, linker in LINKERS:
        same_front = same_lex = same_cf = same_judge = 0
        tot = 0
        for name in TARGETS:
            for sd in SEEDS:
                tot += 1
                mem = prod[(name, sd)]
                fa = set(nondominated(mem, lambda s: v_construct(s, CARGO_A, linker)))
                fb = set(nondominated(mem, lambda s: v_construct(s, CARGO_B, linker)))
                same_front += (fa == fb)
                same_lex += (min(mem, key=lambda s: v_construct(s, CARGO_A, linker))
                             == min(mem, key=lambda s: v_construct(s, CARGO_B, linker)))
                same_cf += (set(cf[(name, "A", lname, sd)])
                            == set(cf[(name, "B", lname, sd)]))
                same_judge += (judge_best(mem, CARGO_A, linker)[0]
                               == judge_best(mem, CARGO_B, linker)[0])
        for label, v in (("재채점 전선 집합", same_front), ("재채점 사전식 1위", same_lex),
                         ("반사실 탐색 출력 집합", same_cf), ("심사 1위(프로덕션 아카이브)", same_judge)):
            print(_w(lname, 10) + _w(label, 30) + _w(v, 12, True) + _w(tot, 8, True))

    # ── 복원 단언 ──
    assert o.objectives is _PROD_OBJECTIVES
    assert o.OFF_BAND == 0.05, o.OFF_BAND
    print(f"\n복원 확인: objectives={o.objectives.__name__} · OFF_BAND={o.OFF_BAND}")


if __name__ == "__main__":
    main()
