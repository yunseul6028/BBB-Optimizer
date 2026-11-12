"""두 번째 화물(cargo B) — **결론이 화물 하나를 넘어서 성립하나**.

왜 돌렸나 (대기열 #2)
---------------------
공개된 모든 수치는 화물 **하나**(`core.config.DEFAULT_CARGO = GSNKGAIIGLM`,
Aβ(25-35))로 측정됐다. `docs/protocol.md` §5 보고 규약은 이 사실을 명시적인 한계로
적는다 — *"벤치마크 전체가 화물 하나로 측정됐고, 그 사실이 축 하나를 구조적으로
측정 불가로 만들었다"*(§23-6). 이 스크립트는 그 한계를 **측정으로** 닫는다.

세 가지를 잰다
--------------
⑴ **전선 구조** — `reachable_front` 의 전선·`front_top` 크기와 **살아남는 셔틀**.
   대표 셔틀이 바뀌나.
⑵ **4축이 애초에 화물에 민감한가** — 전수 조합(링커 13 × 셔틀 12 = 156)에서 축별로
   얼마나 움직이는지. 🔴 **움직이지 않는 축은 "움직이지 않았다"고 적는다** — 그 자체가
   화물이 무엇에 영향을 주고 무엇에 주지 않는지에 대한 발견이다.
⑶ **§29 「공짜 vs 유료」 경계가 이동하나** — 보수성 스윕의 작은 절단면을 cargo B 로
   다시 돌린다(밴드 none/±0.05/±0.02 × 동결 on/off × floor 0.85 × 셔틀 3종 ×
   시드 0~4 · `rounds=40, pop=30` = §29 와 같은 설정).

하네스 검증 (필수)
------------------
cargo B 를 해석하기 **전에** §29/§26-5 의 공개 수치를 내 하네스로 재현한다
(Angiopep-2 · 밴드 ±0.05 · 동결 on · floor 0.85 → 전선 **35** · R11 삭제 **20.0%** ·
사전식 1위 **(1, 0.161, 38.9)** · off 범위 **0.05–0.12**). 재현되지 않으면 하네스가
틀린 것이고 cargo B 해석은 무효다.

제약
----
`core/` 는 **건드리지 않는다**. 파라미터 변주는 모듈 속성 몽키패치로만 하고
(`scripts/sweep_conservatism.run_cell` 과 같은 패턴) `finally` 에서 복원하고
`main()` 끝에서 프로덕션 기본값과 같은지 **단언**한다.

실행: `/usr/bin/python3 scripts/second_cargo.py`   (⚠️ `python3 -I` 금지)
"""

from __future__ import annotations

import argparse
import statistics
import sys
import unicodedata
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import core.config as cfg
import core.developability_opt as dopt
from core.config import LINKER_LIBRARY, SHUTTLES, get_settings
from core.developability import assess_developability, require_protein_analysis
from core.optimizer_agent import OptimizationAgent, rank_key, risk_axes
from scripts.sweep_conservatism import band_label, run_cell

# ---------------------------------------------------------------------------
# 시드 — 이 스크립트의 **모든** 난수 출처를 여기 적는다
# ---------------------------------------------------------------------------
#: `pareto_scrub` 시드 집합. §25-2·§26-5·§29 와 **같은 집합**이라 수치를 나란히 둘 수 있다.
SCRUB_SEEDS: tuple[int, ...] = (0, 1, 2, 3, 4)
#: 선별 유지율의 대응 무작위 변이체 RNG 은 `sweep_conservatism._rng_for` 가
#: (셔틀·밴드·동결·floor·인덱스·서열) → crc32 로 **결정론 파생**한다(세션 불변).
MATCHED_RNG = "crc32(shuttle|band|freeze|floor|idx|seq) — 결정론·세션 불변"
#: 멤버당 대응 무작위 변이체 수(§29 와 동일)
N_MATCHED = 40

# ---------------------------------------------------------------------------
# 화물 — A(기본) vs B(두 번째)
# ---------------------------------------------------------------------------
#: 🔴 cargo B 는 **측정 전에** 골랐다(사후 선택 금지). 선정 기준은 두 가지다:
#:   ① **실재하는 명명된 CNS 펩타이드**여야 한다(무작위 문자열 금지).
#:   ② `assess_developability` 실측에서 cargo A 와 **물리화학적으로 반대쪽**이어야 한다.
#: DSIP(delta sleep-inducing peptide, Trp-Ala-Gly-Gly-Asp-Ala-Ser-Gly-Glu)는
#: 1977년 토끼 뇌 투석액에서 분리된 내인성 나노펩타이드이고 수면·금단·통증 맥락에서
#: 임상적으로 다뤄진 서열이다. 길이가 9aa 로 cargo A(11aa)와 **거의 맞아** 길이 효과를
#: 통제한 채 전하·소수성·응집을 전부 뒤집는다 — 아래 표의 수치가 근거다.
CARGO_A = cfg.DEFAULT_CARGO                 # Aβ(25-35) · 11aa
CARGO_B = "WAGGDASGE"                       # DSIP · 9aa
CARGO_NAMES = {CARGO_A: "Aβ(25-35) (기본)", CARGO_B: "DSIP (두 번째)"}

#: 같이 재서 **고르지 않은 이유**를 남기는 대조군(보고만. 측정에는 쓰지 않는다)
ALSO_CONSIDERED = {
    "GLP-1(7-36)": "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR",
    "Substance P": "RPKPQQFFGLM",
    "davunetide(NAP)": "NAPVSIPQ",
    "PACAP-27": "HSDGIFTDSYSRYRKQMAVKKYLAAVL",
}

#: §29 슬라이스 — 과제 명세대로 **작게** 자른다(밴드 3 × 동결 2 × floor 1 × 셔틀 3)
SLICE_BANDS: tuple[float | None, ...] = (None, 0.05, 0.02)
SLICE_FREEZES = (False, True)
SLICE_FLOOR = 0.85
SLICE_SHUTTLES = ("Angiopep-2", "B6", "RVG29")

#: 🔴 하네스 검증 기준 — `docs/experiments.md` §29 머리말 + §26-5 표의 **공개 수치**
#:   (Angiopep-2 · 동결 on · floor 0.85 · 시드 0~4 · rounds 40 · pop 30)
PUBLISHED = {
    None:  {"front": 55, "r11": 0.145},
    0.05:  {"front": 35, "r11": 0.200},
    0.02:  {"front": 30, "r11": 0.000},
}
PUBLISHED_LEX = {0.05: (1, 0.161, 38.9), None: (1, 0.077, 45.4), 0.02: (1, 0.237, 63.0)}

_AXES = ("dev_liab", "instability", "sel_off", "dev_agg_hot")


# ---------------------------------------------------------------------------
# 표시 보조
# ---------------------------------------------------------------------------
def _w(s: str, width: int, right: bool = False) -> str:
    """표시 폭 패딩 — 한글·이모지는 2칸이라 `f"{s:<16}"` 이 어긋난다."""
    disp = sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in s)
    pad = " " * max(0, width - disp)
    return pad + s if right else s + pad


def _pct(x: float | None) -> str:
    return "—" if x is None else f"{x:.1%}"


def _spearman(xs: list[float], ys: list[float]) -> float | None:
    """순위상관 — 외부 의존성 없이. 동점은 평균 순위."""
    n = len(xs)
    if n < 3:
        return None

    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return None if dx == 0 or dy == 0 else round(num / (dx * dy), 3)


# ---------------------------------------------------------------------------
# ⓐ 화물 특성 표
# ---------------------------------------------------------------------------
def cargo_table() -> None:
    print("\n" + "=" * 118)
    print("ⓐ 화물 특성 — **왜 이 둘이 대조인가** (전부 `assess_developability` 실측, 외부 수치 아님)")
    print("=" * 118)
    print(_w("화물", 26) + _w("서열", 33) + _w("len", 5, True) + _w("순전하", 8, True)
          + _w("전하밀도", 10, True) + _w("GRAVY", 8, True) + _w("pI", 6, True)
          + _w("응집핫스팟", 12, True) + _w("불안정성", 10, True) + _w("liab", 6, True))
    rows = [(CARGO_NAMES[CARGO_A], CARGO_A), (CARGO_NAMES[CARGO_B], CARGO_B)]
    rows += [(f"  (미채택) {k}", v) for k, v in ALSO_CONSIDERED.items()]
    for name, seq in rows:
        d = assess_developability(seq)
        print(_w(name, 26) + _w(seq, 33) + f"{len(seq):>5}{d.net_charge:>8.1f}"
              + f"{d.charge_density:>10.2f}{d.gravy:>8.2f}{d.pI:>6.1f}"
              + f"{d.agg_hotspot:>12.2f}{dopt._instability(seq):>10.1f}{d.n_liabilities:>6}")
    da, db = assess_developability(CARGO_A), assess_developability(CARGO_B)
    print(f"\n  A liability: {da.liabilities}")
    print(f"  B liability: {db.liabilities}")
    print("  → 대조의 내용: 길이는 11 vs 9 로 **거의 맞춰** 길이 효과를 통제하고,"
          f" 순전하 {da.net_charge:+.1f}→{db.net_charge:+.1f}(부호 반전) ·"
          f" GRAVY {da.gravy:+.2f}→{db.gravy:+.2f}(부호 반전) ·"
          f" pI {da.pI:.1f}→{db.pI:.1f} ·"
          f" 응집 핫스팟 {da.agg_hotspot:.2f}→{db.agg_hotspot:.2f}(절반) 로 네 축을 전부 뒤집는다.")
    print("  → 미채택 이유: GLP-1(7-36)·PACAP-27 은 길이가 30/27aa 라 **길이와 화학이"
          " 섞인다**. Substance P 는 11aa 로 길이는 맞지만 C말단 `…GLM` 이 cargo A 와"
          " 같아 **접합부 화학이 동일**해져 대조가 약해진다. davunetide 는 순전하 0 ·"
          " GRAVY −0.06 으로 전하·소수성 대조가 약하다.")


# ---------------------------------------------------------------------------
# ⓑ 전선 구조
# ---------------------------------------------------------------------------
def front_structure(agent: OptimizationAgent) -> dict:
    print("\n" + "=" * 118)
    print("ⓑ 전선 구조 — `reachable_front(cargo, refresh=True)` · 라이브러리 전수"
          f" {len(LINKER_LIBRARY)}링커 × {len(SHUTTLES)}셔틀 = {len(LINKER_LIBRARY)*len(SHUTTLES)}조합")
    print("=" * 118)
    print("🔴 `refresh=True` 필수 — 전선은 `.cache/reachable/` 에 캐시된다.")
    fronts = {}
    for cargo in (CARGO_A, CARGO_B):
        fronts[cargo] = agent.reachable_front(cargo, refresh=True)

    def shuttle_set(rows):
        out = []
        for r in rows:
            nm = next((k for k, v in SHUTTLES.items() if v["seq"] == r["shuttle"]), r["shuttle"])
            if nm not in out:
                out.append(nm)
        return out

    print(_w("항목", 26) + _w(CARGO_NAMES[CARGO_A], 34) + _w(CARGO_NAMES[CARGO_B], 34))
    def row(label, fa, fb):
        print(_w(label, 26) + _w(str(fa), 34) + _w(str(fb), 34))
    A, B = fronts[CARGO_A], fronts[CARGO_B]
    row("채점 조합", A["n_total"], B["n_total"])
    row("비독성(하드 게이트 통과)", A["n_nontoxic"], B["n_nontoxic"])
    row("후보군(근거 A·B·C·非CPP)", A["n_gate"], B["n_gate"])
    row("근거 그룹 분포", A["by_group"], B["by_group"])
    row("최상위 도달 그룹", A["best_group"], B["best_group"])
    row("전선 크기", len(A["front"]), len(B["front"]))
    row("front_top 크기", len(A["front_top"]), len(B["front_top"]))
    row("전선 셔틀", shuttle_set(A["front"]), shuttle_set(B["front"]))
    row("front_top 셔틀", shuttle_set(A["front_top"]), shuttle_set(B["front_top"]))
    rep_a = sorted(A["front_top"], key=rank_key)[0] if A["front_top"] else None
    rep_b = sorted(B["front_top"], key=rank_key)[0] if B["front_top"] else None
    row("대표(rank_key 1위)", rep_a and rep_a["label"], rep_b and rep_b["label"])
    row("  그 축 (liab,II,off,agg)", rep_a and risk_axes(rep_a), rep_b and risk_axes(rep_b))
    row("축별 최선", A["best"], B["best"])

    print("\n  front_top 전원:")
    for tag, f in (("A", A), ("B", B)):
        for r in sorted(f["front_top"], key=rank_key):
            print(f"    [{tag}] {_w(r['label'], 26)} 등급 {_w(str(r['tier']), 4)}"
                  f" liab {r['dev_liab']} · II {r['instability']}"
                  f" · off {r['sel_off']:.2f} · agg {r['dev_agg_hot']:.2f}"
                  f" · BBB(보고) {r['bbb_report']:.3f}")
    return fronts


# ---------------------------------------------------------------------------
# ⓒ 4축의 화물 민감도
# ---------------------------------------------------------------------------
def axis_sensitivity(agent: OptimizationAgent) -> None:
    print("\n" + "=" * 118)
    print("ⓒ 4축이 애초에 **화물에 민감한가** — 전수 조합을 두 화물로 각각 채점해 축별로 비교")
    print("=" * 118)
    combos = [{"label": f"{ln}+{sn}", "linker": LINKER_LIBRARY[ln]["seq"],
               "shuttle": SHUTTLES[sn]["seq"]}
              for ln in LINKER_LIBRARY for sn in SHUTTLES]
    rows = {}
    for cargo in (CARGO_A, CARGO_B):
        acc = []
        for i in range(0, len(combos), 24):
            _, chunk = agent._evaluate(cargo, {"constructs": combos[i:i + 24]})
            acc.extend(chunk)
        rows[cargo] = {r["label"]: r for r in acc}
    labels = [c["label"] for c in combos if c["label"] in rows[CARGO_A] and c["label"] in rows[CARGO_B]]

    print(_w("축", 16) + _w("계산 범위", 22) + _w("A 평균", 10, True) + _w("B 평균", 10, True)
          + _w("변한 조합", 11, True) + _w("평균Δ", 9, True) + _w("평균|Δ|", 10, True)
          + _w("최대|Δ|", 10, True) + _w("ρ(A,B)", 9, True))
    scope = {"dev_liab": "구성체 전체 ⟵화물", "instability": "구성체 전체 ⟵화물",
             "sel_off": "셔틀 단독", "dev_agg_hot": "설계모듈(링커+셔틀)"}
    for ax in _AXES:
        a = [float(rows[CARGO_A][l][ax] or 0.0) for l in labels]
        b = [float(rows[CARGO_B][l][ax] or 0.0) for l in labels]
        d = [y - x for x, y in zip(a, b)]
        nz = sum(1 for v in d if abs(v) > 1e-9)
        print(_w(ax, 16) + _w(scope[ax], 22)
              + f"{statistics.fmean(a):>10.3f}{statistics.fmean(b):>10.3f}"
              + f"{nz:>6}/{len(labels):<5}"
              + f"{statistics.fmean(d):>9.3f}{statistics.fmean(abs(v) for v in d):>10.3f}"
              + f"{max(abs(v) for v in d):>10.3f}"
              + _w(str(_spearman(a, b)), 9, True))
    # 순위에 실제로 들어가는 값 — 불안정성은 **이진**이다
    ba = [risk_axes(rows[CARGO_A][l])[1] for l in labels]
    bb = [risk_axes(rows[CARGO_B][l])[1] for l in labels]
    print(_w("instability(이진)", 16) + _w("40 기준 · 순위 축", 22)
          + f"{statistics.fmean(ba):>10.3f}{statistics.fmean(bb):>10.3f}"
          + f"{sum(1 for x, y in zip(ba, bb) if x != y):>6}/{len(labels):<5}"
          + _w("—", 29, True) + _w(str(_spearman(ba, bb)), 9, True))
    # 하드 게이트(독성)와 보고치(BBB)
    for ax, scp in (("tox", "구성체 전체 ⟵화물"), ("bbb", "bbb_scoring_seq ⟵화물")):
        a = [float(rows[CARGO_A][l][ax]) for l in labels]
        b = [float(rows[CARGO_B][l][ax]) for l in labels]
        d = [y - x for x, y in zip(a, b)]
        print(_w(f"({ax}: 축 아님)", 16) + _w(scp, 22)
              + f"{statistics.fmean(a):>10.3f}{statistics.fmean(b):>10.3f}"
              + f"{sum(1 for v in d if abs(v)>1e-9):>6}/{len(labels):<5}"
              + f"{statistics.fmean(d):>9.3f}{statistics.fmean(abs(v) for v in d):>10.3f}"
              + f"{max(abs(v) for v in d):>10.3f}" + _w(str(_spearman(a, b)), 9, True))
    tfa = sum(1 for l in labels if rows[CARGO_A][l]["toxic"])
    tfb = sum(1 for l in labels if rows[CARGO_B][l]["toxic"])
    flip = sum(1 for l in labels if rows[CARGO_A][l]["toxic"] != rows[CARGO_B][l]["toxic"])
    print(f"\n  독성 하드 게이트 탈락: A {tfa} · B {tfb} · 판정 뒤집힘 {flip}/{len(labels)}")
    # 4축 전체 벡터가 같은 조합 수 = 화물이 순위에 아무 영향도 못 주는 조합
    same = sum(1 for l in labels
               if risk_axes(rows[CARGO_A][l]) == risk_axes(rows[CARGO_B][l]))
    print(f"  4축 벡터가 **완전히 동일**한 조합: {same}/{len(labels)} ({same/len(labels):.1%})")
    print("  🔴 구조적 사실: `sel_off` 는 셔틀 단독, `dev_agg_hot` 는 설계모듈(링커+셔틀) 로"
          " 계산된다(`optimizer_agent._evaluate`) — 즉 **정의상 화물에 불변**이다."
          " 화물이 닿는 순위 축은 `dev_liab` 와 `instability` **둘뿐**이다.")


def agg_scope_check() -> None:
    """🔴 §20-4 의 **근거 자체**가 화물 특이적인가 — 서브프로세스 0회·결정론.

    §20-4(= `risk_axes` 주석 ②)는 응집 순위 축을 '전체 구성체'에서 '설계 모듈'로 좁힌
    근거를 이렇게 적는다: *"화물은 고정이고 이 화물(GSNKGAIIGLM)은 `AIIGLM` 때문에 자체
    핫스팟이 0.80 이다. 전체 구성체로 재면 12종 셔틀 **전부 0.80** 으로 묶여 축이 죽는다."*
    그 측정은 화물 **하나**로 돌았다. 두 번째 화물에서도 축이 죽는지 확인한다.
    """
    print("\n" + "=" * 118)
    print("ⓒ-2 🔴 §20-4 의 근거를 두 번째 화물로 — **'전체 구성체로 재면 축이 죽는다'가"
          " 화물 특이적인가**")
    print("=" * 118)
    combos = [(LINKER_LIBRARY[l]["seq"], SHUTTLES[s]["seq"])
              for l in LINKER_LIBRARY for s in SHUTTLES]
    print(_w("범위", 34) + _w("평균", 9, True) + _w("최소", 8, True) + _w("최대", 8, True)
          + _w("0.80 동점", 11, True) + _w("판정", 16, True))
    for cargo in (CARGO_A, CARGO_B):
        v = [assess_developability(cargo + l + s).agg_hotspot for l, s in combos]
        tie = sum(1 for x in v if abs(x - 0.80) < 1e-9)
        dead = tie >= 0.9 * len(v)
        print(_w(f"전체 구성체 · {CARGO_NAMES[cargo]}", 34)
              + f"{statistics.fmean(v):>9.3f}{min(v):>8.2f}{max(v):>8.2f}"
              + f"{tie:>6}/{len(v):<5}"
              + _w("🔴 축이 죽는다" if dead else "🟢 축이 산다", 16, True))
    m = [assess_developability(l + s).agg_hotspot for l, s in combos if (l + s)]
    print(_w("설계모듈(현행 순위 축·화물 불변)", 34)
          + f"{statistics.fmean(m):>9.3f}{min(m):>8.2f}{max(m):>8.2f}"
          + _w("—", 11, True) + _w("—", 16, True))
    print("  → §20-4 의 '전부 0.80 으로 묶인다'는 **cargo A 에서만** 성립한다. 결정(모듈로"
          " 좁히기)은 길이 해킹 면역·'우리가 바꿀 수 있는 모듈' 논거로 여전히 서지만,"
          " **그 결정의 근거로 인용된 수치는 일반 명제가 아니다.**")


# ---------------------------------------------------------------------------
# ⓓ §29 슬라이스 — 하네스 검증 먼저
# ---------------------------------------------------------------------------
def slice_sweep(rounds: int, pop: int) -> list[dict]:
    print("\n" + "=" * 118)
    print("ⓓ §29 보수성 슬라이스 — 하네스 검증 → cargo B 재실행")
    print("=" * 118)
    print(f"격자: off 밴드 {[band_label(b) for b in SLICE_BANDS]} × 동결 {list(SLICE_FREEZES)}"
          f" × rmt_floor {SLICE_FLOOR} × 셔틀 {list(SLICE_SHUTTLES)}"
          f" · 시드 {list(SCRUB_SEEDS)} · rounds={rounds}, pop={pop} (§29 와 동일)")

    print("\n[검증] §29 머리말 + §26-5 공개 수치와 대조 (Angiopep-2 · 동결 on · floor 0.85)")
    print(_w("밴드", 8) + _w("전선(실측)", 12, True) + _w("전선(공개)", 12, True)
          + _w("R11(실측)", 11, True) + _w("R11(공개)", 11, True)
          + _w("사전식1위(실측)", 22, True) + _w("사전식1위(공개)", 22, True) + _w("판정", 8, True))
    ok = True
    seq_a2 = SHUTTLES["Angiopep-2"]["seq"]
    for band in SLICE_BANDS:
        r = run_cell("Angiopep-2", seq_a2, band, True, SLICE_FLOOR,
                     SCRUB_SEEDS, rounds, pop, N_MATCHED)
        exp = PUBLISHED[band]
        hit = (r["front"] == exp["front"]
               and abs((r["r11_rate"] or 0) - exp["r11"]) < 0.0051
               and r["lex_best"] == PUBLISHED_LEX[band])
        ok = ok and hit
        print(_w(band_label(band), 8) + f"{r['front']:>12}{exp['front']:>12}"
              + _w(_pct(r["r11_rate"]), 11, True) + _w(_pct(exp["r11"]), 11, True)
              + _w(str(r["lex_best"]), 22, True) + _w(str(PUBLISHED_LEX[band]), 22, True)
              + _w("🟢 일치" if hit else "🔴 불일치", 8, True))
    print(f"  → 하네스 검증: {'🟢 **재현됨**' if ok else '🔴 **재현 실패 — cargo B 해석 무효**'}")
    if not ok:
        print("  재현되지 않았다. 여기서 멈춘다(프로토콜: 하네스가 틀렸으면 cargo B 는 읽지 않는다).")
        return []

    print("\n🔴 **구조적 발견 — 이 스윕은 화물을 보지 않는다.**")
    print("   `run_cell` → `dopt.pareto_scrub(셔틀서열, ...)` 이고, 그 안의"
          " `objectives()`·`off_target()`·`shuttle_similarity()`·`_frozen_idx()` 는"
          " **전부 셔틀 서열만** 받는다. 화물은 어떤 인자에도 들어가지 않는다"
          "(`grep DEFAULT_CARGO core/*.py` → `config.py` 정의 한 줄뿐).")
    print("   따라서 §29 의 「공짜 vs 유료」 경계는 화물이 바뀌어도 **이동할 수 없다**."
          " '변화를 못 찾았다'가 아니라 **변할 것이 없다**. 아래에서 수치로 확인한다.")

    rows_a, rows_b = [], []
    for name in SLICE_SHUTTLES:
        seq = SHUTTLES[name]["seq"]
        for band in SLICE_BANDS:
            for freeze in SLICE_FREEZES:
                rows_a.append(run_cell(name, seq, band, freeze, SLICE_FLOOR,
                                       SCRUB_SEEDS, rounds, pop, N_MATCHED))
                # cargo B 조건 — 화물을 바꿔서 **같은 함수를 다시** 돌린다
                orig = cfg.DEFAULT_CARGO
                try:
                    cfg.DEFAULT_CARGO = CARGO_B
                    rows_b.append(run_cell(name, seq, band, freeze, SLICE_FLOOR,
                                           SCRUB_SEEDS, rounds, pop, N_MATCHED))
                finally:
                    cfg.DEFAULT_CARGO = orig

    hdr = (_w("셔틀", 12) + _w("밴드", 7) + _w("동결", 5)
           + _w("전선 A", 8, True) + _w("전선 B", 8, True)
           + _w("유지율 A", 10, True) + _w("유지율 B", 10, True)
           + _w("전하삭제 A", 12, True) + _w("전하삭제 B", 12, True)
           + _w("R11 A", 8, True) + _w("R11 B", 8, True)
           + _w("동일?", 7, True))
    print("\n[슬라이스] cargo A vs cargo B — 같은 셀을 두 화물로")
    print(hdr)
    identical = True
    for ra, rb in zip(rows_a, rows_b):
        same = all(ra[k] == rb[k] for k in
                   ("front", "edited", "retention", "cd_per_cand", "r11_rate",
                    "lex_best", "best_each", "off_lo", "off_hi", "tier_fail"))
        identical = identical and same
        print(_w(ra["shuttle"], 12) + _w(band_label(ra["band"]), 7)
              + _w("on" if ra["freeze"] else "off", 5)
              + f"{ra['front']:>8}{rb['front']:>8}"
              + _w(_pct(ra["retention"]), 10, True) + _w(_pct(rb["retention"]), 10, True)
              + f"{ra['cd_per_cand']:>12.2f}{rb['cd_per_cand']:>12.2f}"
              + _w(_pct(ra["r11_rate"]), 8, True) + _w(_pct(rb["r11_rate"]), 8, True)
              + _w("=" if same else "🔴≠", 7, True))
    print(f"  → 전 {len(rows_a)}셀에서 모든 측정치가 "
          f"{'**비트 단위로 동일**' if identical else '🔴 **달랐다** — 구조 분석이 틀렸다'}.")
    return rows_a


# ---------------------------------------------------------------------------
# ⓔ 전선 멤버를 **두 화물로 다시 심사**
# ---------------------------------------------------------------------------
def rejudge_front(agent: OptimizationAgent, rounds: int, pop: int) -> None:
    """탐색은 화물을 못 봐도 **심사(4축)**는 본다. 같은 설계체를 두 화물로 채점한다."""
    print("\n" + "=" * 118)
    print("ⓔ 탐색이 화물을 못 보는 대가 — **같은 설계체를 두 화물로 다시 심사**")
    print("=" * 118)
    print("프로덕션 셀(밴드 ±0.05 · 동결 on · floor 0.85)의 전선 멤버를"
          " 직접융합 구성체로 만들어 4축을 두 화물로 채점한다.")
    print("🔴 **n = 시드 1개(seed=0)**다 — ⓓ 와 달리 시드 합산이 아니다. 멤버 수가"
          " 7/3/11 뿐이므로 아래 비율은 **참고치**이고 추세로만 읽는다.")
    print(_w("셔틀", 12) + _w("멤버", 6, True) + _w("liab A→B", 12, True)
          + _w("II A→B", 14, True) + _w("II이진 A→B", 13, True)
          + _w("4축 동일", 10, True) + _w("A내 지배관계", 14, True)
          + _w("B내 지배관계", 14, True) + _w("뒤집힘", 8, True))
    from core.optimizer_agent import dominates
    for name in SLICE_SHUTTLES:
        seq = SHUTTLES[name]["seq"]
        front = dopt.pareto_scrub(seq, rounds=rounds, pop=pop,
                                  rmt_floor=SLICE_FLOOR, seed=0)
        combos = [{"label": f"m{i}", "linker": "", "shuttle": m["seq"]}
                  for i, m in enumerate(front)]
        got = {}
        for cargo in (CARGO_A, CARGO_B):
            acc = []
            for i in range(0, len(combos), 24):
                _, ch = agent._evaluate(cargo, {"constructs": combos[i:i + 24]})
                acc.extend(ch)
            got[cargo] = acc
        ra, rb = got[CARGO_A], got[CARGO_B]
        la = [r["dev_liab"] for r in ra]; lb = [r["dev_liab"] for r in rb]
        ia = [r["instability"] for r in ra]; ib = [r["instability"] for r in rb]
        bi_a = [risk_axes(r)[1] for r in ra]; bi_b = [risk_axes(r)[1] for r in rb]
        same4 = sum(1 for x, y in zip(ra, rb) if risk_axes(x) == risk_axes(y))
        pa = sum(1 for i in range(len(ra)) for j in range(len(ra))
                 if i != j and dominates(ra[i], ra[j]))
        pb = sum(1 for i in range(len(rb)) for j in range(len(rb))
                 if i != j and dominates(rb[i], rb[j]))
        flip = sum(1 for i in range(len(ra)) for j in range(len(ra))
                   if i != j and dominates(ra[i], ra[j]) != dominates(rb[i], rb[j]))
        print(_w(name, 12) + f"{len(front):>6}"
              + _w(f"{statistics.fmean(la):.2f}→{statistics.fmean(lb):.2f}", 12, True)
              + _w(f"{statistics.fmean(ia):.1f}→{statistics.fmean(ib):.1f}", 14, True)
              + _w(f"{statistics.fmean(bi_a):.2f}→{statistics.fmean(bi_b):.2f}", 13, True)
              + _w(f"{same4}/{len(ra)}", 10, True) + f"{pa:>14}{pb:>14}" + f"{flip:>8}")
    print("  ⚠️ 링커는 **직접융합 고정**이다(링커 축을 섞지 않기 위해). 'A내/B내 지배관계'는"
          " 전선 멤버쌍 중 한쪽이 다른 쪽을 지배하는 쌍의 수 — 탐색의 3목적에서는 전부"
          " 비지배였는데, 화물을 붙여 4축으로 심사하면 지배가 **생긴다**.")


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description="두 번째 화물 — 결론의 화물 일반화 측정")
    ap.add_argument("--rounds", type=int, default=40, help="pareto_scrub rounds (기본 40 = §29)")
    ap.add_argument("--pop", type=int, default=30, help="pareto_scrub pop (기본 30 = §29)")
    ap.add_argument("--skip-front", action="store_true",
                    help="전선·축 민감도(서브프로세스 필요) 생략 — 스윕만")
    args = ap.parse_args()

    # 🔴 수천 번의 호출 **전에** 엔진 무결성 확인(R10 · silent-engine-degradation).
    require_protein_analysis("두 번째 화물(불안정성 지수·GRAVY·off-target)")
    st = get_settings()
    if not args.skip_front and not (st.use_deepb3p_local and st.use_toxinpred3_local):
        raise SystemExit("🔴 로컬 deepB3P/ToxinPred3 가 없다 — 전선 측정은 실측이어야 한다.")

    # 복원 단언용 원본 — 몽키패치 **전에** 잡아 둔다
    prod_band, prod_frozen, prod_cargo = dopt.OFF_BAND, dopt._frozen_idx, cfg.DEFAULT_CARGO

    print("=" * 118)
    print("두 번째 화물(cargo B) — 결론이 화물 하나를 넘어서 성립하나 (대기열 #2)")
    print("=" * 118)
    print(f"python  : {sys.version.split()[0]} @ {sys.executable}")
    try:
        import Bio
        print(f"Biopython: {Bio.__version__}")
    except Exception as e:  # noqa: BLE001
        raise SystemExit(f"🔴 Biopython 없음: {e}")
    print(f"🎲 시드  : pareto_scrub seeds = {list(SCRUB_SEEDS)}"
          f" · 대응 널 RNG = {MATCHED_RNG} · 멤버당 변이체 {N_MATCHED}개")
    print(f"화물 A  : {CARGO_A}  ({CARGO_NAMES[CARGO_A]})")
    print(f"화물 B  : {CARGO_B}  ({CARGO_NAMES[CARGO_B]})")
    print(f"엔진    : deepB3P(local)={st.use_deepb3p_local} · ToxinPred3(local)={st.use_toxinpred3_local}")
    print(f"프로덕션: OFF_BAND={prod_band} · DEFAULT_CARGO={prod_cargo}")

    cargo_table()

    agent = None
    if not args.skip_front:
        agent = OptimizationAgent(st)
        front_structure(agent)
        axis_sensitivity(agent)
    agg_scope_check()   # 서브프로세스 0회 — `--skip-front` 에서도 돈다

    slice_sweep(args.rounds, args.pop)

    if agent is not None:
        rejudge_front(agent, args.rounds, args.pop)

    # ── 🔴 프로덕션 기본값 복원 단언 ─────────────────────────────────────────
    assert dopt.OFF_BAND == prod_band == 0.05, f"OFF_BAND 오염: {dopt.OFF_BAND}"
    assert dopt._frozen_idx is prod_frozen, "_frozen_idx 몽키패치가 복원되지 않았다"
    assert cfg.DEFAULT_CARGO == prod_cargo == "GSNKGAIIGLM", \
        f"DEFAULT_CARGO 오염: {cfg.DEFAULT_CARGO}"
    # 복원이 **실효적**인지도 확인한다(원본은 Angiopep-2 의 K10·K15 → 0-based 9·14).
    assert dopt._frozen_idx(SHUTTLES["Angiopep-2"]["seq"]) == frozenset({9, 14}), \
        "동결 조회가 실효하지 않는다 — 몽키패치 잔존 의심"
    print("\n🟢 프로덕션 기본값 복원 확인: OFF_BAND=0.05 · _frozen_idx 원본 ·"
          f" DEFAULT_CARGO={cfg.DEFAULT_CARGO}")


if __name__ == "__main__":
    main()
