"""보수성(conservatism) 스윕 — **"축을 빼고 제약을 조이는 대가"를 곡선으로** 잰다.

심사 대상 명제 (가정이 아니라 **피고**)
----------------------------------------
이 프로젝트의 최근 결론들은 한 방향을 향한다 — *축을 순위에서 빼고(§19 BBB · §26 off)
제약을 조여도(밴드 ±0.05 · 보존 위치 동결) **대가는 거의 없고 선별 능력은 남는다**.*
그런데 그 주장을 떠받치는 측정은 **1차원 절단면 하나**다: `docs/experiments.md` §26-5 의
off 밴드 폭 표 — 셔틀 **하나**(Angiopep-2) · 5시드 · `rmt_floor` 고정. 이 스크립트는
그 점 하나를 **보수성에 대한 곡선**으로, 그리고 **셔틀 3종**으로 넓힌다.

🔴 **관계는 이미 비단조로 측정돼 있다.** §26-5 에서 가장 보수적인 설정(밴드 ±0.00)은
개발성을 **크게** 망쳤다 — 최선 목적 `(1, 0.077, 45.4)` → `(2, 0.651, 57.1)`, 전선 55→11.
즉 "조여도 공짜"는 **전 구간에서 성립하지 않는다**. 이 스윕의 목적은 그 명제를 **확인**
하는 것이 아니라 **어디까지 공짜이고 어디서 값을 치르기 시작하는지 지점을 찍는 것**이다.
수치가 명제를 반박하면 반박한 대로 적는다 (`paper-ideas-not-prose` · 지표 게이밍 금지).

보수성 격자 (한 점을 고정하지 않고 **쓸어 본다**)
--------------------------------------------------
  · off-target 제약 : 없음 / 밴드 ±0.10 / ±0.05(현재 기본값) / ±0.02 / ±0.00
  · 문헌 보존 위치 동결 : off / on  (`core.evidence.CONSERVED_POSITIONS`)
  · `rmt_floor` : 0.80 / 0.85(기본) / 0.90
  → 5 × 2 × 3 = **30 셀** × 셔틀 3종 × 시드 5개(0~4). 셀별 전선은 **5시드 합산**이다
    (§25-2·§26-5 와 같은 집계라 숫자를 바로 나란히 놓을 수 있다).

세 가지 결과 계열 (한 셀마다 전부 측정)
----------------------------------------
⑴ **해킹률** — 전선 멤버가 Angiopep-2 의 K10/K15 를 바꾼 비율 · 후보당 **전하 삭제 수**
   (부모의 K/R 자리가 K/R 이 아니게 된 개수) · **R11 삭제율**. K10/K15·R11 은 문헌이
   위치를 집어 준 Angiopep-2 에만 적용된다(§25-5). 전하 삭제는 **3종 전부** 보고한다.
⑵ **전선 품질** — 전선 크기, 그리고 세 목적(liability · 응집 핫스팟 · 불안정성)에서
   **도달한 최선값**. ⚠️ 축별 최소값(componentwise min)이라 한 멤버의 벡터가 아니다 —
   §26-5 의 "최선 목적"(사전식 1위 멤버)과 비교할 수 있게 그 벡터도 함께 찍는다.
⑶ 🔴 **선별 유지율** — 이 스윕이 새로 넣는 핵심 측정. 전선 멤버가 `n_edits = k` 라면
   **같은 셔틀 · 같은 k 치환**의 무작위 변이체를 생성해(위치·잔기 무작위, 시드 고정),
   그 멤버가 세 목적에서 **지배하는** 비율을 센다(`_dominates` — `pareto_scrub` 과 같은
   정의). 전선 멤버 평균. 질문이 직접적이다: *조이고 나서도 옵티마이저의 출력이 같은
   크기의 무작위 편집을 여전히 이기나?* (`matched-null-principle`)
   같이 센다: **근거 등급 상속 실패 수** — `core.evidence.inherit` 가 등급을 못 주는
   멤버는 하류에서 탈락한다. **등급을 잃고 살아남은 전선은 선별이 아니다.**

⚠️ 보증하지 않는 것: 무작위 대응 널은 제약(floor·밴드·동결)을 **받지 않는다**. 그게
   의도다 — "제약 없는 같은 크기 편집"이 기준선이다. 따라서 선별 유지율은 "제약 안에서
   최선인가"가 아니라 "**제약을 받고도 무제약 무작위를 이기나**"를 잰다.

⚠️ `core/` 는 **건드리지 않는다.** 밴드·동결은 스크립트에서 모듈 속성을
   몽키패치해 바꾸고(`finally` 복원), 실행 끝에 **원래 값으로 돌아왔는지 단언**한다.
   프로덕션 기본값(`OFF_BAND = 0.05` · `_frozen_idx` 원본)은 그대로다.

재현(서브프로세스 0회·완전 결정론, Biopython 필요 — `python3 -I` 금지):
  $ /usr/bin/python3 scripts/sweep_conservatism.py
빠른 반복용 축소 격자:
  $ /usr/bin/python3 scripts/sweep_conservatism.py --quick
"""
from __future__ import annotations

import argparse
import random
import sys
import unicodedata
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import core.developability_opt as dopt
from core.config import SHUTTLES
from core.developability import require_protein_analysis
from core.evidence import CONSERVED_POSITIONS, INHERIT_SIM_MIN, inherit

_STD = "ACDEFGHIKLMNPQRSTVWY"
_TOL = 1e-9

# ── 심사할 셔틀 3종 — **성격이 다르게** 고른다 ───────────────────────────────
# 선정 기준(측정값, `off_target`/`objectives` 실측):
#   · Angiopep-2 (19aa, off 0.09, liab 3, 불안정 71.5) — 라이브러리에서 **문헌 보존
#     위치가 있는 유일한 셔틀**(K10·K15). 동결 축을 잴 수 있는 곳은 여기뿐이다.
#   · B6 (9aa, off **0.370** = rmt_sim 이 살아 있는 RMT 중 최고, liab 2, 응집 -0.553)
#     — **고위험 off** 쪽 끝. 짧고 응집 핫스팟이 이미 음수라 "내려갈 데가 적은" 셔틀.
#   · RVG29 (29aa, off **0.060** = 최저 비(非)영, liab **6** = 최다) — **저위험 off ·
#     개발성 부채 최다** 쪽 끝. 길어서 편집 여지가 가장 넓다.
# 제외: THR·Tet1 은 `rmt_sim = 0.0`(§26-8) 이라 어떤 floor 에서도 전선이 시드 하나로
#   붕괴한다 — 보수성 효과가 아니라 참조셋 결손을 재게 된다. Leptin30 은 off 가 **정확히
#   0.000**(축의 바닥)이라 밴드가 한쪽만 작동한다 — 바닥 효과가 밴드 효과로 오독된다
#   (`null_developability` 가 같은 함정을 '바닥 제외'로 처리한 것과 같은 이유).
SEED_SHUTTLES = ("Angiopep-2", "B6", "RVG29")

#: off 제약 격자. None = 제약 없음(옛 3축 상태)
BANDS: tuple[float | None, ...] = (None, 0.10, 0.05, 0.02, 0.00)
FREEZES = (False, True)
FLOORS = (0.80, 0.85, 0.90)

#: 프로덕션 기본값 — 이 셀이 현재 파이프라인과 **같은 설정**이다
PROD_CELL = (0.05, True, 0.85)

#: 선별 유지율의 셀·멤버당 대응 무작위 변이체 수
N_MATCHED = 40


def _w(s: str, width: int, right: bool = False) -> str:
    """표시 폭 기준 패딩 — 한글·이모지는 2칸을 차지하므로 `f"{s:<16}"` 이 어긋난다."""
    disp = sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in s)
    pad = " " * max(0, width - disp)
    return pad + s if right else s + pad


def band_label(b: float | None) -> str:
    return "없음" if b is None else f"±{b:.2f}"


def _rng_for(*parts) -> random.Random:
    """문자열 키 → 결정론 난수. ⚠️ `hash()` 는 세션마다 달라지므로 crc32 를 쓴다."""
    key = "|".join(str(p) for p in parts).encode()
    return random.Random(zlib.crc32(key))


# ---------------------------------------------------------------------------
# 해킹률 측정자
# ---------------------------------------------------------------------------
def charge_deleted(parent: str, child: str) -> int:
    """부모의 **K/R 자리**가 K/R 이 아니게 된 개수(= 전하 삭제 치환 수).

    길이가 같을 때만 위치가 1:1 이다(`_mutate` 는 길이 보존이므로 전선에서는 항상 같다).
    """
    if len(parent) != len(child):
        return 0
    return sum(1 for a, b in zip(parent, child) if a in "KR" and b not in "KR")


def mutated_positions(parent: str, child: str, pos1: tuple[int, ...]) -> bool:
    """1-based 위치 중 하나라도 바뀌었나."""
    if len(parent) != len(child):
        return False
    return any(parent[p - 1] != child[p - 1] for p in pos1 if 1 <= p <= len(parent))


def r11_deleted(parent: str, child: str) -> bool:
    """Angiopep-2 전용 — R11(index 10)이 K/R 이 아니게 됐나 (§25-5 의 해킹면)."""
    if len(parent) != len(child) or len(parent) < 11:
        return False
    return parent[10] in "KR" and child[10] not in "KR"


# ---------------------------------------------------------------------------
# 🔴 선별 유지율 — 같은 크기 무작위 편집 대응 널
# ---------------------------------------------------------------------------
def matched_variant(parent: str, k: int, rng: random.Random) -> str:
    """부모에서 **정확히 k 치환**된 무작위 변이체.

    ⚠️ 잔기를 다시 뽑아 **반드시 달라지게** 한다. `rng.choice(_STD)` 를 그냥 쓰면 1/20
    확률로 같은 잔기가 나와 실효 편집 수가 k 보다 작아지고, 그러면 "같은 크기"가 아니다
    (전선 멤버의 `n_edits` 는 **실제로 다른 자리 수**이므로 짝이 어긋난다).
    제약(floor·밴드·동결)은 **일부러 걸지 않는다** — 무제약 무작위가 기준선이다.
    """
    s = list(parent)
    for i in rng.sample(range(len(s)), min(k, len(s))):
        cur = s[i]
        while True:
            c = rng.choice(_STD)
            if c != cur:
                s[i] = c
                break
    return "".join(s)


def retention(parent: str, member: str, k: int, n: int,
              rng: random.Random) -> float:
    """멤버가 **같은 크기 무작위 변이체**를 세 목적에서 지배하는 비율.

    지배 정의는 `dopt._dominates` — `pareto_scrub` 이 전선을 만들 때 쓰는 **그 함수**다.
    (비교 불가(incomparable)는 지배가 아니므로 분자에 들어가지 않는다 → 보수적 측정.)
    """
    mobj = dopt.objectives(member)
    hit = 0
    for _ in range(n):
        v = matched_variant(parent, k, rng)
        if dopt._dominates(mobj, dopt.objectives(v)):
            hit += 1
    return hit / n


# ---------------------------------------------------------------------------
# 한 셀 = (밴드, 동결, floor) × 셔틀 — 5시드 합산
# ---------------------------------------------------------------------------
def run_cell(shuttle_name: str, seq: str, band: float | None, freeze: bool,
             floor: float, seeds: tuple[int, ...], rounds: int, pop: int,
             n_matched: int) -> dict:
    """전선을 만들고 세 결과 계열을 센다. 모듈 속성 몽키패치는 **여기서만** 한다."""
    orig_band, orig_frozen = dopt.OFF_BAND, dopt._frozen_idx
    try:
        # 🔴 프로덕션 기본값은 건드리지 않는다 — 아래 finally 에서 반드시 복원되고,
        #    main() 끝에서 원래 값과 같은지 **단언**한다.
        dopt.OFF_BAND = float("inf") if band is None else band
        if not freeze:
            dopt._frozen_idx = lambda _s: frozenset()
        members: list[dict] = []
        for sd in seeds:
            members += dopt.pareto_scrub(seq, rounds=rounds, pop=pop,
                                         rmt_floor=floor, seed=sd)
    finally:
        dopt.OFF_BAND, dopt._frozen_idx = orig_band, orig_frozen

    edited = [m for m in members if m["n_edits"] > 0]
    cons = CONSERVED_POSITIONS.get(shuttle_name, ())

    # ⑵ 전선 품질 — 축별 최소값 + (§26-5 비교용) 사전식 1위 멤버의 벡터
    objs = [tuple(m["objectives"][k] for k in dopt.OBJECTIVE_NAMES) for m in members]
    best_each = tuple(min(o[i] for o in objs) for i in range(3)) if objs else None
    lex_best = min(objs) if objs else None

    # ⑴ 해킹률
    n_cd = [charge_deleted(seq, m["seq"]) for m in edited]
    cons_hits = sum(1 for m in edited if cons and mutated_positions(seq, m["seq"], cons))
    r11 = sum(1 for m in edited if shuttle_name == "Angiopep-2"
              and r11_deleted(seq, m["seq"]))

    # ⑶ 선별 유지율 + 근거 등급 상속 실패
    rets, zero = [], 0
    for idx, m in enumerate(edited):
        rng = _rng_for(shuttle_name, band, freeze, floor, idx, m["seq"])
        r = retention(seq, m["seq"], m["n_edits"], n_matched, rng)
        rets.append(r)
        zero += r <= _TOL
    tier_fail = sum(1 for m in members if inherit(m["seq"])[0] is None)

    offs = [m["off_target"] for m in members]
    return {
        "shuttle": shuttle_name, "band": band, "freeze": freeze, "floor": floor,
        "front": len(members), "edited": len(edited),
        "best_each": best_each, "lex_best": lex_best,
        "cd_per_cand": (sum(n_cd) / len(n_cd)) if n_cd else 0.0,
        "cons_rate": (cons_hits / len(edited)) if edited and cons else None,
        "r11_rate": (r11 / len(edited)) if edited and shuttle_name == "Angiopep-2" else None,
        "retention": (sum(rets) / len(rets)) if rets else None,
        "ret_zero": zero,
        "tier_fail": tier_fail,
        "off_lo": min(offs) if offs else None, "off_hi": max(offs) if offs else None,
        "n_edits_mean": (sum(m["n_edits"] for m in edited) / len(edited)) if edited else 0.0,
    }


# ---------------------------------------------------------------------------
# 출력
# ---------------------------------------------------------------------------
_MAIN_HDR = (_w("밴드", 7) + _w("동결", 5) + _w("floor", 6, True)
             + _w("전선", 6, True) + _w("편집", 6, True)
             + _w("최선liab", 9, True) + _w("최선응집", 10, True)
             + _w("최선불안정", 12, True)
             + _w("전하삭제/후보", 15, True) + _w("K10/15", 8, True)
             + _w("R11삭제", 9, True) + _w("🔴선별유지율", 14, True)
             + _w("유지0", 7, True) + _w("등급탈락", 10, True)
             + _w("off 범위", 13, True))


def _pct(x: float | None) -> str:
    return "—" if x is None else f"{x:.1%}"


def main_row(r: dict) -> str:
    be = r["best_each"]
    mark = " ◀기본" if (r["band"], r["freeze"], r["floor"]) == PROD_CELL else ""
    return (_w(band_label(r["band"]), 7) + _w("on" if r["freeze"] else "off", 5)
            + f"{r['floor']:>6.2f}{r['front']:>6}{r['edited']:>6}"
            + (f"{be[0]:>9}{be[1]:>10.3f}{be[2]:>12.1f}" if be else _w("—", 31, True))
            + f"{r['cd_per_cand']:>15.2f}" + _w(_pct(r["cons_rate"]), 8, True)
            + _w(_pct(r["r11_rate"]), 9, True) + _w(_pct(r["retention"]), 14, True)
            + f"{r['ret_zero']:>7}" + f"{r['tier_fail']:>10}"
            + _w(f"{r['off_lo']:.2f}–{r['off_hi']:.2f}" if r["off_lo"] is not None
                 else "—", 13, True) + mark)


def marginal(rows: list[dict], key: str, label: str, levels, fmt,
             note: str = "") -> None:
    """한 인자만 남기고 **나머지를 접은** 요약 — 각 인자의 모양을 따로 본다."""
    print(f"\n[주변 요약 {label}] — 나머지 인자(와 셔틀)를 모두 접었다")
    print(_w(label, 10) + _w("셀", 5, True) + _w("전선합", 8, True)
          + _w("전선/셀", 9, True) + _w("최선liab", 9, True)
          + _w("최선응집", 10, True) + _w("최선불안정", 12, True)
          + _w("전하삭제/후보", 15, True) + _w("R11삭제", 9, True)
          + _w("🔴선별유지율", 14, True) + _w("등급탈락", 10, True))
    for lv in levels:
        sub = [r for r in rows if r[key] == lv]
        if not sub:
            continue
        fr = sum(r["front"] for r in sub)
        be = [r["best_each"] for r in sub if r["best_each"]]
        cd = [r["cd_per_cand"] for r in sub if r["edited"]]
        r11 = [r["r11_rate"] for r in sub if r["r11_rate"] is not None]
        ret = [r["retention"] for r in sub if r["retention"] is not None]
        print(_w(fmt(lv), 10) + f"{len(sub):>5}{fr:>8}{fr/len(sub):>9.1f}"
              + f"{min(b[0] for b in be):>9}{min(b[1] for b in be):>10.3f}"
              + f"{min(b[2] for b in be):>12.1f}"
              + f"{sum(cd)/len(cd) if cd else 0:>15.2f}"
              + _w(_pct(sum(r11)/len(r11)) if r11 else "—", 9, True)
              + _w(_pct(sum(ret)/len(ret)) if ret else "—", 14, True)
              + f"{sum(r['tier_fail'] for r in sub):>10}")
    if note:
        print(note)


def per_shuttle_matrix(rows: list[dict], key: str, label: str, levels, fmt,
                       shuttle_names: tuple[str, ...]) -> None:
    """🔴 **셔틀별로 답이 다른지**를 보이는 행렬 — 접으면 사라지는 신호다."""
    print(f"\n  ·· {label} × 셔틀 (선별유지율 / 전선/셀) — 접으면 사라지는 차이")
    print("  " + _w(label, 10) + "".join(_w(n, 22, True) for n in shuttle_names))
    for lv in levels:
        cells = []
        for n in shuttle_names:
            sub = [r for r in rows if r[key] == lv and r["shuttle"] == n]
            ret = [r["retention"] for r in sub if r["retention"] is not None]
            fr = sum(r["front"] for r in sub) / len(sub) if sub else 0.0
            cells.append(_w(f"{_pct(sum(ret)/len(ret)) if ret else '—'} / {fr:.1f}",
                            22, True))
        print("  " + _w(fmt(lv), 10) + "".join(cells))


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="보수성(off 밴드 · 보존 위치 동결 · rmt_floor) 스윕 — "
                    "해킹률 · 전선 품질 · 🔴 선별 유지율")
    ap.add_argument("--seeds", type=str, default="0,1,2,3,4",
                    help="pareto_scrub 시드 목록 (기본 0,1,2,3,4 — §25-2·§26-5 와 동일)")
    ap.add_argument("--quick", action="store_true",
                    help="축소 격자(밴드 3 × 동결 2 × floor 1, 시드 0~2) — 반복용")
    ap.add_argument("--matched", type=int, default=N_MATCHED,
                    help=f"선별 유지율의 멤버당 대응 무작위 변이체 수 (기본 {N_MATCHED})")
    ap.add_argument("--rounds", type=int, default=40, help="pareto_scrub rounds (기본 40)")
    ap.add_argument("--pop", type=int, default=30, help="pareto_scrub pop (기본 30)")
    args = ap.parse_args()

    # 🔴 수천 번의 호출 **전에** 의존성을 한 번 확인한다. Biopython 이 없으면 불안정성·
    #    GRAVY 가 조용히 대체값으로 떨어져 "보수성은 공짜"라는 **거짓 결론**이 나온다
    #    (`silent-engine-degradation`). 메시지는 require_protein_analysis 단일 출처.
    require_protein_analysis("보수성 스윕(불안정성 지수·GRAVY·off-target)")

    seeds = tuple(int(x) for x in args.seeds.replace(" ", "").split(",") if x != "")
    bands, freezes, floors = BANDS, FREEZES, FLOORS
    if args.quick:
        bands, floors = (None, 0.05, 0.00), (0.85,)
        if args.seeds == "0,1,2,3,4":
            seeds = (0, 1, 2)

    # 🔴 복원 단언용 원본 — 몽키패치 전에 잡아 둔다
    prod_band, prod_frozen = dopt.OFF_BAND, dopt._frozen_idx

    print("=" * 150)
    print("보수성 스윕 — 축 제거·제약 강화의 대가를 **곡선**으로 (§26-5 의 1차원 절단면 확장)")
    print(f"격자: off 밴드 {[band_label(b) for b in bands]} × 보존위치동결 {list(freezes)}"
          f" × rmt_floor {list(floors)} = {len(bands)*len(freezes)*len(floors)} 셀")
    print(f"셔틀 {len(SEED_SHUTTLES)}종: " + " · ".join(
        f"{n}(off {dopt.off_target(SHUTTLES[n]['seq']):.3f} · "
        f"{len(SHUTTLES[n]['seq'])}aa · liab {dopt.objectives(SHUTTLES[n]['seq'])[0]})"
        for n in SEED_SHUTTLES))
    print(f"시드 **{list(seeds)}** ({len(seeds)}회, 전선은 시드 합산) · "
          f"pareto_scrub(rounds={args.rounds}, pop={args.pop}, max_archive=40 기본값"
          f"{' — 축소 없음' if (args.rounds, args.pop) == (40, 30) else ' ⚠️ 전 셀에 동일하게 축소'})")
    print(f"선별 유지율: 멤버당 대응 무작위 변이체 {args.matched}개(같은 셔틀·같은 k 치환·무제약)"
          f" · 지배 정의 = pareto_scrub 의 _dominates")
    print(f"프로덕션 기본값 = 밴드 ±{prod_band:.2f} · 동결 on · floor 0.85 (표에서 ◀기본)")
    print("⚠️ 해킹률·선별유지율은 **실제로 편집된 멤버**(n_edits≥1)만의 통계다 — 시드 자신은"
          " 전선에 들어가지만 편집이 0이라 분모에서 뺀다('편집' 열이 그 분모).")
    print("⚠️ 전선은 **시드 합산**이라 시드 서열 자신이 시드 수만큼 중복 계수된다"
          "(§26-5 와 같은 집계). '편집' 열이 0 이면 전선이 **시드 하나로 붕괴**한 것이다.")
    print("⚠️ 보존 위치가 없는 셔틀(B6·RVG29)은 동결 on/off 가 **같은 실행**이다 — 행이"
          " 동일한 것은 버그가 아니라 문헌 간극이고, 동결 주변 요약을 2/3 희석한다.")
    print("⚠️ '최선liab/응집/불안정'은 **축별 최소값**(componentwise min)이고 한 멤버의 벡터가"
          " 아니다 — §26-5 와 같은 '사전식 1위 멤버' 벡터는 셔틀별 꼬리에 따로 찍는다.")
    print("=" * 150)

    rows: list[dict] = []
    for name in SEED_SHUTTLES:
        seq = SHUTTLES[name]["seq"]
        print(f"\n{'─' * 150}")
        print(f"■ {name}  `{seq}`  ({len(seq)}aa · 시드 목적 {dopt.objectives(seq)}"
              f" · 시드 off {dopt.off_target(seq):.3f}"
              f" · 문헌 보존 위치 {CONSERVED_POSITIONS.get(name) or '없음(= 어디가 중요한지 모름)'})")
        print("─" * 150)
        print(_MAIN_HDR)
        for band in bands:
            for freeze in freezes:
                for floor in floors:
                    r = run_cell(name, seq, band, freeze, floor, seeds,
                                 args.rounds, args.pop, args.matched)
                    rows.append(r)
                    print(main_row(r))
        sub = [r for r in rows if r["shuttle"] == name]
        print("  사전식 1위 멤버 벡터(§26-5 '최선 목적'과 같은 정의): "
              + " · ".join(f"{band_label(b)}/{'on' if f else 'off'}/{fl:.2f}="
                           f"{r['lex_best']}"
                           for b, f, fl in ((PROD_CELL[0], True, 0.85),
                                            (None, True, 0.85), (0.00, True, 0.85))
                           for r in sub
                           if (r["band"], r["freeze"], r["floor"]) == (b, f, fl)))

    # ── 세 개의 주변 요약 ────────────────────────────────────────────────────
    print("\n" + "=" * 150)
    print("주변 요약 — 각 인자의 **모양을 따로** 본다 (나머지 인자는 접음)")
    print("=" * 150)
    marginal(rows, "band", "off 밴드", bands, band_label)
    per_shuttle_matrix(rows, "band", "off 밴드", bands, band_label, SEED_SHUTTLES)
    marginal(rows, "freeze", "보존위치동결", freezes, lambda v: "on" if v else "off",
             note="🔴 이 요약은 **1/3 만 유효하다** — 문헌 보존 위치가 있는 셔틀은"
                  " Angiopep-2 하나이고 나머지 2종은 on/off 가 동일 실행이다."
                  " 동결의 효과는 아래 행렬의 **Angiopep-2 열에서만** 읽을 것(§25-5).")
    per_shuttle_matrix(rows, "freeze", "보존위치동결", freezes,
                       lambda v: "on" if v else "off", SEED_SHUTTLES)
    marginal(rows, "floor", "rmt_floor", floors, lambda v: f"{v:.2f}")
    per_shuttle_matrix(rows, "floor", "rmt_floor", floors, lambda v: f"{v:.2f}",
                       SEED_SHUTTLES)

    # ── 🔴 복원 단언: 프로덕션 기본값이 그대로인가 ──────────────────────────
    assert dopt.OFF_BAND == prod_band, f"OFF_BAND 복원 실패: {dopt.OFF_BAND} != {prod_band}"
    assert dopt._frozen_idx is prod_frozen, "_frozen_idx 복원 실패 (몽키패치가 남았다)"
    assert dopt.OFF_BAND == 0.05, f"OFF_BAND 가 프로덕션 기본값 0.05 가 아니다: {dopt.OFF_BAND}"
    assert dopt._frozen_idx(SHUTTLES["Angiopep-2"]["seq"]) == frozenset({9, 14}), \
        "_frozen_idx 가 Angiopep-2 의 K10/K15 를 더 이상 동결하지 않는다"
    print("\n✅ 복원 단언 통과 — OFF_BAND = 0.05 · `_frozen_idx` 원본 "
          "(Angiopep-2 → {9,14} 동결) · `core/` 파일 **무수정**")

    interpret(rows, bands, freezes, floors, seeds)


# ---------------------------------------------------------------------------
def interpret(rows, bands, freezes, floors, seeds) -> None:
    """숫자에서만 읽는 해석 — 명제를 확인하러 온 게 아니라 **지점을 찍으러** 왔다."""
    def pick(**kw):
        return [r for r in rows if all(r[k] == v for k, v in kw.items())]

    def avg(sub, key):
        xs = [r[key] for r in sub if r[key] is not None]
        return sum(xs) / len(xs) if xs else float("nan")

    base = pick(band=PROD_CELL[0], freeze=True, floor=0.85)
    none = pick(band=None, freeze=True, floor=0.85)
    tight = pick(band=0.00, freeze=True, floor=0.85)

    print("\n" + "=" * 150)
    print("해석 — **측정된 숫자에서만**. 🔴 §26-5 가 이미 비단조를 보였으므로 '조이면 공짜'를"
          " 확인하려 들지 않는다")
    print("=" * 150)
    # ── 하네스 검정: §26-5 의 **발표된 행**을 재현하나 ────────────────────────
    # 이 스윕은 §26-5 를 확장한 것이므로, 겹치는 절단면(Angiopep-2 · 동결 on ·
    # floor 0.85 · 시드 0~4)에서 **같은 숫자**가 나와야 한다. 안 나오면 확장이 아니라
    # 다른 측정이다 — 그래서 해석보다 먼저 찍는다.
    pub = {None: (55, 0.15, (1, 0.077, 45.4)), 0.05: (35, 0.20, (1, 0.161, 38.9)),
           0.02: (30, 0.00, (1, 0.237, 63.0)), 0.00: (11, 0.00, (2, 0.651, 57.1))}
    if tuple(seeds) != (0, 1, 2, 3, 4):
        # 발표 수치는 **시드 0~4 합산**이다 — 다른 시드 집합과 비교하면 안 된다.
        pub = {}
        print(f"  [하네스 검정 생략] §26-5 발표 행은 **시드 0~4 합산**이고 지금은"
              f" {list(seeds)} 다 — 비교 불가(축소 격자의 수치를 §26-5 와 나란히 쓰지 말 것).")
    else:
        print("  [하네스 검정] §26-5 발표 행 재현"
              " (Angiopep-2 · 동결 on · floor 0.85 · 시드 0~4)")
    if pub:
        print("  " + _w("밴드", 8) + _w("전선(발표→측정)", 18, True)
              + _w("R11(발표→측정)", 18, True) + _w("사전식최선(발표)", 24, True)
              + _w("(측정)", 22, True) + _w("일치", 8, True))
    for b, (fr, r11, lex) in pub.items():
        got = pick(band=b, freeze=True, floor=0.85, shuttle="Angiopep-2")
        if not got:
            continue
        g = got[0]
        ok = (g["front"] == fr and abs((g["r11_rate"] or 0) - r11) <= 0.01
              and g["lex_best"] == lex)
        print("  " + _w(band_label(b), 8) + _w(f"{fr}→{g['front']}", 18, True)
              + _w(f"{r11:.0%}→{_pct(g['r11_rate'])}", 18, True)
              + _w(str(lex), 24, True) + _w(str(g["lex_best"]), 22, True)
              + _w("✅" if ok else "⚠️ 불일치", 8, True))

    print("   (⑴⑵ 는 **동결 on · floor 0.85 · 셔틀 3종 평균** 절단면이다 — 주변 요약표는"
          " 동결·floor 까지 접은 값이라 수치가 다르다.)")
    print(f"⑴ **어디까지 공짜인가** — 밴드를 없음 → ±0.05(기본값)로 조일 때, 전선/셀 "
          f"{sum(r['front'] for r in none)/max(1,len(none)):.1f} → "
          f"{sum(r['front'] for r in base)/max(1,len(base)):.1f}, "
          f"선별유지율 {avg(none,'retention'):.1%} → {avg(base,'retention'):.1%}, "
          f"전하삭제/후보 {avg(none,'cd_per_cand'):.2f} → {avg(base,'cd_per_cand'):.2f}.")
    print(f"⑵ **어디서 값을 치르기 시작하나** — 밴드 ±0.00(최대 보수)에서 전선/셀 "
          f"{sum(r['front'] for r in tight)/max(1,len(tight)):.1f}, "
          f"선별유지율 {avg(tight,'retention'):.1%}, "
          f"전하삭제/후보 {avg(tight,'cd_per_cand'):.2f} — "
          "위 주변 요약의 '최선응집·최선불안정' 열과 함께 읽을 것(해킹률이 0 이 되는 대신"
          " 개발성 축에서 값을 치르는지가 거기 보인다).")
    # 🔴 평균이 지우는 것 — **셔틀별로 답이 다르다.** 같은 절단면을 셔틀로 쪼개 찍는다.
    print("  🔴 셔틀별 분기(동결 on · floor 0.85) — 평균 한 줄로 말하면 안 되는 이유")
    print("  " + _w("셔틀", 14) + _w("밴드", 8) + _w("전선", 6, True)
          + _w("편집", 6, True) + _w("최선liab", 9, True) + _w("최선응집", 10, True)
          + _w("최선불안정", 12, True) + _w("선별유지율", 12, True))
    for nm in SEED_SHUTTLES:
        for b in (None, 0.05, 0.00):
            for r in pick(band=b, freeze=True, floor=0.85, shuttle=nm):
                be = r["best_each"]
                print("  " + _w(nm, 14) + _w(band_label(b), 8)
                      + f"{r['front']:>6}{r['edited']:>6}"
                      + (f"{be[0]:>9}{be[1]:>10.3f}{be[2]:>12.1f}" if be
                         else _w("—", 31, True))
                      + _w(_pct(r["retention"]), 12, True))

    allret = [r["retention"] for r in rows if r["retention"] is not None]
    worst = min((r for r in rows if r["retention"] is not None),
                key=lambda r: r["retention"], default=None)
    print(f"⑶ **선별 유지율이 붕괴하는가** — 전 셀 범위 "
          f"{min(allret):.1%} ~ {max(allret):.1%}"
          + (f", 최저는 {worst['shuttle']} / 밴드 {band_label(worst['band'])} / 동결 "
             f"{'on' if worst['freeze'] else 'off'} / floor {worst['floor']:.2f}"
             f" ({worst['retention']:.1%}, 유지 0 멤버 {worst['ret_zero']}/{worst['edited']})"
             if worst else "")
          + ". 선별유지율 0% 셀 "
          + f"{sum(1 for r in rows if r['retention'] is not None and r['retention'] <= _TOL)}"
          + f"/{len(allret)} · 전선이 시드 하나로 붕괴한 셀(편집 0) "
          + f"{sum(1 for r in rows if r['edited'] == 0)}/{len(rows)}.")
    tf = sum(r["tier_fail"] for r in rows)
    print(f"⑷ **등급을 잃고 살아남은 전선인가** — 근거 등급 상속 실패 멤버 총 {tf}개"
          f" (상속 하한 {INHERIT_SIM_MIN}). 0 이면 전선은 '등급을 잃어서' 살아남은 것이"
          " 아니다 — 선별로 읽어도 되는 조건 하나가 충족된다.")
    print("⚠️ 이 스윕이 **말하지 않는 것**: ⓐ 선별유지율은 세 개발성 축에 대한 상대 비교이고"
          " 실제 제조가능성 라벨이 아니다(외부 라벨은 §27 에서 깨졌다). ⓑ 보존 위치 동결의"
          " 효과는 Angiopep-2 **한 종**에서만 측정된다 — 문헌이 위치를 집어 준 셔틀이"
          " 하나뿐이기 때문이고, 이 간극은 성과가 아니라 **한계**다(§25-5).")
    print(f"재현: /usr/bin/python3 scripts/sweep_conservatism.py --seeds {','.join(map(str, seeds))}")


if __name__ == "__main__":
    main()
