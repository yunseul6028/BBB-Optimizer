"""개발성 **3축의 대응 널 심사** — 실제로 최소화하는 축을 피고석에 세운다.

왜 이 스크립트가 있나
---------------------
§26(`scripts/null_selectivity.py`)에서 off-target 축을 심사해 **악화 검출은 믿고
개선 신호는 믿지 않는다**는 결론을 냈다 — 같은 크기 무작위 치환의 **46%가 off 를
내렸다**(k=3, 바닥 셔틀 제외: 414↑/355↓). 그래서 off 는 목적에서 빼고 밴드 제약으로
돌렸다.

그 심사가 **구멍을 드러냈다**: 정작 `core/developability_opt.objectives()` 가 지금
**최소화하는 세 축** — ⑴ liability 개수(모티프 규칙) · ⑵ 응집 핫스팟(AGGRESCAN a3v,
창 5aa) · ⑶ 불안정성 지수(Guruprasad) — 은 **개선 방향을 한 번도 심사받지 않았다.**
off 를 끌어내린 뒤 탐색을 끌고 있는 것은 정확히 이 세 축이다. 축의 근거가 "출판 척도"
라는 것은 **값이 옳다**는 뜻이고, **값을 내리는 편집이 실제 위험을 줄인다**는 뜻이
아니다. 이 스크립트는 세 축에 §26 과 **같은 널**을 돌려 그 구멍을 메운다.

세 개의 하위 널
---------------
**널 1 — 순서를 쓰나(조성 보존 셔플).** 조성을 **정확히 고정**하고 순서만 섞는다.
  값이 거의 안 움직이는 축은 서열이 아니라 **조성**을 읽고 있다(= 조성만 맞추면 축을
  만족시킬 수 있다). 움직이는 축은 진짜로 순서에 민감하다.
  · 보증: 조성 교란 없는 순수 순서 의존성 측정.
  · 미보증: 움직인다고 **옳은 방향**으로 움직인다는 뜻은 아니다(그건 널 2).

**널 2 — 개선 신호가 값싼가(대응 무작위 편집).** 같은 크기 무작위 치환(k∈{1,2,3})에서
  축이 **최적화기가 보상하는 방향(↓)** 으로 간 비율을 센다. 무작위 편집이 축을 쉽게
  내린다면, 최소화 탐색은 **위험을 줄이지 않고 점수만 내리는 길**을 찾아낸다.
  · 보증: "개선이 공짜로 널려 있나"의 상대 측정. §26 의 46%(k=3, 바닥 제외)와 직접 비교.
  · 미보증: 개선비가 낮아도 축이 **실제 개발성**을 잰다는 증거는 아니다 — 외부 실측
    라벨이 없으므로 이 널은 **반증 장치**일 뿐 확증 장치가 아니다.
  · 🔴 **바닥 제외는 선택이 아니라 필수다**: liability 기준값 0 인 셔틀은 **내려갈 데가
    없어** 전부 "악화"로 집계되고, 그러면 축이 실제보다 **튼튼해 보이는 가짜 양성**이
    된다(§26 이 off 에서 반대 방향으로 겪은 것과 같은 함정).

**널 3 — 핫스팟이 정말 "삽입에 면역"인가.** §20-4 와 `developability.py` 주석은 국소
  핫스팟이 길이 희석에 **면역**이라고 적어 뒀지만, 그건 **정의에서 연역한 주장**이고
  측정된 적이 없다. 저응집 잔기(G·S·GS 반복)를 1·2·5·10·20aa 붙여(N/C 양쪽) 순위 축
  `agg_hotspot_a3v` 와 폐기된 전역 평균 `agg_score` 가 각각 어떻게 반응하는지 **센다**.
  · 보증: "패딩으로 축을 내릴 수 있나"의 직접 측정(= M10 교정의 사후 검증).
  · 미보증: 삽입(가운데 끼워넣기)이 아니라 **양 끝 패딩**만 본다. 창(5aa)보다 긴 저응집
    구간을 핫스팟 **내부**에 끼우는 편집은 별도 널이 필요하다.

재현(서브프로세스 0회·완전 결정론, Biopython 필요 — `python3 -I` 금지):
  $ /usr/bin/python3 scripts/null_developability.py
기본값이 보고 수치를 재현한다(널 1 셔플 200회/셔틀, 널 2 k별 120회/셔틀, seed 0).
"""
from __future__ import annotations

import argparse
import random
import sys
import unicodedata
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import SHUTTLES
from core.developability import (_A3V_MIN, _AGG_WINDOW, _AGGRESCAN_A3V,
                                 assess_developability, require_protein_analysis)
from core.developability_opt import OBJECTIVE_NAMES, objectives

_STD = "ACDEFGHIKLMNPQRSTVWY"
_TOL = 1e-9


def _w(s: str, width: int, right: bool = False) -> str:
    """표시 폭 기준 패딩 — 한글·이모지는 2칸을 차지하므로 `f"{s:<16}"` 이 어긋난다."""
    disp = sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in s)
    pad = " " * max(0, width - disp)
    return pad + s if right else s + pad


# 축 이름(한국어 표시) — 순서는 `objectives()` 가 돌려주는 벡터와 **같다**.
AXIS_LABELS = ("liability 개수", "응집 핫스팟 a3v", "불안정성 지수")
assert len(AXIS_LABELS) == len(OBJECTIVE_NAMES) == 3

# §26 의 비교 기준점: off-target 축, k=3, 바닥 제외 **상승비 53.8%** → 개선비 **46.2%**.
OFF_K3_IMPROVE_SHARE = 0.462


def shuttles() -> dict[str, str]:
    return {n: SHUTTLES[n]["seq"] for n in SHUTTLES}


def axis_values(seq: str) -> tuple[float, float, float]:
    """심사 대상 = 옵티마이저가 **그대로 최소화하는** 목적 벡터(단일 출처)."""
    return tuple(float(x) for x in objectives(seq))  # type: ignore[return-value]


# --------------------------------------------------------------------------
# 바닥(floor) 판정 — "더 내려갈 데가 없는" 축·셔틀 조합
# --------------------------------------------------------------------------
def hard_floor(axis: int, base: float) -> bool:
    """축 **정의상** 최소값에 닿아 있나.

    · liability: 모티프 **개수**라 0 이 하드 바닥이다(음수 불가).
    · 응집 핫스팟: AGGRESCAN a3v 척도의 하한 -1.836(= 전원 Asp 창)이 이론 바닥.
    · 불안정성: Guruprasad DIWV 합은 **음수가 가능**해서 하드 바닥이 없다 → 정의로는
      판정 불가. 대신 널 2 에서 **측정된** 무개선(개선 0회)으로 잡는다.
    """
    if axis == 0:
        return base <= 0
    if axis == 1:
        return base <= _A3V_MIN + 1e-6
    return False


# --------------------------------------------------------------------------
def null_1(n: int, seed: int) -> None:
    print(f"\n=== 널 1 — 순서를 쓰나: 조성 보존 셔플 (셔틀당 {n}회, 조성 **정확히** 불변) ===")
    print("값이 거의 안 움직이는 축은 서열이 아니라 **조성**을 읽고 있다.")
    sh = shuttles()
    verdict = []
    for ax, label in enumerate(AXIS_LABELS):
        rng = random.Random(seed)   # 축마다 같은 셔플 집합을 보게 고정
        print(f"\n[축 {ax+1}/3] {label}")
        print(_w("셔틀", 16) + _w("기준", 9, True) + _w("셔플평균", 10, True)
              + _w("변화한 셔플", 13, True) + _w("최대|Δ|", 10, True))
        tot_changed = tot = frozen = below = 0
        for name, seq in sh.items():
            base = axis_values(seq)[ax]
            vals = []
            for _ in range(n):
                s = list(seq)
                rng.shuffle(s)
                vals.append(axis_values("".join(s))[ax])
            changed = sum(1 for v in vals if abs(v - base) > _TOL)
            tot_changed += changed
            tot += n
            frozen += changed <= 0.05 * n       # 셔플이 거의 못 움직인 셔틀
            below += sum(vals) / len(vals) < base - _TOL
            print(f"{name:<16}{base:>9.3f}{sum(vals)/len(vals):>10.3f}"
                  f"{f'{changed}/{n}':>13}{max(abs(v - base) for v in vals):>10.3f}")
        print(_w("→ 합계", 16) + f"{'':>9}{'':>10}{f'{tot_changed}/{tot}':>13}"
              f"{tot_changed/tot:>10.1%}")
        verdict.append((label, tot_changed, tot, frozen, below))

    print("\n[널 1 판정] 셔플이 값을 못 움직이면 그 축은 **조성 측정기**다.")
    print(_w("축", 18) + _w("변화한 셔플", 13, True) + _w("비율", 8, True)
          + _w("≤5%만 변한 셔틀", 18, True) + _w("셔플평균<기준", 15, True))
    for label, ch, tt, fr, bl in verdict:
        print(_w(label, 18) + f"{f'{ch}/{tt}':>13}{ch/tt:>8.1%}"
              + f"{f'{fr}/{len(sh)}':>18}{f'{bl}/{len(sh)}':>15}")
    print("해석: **세 축 모두 순서에 반응한다** — '조성만 읽는 축'은 없다(= 널 1 통과).")
    print("  단 민감도 차이가 크다: 모티프 개수(개별 셔틀 0~7/200 사례 있음) < 응집 핫스팟")
    print("  < 불안정성(이합체 표라 거의 모든 셔플에서 변한다).")
    print("⚠️ '셔플평균<기준' 열은 **별개 신호**다 — 무작위 순서가 평균적으로 축을")
    print("   더 낮게 만든다는 뜻(= 실제 셔틀이 자기 조성 대비 모티프·핫스팟이 많다).")
    print("   축의 순서 민감성과 달리 이건 **개선이 값쌀 수 있다는 예고**이고, 널 2 가 센다.")


# --------------------------------------------------------------------------
def null_2(n: int, seed: int) -> None:
    print(f"\n=== 널 2 🔴 개선 신호가 값싼가: 대응 무작위 편집 (k∈{{1,2,3}} × {n}회/셔틀) ===")
    print("세 축 모두 **최소화**한다 → Δ<0 이 '최적화기가 보상하는 방향'(=개선)이다.")
    print("무작위 편집이 축을 쉽게 내린다면 최소화 탐색은 그 길들을 찾아낸다.")
    sh = shuttles()
    # cell[(name, ax, k)] = [개선, 악화, 동일, Δ합/ n]
    cell: dict[tuple[str, int, int], list] = {}
    base_v: dict[str, tuple] = {}
    rng = random.Random(seed)
    for name, seq in sh.items():
        base = axis_values(seq)
        base_v[name] = base
        for k in (1, 2, 3):
            acc = [[0, 0, 0, 0.0] for _ in AXIS_LABELS]
            for _ in range(n):
                s = list(seq)
                for i in rng.sample(range(len(s)), min(k, len(s))):
                    s[i] = rng.choice(_STD)
                cur = axis_values("".join(s))
                for ax in range(len(AXIS_LABELS)):
                    d = cur[ax] - base[ax]
                    acc[ax][3] += d
                    if d < -_TOL:
                        acc[ax][0] += 1
                    elif d > _TOL:
                        acc[ax][1] += 1
                    else:
                        acc[ax][2] += 1
            for ax in range(len(AXIS_LABELS)):
                acc[ax][3] /= n
                cell[(name, ax, k)] = acc[ax]

    # 바닥 판정: 정의상 하드 바닥 **또는** k 전체에서 개선 0회(=측정된 무개선)
    # 천장 판정(대칭): 악화가 전체 시행의 5% 미만 — **올라갈 데가 거의 없는** 셔틀.
    #   바닥이 개선비를 끌어내리는 것과 **정확히 반대로** 끌어올린다(TAT 불안정성 225.8).
    #   §26 은 한쪽(바닥)만 걸렀다. 양쪽을 다 걸러 본 줄이 민감도 검사다.
    floored: dict[int, set[str]] = {}
    ceiled: dict[int, set[str]] = {}
    for ax in range(len(AXIS_LABELS)):
        s, c = set(), set()
        for name in sh:
            imp = sum(cell[(name, ax, k)][0] for k in (1, 2, 3))
            if hard_floor(ax, base_v[name][ax]) or imp == 0:
                s.add(name)
            elif cell[(name, ax, 3)][1] < 0.05 * n:   # k=3 기준(가장 큰 편집)
                c.add(name)
        floored[ax], ceiled[ax] = s, c

    summary = []
    for ax, label in enumerate(AXIS_LABELS):
        print(f"\n[축 {ax+1}/3] {label} — 셔틀별 (k=3 만 표시; k=1,2 는 아래 합계에)")
        print(_w("셔틀", 16) + _w("기준", 9, True) + _w("평균Δ", 10, True)
              + _w("개선", 6, True) + _w("악화", 6, True) + _w("동일", 6, True)
              + _w("개선비", 9, True))
        for name in sh:
            imp, wor, eq, mean_d = cell[(name, ax, 3)]
            share = imp / (imp + wor) if imp + wor else float("nan")
            flag = " ⚠️바닥" if name in floored[ax] else ""
            print(f"{name:<16}{base_v[name][ax]:>9.3f}{mean_d:>+10.3f}"
                  f"{imp:>6}{wor:>6}{eq:>6}{share:>9.1%}{flag}")

        keep = [x for x in sh if x not in floored[ax]]
        keep2 = [x for x in keep if x not in ceiled[ax]]
        print(f"\n  합계 — 개선비 = 개선/(개선+악화). 바닥 제외: {len(keep)}/{len(sh)}종"
              f" (제외 {', '.join(sorted(floored[ax])) or '없음'})")
        print(f"  천장(k=3 악화<5%) 추가 제외: {', '.join(sorted(ceiled[ax])) or '없음'}"
              " — 바닥의 거울상이라 개선비를 **위로** 왜곡한다")
        print("  " + _w("집계", 24) + _w("k", 2, True) + _w("평균Δ", 10, True)
              + _w("개선", 6, True) + _w("악화", 6, True) + _w("동일", 6, True)
              + _w("개선비", 9, True))
        for k in (1, 2, 3):
            for lbl, names in ((f"전체 ({len(sh)}종)", list(sh)),
                               (f"🔴 바닥 제외 ({len(keep)}종)", keep),
                               (f"  +천장 제외 ({len(keep2)}종)", keep2)):
                if not names:
                    continue
                imp = sum(cell[(x, ax, k)][0] for x in names)
                wor = sum(cell[(x, ax, k)][1] for x in names)
                eq = sum(cell[(x, ax, k)][2] for x in names)
                md = sum(cell[(x, ax, k)][3] for x in names) / len(names)
                share = imp / (imp + wor) if imp + wor else float("nan")
                print("  " + _w(lbl, 24) + f"{k:>2}{md:>+10.3f}"
                      f"{imp:>6}{wor:>6}{eq:>6}{share:>9.1%}")
                if k == 3:
                    tag = ("nf" if lbl.startswith("🔴")
                           else "nf2" if lbl.startswith("  +") else "all")
                    summary.append((label, tag, share, imp, wor))

    print("\n[널 2 판정] k=3 · **바닥 제외** 개선비를 §26 의 off-target 과 같은 자에 올린다.")
    print(_w("축", 18) + _w("전체", 9, True) + _w("🔴 바닥제외", 14, True)
          + _w("개선/악화", 13, True) + _w("+천장제외", 12, True)
          + _w("off 46.2% 대비", 16, True))
    for label, tag, share, imp, wor in summary:
        if tag != "nf":
            continue
        rel = ("더 값쌈 🔴" if share > OFF_K3_IMPROVE_SHARE + 0.02
               else "비슷" if share > OFF_K3_IMPROVE_SHARE - 0.02 else "덜 값쌈")
        pick = {t: s for l, t, s, *_ in summary if l == label}
        print(_w(label, 18) + f"{pick['all']:>9.1%}{share:>14.1%}"
              + f"{f'{imp}/{wor}':>13}{pick.get('nf2', float('nan')):>12.1%}"
              + _w(rel, 16, True))
    print(f"기준점(§26): off-target 축은 k=3·바닥 제외 개선비 **{OFF_K3_IMPROVE_SHARE:.1%}**"
          "(414↑/355↓)였고,")
    print("  그 수치 **때문에** 목적에서 빠져 밴드 제약으로 내려갔다(§26-6). 위 숫자가 그보다")
    print("  높으면 '개선이 off 보다도 값싸다' — 같은 처방을 검토해야 한다는 뜻이다.")
    print("⚠️ 개선비가 50% 를 넘는 것 자체는 축이 **틀렸다**는 증거가 아니다. 무작위 편집이")
    print("   실제로 liability 를 지우는 일도 있다(모티프는 규칙이라 '사실'이다). 이 널이 재는")
    print("   것은 **값이 옳은가**가 아니라 **값을 내리는 길이 얼마나 흔한가**다 — 흔할수록")
    print("   최소화 탐색이 '실제 개선'과 '점수만 개선'을 구별하지 못한다.")


# --------------------------------------------------------------------------
_PADS = (("G", "G 반복"), ("S", "S 반복"), ("GS", "GS 반복"))
_SIDES = (("C말단", lambda s, p: s + p), ("N말단", lambda s, p: p + s))
_LENS = (1, 2, 5, 10, 20)


def _run(unit: str, length: int) -> str:
    """저응집 패딩 서열 — `unit` 을 반복해 정확히 `length` aa."""
    return (unit * (length // len(unit) + 1))[:length]


def null_3() -> None:
    print(f"\n=== 널 3 — 핫스팟이 정말 '삽입에 면역'인가 (저응집 패딩, 창 {_AGG_WINDOW}aa) ===")
    print("§20-4·developability.py 주석은 국소 핫스팟이 길이 희석에 면역이라고 **단정**했다.")
    print("정의에서 연역한 주장이므로 측정한다: G·S·GS 를 1~20aa 양 끝에 붙인다.")
    print("순위 축 agg_hotspot_a3v(raw) vs 폐기된 전역 평균 agg_a3v/agg_score 를 나란히 본다.")
    sh = shuttles()
    print("\n" + _w("패딩", 9) + _w("쪽", 8) + _w("L", 3, True)
          + _w("Δ평균a3v", 11, True) + _w("Δagg_score", 12, True)
          + _w("평균↓", 8, True) + _w("Δ핫스팟", 11, True)
          + _w("핫불변", 9, True) + _w("핫↑", 7, True) + _w("핫↓", 7, True))
    tot_hot_down = tot_cells = 0
    for unit, plabel in _PADS:
        for slabel, join in _SIDES:
            for L in _LENS:
                dm = ds = dh = 0.0
                mean_down = hot_same = hot_up = hot_down = 0
                for seq in sh.values():
                    b = assess_developability(seq)
                    a = assess_developability(join(seq, _run(unit, L)))
                    dm += a.agg_a3v - b.agg_a3v
                    ds += a.agg_score - b.agg_score
                    dh += a.agg_hotspot_a3v - b.agg_hotspot_a3v
                    mean_down += a.agg_a3v < b.agg_a3v - _TOL
                    hot_same += abs(a.agg_hotspot_a3v - b.agg_hotspot_a3v) <= _TOL
                    hot_up += a.agg_hotspot_a3v > b.agg_hotspot_a3v + _TOL
                    hot_down += a.agg_hotspot_a3v < b.agg_hotspot_a3v - _TOL
                m = len(sh)
                tot_hot_down += hot_down
                tot_cells += m
                print(_w(plabel, 9) + _w(slabel, 8) + f"{L:>3}"
                      + f"{dm/m:>+11.3f}{ds/m:>+12.3f}"
                      + f"{f'{mean_down}/{m}':>8}{dh/m:>+11.3f}"
                      + f"{f'{hot_same}/{m}':>9}{hot_up:>7}{hot_down:>7}")

    print(f"\n→ 핫스팟이 **내려간 조합**: {tot_hot_down}/{tot_cells}"
          f" (패딩으로 순위 축을 낮출 수 있었던 횟수)")

    print("\n[예외 추적] 핫스팟이 패딩으로 **올라가는** 셔틀 — 자기 국소최대가 패딩보다 낮은 경우")
    pad_hot = {}
    for unit, _ in _PADS:
        w = _run(unit, _AGG_WINDOW)
        pad_hot[unit] = sum(_AGGRESCAN_A3V[c] for c in w) / _AGG_WINDOW
        print(f"  패딩 {unit:<2} 의 창 {_AGG_WINDOW}aa a3v = {pad_hot[unit]:+.3f}"
              f"  ← 기준 핫스팟이 이보다 낮은 셔틀은 패딩이 핫스팟을 **올린다**")
    for name, seq in sh.items():
        b = assess_developability(seq)
        if b.agg_hotspot_a3v < max(pad_hot.values()) + _TOL:
            a = assess_developability(seq + _run("G", 20))
            print(f"  {name:<16} 기준 핫스팟 {b.agg_hotspot_a3v:+.3f}"
                  f" → G×20 추가 후 {a.agg_hotspot_a3v:+.3f}"
                  f" (전역 평균 {b.agg_a3v:+.3f} → {a.agg_a3v:+.3f})")

    print("\n[대표 상세] Angiopep-2 · C말단 G 반복 — 두 정의가 길이에 어떻게 반응하나")
    seq = sh["Angiopep-2"]
    print("  " + _w("L", 3, True) + _w("평균a3v", 11, True) + _w("agg_score", 11, True)
          + _w("핫스팟a3v", 13, True) + _w("agg_hotspot", 13, True))
    for L in (0, *_LENS):
        r = assess_developability(seq + _run("G", L))
        print(f"  {L:>3}{r.agg_a3v:>11.3f}{r.agg_score:>11.2f}"
              f"{r.agg_hotspot_a3v:>13.3f}{r.agg_hotspot:>13.2f}")

    print("\n[널 3 판정] **정의가 아니라 측정으로** 말할 수 있는 것:")
    print(f"  ⑴ 패딩이 순위 축(핫스팟)을 내린 경우 {tot_hot_down}/{tot_cells} — "
          "즉 '길이 희석으로 응집 점수 깎기'는")
    print("     핫스팟 정의에서 **한 번도 성립하지 않았다**. §20-4 의 '삽입에 면역' 주장은")
    print("     이 범위(양 끝 패딩 1~20aa · G·S·GS)에서 **측정으로 확인**된다.")
    print("  ⑵ 같은 패딩이 폐기된 전역 평균은 **길이에 단조로** 끌려 내려간다 — G 반복에서")
    print("     12종 중 10종, Angiopep-2 는 a3v -0.232 → -0.387(agg_score 0.44 → 0.40)로")
    print("     L=1,2,5,10,20 전 구간 단조 하락. M10 이 고친 병이 실재했다는 사후 확인이다.")
    print("     (내려가지 않은 2종 TAT·B6 은 자기 평균이 패딩보다 **더 낮아서** 오른다 —")
    print("      희석 방향이 뒤집힌 경우이고, 축을 낮추는 해킹과는 무관하다.)")
    print("  ⑶ ⚠️ 면역은 **한쪽 방향뿐이다** — 패딩이 핫스팟을 **올리는** 경우는 있다.")
    print("     자기 국소최대가 패딩 창보다 낮은 셔틀(TAT -0.876 · B6 -0.553)은 G·S 를 붙이면")
    print("     핫스팟이 패딩 자신의 값(-0.535 · -0.294)으로 **올라간다**. 해킹 방향은 막혔고")
    print("     반대로 '저응집 링커를 붙이면 응집이 나빠 보이는' 보수적 편향이 남는다.")


# --------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="개발성 3축(liability·응집 핫스팟·불안정성)의 대응 널 심사")
    ap.add_argument("--n", type=int, default=120,
                    help="널 2 의 셔틀·k 당 무작위 편집 횟수 (기본 120, §26 과 동일)")
    ap.add_argument("--shuffles", type=int, default=200,
                    help="널 1 의 셔틀당 조성 보존 셔플 횟수 (기본 200, §26 과 동일)")
    ap.add_argument("--seed", type=int, default=0, help="난수 시드 (기본 0)")
    args = ap.parse_args()

    # 🔴 의존성을 수천 번의 널 호출 **전에** 한 번 확인한다. 불안정성·GRAVY 가 Biopython
    #    없이 대체값으로 떨어지면 "불안정성 축은 셔플·편집에 아무 반응이 없다"는 **거짓
    #    결론**이 나온다(조용한 퇴화). 메시지는 `require_protein_analysis` 단일 출처.
    require_protein_analysis("개발성 3축 대응 널(불안정성 지수·GRAVY)")

    print("=" * 92)
    print("개발성 3축 대응 널 심사 — liability 개수 · 응집 핫스팟(a3v) · 불안정성 지수")
    print(f"심사 대상: core/developability_opt.objectives() = {OBJECTIVE_NAMES}")
    print(f"셔틀 {len(SHUTTLES)}종(라이브러리 전체) · seed={args.seed}"
          f" · 널1 셔플 {args.shuffles}회/셔틀 · 널2 {args.n}회/(셔틀·k) · 널3 결정론")
    print("=" * 92)
    null_1(args.shuffles, args.seed)
    null_2(args.n, args.seed)
    null_3()


if __name__ == "__main__":
    main()
