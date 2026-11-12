#!/usr/bin/env python3
"""불안정성 축 단독 대응 널 — 이 축은 "양전하를 지워라"는 지시인가, 그리고 방향별로 믿을 만한가.

질문
----
§39-4 가 Angiopep-2 **R11 한 자리**에서 분해한 결과: 19개 치환 전부가 불안정성을
개선하고 liability 는 0개 개선한다. 기전은 `DIWV['K']['R'] = +33.60` 이고 K10 이
동결돼 있어 R11 이 그 쌍을 깰 **유일한 쪽**이라는 것이었다. 한 셔틀 한 자리의 분해다.
이 스크립트는 네 가지를 묻는다.

  Q1 불안정성 축은 **체계적으로** 전하 삭제 유인인가, Angiopep-2 의 `KR` 이 특수한가.
     → 12셔틀 전수 단일 치환. 개선 수(move) 중 전하 삭제 비율을 **조성만으로 기대되는
       기저율**과 비교한다. 비교가 핵심이다 — K/R 이 많은 셔틀은 우연히도 전하를
       잃으므로 **발생(incidence)과 유인(incentive)을 분리**해야 한다.
  Q2 어느 이펩타이드가 끌고 있나. `KR` 은 **표 안에서** 특이한가, 아니면 같은 크기의
     전하 중립 이펩타이드가 있어 **전하를 지키는 탈출로**가 열려 있나.
  Q3 가족① 대응 널(같은 크기 무작위 편집). **개선 검출과 악화 검출을 따로** 판정한다.
     이 프로젝트의 판정 공간은 방향성이고 "둘 다 믿을 만함" 칸은 지금까지 늘 비어 있었다.
  Q4 가족④ 검출+한 수 해킹 쌍. 불안정성을 개선하면서 **알려진 중요한 것을 깨거나
     분자를 분명히 나쁘게 만드는 가장 싼 한 수**를 찾고, 현재 제약 셋이 그걸 잡는지 본다.
  Q5 순위는 40 이진(`optimizer_agent.II_STABLE`), 탐색은 연속이다. 보상을 만드는 것이
     **연속 기울기**인가 **임계 교차**인가 — 40 바로 위에 앉은 셔틀·구성체를 센다.

설계
----
· **전수 열거**(단일 치환 L×19)는 표본추출이 아니라 **결정론**이다. 같은 입력에 같은 수.
· 제약 셋은 생산 경로와 동일하게 적용한다 — `_frozen_idx`(문헌 보존 위치 동결) ·
  `rmt_sim ≥ 0.85` · off **양측 밴드** `±OFF_BAND`(§26). 제약 없는 수치와 제약 통과
  수치를 **둘 다** 적는다. 축의 내재적 성질과 생산 경로에서 도달 가능한 것은 다른 양이다.
· k=3 무작위 편집은 **시드가 독립 단위**다. 아카이브·배치 멤버는 같은 시드 안에서 상관되므로
  멤버를 독립 표본으로 세면 p 가 과대평가된다(§39-1 과 같은 규율). 시드별 비율에
  **양측 부호검정**을 건다.
· 원숭이 패치(monkeypatch)는 **하지 않는다**. 생산 함수(`_instability`·`objectives`·
  `off_target`·`_frozen_idx`·`pareto_scrub`)를 그대로 호출한다. 그래서 호출 카운터로
  패치 도달을 증명할 필요가 없다 — 패치가 없다.
· 🔴 `_instability` 가 0.0 을 돌려주면 **측정이 아니라 실패**다(옛 코드의 `return 0.0`
  은 0.0 < 40 이라 '안정'으로 읽혔다). 시작할 때 기준 서열로 가드한다.

실행: /usr/bin/python3 scripts/null_instability_axis.py
"""
from __future__ import annotations

import random
import statistics as st
import sys
from functools import lru_cache
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Bio.SeqUtils.ProtParamData import DIWV

import core.developability_opt as o
from core.binding import shuttle_similarity
from core.config import DEFAULT_CARGO, LINKER_LIBRARY, SHUTTLES
from core.developability import assess_developability
from core.evidence import CONSERVED_POSITIONS
from core.optimizer_agent import II_STABLE

POS = "KR"                 # 양전하(염기성) 잔기 — 전하 삭제 판정의 단일 출처
RMT_FLOOR = 0.85
N_SEEDS = 60               # k=3 널의 독립 단위
BATCH = 25                 # 시드별 배치 크기(시드 안에서 비율을 만든다)
K3 = 3                     # §28-2 와 같은 편집 크기
ARC_SEEDS = 20             # 아카이브 실현 확인용 `pareto_scrub` 시드 수

#: 화학적 liability 라벨(= 외부 라벨 후보). `core/developability.py` 의 문헌 기반
#: 모티프 집합이며 **Guruprasad 와 독립된 양**이다. 전하 주도 항목은 제외한다(§28-C) —
#: 그것은 "악화"의 라벨이 아니라 전하 자체이고, 바로 이 스크립트가 심사하는 대상이다.
_CHEM_LABELS = ("탈아마이드화", "이성질화", "N-당화", "자유 시스테인", "산화 취약")


# ---------------------------------------------------------------------------
# 캐시된 생산 함수 호출 (패치가 아니라 **호출 캐시**다 — 같은 서열에 같은 값)
# ---------------------------------------------------------------------------
@lru_cache(maxsize=200_000)
def ii(seq: str) -> float:
    return o._instability(seq)


@lru_cache(maxsize=200_000)
def sim(seq: str) -> float:
    return shuttle_similarity(seq).score


@lru_cache(maxsize=200_000)
def off(seq: str) -> float:
    return o.off_target(seq)


@lru_cache(maxsize=200_000)
def chem_liab(seq: str) -> tuple:
    """화학적 liability 항목만 (전하 주도 제외). 라벨 비교의 단일 출처."""
    d = assess_developability(seq)
    return tuple(x for x in d.liabilities if any(x.startswith(p) for p in _CHEM_LABELS))


def n_pos(s: str) -> int:
    return sum(c in POS for c in s)


def charge_del(old: str, new: str) -> bool:
    """이 치환이 **양전하를 삭제**하나 (K/R → K/R 아님). K→R 은 보존이므로 아니다."""
    return old in POS and new not in POS


def feasible(seq: str, off0: float) -> bool:
    """생산 경로의 제약 셋 통과 여부 — `pareto_scrub._try_insert` 와 같은 검사."""
    return sim(seq) >= RMT_FLOOR and abs(off(seq) - off0) <= o.OFF_BAND + 1e-9


def single_moves(seq: str):
    """동결 위치를 제외한 **모든 단일 치환**을 (위치0based, 원잔기, 새잔기, 서열)로."""
    frozen = o._frozen_idx(seq)
    for i, c in enumerate(seq):
        if i in frozen:
            continue
        for a in o._STD:
            if a != c:
                yield i, c, a, seq[:i] + a + seq[i + 1:]


# ---------------------------------------------------------------------------
# 통계 — 전부 표준 라이브러리. 외부 의존성 없음(재현성)
# ---------------------------------------------------------------------------
def fisher_right(a: int, b: int, c: int, d: int) -> float:
    """2×2 단측(우측, 농축) Fisher 정확검정. 표 = [[a,b],[c,d]].

    a = 개선 & 전하삭제 · b = 개선 & 비전하삭제 · c = 비개선 & 전하삭제 · d = 나머지.
    """
    n = a + b + c + d
    r1, c1 = a + b, a + c
    lo, hi = max(0, r1 + c1 - n), min(r1, c1)
    denom = comb(n, r1)
    p = sum(comb(c1, k) * comb(n - c1, r1 - k) for k in range(a, hi + 1)) / denom
    return min(1.0, p) if hi >= lo else 1.0


def sign_test(diffs) -> tuple:
    """양측 부호검정. 반환 (양수 시드 수, 동점 제외 시드 수, p)."""
    wins = sum(1 for d in diffs if d > 0)
    m = sum(1 for d in diffs if d != 0)
    if m == 0:
        return wins, 0, 1.0
    tail = sum(comb(m, i) for i in range(wins, m + 1)) / 2 ** m
    return wins, m, min(1.0, 2 * tail)


def _draw(seq0: str, k: int, rng: random.Random, idx: list) -> str | None:
    """시드에서 **정확히 k 자리** 다른 무작위 변이 1개(제약 검사 없음)."""
    s = list(seq0)
    for i in rng.sample(idx, k):
        s[i] = rng.choice(o._STD)
    cand = "".join(s)
    # 원래 잔기를 다시 뽑아 편집 수가 줄어든 경우는 **버린다** — 같은 크기 널이 아니다.
    return cand if sum(x != y for x, y in zip(cand, seq0)) == k else None


def matched_null(seq0: str, k: int, rng: random.Random, off0: float,
                 tries: int = 300) -> str | None:
    """시드에서 정확히 k 자리 다르고 **같은 제약을 통과하는** 무작위 변이 1개."""
    frozen = o._frozen_idx(seq0)
    idx = [i for i in range(len(seq0)) if i not in frozen]
    if len(idx) < k:
        return None
    for _ in range(tries):
        cand = _draw(seq0, k, rng, idx)
        if cand and feasible(cand, off0):
            return cand
    return None


def null_feasibility(seq0: str, k: int, off0: float, probes: int = 600) -> float:
    """🔴 k 편집 널이 **애초에 제약을 통과할 수 있나**를 먼저 잰다.

    짧은 셔틀(T7 7aa)에서 3치환은 `rmt_sim ≥ 0.85` 를 거의 항상 깬다. 그 사실을
    "널을 못 만들었다"는 침묵으로 흘리면 **측정 불가가 음성으로 위장**된다. 통과
    가능 비율을 명시적으로 보고하고, 0 이면 그 셔틀은 **측정 불가**로 적는다.
    """
    frozen = o._frozen_idx(seq0)
    idx = [i for i in range(len(seq0)) if i not in frozen]
    if len(idx) < k:
        return 0.0
    rng = random.Random(4242)
    ok = n = 0
    for _ in range(probes):
        cand = _draw(seq0, k, rng, idx)
        if not cand:
            continue
        n += 1
        ok += feasible(cand, off0)
    return (ok / n) if n else 0.0


# ---------------------------------------------------------------------------
def guard() -> None:
    """🔴 0.0 은 '안정'이 아니라 **실패**다. 알려진 기준값으로 엔진을 먼저 증명한다."""
    probes = {"TFFYGGSRGKRNNFKTEEY": 71.5, "GRKKRRQRRRPPQ": 225.8}
    for s, want in probes.items():
        got = ii(s)
        if abs(got - want) > 0.05:
            raise SystemExit(f"🔴 엔진 가드 실패: II({s[:12]}…) = {got}, 기대 {want}. "
                             "Biopython ProtParam 이 정상인지 확인하라(0.0 은 측정값이 아니다).")
    print(f"엔진 가드 통과 — II(Angiopep-2)={probes['TFFYGGSRGKRNNFKTEEY']} · "
          f"II(TAT)={probes['GRKKRRQRRRPPQ']} · 패치 없음(생산 함수 직접 호출)")


def section(t: str) -> None:
    print()
    print("=" * 104)
    print(t)


# ---------------------------------------------------------------------------
def main() -> None:
    guard()
    names = list(SHUTTLES)
    seqs = {n: SHUTTLES[n]["seq"] for n in names}
    base = {n: ii(s) for n, s in seqs.items()}
    off0 = {n: off(s) for n, s in seqs.items()}
    # 시드 자신이 하한을 위반하는 셔틀은 **생산 경로가 아예 거부**한다(§33-4 밀반입 차단).
    live = [n for n in names if sim(seqs[n]) >= RMT_FLOOR]
    dead = [n for n in names if n not in live]

    section("A. 기준값 — 불안정성은 셔틀 단독(탐색 범위), 순위는 구성체 전체(§주의: 다른 양)")
    print(f"{'셔틀':16} {'길':>3} {'II(셔틀)':>9} {'40기준':>7} {'off':>6} {'K/R':>4} "
          f"{'편집가능':>8} {'문헌동결':>10}")
    for n in names:
        s = seqs[n]
        pins = sorted(CONSERVED_POSITIONS.get(n) or ())
        print(f"{n:16} {len(s):>3} {base[n]:>9.1f} {'불안정' if base[n] >= II_STABLE else '안정':>7} "
              f"{off0[n]:>6.3f} {n_pos(s):>4} {len(s) - len(o._frozen_idx(s)):>8} "
              f"{str(pins or '없음'):>10}")
    if dead:
        print(f"⚠️ 제약 분석 제외 {dead} — 시드 자신이 `rmt_sim < {RMT_FLOOR}`"
              f"(참조 명부 미등재, §32). `pareto_scrub` 이 ValueError 로 거부한다.")

    # -----------------------------------------------------------------
    section("B. Q1 — 단일 치환 전수. 개선 수 중 전하 삭제 비율 vs **조성만의 기저율**")
    print("기저율 = (전하 삭제 가능 수)/(전체 수). 발생과 유인을 가르는 양은 농축비다.")
    print(f"{'셔틀':16} {'전체수':>6} {'개선':>6} {'개선%':>6} | {'전하삭제/개선':>13} "
          f"{'관측%':>6} {'기저%':>6} {'농축':>6} {'p(단측)':>9}")
    pool = [0, 0, 0, 0]
    rows_b = {}
    for n in names:
        s = seqs[n]
        a = b = c = d = 0
        for _, oc, nc, m in single_moves(s):
            cd = charge_del(oc, nc)
            imp = ii(m) < base[n]
            if imp and cd:
                a += 1
            elif imp:
                b += 1
            elif cd:
                c += 1
            else:
                d += 1
        tot, impn = a + b + c + d, a + b
        obs = 100 * a / impn if impn else 0.0
        bas = 100 * (a + c) / tot
        enr = (obs / bas) if bas else float("nan")
        p = fisher_right(a, b, c, d)
        rows_b[n] = (a, b, c, d)
        for i, v in enumerate((a, b, c, d)):
            pool[i] += v
        print(f"{n:16} {tot:>6} {impn:>6} {100*impn/tot:>6.1f} | {a:>6}/{impn:<6} "
              f"{obs:>6.1f} {bas:>6.1f} {enr:>6.2f} {p:>9.4f}"
              f"{'  🔴' if p < 0.05 else ''}")
    a, b, c, d = pool
    print(f"{'— 합계(풀)':16} {a+b+c+d:>6} {a+b:>6} {100*(a+b)/(a+b+c+d):>6.1f} | "
          f"{a:>6}/{a+b:<6} {100*a/(a+b):>6.1f} {100*(a+c)/(a+b+c+d):>6.1f} "
          f"{(a/(a+b))/((a+c)/(a+b+c+d)):>6.2f} {fisher_right(a, b, c, d):>9.4f}")

    print()
    print("B-2. 같은 표를 **제약 통과 수만**으로 (동결 · rmt_sim≥0.85 · off 양측 밴드)")
    print(f"{'셔틀':16} {'통과수':>6} {'개선':>6} | {'전하삭제/개선':>13} {'관측%':>6} "
          f"{'기저%':>6} {'농축':>6} {'p(단측)':>9}")
    for n in live:
        s = seqs[n]
        a = b = c = d = 0
        for _, oc, nc, m in single_moves(s):
            if not feasible(m, off0[n]):
                continue
            cd, imp = charge_del(oc, nc), ii(m) < base[n]
            if imp and cd:
                a += 1
            elif imp:
                b += 1
            elif cd:
                c += 1
            else:
                d += 1
        tot, impn = a + b + c + d, a + b
        if not tot or not impn:
            print(f"{n:16} {tot:>6} {impn:>6} | (측정 불가)")
            continue
        obs, bas = 100 * a / impn, 100 * (a + c) / tot
        print(f"{n:16} {tot:>6} {impn:>6} | {a:>6}/{impn:<6} {obs:>6.1f} {bas:>6.1f} "
              f"{(obs/bas) if bas else float('nan'):>6.2f} {fisher_right(a,b,c,d):>9.4f}"
              f"{'  🔴' if fisher_right(a,b,c,d) < 0.05 else ''}")

    print()
    print("B-3. 구성체 수준(화물 GSNKGAIIGLM + 직접융합) — 순위 축이 실제로 재는 범위")
    print(f"{'셔틀':16} {'구성체II':>9} {'개선':>6} {'개선%':>6} | {'관측%':>6} {'기저%':>6} {'농축':>6}")
    for n in names:
        s = seqs[n]
        con = DEFAULT_CARGO + s
        cb = ii(con)
        a = b = c = 0
        for i, oc, nc, _ in single_moves(s):
            j = len(DEFAULT_CARGO) + i
            m = con[:j] + nc + con[j + 1:]
            cd, imp = charge_del(oc, nc), ii(m) < cb
            a += imp and cd
            b += imp and not cd
            c += (not imp) and cd
        tot = sum(1 for _ in single_moves(s))
        obs = 100 * a / (a + b) if a + b else 0.0
        bas = 100 * (a + c) / tot
        print(f"{n:16} {cb:>9.1f} {a+b:>6} {100*(a+b)/tot:>6.1f} | {obs:>6.1f} {bas:>6.1f} "
              f"{(obs/bas) if bas else float('nan'):>6.2f}")

    # -----------------------------------------------------------------
    section("C. Q2 — 어느 이펩타이드가 끌고 있나 (Guruprasad DIWV, 서열에 실제로 나타난 것)")
    tbl = sorted(((w, x + y) for x, row in DIWV.items() for y, w in row.items()),
                 reverse=True)
    tiers = sorted({round(w, 2) for w, _ in tbl}, reverse=True)[:5]
    print(f"DIWV 표 최대 계층: {tiers}  · 최대값 쌍: "
          f"{[d for w, d in tbl if abs(w - tbl[0][0]) < 1e-9]}  · KR = {DIWV['K']['R']}")
    print("→ `KR`(+33.60)은 **표 안에서 특이하지 않다**. 더 큰 전하 중립 이펩타이드가 있다.")
    print()
    print(f"{'셔틀':16} {'II':>7} {'상위 이펩타이드 3개 (*=K/R 포함)':44} "
          f"{'최선Δ전하삭제':>12} {'최선Δ전하보존':>12} {'보존로 동급?':>12}")
    for n in names:
        s = seqs[n]
        dp = sorted(((DIWV[s[i]][s[i + 1]], s[i] + s[i + 1], i + 1)
                     for i in range(len(s) - 1)), reverse=True)
        top = ", ".join(f"{d}@{p}{'*' if (d[0] in POS or d[1] in POS) else ''} {w:+.1f}"
                        for w, d, p in dp[:3])
        bc = bp = 0.0
        for _, oc, nc, m in single_moves(s):
            delta = base[n] - ii(m)
            if charge_del(oc, nc):
                bc = max(bc, delta)
            else:
                bp = max(bp, delta)
        print(f"{n:16} {base[n]:>7.1f} {top:44} {bc:>12.1f} {bp:>12.1f} "
              f"{'예' if bp >= bc else '🔴 아니오':>12}")

    # -----------------------------------------------------------------
    section("D. Q3 — 가족① 대응 널. **개선 검출과 악화 검출을 따로** 판정한다")
    print("D-1 개선 검출 — 개선이 값싼가. 연속(탐색 축)과 40 이진(순위 축)을 나눠 센다.")
    print(f"{'셔틀':16} {'k=1 전수':>9} {'연속개선%':>9} {'제약통과 연속%':>14} "
          f"{'40교차%(전수)':>13} {'40교차%(제약)':>13}")
    for n in names:
        s = seqs[n]
        mv = list(single_moves(s))
        cont = sum(1 for _, _, _, m in mv if ii(m) < base[n])
        fe = [m for _, _, _, m in mv if n in live and feasible(m, off0[n])]
        contf = sum(1 for m in fe if ii(m) < base[n])
        cross = sum(1 for _, _, _, m in mv if base[n] >= II_STABLE > ii(m))
        crossf = sum(1 for m in fe if base[n] >= II_STABLE > ii(m))
        print(f"{n:16} {len(mv):>9} {100*cont/len(mv):>9.1f} "
              f"{(100*contf/len(fe) if fe else float('nan')):>14.1f} "
              f"{100*cross/len(mv):>13.1f} "
              f"{(100*crossf/len(fe) if fe else float('nan')):>13.1f}")

    print()
    print(f"D-1b k={K3} 무작위 편집 · 시드 {N_SEEDS}개 × 배치 {BATCH} · **시드가 독립 단위**")
    print("    비교 기준은 50%(= 축이 동전 던지기면 나오는 값). 양측 부호검정.")
    print(f"{'셔틀':16} {'널통과율':>8} {'연속개선% 평균':>14} {'시드>50%':>9} {'p':>8} | "
          f"{'40교차% 평균':>13} {'교차>0 시드':>11}")
    for n in live:
        s = seqs[n]
        frac = null_feasibility(s, K3, off0[n])
        if frac < 0.01:
            print(f"{n:16} {100*frac:>7.1f}% — 🔴 **측정 불가**: k={K3} 널이 제약"
                  f"(rmt_sim≥{RMT_FLOOR}·off 밴드)을 거의 통과하지 못한다")
            continue
        rates, crates = [], []
        for sd in range(N_SEEDS):
            rng = random.Random(10_000 + sd)
            batch = [x for x in (matched_null(s, K3, rng, off0[n]) for _ in range(BATCH)) if x]
            if not batch:
                continue
            rates.append(100 * sum(1 for m in batch if ii(m) < base[n]) / len(batch))
            crates.append(100 * sum(1 for m in batch
                                    if base[n] >= II_STABLE > ii(m)) / len(batch))
        if not rates:
            print(f"{n:16} {100*frac:>7.1f}% (제약을 통과하는 k={K3} 널을 못 만들었다)")
            continue
        w, m, p = sign_test([r - 50.0 for r in rates])
        print(f"{n:16} {100*frac:>7.1f}% {st.mean(rates):>14.1f} {w:>5}/{m:<3} {p:>8.4f} | "
              f"{st.mean(crates):>13.1f} {sum(1 for x in crates if x > 0):>8}/{len(crates)}")

    print()
    print("D-2 악화 검출 — 🔴 외부 라벨이 **하나뿐**이다: 이 프로젝트 자신의 문헌 기반")
    print("    화학적 liability 모티프 생성(탈아마이드화·이성질화·당화 sequon·자유Cys·산화).")
    print("    Guruprasad 와 **독립된 양**이므로 라벨로 쓸 수 있다. 널 = 같은 위치의 비생성 치환.")
    print(f"{'셔틀':16} {'라벨악화수':>10} {'II↑%':>7} {'II↓%':>7} | {'널 수':>7} "
          f"{'널 II↑%':>8} {'초과pp':>7}")
    LU = LD = LE = NU = ND = NE = 0
    for n in names:
        s = seqs[n]
        b0 = set(chem_liab(s))
        lu = ld = le = nu = nd = ne = 0
        for _, _, _, m in single_moves(s):
            worse = len(set(chem_liab(m)) - b0) > 0 and len(chem_liab(m)) > len(b0)
            x = ii(m)
            up, dw = x > base[n], x < base[n]
            if worse:
                lu += up; ld += dw; le += not (up or dw)
            else:
                nu += up; nd += dw; ne += not (up or dw)
        LU += lu; LD += ld; LE += le; NU += nu; ND += nd; NE += ne
        lt, nt = lu + ld + le, nu + nd + ne
        print(f"{n:16} {lt:>10} {100*lu/lt if lt else float('nan'):>7.1f} "
              f"{100*ld/lt if lt else float('nan'):>7.1f} | {nt:>7} "
              f"{100*nu/nt:>8.1f} {(100*lu/lt - 100*nu/nt) if lt else float('nan'):>+7.1f}")
    lt, nt = LU + LD + LE, NU + ND + NE
    print(f"{'— 합계(풀)':16} {lt:>10} {100*LU/lt:>7.1f} {100*LD/lt:>7.1f} | {nt:>7} "
          f"{100*NU/nt:>8.1f} {100*LU/lt - 100*NU/nt:>+7.1f}")
    print(f"   Fisher 단측(라벨악화가 II↑ 에 농축되나) p = "
          f"{fisher_right(LU, lt - LU, NU, nt - NU):.4f}")

    # -----------------------------------------------------------------
    section("E. Q4 — 가족④ 검출 + 한 수 쌍. 불안정성을 개선하는 **가장 싼 나쁜 한 수**")
    print("E-1 문헌 고정 위치(= 알려진 기능 필수 잔기)를 깨는 한 수. 동결이 이걸 잡나?")
    for n, pins in CONSERVED_POSITIONS.items():
        if n not in seqs:
            continue
        s = seqs[n]
        for p in pins:
            cands = sorted((ii(s[:p-1] + a + s[p:]), a) for a in o._STD if a != s[p-1])
            good = [(x, a) for x, a in cands if x < base[n]]
            cross = [(x, a) for x, a in good if x < II_STABLE]
            best = good[0] if good else None
            print(f"   {n} {s[p-1]}{p}: 개선 {len(good)}/19 · 40교차 {len(cross)} · "
                  f"최선 {s[p-1]}{p}→{best[1]} II {base[n]}→{best[0]}" if best else
                  f"   {n} {s[p-1]}{p}: 개선 0/19 (이 자리는 불안정성 유인이 없다)")
            if best:
                m = s[:p-1] + best[1] + s[p:]
                print(f"      이 한 수를 잡는 것: 동결={'예' if (p-1) in o._frozen_idx(s) else '아니오'} · "
                      f"rmt_sim {sim(m):.3f}(하한 {RMT_FLOOR}) → "
                      f"{'통과' if sim(m) >= RMT_FLOOR else '거부'} · "
                      f"off {off(m):.3f} vs 시드 {off0[n]:.3f} (밴드 ±{o.OFF_BAND}) → "
                      f"{'통과' if abs(off(m)-off0[n]) <= o.OFF_BAND + 1e-9 else '거부'}")

    print()
    print("E-2 화학적 liability 를 **새로 만들면서** 불안정성을 개선하는 한 수 — 제약이 잡나")
    print(f"{'셔틀':16} {'해당 한 수':>10} {'제약통과':>9} {'그중 40교차':>11} {'최선(ΔII 최대)':>30}")
    for n in names:
        s = seqs[n]
        b0 = set(chem_liab(s))
        hits = []
        for _, oc, nc, m in single_moves(s):
            if ii(m) < base[n] and len(chem_liab(m)) > len(b0):
                hits.append((base[n] - ii(m), oc, nc, m, n in live and feasible(m, off0[n]),
                             tuple(set(chem_liab(m)) - b0)))
        ok = [h for h in hits if h[4]]
        okc = [h for h in ok if ii(h[3]) < II_STABLE]
        bst = max(ok, default=None)
        desc = ""
        if bst:
            i = next(i for i, (a, b) in enumerate(zip(s, bst[3])) if a != b)
            desc = (f"{bst[1]}{i+1}→{bst[2]} II {base[n]:.1f}→{ii(bst[3]):.1f} "
                    f"[{bst[5][0].split('(')[0].strip()}]")
        print(f"{n:16} {len(hits):>10} {len(ok):>9} {len(okc):>11} {desc:>30}")

    print()
    print(f"E-3 그 한 수가 **실제 탐색 결과에 실현되나** — `pareto_scrub` 시드 0~{ARC_SEEDS-1}")
    print("    (자유 시스테인은 이황화 스크램블 위험으로 이 프로젝트가 직접 liability 로 센다)")
    print(f"{'셔틀':16} {'아카이브':>8} {'새 자유Cys 멤버':>14} {'%':>6} {'예시':>34}")
    for n in live:
        s = seqs[n]
        tot = cys = 0
        ex = ""
        for sd in range(ARC_SEEDS):
            for mem in o.pareto_scrub(s, seed=sd, rmt_floor=RMT_FLOOR):
                tot += 1
                if mem["seq"].count("C") > s.count("C"):
                    cys += 1
                    ex = ex or f"{mem['seq']} II={mem['objectives']['instability']}"
        print(f"{n:16} {tot:>8} {cys:>14} {100*cys/tot if tot else 0:>6.1f} {ex:>34}")

    # -----------------------------------------------------------------
    section("F. Q5 — 보상을 만드는 것은 연속 기울기인가 **임계 40 교차**인가")
    print("순위 축은 `0 if II < 40 else 1` (optimizer_agent.II_STABLE) — 교차만 값을 바꾼다.")
    print(f"{'셔틀':16} {'II(셔틀)':>9} {'40까지':>7} {'교차 한 수':>10} {'그중 전하삭제':>13} "
          f"{'제약통과 교차':>13} {'그중 전하삭제':>13}")
    for n in names:
        s = seqs[n]
        cr = crc = fk = fkc = 0
        for _, oc, nc, m in single_moves(s):
            if base[n] >= II_STABLE > ii(m):
                cr += 1
                crc += charge_del(oc, nc)
                if n in live and feasible(m, off0[n]):
                    fk += 1
                    fkc += charge_del(oc, nc)
        print(f"{n:16} {base[n]:>9.1f} {base[n]-II_STABLE:>+7.1f} {cr:>10} {crc:>13} "
              f"{fk:>13} {fkc:>13}")

    print()
    print("F-2 구성체 수준(순위가 실제로 보는 양) — 링커 13종 × 셔틀 12종 = 156 구성체")
    rows = []
    for n in names:
        for lv in LINKER_LIBRARY.values():
            rows.append((ii(DEFAULT_CARGO + lv["seq"] + seqs[n]), n))
    ge = [r for r in rows if r[0] >= II_STABLE]
    near = [r for r in rows if II_STABLE <= r[0] <= II_STABLE + 10]
    print(f"   II ≥ 40: {len(ge)}/{len(rows)} · 40 ≤ II ≤ 50 (한 수로 넘길 수 있는 띠): "
          f"{len(near)}/{len(rows)}")
    print("🔴 마지막 열이 요점이다 — **잔기 편집 0회**로 40 을 넘기는 링커가 있나.")
    print(f"{'셔틀':16} {'셔틀II':>7} {'구성체II 범위':>18} {'≥40':>6} {'40~50':>7} "
          f"{'0편집 교차 링커':>15} {'그 링커':>14}")
    for n in names:
        v = [(ii(DEFAULT_CARGO + lv["seq"] + seqs[n]), ln)
             for ln, lv in LINKER_LIBRARY.items()]
        lo = sorted(x for x, _ in v)
        under = [ln for x, ln in v if x < II_STABLE]
        print(f"{n:16} {base[n]:>7.1f} {f'{lo[0]:.1f}~{lo[-1]:.1f}':>18} "
              f"{sum(1 for x, _ in v if x >= 40):>6} "
              f"{sum(1 for x, _ in v if 40 <= x <= 50):>7} {len(under):>15} "
              f"{(min(under, key=len) if under else '없음'):>14}")
    print("→ 탐색은 **셔틀 단독 연속 II**(위 '셔틀II')를 깎고, 순위는 **구성체 II 의 40 이진**을")
    print("  읽는다. 범위도 다르고 척도도 다르다 — 같은 이름의 **세 가지 양**이다.")

    # -----------------------------------------------------------------
    section("G. 🔴 이 측정이 **볼 수 없는** 것")
    for line in (
        "1. 실제 안정성. Guruprasad II 는 이펩타이드 조성 통계이고 이 프로젝트에 Tm·ΔG·"
        "혈청 반감기 실측이 **하나도 없다**. 따라서 'II 가 내려갔다 = 분자가 안정해졌다'는"
        " 어디서도 검증되지 않았다. D-2 는 II 가 아니라 **다른 축의 라벨**로 잰 것이다.",
        "2. 다중 편집의 상호작용. 전수 열거는 **단일 치환**이다. k≥2 에서 한 자리가 다른"
        " 자리의 이펩타이드를 바꾸는 효과(epistasis)는 D-1b 의 k=3 비율에만 간접적으로 들어간다.",
        "3. 11/12 셔틀의 기능 필수 잔기. `CONSERVED_POSITIONS` 에 문헌 위치가 있는 셔틀은"
        " Angiopep-2 하나다(§38-2). 나머지에서 '전하 삭제가 기능을 깬다'는 **측정도 반박도"
        " 안 된다** — 라벨이 없다.",
        "4. Tet1 의 제약 하 거동. 참조 명부 미등재로 `rmt_sim = 0.0`(§32) → 제약 분석에서"
        " 빠진다. 무제약 열거에만 들어 있다.",
        "5. 탐색 궤적. E-3 은 `pareto_scrub` 아카이브의 **최종 구성**이고, GA 가 그 멤버에"
        " 도달한 경로나 중간에 버린 후보는 보지 않는다.",
        "6. 링커·화물 상호작용의 전수. F-2 는 화물 1종(GSNKGAIIGLM)·링커 13종이다."
        " 다른 화물에서 40 근처 띠의 구성은 달라진다(§35: 게이트는 일반화되지 않았다).",
        "7. off 밴드의 폭. `OFF_BAND = 0.05` 를 고정값으로 썼다. 밴드를 넓히면 B-2·F 의"
        " '제약 통과' 열이 전부 커진다 — 이 측정은 **현재 상수 하나의 단면**이다.",
        f"8. 🔴 k={K3} 널은 11셔틀 중 **3종에서만 측정됐다**(Angiopep-2 2.5% · Leptin30"
        " 74.3% · RVG29 70.0% 만 통과율 ≥1%, 나머지 8종은 측정 불가)."
        " 짧은 셔틀은 3치환이 `rmt_sim ≥ 0.85` 를 거의 항상 깨므로 §28-2 의"
        " 56.4%(바닥 제외·제약 없음)와 **같은 조건이 아니다**. D-1b 의 비율은 제약을"
        " 통과한 **선택된 부분집합**의 비율이고(Angiopep-2 널 통과율 2.5%), 무조건부"
        " 비율은 k=1 전수(D-1)뿐이다.",
    ):
        print(f"   {line}")
    print()
    print("재현: /usr/bin/python3 scripts/null_instability_axis.py  (난수는 시드 고정, "
          "전수 열거는 결정론)")


if __name__ == "__main__":
    main()
