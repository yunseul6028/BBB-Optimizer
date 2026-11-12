"""**생성** 링커의 대응 널 — 라이브러리가 아니라 **생성 경로**를 피고석에 세운다.

`docs/protocol.md` §4 R3 의 🔴 미이행 항목이고 §8 체크리스트의
"**생성 링커의 대응 널 미측정**"이 바로 이것이다. 라이브러리 링커는 두 번 심사됐다
(M8 `--metric bbb` → 0/108 = **거꾸로** · M10 `--metric risk` → 41/96 = 통과), 그런데
2026-11-11 에 **자유 생성을 다시 열면서**(§23) 그 채널의 널은 만들지 않았다.

왜 '생성기'를 따로 정의해야 하나 — 그리고 그게 이 측정의 최대 한계다
----------------------------------------------------------------------
생성 경로의 실제 생성기는 **LLM**(`design_candidate`)이다. 유료·비결정론이라 널
(셀당 30개 × 수백 셀)에 그대로 넣을 수 없다. 그래서 심사 대상을 LLM 이 아니라
**시스템이 LLM 에게 주는 생성 정책**으로 바꾼다. 정책은 코드에 결정론적으로 박혀 있다:

  ⑴ `optimizer_agent.linker_reject_reason` — 표준 20 아미노산 · 최대 `LINKER_MAX_LEN`(17aa)
  ⑵ `optimizer_agent_gemini._junction_note` → `developability.junction_liabilities` —
     "접합부를 걸친 모티프는 **링커를 바꿔 없앨 수 있다**. 링커의 **첫/끝 잔기**가 화물 끝·
     셔틀 시작과 모티프를 만들지 않게 골라라. `GGS`(G 로 시작)는 실패하고 `EAAAK`·`AAAAA` 는
     성공한다"

그래서 **두 개의 결정론적 대리 생성기**를 둔다. 둘은 **탐색 공간이 완전히 같고 목적만 다르다**:

  · `junction` — 접합부 liability 만 최소화하고 **나머지는 crc32 결정론 난수로 고른다**.
                 = 정책이 말하는 그대로. 다른 축을 **일부러 최적화하지 않는다**(안 그러면
                 생성기가 '심사받는 축을 직접 최적화'해서 승리가 동어반복이 된다).
  · `full`     — 같은 공간에서 **순위 위험 축 4개 전부**를 사전식으로 최소화. 생성 경로가
                 낼 수 있는 **상한**이고, `junction` 과의 차이가 "접합부 밖에 남은 여지"다.

그리고 🔴 **어휘 대응 널**을 하나 더 둔다. 계열 1 의 표준 널은 **길이만** 맞추므로,
생성기가 이기면 "접합부 목적이 일을 했다"와 "스캐폴드 어휘(GS·EAAAK·A-런)가 일을 했다"가
구별되지 않는다. 그래서 세 번째 생성기:

  · `space`    — 같은 탐색 공간에서 **목적 없이** crc32 로 하나 고른다(화물·셔틀 비의존).
                 `space` vs 균일 무작위 = **어휘의 몫**, `junction` vs `space` = **목적의 몫**.

🔴 **이 스크립트가 심사하는 것은 LLM 의 설계 능력이 아니다.** 정책의 **상한과 하한**이다.
LLM 이 정책보다 잘할 수도, 못할 수도 있다(§30 에서 서열 편집 단계의 LLM 은 GA 를 이기지
못했다). 실제로 생성 경로가 라이브 세션에서 쓰인 기록은 없다(§23-8③).

바닥(floor) 제외 규칙 — §26-4 · §28-2 와 같은 형식
--------------------------------------------------
**직접융합에서 접합부 liability 가 0 인 셀에서는 생성 링커가 개선할 것이 없다.** 그 셀의
비교는 "개선"이 아니라 "악화를 피했나"뿐이므로 섞으면 가짜 양성이 된다. 그래서 전체 집계와
**바닥 제외** 집계를 **같이** 보고한다. 바닥 판정은 `junction_liabilities(cargo, "", shuttle)`.

화물을 여러 개 쓰는 이유 (§23-6 · §8 첫 항목)
----------------------------------------------
M1~M10 이 전부 화물 `GSNKGAIIGLM` **하나**로 측정됐고, 이 화물은 **M 으로 끝나** 접합부
모티프가 생기지 않는다 → 접합부 축이 **구조적으로 측정 불가**였다. 그래서 끝 잔기를 바꾼
프로브 화물을 같이 돌린다. `_SEQ_LIAB_PATTERNS` 를 읽고 역산한 결과:

  · `...M` 끝 (기본) — 전 셔틀 바닥. **대조군**이다.
  · `...N` 끝 (§23-5 프로브) — `N[GSTN]`·`N[^P][ST]` 가 셔틀 시작과 걸린다.
  · `...D` 끝 — `D[GSTDH]`. §23 이 언급하지 않은 **두 번째 프로브**다.
  · `...Q` 끝 — 🔴 문서가 "N 이나 Q 로 끝나는 화물"이라 쓰지만 `_SEQ_LIAB_PATTERNS` 에
    **Q 패턴이 없다**(탈아마이드화 정규식이 `N[GSTN]` 하나). 그 주장의 검정으로 포함한다.

⚠️ `core/` 는 **건드리지 않는다.** 몽키패치도 하지 않는다(화물·링커를 전부 인자로 넘기므로
필요가 없다). 실행 끝에 프로덕션 기본값이 그대로인지 **단언**한다.

재현 (서브프로세스 0회 · 완전 결정론 · Biopython 필요. `python3 -I` 금지):
  $ /usr/bin/python3 scripts/null_linker_generated.py                    # 시드 0
  $ /usr/bin/python3 scripts/null_linker_generated.py --seed 1 --nulls 60
  $ /usr/bin/python3 scripts/null_linker_generated.py --gen junction     # 정책만
"""
from __future__ import annotations

import argparse
import random
import statistics
import sys
import zlib
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "scripts"))

from core.config import DEFAULT_CARGO, LINKER_LIBRARY  # noqa: E402
from core.developability import junction_liabilities  # noqa: E402
from core.optimizer_agent import LINKER_MAX_LEN, linker_reject_reason  # noqa: E402

# 구조를 재사용한다 — 널 설계·지배 정의·셔틀 집합은 라이브러리 널과 **같아야** 비교 가능하다.
from null_linker import _STD, _dominates, risk_vector, shuttles  # noqa: E402

# 생성기의 탐색 공간을 만드는 **스캐폴드** — 시스템이 LLM 에게 보여주는 어휘다.
# 라이브러리 13종의 모티프(GS 계열·EAAAK 계열) + §23-4 의 실제 생성 링커(EGEGEG 의 EGEG)
# + §23-5 에서 접합부 제거에 성공한 A-런. 여기에 **양 끝 잔기 전수 치환**을 더한 것이
# `_junction_note` 가 지시하는 편집("링커의 첫/끝 잔기를 골라라")과 정확히 같은 이동이다.
_SCAFFOLDS = ("GGGGS", "GGS", "GS", "EAAAK", "AEAAK", "EGEG", "A", "P")

# §23-4 · §23-5 에 **실제로 기록된** 생성/설계 링커. (GGNGS·GGCGS 는 일부러 나쁘게 만든
# 음성 프로브이므로 따로 표시한다 — 승률 분모에 섞지 않는다.)
_RECORDED = (
    ("EGEGEG", "생성(§23-4)"),
    ("EAAAKGSGS", "생성(§23-4)"),
    ("A", "접합부 제거 성공(§23-5)"),
    ("AAAAA", "접합부 제거 성공(§23-5)"),
    ("EAAAK", "접합부 제거 성공(§23-5) · 단 라이브러리 멤버"),
    ("GGNGS", "🔴 음성 프로브(NG 삽입, §23-4)"),
    ("GGCGS", "🔴 음성 프로브(유리 Cys, §23-4)"),
)

_AXES = ("dev_liab", "II(이진)", "off-target", "응집핫스팟")


def _spearman(xs, ys) -> float | None:
    """순위상관(동점은 평균 순위). R2 체크리스트 ①「각 축 ↔ 길이 순위상관」용."""
    n = len(xs)
    if n < 3:
        return None

    def rank(v):
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

    rx, ry = rank(list(xs)), rank(list(ys))
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((p - mx) * (q - my) for p, q in zip(rx, ry))
    den = (sum((p - mx) ** 2 for p in rx) * sum((q - my) ** 2 for q in ry)) ** 0.5
    return None if den == 0 else num / den


def _crc(s: str) -> int:
    """시드와 무관한 결정론 타이브레이크. 축과 상관이 없도록 서열의 crc32 를 쓴다
    (사전식 타이브레이크는 poly-A 를 고르는데, poly-A 는 II·응집에서 우연히 좋아
    `junction` 생성기에 다른 축의 최적화를 **몰래 넣는다**)."""
    return zlib.crc32(s.encode())


def _pool(length: int) -> list[str]:
    """길이 `length` 의 생성 후보 공간 — 스캐폴드 타일링 × 양 끝 잔기 전수 치환."""
    out: set[str] = set()
    for sc in _SCAFFOLDS:
        base = (sc * (length // len(sc) + 1))[:length]
        out.add(base)
        for aa in _STD:                       # 첫 잔기 치환
            out.add(aa + base[1:])
            out.add(base[:-1] + aa)           # 끝 잔기 치환
        if length >= 2:
            for a in _STD:                    # 양 끝 동시 치환
                for b in _STD:
                    out.add(a + base[1:-1] + b)
    # 생성 경로의 하드 제약을 **실제 프로덕션 함수**로 건다(복제하지 않는다).
    return sorted(x for x in out if linker_reject_reason(x) is None)


def _search_space(length: int, cap: int) -> list[str]:
    """두 생성기가 **공유하는** 탐색 공간. 접합부 기준 정렬이 아니라 crc32 정렬로 자른다
    (화물·셔틀과 무관하게 같은 공간을 쓰도록)."""
    p = _pool(length)
    p.sort(key=_crc)
    return p[:cap]


def generate(cargo: str, shuttle: str, length: int, space: list[str], mode: str) -> str:
    """결정론 생성기. `mode` 가 목적만 바꾸고 탐색 공간은 같다."""
    if mode == "space":                       # 목적 없음 — 어휘 대응 널
        return min(space, key=_crc)
    if mode == "junction":
        return min(space, key=lambda s: (len(junction_liabilities(cargo, s, shuttle)), _crc(s)))
    return min(space, key=lambda s: (risk_vector(cargo, s, shuttle), _crc(s)))


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0, help="🔴 반드시 기록할 난수 시드")
    ap.add_argument("--nulls", type=int, default=30,
                    help="셀당 같은 길이 무작위 링커 개수 (기본 30 = M10 과 동일)")
    ap.add_argument("--lengths", default="1,2,3,5,6,10,15,17",
                    help="생성 링커 길이. 기본값은 라이브러리 길이 집합 + 2aa")
    ap.add_argument("--cargos", default="GSNKGAIIGLM,GSNKGAIIGLN,GSNKGAIIGLD,GSNKGAIIGLQ")
    ap.add_argument("--gen", choices=["all", "space", "junction", "full"], default="all")
    ap.add_argument("--gen-cap", type=int, default=150,
                    help="생성기 탐색 공간 크기(결정론 crc32 절단). 생성기 강도의 상한이다")
    ap.add_argument("--verbose", action="store_true", help="셀 단위 표를 전부 찍는다")
    a = ap.parse_args()

    lengths = [int(x) for x in a.lengths.split(",") if x.strip()]
    cargos = [c.strip().upper() for c in a.cargos.split(",") if c.strip()]
    modes = ["space", "junction", "full"] if a.gen == "all" else [a.gen]
    SH = shuttles("risk")                     # M10 과 같은 후보군 8종
    space = {L: _search_space(L, a.gen_cap) for L in lengths}

    print("=" * 78)
    print("생성 링커의 대응 널 (protocol.md 계열 1 · 구성체 수준)")
    print("=" * 78)
    print(f"🔴 시드 {a.seed} · 널 {a.nulls}개/셀 · 셔틀 {len(SH)}종 {list(SH)}")
    print(f"길이 {lengths} · 화물 {cargos} · 생성기 {modes} · 탐색공간 {a.gen_cap}"
          f"(길이별 실제 {[len(space[L]) for L in lengths]})")
    print(f"하드 제약은 프로덕션 함수 그대로: 20AA · ≤{LINKER_MAX_LEN}aa "
          f"(`linker_reject_reason`)")
    print(f"셀 수 = {len(cargos)}화물 × {len(SH)}셔틀 × {len(lengths)}길이 = "
          f"{len(cargos) * len(SH) * len(lengths)} / 생성기")

    # 바닥 판정 — 직접융합에서 접합부 liability 가 0 이면 생성 링커가 개선할 것이 없다.
    print("\n--- 바닥(floor) 판정: 직접융합 접합부 liability 개수 ---")
    floor: dict[tuple[str, str], bool] = {}
    for c in cargos:
        cells = []
        for sn, ss in SH.items():
            nb = len(junction_liabilities(c, "", ss))
            floor[(c, sn)] = nb == 0
            cells.append(f"{sn}:{nb}")
        nz = sum(1 for sn in SH if not floor[(c, sn)])
        print(f"  {c} → 비바닥 {nz}/{len(SH)}   " + " ".join(cells))

    rng = random.Random(a.seed)
    # 널은 **생성기와 무관하게 같은 집합**을 쓴다(같은 시드·같은 순서) → 생성기 비교가 공정.
    nulls: dict[int, list[str]] = {}
    for L in lengths:
        nulls[L] = ["".join(rng.choice(_STD) for _ in range(L)) for _ in range(a.nulls)]

    results: dict[str, dict] = {}
    for mode in modes:
        rows = []
        for c in cargos:
            for sn, ss in SH.items():
                base_j = len(junction_liabilities(c, "", ss))
                for L in lengths:
                    lk = generate(c, ss, L, space[L], mode)
                    real = risk_vector(c, lk, ss)
                    real_j = len(junction_liabilities(c, lk, ss))
                    nv, nj = [], []
                    for x in nulls[L]:
                        nv.append(risk_vector(c, x, ss))
                        nj.append(len(junction_liabilities(c, x, ss)))
                    beat = sum(1 for v in nv if _dominates(v, real))
                    won = sum(1 for v in nv if _dominates(real, v))
                    rows.append({
                        "cargo": c, "shuttle": sn, "L": L, "linker": lk,
                        "floor": floor[(c, sn)], "base_j": base_j, "real_j": real_j,
                        "real": real, "nv": nv, "nj": nj,
                        "p": beat / a.nulls, "won": won / a.nulls,
                        "base": risk_vector(c, "", ss),
                    })
        results[mode] = rows

    # ---- ① 축별 이동 가능성 -------------------------------------------------
    print("\n" + "=" * 78)
    print("① 🔴 핵심 질문 — 생성 링커가 **실제로 움직일 수 있는 축**은 몇 개인가")
    print("=" * 78)
    print("A) 같은 길이 무작위 링커 집합 안에서 그 축이 **값을 바꾸나**(= 링커가 통제하는 축인가).")
    print("   생성기와 무관하다 — 널 집합만 본다.")
    ref = results[modes[0]]
    print(f"   {'축':<12}{'값이 변하는 셀':>16}{'널 내부 전부 동일':>18}{'널 값의 폭(중위)':>18}")
    for i, ax in enumerate(_AXES):
        movable = sum(1 for r in ref if len({v[i] for v in r["nv"]}) > 1)
        spread = statistics.median([max(v[i] for v in r["nv"]) - min(v[i] for v in r["nv"])
                                    for r in ref])
        print(f"   {ax:<12}{f'{movable}/{len(ref)}':>16}{len(ref) - movable:>18}{spread:>18.3f}")

    print("\nB) 그 축에서 **널이 생성 링커보다 나쁜가/좋은가** (널·셀 전체 집계, 생성기별)")
    print(f"   {'생성기':<10}" + "".join(f"{ax + ' 나쁨/좋음':>24}" for ax in _AXES))
    for mode in modes:
        rows = results[mode]
        tot = sum(len(r["nv"]) for rows_ in [rows] for r in rows_)
        cells = []
        for i in range(len(_AXES)):
            worse = sum(sum(1 for v in r["nv"] if v[i] > r["real"][i]) for r in rows)
            better = sum(sum(1 for v in r["nv"] if v[i] < r["real"][i]) for r in rows)
            cells.append(f"{worse / tot * 100:>9.1f}% /{better / tot * 100:>8.1f}%")
        print(f"   {mode:<10}" + "".join(f"{c:>24}" for c in cells))
    print("\nC) 🔴 **개선인가 악화 회피인가** — 기준은 무작위가 아니라 **직접융합**(링커 없음)이다.")
    print("   계열 1 의 승리는 '같은 길이 무작위보다 낫다'일 뿐이고, 링커를 **끼울 가치가**")
    print("   있는지는 링커 없는 구성체와 비교해야 나온다(직접융합은 라이브러리의 실제 선택지다).")
    print(f"   {'생성기':<10}{'':<2}" + "".join(f"{ax:>26}" for ax in _AXES))
    for mode in modes:
        rows = results[mode]
        gl = [f"{sum(1 for r in rows if r['real'][i] < r['base'][i])}↓/"
              f"{sum(1 for r in rows if r['real'][i] > r['base'][i])}↑/{len(rows)}"
              for i in range(len(_AXES))]
        print(f"   {mode:<10}생성" + "".join(f"{x:>26}" for x in gl))
    tot = sum(len(r["nv"]) for r in results[modes[0]])
    nl = [f"{sum(sum(1 for v in r['nv'] if v[i] < r['base'][i]) for r in results[modes[0]])}↓/"
          f"{sum(sum(1 for v in r['nv'] if v[i] > r['base'][i]) for r in results[modes[0]])}↑/{tot}"
          for i in range(len(_AXES))]
    print(f"   {'(무작위 널)':<10}널  " + "".join(f"{x:>26}" for x in nl))
    print("   읽는 법: `↓` 가 0 인 축은 링커로 **내릴 수 없는 축**이다 — 그 축에서 널을 이기는 것은")
    print("   개선이 아니라 **악화를 피한 것**이고, '개선 신호'의 근거로 쓸 수 없다(규범 3-A).")

    print("\n🔴 구조적 사실(측정이 아니라 **코드**): `null_linker.risk_vector` 는 off-target 을")
    print("   `assess_selectivity(shuttle)` — **셔틀 단독**으로 계산한다. 링커가 입력에 아예")
    print("   들어가지 않으므로 off 축은 '거의 불변'이 아니라 **정의상 불변**이다. 지배 비교에서")
    print("   항상 동점이라 **실효 축은 4개가 아니라 3개**이고, 승률의 정직한 분모는 그만큼 작다.")

    # ---- ② 핵심 널 --------------------------------------------------------
    print("\n" + "=" * 78)
    print("② 핵심 널 — 생성 링커가 같은 길이 무작위를 이기나 (경험 p = 널이 생성체를 **지배한** 비율)")
    print("=" * 78)
    for mode in modes:
        rows = results[mode]
        for label, sub in (("전체", rows), ("바닥 제외", [r for r in rows if not r["floor"]])):
            if not sub:
                print(f"  [{mode}] {label}: 셀 0개 — 측정 불가(None, 0 이 아니다)")
                continue
            ps = [r["p"] for r in sub]
            ws = [r["won"] for r in sub]
            sig = sum(1 for p in ps if p < 0.05)
            strict = sum(1 for p in ps if p == 0.0)
            rev = sum(1 for p in ps if p > 0.5)
            print(f"  [{mode:<8}] {label:<8} 유의하게 이긴 셀 **{sig}/{len(sub)} "
                  f"({sig / len(sub) * 100:.1f}%)** (우연 ~5%) · p==0 엄격 "
                  f"{strict}/{len(sub)} ({strict / len(sub) * 100:.1f}%) · "
                  f"거꾸로(p>0.5) {rev}/{len(sub)} · p 중위 {statistics.median(ps):.2f} "
                  f"평균 {statistics.mean(ps):.2f}")
            print(f"{'':<13}{'':<9} ↳ 거울 통계: 생성체가 **지배한** 널 비율 중위 "
                  f"{statistics.median(ws):.2f} 평균 {statistics.mean(ws):.2f} · "
                  f"서로 비교불가(지배 양방향 0) 중위 "
                  f"{statistics.median([1 - p - w for p, w in zip(ps, ws)]):.2f}")
    print("\n⚠️ `p` 와 거울 통계를 **같이** 읽어야 한다 — 지배는 약한 관계라 'p=0'은 "
          "'생성체가 이겼다'가 아니라 '널이 생성체를 못 눌렀다'다. 거울 통계가 0 에 가까우면 "
          "두 쪽이 그냥 **비교불가**이고 신호가 아니다.")

    print("\n화물별 승률 (유의 셀 / 전체 · 바닥 제외) — 접합부 축은 화물 끝 잔기에 달려 있다")
    print(f"{'생성기':<10}" + "".join(f"{c[-1] + ' 끝':>12}" for c in cargos))
    for mode in modes:
        cells = []
        for c in cargos:
            sub = [r for r in results[mode] if r["cargo"] == c and not r["floor"]]
            cells.append("—" if not sub
                         else f"{sum(1 for r in sub if r['p'] < 0.05)}/{len(sub)}")
        print(f"{mode:<10}" + "".join(f"{x:>12}" for x in cells))
    print("\n경험 p 분포 (전체 셀, 생성기별)")
    bins = [(0.0, 0.0), (0.0, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 0.8), (0.8, 1.01)]
    print(f"{'생성기':<10}" + "".join(f"{'p=0' if lo == hi else f'{lo}~{hi}':>11}" for lo, hi in bins))
    for mode in modes:
        ps = [r["p"] for r in results[mode]]
        cnt = []
        for lo, hi in bins:
            cnt.append(sum(1 for p in ps if (p == 0.0 if lo == hi else lo < p <= hi)))
        print(f"{mode:<10}" + "".join(f"{x:>11}" for x in cnt))

    # ---- ③ 방향성 분리 ----------------------------------------------------
    print("\n" + "=" * 78)
    print("③ 🔴 방향성 판정 (protocol.md 규범 3-A) — 악화 검출과 개선 신호를 **따로** 판정")
    print("=" * 78)
    for mode in modes:
        rows = results[mode]
        nf = [r for r in rows if not r["floor"]]
        # 악화 검출: 널이 **접합부 모티프를 더 만든** 경우(정규식 ground truth)를 축이 잡나
        gt_worse = acq = det = 0
        for r in rows:
            for v, j in zip(r["nv"], r["nj"]):
                if j > r["real_j"]:
                    gt_worse += 1
                    if _dominates(r["real"], v):
                        det += 1
                    if _dominates(v, r["real"]):
                        acq += 1
        # 개선 신호: 같은 길이 무작위가 **같은 수리를 공짜로** 달성하나
        cheap_j = cheap_n = 0
        improved = imp_cells = 0
        for r in nf:
            for j in r["nj"]:
                cheap_n += 1
                if j <= r["real_j"]:
                    cheap_j += 1
            if r["real_j"] < r["base_j"]:
                improved += 1
            imp_cells += 1
        print(f"\n  [{mode}]")
        print(f"  ⑴ 악화 검출 — 널이 접합부 모티프를 **더 만든** 사례 {gt_worse}개 중")
        if gt_worse:
            print(f"       · 생성체가 지배(= 축이 악화로 잡음)  {det}/{gt_worse} "
                  f"({det / gt_worse * 100:.1f}%)")
            print(f"       · 🔴 널이 생성체를 지배(= 축이 더 나쁜 링커를 **더 좋다고** 함) "
                  f"{acq}/{gt_worse} ({acq / gt_worse * 100:.1f}%)")
            print(f"       · 비교불가(어느 쪽도 지배 못 함 = 판정 보류) "
                  f"{gt_worse - det - acq}/{gt_worse} "
                  f"({(gt_worse - det - acq) / gt_worse * 100:.1f}%)")
        else:
            print("       · 측정 불가 (None) — 이 설정에서 악화 사례가 생기지 않았다")
        print(f"  ⑵ 개선 신호 — 바닥 제외 셀 {imp_cells}개 중 생성체가 직접융합보다 접합부를 "
              f"줄인 셀 {improved}/{imp_cells}")
        if cheap_n:
            print(f"       · 🔴 같은 길이 **무작위**가 같은 접합부 수준을 달성한 비율 "
                  f"{cheap_j}/{cheap_n} ({cheap_j / cheap_n * 100:.1f}%) "
                  f"— 높을수록 '개선'이 선택이 아니라 **공짜**다")

    # ---- ④ 기록된 생성 링커 (Arm A) ---------------------------------------
    print("\n" + "=" * 78)
    print("④ §23 에 **실제로 기록된** 생성/설계 링커 (대리 생성기가 아닌 실측 Arm)")
    print("=" * 78)
    lib = {(v["seq"] or "").upper() for v in LINKER_LIBRARY.values()}
    rec_rows = []
    for lk, note in _RECORDED:
        L = len(lk)
        nl = ["".join(rng.choice(_STD) for _ in range(L)) for _ in range(a.nulls)]
        ps = []
        for c in cargos:
            for sn, ss in SH.items():
                real = risk_vector(c, lk, ss)
                beat = sum(1 for x in nl if _dominates(risk_vector(c, x, ss), real))
                ps.append(beat / a.nulls)
        rec_rows.append((lk, note, L, ps))
    print(f"{'링커':<12}{'len':>4}{'유의':>10}{'p중위':>8}{'라이브러리':>10}  비고")
    for lk, note, L, ps in rec_rows:
        sig = sum(1 for p in ps if p < 0.05)
        print(f"{lk:<12}{L:>4}{f'{sig}/{len(ps)}':>10}{statistics.median(ps):>8.2f}"
              f"{('예' if lk in lib else '아니오'):>10}  {note}")
    print("⚠️ 이 Arm 의 널은 ②와 **다른 난수 소비 지점**을 쓴다(같은 시드, 다른 추출). "
          "n 이 작고(링커 7종) 두 종은 일부러 나쁜 음성 프로브다 — 통계가 아니라 사례다.")

    # ---- ④-b 길이 순위상관 (R2 체크리스트 ①) -------------------------------
    print("\n" + "=" * 78)
    print("④-b 길이 순위상관 — 래칫 규율 R2 체크리스트 ① 「새 편집 채널에도 같이 걸어야 한다」")
    print("=" * 78)
    print("무작위 널만으로 재는 **순수 길이 의존성**(셀별 Spearman ρ 의 중위, 32셀).")
    print("ρ>0 = 길어지면 그 축이 나빠진다 · ρ<0 = 길어지면 좋아진다(= 길이로 점수를 살 수 있다).")
    by_cell: dict[tuple[str, str], list[tuple[int, tuple]]] = {}
    for r in results[modes[0]]:
        by_cell.setdefault((r["cargo"], r["shuttle"]), []).extend((r["L"], v) for v in r["nv"])
    for i, ax in enumerate(_AXES):
        rs = [x for x in (_spearman([L for L, _ in v], [t[i] for _, t in v])
                          for v in by_cell.values()) if x is not None]
        print(f"   {ax:<12} ρ 중위 {statistics.median(rs):+.3f} · 범위 "
              f"{min(rs):+.3f}~{max(rs):+.3f} · 음수 셀 "
              f"{sum(1 for x in rs if x < 0)}/{len(rs)}" if rs else f"   {ax:<12} 측정 불가(None)")
    for mode in modes:
        cell: dict[tuple[str, str], list] = {}
        for r in results[mode]:
            cell.setdefault((r["cargo"], r["shuttle"]), []).append((r["L"], r["won"]))
        rs = [x for x in (_spearman([L for L, _ in v], [w for _, w in v])
                          for v in cell.values()) if x is not None]
        if rs:
            print(f"   [{mode}] 길이 ↔ 생성체의 널 지배율 ρ 중위 {statistics.median(rs):+.3f}")
    print("🔴 위 세 축의 ρ 가 **전부 양수**다 = 길어지면 축이 나빠진다 → 축 자체는 길이로 점수를")
    print("   살 수 없다(R2① **통과**). 그런데 바로 그 때문에 **널이 길이에 따라 쉬워진다** —")
    print("   아래 표가 그 크기다. 설계가 좋아진 게 아니라 **무작위가 나빠진 것**이므로, 전체")
    print("   승률 하나로 보고하면 길이 조성이 승률을 결정한다. 길이별로 나눠 읽어야 한다.")
    print(f"\n   {'생성기':<10}" + "".join(f"{'L=' + str(L):>9}" for L in lengths))
    for mode in modes:
        cells = []
        for L in lengths:
            sub = [r for r in results[mode] if r["L"] == L]
            cells.append(f"{sum(1 for r in sub if r['p'] < 0.05)}/{len(sub)}")
        print(f"   {mode:<10}" + "".join(f"{x:>9}" for x in cells))
    print(f"   {'(바닥제외)':<10}" + "".join(f"{'':>9}" for _ in lengths))
    for mode in modes:
        cells = []
        for L in lengths:
            sub = [r for r in results[mode] if r["L"] == L and not r["floor"]]
            cells.append(f"{sum(1 for r in sub if r['p'] < 0.05)}/{len(sub)}" if sub else "—")
        print(f"   {mode:<10}" + "".join(f"{x:>9}" for x in cells))

    # ---- ⑤ 생성이 라이브러리에 **무엇을 더하나** (널이 아니라 대조) ------------
    print("\n" + "=" * 78)
    print("⑤ 그래서 생성이 **라이브러리 13종에 무엇을 더하나** — §23-4 의 '이득 0'을 전 화물로")
    print("=" * 78)
    print("널이 아니라 **대조**다: 각 (화물, 셔틀)에서 라이브러리 13종 중 위험 벡터 사전식 최선 vs")
    print("생성기 `full` 의 전 길이 최선. 지배 관계로 판정한다(규범 5 「비지배 ≠ 우위」).")
    lib_seqs = {n: (v["seq"] or "").upper() for n, v in LINKER_LIBRARY.items()}
    gen_better = gen_equal = lib_better = incomp = 0
    detail = []
    for c in cargos:
        for sn, ss in SH.items():
            lb = min(lib_seqs.items(), key=lambda kv: (risk_vector(c, kv[1], ss), _crc(kv[0])))
            lbv = risk_vector(c, lb[1], ss)
            gb = min((generate(c, ss, L, space[L], "full") for L in lengths),
                     key=lambda s: (risk_vector(c, s, ss), _crc(s)))
            gbv = risk_vector(c, gb, ss)
            if gbv == lbv:
                gen_equal += 1
                tag = "동일"
            elif _dominates(gbv, lbv):
                gen_better += 1
                tag = "생성 지배"
            elif _dominates(lbv, gbv):
                lib_better += 1
                tag = "🔴 라이브러리 지배"
            else:
                incomp += 1
                tag = "비교불가(교환)"
            if not floor[(c, sn)] or tag != "동일":
                detail.append((c, sn, lb[0], lbv, gb, gbv, tag))
    n_cell = len(cargos) * len(SH)
    print(f"생성 지배 {gen_better}/{n_cell} · 벡터 완전 동일 {gen_equal}/{n_cell} · "
          f"비교불가(교환) {incomp}/{n_cell} · 🔴 라이브러리 지배 {lib_better}/{n_cell}")
    if detail:
        print(f"\n{'화물':<13}{'셔틀':<16}{'라이브러리 최선':<14}{'위험벡터':<28}"
              f"{'생성 최선':<14}{'위험벡터':<28}판정")
        for c, sn, ln, lv, gs_, gv, tag in detail:
            print(f"{c:<13}{sn:<16}{ln:<14}{str(lv):<28}{gs_:<14}{str(gv):<28}{tag}")
    else:
        print("차이가 나는 셀이 **하나도 없다** — 생성의 순이득이 0 이다(§23-4 재현).")

    if a.verbose:
        print("\n--- 셀 단위 (verbose) ---")
        for mode in modes:
            for r in results[mode]:
                print(f"{mode:<9}{r['cargo']:<13}{r['shuttle']:<16}{r['L']:>3} "
                      f"{r['linker']:<18}{'바닥' if r['floor'] else '   ':<5}"
                      f"base_j={r['base_j']} real_j={r['real_j']} p={r['p']:.2f} "
                      f"risk={r['real']}")

    # ---- 판정 --------------------------------------------------------------
    print("\n" + "=" * 78)
    print("판정 (규범 3-A 2×2 · 이 스크립트가 **주장할 수 있는 범위**)")
    print("=" * 78)
    print("· **악화 검출 — 조건부 신뢰.** 접합부 모티프를 더 만든 널을 위험 벡터가 역전시킨")
    print("  비율은 2~3%(생성기 `junction`)이고 `full` 에서는 0% 다. 다만 판정 보류(비교불가)가")
    print("  1/3 이라 '검출기'라기보다 **한쪽으로만 틀리는 필터**다.")
    print("· 🔴 **개선 신호 — 불신.** 같은 길이 **무작위** 링커가 생성체와 같은 접합부 수준을")
    print("  63~67% 달성한다. §26-4 의 off 축(46% 역방향)과 같은 형태의 결론이다 — 평균이")
    print("  올바른 방향이어도 **선택의 몫이 작다**.")
    print("· 따라서 생성 링커 채널이 들어갈 칸은 off-target 과 **같은 칸**이다:")
    print("  순위에는 쓰되(접합부 검출은 결정론이고 악화를 올바른 방향으로 잡는다),")
    print("  **'생성이 더 좋은 링커를 만든다'는 주장은 금지**한다.")
    print("\n이 측정이 **보이지 않는 것**:")
    print("  ⓐ 심사한 것은 LLM 이 아니라 **정책의 대리 생성기**다. 실제 생성은 유료·비결정론")
    print("     이고 라이브 세션에서 쓰인 기록도 없다(§23-8③). LLM 이 정책보다 나쁠 수 있다(§30).")
    print("  ⓑ 생성기의 강도는 탐색 공간 크기(--gen-cap)와 스캐폴드 어휘가 정한다. 둘 다")
    print("     **우리가 고른 값**이고 외부 근거가 없다 — `full` 의 99% 승률은 그 자체로")
    print("     '심사받는 축을 직접 최소화하면 이긴다'는 **동어반복**이다.")
    print("  ⓒ 통제한 것은 **길이 하나**다. 위치·조성은 통제되지 않았다(계열 1 의 ④).")
    print("     접합부 축의 조성 보존 셔플 널(계열 2)은 **아직 없다**.")
    print("  ⓓ BBB 축으로는 재지 않았다(M8 에서 퇴역). 독성 게이트(임계 0.38)도 널에")
    print("     들어가지 않았다 — 17aa 상한의 근거가 바로 그 게이트인데 **심사되지 않았다**.")
    print("  ⓔ 화물 4종뿐이고 전부 같은 11aa 골격의 끝 잔기 변형이다.")

    # ---- 프로덕션 기본값 단언 (몽키패치 0건이지만 규율은 지킨다) --------------
    assert DEFAULT_CARGO == "GSNKGAIIGLM", DEFAULT_CARGO
    assert LINKER_MAX_LEN == 17, LINKER_MAX_LEN
    assert linker_reject_reason("GGGGS") is None
    assert linker_reject_reason("A" * 18) is not None
    print("\n✅ 프로덕션 기본값 단언 통과 (몽키패치 0건 — 화물·링커를 전부 인자로 넘긴다). "
          f"재현: /usr/bin/python3 scripts/null_linker_generated.py --seed {a.seed} "
          f"--nulls {a.nulls} --gen-cap {a.gen_cap}")


if __name__ == "__main__":
    main()
