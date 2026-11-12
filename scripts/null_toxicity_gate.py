#!/usr/bin/env python3
"""§34-10ⓓ — 독성 게이트(임계 0.38)의 대응 널. **17aa 상한의 근거를 피고석에**.

왜 돌렸나
---------
`core/optimizer_agent.py:LINKER_MAX_LEN = 17` 의 **유일한** 근거는 독성 게이트다.
주석이 적는 기전 주장은 이것이다 — *"링커를 늘리면 구성체 독성이 단조 감소한다
(GGS 0.27 → (GGGGS)3 0.20). 독성은 순위 축이 아니라 하드 게이트(0.38)라서 길이로
희석하면 독성 셔틀이 게이트를 **우회**할 수 있다."*

즉 상한은 "독성을 낮추려고" 둔 게 아니라 **희석을 구멍으로 보고 막으려고** 둔 것이다.
그런데 §34-8 에서 `full` 생성기의 전 길이 최선이 32셀 중 23셀에서 15~17aa, **14셀이
정확히 17aa 상한**임이 확인됐다 — 상한이 실제로 조여 있다(binding). **지금 가장 조여
있는 제약이 널에 한 번도 들어가지 않았다.**

이 스크립트가 묻는 네 가지
--------------------------
Q1 **기전이 존재하나.** 구성체 = `cargo + linker + shuttle` 에서 ToxinPred3 ML Score 가
   링커 길이 L 의 **단조 감소 함수**인가. L = 0..25 로 **상한을 넘겨서** 쓸어 본다.
   링커 어휘 4계열(GS · EAAAK · poly-G · poly-K)을 같은 길이 격자에 올려 **길이와 조성을
   분리**한다. 어디서 0.38 을 가로지르나.
Q2 **상한이 자기 근거와 싸우고 있나.** 독성이 길이에 단조 감소라면 17 상한은 **더 이상의
   희석을 막는다**. (화물, 셔틀) 쌍 중 L≤17 에서 전멸이지만 L≤25 에서 통과하는 쌍이
   몇 개인가. 그리고 그 반대 — **상한 안에서 이미 우회가 되는** 쌍은 몇 개인가.
Q3 **게이트의 대응 널**(계열 1). 같은 길이 **무작위** 링커 vs 라이브러리 어휘 링커.
   게이트의 통과/탈락 판정이 **길이를 넘는 정보**를 담는가. 같은 길이에서 통과율이 같다면
   게이트는 서열이 아니라 **길이를 재고 있다**.
   보조로 계열 2(조성 보존 셔플)를 라이브러리 최장 링커 `A(EAAAK)3A`(= 상한을 정의하는
   그 서열)에 돌린다 — ToxinPred3 Model 1 은 AAC(순서 불변) + DPC(순서 민감)이므로
   셔플 폭이 **DPC 성분의 크기**를 준다.
Q4 **17 이 특별한가.** 15·17·20 에서 데이터가 달라지는 게 있는지 본다.

표본 수 (명시)
--------------
화물 2종(`GSNKGAIIGLM` Aβ(25-35) 기본 · `WAGGDASGE` DSIP — §35 와 같은 cargo B)
셔틀 7종(고전하 TAT·Penetratin·B6·RVG29 / 저전하 ApoE(159-167)2·T7·Angiopep-2)

| 블록 | 구성 | 서열 수 |
|---|---|---|
| Q1 길이 격자 | 2 화물 × 7 셔틀 × L 0..25(26) × 어휘 4계열 | 1,456 |
| Q3 계열 1 널 | 2 화물 × 7 셔틀 × L 1..25(25) × 무작위 20개 | 7,000 |
| Q3 계열 2 널 | 2 화물 × 7 셔틀 × `A(EAAAK)3A` 순열 20개 | 280 |
| 참조: 실제 라이브러리 | 2 화물 × 7 셔틀 × 링커 13종 | 182 |
| 결정론 재측정 | 위 중 48개를 두 번째 배치로 재호출 | 48 |

총 **약 8,966 서열**. 실측 처리량 ≈ 41 ms/서열(배치당 고정비 ≈ 2s) → **약 6~7분**.
🔴 ToxinPred3 는 **직렬화**되어 있다(`core/toxicity.py` 의 flock). 병렬화 금지 —
공유 cwd 에 고정 이름 중간파일을 쓰기 때문에 동시 호출이 **점수를 조용히 오염**시킨다.

이 측정이 **볼 수 없는 것** (결론과 함께 반드시 읽을 것)
-------------------------------------------------------
① **ToxinPred3 ML Score 가 독성인지는 심사하지 않는다.** 여기서 재는 양은 오직
   *"프로덕션 게이트가 쓰는 그 숫자"* 다. 외부 라벨 검증(계열 3)은 **없다** — 즉
   "길이가 늘면 독성이 준다"가 아니라 **"길이가 늘면 ToxinPred3 점수가 준다"** 가 결론이다.
   조성 기반 분류기에서 그건 **기전이 아니라 산술일 수 있다**(AAC 분모가 커진다).
② **화물 2종·셔틀 7종**이다. 라이브러리 셔틀 12종 중 5종(Leptin30·CDX·THR·Tet1·SynB1)은
   들어 있지 않다.
③ 무작위 널은 **균일 20AA** 다. 조성이 라이브러리 어휘(G/S/E/A/K)와 전혀 다르므로
   "길이 통제"만 되고 조성은 통제되지 않는다(계열 1 의 ④ 참조).
④ 임계 0.38 **자체의 근거**는 심사 범위 밖이다(ToxinPred3 기본값을 그대로 받았다).
⑤ 구성체 **접힘·제조성·실제 세포 독성**은 전혀 재지 않는다.

실행: `/usr/bin/python3 scripts/null_toxicity_gate.py`   (⚠️ `python3 -I` 금지)
"""
from __future__ import annotations

import random
import statistics as st
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import DEFAULT_CARGO, LINKER_LIBRARY, SHUTTLES, get_settings  # noqa: E402
from core.optimizer_agent import LINKER_MAX_LEN  # noqa: E402
from core.toxicity import get_toxicity_predictor  # noqa: E402

# ---------------------------------------------------------------------------
# 파라미터 — 이 스크립트의 **모든** 난수·표본 출처를 여기 적는다
# ---------------------------------------------------------------------------
CARGO_A = DEFAULT_CARGO              # Aβ(25-35) · 11aa (프로덕션 기본)
CARGO_B = "WAGGDASGE"                # DSIP · 9aa  (§35 의 cargo B 와 동일)
CARGOS = {"Aβ25-35": CARGO_A, "DSIP": CARGO_B}

#: 고전하 4종 + 저전하 3종. 저전하 쪽에 등급 A(Angiopep-2)를 넣어 "근거 있는 셔틀이
#: 게이트에서 어떻게 되는가"가 같이 보이게 한다(폐지된 MIN_BBB_SIGNAL 의 교훈).
TARGETS = ("TAT", "Penetratin", "B6", "RVG29", "ApoE(159-167)2", "T7", "Angiopep-2")

L_MAX = 25                           # 🔴 상한(17)을 **넘겨서** 쓸어본다
LENGTHS = tuple(range(0, L_MAX + 1))
N_RAND = 20                          # 계열 1: 셀당 같은 길이 무작위 링커 수
N_SHUF = 20                          # 계열 2: 셀당 조성 보존 순열 수
RNG_SEED = 20261112                  # 전역 시드(셀별로 결정론 파생)
AA20 = "ACDEFGHIKLMNPQRSTVWY"
POS, NEG = "KR", "DE"

#: 길이 격자에 올리는 링커 어휘 — 같은 L 에서 **조성만** 다르다.
MOTIFS = {"GS": "GGGGS", "EAAAK": "EAAAK", "polyG": "G", "polyK": "K"}
CAP_LINKER = LINKER_LIBRARY["A(EAAAK)3A"]["seq"]   # 17aa — 상한을 **정의하는** 서열


def motif_linker(motif: str, length: int) -> str:
    """모티프를 순환 반복해 정확히 `length` aa 로 자른다. length=0 → 직접융합."""
    if length == 0:
        return ""
    m = MOTIFS[motif]
    return (m * (length // len(m) + 1))[:length]


def _rng(*parts) -> random.Random:
    """셀 좌표 → 결정론 RNG. 세션·실행 순서에 불변."""
    return random.Random(f"{RNG_SEED}:" + ":".join(str(p) for p in parts))


def net_charge(s: str) -> int:
    return sum(c in POS for c in s) - sum(c in NEG for c in s)


# ---------------------------------------------------------------------------
# 채점 — 프로덕션 경로 그대로. 서브프로세스 호출 수를 **센다**(계측의 계측).
# ---------------------------------------------------------------------------
class Scorer:
    """`ToxinPred3LocalPredictor.predict_many` 를 배치로 호출하고 캐시한다.

    🔴 `predict_many` 를 **직접** 쓴다(재구현하지 않는다). `is_toxic` 판정도
    예측기가 내린 것을 그대로 읽는다 — 임계값을 여기서 다시 쓰면 드리프트한다.
    """

    def __init__(self):
        self.pred = get_toxicity_predictor(get_settings())
        if self.pred.source != "toxinpred3-local":
            raise SystemExit(
                f"🔴 실측 예측기가 아니다(source={self.pred.source!r}). "
                "폴백 상수로는 이 측정이 무의미하다 — vendor/toxinpred3 연결을 확인할 것.")
        self.cache: dict[str, tuple[float, bool]] = {}
        self.n_batches = 0           # 서브프로세스 호출 수
        self.n_scored = 0            # 실제로 채점된 서열 수(캐시 미스)
        self.seconds = 0.0

    def score(self, seqs: list[str], batch: int = 1200) -> None:
        todo = sorted({s for s in seqs if s not in self.cache})
        for i in range(0, len(todo), batch):
            chunk = todo[i:i + batch]
            t0 = time.time()
            res = self.pred.predict_many(chunk)      # ← 프로덕션과 동일 호출
            self.seconds += time.time() - t0
            self.n_batches += 1
            self.n_scored += len(chunk)
            for s, r in zip(chunk, res):
                self.cache[s] = (r.risk, r.is_toxic)
            print(f"   ... 채점 {min(i + batch, len(todo)):>5}/{len(todo)} "
                  f"(누적 {self.seconds:.0f}s)", flush=True)

    def risk(self, s: str) -> float:
        return self.cache[s][0]

    def toxic(self, s: str) -> bool:
        return self.cache[s][1]


def construct(cargo: str, linker: str, shuttle: str) -> str:
    """프로덕션 조립식(`optimizer_agent._enumerate`: cargo + lk + sh)."""
    return cargo + linker + shuttle


def _w(s, n):   # 한글 폭 보정 없이 단순 좌정렬(표 비교용)
    return str(s).ljust(n)


# ---------------------------------------------------------------------------
def main() -> None:
    sc = Scorer()
    thr = sc.pred.threshold
    print(f"예측기 {sc.pred.source} · 임계 {thr} · 상한 LINKER_MAX_LEN = {LINKER_MAX_LEN}")
    print(f"화물 {len(CARGOS)}종 · 셔틀 {len(TARGETS)}종 · L = 0..{L_MAX} · "
          f"어휘 {len(MOTIFS)}계열 · 무작위 {N_RAND}/셀 · 셔플 {N_SHUF}/셀\n")

    # ── 서열 조립 (전부 한 번에 모아 배치 채점) ────────────────────────────
    grid: dict = {}      # (cargo, shuttle, motif, L) -> seq
    nulls: dict = {}     # (cargo, shuttle, L) -> [seq] (계열 1)
    shufs: dict = {}     # (cargo, shuttle) -> [seq]    (계열 2)
    libs: dict = {}      # (cargo, shuttle, 링커명) -> seq
    for cn, cargo in CARGOS.items():
        for sn in TARGETS:
            sh = SHUTTLES[sn]["seq"]
            for mo in MOTIFS:
                for L in LENGTHS:
                    grid[(cn, sn, mo, L)] = construct(cargo, motif_linker(mo, L), sh)
            for L in LENGTHS[1:]:
                r = _rng(cn, sn, L)
                nulls[(cn, sn, L)] = [
                    construct(cargo, "".join(r.choice(AA20) for _ in range(L)), sh)
                    for _ in range(N_RAND)]
            r = _rng("shuf", cn, sn)
            shufs[(cn, sn)] = [
                construct(cargo, "".join(r.sample(CAP_LINKER, len(CAP_LINKER))), sh)
                for _ in range(N_SHUF)]
            for ln, v in LINKER_LIBRARY.items():
                libs[(cn, sn, ln)] = construct(cargo, v["seq"] or "", sh)

    allseq = (list(grid.values()) + [s for v in nulls.values() for s in v]
              + [s for v in shufs.values() for s in v] + list(libs.values()))
    print(f"조립 완료: 고유 서열 {len(set(allseq))}개 (중복 포함 {len(allseq)})")
    sc.score(allseq)
    print(f"✅ 채점 {sc.n_scored}서열 / 서브프로세스 {sc.n_batches}회 / {sc.seconds:.0f}s\n")

    # ── 하네스 검증 ①: 주석이 인용한 실측치 재현 (GGS 0.27 → (GGGGS)3 0.20) ──
    c = construct(CARGO_A, "GGS", SHUTTLES["T7"]["seq"])
    d = construct(CARGO_A, "GGGGSGGGGSGGGGS", SHUTTLES["T7"]["seq"])
    print("하네스 검증 ① — optimizer_agent.py 주석의 인용 실측치")
    print(f"   화물 Aβ25-35 + T7:  GGS(3aa) {sc.risk(c):.3f}  →  "
          f"(GGGGS)3(15aa) {sc.risk(d):.3f}   [주석: 0.27 → 0.20]")

    # ── 하네스 검증 ②: 결정론 (같은 서열 두 번째 배치 = 같은 점수) ─────────
    probe = sorted(set(allseq))[::187][:48]
    before = {s: sc.risk(s) for s in probe}
    sc.cache = {k: v for k, v in sc.cache.items() if k not in before}
    sc.score(probe)
    bad = [s for s in probe if abs(sc.risk(s) - before[s]) > 1e-12]
    print(f"하네스 검증 ② — 재측정 결정론: {len(probe)}서열 중 불일치 {len(bad)}개 "
          f"{'✅' if not bad else '🔴 ' + str(bad[:2])}\n")

    # ===================================================================
    print("=" * 108)
    print("Q1 기전 — ToxinPred3 ML Score vs 링커 길이 L. 🔴 재는 양은 **점수**이고 독성이 아니다")
    print("   (어휘 GS = 라이브러리 주력. ⚠️ 는 점수 > 임계 0.38 = 게이트 탈락)")
    for cn in CARGOS:
        print(f"\n── 화물 {cn} ──")
        print(_w("셔틀", 16) + _w("순전하", 7) + _w("어휘", 7)
              + "".join(f"L{L:<5}" for L in (0, 1, 3, 5, 10, 15, 17, 20, 25))
              + "  단조?  0.38교차")
        for sn in TARGETS:
            q = net_charge(SHUTTLES[sn]["seq"])
            for mo in MOTIFS:
                ser = [sc.risk(grid[(cn, sn, mo, L)]) for L in LENGTHS]
                mono = all(b <= a + 1e-9 for a, b in zip(ser, ser[1:]))
                # 0.38 을 **아래로** 가로지르는 첫 L
                cross = next((L for L in LENGTHS
                              if not sc.toxic(grid[(cn, sn, mo, L)])
                              and sc.toxic(grid[(cn, sn, mo, 0)])), None)
                cells = "".join(
                    f"{sc.risk(grid[(cn, sn, mo, L)]):.3f}{'⚠' if sc.toxic(grid[(cn, sn, mo, L)]) else ' '}"
                    for L in (0, 1, 3, 5, 10, 15, 17, 20, 25))
                print(_w(sn if mo == "GS" else "", 16) + _w(f"{q:+d}" if mo == "GS" else "", 7)
                      + _w(mo, 7) + cells
                      + f"  {'단조↓' if mono else '  -  '}  "
                      + (f"L={cross}" if cross is not None else "—"))

    print("\n   길이 순위상관(Spearman 대용: 인접 감소 비율) · 어휘별 집계")
    for mo in MOTIFS:
        drops, rises = 0, 0
        for cn in CARGOS:
            for sn in TARGETS:
                ser = [sc.risk(grid[(cn, sn, mo, L)]) for L in LENGTHS]
                for a, b in zip(ser, ser[1:]):
                    if b < a - 1e-9:
                        drops += 1
                    elif b > a + 1e-9:
                        rises += 1
        tot = drops + rises
        print(f"   {_w(mo, 7)} 인접 쌍 {tot:>4}개 중 감소 {drops:>4} ({100*drops/max(tot,1):5.1f}%) · "
              f"증가 {rises:>4} ({100*rises/max(tot,1):5.1f}%)")

    # ── Q1b 기전 주장을 **문장 그대로** 심사한다 ────────────────────────────
    # 인용된 주장(README:226 · 기획서:405 · §23-3 · optimizer_agent.py 주석):
    #   "링커를 늘리면 구성체 독성이 **단조 감소**한다 (0.27 → 0.20)"
    # 근거는 (화물 1, 셔틀 1, 링커 2점)이다. 아래는 같은 문장을 56계열로 심사한다.
    def spearman(xs, ys) -> float:
        def rank(v):
            s = sorted(range(len(v)), key=lambda i: v[i])
            r = [0.0] * len(v)
            i = 0
            while i < len(s):
                j = i
                while j + 1 < len(s) and abs(v[s[j + 1]] - v[s[i]]) < 1e-12:
                    j += 1
                avg = (i + j) / 2 + 1
                for k in range(i, j + 1):
                    r[s[k]] = avg
                i = j + 1
            return r
        rx, ry = rank(xs), rank(ys)
        mx, my = st.mean(rx), st.mean(ry)
        num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
        den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
        return num / den if den else 0.0

    print("\n   Q1b 「단조 감소」를 문장 그대로 심사 — 56계열(화물2 × 셔틀7 × 어휘4), L = 0..25")
    print(_w("어휘", 8) + _w("단조↓ 계열", 12) + _w("ρ(L,점수) 중위", 15) + _w("ρ<0 계열", 10)
          + _w("Δ(L17−L0) 중위", 16) + _w("Δ 범위", 20) + "L0통과→L≤17탈락")
    rev_total, rev_rows = 0, []
    for mo in MOTIFS:
        mono, rhos, deltas, rev = 0, [], [], 0
        for cn in CARGOS:
            for sn in TARGETS:
                ser = [sc.risk(grid[(cn, sn, mo, L)]) for L in LENGTHS]
                if all(b <= a + 1e-9 for a, b in zip(ser, ser[1:])):
                    mono += 1
                rhos.append(spearman(list(LENGTHS), ser))
                deltas.append(ser[LENGTHS.index(17)] - ser[0])
                # 🔴 **역방향 사례** — 희석 주장의 반례: 짧으면 통과인데 길면 탈락
                if not sc.toxic(grid[(cn, sn, mo, 0)]):
                    hits = [L for L in LENGTHS if L <= LINKER_MAX_LEN
                            and sc.toxic(grid[(cn, sn, mo, L)])]
                    if hits:
                        rev += 1
                        rev_rows.append((cn, sn, mo, hits[0], ser[0],
                                         max(ser[:LINKER_MAX_LEN + 1])))
        rev_total += rev
        print(_w(mo, 8) + _w(f"{mono}/14", 12) + _w(f"{st.median(rhos):+.3f}", 15)
              + _w(f"{sum(1 for r in rhos if r < 0)}/14", 10)
              + _w(f"{st.median(deltas):+.3f}", 16)
              + _w(f"{min(deltas):+.3f} … {max(deltas):+.3f}", 20) + f"{rev}/14")
    print(f"\n   🔴 역방향(짧으면 게이트 통과 · 길면 탈락) 계열 **{rev_total}/56** — "
          f"상한 17 은 이 방향을 전혀 막지 않는다(전부 L≤17 안에서 일어난다):")
    for cn, sn, mo, L1, r0, rmax in rev_rows:
        print(f"      {_w(cn, 10)}{_w(sn, 16)}{_w(mo, 7)} L=0 {r0:.3f} 통과 → "
              f"L={L1} 에서 첫 탈락 · L≤17 최대 {rmax:.3f}")

    # ===================================================================
    print("\n" + "=" * 108)
    print("Q2 상한이 자기 근거와 싸우고 있나 — (화물, 셔틀) 쌍별 **달성 가능한 최소 점수**")
    print("   후보 집합 = 어휘 4계열 + 같은 길이 무작위 20개 (그 길이에서 우리가 찾을 수 있는 최선)")
    print()
    print(_w("화물", 10) + _w("셔틀", 16) + _w("L=0", 9) + _w("min L≤17", 11)
          + _w("min L≤25", 11) + _w("최선 L(≤17)", 12) + _w("최선 L(≤25)", 12) + "판정")
    n_cap_blocks, n_already, n_hopeless, n_fine = 0, 0, 0, 0
    for cn in CARGOS:
        for sn in TARGETS:
            def pool(L):
                xs = [grid[(cn, sn, mo, L)] for mo in MOTIFS]
                xs += nulls.get((cn, sn, L), [])
                return xs
            best = {L: min(pool(L), key=sc.risk) for L in LENGTHS}
            b17 = min((L for L in LENGTHS if L <= LINKER_MAX_LEN), key=lambda L: sc.risk(best[L]))
            b25 = min(LENGTHS, key=lambda L: sc.risk(best[L]))
            r0, r17, r25 = sc.risk(best[0]), sc.risk(best[b17]), sc.risk(best[b25])
            t0, t17, t25 = sc.toxic(best[0]), sc.toxic(best[b17]), sc.toxic(best[b25])
            if t17 and not t25:
                verdict, n_cap_blocks = "🔴 상한이 통과를 막는다", n_cap_blocks + 1
            elif t0 and not t17:
                verdict, n_already = "⚠️ 상한 **안에서** 이미 우회", n_already + 1
            elif t17 and t25:
                verdict, n_hopeless = "   길이로 안 됨", n_hopeless + 1
            else:
                verdict, n_fine = "   L=0 부터 통과", n_fine + 1
            print(_w(cn, 10) + _w(sn, 16)
                  + _w(f"{r0:.3f}{'⚠' if t0 else ''}", 9)
                  + _w(f"{r17:.3f}{'⚠' if t17 else ''}", 11)
                  + _w(f"{r25:.3f}{'⚠' if t25 else ''}", 11)
                  + _w(b17, 12) + _w(b25, 12) + verdict)
    print(f"\n   요약 {len(CARGOS)*len(TARGETS)}쌍: 상한이 통과를 막는 쌍 **{n_cap_blocks}** · "
          f"상한 안에서 이미 우회되는 쌍 **{n_already}** · "
          f"길이로 해결 안 되는 쌍 {n_hopeless} · L=0 부터 통과 {n_fine}")

    # ===================================================================
    print("\n" + "=" * 108)
    print("Q3 계열 1 대응 널 — 같은 길이 무작위 링커 vs 라이브러리 어휘(GS). 게이트가 서열을 보나 길이를 보나")
    print(f"   셀 = (화물 2 × 셔틀 7) = 14. 널 {N_RAND}/셀 → L 당 {14*N_RAND}관측")
    print()
    print(_w("L", 4) + _w("GS 통과셀", 11) + _w("널 통과율", 11) + _w("널 점수 평균", 14)
          + _w("GS 점수 평균", 14) + _w("경험 p 중위", 13) + "GS가 널 하위5%인 셀")
    for L in LENGTHS[1:]:
        gs_pass = sum(not sc.toxic(grid[(cn, sn, "GS", L)]) for cn in CARGOS for sn in TARGETS)
        npass, ntot, nr, gr, ps, winlo = 0, 0, [], [], [], 0
        for cn in CARGOS:
            for sn in TARGETS:
                ns = nulls[(cn, sn, L)]
                npass += sum(not sc.toxic(x) for x in ns)
                ntot += len(ns)
                nr += [sc.risk(x) for x in ns]
                g = sc.risk(grid[(cn, sn, "GS", L)])
                gr.append(g)
                # 경험 p = (GS 이하인 널 수 + 1)/(n+1) — 낮을수록 GS 가 널보다 저독성
                p = (sum(1 for x in ns if sc.risk(x) <= g) + 1) / (len(ns) + 1)
                ps.append(p)
                if p <= 0.05:
                    winlo += 1
        print(_w(L, 4) + _w(f"{gs_pass}/14", 11) + _w(f"{100*npass/ntot:5.1f}%", 11)
              + _w(f"{st.mean(nr):.3f}", 14) + _w(f"{st.mean(gr):.3f}", 14)
              + _w(f"{st.median(ps):.3f}", 13) + f"{winlo}/14")

    print("\n   계열 2 조성 보존 셔플 — 상한을 정의하는 서열 A(EAAAK)3A(17aa)의 순열")
    print("   (AAC 는 순서 불변이므로 폭 = DPC 성분. 폭이 0 이면 게이트는 조성만 본다)")
    print(_w("화물", 10) + _w("셔틀", 16) + _w("원본", 9) + _w("셔플 평균", 11)
          + _w("셔플 최소", 11) + _w("셔플 최대", 11) + "판정 뒤집힘")
    flips = 0
    for cn in CARGOS:
        for sn in TARGETS:
            orig = construct(CARGOS[cn], CAP_LINKER, SHUTTLES[sn]["seq"])
            xs = [sc.risk(x) for x in shufs[(cn, sn)]]
            f = sum(sc.toxic(x) != sc.toxic(orig) for x in shufs[(cn, sn)])
            flips += f
            print(_w(cn, 10) + _w(sn, 16) + _w(f"{sc.risk(orig):.3f}", 9)
                  + _w(f"{st.mean(xs):.3f}", 11) + _w(f"{min(xs):.3f}", 11)
                  + _w(f"{max(xs):.3f}", 11) + f"{f}/{N_SHUF}")
    print(f"   총 판정 뒤집힘 {flips}/{len(CARGOS)*len(TARGETS)*N_SHUF}")

    # ===================================================================
    print("\n" + "=" * 108)
    print("Q4 17 이 특별한가 — 15 / 17 / 20 을 나란히")
    print(_w("양", 34) + _w("L=15", 12) + _w("L=17", 12) + _w("L=20", 12) + "17이 특별?")
    rows = []
    for L in (15, 17, 20):
        gs = [sc.risk(grid[(cn, sn, "GS", L)]) for cn in CARGOS for sn in TARGETS]
        gp = sum(not sc.toxic(grid[(cn, sn, "GS", L)]) for cn in CARGOS for sn in TARGETS)
        ns = [sc.risk(x) for cn in CARGOS for sn in TARGETS for x in nulls[(cn, sn, L)]]
        nq = sum(not sc.toxic(x) for cn in CARGOS for sn in TARGETS for x in nulls[(cn, sn, L)])
        rows.append((st.mean(gs), gp, st.mean(ns), 100 * nq / len(ns)))
    for i, lab in enumerate(["GS 점수 평균", "GS 통과셀 /14", "무작위 점수 평균", "무작위 통과율 %"]):
        v = [r[i] for r in rows]
        fmt = (lambda x: f"{x:.3f}") if i in (0, 2) else (lambda x: f"{x:.1f}")
        print(_w(lab, 34) + "".join(_w(fmt(x), 12) for x in v)
              + ("단조, 17 에 불연속 없음" if (v[0] - v[1]) * (v[1] - v[2]) > 0 else "비단조"))
    print("\n   ⚠️ 17 은 데이터에서 온 값이 아니다 — "
          "`LINKER_MAX_LEN = max(len(라이브러리 링커))` 로 **A(EAAAK)3A 의 길이**다.")

    # ── 종결: 프로덕션 상수가 변하지 않았음을 단언 ─────────────────────────
    assert LINKER_MAX_LEN == 17, LINKER_MAX_LEN
    assert get_settings().toxicity_threshold == 0.38, get_settings().toxicity_threshold
    assert DEFAULT_CARGO == "GSNKGAIIGLM", DEFAULT_CARGO
    print(f"\n단언 통과 (상한 17 · 임계 0.38 · 기본 화물 불변). "
          f"ToxinPred3 서브프로세스 {sc.n_batches}회 · {sc.seconds:.0f}s · "
          f"채점 {sc.n_scored}서열.")


if __name__ == "__main__":
    main()
