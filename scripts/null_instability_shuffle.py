#!/usr/bin/env python3
"""§42 공백 ⓐ — 불안정성 축의 **계열 2**(조성 보존 셔플) 대응 널.

왜 돌리나
--------
§42 는 불안정성 축(①셔틀 단독 연속 II)에 계열 1·3·4 를 적용하고 **계열 2 를 비워 뒀다**
(`docs/protocol.md` §6 의 새 행이 "계열 2 미적용"으로 적는다). 이 절이 그 칸을 채운다.

계열 2 가 이 축에 **정확히 맞는** 이유가 있다. Guruprasad 불안정성 지수(II)는
**이펩타이드 가중치의 합**이다. 조성 보존 셔플은

  · 아미노산 조성(AAC)을 **소수점까지** 고정하고 → **순전하·GRAVY·분자량 전부 불변**
  · **이펩타이드 조성만** 파괴한다

따라서 셔플은 "II 가 조성에서 오나 **순서**에서 오나"를 가르는 유일한 깨끗한 칼이다.
그리고 이 축에서 그 질문은 특별히 비싸다 — §42-1 이 측정한 전하 삭제 유인이
**순전하를 전혀 건드리지 않고** 해소될 수 있는지를 직접 묻기 때문이다.

묻는 네 가지
-----------
Q1 **방향.** 문헌 셔틀의 II 는 같은 조성의 무작위 순열보다 높은가 낮은가.
   경험 p = (셔플 II ≤ 원본 II 인 개수 + 1)/(n+1). 🔴 **낮은 p = 셔플이 더 "안정"** 이고,
   그것은 축이 **순서 파괴를 보상한다**는 뜻이다 — 셔플은 수용체 에피토프가 죽은 서열이다.

Q2 🔴 **공짜 길이 또 있나.** §42-5 는 "잔기 편집 0회로 임계를 넘기는 링커가 10/12" 를
   찾았다. 여기서 묻는 것은 그 자매 질문이다 — **조성 편집 0회**(= 순전하 삭제 0회)로
   ③(구성체 II 의 40 이진)을 넘기는 **순열**이 있나. 있으면 §42-1 의 전하 삭제 유인은
   두 번째로 불필요해진다.

Q3 **악화 검출.** 셔플은 기능이 죽은 음성 대조군이다(순서 파괴 = 에피토프 파괴).
   off 축은 이 검사를 **200/200** 으로 통과했다(§26-3: 셔플하면 off 가 반드시 오른다).
   불안정성 축은 같은 검사에서 어떻게 되나. "셔플이 더 안정" 비율이 높으면
   **악화 검출 실패**다 — 단 II 는 off 와 달리 애초에 기능 축이 아니므로,
   그 실패는 "축이 고장"이 아니라 **"축이 혼자서는 기능을 지키지 못한다"** 로만 읽는다.

Q4 **제약이 잡나.** II 를 개선한 셔플 중 `rmt_sim ≥ 0.85`(기능 보존 프록시)를
   통과하는 것이 몇 개인가. 0 이면 제약이 이 경로를 **구조적으로** 막는다는 뜻이고,
   그것이 이 축을 쓰는 현재 설계의 안전성 근거가 된다.

표본·결정론
----------
셔틀 **12종**(라이브러리 전수) × 순열 **200개/셔틀** = 2,400 서열.
시드 **0** 고정(`random.Random(0)`, 셔틀별로 같은 시드에서 파생) · 외부 호출 없음
(Biopython ProtParam + 로컬 정렬기) → **결정론**. 두 번 돌리면 비트 단위로 같다.
구성체 블록: 화물 `GSNKGAIIGLM` × 링커 3종(직접융합 · `GGGGS` · `A(EAAAK)3A`) × 12종 ×
200 순열 = 7,200 평가. ProtParam 은 서열당 ~0.2ms 라 전체 **1분 미만**.

🔴 이 측정이 **볼 수 없는 것**
-----------------------------
① **실제 안정성이 아니다.** §42-8ⓐ 와 같다 — Tm·ΔG·반감기 실측이 **0건**이므로
   "II 가 낮다 = 안정하다"는 이 프로젝트 어디서도 검증되지 않았다. 여기서 재는 양은
   오직 **"순위·탐색 축이 쓰는 그 숫자"** 다.
② **셔플이 기능을 잃는다는 것도 측정이 아니라 전제다.** 근거는 메커니즘(RMT 는 서열
   순서 의존 결합)과 §26-3 의 off 축 거동이고, 결합 실측은 없다. `rmt_sim` 은 프록시다.
③ **순열은 동일잔기 반복 서열에서 유효 표본이 줄어든다**(프로토콜 §2.2②). 12종에
   homopolymer 는 없지만 중복 순열은 생길 수 있어 **고유 순열 수를 함께 보고**한다.
④ 구성체 블록의 링커는 **3종**이고 화물은 **1종**이다. §42-5 의 13종 전수가 아니다.
⑤ ③(40 이진)의 교차 여부는 화물에 의존한다(§35) — 다른 화물에서 숫자가 달라진다.
"""

from __future__ import annotations

import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.binding import shuttle_similarity  # noqa: E402
from core.config import LINKER_LIBRARY, SHUTTLES  # noqa: E402
from core.developability_opt import _instability  # noqa: E402  (생산 함수 직접 호출)

N_SHUFFLE = 200
SEED = 0
RMT_FLOOR = 0.85          # pareto_scrub 의 기능 보존 제약과 동일
II_STABLE = 40.0          # ③ 이진 임계 (core/developability_opt)
CARGO = "GSNKGAIIGLM"     # 벤치마크 기본 화물 (Aβ 25-35)
LINKERS = ["", "GGGGS", "AEAAAKEAAAKEAAAKA"]
LINKER_NAMES = ["직접융합", "GGGGS", "A(EAAAK)3A"]


def ii(seq: str) -> float:
    """① 과 같은 양을 같은 방법으로 — `developability_opt._instability` 와 동일 반올림.

    🔴 생산 함수를 **직접 쓴다**(복제하지 않는다). 복제하면 반올림이 갈려
    §36 계열의 "같은 이름 다른 양"이 하나 더 생긴다."""
    return _instability(seq)


def shuffles(seq: str, n: int, seed: int) -> list[str]:
    """조성 보존 순열 n 개. 원본과 같은 서열은 버리지 않고 그대로 센다 —
    버리면 '원본이 순열 중 하나일 확률'이 사라져 경험 p 가 편향된다."""
    rng = random.Random(seed)
    out = []
    chars = list(seq)
    for _ in range(n):
        c = chars[:]
        rng.shuffle(c)
        out.append("".join(c))
    return out


def main() -> None:
    print("=" * 78)
    print("§42 공백 ⓐ — 불안정성 축 계열 2(조성 보존 셔플) 대응 널")
    print(f"셔틀 12종 × 순열 {N_SHUFFLE}개 · 시드 {SEED} · 외부 호출 없음(결정론)")
    print("🔴 재는 양 = ① 셔틀 단독 II(탐색 축). 실제 안정성이 아니다.")
    print("=" * 78)

    # ── Q1·Q3: 셔틀 단독 ① ───────────────────────────────────────────────
    print("\n### Q1·Q3 셔틀 단독 ① — 문헌 셔틀 vs 같은 조성의 순열\n")
    print("| 셔틀 | ①II | 셔플 평균 | 셔플 최소 | 셔플 최대 | **셔플이 더 안정** | 경험 p | 고유순열 |")
    print("|---|---|---|---|---|---|---|---|")
    rows = []
    for name, meta in SHUTTLES.items():
        seq = meta["seq"]
        base = ii(seq)
        sh = shuffles(seq, N_SHUFFLE, SEED)
        vals = [ii(s) for s in sh]
        better = sum(1 for v in vals if v <= base)      # 더 '안정'(II 낮음)
        p = (better + 1) / (N_SHUFFLE + 1)
        uniq = len(set(sh))
        rows.append((name, seq, base, vals, better, p))
        flag = " 🔴" if better > N_SHUFFLE * 0.5 else ""
        print(f"| {name} | {base} | {statistics.mean(vals):.1f} | {min(vals):.1f} | "
              f"{max(vals):.1f} | **{better}/{N_SHUFFLE} ({100*better/N_SHUFFLE:.1f}%)**{flag} | "
              f"{p:.4f} | {uniq} |")

    med = statistics.median([r[5] for r in rows])
    n_sig = sum(1 for r in rows if r[5] < 0.05)
    n_half = sum(1 for r in rows if r[4] > N_SHUFFLE * 0.5)
    print(f"\n**경험 p 중위 {med:.3f} · 유의(p<0.05) {n_sig}/12 · "
          f"셔플 과반이 더 안정인 셔틀 {n_half}/12**")

    # ── Q2: 구성체 ③ 임계를 순열만으로 넘나 ──────────────────────────────
    print("\n### Q2 🔴 조성 편집 0회로 ③(구성체 II<40)을 넘나\n")
    print(f"화물 `{CARGO}` × 링커 3종. 원본이 **이미 넘은** 셀은 `—` 로 적는다.\n")
    print("| 셔틀 | 링커 | 원본 ②II | ③ | 순열 중 교차 | 최소 ②II(순열) |")
    print("|---|---|---|---|---|---|")
    cross_rows = []
    for name, meta in SHUTTLES.items():
        seq = meta["seq"]
        sh = shuffles(seq, N_SHUFFLE, SEED)
        for lk, lname in zip(LINKERS, LINKER_NAMES):
            base_c = ii(CARGO + lk + seq)
            if base_c < II_STABLE:
                continue  # 원본이 이미 안정 → Q2 의 질문 대상이 아니다
            vals = [ii(CARGO + lk + s) for s in sh]
            cross = sum(1 for v in vals if v < II_STABLE)
            cross_rows.append((name, lname, base_c, cross, min(vals)))
            mark = f"**{cross}/{N_SHUFFLE}**" + (" 🔴" if cross else "")
            print(f"| {name} | {lname} | {base_c} | 불안정 | {mark} | {min(vals):.1f} |")

    n_cells = len(cross_rows)
    n_crossed = sum(1 for r in cross_rows if r[3] > 0)
    print(f"\n**순열로 교차 가능한 셀 {n_crossed}/{n_cells}** "
          f"(원본이 이미 안정한 셀은 분모에서 제외)")

    # ── Q4: 제약이 그 경로를 잡나 ────────────────────────────────────────
    print("\n### Q4 제약(`rmt_sim ≥ 0.85`)이 '셔플로 안정화' 경로를 잡나\n")
    print("| 셔틀 | 원본 rmt_sim | ①개선 순열 | 그중 rmt_sim≥0.85 | 순열 rmt_sim 최대 |")
    print("|---|---|---|---|---|")
    total_pass = 0
    for name, seq, base, vals, better, _p in rows:
        sh = shuffles(seq, N_SHUFFLE, SEED)
        base_sim = shuttle_similarity(seq).score
        improved = [s for s, v in zip(sh, vals) if v < base]
        sims = [shuttle_similarity(s).score for s in improved] if improved else []
        passed = sum(1 for x in sims if x >= RMT_FLOOR)
        total_pass += passed
        mx = f"{max(sims):.3f}" if sims else "—"
        print(f"| {name} | {base_sim:.3f} | {len(improved)} | **{passed}** | {mx} |")
    print(f"\n**①을 개선한 순열 중 기능 제약을 통과한 것: 총 {total_pass}개**")

    print("\n" + "=" * 78)
    print("판정은 docs/experiments.md §46 에 적는다. 이 스크립트는 수치만 낸다.")
    print("=" * 78)


if __name__ == "__main__":
    main()
