"""링커 축의 **대응 널** — 라이브러리 링커가 같은 길이 무작위 링커를 이기나 (M8 · M10 재현).

`docs/protocol.md` 계열 1(대응 무작위 편집)의 **구성체 수준 변종**이다. 편집 크기를
"같은 길이"로 고정하고, 라이브러리가 고른 링커를 **같은 길이 무작위 링커**와 겨룬다.

왜 두 지표를 한 스크립트에 두나
------------------------------
같은 널 설계에 **지표만 바꿔** 넣으면 "축이 떨어지는가"가 축의 문제인지 널의 문제인지
갈린다. 실제로 그렇게 갈렸다:

  · `--metric bbb`  (M8, §17) → 라이브러리가 이긴 비율 **0/108 (0.0%)**, 경험 p 중위 ~0.86
                                = **거꾸로**. BBB 로는 링커를 고를 수 없다.
  · `--metric risk` (M10, §20) → **41/96 (42.7%)**, 무작위가 이긴 비율 중위 0.03~0.10
                                = **통과**. 링커 선택의 근거는 개발성이었다.

즉 **같은 널에서 한 축은 떨어지고 다른 축은 통과한다** — 널이 축을 구별한다는 직접 증거다.

실행 (프로젝트 루트, Biopython 설치된 파이썬):
  $ /usr/bin/python3 scripts/null_linker.py --metric risk          # M10 (로컬·무료·빠름)
  $ /usr/bin/python3 scripts/null_linker.py --metric bbb           # M8 (deepB3P 서브프로세스)

⚠️ 시드는 **반드시 기록**한다 — 옛 M8·M10 측정은 시드가 기록되지 않아 재현 불가였다.
"""
from __future__ import annotations

import argparse
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import DEFAULT_CARGO, LINKER_LIBRARY, SHUTTLES
from core.developability import assess_developability
from core.evidence import group_rank, inherit
from core.optimizer_agent import risk_axes
from core.selectivity import assess_selectivity
from core.stability import assess_stability

_STD = "ACDEFGHIKLMNPQRSTVWY"


def _seq(d) -> str:
    return d["seq"] if isinstance(d, dict) else d


def linkers(include_direct: bool = False) -> dict[str, str]:
    """널에 쓰는 링커 집합.

    ⚠️ **직접융합(링커 없음)은 제외**한다 — 길이 0 의 무작위 링커를 만들 수 없어 대응 널이
    성립하지 않는다. 그래서 라이브러리 13종 중 **12종**이 분모다(옛 기록의 '13종'은 오기).
    """
    out = {n: _seq(d) for n, d in LINKER_LIBRARY.items()}
    if not include_direct:
        out = {n: s for n, s in out.items() if s}
    return out


def shuttles(metric: str) -> dict[str, str]:
    """셔틀 집합 — 지표에 따라 다르다(옛 측정의 설계를 그대로 따른다).

    · `bbb`  (M8)  = **RMT 셔틀 9종** 전부
    · `risk` (M10) = **후보군 8종** — RMT 9종에서 **근거 등급 D** 를 뺀 것.
      위험 축은 전부 감점 방향이라 근거 그룹이 없는 후보를 섞으면 축이 의미를 잃는다.
    """
    rmt = {n: _seq(d) for n, d in SHUTTLES.items() if "RMT" in SHUTTLES[n].get("note", "")}
    if metric == "bbb":
        return rmt
    return {n: s for n, s in rmt.items() if group_rank(inherit(s)[0]) < group_rank("D")}


# --------------------------------------------------------------------------
def risk_vector(cargo: str, linker: str, shuttle: str) -> tuple:
    """순위 위험 축 4개 — `optimizer_agent._evaluate` 와 **같은 범위 선택**으로 계산.

    deepB3P 를 **쓰지 않는다**(위험 축에 BBB 가 없다) → 서브프로세스 0회·전부 결정론.
    """
    full = cargo + linker + shuttle
    mod = linker + shuttle
    return risk_axes({
        "dev_liab": assess_developability(full).n_liabilities,
        "instability": assess_stability(full).instability_index,
        "sel_off": assess_selectivity(shuttle).off_target_risk,
        "dev_agg_hot": assess_developability(mod).agg_hotspot if mod else 0.0,
    })


def _dominates(a: tuple, b: tuple) -> bool:
    return all(x <= y for x, y in zip(a, b)) and a != b


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metric", choices=["risk", "bbb"], default="risk",
                    help="risk = M10(위험 축, 로컬) · bbb = M8(deepB3P)")
    ap.add_argument("--nulls", type=int, default=None,
                    help="조합당 무작위 링커 개수 (기본: risk 30 · bbb 20 — 옛 측정 설계)")
    ap.add_argument("--seed", type=int, default=0, help="🔴 반드시 기록할 난수 시드")
    ap.add_argument("--cargo", default=DEFAULT_CARGO)
    a = ap.parse_args()
    n_null = a.nulls if a.nulls is not None else (30 if a.metric == "risk" else 20)
    rng = random.Random(a.seed)
    LK, SH = linkers(), shuttles(a.metric)

    print(f"지표 {a.metric} · 셔틀 {len(SH)}종 × 링커 {len(LK)}종 = {len(SH) * len(LK)}조합 "
          f"· 널 {n_null}개/조합 · 시드 {a.seed} · 화물 {a.cargo}")

    pred = None
    if a.metric == "bbb":
        from core.predictors import DeepB3PLocalPredictor
        from core.config import bbb_scoring_seq
        pred = DeepB3PLocalPredictor()

    rows, wins, pvals = [], 0, []
    for sname, sseq in SH.items():
        for lname, lseq in LK.items():
            nulls = ["".join(rng.choice(_STD) for _ in range(len(lseq)))
                     for _ in range(n_null)]
            if a.metric == "risk":
                real = risk_vector(a.cargo, lseq, sseq)
                nv = [risk_vector(a.cargo, x, sseq) for x in nulls]
                # 경험 p = 널이 **실제 링커를 지배한** 비율 (작을수록 라이브러리 우세)
                beat = sum(1 for v in nv if _dominates(v, real))
            else:
                seqs = [bbb_scoring_seq(a.cargo, x, sseq) for x in [lseq] + nulls]
                sc = [r.bbb_permeability for r in pred.predict_many(seqs)]
                real, nv = sc[0], sc[1:]
                beat = sum(1 for v in nv if v > real)
            p = beat / n_null
            pvals.append(p)
            sig = p < 0.05
            wins += sig
            rows.append((sname, lname, len(lseq), p, sig))

    print(f"\n{'셔틀':<16}{'링커':<14}{'len':>4}{'경험 p':>9}{'유의':>6}")
    for sname, lname, ln, p, sig in rows:
        print(f"{sname:<16}{lname:<14}{ln:>4}{p:>9.2f}{'✓' if sig else '':>6}")

    print(f"\n=== 요약 (지표 {a.metric} · 시드 {a.seed}) ===")
    print(f"라이브러리 링커가 유의하게 이긴 조합: **{wins}/{len(rows)} "
          f"({wins / len(rows) * 100:.1f}%)**   (우연 기대 ~5%)")
    print(f"경험 p 중위: {statistics.median(pvals):.2f} · 평균 {statistics.mean(pvals):.2f}")
    print("해석: p 가 0 에 가까우면 라이브러리 우세, 0.5 이상이면 무작위와 구별 불가, "
          "1 에 가까우면 **거꾸로**(무작위가 더 좋다).")


if __name__ == "__main__":
    main()
