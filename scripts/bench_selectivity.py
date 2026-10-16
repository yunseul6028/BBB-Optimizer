"""선택성 축 vs deepB3P 단독 — RMT vs CPP 판별력 재현 벤치 (docs/experiments.md §1).

우리 메커니즘(선택성 축, 규칙 기반)이 deepB3P 단독보다 RMT 리간드와 CPP를 잘 분리하는지를
AUC(RMT>CPP)로 비교한다. 보상 해킹(deepB3P가 양이온 CPP를 과대평가) 실증.

평가셋:
  RMT = 라이브러리 RMT 셔틀 전체(note에 "RMT")
  CPP = 라이브러리 CPP + 홀드아웃(R9·Transportan·MAP·pVEC) — 라이브러리 밖 = 비순환성 검증

실행: 프로젝트 루트에서 Biopython이 설치된 파이썬으로 (deepB3P는 vendor/.venv-deepb3p subprocess).
  $ python scripts/bench_selectivity.py
결과(시드 고정·결정론): 선택성 AUC 1.0 · deepB3P AUC ~0.40 (둘 다 <0.5 = CPP를 더 높게 = 보상해킹).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import SHUTTLES
from core.selectivity import assess_selectivity
from core.predictors import DeepB3PLocalPredictor

HOLDOUT_CPP = {               # 라이브러리 밖 = 비순환성(암기 아님) 검증용
    "R9":          "RRRRRRRRR",
    "Transportan": "GWTLNSAGYLLGKINLKALAALAKKIL",
    "MAP":         "KLALKLALKALKAALKLA",
    "pVEC":        "LLIILRRRIRKQAHAHSK",
}


def auc(pos, neg) -> float:
    """순위 기반 AUC = P(pos > neg) (Mann-Whitney)."""
    c = t = 0
    for p in pos:
        for n in neg:
            t += 1
            c += 1 if p > n else (0.5 if p == n else 0)
    return c / t if t else 0.5


def main() -> None:
    rmt = {n: SHUTTLES[n]["seq"] for n in SHUTTLES if "RMT" in SHUTTLES[n]["note"]}
    lib_cpp = {n: SHUTTLES[n]["seq"] for n in SHUTTLES if "CPP" in SHUTTLES[n]["note"]}
    cpp = {**lib_cpp, **HOLDOUT_CPP}
    print(f"RMT {len(rmt)}종: {list(rmt)}")
    print(f"CPP {len(cpp)}종: 라이브러리 {list(lib_cpp)} + 홀드아웃 {list(HOLDOUT_CPP)}\n")

    names = list(rmt) + list(cpp)
    seqs = [rmt[n] for n in rmt] + [cpp[n] for n in cpp]
    sel = [assess_selectivity(s).selectivity for s in seqs]
    bbb = [r.bbb_permeability for r in DeepB3PLocalPredictor().predict_many(seqs)]

    nR = len(rmt)
    print(f"{'셔틀':<16}{'군':>5}{'선택성':>9}{'deepB3P':>10}")
    for i, n in enumerate(names):
        print(f"{n:<16}{'RMT' if i < nR else 'CPP':>5}{sel[i]:>9.3f}{bbb[i]:>10.3f}")

    print("\n=== 요약 ===")
    print(f"{'축':<12}{'RMT평균':>9}{'CPP평균':>9}{'AUC(RMT>CPP)':>14}")
    print(f"{'선택성':<12}{sum(sel[:nR])/nR:>9.3f}{sum(sel[nR:])/(len(sel)-nR):>9.3f}"
          f"{auc(sel[:nR], sel[nR:]):>14.3f}")
    print(f"{'deepB3P':<12}{sum(bbb[:nR])/nR:>9.3f}{sum(bbb[nR:])/(len(bbb)-nR):>9.3f}"
          f"{auc(bbb[:nR], bbb[nR:]):>14.3f}")


if __name__ == "__main__":
    main()
