"""융합체 자체 안정성 — Biopython ProtParam 기반 서열 지표.

'결합(융합체) 자체'가 분해·응집 없이 얼마나 안정한지를 **서열 조성만으로** 근사한다.
deepB3P(투과)·ToxinPred3(독성)가 못 보는 축 — 이 분자가 시험관/체내에서 형태를
유지할 가능성 — 을 보완한다.

주요 지표:
  - 불안정성 지수(instability index, Guruprasad 1990): **< 40 이면 안정**으로 예측.
  - 지방족 지수(aliphatic index, Ikai 1980): 높을수록 열안정성↑.
  - GRAVY(평균 소수성): 양수=소수성(응집 경향↑), 음수=친수성.
  - 방향족성(aromaticity).

⚠️ in-vitro 실측(Tm·ΔG·응집 실험)이 아니라 **조성 기반 경험적 지표**다. 상대 비교·
   스크리닝용이며, 실제 개발 전 실험 검증이 필요하다.
"""

from __future__ import annotations

from dataclasses import dataclass

_STD = set("ACDEFGHIKLMNPQRSTVWY")  # ProtParam은 표준 20종 아미노산만 처리


@dataclass
class StabilityResult:
    instability_index: float = 0.0  # < 40 = 안정 예측
    aliphatic_index: float = 0.0    # 높을수록 열안정
    gravy: float = 0.0              # 평균 소수성(양수=소수성)
    aromaticity: float = 0.0
    stable: bool = False
    verdict: str = ""
    error: str = ""


def assess_stability(sequence: str) -> StabilityResult:
    """융합체 전체 서열의 안정성 지표를 계산한다."""
    seq = "".join(ch for ch in (sequence or "").upper() if ch in _STD)
    if len(seq) < 5:
        return StabilityResult(error="서열이 너무 짧아 안정성 계산 불가")
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis

        pa = ProteinAnalysis(seq)
        ii = float(pa.instability_index())
        gravy = float(pa.gravy())
        arom = float(pa.aromaticity())
        perc = pa.amino_acids_percent  # {AA: 0~1 fraction}
    except Exception as exc:  # noqa: BLE001
        return StabilityResult(error=f"안정성 계산 오류: {exc}")

    # 지방족 지수 (Ikai 1980): mole% 기준 상대 부피
    ai = (perc.get("A", 0) * 100
          + 2.9 * perc.get("V", 0) * 100
          + 3.9 * (perc.get("I", 0) + perc.get("L", 0)) * 100)

    stable = ii < 40.0
    r = StabilityResult(
        instability_index=round(ii, 1),
        aliphatic_index=round(ai, 1),
        gravy=round(gravy, 2),
        aromaticity=round(arom, 3),
        stable=stable,
    )
    if ii < 40.0:
        r.verdict = f"안정 예측 (불안정성 지수 {ii:.1f} < 40)"
    else:
        extra = " (특히 높음)" if ii >= 55 else ""
        r.verdict = f"불안정 경향 (불안정성 지수 {ii:.1f} ≥ 40){extra}"
    return r
