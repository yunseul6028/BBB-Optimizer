"""용해도(solubility) — 서열 기반 경량 근사.

용해도에 유리한 요인: (1) 친수성↑(소수성 GRAVY↓), (2) 순전하↑(하전 잔기가 물과 친화),
(3) 응집 경향↓. CamSol/Protein-Sol이 쓰는 원리(소수성+전하+패터닝)의 경량 휴리스틱.

⚠️ 실측 용해도(mg/mL)가 아니라 **상대 지표**(0~1, 높을수록 잘 녹음). 개발성 판단의 보조 축.
"""

from __future__ import annotations

from dataclasses import dataclass

_STD = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class SolubilityResult:
    score: float = 0.5          # 0~1 (높을수록 용해도 좋음)
    gravy: float = 0.0          # 소수성(양수=소수성)
    charge_density: float = 0.0
    level: str = "보통"          # 낮음 | 보통 | 높음
    verdict: str = ""
    error: str = ""


def _clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def assess_solubility(sequence: str) -> SolubilityResult:
    """융합체 서열의 용해도 근사 점수(0~1)."""
    seq = "".join(c for c in (sequence or "").upper() if c in _STD)
    if len(seq) < 5:
        return SolubilityResult(error="서열이 너무 짧아 용해도 계산 불가")

    from .developability import _HYDROPHOBIC, _longest_hydrophobic_run, _net_charge_ph74
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        gravy = round(float(ProteinAnalysis(seq).gravy()), 2)
    except Exception:  # noqa: BLE001
        gravy = 0.0
    dens = round(_net_charge_ph74(seq) / len(seq), 2)
    hp_frac = sum(1 for c in seq if c in _HYDROPHOBIC) / len(seq)
    agg = _clamp(0.7 * hp_frac + 0.1 * _longest_hydrophobic_run(seq))

    hydrophil = _clamp((-gravy + 0.8) / 1.5)      # 친수성↑ = 용해도↑
    charge_bonus = _clamp(abs(dens) / 0.3)        # 하전 잔기 = 용해도↑(보너스)
    score = round(_clamp(0.6 * hydrophil + 0.2 * charge_bonus + 0.2 * (1 - agg)), 2)

    level = "높음" if score >= 0.6 else "보통" if score >= 0.35 else "낮음"
    r = SolubilityResult(score=score, gravy=gravy, charge_density=dens, level=level)
    r.verdict = f"용해도 {level}(점수 {score}) — GRAVY {gravy}, 순전하밀도 {dens}"
    return r
