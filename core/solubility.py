"""용해도(solubility) — 서열 기반 경량 지표 (출판 척도 + 설계 가중치).

용해도에 유리한 요인: (1) 친수성↑(GRAVY↓), (2) 순전하↑(하전 잔기가 물과 친화),
(3) 응집 경향↓. 세 항의 가중합(0~1). 각 항의 **계산은 출판된 아미노산 척도**를 쓴다:

  · 친수성  : Kyte-Doolittle GRAVY (ProtParam)
  · 전하    : pH 7.4 순전하 밀도
  · 응집    : **AGGRESCAN a3v** 아미노산별 응집 성향 척도
             (Conchillo-Solé et al., BMC Bioinformatics 2007; Sánchez de Groot 2005)
             — in vivo 실험 유래, **펩타이드에 적합**(전체 단백질 전용 아님).

가중치(0.6/0.2/0.2)와 GRAVY·전하의 [0,1] 정규화 상수만 **미피팅 설계값**이다. 용해도는
랭킹(유효점수)에 쓰지 않는 **표시·보조 축**이라 엄밀 피팅보다 가벼움·투명성을 택했다.

※ 출판 닫힌 식 Wilkinson-Harrison(1991)도 검토했으나, 대장균 재조합 **전체 단백질**로
  유도되어 짧은 펩타이드에서 비상식적 결과(GS 링커→불용성, 순수 소수성→가용성; 소수성
  항 부재·turn 논리의 단백질 특이성)를 내어 **기각**했다. 대신 응집항을 **펩타이드 유래
  AGGRESCAN 척도**로 두어 모달리티를 맞췄다. (펩타이드 전용 CamSol-PTM(2023)은 웹서버
  의존이라 로컬 엔진엔 미채택 — 향후 후향 검증 레퍼런스로 활용 가능.)

⚠️ 실측 용해도(mg/mL)가 아니라 **상대 지표**(0~1, 높을수록 잘 녹음). 개발성 판단의 보조 축.
"""

from __future__ import annotations

from dataclasses import dataclass

_STD = set("ACDEFGHIKLMNPQRSTVWY")

# AGGRESCAN a3v — 아미노산별 내재 응집 성향(높을수록 응집↑). 출판 척도(자유 파라미터 아님).
_AGGRESCAN_A3V = {
    "I": 1.822, "F": 1.754, "V": 1.594, "L": 1.380, "Y": 1.159, "W": 1.037,
    "M": 0.910, "C": 0.604, "A": -0.036, "T": -0.159, "S": -0.294, "P": -0.334,
    "G": -0.535, "K": -0.931, "H": -1.033, "Q": -1.231, "R": -1.240, "N": -1.302,
    "E": -1.412, "D": -1.836,
}
_A3V_MIN, _A3V_MAX = -1.836, 1.822     # 척도 자체의 극값(정규화용 — 설계가 아니라 척도 고유)


@dataclass
class SolubilityResult:
    score: float = 0.5          # 0~1 (높을수록 용해도 좋음)
    gravy: float = 0.0          # 소수성(양수=소수성)
    charge_density: float = 0.0
    agg_a3v: float = 0.0        # 평균 AGGRESCAN a3v (높을수록 응집↑)
    level: str = "보통"          # 낮음 | 보통 | 높음
    verdict: str = ""
    error: str = ""


def _clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def assess_solubility(sequence: str) -> SolubilityResult:
    """융합체 서열의 용해도 근사 점수(0~1). 응집항은 AGGRESCAN a3v 척도 기반."""
    seq = "".join(c for c in (sequence or "").upper() if c in _STD)
    if len(seq) < 5:
        return SolubilityResult(error="서열이 너무 짧아 용해도 계산 불가")

    from .developability import _net_charge_ph74
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        gravy = round(float(ProteinAnalysis(seq).gravy()), 2)
    except Exception:  # noqa: BLE001
        gravy = 0.0
    dens = round(_net_charge_ph74(seq) / len(seq), 2)

    # 응집: AGGRESCAN a3v 평균 → 척도 극값으로 [0,1] 정규화(높을수록 응집↑)
    mean_a3v = sum(_AGGRESCAN_A3V.get(c, 0.0) for c in seq) / len(seq)
    agg = _clamp((mean_a3v - _A3V_MIN) / (_A3V_MAX - _A3V_MIN))

    hydrophil = _clamp((-gravy + 0.8) / 1.5)      # 친수성↑ = 용해도↑
    charge_bonus = _clamp(abs(dens) / 0.3)        # 하전 잔기 = 용해도↑(보너스)
    score = round(_clamp(0.6 * hydrophil + 0.2 * charge_bonus + 0.2 * (1 - agg)), 2)

    level = "높음" if score >= 0.6 else "보통" if score >= 0.35 else "낮음"
    r = SolubilityResult(score=score, gravy=gravy, charge_density=dens,
                         agg_a3v=round(mean_a3v, 3), level=level)
    r.verdict = (f"용해도 {level}(점수 {score}) — GRAVY {gravy}, 순전하밀도 {dens}, "
                 f"AGGRESCAN a3v {round(mean_a3v, 2)}")
    return r
