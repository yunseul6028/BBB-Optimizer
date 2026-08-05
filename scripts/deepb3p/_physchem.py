# -- coding: utf-8 --
"""생성 루프용 물성 가드레일 — 양전하(양이온) 편향 억제.

deepB3P 는 양전하 펩타이드를 BBB 고투과로 과적합하는 경향이 있다(TAT/Penetratin 등).
생성 적합도에서 BBB 점수만 최대화하면 생성기가 이 편향을 착취해 '무지성 양이온' 덩어리를
뽑는다. 여기서 순전하(pH 7.4)와 Eisenberg 소수성 모멘트(μH, 양친매성)를 결합한 페널티를
적합도에 감점항으로 더해 이를 상쇄한다.

  penalty = excess_pos_charge · (W_CHARGE + W_AMPHI · max(0, AMPHI_REF − μH))
    · excess_pos_charge = max(0, netcharge − CHARGE_TOL)   # 허용치 초과 '양전하'만 감점
    · 양친매성(μH)이 낮을수록(=무지성 양이온) 같은 전하라도 더 크게 감점

순전하 공식은 core/developability._net_charge_ph74 와 동일(앱 전체 일관성).
"""
import math

# --- 튜닝 파라미터 (편향 억제 강도) ---
CHARGE_TOL = 3.0     # 허용 순전하(+3까지는 RMT 셔틀에 흔함, 초과분만 감점)
W_CHARGE = 0.04      # 초과 전하 단위당 기본 감점
W_AMPHI = 0.15       # 비양친매성 가중(무지성 양이온 추가 감점)
AMPHI_REF = 0.35     # 이 μH 미만이면 '비양친매' 취급

# Eisenberg consensus 소수성 스케일
_EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "Q": -0.85,
    "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50,
    "M": 0.64, "F": 1.19, "P": 0.12, "S": -0.18, "T": -0.05, "W": 0.81,
    "Y": 0.26, "V": 1.08,
}


def net_charge_ph74(seq: str) -> float:
    """pH 7.4 순전하 — K,R(+1)·H(+0.1) − D,E(−1). core/developability 와 동일."""
    pos = seq.count("K") + seq.count("R") + 0.1 * seq.count("H")
    neg = seq.count("D") + seq.count("E")
    return round(pos - neg, 1)


def hydrophobic_moment(seq: str, angle: float = 100.0) -> float:
    """Eisenberg 소수성 모멘트 μH (α-helix δ=100°), 잔기수로 정규화."""
    s = [c for c in seq if c in _EISENBERG]
    if not s:
        return 0.0
    rad = math.radians(angle)
    sx = sum(_EISENBERG[c] * math.sin(i * rad) for i, c in enumerate(s))
    sy = sum(_EISENBERG[c] * math.cos(i * rad) for i, c in enumerate(s))
    return math.sqrt(sx * sx + sy * sy) / len(s)


def charge_guardrail(seq: str):
    """(penalty, net_charge, muH) — 적합도 감점항 + 표시용 물성."""
    q = net_charge_ph74(seq)
    excess = max(0.0, q - CHARGE_TOL)
    muH = hydrophobic_moment(seq)
    penalty = excess * (W_CHARGE + W_AMPHI * max(0.0, AMPHI_REF - muH))
    return round(penalty, 4), q, round(muH, 3)
