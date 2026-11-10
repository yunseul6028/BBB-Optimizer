"""개발성(developability) 평가 — 서열 기반 liability·응집·용해도·전하.

계산이 유효한 영역: 항체·펩타이드 치료제 개발에서 **표준으로 쓰는 서열 수준 위험 신호**.
실험 전에 "만들기 어렵다/불안정하다/비특이 결합·독성 위험/잘 안 녹는다"를 미리 걸러낸다.
100% 서열 기반이라 즉시 계산 — LLM·실험 불필요.

평가 항목:
  1) 서열 liability 모티프
     - 탈아마이드화(deamidation): N-[G/S/T/N] (NG가 최고 위험)
     - 이성질화(isomerization): D-[G/S/T/D/H]
     - N-당화 sequon: N-X-[S/T] (X≠P)
     - 자유 시스테인(free Cys): 이황화 스크램블·응집 위험
     - 산화 취약(oxidation): Met·Trp
     - 프로테아제 절단부위(dibasic): RR/KR/RK/KK (furin·트립신) — 체내 안정성↓
  2) 응집 경향(aggregation): **AGGRESCAN a3v** 아미노산별 응집 성향 척도
     (Conchillo-Solé et al., BMC Bioinformatics 2007; Sánchez de Groot 2005)
     — in vivo 실험 유래, **펩타이드에 적합**.
  3) 용해도(solubility): 친수성(GRAVY)·전하·(1−응집)의 가중합 **상대지표**(0~1).
  4) 전하: pH 7.4 순전하·전하밀도·pI — **과도한 양전하 = 비특이 결합·빠른 청소·독성 위험**
     (BBB 셔틀은 양이온성이 많아 특히 중요)

※ 응집·용해도는 예전 **별도 '용해도 축'**이었으나, 같은 물성(응집)을 두 축이 **중복 신호**하던
  것을 정리해 **개발성 축으로 통합**했다. 응집은 이제 거친 휴리스틱이 아니라 출판 AGGRESCAN
  a3v 척도 **하나로만**(단일 출처) 계산하고, 용해도 점수는 개발성의 보조 지표로 함께 보고한다.

⚠️ 서열 규칙 기반 경험 지표다(TAP/Therapeutic Antibody Profiler 개념의 펩타이드판).
   용해도 가중치(0.6/0.2/0.2)·정규화 상수는 랭킹에 쓰지 않는 **표시·보조용 미피팅 설계값**이다.
   실제 개발성·실측 용해도(mg/mL)는 발현·정제·가속 안정성 실험으로 확정.
   (출판 닫힌 식 Wilkinson-Harrison(1991)은 재조합 전체 단백질 유도라 짧은 펩타이드에서 비상식적
    결과 → 기각. 펩타이드 전용 CamSol-PTM(2023)은 웹서버 의존이라 로컬 미채택 — 후향검증 레퍼런스.)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

_STD = set("ACDEFGHIKLMNPQRSTVWY")

# AGGRESCAN a3v — 아미노산별 내재 응집 성향(높을수록 응집↑). 출판 척도(자유 파라미터 아님).
_AGGRESCAN_A3V = {
    "I": 1.822, "F": 1.754, "V": 1.594, "L": 1.380, "Y": 1.159, "W": 1.037,
    "M": 0.910, "C": 0.604, "A": -0.036, "T": -0.159, "S": -0.294, "P": -0.334,
    "G": -0.535, "K": -0.931, "H": -1.033, "Q": -1.231, "R": -1.240, "N": -1.302,
    "E": -1.412, "D": -1.836,
}
_A3V_MIN, _A3V_MAX = -1.836, 1.822     # 척도 자체의 극값(정규화용 — 설계가 아니라 척도 고유)

# 응집 **핫스팟** 창 길이(aa). AGGRESCAN 이 응집 경향을 다루는 방식은 전역 평균이 아니라
# **연속 구간(hot spot)** 이다. 5 는 AGGRESCAN 의 a4v/HSA 창 길이를 따른다(자유 파라미터 아님).
#
# ⚠️ 왜 전역 평균을 순위 축에서 뺐나 (M10 실측, docs/experiments.md §20)
# 전역 평균 a3v 는 **길이로 희석된다** — 저응집 잔기(G·S·E·K·D)를 붙이면 분모만 커져
# 점수가 내려간다. 실측: 후보군 셔틀 8종 전부에서 링커 길이↔응집 점수 순위상관이 음수
# (평균 **-0.697**, 8/8)이고, 그 결과 도달 가능 전선이 **가장 긴 링커 A(EAAAK)3A 로
# 전원 수렴**했다. 이것은 M8 이 BBB 에서 찾은 길이 아티팩트와 **구조가 같다** — 지표를
# 길이로 정규화해 두면 "무해한 잔기를 붙이기"가 항상 이긴다.
# 핫스팟(국소 최대)은 이 희석에 면역이다: 링커를 끼워도 셔틀·화물의 국소 최대는 그대로다.
_AGG_WINDOW = 5


@dataclass
class DevelopabilityResult:
    liabilities: list = field(default_factory=list)  # 사람이 읽는 위험 신호
    n_liabilities: int = 0
    free_cys: int = 0
    net_charge: float = 0.0        # pH 7.4 순전하
    charge_density: float = 0.0    # 순전하 / 길이
    pI: float = 0.0
    agg_a3v: float = 0.0           # 평균 AGGRESCAN a3v (raw, 높을수록 응집↑)
    agg_score: float = 0.0         # 0~1 정규화 **전역 평균** 응집 — ⚠️ 보고용(길이로 희석됨)
    agg_hotspot_a3v: float = 0.0   # 창 _AGG_WINDOW 슬라이딩 **최대** a3v (raw)
    agg_hotspot: float = 0.0       # 0~1 정규화 국소 최대 — **순위·Pareto 축은 이것**
    gravy: float = 0.0             # 소수성(양수=소수성)
    sol_score: float = 0.5         # 0~1 용해도 상대지표(높을수록 가용)
    sol_level: str = "보통"         # 낮음 | 보통 | 높음
    risk_level: str = "낮음"        # 낮음 | 보통 | 높음
    verdict: str = ""
    error: str = ""


def _clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def _net_charge_ph74(seq: str) -> float:
    pos = seq.count("K") + seq.count("R") + 0.1 * seq.count("H")  # His ~10% 양성@7.4
    neg = seq.count("D") + seq.count("E")
    return round(pos - neg, 1)


def assess_developability(sequence: str) -> DevelopabilityResult:
    seq = "".join(ch for ch in (sequence or "").upper() if ch in _STD)
    if len(seq) < 4:
        return DevelopabilityResult(error="서열이 너무 짧아 개발성 계산 불가")

    liab = []
    # 1) 모티프 liability
    for m in re.finditer(r"N[GSTN]", seq):
        liab.append(f"탈아마이드화 {m.group()}@{m.start()+1}")
    for m in re.finditer(r"D[GSTDH]", seq):
        liab.append(f"이성질화 {m.group()}@{m.start()+1}")
    for m in re.finditer(r"N[^P][ST]", seq):
        liab.append(f"N-당화 sequon {m.group()}@{m.start()+1}")
    for m in re.finditer(r"(RR|KR|RK|KK)", seq):
        liab.append(f"프로테아제 절단부위 {m.group()}@{m.start()+1}")
    free_cys = seq.count("C")
    if free_cys:
        liab.append(f"자유 시스테인 {free_cys}개(이황화 스크램블 위험)")
    n_ox = seq.count("M") + seq.count("W")
    if n_ox >= 2:
        liab.append(f"산화 취약 잔기(Met/Trp) {n_ox}개")

    # 2) 전하
    net = _net_charge_ph74(seq)
    density = round(net / len(seq), 2)
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        gravy = round(float(ProteinAnalysis(seq).gravy()), 2)
        pI = round(float(ProteinAnalysis(seq).isoelectric_point()), 1)
    except Exception:  # noqa: BLE001
        gravy, pI = 0.0, 0.0
    if density >= 0.30:
        liab.append(f"과도한 양전하(밀도 {density}) — 비특이 결합·독성·빠른 청소 위험")

    # 3) 응집 경향 — AGGRESCAN a3v
    #    · agg_score  = 전역 평균 (보고용. **순위에 쓰지 말 것** — 길이로 희석된다)
    #    · agg_hotspot = 길이 _AGG_WINDOW 슬라이딩 창의 **최대** 평균 (순위 축)
    vals = [_AGGRESCAN_A3V.get(c, 0.0) for c in seq]
    mean_a3v = sum(vals) / len(vals)
    agg = round(_clamp((mean_a3v - _A3V_MIN) / (_A3V_MAX - _A3V_MIN)), 2)
    w = min(_AGG_WINDOW, len(vals))
    run = sum(vals[:w])
    hot_raw = run / w
    for i in range(w, len(vals)):
        run += vals[i] - vals[i - w]
        hot_raw = max(hot_raw, run / w)
    hot = round(_clamp((hot_raw - _A3V_MIN) / (_A3V_MAX - _A3V_MIN)), 2)
    # liability 판정은 **핫스팟**으로 한다 — 응집은 국소 현상이고, 전역 평균으로 보면
    # 긴 구성체의 실제 응집 위험 구간이 희석돼 사라진다.
    if hot >= 0.70:
        liab.append(f"응집 핫스팟(AGGRESCAN a3v 국소최대 {hot_raw:.2f}, 창 {w}aa)")

    # 4) 용해도 상대지표(0~1) — 친수성 + 전하보너스 + (1−응집)
    hydrophil = _clamp((-gravy + 0.8) / 1.5)      # 친수성↑ = 용해도↑
    charge_bonus = _clamp(abs(density) / 0.3)     # 하전 잔기 = 용해도↑(보너스)
    sol_score = round(_clamp(0.6 * hydrophil + 0.2 * charge_bonus + 0.2 * (1 - agg)), 2)
    sol_level = "높음" if sol_score >= 0.6 else "보통" if sol_score >= 0.35 else "낮음"

    r = DevelopabilityResult(
        liabilities=liab, n_liabilities=len(liab), free_cys=free_cys,
        net_charge=net, charge_density=density, pI=pI,
        agg_a3v=round(mean_a3v, 3), agg_score=agg,
        agg_hotspot_a3v=round(hot_raw, 3), agg_hotspot=hot,
        gravy=gravy, sol_score=sol_score, sol_level=sol_level,
    )
    # 위험 등급: liability 수 + 강한 신호 가중
    strong = sum(1 for x in liab if "과도한 양전하" in x or "응집 핫스팟" in x or "자유 시스테인" in x)
    if len(liab) >= 4 or strong >= 2:
        r.risk_level = "높음"
    elif len(liab) >= 2 or strong >= 1:
        r.risk_level = "보통"
    else:
        r.risk_level = "낮음"
    r.verdict = (f"개발성 위험 {r.risk_level} — liability {len(liab)}개, 순전하 {net}"
                 f"(밀도 {density}), 응집 a3v {mean_a3v:.2f}, 용해도 {sol_level}({sol_score})")
    return r
