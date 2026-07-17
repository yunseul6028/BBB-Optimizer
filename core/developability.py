"""개발성(developability) 평가 — 서열 기반 liability·응집·전하.

계산이 유효한 영역: 항체·펩타이드 치료제 개발에서 **표준으로 쓰는 서열 수준 위험 신호**.
실험 전에 "만들기 어렵다/불안정하다/비특이 결합·독성 위험"을 미리 걸러낸다.
100% 서열 기반이라 즉시 계산 — LLM·실험 불필요.

평가 항목:
  1) 서열 liability 모티프
     - 탈아마이드화(deamidation): N-[G/S/T/N] (NG가 최고 위험)
     - 이성질화(isomerization): D-[G/S/T/D/H]
     - N-당화 sequon: N-X-[S/T] (X≠P)
     - 자유 시스테인(free Cys): 이황화 스크램블·응집 위험
     - 산화 취약(oxidation): Met·Trp
     - 프로테아제 절단부위(dibasic): RR/KR/RK/KK (furin·트립신) — 체내 안정성↓
  2) 응집 경향(aggregation): 소수성·방향족 비율 + 최장 소수성 연쇄 (휴리스틱)
  3) 전하: pH 7.4 순전하·전하밀도·pI — **과도한 양전하 = 비특이 결합·빠른 청소·독성 위험**
     (BBB 셔틀은 양이온성이 많아 특히 중요)

⚠️ 서열 규칙 기반 경험 지표다(TAP/Therapeutic Antibody Profiler 개념의 펩타이드판).
   실제 개발성은 발현·정제·가속 안정성 실험으로 확정.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

_STD = set("ACDEFGHIKLMNPQRSTVWY")
_HYDROPHOBIC = set("FILVWYM")  # 응집 유발 소수성·방향족


@dataclass
class DevelopabilityResult:
    liabilities: list = field(default_factory=list)  # 사람이 읽는 위험 신호
    n_liabilities: int = 0
    free_cys: int = 0
    net_charge: float = 0.0        # pH 7.4 순전하
    charge_density: float = 0.0    # 순전하 / 길이
    pI: float = 0.0
    agg_score: float = 0.0         # 0~1 응집 경향(휴리스틱)
    risk_level: str = "낮음"        # 낮음 | 보통 | 높음
    verdict: str = ""
    error: str = ""


def _net_charge_ph74(seq: str) -> float:
    pos = seq.count("K") + seq.count("R") + 0.1 * seq.count("H")  # His ~10% 양성@7.4
    neg = seq.count("D") + seq.count("E")
    return round(pos - neg, 1)


def _longest_hydrophobic_run(seq: str) -> int:
    best = cur = 0
    for ch in seq:
        cur = cur + 1 if ch in _HYDROPHOBIC else 0
        best = max(best, cur)
    return best


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
        pI = round(float(ProteinAnalysis(seq).isoelectric_point()), 1)
    except Exception:  # noqa: BLE001
        pI = 0.0
    if density >= 0.30:
        liab.append(f"과도한 양전하(밀도 {density}) — 비특이 결합·독성·빠른 청소 위험")

    # 3) 응집 경향(휴리스틱): 소수성·방향족 비율 + 최장 연쇄
    hp_frac = sum(1 for c in seq if c in _HYDROPHOBIC) / len(seq)
    run = _longest_hydrophobic_run(seq)
    agg = round(min(1.0, 0.7 * hp_frac + 0.1 * run), 2)
    if agg >= 0.45:
        liab.append(f"응집 경향 높음(소수성 {hp_frac:.0%}, 최장연쇄 {run})")

    r = DevelopabilityResult(
        liabilities=liab, n_liabilities=len(liab), free_cys=free_cys,
        net_charge=net, charge_density=density, pI=pI, agg_score=agg,
    )
    # 위험 등급: liability 수 + 강한 신호 가중
    strong = sum(1 for x in liab if "과도한 양전하" in x or "응집 경향" in x or "자유 시스테인" in x)
    if len(liab) >= 4 or strong >= 2:
        r.risk_level = "높음"
    elif len(liab) >= 2 or strong >= 1:
        r.risk_level = "보통"
    else:
        r.risk_level = "낮음"
    r.verdict = (f"개발성 위험 {r.risk_level} — liability {len(liab)}개, 순전하 {net}"
                 f"(밀도 {density}), 응집 {agg}")
    return r
