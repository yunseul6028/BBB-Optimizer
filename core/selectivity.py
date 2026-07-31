"""선택성 / off-target 위험 — 비특이 조직 흡수 가능성.

양전하·친유성으로 세포막을 직접 뚫는 방식(CPP)은 **여러 조직에 비특이적으로 흡수**되어
off-target 효과·독성을 낸다. 반면 **수용체 매개 수송(RMT, 예: Angiopep-2 → LRP1)** 은
특정 수용체가 있는 조직(BBB)만 노려 **선택적**이다.

이 모듈은 셔틀의 (1) 양전하 밀도, (2) 친유성(GRAVY), (3) 수송 메커니즘(RMT vs CPP)을
종합해 **off-target 위험 / 선택성**을 근사한다. 핵심 논리:
  - 양전하↑ 또는 친유성↑ → 비특이 흡수 위험↑ (둘 중 하나만 높아도 위험 = 확률적 OR)
  - RMT 수용체 표적성↑ → off-target 위험↓ (특정 수용체만 노리므로)

⚠️ 진짜 off-target은 생체분포(biodistribution) 실험으로 확정한다. 서열 기반 위험 프록시다.
"""

from __future__ import annotations

from dataclasses import dataclass, field

_STD = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class SelectivityResult:
    off_target_risk: float = 0.0   # 0~1 (높을수록 비특이 흡수 위험)
    selectivity: float = 1.0       # 1 - off_target_risk
    charge_density: float = 0.0    # pH7.4 순전하 / 길이
    gravy: float = 0.0             # 친유성(양수=소수성)
    mechanism: str = ""            # RMT형(선택적) / CPP형(비특이) / 불명확
    rmt_similarity: float = 0.0
    risk_level: str = "낮음"        # 낮음 | 보통 | 높음
    drivers: list = field(default_factory=list)
    verdict: str = ""
    error: str = ""


def _clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def assess_selectivity(shuttle: str) -> SelectivityResult:
    """셔틀의 off-target(비특이 흡수) 위험과 선택성을 평가한다."""
    seq = "".join(c for c in (shuttle or "").upper() if c in _STD)
    if len(seq) < 4:
        return SelectivityResult(error="셔틀이 너무 짧아 선택성 계산 불가")

    from .developability import _net_charge_ph74
    net = _net_charge_ph74(seq)
    dens = round(net / len(seq), 2)
    cationicity = _clamp(dens / 0.5) if dens > 0 else 0.0   # 밀도 +0.5 = 최대 위험

    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        gravy = round(float(ProteinAnalysis(seq).gravy()), 2)
    except Exception:  # noqa: BLE001
        gravy = 0.0
    lipophilicity = _clamp(gravy / 0.8) if gravy > 0 else 0.0  # GRAVY +0.8 = 최대

    from .binding import REFERENCE_SHUTTLES, shuttle_similarity
    bind = shuttle_similarity(seq)
    rmt_refs = {n for n, r in REFERENCE_SHUTTLES.items()
                if "RMT" in r["mech"] or "수용체" in r["mech"]}
    per = {d["name"]: d["score"] for d in (bind.per_ref or [])}
    rmt_sim = max((per.get(n, 0.0) for n in rmt_refs), default=0.0)
    best_is_rmt = bind.best_ref in rmt_refs and bind.score >= 0.4

    # 비특이 raw: 양전하 또는 친유성 중 하나만 높아도 위험 (확률적 OR)
    raw = 1 - (1 - cationicity) * (1 - lipophilicity)
    # RMT 수용체 표적성은 off-target 위험을 낮춤
    off = round(_clamp(raw * (1 - 0.6 * rmt_sim)), 2)
    sel = round(1 - off, 2)

    drivers = []
    if cationicity >= 0.4:
        drivers.append(f"양전하 밀도 {dens} — 비특이 정전 흡수 위험")
    if lipophilicity >= 0.35:
        drivers.append(f"친유성 GRAVY {gravy} — 비특이 막삽입 위험")
    cpp_scores = {n: per.get(n, 0.0) for n in REFERENCE_SHUTTLES if n not in rmt_refs}
    best_cpp = max(cpp_scores, key=cpp_scores.get, default="")
    if best_is_rmt:
        mechanism = f"RMT형({bind.best_ref} 유사·수용체 선택적)"
        drivers.append("RMT 수용체 표적 → 선택성↑")
    elif best_cpp and cpp_scores[best_cpp] >= 0.4:
        mechanism = f"CPP형({best_cpp} 유사·비특이 막투과)"
    else:
        mechanism = "비특이 막투과형(전하/친유성 주도, 뚜렷한 수용체 유사성 낮음)"

    level = "높음" if off >= 0.6 else "보통" if off >= 0.33 else "낮음"
    r = SelectivityResult(
        off_target_risk=off, selectivity=sel, charge_density=dens, gravy=gravy,
        mechanism=mechanism, rmt_similarity=round(rmt_sim, 2), risk_level=level,
        drivers=drivers)
    r.verdict = f"off-target 위험 {level}(선택성 {sel}) — {mechanism}"
    return r
