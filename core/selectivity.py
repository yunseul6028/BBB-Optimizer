"""전하 기반 비특이 흡수 위험(off-target) — 🔴 **선택성 모듈이 아니다.**

이 파일의 옛 제목은 "선택성 / off-target 위험"이었고, 문헌 논리(CPP 는 비특이, RMT 는
수용체가 있는 조직만 노려 선택적)를 그대로 모듈 이름으로 삼았다. 그 이름이 **측정하지
않은 것을 주장한다**는 것이 세 번의 실측으로 드러났다:

  · §50-5 — 친유성항은 `GRAVY > 0` 에서만 켜지는데 라이브러리 **11/12 가 GRAVY ≤ 0** 이다.
    즉 이 축은 사실상 **양전하 밀도 하나**이고, 전하만으로 계산한 값과 정확히 같다.
  · §50-3 — 위험을 실제로 정하는 **표적 수용체의 조직 분포를 이 축은 본 적이 없다.**
    외부 발현 대조에서 축이 **낙관한 셔틀 4종 · 반대 방향 0종**이었다.
  · §26 — 대응 널 심사에서 **악화 검출만 신뢰할 수 있고 개선 신호는 신뢰할 수 없다.**
    그래서 탐색에서는 단측 상한이 아니라 **밴드 제약**(`OFF_BAND`)으로만 쓴다.

계산: 셔틀의 (1) 양전하 밀도, (2) 친유성(GRAVY) 을 확률적 OR 로 합치고,
(3) RMT 참조 서열과의 유사도로 그 값을 **할인**한다(`off = raw * (1 − 0.6 * rmt_sim)`).
할인은 수용체 선택성의 측정이 아니라 **문헌 셔틀을 덜 벌주기 위한 보정**이다.

반환 필드 `selectivity`(= 1 − off)는 하위 호환 때문에 이름이 남아 있을 뿐이다 —
**사용자·LLM 에게 "선택성"이라는 말로 보고하지 않는다.**

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
    """셔틀의 **전하 기반 비특이 흡수 위험**(off-target)을 평가한다.

    🔴 이름 주의 — 이것은 **수용체 선택성이 아니다.** ⑴ 친유성항은 `GRAVY > 0` 에서만
    켜지는데 라이브러리 **11/12 가 GRAVY ≤ 0** 이라, off 가 **전하항만으로 계산한 값과
    정확히 같다**(§50-5). ⑵ 위험을 실제로 정하는 **표적 수용체의 조직 분포를 이 축은 본
    적이 없다** — 외부 발현 대조에서 축이 **낙관한 셔틀 4종 · 반대 방향 0종**이었다(§50-3).
    ⑶ 대응 널 심사 판정은 **악화 검출만 신뢰할 수 있고 개선 신호는 신뢰할 수 없다**(§26).
    그래서 반환 필드 `selectivity`(= 1 − off)를 **"선택성"이라는 말로 보고하지 않는다.**
    """
    seq = "".join(c for c in (shuttle or "").upper() if c in _STD)
    # 정당하게 계산 불가한 입력(의존성 문제가 아님): 4aa 미만 셔틀은 전하밀도·GRAVY·
    # 레퍼런스 정렬이 모두 무의미하다 → 예외가 아니라 `error` 를 채운 정의된 결과.
    if len(seq) < 4:
        return SelectivityResult(error="셔틀이 너무 짧아 off-target 축 계산 불가")

    from .developability import _net_charge_ph74, require_protein_analysis
    net = _net_charge_ph74(seq)
    dens = round(net / len(seq), 2)
    cationicity = _clamp(dens / 0.5) if dens > 0 else 0.0   # 밀도 +0.5 = 최대 위험

    # ⚠️ 여기가 실측된 조용한 퇴화 지점이다: Biopython 이 없으면 gravy 가 0.0 으로
    #    떨어져 친유성 항이 사라지고 T7 의 off-target 위험이 0.140 → 0.340 으로
    #    **경고 없이** 움직였다. 의존성 부재는 배포 오류 → 즉시 중단.
    ProteinAnalysis = require_protein_analysis("off-target 축의 친유성항(GRAVY)")
    # 감싸지 않는 이유: seq 는 이미 표준 20종으로 걸러진 len>=4 서열이라 gravy() 는
    # 실패할 수 없다. 그래도 실패하면 진짜 버그이므로 예외를 그대로 올린다.
    gravy = round(float(ProteinAnalysis(seq).gravy()), 2)
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
        # 🔴 옛 라벨은 "RMT형(… 유사·수용체 **선택적**)" / "RMT 수용체 표적 → **선택성↑**"
        #    이었다. 두 문구 다 **측정하지 않은 것을 주장한다** — 이 축은 표적 수용체의
        #    조직 분포를 본 적이 없고(§50-3), 하는 일은 RMT 참조 서열과의 **유사도로
        #    전하항을 할인**하는 것뿐이다(`off = raw * (1 − 0.6 * rmt_sim)`).
        #    이 문자열은 `sel_mech` 로 **심사 LLM 에 그대로 들어간다**(§67).
        mechanism = f"RMT형({bind.best_ref} 서열 유사)"
        drivers.append("RMT 참조 서열과 유사 → 전하 기반 위험항을 할인했다 "
                       "(🔴 수용체 선택성을 측정한 것이 아니다)")
    elif best_cpp and cpp_scores[best_cpp] >= 0.4:
        mechanism = f"CPP형({best_cpp} 유사·비특이 막투과)"
    else:
        mechanism = "비특이 막투과형(전하/친유성 주도, 뚜렷한 수용체 유사성 낮음)"

    level = "높음" if off >= 0.6 else "보통" if off >= 0.33 else "낮음"
    r = SelectivityResult(
        off_target_risk=off, selectivity=sel, charge_density=dens, gravy=gravy,
        mechanism=mechanism, rmt_similarity=round(rmt_sim, 2), risk_level=level,
        drivers=drivers)
    # 🔴 옛 문구는 "(**선택성** {sel})" 이었다. 잰 것은 전하(+친유성)이고 수용체
    #    선택성이 아니다(§50). 지금은 LLM 입력에 쓰이지 않지만 되살아날 때를 대비한다.
    r.verdict = f"비특이 흡수 위험 {level}(축값 1−off = {sel}) — {mechanism}"
    return r
