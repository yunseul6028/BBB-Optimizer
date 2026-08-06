"""BBB 전달 축 분해 (RMT-aware) — deepB3P 한 숫자의 개념 혼동을 제거한다.

배경(왜 필요한가):
    deepB3P는 **짧은 펩타이드의 BBB 투과 확률** 분류기다. 이 점수를 융합체 하나의
    "투과율"처럼 제시하면 → *펩타이드 투과 점수 ↔ 항체 RMT 전달* 개념 혼동이 생긴다.
    실제로 큰 항체 융합체가 뇌로 들어가는 병목은 확산 투과가 아니라 **수용체 매개
    트랜스시토시스(RMT)** — 즉 수용체 결합·avidity가 결정한다. 그래서 deepB3P의 유효성은
    셔틀의 **수송 메커니즘에 따라 다르다**:
        · CPP(막 직접투과) 셔틀 → deepB3P가 학습한 현상에 가까움 → 점수가 비교적 유효.
        · RMT(수용체 매개) 셔틀 → 병목이 수용체 결합이므로 deepB3P는 **약한 프록시**.
                                    구조 노출·친화도(avidity)로 보완해야 한다.

이 모듈은 deepB3P 점수를 근거 있는 하위 축으로 **분해**한다(오라클을 새로 지어내지 않는다):
    ① 셔틀 내재 투과   = deepB3P(셔틀 단독)     — 셔틀 자체의 BBB 투과 성향(펩타이드→타당)
    ② 융합 보존        = 융합 점수 / 셔틀 내재   — 화물·링커를 붙여도 셔틀 신호가 보존되나
    ③ 메커니즘 타당도  = RMT/CPP 분류 → deepB3P 신뢰 수준을 정직하게 라벨
    ④ avidity          = 친화도 sweet-spot(KD 있을 때만; 펩타이드는 N/A → 항체 Tier-3)

⚠️ 전부 서열/휴리스틱 **프록시**다. 실측(in vitro BBB transwell·뇌흡수) 대체가 아니다.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .antibody import transcytosis_sweetspot

# 메커니즘 분류용 키워드 (target/mechanism 문자열에서 탐지).
_CPP_HINTS = ("CPP", "세포투과", "막 직접투과", "막직접투과", "막 투과", "직접투과")
_RMT_HINTS = ("RMT", "수용체매개", "수용체 매개", "LRP", "LDLR", "TFR", "트랜스페린",
              "렙틴", "LEPR", "IGF", "INSR", "CD98", "수용체")

# 셔틀 내재 deepB3P 신호가 이 값 미만이면 '보존도' 비율이 무의미(0 나눗셈·노이즈) → 판단 보류.
_INTRINSIC_FLOOR = 0.05


@dataclass
class DeliveryAxis:
    """deepB3P 점수를 분해한 전달 축. 모든 값은 0~1 프록시(또는 None=계산 불가)."""
    shuttle_intrinsic: float          # ① deepB3P(셔틀 단독)
    fusion_bbb: float                 # 융합 연결부위 deepB3P (= 기존 'bbb')
    mechanism: str                    # "RMT" | "CPP" | "불확실"
    is_rmt: bool
    preservation: float | None = None # ② 융합 보존도 (None=셔틀 신호 미약해 판단 불가)
    avidity: float | None = None      # ④ 친화도 sweet-spot (KD 없으면 None)
    kd_nM: float | None = None
    deepb3p_validity: str = ""        # ③ 이 메커니즘에서 deepB3P 점수를 얼마나 믿을지
    basis: str = ""                   # 한 줄 종합 해석
    caveats: list = field(default_factory=list)


def classify_mechanism(target: str = "", mechanism: str = "") -> tuple[str, bool]:
    """(target, mechanism) 문자열 → ('RMT'|'CPP'|'불확실', is_rmt)."""
    text = f"{target} {mechanism}".upper()
    # CPP를 먼저 본다(일부 CPP 설명에 '수용체' 무관 단어가 섞일 수 있어 우선순위 부여).
    if any(h.upper() in text for h in _CPP_HINTS):
        return "CPP", False
    if any(h.upper() in text for h in _RMT_HINTS):
        return "RMT", True
    return "불확실", False


def _validity_note(mechanism: str) -> str:
    if mechanism == "CPP":
        return "CPP(막투과) — deepB3P 학습 현상에 부합, 점수 비교적 유효"
    if mechanism == "RMT":
        return ("RMT(수용체 매개) — 병목은 수용체 결합이므로 deepB3P는 약한 프록시. "
                "구조 노출·avidity로 보완 필요")
    return "메커니즘 불확실 — deepB3P 절대값 신뢰 낮음, 상대 비교로만 해석"


def assess_delivery(shuttle_intrinsic: float, fusion_bbb: float,
                    target: str = "", mechanism: str = "",
                    *, kd_nM: float | None = None,
                    valency: int = 1) -> DeliveryAxis:
    """deepB3P 점수(셔틀 단독 + 융합)를 전달 축으로 분해한다.

    shuttle_intrinsic : deepB3P(셔틀 펩타이드 단독) 0~1
    fusion_bbb        : deepB3P(융합 연결부위) 0~1  ← 기존 파이프라인의 'bbb'
    target, mechanism : 셔틀의 표적·수송 메커니즘 라벨(binding.shuttle_similarity 등에서)
    kd_nM, valency    : 항체 모달리티에서 KD/결합가 있을 때만 avidity 계산.
    """
    si = max(0.0, min(1.0, float(shuttle_intrinsic)))
    fb = max(0.0, min(1.0, float(fusion_bbb)))
    mech, is_rmt = classify_mechanism(target, mechanism)
    caveats: list[str] = []

    # ② 융합 보존도 — 셔틀 신호가 충분할 때만 의미. (셔틀 자체가 0점이면 비율 무의미)
    if si >= _INTRINSIC_FLOOR:
        preservation = round(max(0.0, min(1.0, fb / si)), 3)
    else:
        preservation = None
        caveats.append("셔틀 단독 deepB3P 신호 미약 → 보존도·절대값 해석 불가(상대 비교만)")

    # ④ avidity — KD 있을 때만(펩타이드 셔틀은 보통 N/A).
    avidity = transcytosis_sweetspot(kd_nM)
    if avidity is None:
        caveats.append("펩타이드 셔틀·KD 미측정 → avidity(친화도 sweet-spot)는 항체 모달리티(Tier-3)서 계산")
    if valency >= 2:
        caveats.append("bivalent — TfR류 RMT는 monovalent가 수송에 유리(가교→분해경로 위험)")

    validity = _validity_note(mech)

    # 한 줄 종합 해석 — 오라클 숫자 대신, 근거로 프레이밍.
    if mech == "RMT":
        basis = (f"RMT 셔틀 → 실제 전달은 수용체 결합/avidity가 결정. deepB3P {fb:.2f}는 참고치. "
                 + (f"보존 {preservation:.2f}" if preservation is not None else "보존도 판단 불가")
                 + (f" · avidity {avidity:.2f}" if avidity is not None else " · avidity 미산정"))
    elif mech == "CPP":
        basis = (f"CPP 셔틀 → deepB3P {fb:.2f} 비교적 유효(막투과). "
                 + (f"셔틀 내재 {si:.2f}×보존 {preservation:.2f}" if preservation is not None
                    else "셔틀 내재 신호 미약") + " · 단 CPP는 비특이(off-target) 주의")
    else:
        basis = f"메커니즘 불확실 → deepB3P {fb:.2f}는 상대 비교용으로만."

    return DeliveryAxis(
        shuttle_intrinsic=round(si, 4), fusion_bbb=round(fb, 4),
        mechanism=mech, is_rmt=is_rmt, preservation=preservation,
        avidity=avidity, kd_nM=kd_nM, deepb3p_validity=validity,
        basis=basis, caveats=caveats,
    )
