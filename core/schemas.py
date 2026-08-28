"""계층 간 데이터 구조.

예측기 → 에이전트 → app(UI) 사이에서 주고받는 값 객체들.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any
from enum import Enum


class Verdict(str, Enum):
    ACCEPTED = "accepted"              # 최종 채택
    SUBOPTIMAL = "suboptimal"          # 생존했으나 후순위
    REJECTED_TOXIC = "rejected_toxic"  # 독성 임계값 초과 탈락
    REFERENCE = "reference"            # 화물 단독(셔틀 없음) 기준선


@dataclass
class PredictionResult:
    """BBB 예측기가 반환하는 스코어. 로컬/원격/폴백 공통 형태."""
    bbb_permeability: float        # 0.0 ~ 1.0
    toxicity_risk: float           # 0.0 ~ 1.0 (현재 placeholder)
    source: str = "mock"           # "mock" | "deepb3p-local" | "deep-b3-api"
    note: str = ""
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass
class ToxicityResult:
    """독성 예측기가 반환하는 결과."""
    risk: float                    # 0.0 ~ 1.0 (높을수록 독성)
    is_toxic: bool = False
    source: str = "placeholder"    # "toxinpred3-local" | "placeholder"
    note: str = ""
