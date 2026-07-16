"""계층 간 데이터 구조.

predictors → agent → app(UI) 사이에서 주고받는 값 객체들.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class StepLevel(str, Enum):
    INFO = "info"
    WARN = "warn"
    SUCCESS = "success"
    ERROR = "error"


class Verdict(str, Enum):
    ACCEPTED = "accepted"              # 최종 채택
    SUBOPTIMAL = "suboptimal"          # 생존했으나 후순위
    REJECTED_TOXIC = "rejected_toxic"  # 독성 임계값 초과 탈락
    REFERENCE = "reference"            # 화물 단독(셔틀 없음) 기준선


@dataclass
class Construct:
    """화물 + 링커 + 셔틀로 조립된 융합체 후보 하나."""
    cargo: str                     # 화물(payload) 펩타이드
    linker_name: str               # 링커 라이브러리 이름 (예: "YGGGGS"); 기준선은 ""
    linker: str                    # 링커 서열; 셔틀 없는 기준선은 ""
    shuttle_name: str              # 셔틀 이름 (예: "Angiopep-2"); 없으면 "(없음)"
    shuttle_seq: str               # 셔틀 서열
    sequence: str                  # 조립된 전체 서열 (cargo+linker+shuttle) — 독성·표시용
    bbb_sequence: str = ""         # BBB 계산용 서열 (초과 시 연결부위 윈도우)
    is_standard_linker: bool = False   # 표준 GGGGGS 링커 여부
    truncated: bool = False        # 50aa 초과 → BBB는 연결부위 윈도우로 계산

    @property
    def label(self) -> str:
        if not self.linker and self.shuttle_name == "(없음)":
            return "화물 단독"
        return f"{self.shuttle_name} · {self.linker_name}"


@dataclass
class PredictionResult:
    """예측기가 반환하는 스코어. mock이든 실제 모델이든 이 형태로 통일."""
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


@dataclass
class EvaluatedCandidate:
    """융합체 + 예측결과 + 판정."""
    construct: Construct
    prediction: PredictionResult
    verdict: Verdict = Verdict.SUBOPTIMAL


@dataclass
class AgentStep:
    """에이전트 추론 로그 한 줄."""
    stage: str
    message: str
    level: StepLevel = StepLevel.INFO


@dataclass
class OptimizationResult:
    """최종 산출물."""
    cargo: str
    cargo_only: EvaluatedCandidate | None            # 화물 단독 기준선
    winner: EvaluatedCandidate | None                # 최적 (셔틀,링커) 조합
    candidates: list[EvaluatedCandidate] = field(default_factory=list)
    n_linkers: int = 0
    n_shuttles: int = 0
    truncated_any: bool = False
