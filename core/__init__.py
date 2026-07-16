"""BBB-Optimize AI Agent — core package.

레이어 구성:
    config      : 상수 + 환경변수(API 키/엔드포인트) 로딩, mock/real 모드 판별
    schemas     : 계층 간 주고받는 데이터 구조(dataclass)
    predictors  : BBB 투과율/독성 예측기 (인터페이스 + Mock + 실제 API 스텁)
    agent       : 최적화 오케스트레이션 (변이 전략 + 예측 + 필터링)

UI(app.py)는 오직 `agent`와 `config`만 알면 되고,
예측기가 mock인지 실제 API인지는 알 필요가 없다(의존성 주입).
"""

from .agent import OptimizerAgent, build_agent
from .config import Settings, get_settings

__all__ = ["OptimizerAgent", "build_agent", "Settings", "get_settings"]
