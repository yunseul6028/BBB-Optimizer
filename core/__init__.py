"""BBB-Optimize AI Agent — core package.

레이어 구성:
    config              : 상수 + 환경변수(API 키/엔드포인트) 로딩, 로컬/원격 엔진 판별
    schemas             : 계층 간 주고받는 데이터 구조(dataclass)
    predictors/toxicity : BBB 투과율·독성 예측기 (로컬 실측 > 원격 API > 폴백)
    optimizer_agent*    : 자율 설계 에이전트(Designer–Critic tool-use 루프)
    delivery/binding/…  : 8축 평가 엔진(선택성·delivery 분해·개발성·구조·용해도 등)

UI(app.py)는 오직 에이전트와 `config`만 알면 되고, 예측기가 로컬인지 원격인지는
알 필요가 없다(의존성 주입).
"""

from .config import Settings, get_settings

__all__ = ["Settings", "get_settings"]
