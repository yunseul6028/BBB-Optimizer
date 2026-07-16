"""BBB 투과율 예측기.

인터페이스(BBBPredictor) 아래에 세 구현:
  - MockBBBPredictor       : 간이 휴리스틱 (오프라인 폴백)
  - DeepB3PLocalPredictor  : 로컬 deepB3P 딥러닝 추론 (실측 BBB) ← 현재 기본
  - DeepB3ApiPredictor     : 원격 예측 서버 (자체 호스팅 시)

⚠️ 이 레이어는 **BBB 투과(효능)만** 예측한다. 독성(부작용)은 아직 연결된 도구가
   없어 에이전트가 placeholder(LINKER_TOX_TABLE)로 채운다. [[toxicity-engine-todo]]
   따라서 여기서 반환하는 toxicity_risk는 중립 기본값이며 에이전트가 덮어쓴다.
"""

from __future__ import annotations

import csv
import subprocess
import tempfile
from abc import ABC, abstractmethod
from pathlib import Path

from .config import DEEPB3P_PYTHON, DEEPB3P_REPO, DEEPB3P_RUNNER, Settings, junction_window
from .schemas import PredictionResult

_DEFAULT_TOX = 0.20  # 에이전트가 링커 잔기 기준으로 덮어씀


class BBBPredictor(ABC):
    source: str = "base"

    @abstractmethod
    def predict(self, sequence: str) -> PredictionResult:
        raise NotImplementedError

    def predict_many(self, sequences: list[str]) -> list[PredictionResult]:
        """여러 서열 일괄 예측. 기본은 단건 반복 (배치 이득 구현은 오버라이드)."""
        return [self.predict(s) for s in sequences]


# ---------------------------------------------------------------------------
# Mock — 오프라인 폴백. 소수성/방향족 잔기 비율 기반 간이 휴리스틱.
# ---------------------------------------------------------------------------
class MockBBBPredictor(BBBPredictor):
    source = "mock"
    _HYDROPHOBIC = set("AVLIMFWYC")
    _AROMATIC = set("FWY")

    def predict(self, sequence: str) -> PredictionResult:
        seq = sequence.upper()
        if not seq:
            return PredictionResult(0.0, _DEFAULT_TOX, self.source, "빈 서열")
        hyd = sum(c in self._HYDROPHOBIC for c in seq) / len(seq)
        aro = sum(c in self._AROMATIC for c in seq) / len(seq)
        bbb = max(0.0, min(1.0, 0.15 + 0.7 * hyd + 0.3 * aro))
        return PredictionResult(
            bbb_permeability=round(bbb, 4),
            toxicity_risk=_DEFAULT_TOX,
            source=self.source,
            note="휴리스틱 근사(소수성/방향족 비율)",
            raw={"hydrophobic_frac": round(hyd, 3), "aromatic_frac": round(aro, 3)},
        )


# ---------------------------------------------------------------------------
# 로컬 deepB3P 추론 — 실측 BBB.
# ---------------------------------------------------------------------------
class DeepB3PLocalPredictor(BBBPredictor):
    """vendor/deepB3P 의 5-fold 앙상블을 subprocess로 호출해 실측 BBB 확률을 얻는다."""

    source = "deepb3p-local"

    def __init__(self, python: Path = DEEPB3P_PYTHON, repo: Path = DEEPB3P_REPO,
                 runner: Path = DEEPB3P_RUNNER, timeout: float = 300.0):
        self.python = str(python)
        self.repo = str(repo)
        self.runner = str(runner)
        self.timeout = timeout

    def predict(self, sequence: str) -> PredictionResult:
        return self.predict_many([sequence])[0]

    def predict_many(self, sequences: list[str]) -> list[PredictionResult]:
        if not sequences:
            return []
        probs = self._infer(sequences)
        return [
            PredictionResult(
                bbb_permeability=prob,
                toxicity_risk=_DEFAULT_TOX,
                source=self.source,
                note="deepB3P 5-fold 앙상블 실측",
                raw={"prob": prob},
            )
            for prob in probs
        ]

    def _infer(self, sequences: list[str]) -> list[float]:
        """FASTA 작성 → subprocess 추론 → 확률을 입력 순서대로 반환."""
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "in.fasta"
            out = Path(tmp) / "out.csv"
            # 50aa 초과는 연결부위(링커+셔틀=C말단) 윈도우로 계산 (앞자르기 방지)
            with open(fasta, "w") as f:
                for i, seq in enumerate(sequences):
                    f.write(f">{i}\n{junction_window(seq)}\n")

            proc = subprocess.run(
                [self.python, self.runner, str(fasta), str(out)],
                cwd=self.repo, capture_output=True, text=True, timeout=self.timeout,
            )
            if proc.returncode != 0 or not out.exists():
                raise RuntimeError(
                    f"deepB3P 추론 실패 (rc={proc.returncode}).\n"
                    f"stderr(마지막 500자): {proc.stderr[-500:]}"
                )
            with open(out) as f:
                rows = list(csv.DictReader(f))
        return [float(r["prob"]) for r in rows]  # 입력 순서 = 출력 순서


# ---------------------------------------------------------------------------
# 원격 API 구현 — 자체 호스팅 시.
# ---------------------------------------------------------------------------
class DeepB3ApiPredictor(BBBPredictor):
    source = "deep-b3-api"

    def __init__(self, api_url: str, api_key: str, timeout: float = 30.0):
        self.api_url = api_url.rstrip("/")
        self.api_key = api_key
        self.timeout = timeout

    def predict(self, sequence: str) -> PredictionResult:
        # TODO(원격 연동): requests.post 골격 채우기.
        raise NotImplementedError("DeepB3ApiPredictor 미구현 — 자체 서버가 있을 때만 사용.")


# ---------------------------------------------------------------------------
# 팩토리 — 우선순위: 로컬 deepB3P > 원격 API > Mock
# ---------------------------------------------------------------------------
def get_predictor(settings: Settings) -> BBBPredictor:
    if settings.use_deepb3p_local:
        return DeepB3PLocalPredictor()
    if settings.use_real_predictor:
        return DeepB3ApiPredictor(settings.deepb3_api_url, settings.deepb3_api_key)  # type: ignore[arg-type]
    return MockBBBPredictor()
