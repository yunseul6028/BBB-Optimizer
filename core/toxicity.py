"""독성(부작용) 예측기.

BBBPredictor와 대칭되는 독성 축. 인터페이스 아래에:
  - ToxinPred3LocalPredictor : 로컬 ToxinPred3 (ML Score) — 실측 ← 기본
  - PlaceholderToxicityPredictor : 링커 잔기 기준 상수 (폴백)

독성은 **전체 융합체 서열**에 대해 평가한다(링커 잔기만이 아니라 실제 분자 전체).
ToxinPred3 Model 1(-m 1)은 AAC+DPC 조성 기반 Extra-Tree로 외부도구(BLAST/MERCI) 불필요.
"""

from __future__ import annotations

import csv
import fcntl
import subprocess
import tempfile
from abc import ABC, abstractmethod
from pathlib import Path

from .config import (
    TOXINPRED3_MODEL_ARG,
    TOXINPRED3_PYTHON,
    TOXINPRED3_REPO,
    TOXINPRED3_SCRIPT,
    Settings,
)
from .schemas import ToxicityResult


class ToxicityPredictor(ABC):
    source: str = "base"

    def __init__(self, threshold: float):
        self.threshold = threshold

    @abstractmethod
    def predict_many(self, sequences: list[str]) -> list[ToxicityResult]:
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Placeholder — 링커 첫 잔기 기준 상수 (도구 미설치 시 폴백).
# ---------------------------------------------------------------------------
class PlaceholderToxicityPredictor(ToxicityPredictor):
    source = "placeholder"

    def predict_one(self, sequence: str) -> ToxicityResult:
        # 링커 첫 잔기 = cargo 다음 위치를 모르므로, 서열 내 잔기 조성으로 근사 불가.
        # 폴백은 보수적으로 중립값 사용.
        risk = 0.20
        return ToxicityResult(risk, risk > self.threshold, self.source,
                              "placeholder(중립값, 도구 미연결)")

    def predict_many(self, sequences: list[str]) -> list[ToxicityResult]:
        return [self.predict_one(s) for s in sequences]


# ---------------------------------------------------------------------------
# 로컬 ToxinPred3 — 실측.
# ---------------------------------------------------------------------------
class ToxinPred3LocalPredictor(ToxicityPredictor):
    """vendor/toxinpred3 의 Model 1(ML)을 subprocess로 호출해 독성 확률(ML Score)을 얻는다."""

    source = "toxinpred3-local"

    def __init__(self, threshold: float, python: Path = TOXINPRED3_PYTHON,
                 repo: Path = TOXINPRED3_REPO, script: Path = TOXINPRED3_SCRIPT,
                 model_arg: int = TOXINPRED3_MODEL_ARG, timeout: float = 300.0):
        super().__init__(threshold)
        self.python = str(python)
        self.repo = str(repo)
        self.script = str(script)
        self.model_arg = model_arg
        self.timeout = timeout

    def predict_many(self, sequences: list[str]) -> list[ToxicityResult]:
        if not sequences:
            return []
        scores = self._infer(sequences)
        return [
            ToxicityResult(
                risk=s,
                is_toxic=s > self.threshold,
                source=self.source,
                note=f"ToxinPred3 ML Score (Model {self.model_arg})",
            )
            for s in scores
        ]

    def _infer(self, sequences: list[str]) -> list[float]:
        # ToxinPred3는 서열 1개면 np.loadtxt가 1D를 반환해 crash → 더미로 패딩(결과는 버림)
        n = len(sequences)
        write_seqs = sequences + ["AAAAAAAA"] if n == 1 else sequences
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "in.fa"
            out = Path(tmp) / "out.csv"
            with open(fasta, "w") as f:
                for i, seq in enumerate(write_seqs):
                    f.write(f">{i}\n{seq.upper()}\n")

            # 🔴 **벤더 작업 디렉터리 배타 잠금** (2026-11-12)
            # ToxinPred3 는 중간 산물(`seq.aac`·`seq.dpc`)을 **cwd 에 고정된 이름으로**
            # 쓴다. cwd 가 공유 벤더 디렉터리이므로 **동시 호출이 서로의 중간 파일을
            # 덮어쓴다.** 실측 피해: 같은 입력에 `reachable_front(refresh=True)` 를 세 번
            # 돌렸을 때 후보군이 **89 / 91 / 99** 로 갈렸다(전선은 같았다). 즉 크래시로
            # 끝나지 않고 **점수가 조용히 오염된다** — 측정 인프라가 측정을 바꾼 사례다.
            # 처방: 파일 잠금으로 **직렬화**한다. 동시 실행이 아니면 동작은 완전히 동일하고,
            # 동시 실행이면 느려지는 대신 정확해진다(느린 쪽이 맞는 쪽이다).
            lock_path = Path(self.repo) / ".toxinpred3.lock"
            with open(lock_path, "w") as lock:
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
                try:
                    proc = subprocess.run(
                        [self.python, self.script, "-i", str(fasta), "-o", str(out),
                         "-m", str(self.model_arg), "-d", "2"],
                        cwd=self.repo, capture_output=True, text=True,
                        timeout=self.timeout,
                    )
                finally:
                    fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
            if proc.returncode != 0 or not out.exists():
                raise RuntimeError(
                    f"ToxinPred3 추론 실패 (rc={proc.returncode}).\n"
                    f"stderr(마지막 500자): {proc.stderr[-500:]}"
                )
            with open(out) as f:
                rows = {r["ID"]: r for r in csv.DictReader(f)}
        # 입력 순서(ID=인덱스) 복원 — ToxinPred3는 -d 2라도 순서를 보장하지 않을 수 있음
        return [float(rows[str(i)]["ML Score"]) for i in range(len(sequences))]


# ---------------------------------------------------------------------------
# 팩토리
# ---------------------------------------------------------------------------
def get_toxicity_predictor(settings: Settings) -> ToxicityPredictor:
    if settings.use_toxinpred3_local:
        return ToxinPred3LocalPredictor(threshold=settings.toxicity_threshold)
    return PlaceholderToxicityPredictor(threshold=settings.toxicity_threshold)
