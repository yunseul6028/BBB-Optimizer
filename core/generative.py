"""FBGAN 생성 최적화 (잠재공간 진화 피드백).

사전학습 생성기로 novel 셔틀 후보를 생성하고, deepB3P(BBB)·ToxinPred3(독성)를
적합도로 삼아 latent z 를 진화시킨다. 실제 루프는 vendor/deepB3P/_run_fbgan.py 가
deepB3P venv 에서 돌고, 여기서는 subprocess로 호출·파싱만 한다.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from .config import (
    DEEPB3P_PYTHON,
    DEEPB3P_REPO,
    FBGAN_RUNNER,
    TOXINPRED3_PYTHON,
    TOXINPRED3_REPO,
    Settings,
)


@dataclass
class FBGANResult:
    history: list[dict] = field(default_factory=list)   # [{round, mean_bbb, best_bbb, n_safe}]
    best: list[dict] = field(default_factory=list)       # [{shuttle, sequence, bbb, tox, len}]
    cargo: str = ""
    linker: str = ""


class FBGANOptimizer:
    def __init__(self, timeout: float = 900.0):
        self.timeout = timeout

    def run(self, cargo: str, linker: str, rounds: int = 4, pop: int = 32,
            elite: int = 6, tox_threshold: float = 0.38) -> FBGANResult:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "fbgan.json"
            cmd = [
                str(DEEPB3P_PYTHON), str(FBGAN_RUNNER),
                cargo.upper(), linker.upper(), str(rounds), str(pop), str(elite),
                str(tox_threshold),
                str(TOXINPRED3_PYTHON), str(TOXINPRED3_REPO), "toxinpred3.py",
                str(out),
            ]
            proc = subprocess.run(cmd, cwd=str(DEEPB3P_REPO), capture_output=True,
                                  text=True, timeout=self.timeout)
            if proc.returncode != 0 or not out.exists():
                raise RuntimeError(
                    f"FBGAN 생성 최적화 실패 (rc={proc.returncode}).\n"
                    f"stderr(마지막 600자): {proc.stderr[-600:]}"
                )
            data = json.loads(out.read_text())
        return FBGANResult(history=data.get("history", []), best=data.get("best", []),
                           cargo=data.get("cargo", cargo), linker=data.get("linker", linker))


def get_fbgan(settings: Settings) -> FBGANOptimizer | None:
    return FBGANOptimizer() if settings.use_fbgan_local else None
