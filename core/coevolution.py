"""링커·셔틀 co-evolution 오케스트레이터 (멀티모듈 동시 진화 — 서열 directed evolution).

링커(펩타이드 융합 링커)와 셔틀을 함께 진화시킨다. 셔틀은 **de-novo 생성이 아니라 검증
라이브러리 리간드에서 시드**해 보존적으로 진화한다(SHUTTLES 서열을 러너에 시드로 전달).
실제 루프는 vendor/deepB3P/_run_coevo.py 가 deepB3P venv 에서 돌고, 여기서는 subprocess
호출·파싱만 한다.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from .config import (
    COEVO_RUNNER,
    DEEPB3P_PYTHON,
    DEEPB3P_REPO,
    SHUTTLES,
    TOXINPRED3_PYTHON,
    TOXINPRED3_REPO,
    Settings,
)

# 셔틀 진화 시드 — 검증된 라이브러리 리간드 서열(de-novo 생성 금지, 여기서 시드해 진화).
_SHUTTLE_SEEDS = ",".join(v["seq"] for v in SHUTTLES.values() if v.get("seq"))


@dataclass
class CoevoResult:
    history: list[dict] = field(default_factory=list)   # [{round, mean_bbb, best_bbb, n_safe, n_uniq_linker, best_linker}]
    best: list[dict] = field(default_factory=list)       # [{shuttle, linker, sequence, bbb, tox, len, linker_len}]
    cargo: str = ""


class CoevolutionOptimizer:
    def __init__(self, timeout: float = 1200.0):
        self.timeout = timeout

    def run(self, cargo: str, rounds: int = 4, pop: int = 32,
            elite: int = 6, tox_threshold: float = 0.38) -> CoevoResult:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "coevo.json"
            cmd = [
                str(DEEPB3P_PYTHON), str(COEVO_RUNNER),
                cargo.upper(), str(rounds), str(pop), str(elite), str(tox_threshold),
                str(TOXINPRED3_PYTHON), str(TOXINPRED3_REPO), "toxinpred3.py",
                str(out), _SHUTTLE_SEEDS,
            ]
            proc = subprocess.run(cmd, cwd=str(DEEPB3P_REPO), capture_output=True,
                                  text=True, timeout=self.timeout)
            if proc.returncode != 0 or not out.exists():
                raise RuntimeError(
                    f"링커·셔틀 co-evolution 실패 (rc={proc.returncode}).\n"
                    f"stderr(마지막 600자): {proc.stderr[-600:]}"
                )
            data = json.loads(out.read_text())
        return CoevoResult(history=data.get("history", []), best=data.get("best", []),
                           cargo=data.get("cargo", cargo))


def get_coevolution(settings: Settings) -> CoevolutionOptimizer | None:
    return CoevolutionOptimizer() if settings.use_coevo_local else None
