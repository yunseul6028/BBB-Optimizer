"""모듈별 최적화 → 조립 재순위 오케스트레이터 (decoupled optimization + assembly).

co-evolution(동시 진화)의 대안. 셔틀은 BBB·전하 공간에서, 링커는 개발성 공간에서 각각
최적화한 뒤 N×M 조립하고 결합 재채점으로 재순위한다. 실제 루프는
vendor/deepB3P/_run_modular.py 가 deepB3P venv 에서 돌고, 여기서는 subprocess 호출·파싱만.
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
    MODULAR_RUNNER,
    TOXINPRED3_PYTHON,
    TOXINPRED3_REPO,
    Settings,
)


@dataclass
class ModularResult:
    shuttles: list[dict] = field(default_factory=list)  # ① [{seq, bbb, charge, muH}]
    linkers: list[dict] = field(default_factory=list)   # ② [{seq, dev, charge, len}]
    best: list[dict] = field(default_factory=list)       # ③ [{shuttle, linker, sequence, bbb, tox, charge, muH, fit, len, linker_len}]
    cargo: str = ""
    n_grid: int = 0


class ModularOptimizer:
    def __init__(self, timeout: float = 1200.0):
        self.timeout = timeout

    def run(self, cargo: str, s_rounds: int = 4, pop: int = 32, elite: int = 6,
            top_shuttles: int = 8, l_rounds: int = 8, top_linkers: int = 6,
            tox_threshold: float = 0.38) -> ModularResult:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "modular.json"
            cmd = [
                str(DEEPB3P_PYTHON), str(MODULAR_RUNNER),
                cargo.upper(), str(s_rounds), str(pop), str(elite), str(top_shuttles),
                str(l_rounds), str(top_linkers), str(tox_threshold),
                str(TOXINPRED3_PYTHON), str(TOXINPRED3_REPO), "toxinpred3.py",
                str(out),
            ]
            proc = subprocess.run(cmd, cwd=str(DEEPB3P_REPO), capture_output=True,
                                  text=True, timeout=self.timeout)
            if proc.returncode != 0 or not out.exists():
                raise RuntimeError(
                    f"모듈별 최적화→조립 실패 (rc={proc.returncode}).\n"
                    f"stderr(마지막 600자): {proc.stderr[-600:]}"
                )
            data = json.loads(out.read_text())
        return ModularResult(shuttles=data.get("shuttles", []), linkers=data.get("linkers", []),
                             best=data.get("best", []), cargo=data.get("cargo", cargo),
                             n_grid=data.get("n_grid", 0))


def get_modular(settings: Settings) -> ModularOptimizer | None:
    return ModularOptimizer() if settings.use_modular_local else None
