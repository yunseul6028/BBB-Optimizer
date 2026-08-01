"""자율 설계 에이전트 — 공용 도구 백엔드 베이스.

에이전트의 **실제 생물학 계산 도구**(evaluate/structure/generate)를 담는 베이스 클래스.
브레인(LLM) 루프는 `optimizer_agent_gemini.py`의 `GeminiOptimizationAgent`가 이 클래스를
상속해 구현한다. (Claude 브레인 버전은 `with-claude` 브랜치에 보존 — main은 Gemini 전용.)

  - evaluate  : deepB3P(BBB) + ToxinPred3(독성) + 안정성 + 수용체유사 + 개발성 + 선택성 + 용해도 배치 채점
  - structure : ESMFold 폴딩 → 셔틀 구조 노출도 (느림 ~10초, 최종 후보에만)
  - generate  : FBGAN으로 라이브러리 밖 셔틀 생성 (느림, 필요 시)
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .config import (
    LINKER_LIBRARY,
    STANDARD_LINKER_NAME,
    Settings,
    bbb_scoring_seq,
)
from .predictors import get_predictor
from .toxicity import get_toxicity_predictor

MAX_CANDIDATES_PER_STEP = 24


@dataclass
class AgentEvent:
    kind: str  # reasoning|text|plan|evaluation|structure|generation|progress|choice|critique|final|optimum|error
    text: str = ""
    data: dict = field(default_factory=dict)


class OptimizationAgent:
    """도구 백엔드(공용). LLM 루프는 서브클래스(GeminiOptimizationAgent)가 구현."""

    def __init__(self, settings: Settings, max_rounds: int = 8):
        self.settings = settings
        self.max_rounds = max_rounds
        self.predictor = get_predictor(settings)
        self.tox_predictor = get_toxicity_predictor(settings)

    # ---- 도구 백엔드 (실제 생물학 계산) ------------------------------------
    def _evaluate(self, cargo, tool_input):
        items = (tool_input.get("assemblies") or tool_input.get("constructs") or [])[:MAX_CANDIDATES_PER_STEP]
        clean = []
        for c in items:
            lk = "".join(ch for ch in (c.get("prefix", c.get("linker", "")) or "").upper() if ch.isalpha())
            sh = "".join(ch for ch in (c.get("suffix", c.get("shuttle", "")) or "").upper() if ch.isalpha())
            clean.append({"label": c.get("label", "?"), "linker": lk, "shuttle": sh,
                          "sequence": cargo + lk + sh})
        from .binding import shuttle_similarity
        from .stability import assess_stability
        from .developability import assess_developability
        from .selectivity import assess_selectivity
        from .solubility import assess_solubility
        bbb = self.predictor.predict_many([bbb_scoring_seq(cargo, c["linker"], c["shuttle"]) for c in clean])
        tox = self.tox_predictor.predict_many([c["sequence"] for c in clean])
        thr = self.settings.toxicity_threshold
        rows, lines = [], ["label | primary | penalty | index | match | status | assembly"]
        for c, p, t in zip(clean, bbb, tox):
            toxic = t.risk > thr
            bind = shuttle_similarity(c["shuttle"])
            stab = assess_stability(c["sequence"])
            dev = assess_developability(c["sequence"])
            selr = assess_selectivity(c["shuttle"])   # off-target은 셔틀이 주도
            solr = assess_solubility(c["sequence"])   # 용해도는 전체 융합체
            rows.append({"label": c["label"], "linker": c["linker"], "shuttle": c["shuttle"],
                         "sequence": c["sequence"], "bbb": round(p.bbb_permeability, 4),
                         "tox": round(t.risk, 4), "toxic": toxic,
                         "bind_ref": bind.best_ref, "bind_score": bind.score,
                         "instability": stab.instability_index, "stable": stab.stable,
                         "dev_risk": dev.risk_level, "dev_liab": dev.n_liabilities,
                         "dev_charge": dev.net_charge, "dev_agg": dev.agg_score,
                         "dev_liabilities": dev.liabilities,
                         "sel_off": selr.off_target_risk, "selectivity": selr.selectivity,
                         "sel_level": selr.risk_level, "sel_mech": selr.mechanism,
                         "sel_drivers": selr.drivers,
                         "sol_score": solr.score, "sol_level": solr.level})
            lines.append(f"{c['label']} | {p.bbb_permeability:.3f} | {t.risk:.3f} | "
                         f"{stab.instability_index} | {bind.score:.2f} | "
                         f"{'FAIL(penalty)' if toxic else 'ok'} | {c['sequence']}")
        return ("\n".join(lines)
                + f"\n(primary: maximize 0-1 | penalty: ≤{thr:.2f} else FAIL | "
                  "index: lower better, <40 good | match: higher better 0-1)"), rows

    def _structure(self, cargo, tool_input):
        from .structure import analyze_construct
        lk = "".join(ch for ch in (tool_input.get("prefix", tool_input.get("linker", "")) or "").upper() if ch.isalpha())
        sh = "".join(ch for ch in (tool_input.get("suffix", tool_input.get("shuttle", "")) or "").upper() if ch.isalpha())
        sr = analyze_construct(cargo, lk, sh)
        if sr.error:
            return f"inspect failed: {sr.error}", {"error": sr.error}
        text = (f"exposure={sr.shuttle_exposure} (>=0.35 = exposed), "
                f"confidence(suffix/all)={sr.shuttle_plddt}/{sr.mean_plddt}. "
                f"verdict: {'exposed' if sr.exposed else 'buried/low-confidence'}")
        return text, {"linker": lk, "shuttle": sh, "exposure": sr.shuttle_exposure,
                      "shuttle_plddt": sr.shuttle_plddt, "mean_plddt": sr.mean_plddt,
                      "verdict": sr.verdict, "exposed": sr.exposed}

    def _generate(self, cargo, tool_input):
        from .generative import get_fbgan
        fb = get_fbgan(self.settings)
        if fb is None:
            return "propose unavailable (no local generator).", {"novel": []}
        rounds = max(2, min(4, int(tool_input.get("rounds", 3))))
        linker = LINKER_LIBRARY[STANDARD_LINKER_NAME]["seq"]
        fres = fb.run(cargo, linker, rounds=rounds, tox_threshold=self.settings.toxicity_threshold)
        best = fres.best[:5]
        lines = ["novel suffixes (with standard prefix):", "suffix | primary | penalty"]
        for b in best:
            lines.append(f"{b['shuttle']} | {b['bbb']:.3f} | {b['tox']:.3f}")
        return "\n".join(lines), {"novel": best}
