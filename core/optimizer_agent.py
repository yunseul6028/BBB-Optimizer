"""자율 설계 에이전트 — 공용 도구 백엔드 베이스.

에이전트의 **실제 생물학 계산 도구**(evaluate/structure)를 담는 베이스 클래스.
브레인(LLM) 루프는 `optimizer_agent_gemini.py`의 `GeminiOptimizationAgent`가 이 클래스를
상속해 구현한다. (Claude 브레인 버전은 `with-claude` 브랜치에 보존 — main은 Gemini 전용.)

  - evaluate  : deepB3P(BBB) + ToxinPred3(독성) + 안정성 + 수용체유사 + 개발성 + 선택성 + 용해도 배치 채점
  - structure : ESMFold 폴딩 → 셔틀 구조 노출도 (느림 ~10초, 최종 후보에만)

셔틀은 **de-novo 생성하지 않는다**(검증 리간드만). 라이브러리 밖 탐색은 잔기 편집(design_candidate)
또는 링커·셔틀 co-evolution(라이브러리 시드 진화)로 수행한다.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .config import (
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
        from .delivery import assess_delivery
        # 융합 연결부위 deepB3P + 셔틀 단독 deepB3P를 **한 배치**로 추론(전달 축 분해용).
        junction_seqs = [bbb_scoring_seq(cargo, c["linker"], c["shuttle"]) for c in clean]
        uniq_shuttles: list[str] = []
        for c in clean:
            if c["shuttle"] and c["shuttle"] not in uniq_shuttles:
                uniq_shuttles.append(c["shuttle"])
        preds = self.predictor.predict_many(junction_seqs + uniq_shuttles)
        bbb = preds[: len(clean)]
        shuttle_bbb = {s: preds[len(clean) + i].bbb_permeability for i, s in enumerate(uniq_shuttles)}
        tox = self.tox_predictor.predict_many([c["sequence"] for c in clean])
        thr = self.settings.toxicity_threshold
        rows, lines = [], ["label | primary | shuttleBBB | mech | preserv | penalty | index | match | status | assembly"]
        for c, p, t in zip(clean, bbb, tox):
            toxic = t.risk > thr
            bind = shuttle_similarity(c["shuttle"])
            stab = assess_stability(c["sequence"])
            dev = assess_developability(c["sequence"])
            selr = assess_selectivity(c["shuttle"])   # off-target은 셔틀이 주도
            solr = assess_solubility(c["sequence"])   # 용해도는 전체 융합체
            # 전달 축 분해 — deepB3P를 '셔틀 내재 × 융합 보존 × 메커니즘 타당도'로 나눈다.
            si = shuttle_bbb.get(c["shuttle"], 0.0)
            dax = assess_delivery(si, p.bbb_permeability, target=bind.target, mechanism=bind.mechanism)
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
                         "sol_score": solr.score, "sol_level": solr.level,
                         "shuttle_bbb": round(si, 4), "mechanism": dax.mechanism,
                         "is_rmt": dax.is_rmt, "preservation": dax.preservation,
                         "avidity": dax.avidity, "deepb3p_valid": dax.deepb3p_validity,
                         "delivery_basis": dax.basis})
            _pres = f"{dax.preservation:.2f}" if dax.preservation is not None else "n/a"
            lines.append(f"{c['label']} | {p.bbb_permeability:.3f} | {si:.3f} | {dax.mechanism} | "
                         f"{_pres} | {t.risk:.3f} | {stab.instability_index} | {bind.score:.2f} | "
                         f"{'FAIL(penalty)' if toxic else 'ok'} | {c['sequence']}")
        return ("\n".join(lines)
                + f"\n(primary: deepB3P 융합점수 0-1 — **RMT 셔틀엔 약한 프록시**(수용체 결합이 병목) | "
                  "shuttleBBB: 셔틀 단독 deepB3P | mech: RMT(수용체매개)/CPP(막투과) | "
                  f"preserv: 융합 보존도(융합/셔틀) | penalty: ≤{thr:.2f} else FAIL | "
                  "index: lower better, <40 good | match: 수용체유사도 0-1 | "
                  "→ RMT는 primary 절대값보다 mech·preserv·구조노출로 판단)"), rows

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

    def _coevolve(self, cargo, tool_input):
        """검증 라이브러리 셔틀·링커를 **서열 directed evolution**(라이브러리 시드 point-mutation/
        crossover)으로 진화시켜 라이브러리 밖 후보를 제안하고, 그 결과를 **8축 평가 파이프라인으로
        재채점**한다(선택성·delivery 분해·eff_bbb 랭킹에 그대로 합류). 셔틀은 de-novo 생성이 아니라
        검증 리간드에서 파생된 변이체다."""
        if not self.settings.use_coevo_local:
            return "evolve unavailable (no local co-evolution engine).", []
        from .coevolution import CoevolutionOptimizer
        rounds = max(2, min(3, int(tool_input.get("rounds", 2))))
        try:
            # 라이브 세션용 경량 파라미터(짧은 탐색 버스트) + 짧은 타임아웃
            res = CoevolutionOptimizer(timeout=300.0).run(
                cargo, rounds=rounds, pop=16, elite=4,
                tox_threshold=self.settings.toxicity_threshold)
        except Exception as exc:  # noqa: BLE001
            return f"evolve failed: {type(exc).__name__}: {exc}", []
        pairs = res.best[:6]
        if not pairs:
            return "co-evolution produced no non-toxic candidates.", []
        # 진화된 (링커, 셔틀) 쌍을 8축 평가로 재채점 → 후보 풀에 합류
        constructs = [{"label": f"진화#{i + 1}", "linker": p.get("linker", ""),
                       "shuttle": p.get("shuttle", "")} for i, p in enumerate(pairs)]
        text, rows = self._evaluate(cargo, {"constructs": constructs})
        return ("directed evolution (라이브러리 시드) → 8축 재평가:\n" + text), rows
