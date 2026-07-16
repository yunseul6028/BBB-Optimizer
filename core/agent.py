"""최적화 에이전트 (오케스트레이션 레이어).

책임:
    1. 화물(cargo) 펩타이드 수령
    2. 링커 라이브러리 × 셔틀 라이브러리를 **전수 조합**해 융합체(construct) 조립
    3. 화물 단독 기준선 + 모든 조합을 BBB·독성 예측기로 배치 스코어링
    4. 독성 임계값 필터링 후, 투과율 높고 비독성인 상위 조합 선정
    5. 추론 과정을 AgentStep 스트림으로 방출

construct = cargo + linker + shuttle  (N→C 방향).
UI는 run()이 yield하는 AgentStep을 렌더링하고, 종료 시 StopIteration.value로
최종 OptimizationResult를 받는다.
"""

from __future__ import annotations

from typing import Generator

from .config import (
    LINKER_LIBRARY,
    MODEL_MAX_LEN,
    SHUTTLES,
    STANDARD_LINKER_NAME,
    Settings,
    bbb_scoring_seq,
    junction_window,
)
from .predictors import BBBPredictor, get_predictor
from .toxicity import ToxicityPredictor, get_toxicity_predictor
from .schemas import (
    AgentStep,
    Construct,
    EvaluatedCandidate,
    OptimizationResult,
    StepLevel,
    Verdict,
)


def assemble(cargo: str, linker: str, shuttle_seq: str) -> str:
    """융합체 = 화물 + 링커 + 셔틀 (N→C)."""
    return cargo + linker + shuttle_seq


def build_construct(cargo: str, linker_name: str, linker_seq: str,
                    shuttle_name: str, shuttle_seq: str) -> Construct:
    seq = assemble(cargo, linker_seq, shuttle_seq)
    return Construct(
        cargo=cargo,
        linker_name=linker_name,
        linker=linker_seq,
        shuttle_name=shuttle_name,
        shuttle_seq=shuttle_seq,
        sequence=seq,
        bbb_sequence=bbb_scoring_seq(cargo, linker_seq, shuttle_seq),
        is_standard_linker=(linker_name == STANDARD_LINKER_NAME),
        truncated=len(seq) > MODEL_MAX_LEN,
    )


class OptimizerAgent:
    def __init__(
        self,
        predictor: BBBPredictor,
        tox_predictor: ToxicityPredictor,
        toxicity_threshold: float,
        linkers: list[tuple[str, str]],    # [(name, seq), ...]
        shuttles: list[tuple[str, str]],   # [(name, seq), ...]
        top_n: int = 3,
    ):
        self.predictor = predictor
        self.tox_predictor = tox_predictor
        self.toxicity_threshold = toxicity_threshold
        self.linkers = linkers
        self.shuttles = shuttles
        self.top_n = top_n

    def run(self, cargo: str) -> Generator[AgentStep, None, OptimizationResult]:
        cargo = cargo.strip().upper()

        # --- 1단계: 화물 수령 -----------------------------------------------
        yield AgentStep("1단계",
                        f"화물(cargo) 펩타이드 수령 ({len(cargo)}aa). 라이브러리 전수 조립 준비.",
                        StepLevel.INFO)

        # --- 조립: 링커 × 셔틀 전수 매트릭스 --------------------------------
        constructs = [
            build_construct(cargo, ln, ls, sn, ss)
            for sn, ss in self.shuttles
            for ln, ls in self.linkers
        ]
        cargo_only_c = Construct(
            cargo=cargo, linker_name="", linker="", shuttle_name="(없음)",
            shuttle_seq="", sequence=cargo, bbb_sequence=junction_window(cargo),
            truncated=len(cargo) > MODEL_MAX_LEN,
        )
        yield AgentStep(
            "2단계",
            f"라이브러리 전수 조립: 링커 {len(self.linkers)}종 × 셔틀 {len(self.shuttles)}종 "
            f"= 융합체 {len(constructs)}개 (+화물 단독 기준선).",
            StepLevel.INFO,
        )

        # --- 배치 예측: BBB + 독성 (기준선 + 모든 조합 한 번에) --------------
        all_c = [cargo_only_c] + constructs
        # BBB는 연결부위 윈도우(bbb_sequence), 독성은 전체 서열로 평가
        preds = self.predictor.predict_many([c.bbb_sequence for c in all_c])
        tox_results = self.tox_predictor.predict_many([c.sequence for c in all_c])
        cargo_pred, construct_preds = preds[0], preds[1:]
        cargo_tox, construct_tox = tox_results[0], tox_results[1:]

        cargo_pred.toxicity_risk = cargo_tox.risk
        yield AgentStep(
            "3단계",
            f"deepB3P(BBB) + {self.tox_predictor.source}(독성)로 {len(all_c)}개 서열 배치 실측 완료. "
            f"화물 단독 BBB {cargo_pred.bbb_permeability * 100:.1f}% · "
            f"독성 {cargo_tox.risk * 100:.0f}%{' ⚠️Toxin' if cargo_tox.is_toxic else ''}.",
            StepLevel.INFO,
        )
        cargo_only = EvaluatedCandidate(cargo_only_c, cargo_pred, Verdict.REFERENCE)

        truncated = [c for c in constructs if c.truncated]
        if truncated:
            yield AgentStep(
                "3단계",
                f"ℹ️ {len(truncated)}개 융합체가 {MODEL_MAX_LEN}aa 초과 → BBB는 **연결부위"
                f"(링커+셔틀 = BBB 모듈)** 만으로 계산 (화물 제외, 독성은 전체).",
                StepLevel.INFO,
            )

        # --- 4단계: 독성 필터링 (요약 보고) ---------------------------------
        evaluated: list[EvaluatedCandidate] = []
        rejected: list[str] = []
        for c, pred, tox in zip(constructs, construct_preds, construct_tox):
            pred.toxicity_risk = tox.risk
            pred.note = f"{pred.note} · 독성 {tox.source}"
            cand = EvaluatedCandidate(c, pred)
            if tox.is_toxic:
                cand.verdict = Verdict.REJECTED_TOXIC
                rejected.append(f"{c.label}({tox.risk:.2f})")
            evaluated.append(cand)

        n_survive = len(evaluated) - len(rejected)
        yield AgentStep(
            "4단계",
            f"부작용 스크리닝({self.tox_predictor.source}) 완료 — 임계값 "
            f"{self.toxicity_threshold:.2f} 초과 {len(rejected)}개 탈락, {n_survive}개 생존."
            + (f" 탈락: {', '.join(rejected[:6])}{' …' if len(rejected) > 6 else ''}"
               if rejected else ""),
            StepLevel.WARN if rejected else StepLevel.INFO,
        )

        # --- 최종: 생존 조합 중 BBB 상위 N개 -------------------------------
        survivors = sorted(
            (c for c in evaluated if c.verdict != Verdict.REJECTED_TOXIC),
            key=lambda c: c.prediction.bbb_permeability, reverse=True,
        )
        top = survivors[: self.top_n]
        winner = top[0] if top else None
        if winner:
            winner.verdict = Verdict.ACCEPTED
            names = " > ".join(f"{c.construct.label}({c.prediction.bbb_permeability*100:.0f}%)"
                               for c in top)
            yield AgentStep("최종", f"🏆 베스트 {len(top)}: {names}", StepLevel.SUCCESS)
        else:
            yield AgentStep("최종", "생존 조합 없음 — 임계값/라이브러리 재검토 필요.",
                            StepLevel.ERROR)

        return OptimizationResult(
            cargo=cargo,
            cargo_only=cargo_only,
            winner=winner,
            candidates=evaluated,
            n_linkers=len(self.linkers),
            n_shuttles=len(self.shuttles),
            truncated_any=bool(truncated),
        )


# ---------------------------------------------------------------------------
# 조립(팩토리) — 라이브러리 전체를 주입.
# ---------------------------------------------------------------------------
def build_agent(settings: Settings) -> OptimizerAgent:
    linkers = [(name, spec["seq"]) for name, spec in LINKER_LIBRARY.items()]
    shuttles = [(name, spec["seq"]) for name, spec in SHUTTLES.items()]
    return OptimizerAgent(
        predictor=get_predictor(settings),
        tox_predictor=get_toxicity_predictor(settings),
        toxicity_threshold=settings.toxicity_threshold,
        linkers=linkers,
        shuttles=shuttles,
        top_n=settings.top_n,
    )
