"""개발성 기반 잔기 편집 (liability 스크러빙) — directed evolution.

잔기 편집의 **정당한 쓰임**: BBB(오라클 약함·게임 가능)가 아니라 **개발성**(규칙 기반
liability·AGGRESCAN — 거의 "사실"이라 게임 불가)을 목표로, **검증된 기능을 보존하는
보존적 치환**으로 liability를 최소화한다.

설계 원칙:
  · 목표(오라클): developability liability 개수(사실) + AGGRESCAN 응집(출판 척도)
  · 기능 보존 제약: rmt_sim ≥ floor → 검증 리간드에서 멀어지지 않게(Angiopep-7 교훈:
    서열 유사≠기능보장이므로 '보수적 범위'로만 편집). BBB는 라이브러리 스캐폴드가 담당.
  · 적합도 = (liability 개수, 응집) **사전식(lexicographic) 최소화** — 가중치 없음.

⚠️ rmt_sim은 기능 보존의 **프록시**(유사도≠결합보장). 실제 결합 보존은 구조/실험 확정.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field

from .binding import shuttle_similarity
from .developability import assess_developability
from .solubility import assess_solubility

_STD = "ACDEFGHIKLMNPQRSTVWY"


@dataclass
class ScrubResult:
    start_seq: str = ""
    best_seq: str = ""
    start_n: int = 0
    best_n: int = 0
    start_liabilities: list = field(default_factory=list)
    best_liabilities: list = field(default_factory=list)
    start_agg: float = 0.0
    best_agg: float = 0.0
    start_rmt: float = 0.0
    best_rmt: float = 0.0
    n_edits: int = 0


def _fitness(seq: str, rmt_floor: float):
    """(liability 개수, 응집 a3v) 최소화. rmt_sim<floor면 무효(큰 값). (fit, rmt_sim) 반환."""
    sim = shuttle_similarity(seq).score
    if sim < rmt_floor:
        return (999, 999.0), sim
    dev = assess_developability(seq)
    agg = assess_solubility(seq).agg_a3v
    return (dev.n_liabilities, agg), sim


def _mutate(seq: str, rng: random.Random) -> str:
    """단일 치환(길이 보존 → 정렬·기능 안정)."""
    s = list(seq)
    i = rng.randrange(len(s))
    s[i] = rng.choice(_STD)
    return "".join(s)


def scrub_liabilities(shuttle: str, rounds: int = 25, pop: int = 24,
                      elite_frac: float = 0.25, rmt_floor: float = 0.85,
                      seed: int = 0) -> ScrubResult:
    """보존적 치환으로 개발성 liability를 줄인다(기능=rmt_sim 보존)."""
    rng = random.Random(seed)
    dev0 = assess_developability(shuttle)
    agg0 = assess_solubility(shuttle).agg_a3v
    sim0 = shuttle_similarity(shuttle).score

    best_seq, best_fit, best_sim = shuttle, (dev0.n_liabilities, agg0), sim0
    population = [shuttle] * pop
    n_elite = max(2, int(pop * elite_frac))

    for _ in range(rounds):
        scored = []
        for s in population:
            fit, sim = _fitness(s, rmt_floor)
            scored.append((fit, s, sim))
        scored.sort(key=lambda x: x[0])
        if scored[0][0] < best_fit:
            best_fit, best_seq, best_sim = scored[0][0], scored[0][1], scored[0][2]
        elites = [s for _, s, _ in scored[:n_elite]]
        population = elites + [_mutate(rng.choice(elites), rng)
                               for _ in range(pop - n_elite)]

    devb = assess_developability(best_seq)
    n_edits = sum(1 for a, b in zip(shuttle, best_seq) if a != b)
    return ScrubResult(
        start_seq=shuttle, best_seq=best_seq,
        start_n=dev0.n_liabilities, best_n=devb.n_liabilities,
        start_liabilities=dev0.liabilities, best_liabilities=devb.liabilities,
        start_agg=agg0, best_agg=assess_solubility(best_seq).agg_a3v,
        start_rmt=sim0, best_rmt=best_sim, n_edits=n_edits)
