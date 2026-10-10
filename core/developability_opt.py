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


# ---------------------------------------------------------------------------
# 다목적(Pareto) 개발성 최적화 — "동시 최적화"를 진짜로 (상충하는 신뢰 오라클 축)
# ---------------------------------------------------------------------------
# 목표(모두 최소화·전부 신뢰 가능한 오라클): liability 개수(규칙=사실) · 응집(AGGRESCAN)
# · 불안정성(Guruprasad). 서로 상충하므로 단일 최적해가 아니라 Pareto 전선을 찾는다.
OBJECTIVE_NAMES = ("liability", "aggregation", "instability")


def _instability(seq: str) -> float:
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        return round(float(ProteinAnalysis(seq).instability_index()), 1)
    except Exception:  # noqa: BLE001
        return 0.0


def objectives(seq: str) -> tuple:
    """개발성 다목적 벡터 (모두 ↓ 좋음): (liability 개수, 응집 a3v, 불안정성지수)."""
    dev = assess_developability(seq)
    agg = round(assess_solubility(seq).agg_a3v, 3)
    return (dev.n_liabilities, agg, _instability(seq))


def _dominates(a: tuple, b: tuple) -> bool:
    """a가 b를 지배: 모든 목표에서 a ≤ b, 적어도 하나는 a < b."""
    return all(x <= y for x, y in zip(a, b)) and any(x < y for x, y in zip(a, b))


def pareto_scrub(shuttle: str, rounds: int = 40, pop: int = 30,
                 rmt_floor: float = 0.85, seed: int = 0, max_archive: int = 40):
    """다목적 개발성 최적화 — 기능(rmt_sim≥floor) 보존하며 Pareto 전선을 반환.

    반환: [{seq, objectives, rmt_sim, n_edits}] — 비지배(non-dominated) 후보 집합.
    상충하는 목표(liability↓가 응집↑을 부를 수 있음)를 동시에 다루므로 단일 '정답'이
    아니라 **트레이드오프 전선**을 제시한다.
    """
    rng = random.Random(seed)
    sim0 = shuttle_similarity(shuttle).score
    # (seq, objvec, rmt_sim)
    archive = [(shuttle, objectives(shuttle), sim0)]

    def _try_insert(s):
        sim = shuttle_similarity(s).score
        if sim < rmt_floor:
            return
        obj = objectives(s)
        if any(_dominates(a[1], obj) for a in archive):
            return                                   # 기존에 지배당함
        archive[:] = [a for a in archive if not _dominates(obj, a[1])]  # 지배되는 것 제거
        if all(a[0] != s for a in archive):
            archive.append((s, obj, sim))

    population = [shuttle] * pop
    for _ in range(rounds):
        for s in population:
            _try_insert(s)
        parents = [a[0] for a in archive] or [shuttle]
        population = [_mutate(rng.choice(parents), rng) for _ in range(pop)]
        if len(archive) > max_archive:               # 과대 시 다양성 유지하며 절삭
            archive.sort(key=lambda a: a[1])
            del archive[max_archive:]

    front = sorted(archive, key=lambda a: a[1])
    return [{"seq": s, "objectives": dict(zip(OBJECTIVE_NAMES, obj)),
             "rmt_sim": round(sim, 2),
             "n_edits": sum(1 for a, b in zip(shuttle, s) if a != b)}
            for s, obj, sim in front]
