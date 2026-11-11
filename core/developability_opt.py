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
    """(liability 개수, 응집 **핫스팟** a3v) 최소화. rmt_sim<floor면 무효. (fit, rmt_sim) 반환.

    ⚠️ 응집 목표를 **전역 평균 → 국소 최대(핫스팟)** 로 바꿨다(2026-11-10). 전역 평균은
    길이로 희석되고 이 최적화기에는 **삽입 변이**(`_try_insert`)가 있어서, 평균을 목표로
    두면 "저응집 잔기를 끼워넣기"가 가장 싼 해가 된다 — 응집 위험 구간은 그대로 둔 채
    점수만 내려간다. M10 실측의 길이 아티팩트(링커 길이↔응집 순위상관 평균 -0.697)와
    **같은 함정**이다. 핫스팟은 국소 최대라 삽입으로 낮출 수 없다.
    """
    sim = shuttle_similarity(seq).score
    if sim < rmt_floor:
        return (999, 999.0), sim
    dev = assess_developability(seq)
    agg = dev.agg_hotspot_a3v
    return (dev.n_liabilities, agg), sim


def _mutate(seq: str, rng: random.Random, frozen: frozenset[int] = frozenset()) -> str:
    """단일 치환(길이 보존 → 정렬·기능 안정). `frozen` = 건드리지 않을 0-based 인덱스.

    `frozen` 이 왜 필요한가: `rmt_sim ≥ floor` 는 **전역 평균**이라 치명적 치환 2~3개가
    동일 잔기 뒤에 희석된다. 실측으로 옵티마이저는 off-target 을 낮추려고 Angiopep-2 의
    K10·K15(문헌상 기능 필수)를 지우는 법을 배웠다 — 자세한 근거는
    `core.evidence.CONSERVED_POSITIONS` 주석.
    """
    s = list(seq)
    choices = [i for i in range(len(s)) if i not in frozen]
    if not choices:
        return seq
    i = rng.choice(choices)
    s[i] = rng.choice(_STD)
    return "".join(s)


def _frozen_idx(shuttle: str) -> frozenset[int]:
    """이 셔틀에서 **동결할 0-based 인덱스** (문헌이 집어 준 곳만)."""
    try:
        from .evidence import conserved_positions
        return frozenset(p - 1 for p in conserved_positions(shuttle)
                         if 1 <= p <= len(shuttle))
    except Exception:  # noqa: BLE001
        return frozenset()


def scrub_liabilities(shuttle: str, rounds: int = 25, pop: int = 24,
                      elite_frac: float = 0.25, rmt_floor: float = 0.85,
                      seed: int = 0) -> ScrubResult:
    """보존적 치환으로 개발성 liability를 줄인다(기능=rmt_sim 보존)."""
    rng = random.Random(seed)
    dev0 = assess_developability(shuttle)
    agg0 = dev0.agg_hotspot_a3v
    sim0 = shuttle_similarity(shuttle).score

    frozen = _frozen_idx(shuttle)
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
        population = elites + [_mutate(rng.choice(elites), rng, frozen)
                               for _ in range(pop - n_elite)]

    devb = assess_developability(best_seq)
    n_edits = sum(1 for a, b in zip(shuttle, best_seq) if a != b)
    return ScrubResult(
        start_seq=shuttle, best_seq=best_seq,
        start_n=dev0.n_liabilities, best_n=devb.n_liabilities,
        start_liabilities=dev0.liabilities, best_liabilities=devb.liabilities,
        start_agg=agg0, best_agg=devb.agg_hotspot_a3v,
        start_rmt=sim0, best_rmt=best_sim, n_edits=n_edits)


# ---------------------------------------------------------------------------
# 다목적(Pareto) 개발성 최적화 — "동시 최적화"를 진짜로 (상충하는 신뢰 오라클 축)
# ---------------------------------------------------------------------------
# 목표(모두 최소화·전부 신뢰 가능한 오라클): liability 개수(규칙=사실) · 응집(AGGRESCAN)
# · 불안정성(Guruprasad). 서로 상충하므로 단일 최적해가 아니라 Pareto 전선을 찾는다.
#: 🔴 off-target 은 **목적이 아니라 제약**이다(2026-11-11 대응 널 심사, §26).
#: 축을 심사해 보니 *악화는 믿을 만하지만 개선은 못 믿는다*:
#:   · 조성 보존 셔플(전하·GRAVY 정확히 불변, 결합만 파괴) → off 가 **올바르게 상승**
#:     (Angiopep-2 0.09→0.214 · B6 0.37→0.853, 200/200) ⇒ **실제 악화는 검출한다**
#:   · 같은 크기 무작위 치환 → k=3 에서 **46% 가 off 를 내린다**(바닥 셔틀 제외 414↑/355↓)
#:     ⇒ **위험이 줄지 않았는데 점수가 내려가는 길이 널려 있다**
#:   · 전하만 지운 CPP → TAT 0.95→**0.00** · R9 0.95→**0.00** ⇒ **해킹 1수로 뚫린다**
#: 그래서 off 는 `objectives()` 에서 빼고 **양측 밴드 제약**으로 쓴다 — 시드 off 에서
#: `OFF_BAND` 이상 **벗어나면**(올라가든 내려가든) 거부. 즉 축을 "검증된 셔틀의 프로파일에서
#: 얼마나 벗어났나"의 **변화 검출기**로만 쓰고, 어느 방향으로도 **보상으로 쓰지 않는다.**
#:
#: ⚠️ **단측 상한(`off ≤ 시드`)은 실패했다** — 실측하니 오히려 전하 삭제를 **선택**했다
#: (R11 삭제 54%: 제약을 만족하는 가장 쉬운 길이 전하 삭제이므로 생존자가 그쪽으로 쏠린다).
#: 밴드 폭은 **측정으로 골랐다**(§26-5, Angiopep-2 · 시드 off 0.09 · 5시드):
#:
#:   | 제약 | 전선 | R11 삭제 | off 범위 | 최선 목적 |
#:   |---|---|---|---|---|
#:   | 없음(옛) | 55 | 15% | 0.05–**0.25** | (1, 0.077, 45.4) |
#:   | 단측 `≤` | 26 | **54%** | 0.00–0.06 | (1, 0.077, 45.4) |
#:   | **밴드 ±0.05** | 35 | **20%** | 0.05–0.12 | (1, 0.161, 38.9) |
#:   | 밴드 ±0.02 | 30 | 0% | 0.10–0.11 | (1, 0.237, 63.0) |
#:   | 밴드 ±0.00 | 11 | 0% | 0.09 | (2, 0.651, 57.1) |
#:
#: ±0.05 를 골랐다: M4 의 실패(off 0.09→0.25)를 막으면서 전하 삭제가 **압력 없는 기저선
#: (15%)** 수준으로 돌아오고, 개발성 이득을 대부분 지킨다. 더 좁히면 해킹은 0% 가 되지만
#: 응집·불안정성 대가가 급격히 커진다 — **판단이고, 상수 하나로 되돌릴 수 있다.**
OFF_BAND = 0.05
OBJECTIVE_NAMES = ("liability", "aggregation", "instability")


def _instability(seq: str) -> float:
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        return round(float(ProteinAnalysis(seq).instability_index()), 1)
    except Exception:  # noqa: BLE001
        return 0.0


def objectives(seq: str) -> tuple:
    """다목적 벡터 (모두 ↓ 좋음): (liability 개수, 응집 핫스팟 a3v, 불안정성지수).

    off-target 은 왜 **여기 없고 제약인가** (§25 → §26)
    -----------------------------------------------------------------
    옛 벡터는 3축(liability·응집·불안정성)이었고 **off-target 이 빠져 있었다**. 그런데
    순위(`optimizer_agent.risk_axes`)는 **4축**이다. 이 어긋남에 측정된 대가가 있었다 —
    M4 콜д 3세션의 설계체가 **전부** off 0.09 → 0.20~0.23 으로 악화했고(원장 기록 5/5),
    그래서 "전선 멤버를 지배하는 설계체 0개"가 됐다. 옵티마이저가 **못 보는 축으로
    값을 치르고 있었던 것**이다. 원칙: **탐색의 목적함수는 심사의 축과 같아야 한다.**

    1차 처방은 off 를 **목적에 넣는 것**이었고, 그러자 옵티마이저가 off 식의 양전하 항을
    타고 **문헌이 기능이라 말한 라이신(K10·K15)을 지우는 법을 배웠다**(§25). 보존 위치를
    동결해 막으니 **R11(문헌 미기재)로 옮겨 갔다.** 그래서 축 자체를 대응 널로 심사했고
    (§26), **악화 검출은 신뢰할 수 있지만 개선 신호는 신뢰할 수 없다**는 결론이 나왔다
    (셔플 널 통과 · 무작위 편집 46% 가 off 를 내림 · 전하 제거 1수로 TAT 0.95→0.00).
    따라서 off 는 **양측 밴드 제약**으로만 쓴다(`OFF_BAND`) — 시드에서 벗어나면 어느
    방향이든 거부. 이게 M4 의 실패(설계체 off 0.09→0.20~0.25, 5/5)를 막으면서 해킹
    보상도 없앤다. 단측 상한은 **오히려 전하 삭제를 선택**해서 실패했다(§26-5).
    ⚠️ BBB 는 목적도 제약도 아니다 — 순위에서 뺀 축이고(M6·M8), 서브프로세스가 돈다.

    ⚠️ 불안정성은 **연속값**을 유지한다(순위는 임계 40 이진). 의도된 불일치다 —
    **탐색에는 기울기가 필요하고**(40 아래로 내려가는 길을 찾아야 한다), 이진화는
    길이 해킹을 막는 **판정 장치**다. 대신 전선을 에이전트에게 줄 때 순위 축으로
    중복을 제거해, 순위 관점에서 같은 점을 여러 개 제시하지 않는다.
    """
    dev = assess_developability(seq)
    agg = round(dev.agg_hotspot_a3v, 3)
    return (dev.n_liabilities, agg, _instability(seq))


def off_target(seq: str) -> float:
    """셔틀의 off-target 위험 — `_evaluate` 와 **같은 호출**(셔틀만 본다 → 링커 불변).

    전하·GRAVY·BLOSUM 의 로컬 결정론 계산이라 서브프로세스가 **0회**다.
    """
    try:
        from .selectivity import assess_selectivity
        return round(float(assess_selectivity(seq).off_target_risk), 3)
    except Exception:  # noqa: BLE001
        # 🔴 조용히 0.0 으로 떨어지면 "off 가 완벽하다"는 **거짓 신호**가 되어 전선을
        #    오염시킨다. 최악값으로 떨어뜨려 그 후보가 전선에 못 들어가게 한다.
        return 1.0


def _dominates(a: tuple, b: tuple) -> bool:
    """a가 b를 지배: 모든 목표에서 a ≤ b, 적어도 하나는 a < b."""
    return all(x <= y for x, y in zip(a, b)) and any(x < y for x, y in zip(a, b))


def pareto_scrub(shuttle: str, rounds: int = 40, pop: int = 30,
                 rmt_floor: float = 0.85, seed: int = 0, max_archive: int = 40):
    """다목적 개발성 최적화 — 기능(rmt_sim≥floor) 보존하며 Pareto 전선을 반환.

    반환: [{seq, objectives, off_target, off_seed, rmt_sim, n_edits}] — 비지배 후보 집합.
    제약 셋: `rmt_sim ≥ rmt_floor` · 문헌 보존 위치 동결 · **|off − 시드 off| ≤ OFF_BAND**.
    상충하는 목표(liability↓가 응집↑을 부를 수 있음)를 동시에 다루므로 단일 '정답'이
    아니라 **트레이드오프 전선**을 제시한다.
    """
    rng = random.Random(seed)
    frozen = _frozen_idx(shuttle)
    sim0 = shuttle_similarity(shuttle).score
    # 🔴 off-target 은 **양측 밴드 제약**이다(§26) — 시드에서 ±OFF_BAND 를 벗어나면 거부.
    #    목적에 넣으면 축의 양전하 항이 해킹 보상이 되고(§25: K10/K15 → 막으면 R11),
    #    단측 상한으로 두면 제약 만족이 곧 전하 삭제라서 **오히려 더 쏠린다**(§26-5).
    off0 = off_target(shuttle)
    # (seq, objvec, rmt_sim)
    archive = [(shuttle, objectives(shuttle), sim0)]

    def _try_insert(s):
        sim = shuttle_similarity(s).score
        if sim < rmt_floor:
            return
        if abs(off_target(s) - off0) > OFF_BAND + 1e-9:   # off 밴드 제약 — §26
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
        population = [_mutate(rng.choice(parents), rng, frozen) for _ in range(pop)]
        if len(archive) > max_archive:               # 과대 시 다양성 유지하며 절삭
            archive.sort(key=lambda a: a[1])
            del archive[max_archive:]

    front = sorted(archive, key=lambda a: a[1])
    return [{"seq": s, "objectives": dict(zip(OBJECTIVE_NAMES, obj)),
             # 목적이 아니라 **충족된 제약**으로 보고한다(시드 ± OFF_BAND).
             "off_target": off_target(s), "off_seed": off0,
             "rmt_sim": round(sim, 2),
             "n_edits": sum(1 for a, b in zip(shuttle, s) if a != b)}
            for s, obj, sim in front]
