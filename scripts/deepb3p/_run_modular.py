# -- coding: utf-8 --
"""모듈별 최적화 → 조립 재순위 (decoupled optimization + assembly).

co-evolution(개체=(셔틀,링커) 동시 진화)의 대안. 각 모듈을 **자기 목적함수가 유효한
자기 공간**에서 따로 최적화한 뒤 조립하고, 마지막에 결합 재채점으로 궁합(상호작용)을 회수한다.

  ① 셔틀 공간: 잠재 z 진화 → BBB↑·전하 억제 상위 N개  (deepB3P in-process, 독성은 조립에서)
  ② 링커 공간: 링커 유전자 진화 → **개발성** 상위 M개
       개발성 = 유연성(G/S/A/P)↑ · 프로테아제부위(KR/RR/RK/KK)↓ · 저전하 · 저응집 · 컴팩트
       (링커는 BBB가 없으므로 deepB3P 아님 — 규칙 기반, 순수 파이썬)
  ③ 조립: N×M 전체 융합체 → BBB(연결부위)·독성(전체)·전하 가드레일로 결합 재채점 → 재순위

효율: 느린 독성 예측(ToxinPred3 subprocess)을 **조립 단계 1회**만 호출(co-evolution은 매 라운드
호출)하고, 링커 최적화는 예측기 없이 즉시. deepB3P는 1회 로딩.

Usage:
  python _run_modular.py <cargo> <s_rounds> <pop> <elite> <topN> <l_rounds> <topM> <tox_thr> \
      <toxpy> <toxrepo> <toxscript> <out.json>
"""
import json
import sys

import numpy as np
import functools
import torch

torch.load = functools.partial(torch.load, map_location="cpu")

from _run_fbgan import HIDDEN, N_CHARS, SEQ_LEN, BBBScorer, bbb_region, decode, score_tox
from _run_coevo import LINK_SEEDS, cross_linker, mutate_linker
from _physchem import charge_guardrail, net_charge_ph74
from fbgan.models import Generator

_HYDROPHOBIC = set("FILVWYM")


# ---------------- ② 링커 개발성 목적함수 (규칙 기반, 예측기 불요) ----------------
def _protease_sites(lk: str) -> int:
    """furin/트립신 유사 이염기 절단부위(KR·RR·RK·KK) 수."""
    return sum(1 for i in range(len(lk) - 1) if lk[i] in "KR" and lk[i + 1] in "KR")


def linker_dev_fitness(lk: str) -> float:
    """링커 개발성 점수(높을수록 좋음) — 유연·비절단·저전하·저응집·컴팩트."""
    n = max(len(lk), 1)
    flex = sum(c in "GSAP" for c in lk) / n            # 유연 잔기 비율
    hydro = sum(c in _HYDROPHOBIC for c in lk) / n     # 응집 위험
    return (1.0 * flex
            - 0.3 * _protease_sites(lk)
            - 0.15 * abs(net_charge_ph74(lk))
            - 0.5 * hydro
            - 0.02 * n)


# ---------------- ① 셔틀 공간 진화 (BBB·전하, 독성은 조립에서) ----------------
def evolve_shuttles(G, bbb, rounds, pop, elite, topN):
    rng = np.random.RandomState(2022)
    z = rng.randn(pop, 128)
    sigma = 0.6
    best = {}   # seq -> (fit, bbb, charge, muH)
    for _ in range(rounds):
        shuttles = decode(G, z)
        bscore = bbb.score(shuttles)                   # 셔틀 단독 BBB
        guards = [charge_guardrail(s) for s in shuttles]
        fit = np.array([b - g[0] for b, g in zip(bscore, guards)])
        for i, s in enumerate(shuttles):
            prev = best.get(s)
            if prev is None or fit[i] > prev[0]:
                best[s] = (float(fit[i]), float(bscore[i]), guards[i][1], guards[i][2])
        order = np.argsort(fit)[::-1]
        elites = z[order[:elite]]
        children = []
        while len(children) < pop - elite - 2:
            children.append(elites[rng.randint(elite)] + sigma * rng.randn(128))
        z = np.vstack([elites, np.array(children), rng.randn(2, 128)])
        sigma *= 0.85
    top = sorted(best.items(), key=lambda kv: kv[1][0], reverse=True)[:topN]
    return [{"seq": s, "bbb": round(v[1], 4), "charge": v[2], "muH": v[3]} for s, v in top]


# ---------------- ② 링커 공간 진화 (개발성) ----------------
def evolve_linkers(rounds, pop, elite, topM):
    rng = np.random.RandomState(7)
    L = [LINK_SEEDS[i % len(LINK_SEEDS)] for i in range(pop)]
    best = {}
    for _ in range(rounds):
        fits = [linker_dev_fitness(lk) for lk in L]
        for lk, f in zip(L, fits):
            if lk not in best or f > best[lk]:
                best[lk] = f
        order = sorted(range(pop), key=lambda i: fits[i], reverse=True)
        elites = [L[i] for i in order[:elite]]
        children = []
        while len(children) < pop - elite - 2:
            a, b = elites[rng.randint(elite)], elites[rng.randint(elite)]
            children.append(mutate_linker(cross_linker(a, b, rng), rng))
        L = elites + children + [LINK_SEEDS[rng.randint(len(LINK_SEEDS))] for _ in range(2)]
    top = sorted(best.items(), key=lambda kv: kv[1], reverse=True)[:topM]
    return [{"seq": lk, "dev": round(v, 4), "charge": net_charge_ph74(lk), "len": len(lk)}
            for lk, v in top]


def main(cargo, s_rounds, pop, elite, topN, l_rounds, topM, tox_thr,
         toxpy, toxrepo, toxscript, out_path):
    G = Generator(n_chars=N_CHARS, seq_len=SEQ_LEN, bs=pop, hidden=HIDDEN)
    G.load_state_dict(torch.load("fbgan/checkpoint/G_weights_1000.pth"))
    bbb = BBBScorer()

    shuttles = evolve_shuttles(G, bbb, s_rounds, pop, elite, topN)   # ①
    linkers = evolve_linkers(l_rounds, pop, elite, topM)             # ②

    # ③ 조립: N×M → 결합 재채점 (BBB 연결부위 · 독성 전체 · 전하 가드레일)
    meta, constructs, regions = [], [], []
    for sh in shuttles:
        for lk in linkers:
            meta.append((sh, lk))
            constructs.append(cargo + lk["seq"] + sh["seq"])
            regions.append(bbb_region(cargo, lk["seq"], sh["seq"]))
    bscore = bbb.score(regions)
    tscore = score_tox(constructs, toxpy, toxrepo, toxscript)        # 느린 예측 1회만
    pool = []
    for (sh, lk), c, b, t, r in zip(meta, constructs, bscore, tscore, regions):
        if t > tox_thr:
            continue
        g = charge_guardrail(r)
        pool.append({"shuttle": sh["seq"], "linker": lk["seq"], "sequence": c,
                     "bbb": round(float(b), 4), "tox": round(float(t), 4),
                     "charge": g[1], "muH": g[2], "fit": round(float(b - g[0]), 4),
                     "len": len(c), "linker_len": len(lk["seq"])})
    seen, best_combos = set(), []
    for item in sorted(pool, key=lambda x: x["fit"], reverse=True):
        if item["sequence"] in seen:
            continue
        seen.add(item["sequence"]); best_combos.append(item)

    json.dump({"cargo": cargo, "shuttles": shuttles, "linkers": linkers,
               "best": best_combos[:5], "n_grid": len(constructs)},
              open(out_path, "w"), ensure_ascii=False)
    print("OK")


if __name__ == "__main__":
    a = sys.argv
    main(a[1], int(a[2]), int(a[3]), int(a[4]), int(a[5]), int(a[6]), int(a[7]),
         float(a[8]), a[9], a[10], a[11], a[12])
