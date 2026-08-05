# -- coding: utf-8 --
"""링커·셔틀 잠재 co-evolution 루프 (멀티모듈 동시 진화).

FBGAN(_run_fbgan.py)이 '셔틀만' latent z 로 진화시키는 것과 달리, 여기서는
개체 = (z_셔틀 ∈ R^128, 링커 펩타이드) 두 모듈을 **함께** 진화시킨다.

  개체     : (z, linker)  — 셔틀은 잠재 z, 링커는 아미노산 서열(펩타이드 융합 링커)
  조립     : cargo + linker + shuttle(=G.decode(z))
  적합도   : 비독성이면 BBB(deepB3P), 독성이면 페널티. 미세 길이 절약항으로 컴팩트 링커 선호.
  선택     : 결합 적합도 상위 elite 개체(z·linker 쌍째로 보존)
  재생산   : child_z = eliteA.z + 변이,  child_linker = crossover(eliteA, eliteB) 후 변이
             → 좋은 셔틀 잠재(A)와 좋은 링커 계통(A×B)이 재조합되며 **공진화**
  이민     : 소량은 새 랜덤 z + 라이브러리 시드 링커 (탐색 유지)

BBB·독성 채점기는 FBGAN 러너와 동일 구현을 재사용한다(중복 제거).

Usage:
  python _run_coevo.py <cargo> <rounds> <pop> <elite> <tox_thr> \
      <toxpy> <toxrepo> <toxscript> <out.json>
"""
import json
import sys

import numpy as np

import torch                          # noqa: E402  (torch.load 패치는 _run_fbgan import 시 적용됨)
from _run_fbgan import (               # 동일 venv·동일 스코어러 재사용
    HIDDEN, N_CHARS, SEQ_LEN,
    BBBScorer, bbb_region, decode, score_tox,
)
from fbgan.models import Generator

from _physchem import charge_guardrail   # 양전하 편향 억제 물성 가드레일

# ---------------- 링커 유전자(펩타이드 융합 링커) ----------------
# 알파벳: 유연(G,S,A,P)·전하/강직(E,K)·방향족(Y,W,F)·소수성(V) — 링커에서 흔한 잔기군
LINK_ALPHABET = "GSAEKPYWFV"
LINK_MIN, LINK_MAX = 3, 18            # 라이브러리 최장(A(EAAAK)3A=17)에 맞춘 상한
LINK_SEEDS = [                        # 초기·이민 개체 시드: 검증된 링커 모티프
    "GGS", "GGGGS", "GGGGGS", "GGGGSGGGGS", "GGGGSGGGGSGGGGS",
    "YGGGGS", "WGGGGS", "FGGGGS", "VGGGGS",
    "EAAAK", "EAAAKEAAAK", "AEAAAKEAAAKEAAAKA",
]
LEN_PARSIMONY = 0.003                 # 미세 길이 절약(동점 시 컴팩트 링커 선호), BBB 신호는 지배 못함


def _clip_len(lk):
    if len(lk) > LINK_MAX:
        lk = lk[:LINK_MAX]
    while len(lk) < LINK_MIN:
        lk += "G"
    return lk


def mutate_linker(lk, rng):
    """점변이(치환) / 삽입 / 결실 중 하나."""
    s = list(lk)
    op = rng.rand()
    if op < 0.6 and s:                                   # 치환
        s[rng.randint(len(s))] = LINK_ALPHABET[rng.randint(len(LINK_ALPHABET))]
    elif op < 0.8 and len(s) < LINK_MAX:                 # 삽입
        s.insert(rng.randint(len(s) + 1), LINK_ALPHABET[rng.randint(len(LINK_ALPHABET))])
    elif len(s) > LINK_MIN:                              # 결실
        del s[rng.randint(len(s))]
    return _clip_len("".join(s))


def cross_linker(a, b, rng):
    """단일점 교차 — 두 부모 링커 계통을 재조합."""
    if not a or not b:
        return _clip_len(a or b or "GGGGS")
    ca = rng.randint(1, len(a)) if len(a) > 1 else 1
    cb = rng.randint(1, len(b)) if len(b) > 1 else 1
    return _clip_len(a[:ca] + b[cb:])


def main(cargo, rounds, pop, elite, tox_thr, toxpy, toxrepo, toxscript, out_path):
    G = Generator(n_chars=N_CHARS, seq_len=SEQ_LEN, bs=pop, hidden=HIDDEN)
    G.load_state_dict(torch.load("fbgan/checkpoint/G_weights_1000.pth"))
    bbb = BBBScorer()

    rng = np.random.RandomState(2022)
    Z = rng.randn(pop, 128)                                        # 셔틀 잠재 모듈
    L = [LINK_SEEDS[i % len(LINK_SEEDS)] for i in range(pop)]       # 링커 모듈(시드)
    sigma = 0.6
    history, best = [], {}

    for rd in range(rounds):
        shuttles = decode(G, Z)
        constructs = [cargo + lk + sh for lk, sh in zip(L, shuttles)]
        regions = [bbb_region(cargo, lk, sh) for lk, sh in zip(L, shuttles)]
        bbb_scores = bbb.score(regions)
        tox_scores = score_tox(constructs, toxpy, toxrepo, toxscript)
        guards = [charge_guardrail(r) for r in regions]   # (penalty, net_charge, muH)

        fitness = np.array([
            (b if t <= tox_thr else b - 1.0) - g[0] - LEN_PARSIMONY * len(lk)
            for b, t, g, lk in zip(bbb_scores, tox_scores, guards, L)
        ])
        order = np.argsort(fitness)[::-1]

        # 전역 베스트 (비독성 중 보정 적합도 최고)
        for i in order:
            if tox_scores[i] <= tox_thr:
                cand = {"shuttle": shuttles[i], "linker": L[i], "sequence": constructs[i],
                        "bbb": round(float(bbb_scores[i]), 4), "tox": round(float(tox_scores[i]), 4),
                        "charge": guards[i][1], "muH": guards[i][2],
                        "fit": round(float(fitness[i]), 4),
                        "len": len(constructs[i]), "linker_len": len(L[i])}
                if not best or cand["fit"] > best["fit"]:
                    best = cand
                break

        safe = [b for b, t in zip(bbb_scores, tox_scores) if t <= tox_thr]
        history.append({"round": rd + 1,
                        "mean_bbb": round(float(np.mean(bbb_scores)), 4),
                        "best_bbb": round(float(np.max(bbb_scores)), 4),
                        "n_safe": len(safe),
                        "n_uniq_linker": len(set(L)),
                        "best_linker": L[int(order[0])],
                        "mean_charge": round(float(np.mean([g[1] for g in guards])), 1)})

        # --- 공진화 재생산: 셔틀 잠재 변이 + 링커 교차·변이 ---
        e_idx = order[:elite]
        childZ, childL = [], []
        n_children = pop - elite - 2
        while len(childZ) < n_children:
            pa = int(e_idx[rng.randint(elite)])
            pb = int(e_idx[rng.randint(elite)])
            childZ.append(Z[pa] + sigma * rng.randn(128))              # 셔틀 잠재(A) 변이
            childL.append(mutate_linker(cross_linker(L[pa], L[pb], rng), rng))  # 링커 A×B 재조합
        # elite 보존 + 자식 + 이민 2(새 z + 시드 링커)
        Z = np.vstack([Z[e_idx], np.array(childZ), rng.randn(2, 128)])
        L = ([L[i] for i in e_idx] + childL +
             [LINK_SEEDS[rng.randint(len(LINK_SEEDS))] for _ in range(2)])
        sigma *= 0.85

    # 마지막 세대 수집
    final = decode(G, Z)
    fc = [cargo + lk + sh for lk, sh in zip(L, final)]
    fregions = [bbb_region(cargo, lk, sh) for lk, sh in zip(L, final)]
    fb = bbb.score(fregions)
    ft = score_tox(fc, toxpy, toxrepo, toxscript)
    fg = [charge_guardrail(r) for r in fregions]
    pool = [{"shuttle": sh, "linker": lk, "sequence": c,
             "bbb": round(float(b), 4), "tox": round(float(t), 4),
             "charge": g[1], "muH": g[2], "fit": round(float(b - g[0]), 4),
             "len": len(c), "linker_len": len(lk)}
            for sh, lk, c, b, t, g in zip(final, L, fc, fb, ft, fg) if t <= tox_thr]
    if best:
        pool.append(best)
    seen, top = set(), []
    for item in sorted(pool, key=lambda x: x["fit"], reverse=True):
        if item["sequence"] in seen:
            continue
        seen.add(item["sequence"]); top.append(item)

    json.dump({"history": history, "best": top[:5], "cargo": cargo},
              open(out_path, "w"), ensure_ascii=False)
    print("OK")


if __name__ == "__main__":
    a = sys.argv
    main(a[1], int(a[2]), int(a[3]), int(a[4]), float(a[5]), a[6], a[7], a[8], a[9])
