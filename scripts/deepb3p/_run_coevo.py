# -- coding: utf-8 --
"""링커·셔틀 co-evolution 루프 (멀티모듈 동시 진화 — 서열 directed evolution).

개체 = (셔틀 서열, 링커 서열) 두 모듈을 **함께** 진화시킨다. 핵심 원칙:
  · 셔틀은 **de-novo 생성 금지** — 검증된 라이브러리 리간드(Angiopep·T7·ApoE·RVG29…)에서
    **시드**해 보존적 point-mutation/crossover 로 진화한다(잠재공간 GAN 생성 폐기, 가짜 셔틀 차단).
  · 링커는 검증 링커 모티프에서 시드해 진화한다.

  개체     : (shuttle, linker)  — 둘 다 아미노산 서열
  조립     : cargo + linker + shuttle
  적합도   : 비독성이면 BBB(deepB3P), 독성이면 페널티. 전하 가드레일·미세 길이 절약항.
  선택     : 결합 적합도 상위 elite 개체(shuttle·linker 쌍째 보존)
  재생산   : child_shuttle = mutate(cross(A.sh, B.sh)),  child_linker = mutate(cross(A.lk, B.lk))
             → 좋은 셔틀 계통과 좋은 링커 계통이 재조합되며 **공진화**
  이민     : 소량은 **라이브러리 시드** 셔틀 + 라이브러리 시드 링커 (탐색 유지, de-novo 아님)

BBB·독성 채점기는 _run_fbgan 구현을 재사용한다(셔틀 생성기 G/latent 는 더는 쓰지 않음).

Usage:
  python _run_coevo.py <cargo> <rounds> <pop> <elite> <tox_thr> \
      <toxpy> <toxrepo> <toxscript> <out.json> <shuttle_seeds_csv>
"""
import json
import sys

import numpy as np

from _run_fbgan import BBBScorer, bbb_region, score_tox   # 동일 venv·동일 스코어러 재사용
from _physchem import charge_guardrail                     # 양전하 편향 억제 물성 가드레일

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


# ---------------- 셔틀 유전자(검증 리간드 서열의 directed evolution) ----------------
# 셔틀은 de-novo 생성 금지 — 검증된 RMT/CPP 리간드 서열에서 시드해 **보존적으로** 진화한다.
# (수용체 결합 근거 없는 무작위 서열이 deepB3P 편향을 게이밍하는 것을 원천 차단)
SH_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"   # 표준 20 아미노산
SH_MIN, SH_MAX = 5, 40                  # 라이브러리 셔틀 길이(7~30) 포괄, deepB3P 50aa 이내


def _clip_sh(sh):
    if len(sh) > SH_MAX:
        sh = sh[:SH_MAX]
    while len(sh) < SH_MIN:
        sh += "G"
    return sh


def mutate_shuttle(sh, rng):
    """점변이(치환) 위주 보존적 진화 / 삽입 / 결실 — 시드 리간드에서 크게 벗어나지 않게."""
    s = list(sh)
    op = rng.rand()
    if op < 0.7 and s:                                   # 치환(보존적)
        s[rng.randint(len(s))] = SH_ALPHABET[rng.randint(len(SH_ALPHABET))]
    elif op < 0.85 and len(s) < SH_MAX:                  # 삽입
        s.insert(rng.randint(len(s) + 1), SH_ALPHABET[rng.randint(len(SH_ALPHABET))])
    elif len(s) > SH_MIN:                                # 결실
        del s[rng.randint(len(s))]
    return _clip_sh("".join(s))


def cross_shuttle(a, b, rng):
    """단일점 교차 — 두 부모 셔틀 계통을 재조합."""
    if not a or not b:
        return _clip_sh(a or b or "GGGGS")
    ca = rng.randint(1, len(a)) if len(a) > 1 else 1
    cb = rng.randint(1, len(b)) if len(b) > 1 else 1
    return _clip_sh(a[:ca] + b[cb:])


def main(cargo, rounds, pop, elite, tox_thr, toxpy, toxrepo, toxscript, out_path, sh_seeds_csv):
    bbb = BBBScorer()
    sh_seeds = [_clip_sh(s) for s in sh_seeds_csv.split(",") if s] or ["GGGGS"]

    rng = np.random.RandomState(2022)
    S = [sh_seeds[i % len(sh_seeds)] for i in range(pop)]           # 셔틀 모듈(라이브러리 시드)
    L = [LINK_SEEDS[i % len(LINK_SEEDS)] for i in range(pop)]       # 링커 모듈(시드)
    history, best = [], {}

    for rd in range(rounds):
        shuttles = [_clip_sh(s) for s in S]
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

        # --- 공진화 재생산: 셔틀·링커 각각 교차·변이 ---
        e_idx = order[:elite]
        childS, childL = [], []
        n_children = pop - elite - 2
        while len(childS) < n_children:
            pa = int(e_idx[rng.randint(elite)])
            pb = int(e_idx[rng.randint(elite)])
            childS.append(mutate_shuttle(cross_shuttle(S[pa], S[pb], rng), rng))   # 셔틀 A×B
            childL.append(mutate_linker(cross_linker(L[pa], L[pb], rng), rng))     # 링커 A×B
        # elite 보존 + 자식 + 이민 2(라이브러리 시드 셔틀·링커, de-novo 아님)
        S = ([S[i] for i in e_idx] + childS +
             [sh_seeds[rng.randint(len(sh_seeds))] for _ in range(2)])
        L = ([L[i] for i in e_idx] + childL +
             [LINK_SEEDS[rng.randint(len(LINK_SEEDS))] for _ in range(2)])

    # 마지막 세대 수집
    final = [_clip_sh(s) for s in S]
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
    main(a[1], int(a[2]), int(a[3]), int(a[4]), float(a[5]), a[6], a[7], a[8], a[9], a[10])
