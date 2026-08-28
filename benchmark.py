"""에이전트 성능 벤치마크 — 평가 세트(evaluation set) + 평가 지표(metrics).

심사 기준 "에이전트 성능 평가를 위한 적절한 evaluation set과 평가 지표를 설정했는가"에 대응.

LLM 없이 로컬 엔진(deepB3P·ToxinPred3)만으로 도는 핵심 지표 — 지금 바로 재현 가능:
  M1) 판별력(Discrimination): 알려진 BBB 투과 펩타이드(양성) vs 스크램블·무작위(음성)를
      deepB3P가 분리하는가 (AUC + 평균 분리). → 점수 엔진이 신호를 잡음.
  M2) 탐색 공간 & 최적성: 라이브러리 전수 탐색(oracle)의 비독성 최적 vs 무작위 탐색 효율 곡선.
      → "얼마나 탐색해야 최적에 근접하나"의 기준선. 에이전트는 이보다 적은 평가로 도달해야 가치.
  M3) 결정론적 가드레일: 독성 펩타이드(멜리틴 등)를 파이프라인이 올바로 탈락시키는가.

LLM이 필요한 지표(옵션, 여기선 미실행 — 쿼터 리셋/크레딧 후):
  M4) 에이전트 효율/최적성 갭: 에이전트가 M2 최적 대비 몇 번의 평가로 도달/초과하나.
  M5) 심사(Critic) 효과: 일부러 나쁜 후보를 주면 Critic이 REVISE로 걸러내나.

실행:  python benchmark.py           # M1~M3 (로컬, 무료)
"""

from __future__ import annotations

import random

from core.config import (
    DEFAULT_CARGO,
    LINKER_LIBRARY,
    SHUTTLES,
    TOXICITY_THRESHOLD,
    bbb_scoring_seq,
    get_settings,
)
from core.predictors import get_predictor
from core.toxicity import get_toxicity_predictor

SEED = 42

# ─────────────────────────────────────────────────────────────────────────────
# 평가 세트 (Evaluation set)
# ─────────────────────────────────────────────────────────────────────────────
# 양성: 문헌에 보고된 BBB 투과/세포투과(CPP)·셔틀 펩타이드
POSITIVES = {
    "TAT": "GRKKRRQRRRPPQ",
    "Penetratin": "RQIKIWFQNRRMKWKK",
    "Angiopep-2": "TFFYGGSRGKRNNFKTEEY",
    "Angiopep-1": "TFFYGGCRGKRNNFKTEEY",
    "SynB1": "RGGRLSYSRRRFSTSTGR",
    "Transportan": "GWTLNSAGYLLGKINLKALAALAKKIL",
    "pVEC": "LLIILRRRIRKQAHAHSK",
    "R8": "RRRRRRRR",
    "MAP": "KLALKLALKALKAALKLA",
    "pep-1": "KETWWETWWTEWSQPKKKRKV",
    "RVG29": "YTIWMPENPRPGTPCDIFTNSRGKRASNG",
    "TP10": "AGYLLGKINLKALAALAKKIL",
}

# 음성(무작위): CPP가 아닌 무작위/친수성 서열 — 투과 신호가 없어야 함
NEGATIVES_RANDOM = {
    "poly-Ser": "SSSSSSSSSSSS",
    "poly-Glu": "EEEEEEEEEE",
    "rand-hydrophil": "SGNTDSQGTNSDGQ",
    "rand-acidic": "DEDEEGSDEGSDE",
    "rand-mix1": "AGSTNQDTGSANMT",
    "rand-mix2": "TGSNADQMTSGANT",
    "albumin-frag": "DAHKSEVAHRFKDLGE",
    "insulinB-frag": "FVNQHLCGSHLVEAL",
}

# 독성 가드레일용: 알려진 독성 펩타이드(양성) vs 양성대조 비독성
TOXIC_KNOWN = {
    "Melittin": "GIGAVLKVLTTGLPALISWIKRKRQQ",   # 벌독 — 강한 세포용해 독성
    "Mastoparan": "INLKALAALAKKIL",              # 말벌독
}
BENIGN_KNOWN = {
    "GS-linker": "GGGGSGGGGS",
    "Angiopep-2": "TFFYGGSRGKRNNFKTEEY",
}


def _scramble(seq: str, rng: random.Random) -> str:
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def _auc(pos, neg) -> float:
    """순위 기반 AUC = 무작위 양성이 무작위 음성보다 높을 확률(Mann-Whitney)."""
    if not pos or not neg:
        return 0.5
    c = t = 0
    for p in pos:
        for n in neg:
            t += 1
            c += 1 if p > n else (0.5 if p == n else 0)
    return c / t


def _bbb_scores(predictor, seqs):
    return [r.bbb_permeability for r in predictor.predict_many(seqs)]


# ─────────────────────────────────────────────────────────────────────────────
# M1) 판별력
# ─────────────────────────────────────────────────────────────────────────────
def m1_discrimination(predictor):
    rng = random.Random(SEED)
    pos_names, pos_seqs = list(POSITIVES), list(POSITIVES.values())
    scr_seqs = [_scramble(s, rng) for s in pos_seqs]          # 조성 동일·순서 파괴
    rnd_seqs = list(NEGATIVES_RANDOM.values())

    pos = _bbb_scores(predictor, pos_seqs)
    scr = _bbb_scores(predictor, scr_seqs)
    rnd = _bbb_scores(predictor, rnd_seqs)

    print("── M1. 판별력 (deepB3P 점수, 0~1) ──")
    print(f"  양성(CPP/셔틀 {len(pos)}개)   평균 {sum(pos)/len(pos):.3f}")
    print(f"  음성-무작위({len(rnd)}개)      평균 {sum(rnd)/len(rnd):.3f}")
    print(f"  음성-스크램블({len(scr)}개)     평균 {sum(scr)/len(scr):.3f} (조성 동일, 순서 파괴)")
    print(f"  AUC 양성 vs 무작위:   {_auc(pos, rnd):.3f}  (1.0=완벽 분리, 0.5=무작위)")
    print(f"  AUC 양성 vs 스크램블: {_auc(pos, scr):.3f}  (순서 신호 포착도)")
    worst = sorted(zip(pos_names, pos), key=lambda x: x[1])[:3]
    print(f"  양성 중 저점(주의): {', '.join(f'{n}={s:.2f}' for n, s in worst)}")


# ─────────────────────────────────────────────────────────────────────────────
# M2) 탐색 공간 & 최적성 + 무작위 탐색 효율
# ─────────────────────────────────────────────────────────────────────────────
def _library_scores(cargo, predictor, tox_pred, thr):
    combos = [(ln, sn) for ln in LINKER_LIBRARY for sn in SHUTTLES]
    bbb_seqs, full_seqs = [], []
    for ln, sn in combos:
        lk, sh = LINKER_LIBRARY[ln]["seq"], SHUTTLES[sn]["seq"]
        bbb_seqs.append(bbb_scoring_seq(cargo, lk, sh))
        full_seqs.append(cargo + lk + sh)
    bbb = _bbb_scores(predictor, bbb_seqs)
    tox = [r.risk for r in tox_pred.predict_many(full_seqs)]
    rows = []
    for (ln, sn), b, t in zip(combos, bbb, tox):
        rows.append({"linker": ln, "shuttle": sn, "bbb": b, "tox": t, "toxic": t > thr})
    return rows


def m2_search(cargo, predictor, tox_pred, thr):
    rows = _library_scores(cargo, predictor, tox_pred, thr)
    valid = [r for r in rows if not r["toxic"]]
    opt = max(valid, key=lambda r: r["bbb"]) if valid else None
    print("\n── M2. 탐색 공간 & 최적성 ──")
    print(f"  화물: {cargo}  |  탐색 공간: 링커 {len(LINKER_LIBRARY)} × 셔틀 {len(SHUTTLES)} = {len(rows)}개 조합")
    print(f"  비독성 후보: {len(valid)}/{len(rows)}")
    if opt is None:
        print("  (비독성 최적 없음)")
        return
    print(f"  🏆 라이브러리 최적(oracle): {opt['linker']}+{opt['shuttle']}  "
          f"BBB={opt['bbb']:.3f}  독성={opt['tox']:.3f}")

    # 무작위 탐색 효율: K개 평가 시 최적의 몇 %에 도달 (여러 시드 평균)
    print("  무작위 탐색 효율(최적 대비 도달률, 시드 20개 평균):")
    for k in (3, 6, 12, 24):
        reach = []
        for seed in range(20):
            rng = random.Random(seed)
            sample = rng.sample(rows, min(k, len(rows)))
            v = [r["bbb"] for r in sample if not r["toxic"]]
            reach.append((max(v) / opt["bbb"]) if v else 0.0)
        print(f"    {k:2d}개 평가 → 최적의 {100*sum(reach)/len(reach):5.1f}%")
    print("  → 에이전트(M4)는 이보다 적은 평가로 최적 근접/초과해야 가치가 증명됨.")


# ─────────────────────────────────────────────────────────────────────────────
# M3) 결정론적 가드레일 (독성 탈락)
# ─────────────────────────────────────────────────────────────────────────────
def m3_guardrails(tox_pred, thr):
    print("\n── M3. 결정론적 가드레일 (독성 탈락, 임계값 {:.2f}) ──".format(thr))
    tox_seqs = list(TOXIC_KNOWN.values())
    ben_seqs = list(BENIGN_KNOWN.values())
    tox_r = [r.risk for r in tox_pred.predict_many(tox_seqs)]
    ben_r = [r.risk for r in tox_pred.predict_many(ben_seqs)]
    hits = 0
    for name, r in zip(TOXIC_KNOWN, tox_r):
        flag = r > thr
        hits += flag
        print(f"  독성 알려짐 {name:12s} risk={r:.2f} → {'✅ 탈락' if flag else '❌ 통과(오탐)'}")
    for name, r in zip(BENIGN_KNOWN, ben_r):
        print(f"  비독성 대조 {name:12s} risk={r:.2f} → {'❌ 탈락(오탐)' if r>thr else '✅ 통과'}")
    print(f"  독성 검출율: {hits}/{len(TOXIC_KNOWN)}")


def main():
    settings = get_settings()
    if not settings.use_deepb3p_local or not settings.use_toxinpred3_local:
        print("⚠️ 로컬 엔진(deepB3P/ToxinPred3) 미감지 — 벤치마크 불가.")
        return
    predictor = get_predictor(settings)
    tox_pred = get_toxicity_predictor(settings)
    thr = TOXICITY_THRESHOLD
    print("=" * 68)
    print("BBB-Optimizer 벤치마크 — 평가 세트 + 지표 (로컬 엔진, LLM 미사용)")
    print("=" * 68)
    m1_discrimination(predictor)
    m2_search(DEFAULT_CARGO, predictor, tox_pred, thr)
    m3_guardrails(tox_pred, thr)
    print("=" * 68)
    print("M4(에이전트 효율)·M5(심사 효과)는 LLM 필요 — benchmark_agent.py 참고(쿼터 리셋 후).")


if __name__ == "__main__":
    main()
