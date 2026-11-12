#!/usr/bin/env python3
"""방향성 쌍 검증 — 문헌이 우열을 아는 쌍을 모델이 맞게 순서 매기나.

## 왜 이게 필요한가

지금까지의 모든 숫자는 deepB3P 의 학습 분포 안이거나 우리 자체 척도다. 순위가 **현실의
무엇과도 대조된 적이 없다**. 그래서 "스크리닝 도구"라고 말할 근거가 없다.

가장 싸고 정직한 외부 대조가 **방향성 쌍**이다 — 문헌이 "A가 B보다 뇌 전달이 좋다"고
보고한 쌍들을 모아, 모델이 그 순서를 맞추는지 본다. 절대값 예측이 아니라 **순위 일치율**만
보므로 척도 차이에 강건하고, 결과가 어느 쪽이든 쓸 수 있는 숫자가 나온다.

## 두 종류의 쌍 — 왜 나눠 보나

**(A) 문헌 점 변이 쌍** — 길이가 같고 잔기 1~2개만 다르다. **가장 진단적**이다. 조성도
길이도 거의 같으므로, 모델이 맞히면 그건 조성·길이 아티팩트가 아니라 **서열을 본** 것이다.
문제는 **수가 극히 적다**는 것이다. BBB 펩타이드에서 점 변이 + 뇌 전달 실측을 함께 보고한
연구가 드물다. 아래 2쌍이 출처가 확실한 전부이고, 1쌍은 문헌이 서로 엇갈린다.

**(B) 스크램블 쌍** — 조성 동일·순서 파괴. 수용체 결합 리간드는 순서가 깨지면 결합을 잃는
것이 **메커니즘상 당연**하므로(에피토프가 사라진다), "원본 > 스크램블"을 방향성 쌍으로 쓴다.
각 쌍이 개별로 실측된 건 아니라서 (A)보다 약한 근거지만, **라이브러리 전체에 적용되어
표본이 확보된다**. (A)의 희소성을 (B)가 보완하는 구조다.

(B)에서는 세 기준을 따로 본다:
  - vs 스크램블 **평균**: 느슨한 기준(평균적인 셔플보다 나은가)
  - vs 스크램블 **최대**: 직관적이지만 **표본 수에 의존**한다(셔플을 많이 뽑으면 최대가
    올라간다) → 참고용으로만 본다
  - **경험적 백분위 p** = (원본 이상인 스크램블 수 + 1) / (n + 1) : **n 에 불변**한 정규
    통계량이다. 원본이 진짜 특별하면 p 가 작다. 주 지표로 이걸 쓴다.
"vs 최대"가 그래도 의미 있는 이유 — 실제 설계에서 에이전트는 **최대값을 노린다**. 원본이
최선의 셔플에 밀리면 그 점수로는 "서열이 좋다"고 말할 수 없다. 다만 통계적 주장은 p 로 한다.

## 출처 (2026-11-06 확인)

**A1. Angiopep-2 > Angiopep-7** (K10R·K15R, 19aa 동일 길이, 2치환)
  Angiopep-7 은 뇌 실질 축적이 **약 85% 낮고**, 여러 연구에서 **음성 대조군**으로 쓰인다.
  원 계열 보고: Demeule et al., J Neurochem 2008 (doi:10.1111/j.1471-4159.2008.05492.x);
  Demeule et al., J Pharmacol Exp Ther 2008; 324(3):1064–1072.
  ⚠️ 85% 수치는 **2차 출처**를 통해 확인했다(원문 수치 직접 확인 전).

**A2. Penetratin > Penetratin W48F** (16aa 동일 길이, 1치환)
  Christiaens et al., Eur J Biochem 2004 (doi:10.1111/j.1432-1033.2004.04022.x):
  W48F 는 세포 결합이 야생형의 **67%**, 내재화가 **48%**(4 °C)로 감소.
  ⚠️ Penetratin 은 CPP 다 — 우리 체계에서 **비교군**이며 RMT 가 아니다. 그래도
  "모델이 1잔기 치환의 방향을 아는가"라는 질문에는 유효한 쌍이다.

**A3(논쟁). Penetratin W56F > W48F**
  같은 Christiaens 2004 에서 W56F 는 야생형의 70%, W48F 는 48%(4 °C) → W56F 가 위.
  ⚠️ 그런데 Fischer et al. 은 **Ala** 치환에서 반대 패턴(W56A 가 더 큰 감소)을 보고했고,
  Derossi 1994 와 Thorén 2003 은 이중 변이에서 서로 반대 결론을 냈다. 세포주·온도·
  측정법이 달라서다. **일치율 집계에서 제외**하고 참고로만 출력한다.

**(B) 스크램블 쌍 근거** — 수용체 매개 전달은 특정 에피토프 인식에 의존한다. 순서를
파괴하면 그 에피토프가 사라지므로 결합·전달이 떨어지는 것이 기대값이다. 이 전제가
깨지는 경우(비특이 흡수가 지배하는 양이온성 CPP)는 **결과에서 그렇게 드러날 것**이고,
그것 자체가 정보다.

사용: /usr/bin/python3 scripts/directional_pairs.py [--cargo GSNKGAIIGLM] [--n 20]
"""
from __future__ import annotations

import argparse
import statistics
import sys
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import (OFF_TARGET_PENALTY, SHUTTLES, bbb_scoring_seq,  # noqa: E402
                         get_settings)
from core.order_specificity import scrambles  # noqa: E402
from core.predictors import get_predictor  # noqa: E402

DEFAULT_CARGO = "GSNKGAIIGLM"

ANGIOPEP2 = "TFFYGGSRGKRNNFKTEEY"
ANGIOPEP7 = "TFFYGGSRGRRNNFRTEEY"      # K10R · K15R
PENETRATIN = "RQIKIWFQNRRMKWKK"
PEN_W48F = "RQIKIFFQNRRMKWKK"          # 6번째 잔기(=Antp W48) W→F
PEN_W56F = "RQIKIWFQNRRMKFKK"          # 14번째 잔기(=Antp W56) W→F

# (라벨, 더 좋다고 보고된 쪽, 덜 좋은 쪽, 근거 요약, 집계 포함 여부)
LIT_PAIRS = [
    ("Angiopep-2 > Angiopep-7", ANGIOPEP2, ANGIOPEP7,
     "뇌 실질 축적 ~85% 감소·음성 대조군 (Demeule 2008; 2차 출처 확인)", True),
    ("Penetratin > W48F", PENETRATIN, PEN_W48F,
     "내재화 야생형의 48%, 결합 67% (Christiaens 2004, 4°C)", True),
    ("Penetratin W56F > W48F", PEN_W56F, PEN_W48F,
     "70% vs 48% (Christiaens 2004) — ⚠️ Fischer/Derossi/Thorén 와 엇갈림", False),
]


def binom_p(k: int, n: int) -> float:
    """H0: 동전던지기(p=0.5). k 이상 맞힐 단측 확률."""
    if n == 0:
        return 1.0
    return sum(comb(n, i) for i in range(k, n + 1)) / 2 ** n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cargo", default=DEFAULT_CARGO)
    # 기본값 50 = **문서에 적힌 수치를 재현하는 값**(§14-1·EVALUATION M6). 옛 기본값 20 은
    # 발표 수치를 재현하지 않아 "재현 스크립트가 재현하지 않는" 상태였다.
    ap.add_argument("--n", type=int, default=50, help="셔틀당 스크램블 개수")
    args = ap.parse_args()

    st = get_settings()
    if not st.use_deepb3p_local:
        print("⚠️ 로컬 deepB3P 미감지 — 측정 불가.")
        return
    pr = get_predictor(st)
    C = args.cargo

    def score(seqs, fuse):
        qs = [bbb_scoring_seq(C, "", s) for s in seqs] if fuse else list(seqs)
        return [p.bbb_permeability for p in pr.predict_many(qs)]

    print("=" * 78)
    print("방향성 쌍 검증 — 문헌이 우열을 아는 쌍을 모델이 맞게 순서 매기나")
    print("=" * 78)
    print(f"화물 {C} · 셔틀당 스크램블 {args.n}개")
    print("순위 일치만 본다(절대값 예측이 아니다). H0 = 동전던지기(p=0.5).")

    # ── (A) 문헌 점 변이 쌍 ──────────────────────────────────────────────
    for fuse, lab in ((False, "셔틀 단독"), (True, "융합체")):
        print(f"\n── (A) 문헌 점 변이 쌍 · {lab} ──")
        hi = score([p[1] for p in LIT_PAIRS], fuse)
        lo = score([p[2] for p in LIT_PAIRS], fuse)
        ok = n = 0
        for (lbl, _a, _b, src, counted), x, y in zip(LIT_PAIRS, hi, lo):
            right = x > y
            if counted:
                n += 1
                ok += right
            mark = "✅" if right else "❌"
            tag = "" if counted else "  (집계 제외)"
            print(f"  {mark} {lbl:26s} {x:.3f} vs {y:.3f}  (Δ {x-y:+.3f}){tag}")
            print(f"       근거: {src}")
        print(f"  → 집계 대상 일치 {ok}/{n}" + (f"  (p={binom_p(ok, n):.3f})" if n else ""))

    # ── (B) 스크램블 쌍 ─────────────────────────────────────────────────
    rmt = {k: v["seq"] for k, v in SHUTTLES.items() if "CPP" not in v["note"]}
    cpp = {k: v["seq"] for k, v in SHUTTLES.items() if "CPP" in v["note"]}

    for fuse, lab in ((False, "셔틀 단독"), (True, "융합체")):
        print(f"\n── (B) 스크램블 쌍 (원본 > 스크램블이어야 함) · {lab} ──")
        print(f"  {'셔틀':16s} {'원본':>6s} {'스크평균':>8s} {'스크최대':>8s} "
              f"{'경험p':>7s}  {'vs평균':>6s} {'vs최대':>6s}")
        tallies = {}
        for grp, items in (("후보군 RMT", rmt), ("비교군 CPP", cpp)):
            w_mean = w_max = sig = tot = 0
            rows = []
            for nm, sq in items.items():
                scr = scrambles(sq, args.n, seed=0)
                if not scr:
                    continue
                vals = score([sq] + scr, fuse)
                real, sv = vals[0], vals[1:]
                m, mx = statistics.fmean(sv), max(sv)
                # 경험적 백분위: 원본 이상인 스크램블이 몇 개인가 (+1 보정, n 에 불변)
                pe = (sum(1 for v in sv if v >= real) + 1) / (len(sv) + 1)
                tot += 1
                w_mean += real > m
                w_max += real > mx
                sig += pe <= 0.05
                rows.append((nm, real, m, mx, pe))
            print(f"  [{grp}]")
            for nm, real, m, mx, pe in sorted(rows, key=lambda r: r[4]):
                print(f"  {nm:16s} {real:6.3f} {m:8.3f} {mx:8.3f} {pe:7.3f}  "
                      f"{'✅' if real > m else '❌':>6s} {'✅' if real > mx else '❌':>6s}")
            tallies[grp] = (w_mean, w_max, sig, tot, [r[4] for r in rows])
        for grp, (wm, wx, sig, tot, pes) in tallies.items():
            print(f"  → {grp}: 경험p 중위 **{statistics.median(pes):.3f}** · "
                  f"p≤0.05 인 셔틀 {sig}/{tot}  |  vs평균 {wm}/{tot}"
                  f"(이항 p={binom_p(wm, tot):.3f}) · vs최대 {wx}/{tot}")

    # ── (C) 시스템이 실제로 쓰는 유효점수로도 본다 ────────────────────────
    # eff = BBB − λ·off_target. 문헌 쌍에서 off_target 축이 방향을 되살리는지 확인.
    print("\n── (C) 문헌 쌍 · 시스템 유효점수(eff = BBB − 0.25·off_target) 기준 ──")
    from core.selectivity import assess_selectivity
    hi = score([p[1] for p in LIT_PAIRS], True)
    lo = score([p[2] for p in LIT_PAIRS], True)
    ok = n = 0
    for (lbl, a, b, _src, counted), x, y in zip(LIT_PAIRS, hi, lo):
        ea = x - OFF_TARGET_PENALTY * assess_selectivity(a).off_target_risk
        eb = y - OFF_TARGET_PENALTY * assess_selectivity(b).off_target_risk
        right = ea > eb
        if counted:
            n += 1
            ok += right
        tag = "" if counted else "  (집계 제외)"
        print(f"  {'✅' if right else '❌'} {lbl:26s} eff {ea:+.3f} vs {eb:+.3f}  "
              f"(Δ {ea-eb:+.3f}){tag}")
    print(f"  → 집계 대상 일치 {ok}/{n}" + (f"  (p={binom_p(ok, n):.3f})" if n else ""))

    print("\n" + "=" * 78)
    print("해석 주의: (A)는 표본이 2개다 — 일치든 불일치든 **통계적 결론이 아니다**.")
    print("(B)는 표본이 확보되지만 각 쌍이 개별 실측된 게 아니라 메커니즘상의 기대값이다.")
    print("두 결과가 **같은 방향**을 가리킬 때만 무게를 둘 수 있다.")


if __name__ == "__main__":
    main()
