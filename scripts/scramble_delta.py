#!/usr/bin/env python3
"""스크램블 델타 측정 — deepB3P 점수가 '순서'를 보나 '조성'을 보나.

왜 이게 필요한가
----------------
M1 실측: 양성 vs 스크램블 AUC = 0.618 (양성 vs 무작위는 0.807). 즉 조성만 같고
**순서를 깨도 점수가 꽤 유지된다** → deepB3P 점수의 상당 부분이 아미노산 조성이다.
그래서 에이전트가 아르기닌을 쌓으면 점수가 오른다. 보상 해킹의 메커니즘이 이것이다.

대응: 후보마다 **자기 스크램블을 대조군으로** 같이 채점하고
    delta = bbb(서열) - 평균 bbb(스크램블 n개)
를 요구한다. 조성 쌓기로는 delta 가 0 에 붙어서 못 통과한다.

이 스크립트는 **게이트로 쓸 수 있는지 먼저 재는** 교정(calibration) 도구다.
라이브러리의 정당한 RMT 셔틀들이 줄줄이 탈락하면 게이트가 아니라 **지표**로만 써야 한다.
그 판단에 필요한 분포를 뽑는 것이 목적이고, 결론을 내리지 않는다.

측정 두 가지
-----------
A) 셔틀 단독       — 셔틀 서열 자체의 순서 특이성
B) 융합체(직접융합) — 에이전트가 실제로 최적화하는 대상(cargo + 셔틀)

스크램블은 **셔틀 잔기만** 섞는다. 융합체 전체를 섞으면 화물까지 깨져서
"셔틀의 순서가 일을 하나"라는 질문에 답하지 못한다.

사용: /usr/bin/python3 scripts/scramble_delta.py [--n 20] [--cargo GSNKGAIIGLM]
"""
from __future__ import annotations

import argparse
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import SHUTTLES, bbb_scoring_seq, get_settings  # noqa: E402
from core.order_specificity import DELTA_MIN, scrambles  # noqa: E402
from core.predictors import get_predictor  # noqa: E402

DEFAULT_CARGO = "GSNKGAIIGLM"

# 에이전트가 실제로 만들어낸/만들 수 있는 **조성 쌓기 의심 서열**을 같이 재서
# "정당한 셔틀"과 분리되는지 본다. 이게 분리되지 않으면 게이트는 쓸모가 없다.
PROBES = {
    "에이전트설계 ApoE변이": "LRVYQAGRRLRVYQAGRR",   # 공진화 세션이 올린 후보(R 치환)
    "조성쌓기 R/K 혼합": "RKRKRKRKRKRK",
    "조성쌓기 Arg-rich": "RRRRRRRRRG",       # 반복서열이지만 섞을 수 있다(G 하나)
    "조성쌓기 Arg+방향족": "RRRRWRRRRWRR",
    # 아래 둘은 **동일잔기 반복**이라 유효 스크램블이 0개 → 델타 측정 불가(n/a).
    # 일부러 남겨둔다: 측정 불가가 n/a 로 찍히는지 확인하는 자기검사 역할.
    "조성쌓기 poly-Arg(측정불가)": "RRRRRRRRR",
    "음성대조 poly-Ser(측정불가)": "SSSSSSSSSSSS",
}


def measure(predictor, items: dict[str, str], cargo: str, n: int, fuse: bool):
    """items(이름→셔틀서열) 각각에 대해 (원본점수, 스크램블 평균/표준편차, delta)."""
    def prep(sh: str) -> str:
        return bbb_scoring_seq(cargo, "", sh) if fuse else sh

    # 모든 서열을 **한 번에** 예측해 서브프로세스 왕복을 줄인다.
    names = list(items)
    scr = {nm: scrambles(items[nm], n) for nm in names}
    seqs = [prep(items[nm]) for nm in names]
    for nm in names:
        seqs.extend(prep(s) for s in scr[nm])
    preds = [p.bbb_permeability for p in predictor.predict_many(seqs)]

    rows, cur = [], len(names)
    for i, nm in enumerate(names):
        real = preds[i]
        sv = preds[cur:cur + len(scr[nm])]
        cur += len(scr[nm])
        if not sv:
            # 유효 스크램블 0개 = 측정 불가. 0 으로 때우지 않는다(위 docstring 참고).
            rows.append({"name": nm, "seq": items[nm], "real": real, "n_scr": 0,
                         "scr_mean": None, "scr_sd": None, "scr_max": None, "delta": None})
            continue
        mean = statistics.fmean(sv)
        rows.append({"name": nm, "seq": items[nm], "real": real, "n_scr": len(sv),
                     "scr_mean": mean,
                     "scr_sd": statistics.pstdev(sv) if len(sv) > 1 else 0.0,
                     "scr_max": max(sv), "delta": real - mean})
    return rows


def _f(v, w=6, p=3):
    return f"{'    n/a':>{w}s}" if v is None else f"{v:{w}.{p}f}"


def show(title: str, rows: list[dict], note: str = "") -> None:
    print(f"\n── {title} ──")
    if note:
        print(f"  {note}")
    print(f"  {'후보':22s} {'원본':>6s} {'스크램블평균':>10s} {'±sd':>6s} "
          f"{'스크램블최대':>10s} {'델타':>7s}  n")
    # 측정 불가(None)는 맨 아래로 내린다 — 정렬에 섞여 순위를 왜곡하지 않게.
    for r in sorted(rows, key=lambda x: (x["delta"] is None, -(x["delta"] or 0))):
        d = "    n/a" if r["delta"] is None else f"{r['delta']:+7.3f}"
        print(f"  {r['name']:22s} {r['real']:6.3f} {_f(r['scr_mean'], 10)} "
              f"{_f(r['scr_sd'])} {_f(r['scr_max'], 10)} {d}  {r['n_scr']}")
    nd = [r["name"] for r in rows if r["delta"] is None]
    if nd:
        print(f"  ℹ️ 측정 불가(유효 스크램블 0개 — 동일잔기 반복서열): {', '.join(nd)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20, help="후보당 스크램블 개수")
    ap.add_argument("--cargo", default=DEFAULT_CARGO)
    args = ap.parse_args()

    settings = get_settings()
    if not settings.use_deepb3p_local:
        print("⚠️ 로컬 deepB3P 미감지 — 측정 불가.")
        return
    predictor = get_predictor(settings)

    rmt = {k: v["seq"] for k, v in SHUTTLES.items() if "CPP" not in v["note"]}
    cpp = {k: v["seq"] for k, v in SHUTTLES.items() if "CPP" in v["note"]}

    print("=" * 74)
    print("스크램블 델타 — deepB3P 점수의 '순서 특이성' 측정")
    print("=" * 74)
    print(f"화물 {args.cargo} · 후보당 스크램블 {args.n}개 (조성 보존·순서 파괴)")
    print("델타 = 원본 BBB − 스크램블 평균 BBB.  델타≈0 이면 그 점수는 **조성**에서 온 것이다.")
    print(f"설계체 하드 게이트 마진 = {DELTA_MIN:+.3f} (core/order_specificity.py)")

    for fuse, label in ((False, "A) 셔틀 단독"), (True, f"B) 융합체 ({args.cargo} + 직접융합 + 셔틀)")):
        r_rmt = measure(predictor, rmt, args.cargo, args.n, fuse)
        r_cpp = measure(predictor, cpp, args.cargo, args.n, fuse)
        r_pb = measure(predictor, PROBES, args.cargo, args.n, fuse)
        show(f"{label} · 후보군 RMT 셔틀 {len(rmt)}종", r_rmt)
        show(f"{label} · 비교군 CPP {len(cpp)}종", r_cpp)
        show(f"{label} · 프로브(조성 쌓기 의심)", r_pb,
             "정당한 셔틀과 **분리돼야** 게이트로 쓸 수 있다.")

        def summ(rows):
            d = [r["delta"] for r in rows if r["delta"] is not None]
            return (f"중위 {statistics.median(d):+.3f} (범위 {min(d):+.3f}~{max(d):+.3f}, "
                    f"n={len(d)})") if d else "측정 가능 항목 없음"

        print(f"\n  요약({label}):")
        print(f"    후보군 RMT  델타 {summ(r_rmt)}")
        print(f"    비교군 CPP  델타 {summ(r_cpp)}")
        print(f"    프로브      델타 {summ(r_pb)}")
        neg = [r["name"] for r in r_rmt if r["delta"] is not None and r["delta"] <= 0]
        if neg:
            print(f"    ⚠️ 델타≤0 인 **정당한** RMT 셔틀 {len(neg)}/{len(r_rmt)}: "
                  f"{', '.join(neg)}")
            print("       → 이만큼 오탈락하면 게이트로 못 쓴다. 지표로만 써야 한다.")


if __name__ == "__main__":
    main()
