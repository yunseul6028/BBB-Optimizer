"""응집 축 외부 레이블 심사 — 문헌 아밀로이드 vs 가용성 펩타이드 (docs/experiments.md §3 재기초).

왜 있나
  응집 축은 최적화기가 **능동적으로 최소화하는 3축 중 하나**인데, 이 축이 옳은 방향을
  가리킨다는 외부 레이블 증거는 §3 의 `AUC(가용>응집) = 0.96` 하나뿐이다.
  🔴 **그런데 그 0.96 을 만든 평가셋 20종이 저장소 어디에도 기록돼 있지 않다.**
  즉 이 스크립트는 **0.96 의 재현이 아니다** — 평가셋을 코드에 박아 넣어 지금부터
  재현 가능하게 만든 **새 측정**이다. 두 숫자는 서로 **교환 불가**하며, 여기 나오는
  숫자를 "0.96 을 재현했다"고 쓰면 안 된다. 값이 다르면 그것은 오류가 아니라 **발견**이다.
  ⚠️ 숫자를 맞추려고 평가셋을 손보는 것은 금지(그 순간 외부 레이블 심사가 아니게 된다).

무엇을 재나
  - 양성(낮아야 함) = 문헌 아밀로이드·응집성 펩타이드 10종 (steric-zipper 코어 중심)
  - 음성(높아야 함) = 가용성·거동 양호 펩타이드 10종 (이 프로젝트의 링커 + 표준 태그)
  - 지표 3개를 같은 셋에서 비교한다:
      · `1 − agg_hotspot`  … 창 5aa 슬라이딩 **국소최대** (현재 순위·Pareto 축)
      · `1 − agg_score`    … **전역 평균** (M10 에서 길이 희석 때문에 순위에서 퇴역)
      · `sol_score`        … 용해도 상대지표(보조 지표, 순위 축 아님)
  📌 왜 부호를 뒤집나: 응집 점수는 **높을수록 응집**이라 방향이 §3 의 `AUC(가용>응집)`
     과 반대다. `1 − agg` 로 "높을수록 가용"으로 맞춰야 §3 표와 같은 방향으로 읽힌다
     (원점수 기준 AUC = 1 − 여기 값). 선형 단조 변환이라 AUC 자체는 부호 뒤집기만 한다.
  📌 왜 `sol_score` 도 같이 재나: §3 의 한계로 기록된 값(Sup35 GNNQQNY 0.74, polyQ 0.77)은
     `agg_*` 가 아니라 **`sol_score` 의 값과 일치**한다. 즉 §3 의 0.96 은 `sol_score`
     계열 양에서 계산된 것으로 보이며, 비교 가능하게 두려면 이 열이 필요하다.
     그러나 `sol_score` 는 **내부에 `agg_score`(전역 평균)를 품고** 있어 순위 축이 아니다.

실행: 프로젝트 루트에서 Biopython 이 있는 파이썬으로 (⚠️ `python3 -I` 금지 — core 임포트 깨짐)
  $ /usr/bin/python3 scripts/bench_aggregation.py
완전 결정론(무작위·하위프로세스 없음) — 시드 불필요.

이 AUC 가 **증명하는 것**: 외부 문헌 레이블(아밀로이드 vs 가용성)에 대해 응집 축의 순위
  방향이 맞다는 것, 그리고 핫스팟 전환이 판별력을 얼마나 바꿨는지.
이 AUC 가 **증명하지 못하는 것**: 대응 널(matched null)이 아니다 — 같은 크기의 무작위
  편집과 비교한 것이 아니므로 "이 축으로 **개선**할 수 있다"는 근거가 전혀 아니다.
  또한 평가셋은 저자(나)가 고른 20종이고 구성체 수준이 아니라 **짧은 펩타이드 수준**이다.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.developability import assess_developability
from scripts.bench_selectivity import auc          # 동일 Mann-Whitney 구현 재사용(복제 금지)

# ── 평가셋 (🔴 여기가 기록이다 — 바꾸면 그 측정은 다른 측정이다) ────────────────
# 양성: 문헌 아밀로이드·응집성 펩타이드 — 응집 점수가 **높아야**(= 1−agg 가 낮아야) 맞다.
AMYLOID = {
    "Aβ(16-22) KLVFFAE":  "KLVFFAE",                # Aβ 중심 소수성 클러스터(CHC), 고전 steric zipper
    "Aβ(33-42)":          "GLMVGGVVIA",             # Aβ42 C-말단 소수성 단편(섬유 코어)
    "Tau PHF6":           "VQIVYK",                 # tau 306-311 육펩타이드, 응집 핵
    "IAPP(22-27)":        "NFGAIL",                 # 인간 아밀린(IAPP) 아밀로이드 코어
    "Sup35 GNNQQNY":      "GNNQQNY",                # 효모 Sup35 극성 steric zipper(§3 의 기지 실패 사례)
    "polyQ Q10":          "QQQQQQQQQQ",             # 폴리글루타민 반복(§3 의 기지 실패 사례)
    "α-syn NAC(68-78)":   "GAVVTGVTAVA",            # α-시누클레인 NAC 영역 응집 코어
    "PrP(106-126)":       "KTNMKHMAGAAAAGAVVGGLG",  # 프리온 단백질 응집성 단편 ⚠️ 서열 미대조
    "Insulin B(11-17)":   "LVEALYL",                # 인슐린 B사슬 아밀로이드 코어(steric zipper 결정구조)
    "STVIIE(설계)":        "STVIIE",                 # de novo 설계 아밀로이드 육펩타이드(Serrano 계열)
}

# 음성: 가용성·거동 양호 펩타이드 — 이 프로젝트 링커 라이브러리(core/config.py)의
#       유연·강직 링커가 핵심이다. **이들의 가용성이 링커 라이브러리 전체의 전제**이므로
#       축이 이들을 응집성으로 매기면 그건 축의 문제다. 나머지는 표준 정제·에피토프 태그.
SOLUBLE = {
    "(GGGGS)2":           "GGGGSGGGGS",             # 프로젝트 유연 링커(10aa) — 표준 GS 링커
    "(GGGGS)3":           "GGGGSGGGGSGGGGS",        # 프로젝트 유연 링커(15aa)
    "GS 반복(10aa)":       "GSGSGSGSGS",             # 일반 GS 반복(무구조·친수)
    "(EAAAK)2":           "EAAAKEAAAK",             # 프로젝트 강직 α-헬릭스 링커(10aa)
    "A(EAAAK)3A":         "AEAAAKEAAAKEAAAKA",      # 프로젝트 최장 강직 링커(17aa)
    "His6 태그":           "HHHHHH",                 # 표준 정제 태그(가용)
    "FLAG 태그":           "DYKDDDDK",               # 표준 에피토프 태그(강한 음전하·가용)
    "Myc 태그":            "EQKLISEEDL",             # 표준 에피토프 태그(가용)
    "Strep-tag II":       "WSHPQFEK",               # 표준 친화 태그(가용)
    "poly-E(E10)":        "EEEEEEEEEE",             # 음전하 용해도 증강 태그
}

# 미스 판정 임계 — **판별 임계가 아니라 보고용**이다(AUC 는 순위 기반이라 임계 무관).
# 두 군 전체 중앙값을 기준으로 삼아 자유 파라미터를 만들지 않는다.
_METRICS = ("1−agg_hotspot", "1−agg_score", "sol_score")


def _metrics(seq: str) -> dict:
    r = assess_developability(seq)
    if r.error:
        raise SystemExit(f"개발성 계산 실패: {seq} — {r.error}")
    return {
        "1−agg_hotspot": round(1.0 - r.agg_hotspot, 3),
        "1−agg_score":   round(1.0 - r.agg_score, 3),
        "sol_score":     r.sol_score,
        "hot_raw":       r.agg_hotspot_a3v,
        "mean_raw":      r.agg_a3v,
    }


def main() -> None:
    print("응집 축 외부 레이블 심사 — 문헌 아밀로이드 vs 가용성 펩타이드")
    print("🔴 docs/experiments.md §3 의 AUC 0.96 은 평가셋이 기록되지 않아 **재현 불가**다.")
    print("   아래는 평가셋을 코드에 고정한 **새 측정**이며 0.96 과 교환 불가하다.\n")
    print(f"아밀로이드 {len(AMYLOID)}종 / 가용성 {len(SOLUBLE)}종 (결정론, 무작위 없음)\n")

    amy = {n: _metrics(s) for n, s in AMYLOID.items()}
    sol = {n: _metrics(s) for n, s in SOLUBLE.items()}

    hdr = (f"{'펩타이드':<20}{'군':>6}{'길이':>5}{'hot_raw':>9}{'1−hot':>8}"
           f"{'mean_raw':>10}{'1−mean':>9}{'sol':>7}")
    print(hdr)
    for group, table, seqs in (("응집", amy, AMYLOID), ("가용", sol, SOLUBLE)):
        for n, m in table.items():
            print(f"{n:<20}{group:>6}{len(seqs[n]):>5}{m['hot_raw']:>9.3f}"
                  f"{m['1−agg_hotspot']:>8.3f}{m['mean_raw']:>10.3f}"
                  f"{m['1−agg_score']:>9.3f}{m['sol_score']:>7.3f}")

    print("\n=== 요약: AUC(가용 > 응집) ===")
    print("(지표를 1−agg 로 뒤집어 '높을수록 가용'으로 맞춤 — 원점수 AUC = 1 − 아래 값)")
    print(f"{'지표':<16}{'가용평균':>10}{'응집평균':>10}{'분리':>8}{'AUC':>8}   비고")
    notes = {
        "1−agg_hotspot": "현재 순위·Pareto 축(창 5aa 국소최대)",
        "1−agg_score":   "퇴역 전역평균(M10: 길이 희석)",
        "sol_score":     "보조 지표 — §3 의 0.96 이 계산된 양으로 추정",
    }
    summary = {}
    for k in _METRICS:
        pos = [m[k] for m in sol.values()]      # 가용 = 높아야 함
        neg = [m[k] for m in amy.values()]      # 응집 = 낮아야 함
        mp, mn = sum(pos) / len(pos), sum(neg) / len(neg)
        a = auc(pos, neg)
        summary[k] = (mp, mn, a)
        print(f"{k:<16}{mp:>10.3f}{mn:>10.3f}{mp - mn:>8.3f}{a:>8.3f}   {notes[k]}")

    d = summary["1−agg_hotspot"][2] - summary["1−agg_score"][2]
    verdict = ("핫스팟 전환이 판별력을 **높였다**" if d > 0.005 else
               "핫스팟 전환이 판별력을 **낮췄다**" if d < -0.005 else
               "핫스팟 전환은 판별력을 **바꾸지 않았다**")
    print(f"\n핫스팟(순위 축) − 전역평균(퇴역) AUC 차 = {d:+.3f} → {verdict}.")
    print("⚠️ 핫스팟을 쓰는 이유는 판별력이 아니라 **길이 희석 면역**(M10)이다 — "
          "이 숫자가 낮게 나와도 전역평균으로 되돌리는 근거가 되지 않는다.")

    # ── 미스 목록 — 군 전체 중앙값 기준(자유 파라미터 아님) ──────────────────
    print("\n=== 미스 목록 (순위 축 `1−agg_hotspot` 기준) ===")
    key = "1−agg_hotspot"
    allv = sorted([m[key] for m in amy.values()] + [m[key] for m in sol.values()])
    mid = (allv[len(allv) // 2 - 1] + allv[len(allv) // 2]) / 2
    print(f"기준선 = 20종 전체 중앙값 {mid:.3f} (판정 임계가 아니라 **가시화용**)")

    amy_miss = [(n, m[key]) for n, m in amy.items() if m[key] > mid]      # 응집인데 가용으로 보임
    sol_miss = [(n, m[key]) for n, m in sol.items() if m[key] < mid]      # 가용인데 응집으로 보임
    print(f"\n[응집인데 가용처럼 높게 나온 것] {len(amy_miss)}/{len(amy)}")
    for n, v in sorted(amy_miss, key=lambda x: -x[1]):
        print(f"  - {n:<20}1−hot {v:.3f}  (sol {amy[n]['sol_score']:.3f})")
    print(f"\n[가용인데 응집처럼 낮게 나온 것] {len(sol_miss)}/{len(sol)}")
    for n, v in sorted(sol_miss, key=lambda x: x[1]):
        print(f"  - {n:<20}1−hot {v:.3f}  (sol {sol[n]['sol_score']:.3f})")
    if not sol_miss:
        print("  - 없음 (링커 라이브러리 전원 기준선 위 = 링커 가용성 전제 유지)")

    # ── §3 의 기록된 한계(극성 아밀로이드)가 아직도 재현되나 ──────────────────
    print("\n=== §3 기록 한계 재검: 극성 아밀로이드(steric zipper) ===")
    print("§3 기록: Sup35 GNNQQNY 0.74 · polyQ 0.77 (= 가용으로 오판) — 조성 기반의 원천 한계")
    polar = ("Sup35 GNNQQNY", "polyQ Q10")
    for n in polar:
        m = amy[n]
        print(f"  {n:<20}1−hot {m['1−agg_hotspot']:.3f} · 1−mean {m['1−agg_score']:.3f} "
              f"· sol {m['sol_score']:.3f} → {'오판(가용쪽)' if m[key] > mid else '정상 검출'}")
    still = all(amy[n][key] > mid for n in polar)
    print(f"  판정: 극성 아밀로이드 한계는 {'**여전히 재현된다**' if still else '재현되지 않는다'}"
          f" — 조성 기반 척도로는 원천 불가(Waltz/TANGO 류 구조 기반 필요).")
    if (amy["Sup35 GNNQQNY"]["sol_score"], amy["polyQ Q10"]["sol_score"]) == (0.74, 0.77):
        print("  📌 §3 에 적힌 0.74·0.77 은 `sol_score` 값과 **소수 둘째자리까지 일치**한다 →"
              " §3 의 0.96 은 `sol_score` 계열에서 계산된 것으로 보인다(순위 축 아님).")


if __name__ == "__main__":
    main()
