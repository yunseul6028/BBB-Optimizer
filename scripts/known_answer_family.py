#!/usr/bin/env /usr/bin/python3
"""정답지 확장 — **Angiopep 계열 13종**을 우리 축으로 줄 세워 문헌 순서와 맞대본다 (§56).

왜 이 스크립트가 있나
---------------------
M6(§14)의 가장 진단적인 부분인 **문헌 점 변이 쌍**이 **n=2** 였다. "통계가 아니다"가
`docs/protocol.md` §8 에 공백으로 적혀 있고, 병목은 코드가 아니라 **문헌**이었다.

그 병목이 풀렸다 — 우리가 이미 인용하던 **그 논문**(Demeule et al. 2008)의 본문 표에
**같은 계열 19-mer 11종의 정량 측정값이 통째로 들어 있었다.** 2차 출처만 보고
Angiopep-2/5/7 세 개만 쓰고 있었다.

정답지 (전부 **원문 표에서 직접 전사**, 2026-11-15 PDF 대조)
  · 표 2 — BBCEC 단층 **in vitro 투과율(%)**, 1시간, 필터 ≥3개, 11종
  · 표 3 — **in situ 뇌 관류** 5분 후 분포용적 `ml/100 g/5 min`,
           총뇌(B)·모세혈관·**실질(P)** 과 **P/B 비율**, 8종
  · 본문 — Angiopep-2(= 67 의 C7S)는 실질 흡수가 Angiopep-1 과 **유사**,
           Angiopep-5(K10R) **−60%**, Angiopep-7(K10R·K15R) **−85%**, 셋 다 순전하 +2

🔴 **이 정답지의 한계(쓰기 전에 읽을 것)**
  ⓐ **표 3 은 마우스 2마리**다(원문: "means obtained after the perfusion of two mice").
     Fig. 3(Angiopep-2/5/7)은 "at least three mice". **n 이 작다.**
  ⓑ **한 논문·한 실험실·한 계열**이다. 독립 재현이 아니다.
  ⓒ 표 2(in vitro)와 표 3(in situ)은 **다른 수준**이고, 아래에서 둘이 **엇갈리는 쌍**이 있다.
  ⓓ 🔴 **deepB3P 의 홀드아웃이 아니다.** Angiopep 계열은 2008년 발표라 학습 분포에
     포함됐을 수 있다. 진짜 외부 홀드아웃은 **학습 시점 이후** 서열이어야 한다(§14-6).
  ⓔ 측정 대상은 **셔틀 단독**(화물 없음)이다. 구성체 수준은 참고로만 같이 낸다.
  ⓕ 🔴 **산포가 크다**(2026-11-15 추가). 원문은 `평균 ± S.D.` 이고 표 2 의 S.D. 는
     평균의 10~30% 다. 특히 **67(15.8 ± 5.1)과 91(15.8 ± 3.3)은 완전 동률**이다.
     `SD_T2`·`SD_T3_P` 와 `tie_report()` 가 구분 불가능한 쌍을 센다. 순위상관을
     인용할 때 이 비율을 **같은 문단에** 적을 것.

실행: /usr/bin/python3 scripts/known_answer_family.py   (유료 LLM 0회 · 결정론)
출처: Demeule M, Régina A, Ché C, Poirier J, Nguyen T, Gabathuler R, Castaigne JP,
      Béliveau R. J Pharmacol Exp Ther. 2008;324(3):1064-1072. DOI 10.1124/jpet.107.131318
      표 2(p.1067) · 표 3(p.1067) · 본문 p.1069 · Fig. 3.
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from core.binding import shuttle_similarity
from core.config import DEFAULT_CARGO, bbb_scoring_seq, get_settings
from core.developability import assess_developability
from core.predictors import get_predictor
from core.selectivity import assess_selectivity
from core.stability import assess_stability

# ── 정답지 (원문 표 전사) ─────────────────────────────────────────────────
# name: (서열, 표2 투과율 %, 표3 실질 Vd, 표3 총뇌 Vd, 표3 모세혈관 Vd, 표3 P/B)
# None = 그 표에 없음
FAMILY: dict[str, tuple] = {
    "5":          ("TFFYGGCRAKRNNFKRAKY",  5.1,  95.0, 312.0, 217.0, 0.30),
    "8":          ("TFFYGGCRGKKNNFKRAKY",  3.9,  46.0, 250.0, 204.0, 0.18),
    "67/Ang-1":   ("TFFYGGCRGKRNNFKTEEY", 15.8,  25.0,  38.0,  13.0, 0.66),
    "75":         ("PFFYGGCRGKRNNFKTEEY",  8.7,  None,  None,  None, None),
    "76":         ("TFFYGGCRGKRNNFKTKEY",  8.4,  24.0,  40.0,  16.0, 0.60),
    "77":         ("TFFYGGKRGKRNNFKTKEY",  6.3,  None,  None,  None, None),
    "78":         ("TFFYGGCRGKRNNFKTKRY",  7.5,  16.0, 198.0, 181.0, 0.08),
    "79":         ("TFFYGGKRGKRNNFKTAEY",  7.8,  18.0,  70.0,  52.0, 0.26),
    "81":         ("TFFYGGKRGKRNNFKREKY",  8.8,  None,  None,  None, None),
    "90":         ("RFKYGGCLGNKNNFLRLKY", 13.4,  11.0,  87.0,  76.0, 0.13),
    "91":         ("RFKYGGCLGNKNNYLRLKY", 15.8,  23.0,  47.0,  24.0, 0.49),
}

# 🔴 **원문의 표준편차**(2026-11-15 추가, §57-3).
# 처음 전사할 때 평균만 가져왔는데, 원문은 `평균 ± S.D.` 로 적혀 있고 **산포가 크다**.
# 이걸 빼놓으면 "순위"가 측정 분해능보다 촘촘해 보인다 — 아래 `tie_report()` 가 센다.
SD_T2: dict[str, float] = {        # 표 2 투과율 % 의 S.D. (필터 ≥3개)
    "5": 0.8, "8": 0.3, "67/Ang-1": 5.1, "75": 0.7, "76": 2.7, "77": 2.3,
    "78": 1.0, "79": 0.7, "81": 0.4, "90": 1.2, "91": 3.3,
}
SD_T3_P: dict[str, float] = {      # 표 3 **실질** Vd 의 S.D. (마우스 2마리)
    "5": 58.0, "8": 26.0, "67/Ang-1": 8.0, "76": 7.0,
    "78": 1.0, "79": 5.0, "90": 7.0, "91": 16.0,
}


def tie_report() -> None:
    """측정 분해능으로 **구분되지 않는 쌍**을 센다.

    기준은 보수적으로 `|a-b| <= sd_a + sd_b` (두 S.D. 의 합). 이 안에 들어오면
    "a 가 b 보다 낫다"를 이 표만 보고 말할 수 없다.
    """
    for label, vals, sds in (("표2 투과율", {n: FAMILY[n][1] for n in FAMILY}, SD_T2),
                             ("표3 실질 Vd", {n: FAMILY[n][2] for n in FAMILY
                                            if FAMILY[n][2] is not None}, SD_T3_P)):
        names = [n for n in vals if n in sds]
        pairs = list(combinations(names, 2))
        tied = [(a, b) for a, b in pairs
                if abs(vals[a] - vals[b]) <= sds[a] + sds[b]]
        print(f"  {label:<12} 쌍 {len(pairs):3d} 중 **구분 안 되는 쌍 {len(tied):3d}**"
              f" ({100 * len(tied) / len(pairs):.0f}%)")
        exact = [(a, b) for a, b in pairs if vals[a] == vals[b]]
        for a, b in exact:
            print(f"      🔴 완전 동률: {a} = {b} = {vals[a]}")
# 본문 서술값(Fig. 3) — Angiopep-1 의 실질 흡수를 100 으로 둔 **상대값**이다.
# 🔴 절대 Vd 가 아니므로 표 3 과 **같은 열에 넣지 않는다**.
RELATIVE = {
    "Angiopep-2": ("TFFYGGSRGKRNNFKTEEY", 100.0, "Angiopep-1 과 유사"),
    "Angiopep-5": ("TFFYGGSRGRRNNFKTEEY",  40.0, "−60%"),
    "Angiopep-7": ("TFFYGGSRGRRNNFRTEEY",  15.0, "−85%"),
}
NET_CHARGE_NOTE = {"5": "+6 (원문 명시)", "8": "+6 (원문 명시)",
                   "Angiopep-2": "+2", "Angiopep-5": "+2", "Angiopep-7": "+2"}


N_PERM = 200_000
PERM_SEED = 20261115


def perm_test(xs, ys) -> tuple[float, float]:
    """스피어만 ρ 의 **양측 순열검정**. n=11 이면 전수(39.9M)가 느려 몬테카를로를 쓴다."""
    import random
    obs = spearman(xs, ys)
    rng = random.Random(PERM_SEED)
    sh = list(ys)
    hit = 0
    for _ in range(N_PERM):
        rng.shuffle(sh)
        if abs(spearman(xs, sh)) >= abs(obs) - 1e-12:
            hit += 1
    return obs, (hit + 1) / (N_PERM + 1)


N_BOOT = 10_000
BOOT_SEED = 20261116


def sd_bootstrap(names, vals, sds, axis_vals):
    """🔴 **문헌값의 산포를 ρ 로 전달한다** (§57-3).

    원문이 주는 것은 점이 아니라 `평균 ± S.D.` 다. 각 펩타이드의 문헌값을
    `N(평균, S.D.)` 에서 다시 뽑아 ρ 를 `N_BOOT` 번 계산한다. 축 값은 결정론이라
    고정이고 **문헌 쪽만** 흔든다. 보고할 것은 세 개다 —
      · ρ 의 중위와 2.5~97.5% 구간
      · **부호가 유지되는 비율**(이게 '방향'이 산포를 견디는가에 대한 답)
      · |ρ| 가 관측치만큼 큰 비율
    """
    import random
    rng = random.Random(BOOT_SEED)
    obs = spearman(axis_vals, vals)
    rhos = []
    for _ in range(N_BOOT):
        drawn = [rng.gauss(vals[i], sds[names[i]]) for i in range(len(names))]
        rhos.append(spearman(axis_vals, drawn))
    rhos.sort()
    lo, hi = rhos[int(0.025 * N_BOOT)], rhos[int(0.975 * N_BOOT) - 1]
    same = sum(1 for r in rhos if (r < 0) == (obs < 0)) / N_BOOT
    return obs, median(rhos), lo, hi, same


def spearman(xs, ys) -> float:
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else 0.0


def axes_for(seq: str, pred) -> dict:
    dev = assess_developability(seq)
    sel = assess_selectivity(seq)
    construct = bbb_scoring_seq(DEFAULT_CARGO, "", seq)
    return {
        "deepB3P(셔틀)": pred.predict(seq).bbb_permeability,
        "deepB3P(구성체)": pred.predict(construct).bbb_permeability,
        "rmt_sim": shuttle_similarity(seq).score or 0.0,
        "off_target": sel.off_target_risk,
        "순전하": dev.net_charge,
        "불안정성": assess_stability(seq).instability_index,
        "응집 핫스팟": dev.agg_hotspot,
        "liability 개수": float(dev.n_liabilities),
    }


def main() -> None:
    pred = get_predictor(get_settings())
    names = list(FAMILY)
    ax = {n: axes_for(FAMILY[n][0], pred) for n in names}
    for n, (s, *_rest) in RELATIVE.items():
        ax[n] = axes_for(s, pred)
    axis_names = list(next(iter(ax.values())))

    print("=" * 100)
    print("§56 정답지 확장 — Angiopep 계열 (Demeule 2008 표 2·표 3 원문 전사)")
    print("=" * 100)

    print("\n[1] 정답지")
    print(f"{'펩타이드':<12}{'서열':<22}{'투과율%':>8}{'실질 Vd':>9}{'총뇌 Vd':>9}{'P/B':>7}  순전하")
    for n in names:
        s, t2, p, b, _c, pb = FAMILY[n]
        f = lambda v, w, d=1: (f"{v:>{w}.{d}f}" if v is not None else f"{'—':>{w}}")
        print(f"{n:<12}{s:<22}{f(t2,8)}{f(p,9)}{f(b,9)}{f(pb,7,2)}  "
              f"{NET_CHARGE_NOTE.get(n, '')}")
    print("  (상대값 — 표 3 과 단위가 다르므로 상관 계산에서 제외)")
    for n, (s, rel, note) in RELATIVE.items():
        print(f"{n:<12}{s:<22}{'—':>8}{'Ang-1=100 기준 ' + str(rel):>18}  {note}")

    print("\n[2] 우리 축 값")
    print(f"{'펩타이드':<12}" + "".join(f"{a:>15}" for a in axis_names))
    for n in names + list(RELATIVE):
        print(f"{n:<12}" + "".join(f"{ax[n][a]:>15.3f}" for a in axis_names))

    print("\n[3] 문헌 지표 ↔ 우리 축 순위상관 (스피어만)")
    targets = {
        "표2 투과율(n=11)":  [(n, FAMILY[n][1]) for n in names],
        "표3 실질 Vd(n=8)":  [(n, FAMILY[n][2]) for n in names if FAMILY[n][2] is not None],
        "표3 P/B 비율(n=8)": [(n, FAMILY[n][5]) for n in names if FAMILY[n][5] is not None],
    }
    print(f"{'축':<16}" + "".join(f"{k:>20}" for k in targets))
    for a in axis_names:
        row = ""
        for _k, pairs in targets.items():
            xs = [ax[n][a] for n, _ in pairs]
            ys = [v for _, v in pairs]
            row += f"{spearman(xs, ys):>+20.3f}"
        print(f"{a:<16}{row}")
    print("  (양수 = 축이 클수록 문헌값도 큼. 'off_target·불안정성·응집·liability' 는")
    print("   작을수록 좋다고 보는 축이므로 **음수가 나와야** 방향이 맞는 것이다.)")

    print(f"\n[3b] 순열검정 (양측 · 몬테카를로 {N_PERM:,}회 · 시드 {PERM_SEED})")
    print(f"{'축':<16}" + "".join(f"{k:>24}" for k in targets))
    for a in axis_names:
        row = ""
        for _k, prs in targets.items():
            rho, pv = perm_test([ax[n][a] for n, _ in prs], [v for _, v in prs])
            row += f"{f'{rho:+.3f} (p={pv:.3f})':>24}"
        print(f"{a:<16}{row}")
    print(f"  🔴 다중비교: {len(axis_names)}축 × {len(targets)}표 = {len(axis_names)*len(targets)}회 검정이다.")
    print(f"     본페로니 기준은 0.05/{len(axis_names)*len(targets)} = {0.05/(len(axis_names)*len(targets)):.4f} 이고,")
    print("     아래 유의 결과는 그 기준을 넘지 못한다. 다만 deepB3P 역전은 §14(M6)가")
    print("     **먼저 예측한 방향**이라 탐색적 검정이 아니라 **확인**이다.")

    print("\n[4] 단일 치환 쌍 — 가장 진단적 (M6 의 n=2 를 대체)")
    base = {n: FAMILY[n][0] for n in names} | {n: RELATIVE[n][0] for n in RELATIVE}
    pairs = []
    for a, b in combinations(base, 2):
        sa, sb = base[a], base[b]
        if len(sa) != len(sb):
            continue
        d = [(i + 1, sa[i], sb[i]) for i in range(len(sa)) if sa[i] != sb[i]]
        if len(d) == 1:
            pairs.append((a, b, d[0]))
    print(f"  단일 치환 쌍 **{len(pairs)}개** (M6 은 2개였다)")
    for metric_i, label in ((1, "표2 투과율"), (2, "표3 실질 Vd")):
        usable = [(a, b, d) for a, b, d in pairs
                  if a in FAMILY and b in FAMILY
                  and FAMILY[a][metric_i] is not None and FAMILY[b][metric_i] is not None]
        if not usable:
            continue
        sdmap = SD_T2 if metric_i == 1 else SD_T3_P
        print(f"\n  ── 정답 = {label} ({len(usable)}쌍) ──")
        print("    🔴 '해상' 칸: 원문 S.D. 로 두 값이 **구분되는가**(|a-b| > sd_a+sd_b).")
        print("       '노이즈'인 쌍은 O/X 자체가 의미 없다 — 정답지가 순서를 정하지 못한다.")
        print(f"    {'쌍':<26}{'치환':>8}{'해상':>8}{'문헌':>18}"
              + "".join(f"{a:>15}" for a in axis_names))
        agree = {a: 0 for a in axis_names}
        agree_res = {a: 0 for a in axis_names}
        n_res = 0
        for a, b, (p, x, y) in usable:
            va, vb = FAMILY[a][metric_i], FAMILY[b][metric_i]
            hi, lo = (a, b) if va > vb else (b, a)
            resolved = abs(va - vb) > sdmap.get(a, 0.0) + sdmap.get(b, 0.0)
            n_res += 1 if resolved else 0
            cells = ""
            for an in axis_names:
                ok = ax[hi][an] > ax[lo][an]
                lower_better = an in ("off_target", "불안정성", "응집 핫스팟",
                                      "liability 개수", "순전하")
                correct = (not ok) if lower_better else ok
                if ax[hi][an] == ax[lo][an]:
                    cells += f"{'=':>15}"
                else:
                    agree[an] += 1 if correct else 0
                    agree_res[an] += 1 if (correct and resolved) else 0
                    cells += f"{'O' if correct else 'X':>15}"
            print(f"    {a + ' vs ' + b:<26}{x+str(p)+y:>8}"
                  f"{'해상' if resolved else '노이즈':>8}{hi + ' 높음':>18}{cells}")
        n_cmp = len(usable)
        print(f"    {'일치 합계(전체)':<26}{'':>8}{'':>8}{'':>18}"
              + "".join(f"{str(agree[a]) + '/' + str(n_cmp):>15}" for a in axis_names))
        print(f"    {'🔴 일치(해상된 쌍만)':<24}{'':>8}{'':>8}{'':>18}"
              + "".join(f"{str(agree_res[a]) + '/' + str(n_res):>15}" for a in axis_names))
        if n_res <= 1:
            print(f"    🔴 해상된 쌍이 {n_res}개뿐이다 — 이 표로 축을 줄 세우면 **안 된다**.")

    print("\n[5] 🔴 문헌 산포를 ρ 로 전달 — 부트스트랩 (§57-3)")
    print(f"    문헌값을 N(평균, S.D.) 에서 {N_BOOT:,}번 다시 뽑아 ρ 를 재계산한다.")
    print(f"    {'축':<18}{'관측 ρ':>9}{'중위':>9}{'95% 구간':>20}{'부호 유지':>10}")
    bn = [n for n in names if n in SD_T2]
    bvals = [FAMILY[n][1] for n in bn]
    bsd = SD_T2
    for an in axis_names:
        av = [ax[n][an] for n in bn]
        o, med, lo, hi, same = sd_bootstrap(bn, bvals, bsd, av)
        print(f"    {an:<18}{o:>9.3f}{med:>9.3f}"
              f"{f'[{lo:+.3f}, {hi:+.3f}]':>20}{100 * same:>9.1f}%")
    print("    → '부호 유지'가 95% 미만이면 **방향조차 산포 안에서 뒤집힌다**는 뜻이다.")

    print("\n[6] 🔴 측정 분해능 — 원문 S.D. 로 **구분되지 않는 쌍** 세기 (§57-3)")
    tie_report()
    print("  → 이 비율이 높으면 '우리 축이 문헌 순서를 맞혔다/틀렸다'의 **분모 자체가 무릅니다**.")
    print("    상관계수는 그대로 두되, 해석에 반드시 같이 적을 것.")

    print("\n🔴 한계는 파일 상단 독스트링 ⓐ~ⓕ 를 그대로 읽을 것. 특히 표 3 은 **마우스 2마리**다.")


if __name__ == "__main__":
    main()
