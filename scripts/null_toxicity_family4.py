#!/usr/bin/env python3
"""§44 공백 ⓖ — 독성 게이트의 **계열 4**(검출 + 한 수 해킹) + 외부 라벨 첫 데이터점.

왜 돌리나
--------
§44 가 계열 1·2 를 닫고 **계열 3·4 를 공백으로 남겼다**(`docs/protocol.md` §6 의 독성
게이트 행이 "🔴 계열 3·4 공백"으로 적는다). 그리고 §44-8 이 확인한 구조적 사실 —
**독성은 이 파이프라인의 유일한 하드 게이트다**(`core/optimizer_agent.py:366-368`) —
때문에 이 공백이 현재 가장 비싸다. 다른 축은 틀리면 순위가 흔들리는데, 이 축이 틀리면
**통과/탈락 자체가 틀린다.**

🔴 **먼저 어느 양인지 못 박는다**: 재는 것은 **ToxinPred3 ML Score** 이고 **독성이 아니다**
(§44 와 같은 규율). 아래 "게이트가 잡는다/놓친다"는 전부 *"프로덕션 게이트가 쓰는 그
숫자가 0.38 을 넘나"* 의 줄임말이다.

세 블록
-------
**A. 외부 라벨 프로브** — 교과서적 독성 펩타이드 vs 임상·내인성 무해 펩타이드에
   게이트를 걸어 본다. 🔴 **이것은 계열 3 이 아니다.** 계열 3(라벨 분해 애블레이션)은
   모델 내부의 라벨 조회 항을 떼는 검사이고, ToxinPred3 는 블랙박스라 그게 불가능하다.
   이 블록은 **외부 라벨의 첫 데이터점**일 뿐이고, 두 가지 이유로 **상한값**이다:
     ⑴ 여기 쓰는 독소들은 ToxinPred3 학습 데이터(독소 DB)에 **거의 확실히 들어 있다** →
        양성 쪽은 암기를 재는 것이고 일반화가 아니다.
     ⑵ n 이 10 단위다.
   🔴 **그래서 이 검사는 비대칭이다 — 실패하면 치명적이고 통과하면 거의 무정보다.**
   단 **음성 쪽은 사정이 다르다**: 우리 설계체는 독소가 아니라 **치료용 융합 펩타이드**다.
   무해 펩타이드가 0.38 을 넘으면 그 위양성률은 **우리 후보가 속한 클래스에서의
   실측**이고, 학습셋 오염으로 설명되지 않는다. 그쪽이 이 블록의 진짜 수확이다.

**B. 계열 4 — 한 수 해킹.** 게이트에 **탈락한** 구성체에서 **단 한 수**로 0.38 아래로
   내려오는 수를 전수 열거한다. 세 종류를 **따로** 센다:
     ① 링커 안 단일 치환  ② 링커 라이브러리 1회 교체  ③ 셔틀 안 단일 치환
   그리고 각 수가 **다른 제약에 걸리는지** 확인한다. 🔴 예상되는 비대칭을 미리 적는다 —
   이 시스템의 제약(`rmt_sim ≥ 0.85` · off 양측 밴드 · 문헌 보존 위치)은 **전부 셔틀에만**
   걸린다. 링커를 건드리는 수는 셔틀을 바꾸지 않으므로 **그 제약 전부를 통과한다**(동점).
   즉 ①②가 게이트를 넘기면 그건 **제약이 못 잡는 경로**다.

**C. 계열 4 — 검출 짝(악화 검출).** 게이트를 **통과한** 구성체의 링커 자리에 **독소 조각**을
   같은 길이로 심는다. 게이트가 잡나. 대조는 §44 계열 1 과 같은 **같은 길이 무작위**다.
   잡으면 "악화 검출 신뢰", 놓치면 🔴 **"유일한 하드 게이트가 독소 삽입을 통과시킨다"** 다.

결정론·위생
----------
외부 호출은 ToxinPred3 서브프로세스뿐이고 `core/toxicity.py` 의 flock 으로 **직렬화**된다
(§37-1 경쟁 조건). 무작위는 전부 `random.Random(SEED)` 파생. 몽키패치 **0건** — 생산
함수(`get_toxicity_predictor`)를 직접 호출한다. 끝에 **같은 배치를 두 번 호출해 비트 단위
동일**을 단언한다(§44 와 같은 검사).

이 측정이 **볼 수 없는 것**
--------------------------
ⓐ 블록 A 의 양성 쪽은 **학습셋 오염 가능성이 높다**(위 ⑴). ⓑ 블록 A 의 펩타이드는
   **손으로 고른 것**이고 체계적 벤치마크가 아니다. ⓒ 블록 B 는 **단일 수**만 본다 —
   두 수 이상의 조합은 안 본다. ⓓ 블록 C 의 독소 조각은 **절단된 조각**이므로 실제
   독성을 갖는지는 별개 문제다. 여기서 묻는 것은 *"게이트가 알아보나"* 뿐이다.
   ⓔ 임계 **0.38 자체의 근거는 여전히 범위 밖**이다(ToxinPred3 기본값 수용, §44-9ⓑ).
   ⓕ 실제 세포·동물 독성은 **전혀 재지 않는다**.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.binding import shuttle_similarity  # noqa: E402
from core.config import LINKER_LIBRARY, SHUTTLES, get_settings  # noqa: E402
from core.evidence import CONSERVED_POSITIONS  # noqa: E402
from core.selectivity import assess_selectivity  # noqa: E402
from core.toxicity import get_toxicity_predictor  # noqa: E402

SEED = 0
AA = "ACDEFGHIKLMNPQRSTVWY"
THRESH = 0.38          # core.config.TOXICITY_THRESHOLD
RMT_FLOOR = 0.85
OFF_BAND = 0.05
CHUNK = 900            # ToxinPred3 배치 크기(§44 와 같은 규모)

CARGOS = {"Aβ25-35": "GSNKGAIIGLM", "DSIP": "WAGGDASGE"}

# ── 블록 A 외부 라벨 ────────────────────────────────────────────────────────
# 🔴 손으로 고른 교과서 펩타이드다. 출처는 모두 표준 약리학·독소학 문헌이고,
#    여기서는 **라벨**만 쓴다(서열의 1차 출처 대조는 이 측정의 범위가 아니다).
TOXIC = {
    "Melittin(꿀벌)":        "GIGAVLKVLTTGLPALISWIKRKRQQ",
    "Mastoparan(말벌)":      "INLKALAALAKKIL",
    "Pardaxin(가자미)":      "GFFALIPKIISSPLFKTLLSAVGSALSSSGGQE",
    "Apamin(꿀벌)":          "CNCKAPETALCARRCQQH",
    "Conotoxin-GI(청자고둥)": "ECCNPACGRHYSC",
    "Bombolitin-I(호박벌)":  "IKITTMLAKLGKVLAHV",
}
BENIGN = {
    "Glucagon":        "HSQGTFTSDYSKYLDSRRAQDFVQWLMNT",
    "GLP-1(7-37)":     "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR",
    "Oxytocin":        "CYIQNCPLG",
    "Vasopressin":     "CYFQNCPRG",
    "Leu-enkephalin":  "YGGFL",
    "Somatostatin-14": "AGCKNFFWKTFTSC",
    "Substance P":     "RPKPQQFFGLM",
    "Angiotensin II":  "DRVYIHPF",
    "Bradykinin":      "RPPGFSPFR",
    "DSIP":            "WAGGDASGE",
}


def score_all(pred, seqs: list[str]) -> list[float]:
    """청크로 나눠 점수화. 🔴 병렬화 금지(공유 cwd, §37-1)."""
    out: list[float] = []
    for i in range(0, len(seqs), CHUNK):
        out.extend(r.risk for r in pred.predict_many(seqs[i:i + CHUNK]))
    return out


def main() -> None:
    settings = get_settings()
    pred = get_toxicity_predictor(settings)
    print("=" * 78)
    print("§44 공백 ⓖ — 독성 게이트 계열 4(검출 + 한 수 해킹) + 외부 라벨 첫 데이터점")
    print(f"예측기: {pred.source} · 임계 {THRESH} · 시드 {SEED}")
    print("🔴 재는 양 = ToxinPred3 ML Score. 독성이 아니다.")
    print("=" * 78)
    n_calls = [0]

    # ══ 블록 A ═══════════════════════════════════════════════════════════
    print("\n## A. 외부 라벨 프로브 — 🔴 계열 3 이 아니다(상한값)\n")
    names = list(TOXIC) + list(BENIGN)
    seqs = list(TOXIC.values()) + list(BENIGN.values())
    labels = [1] * len(TOXIC) + [0] * len(BENIGN)
    sc = score_all(pred, seqs)
    n_calls[0] += len(seqs)
    print("| 펩타이드 | 라벨 | ML Score | 게이트 판정 | 맞나 |")
    print("|---|---|---|---|---|")
    for nm, lb, s in zip(names, labels, sc):
        gate = "탈락(독성)" if s > THRESH else "통과(비독성)"
        ok = (s > THRESH) == bool(lb)
        print(f"| {nm} | {'독성' if lb else '무해'} | {s:.3f} | {gate} | "
              f"{'✅' if ok else '🔴 **틀림**'} |")
    tp = sum(1 for lb, s in zip(labels, sc) if lb and s > THRESH)
    fn = sum(1 for lb, s in zip(labels, sc) if lb and s <= THRESH)
    fp = sum(1 for lb, s in zip(labels, sc) if not lb and s > THRESH)
    tn = sum(1 for lb, s in zip(labels, sc) if not lb and s <= THRESH)
    pos = [s for lb, s in zip(labels, sc) if lb]
    neg = [s for lb, s in zip(labels, sc) if not lb]
    auc = sum((1.0 if a > b else 0.5 if a == b else 0.0) for a in pos for b in neg) / (len(pos) * len(neg))
    print(f"\n**민감도 {tp}/{tp+fn} · 특이도 {tn}/{tn+fp} · AUC {auc:.3f}**")
    print(f"**🔴 무해 펩타이드 위양성 {fp}/{len(neg)}** — 이쪽이 학습셋 오염으로 "
          f"설명되지 않는 부분이다(우리 후보가 속한 클래스).")

    # ══ 탈락/통과 구성체 찾기 ════════════════════════════════════════════
    print("\n## 기준선 — 라이브러리 전수에서 탈락 구성체 찾기\n")
    cells, cseqs = [], []
    for cname, cargo in CARGOS.items():
        for sname, smeta in SHUTTLES.items():
            for lname, lmeta in LINKER_LIBRARY.items():
                lk = lmeta["seq"] or ""
                cells.append((cname, cargo, sname, smeta["seq"], lname, lk))
                cseqs.append(cargo + lk + smeta["seq"])
    base = score_all(pred, cseqs)
    n_calls[0] += len(cseqs)
    fails = [(c, s) for c, s in zip(cells, base) if s > THRESH]
    passes = [(c, s) for c, s in zip(cells, base) if s <= THRESH]
    print(f"전수 {len(cells)}셀(화물 2 × 셔틀 12 × 링커 {len(LINKER_LIBRARY)}) → "
          f"**탈락 {len(fails)}셀 ({100*len(fails)/len(cells):.1f}%)** · 통과 {len(passes)}셀")
    fails.sort(key=lambda x: x[1])
    print("\n임계에 가까운 탈락 셀(= 한 수로 넘길 여지가 가장 큰 곳):\n")
    print("| 화물 | 셔틀 | 링커 | ML Score |")
    print("|---|---|---|---|")
    for (cn, _cg, sn, _ss, ln, _lk), s in fails[:8]:
        print(f"| {cn} | {sn} | {ln} | {s:.3f} |")

    # ══ 블록 B ═══════════════════════════════════════════════════════════
    print("\n## B. 계열 4 — 한 수 해킹 (탈락 → 통과)\n")
    targets = fails[:4]
    if not targets:
        print("탈락 셀이 없다 — 이 블록은 측정 불가.")
    print("| 화물 | 셔틀 | 링커 | 원점수 | ① 링커 단일치환 | ② 링커 교체 | "
          "③ 셔틀 단일치환 | ③ 중 제약통과 |")
    print("|---|---|---|---|---|---|---|---|")
    hack_detail = []
    for (cn, cg, sn, ss, ln, lk), s0 in targets:
        # ① 링커 안 단일 치환 (링커가 없으면 측정 불가)
        cand1, meta1 = [], []
        for i in range(len(lk)):
            for a in AA:
                if a == lk[i]:
                    continue
                nl = lk[:i] + a + lk[i + 1:]
                cand1.append(cg + nl + ss)
                meta1.append((i, lk[i], a))
        # ② 링커 라이브러리 교체
        cand2, meta2 = [], []
        for ln2, lm2 in LINKER_LIBRARY.items():
            if ln2 == ln:
                continue
            l2 = lm2["seq"] or ""
            cand2.append(cg + l2 + ss)
            meta2.append(ln2)
        # ③ 셔틀 안 단일 치환 (문헌 보존 위치 제외)
        frozen = {p - 1 for p in CONSERVED_POSITIONS.get(sn, ())}
        cand3, meta3 = [], []
        for i in range(len(ss)):
            if i in frozen:
                continue
            for a in AA:
                if a == ss[i]:
                    continue
                ns = ss[:i] + a + ss[i + 1:]
                cand3.append(cg + lk + ns)
                meta3.append((i, ss[i], a, ns))
        allc = cand1 + cand2 + cand3
        sc_all = score_all(pred, allc)
        n_calls[0] += len(allc)
        s1 = sc_all[:len(cand1)]
        s2 = sc_all[len(cand1):len(cand1) + len(cand2)]
        s3 = sc_all[len(cand1) + len(cand2):]
        c1 = sum(1 for v in s1 if v <= THRESH)
        c2 = sum(1 for v in s2 if v <= THRESH)
        # ③ 의 교차 수 + 그중 제약을 통과하는 것
        base_off = assess_selectivity(ss).off_target_risk
        c3, c3_pass = 0, 0
        for (i, o, a, ns), v in zip(meta3, s3):
            if v > THRESH:
                continue
            c3 += 1
            sim = shuttle_similarity(ns).score
            off = assess_selectivity(ns).off_target_risk
            if sim >= RMT_FLOOR and abs(off - base_off) <= OFF_BAND + 1e-9:
                c3_pass += 1
        hack_detail.append((cn, sn, ln, s0, meta1, s1, meta2, s2, c1, c2, c3, c3_pass))
        f1 = f"**{c1}/{len(cand1)}**" + (" 🔴" if c1 else "") if cand1 else "— (링커 없음)"
        print(f"| {cn} | {sn} | {ln} | {s0:.3f} | {f1} | **{c2}/{len(cand2)}** | "
              f"{c3}/{len(cand3)} | **{c3_pass}** |")

    print("\n🔴 **①②는 셔틀을 바꾸지 않는다 → `rmt_sim`·off 밴드·문헌 동결이 전부 동점이다.**")
    print("   즉 ①②로 넘긴 수는 **제약 셋이 구조적으로 못 잡는다**(통과율 100%).")
    for cn, sn, ln, s0, meta1, s1, meta2, s2, c1, c2, _c3, _cp in hack_detail:
        best = sorted(zip(meta1, s1), key=lambda x: x[1])[:3] if meta1 else []
        if best:
            txt = " · ".join(f"링커 {i+1}번 {o}→{a} ({v:.3f})" for (i, o, a), v in best)
            print(f"   - {cn}+{sn}+{ln} (원 {s0:.3f}): 가장 싼 링커 한 수 — {txt}")
        bl = sorted(zip(meta2, s2), key=lambda x: x[1])[:2]
        if bl:
            print(f"     링커 교체 최선: " + " · ".join(f"`{n}` ({v:.3f})" for n, v in bl))

    # ══ 블록 C ═══════════════════════════════════════════════════════════
    print("\n## C. 계열 4 — 검출 짝: 통과 구성체의 링커 자리에 **독소 조각**을 심는다\n")
    rng = random.Random(SEED)
    # 통과 셀 중 링커가 라이브러리 최장(17aa)인 것들을 고른다 — 같은 길이로 바꿔 심기 위해.
    cand_cells = [(c, s) for c, s in passes if len(c[5]) == 17]
    cand_cells.sort(key=lambda x: x[1])
    picks = cand_cells[::max(1, len(cand_cells) // 6)][:6]
    tox_frag = {nm: sq[:17] for nm, sq in TOXIC.items() if len(sq) >= 17}
    print(f"대상 {len(picks)}셀(링커 길이 17 = 라이브러리 최장) × 독소 조각 "
          f"{len(tox_frag)}종 + 같은 길이 무작위 20개/셀\n")
    print("| 화물 | 셔틀 | 원 ML Score | 독소조각 삽입 최대 | 잡힌 조각 | 무작위 최대 | "
          "무작위 탈락 |")
    print("|---|---|---|---|---|---|---|")
    det_hit = det_tot = 0
    for (cn, cg, sn, ss, ln, lk), s0 in picks:
        tseqs = [cg + f + ss for f in tox_frag.values()]
        rseqs = [cg + "".join(rng.choice(AA) for _ in range(17)) + ss for _ in range(20)]
        sc2 = score_all(pred, tseqs + rseqs)
        n_calls[0] += len(tseqs) + len(rseqs)
        st, sr = sc2[:len(tseqs)], sc2[len(tseqs):]
        hit = sum(1 for v in st if v > THRESH)
        det_hit += hit
        det_tot += len(st)
        rfail = sum(1 for v in sr if v > THRESH)
        print(f"| {cn} | {sn} | {s0:.3f} | {max(st):.3f} | **{hit}/{len(st)}**"
              f"{' 🔴' if hit < len(st) else ''} | {max(sr):.3f} | {rfail}/20 |")
    print(f"\n**독소 조각 삽입 검출: {det_hit}/{det_tot} "
          f"({100*det_hit/max(1,det_tot):.1f}%)**")

    # ══ 결정론 확인 ══════════════════════════════════════════════════════
    probe = cseqs[:40]
    a = score_all(pred, probe)
    b = score_all(pred, probe)
    n_calls[0] += 2 * len(probe)
    assert a == b, "🔴 결정론 깨짐 — 같은 배치가 다른 점수를 냈다(경쟁 조건 의심)"
    print(f"\n✅ 결정론 확인: 같은 40서열 2회 호출 **비트 단위 동일**")
    print(f"총 점수화 서열 {n_calls[0]}개")
    print("=" * 78)
    print("판정은 docs/experiments.md §47 에 적는다. 이 스크립트는 수치만 낸다.")
    print("=" * 78)


if __name__ == "__main__":
    main()
