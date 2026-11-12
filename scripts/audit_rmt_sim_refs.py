#!/usr/bin/env /usr/bin/python3
"""큐 #13b — THR·Tet1 의 `rmt_sim = 0.0` 진단.

`core/selectivity.py` 의 off-target 축:
    off = (1 − (1−cationicity)(1−lipophilicity)) × (1 − 0.6·rmt_sim)
여기서 `rmt_sim` 은 `binding.REFERENCE_SHUTTLES` 중 **RMT 라벨 참조들**과의 최대
정규화 BLOSUM62 유사도다. THR·Tet1 은 라이브러리의 RMT 셔틀인데 0.0 이 나온다고 한다.

측정:
  B1. 12종 전체의 rawOR·cationicity·lipophilicity·rmt_sim·off 표
  B2. THR·Tet1 의 **참조별** 정렬 점수 — 원시 BLOSUM 점수까지 열어 원인 규명
  B3. 참조 세트 ↔ 라이브러리 대조 (누락/잉여)
  B4. 🔴 반사실: rmt_sim 이 올바르면 off 가 얼마가 되는가 + 순위/전선 멤버십 변화

실행:  /usr/bin/python3 scripts/audit_rmt_sim_refs.py
⚠️ `python3 -I` 금지. 서브프로세스 0회(BBB 호출 없음 — 순위 축에 BBB 가 없다).
시드: 없음(전부 결정론적 계산).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.binding import REFERENCE_SHUTTLES, _aligner, shuttle_similarity  # noqa: E402
from core.config import LINKER_LIBRARY, SHUTTLES  # noqa: E402
from core.developability import assess_developability  # noqa: E402
from core.evidence import SHUTTLE_EVIDENCE, group_rank  # noqa: E402
from core.selectivity import _clamp, assess_selectivity  # noqa: E402
from core.stability import assess_stability  # noqa: E402

RMT_REFS = {n for n, r in REFERENCE_SHUTTLES.items()
            if "RMT" in r["mech"] or "수용체" in r["mech"]}


def hdr(t: str) -> None:
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


def parts(seq: str) -> dict:
    """selectivity.assess_selectivity 의 내부 항을 그대로 재계산(검증용)."""
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    from core.developability import _net_charge_ph74
    net = _net_charge_ph74(seq)
    dens = round(net / len(seq), 2)
    cat = _clamp(dens / 0.5) if dens > 0 else 0.0
    gravy = round(float(ProteinAnalysis(seq).gravy()), 2)
    lip = _clamp(gravy / 0.8) if gravy > 0 else 0.0
    raw = 1 - (1 - cat) * (1 - lip)
    return {"net": round(net, 2), "dens": dens, "cat": round(cat, 3),
            "gravy": gravy, "lip": round(lip, 3), "raw": round(raw, 3)}


# ---------------------------------------------------------------------------
def b1_table() -> None:
    hdr("B1. 라이브러리 12종 — off-target 분해 (참조 세트 소속 여부와 함께)")
    print("셔틀 | 기전 | 참조세트? | 전하밀도 | GRAVY | rawOR | rmt_sim | 할인(1−0.6s) | off")
    print("---|---|---|---|---|---|---|---|---")
    for n, meta in SHUTTLES.items():
        seq = meta["seq"]
        p = parts(seq)
        r = assess_selectivity(seq)
        mech = "CPP" if "CPP" in meta["note"] else "RMT"
        in_ref = "✅" if n in REFERENCE_SHUTTLES else "🔴 없음"
        print(f"{n} | {mech} | {in_ref} | {p['dens']} | {p['gravy']} | {p['raw']} | "
              f"**{r.rmt_similarity}** | {round(1 - 0.6 * r.rmt_similarity, 3)} | "
              f"**{r.off_target_risk}**")


def b2_per_ref() -> None:
    hdr("B2. THR·Tet1 의 참조별 정렬 — 원시 BLOSUM62 점수까지")
    aln = _aligner()
    for n in ("THR", "Tet1"):
        seq = SHUTTLES[n]["seq"]
        print(f"\n{n} = {seq} (L={len(seq)})")
        print("  참조 | L | RMT? | raw score | self score | 정규화 | 클램프 후")
        print("  ---|---|---|---|---|---|---")
        for rn, rr in REFERENCE_SHUTTLES.items():
            ref = rr["seq"]
            raw = aln.score(seq, ref)
            selfs = aln.score(ref, ref)
            norm = raw / selfs if selfs > 0 else 0.0
            print(f"  {rn} | {len(ref)} | {'RMT' if rn in RMT_REFS else 'CPP'} | "
                  f"{raw:.1f} | {selfs:.1f} | {norm:.3f} | "
                  f"{max(0.0, min(1.0, norm)):.3f}")
        b = shuttle_similarity(seq)
        print(f"  → shuttle_similarity: best={b.best_ref} score={b.score} "
              f"(`_fitness`/`pareto_scrub` 의 rmt_floor=0.85 와 비교)")
    print("\n진단 재료: 12-mer 를 19~30-mer 참조에 **전역(global)** 정렬하면 7~18개 갭이")
    print("필요하고 갭 벌점이 open −10.0 / extend −0.5 이므로 raw 가 **음수**로 떨어진다.")
    print("→ `max(0.0, ...)` 클램프가 음수를 전부 0.0 으로 뭉갠다.")


def b3_ref_vs_lib() -> None:
    hdr("B3. 참조 세트 ↔ 라이브러리 대조")
    lib, ref = set(SHUTTLES), set(REFERENCE_SHUTTLES)
    print(f"라이브러리 {len(lib)}종 · 참조 세트 {len(ref)}종")
    print(f"🔴 라이브러리에 있는데 참조에 **없는** 것: {sorted(lib - ref)}")
    print(f"   참조에만 있는 것(라이브러리 밖): {sorted(ref - lib)}")
    print(f"   RMT 라벨 참조: {sorted(RMT_REFS)}")
    print("\n함의: 참조 세트에 **있는** 셔틀은 자기 자신과 유사도 1.0 → 할인 ×0.4 를 받고,")
    print("      참조 세트에 **없는** 셔틀(THR·Tet1)은 할인 ×1.0 → **전액 물리적화학 벌점**.")
    print("      즉 off 축의 2.5배 차이가 **측정이 아니라 명부 등재 여부**에서 온다.")


# ---------------------------------------------------------------------------
def _rows(off_override: dict | None = None) -> list[dict]:
    """(링커 × 셔틀) 전조합의 **순위 축 행**. BBB·독성은 순위 축이 아니라 생략."""
    out = []
    for ln, lm in LINKER_LIBRARY.items():
        for sn, sm in SHUTTLES.items():
            linker, shuttle = lm["seq"], sm["seq"]
            full = ("GSNKGAIIGLM" + linker + shuttle).upper()
            mod = linker + shuttle
            dev = assess_developability(full)
            devm = assess_developability(mod)
            stab = assess_stability(full)
            off = (off_override or {}).get(sn)
            if off is None:
                off = assess_selectivity(shuttle).off_target_risk
            out.append({"label": f"{ln}+{sn}", "shuttle_name": sn, "linker": ln,
                        "tier": SHUTTLE_EVIDENCE[sn].tier,
                        "dev_liab": dev.n_liabilities,
                        "instability": stab.instability_index,
                        "sel_off": off, "dev_agg_hot": devm.agg_hotspot})
    return out


def _front(rows: list[dict]) -> list[dict]:
    from core.optimizer_agent import dominates
    return [r for r in rows if not any(dominates(q, r) for q in rows if q is not r)]


def b4_counterfactual() -> None:
    hdr("B4. 🔴 반사실 — rmt_sim 이 올바르면 off 가 얼마고, 순위/전선이 바뀌나")
    base = {n: assess_selectivity(m["seq"]).off_target_risk for n, m in SHUTTLES.items()}
    # 반사실: THR·Tet1 을 참조 세트의 RMT 로 등재했다면 자기 유사도 = 1.0
    #        → 다른 모든 라이브러리 RMT 셔틀과 **같은 대우**가 된다.
    fixed = dict(base)
    for n in ("THR", "Tet1"):
        p = parts(SHUTTLES[n]["seq"])
        fixed[n] = round(_clamp(p["raw"] * (1 - 0.6 * 1.0)), 2)

    print("셔틀 | rawOR | 현재 rmt_sim | 현재 off | 교정 rmt_sim | 교정 off | 차이")
    print("---|---|---|---|---|---|---")
    for n in ("THR", "Tet1"):
        p = parts(SHUTTLES[n]["seq"])
        print(f"{n} | {p['raw']} | 0.0 | **{base[n]}** | 1.0 | **{fixed[n]}** | "
              f"{fixed[n] - base[n]:+.2f}")
    print("\n비교 — 참조 세트에 등재된 라이브러리 RMT 셔틀들(이미 rmt_sim=1.0):")
    for n in ("Angiopep-2", "ApoE(159-167)2", "T7", "B6", "CDX", "RVG29", "Leptin30"):
        p = parts(SHUTTLES[n]["seq"])
        print(f"   {n}: rawOR {p['raw']} → off {base[n]} (할인 ×0.4 적용됨)")

    for tag, off_map in (("현재(THR·Tet1 rmt_sim=0)", base),
                         ("교정(THR·Tet1 rmt_sim=1.0)", fixed)):
        rows = _rows(off_map)
        front = _front(rows)
        fnames = sorted({r["shuttle_name"] for r in front})
        print(f"\n[{tag}] 전조합 {len(rows)}개 · 비지배 전선 {len(front)}개")
        print(f"   전선에 등장하는 셔틀: {fnames}")
        print(f"   THR 전선 멤버 수 = {sum(1 for r in front if r['shuttle_name'] == 'THR')}"
              f" · Tet1 = {sum(1 for r in front if r['shuttle_name'] == 'Tet1')}")
        # 등급 그룹 내부 순위(THR 은 C=그룹1, Tet1 은 D=그룹2)
        for target in ("THR", "Tet1"):
            g = group_rank(SHUTTLE_EVIDENCE[target].tier)
            peers = [r for r in rows if group_rank(r["tier"]) == g]
            pf = _front(peers)
            print(f"   그룹{g}({SHUTTLE_EVIDENCE[target].tier}) 내부 전선 {len(pf)}개 중 "
                  f"{target} = {sum(1 for r in pf if r['shuttle_name'] == target)}개")

    print("\n주의: 등급 사전식 1순위 때문에 THR(C=주의) · Tet1(D=근거없음) 은 어차피")
    print("      최종 후보 자격(`is_candidate`)이 없다 — off 교정이 바꾸는 것은")
    print("      **같은 등급 그룹 안의 전선 멤버십과 보고되는 수치**다.")


def b5_scrub_seed() -> None:
    hdr("B5. 파생 영향 — rmt_floor 0.85 게이트 (개발성 최적화의 입구)")
    print("`pareto_scrub`/`scrub_liabilities` 는 `shuttle_similarity(seq).score ≥ 0.85` 를")
    print("기능 보존 제약으로 쓴다. 참조 세트에 없는 셔틀은 **시드 자신이 이미 미달**이다.\n")
    print("셔틀 | 참조세트? | shuttle_similarity(시드) | 0.85 통과?")
    print("---|---|---|---")
    for n, m in SHUTTLES.items():
        s = shuttle_similarity(m["seq"]).score
        print(f"{n} | {'✅' if n in REFERENCE_SHUTTLES else '🔴 없음'} | {s} | "
              f"{'통과' if s >= 0.85 else '🔴 미달'}")


if __name__ == "__main__":
    b1_table()
    b2_per_ref()
    b3_ref_vs_lib()
    b4_counterfactual()
    b5_scrub_seed()
