"""애블레이션 축별 사후 채점 — "도구가 쓸모 있나"를 그 도구의 축에서 묻는다.

왜 이 스크립트가 필요한가
--------------------------
1차 집계(ablation_report.py)는 모든 설정을 **유효점수(eff = BBB − λ·off_target)** 하나로
비교했다. 그런데 eff는 **선택(selection) 축**이다. 도구 4종 중 둘은 애초에 선택을 하지 않는다.

  evaluate_candidates     → 라이브러리에서 고른다            ⇒ eff 가 맞는 축
  design_candidate        → 라이브러리 **밖**을 만든다        ⇒ eff 가 맞는 축 (oracle 초과 가능)
  optimize_developability → liability 를 깎는다              ⇒ eff 로는 **보이지 않음**
  analyze_structure       → 셔틀이 가려졌는지 본다            ⇒ eff 로는 **보이지 않음**

따라서 "도구를 빼도 eff가 안 나빠진다"는 결론은 **도구가 쓸모없다는 뜻이 아니라,
측정 축이 그 도구를 못 본다는 뜻**이다. 올바른 애블레이션은 각 조건(arm)의 산출물을
**그 조건이 보지 못한 축까지 포함해 독립 심판이 전부 채점**하는 것이다.

⚠️ 여기서 하지 않는 일: **목적함수(eff)에 개발성·구조를 더하지 않는다.**
   그 도구가 계산한 값을 그 도구의 점수에 넣으면 순환이고(도구를 쓰면 당연히 오른다),
   우리가 논문에서 비판하는 바로 그 보상 해킹이다. 바꾸는 건 **보고 축**이지 목적함수가 아니다.

실행:  /usr/bin/python3 scripts/ablation_axes.py .cache/ablation/run2.jsonl [--no-struct]
       (ESMFold 공개 API를 쓰므로 --no-struct 로 구조 축을 건너뛸 수 있다)
"""

from __future__ import annotations

import json
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import OFF_TARGET_PENALTY, TOXICITY_THRESHOLD, get_settings  # noqa: E402
from core.predictors import get_predictor  # noqa: E402
from core.toxicity import get_toxicity_predictor  # noqa: E402

ORDER = ["full", "design_candidate", "optimize_developability",
         "evaluate_candidates",
         "design_candidate+optimize_developability"]
LABEL = {
    "full": "전체(대조군)",
    "design_candidate": "−design_candidate",
    "optimize_developability": "−optimize_developability",
    "evaluate_candidates": "−evaluate_candidates",
    "design_candidate+optimize_developability": "evaluate만(lib-only)",
}
# 각 설정이 '잃은' 도구의 고유 축 — 그 축에서 퇴행이 나와야 그 도구가 값을 한 것이다.
OWN_AXIS = {
    "design_candidate": "eff",
    "evaluate_candidates": "eff",
    "optimize_developability": "liab",
}


def score_final(agent, cargo, linker, shuttle):
    """최종 선택을 **독립적으로** 다시 채점한다(그 조건이 못 본 축까지 전부)."""
    row = agent._score_choice(cargo, {"chosen_label": "사후채점",
                                      "chosen_linker": linker or "",
                                      "chosen_shuttle": shuttle or ""})
    if row is None:
        return None
    return {
        "eff": round((row.get("bbb") or 0.0)
                     - OFF_TARGET_PENALTY * (row.get("sel_off") or 0.0), 4),
        "bbb": row.get("bbb"),
        "off": row.get("sel_off"),
        "liab": len(row.get("dev_liabilities") or []),
        "inst": row.get("instability"),
        "toxic": bool(row.get("toxic")),
    }


def score_structure(cargo, linker, shuttle):
    """셔틀 노출도(analyze_structure 의 고유 축). 측정 불가면 None.

    ⚠️ `analyze_construct` 는 접기가 실패해도 `shuttle_exposure=0.0` 을 돌려준다.
    그대로 쓰면 **"셔틀이 파묻힘(0.0)"과 "못 쟀음"이 구분되지 않아** 없는 퇴행을
    만들어낸다(실측: 6세션 중 4개가 0.000 = 전부 ESMFold 실패였다).
    그래서 `res.error` 가 있으면 반드시 None 을 돌려준다.
    """
    try:
        from core.structure import analyze_construct
        res = analyze_construct(cargo, linker or "", shuttle or "")
        if getattr(res, "error", None):
            return None
        if not getattr(res, "pdb", ""):
            return None
        return {"expo": res.shuttle_exposure,
                "plddt": res.shuttle_plddt,
                "low_conf": bool(res.low_confidence)}
    except Exception:  # noqa: BLE001  — ESMFold 공개 API는 자주 죽는다
        return None


def agg(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None, None
    return sum(vals) / len(vals), st.median(vals)


def main():
    path = Path(sys.argv[1])
    do_struct = "--no-struct" not in sys.argv
    recs = [json.loads(l) for l in path.read_text().splitlines() if l.strip()]

    settings = get_settings()
    from core.optimizer_agent_gemini import get_gemini_agent
    agent = get_gemini_agent(settings)
    if agent is None:
        print("⚠️ Gemini 에이전트 비활성 — 채점에는 LLM이 필요 없지만 agent 객체가 필요하다.")
        return
    # 채점에 LLM 호출은 없다(_score_choice 는 로컬 오라클만 쓴다).

    # ── 분모(oracle) 교정 ────────────────────────────────────────────────
    # ablation_tools 의 oracle_eff 는 **CPP 포함** 전수 최적이다(직접융합+SynB1 = 0.795).
    # 그런데 CPP는 비교군이라 에이전트가 고를 수 없다 → 그 분모로 나눈 %oracle 은
    # 애초에 도달 불가능한 기준이고 모든 조건을 부당하게 낮게 보이게 한다.
    # 올바른 분모는 **후보군(RMT) 라이브러리 최적** = 도달 가능 전선의 상한이다.
    cargo0 = recs[0]["cargo"]
    front = agent.reachable_front(cargo0)
    orc_lib = front.get("front_eff_max") or 0.0
    orc_all = recs[0].get("oracle_eff") or 0.0
    print(f"분모(후보군 라이브러리 최적 eff) = {orc_lib:.4f}"
          f"   [참고: CPP 포함 전수 최적 {orc_all:.4f} — 비교군이라 선택 불가]")
    print("→ 100%를 넘으면 라이브러리 밖을 **설계해서** 넘어선 것이다(design_candidate 의 몫).\n")

    print(f"사후 채점 {len(recs)}세션 · 구조축 {'포함' if do_struct else '생략'}\n")
    scored = []
    for i, r in enumerate(recs, 1):
        fin = r.get("final") or {}
        lk, sh = fin.get("linker"), fin.get("shuttle")
        out = {"config": r["config"], "rep": r.get("rep"), "approved": r.get("approved"),
               "label": fin.get("label"), "linker": lk, "shuttle": sh}
        if not sh:
            out["axes"] = None
        else:
            out["axes"] = score_final(agent, r["cargo"], lk, sh)
            if do_struct and out["axes"] is not None:
                sres = score_structure(r["cargo"], lk, sh)
                # 측정 실패는 None 으로 남긴다 — 0.0 으로 채우면 허위 퇴행이 된다.
                out["axes"]["expo"] = None if sres is None else sres["expo"]
                out["axes"]["expo_plddt"] = None if sres is None else sres["plddt"]
                out["axes"]["expo_lowconf"] = None if sres is None else sres["low_conf"]
        scored.append(out)
        print(f"  [{i}/{len(recs)}] {r['config']:<40} {lk or '직접융합'}+{sh or '—'} "
              f"{out['axes']}", flush=True)

    by: dict[str, list] = {}
    for s in scored:
        by.setdefault(s["config"], []).append(s)
    names = [n for n in ORDER if n in by] + [n for n in by if n not in ORDER]

    print("\n### 축별 사후 채점 (독립 심판이 모든 조건의 최종 선택을 전부 채점)\n")
    print("| 설정 | n | eff↑ | %후보군oracle | 개발성 liability↓ | 불안정성↓ | 셔틀 노출도↑ (측정성공) | 잃은 도구의 축 |")
    print("|---|---|---|---|---|---|---|---|")
    for n in names:
        rs = [s for s in by[n] if s["axes"]]
        if not rs:
            print(f"| {LABEL.get(n, n)} | {len(by[n])} | — 최종 후보 없음 |||||")
            continue
        e, _ = agg([s["axes"]["eff"] for s in rs])
        l, _ = agg([s["axes"]["liab"] for s in rs])
        i_, _ = agg([s["axes"]["inst"] for s in rs])
        xs = [s["axes"].get("expo") for s in rs]
        x, _ = agg(xs)
        n_x = sum(1 for v in xs if v is not None)
        own = OWN_AXIS.get(n, "—")
        own_ko = {"eff": "eff", "liab": "liability", "expo": "노출도", "—": "— (대조군)"}[own]
        xcell = "측정실패" if x is None else f"{x:.3f} ({n_x}/{len(rs)})"
        pct = f"{100 * e / orc_lib:.0f}%" if orc_lib else "—"
        print(f"| {LABEL.get(n, n)} | {len(rs)} | {e:.3f} | {pct} | {l:.2f} | {i_:.1f} | "
              f"{xcell} | **{own_ko}** |")

    print("\n읽는 법: 각 행에서 **굵게 표시된 축**이 그 설정이 잃은 도구의 고유 축이다.")
    print("그 축에서만 퇴행이 크면 = 그 도구가 그 축을 실제로 지키고 있었다는 뜻.")
    print("eff에서만 비교하면 개발성·구조 도구는 구조적으로 보이지 않는다.")

    out_path = path.with_suffix(".axes.json")
    out_path.write_text(json.dumps(scored, ensure_ascii=False, indent=1))
    print(f"\n원자료: {out_path}")


if __name__ == "__main__":
    main()
