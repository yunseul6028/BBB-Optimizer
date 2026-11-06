"""애블레이션 결과 집계 — ablation_tools.py 가 쓴 JSONL을 표로 요약한다.

실행:  /usr/bin/python3 scripts/ablation_report.py <ablation.jsonl> [oracle.json]
"""

from __future__ import annotations

import json
import statistics as st
import sys
from pathlib import Path

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


def main():
    path = Path(sys.argv[1])
    recs = [json.loads(l) for l in path.read_text().splitlines() if l.strip()]
    orc = recs[0].get("oracle_eff") or 0.0
    if len(sys.argv) > 2:
        orc = json.loads(Path(sys.argv[2]).read_text())["oracle_eff"]

    by: dict[str, list] = {}
    for r in recs:
        by.setdefault(r["config"], []).append(r)

    names = [n for n in ORDER if n in by] + [n for n in by if n not in ORDER]
    print(f"oracle 유효점수(라이브러리 전수) = {orc:.4f}\n")
    print("| 설정 | n | eff 평균 | eff 중위 | eff 범위 | %oracle | 심사승인 | 라운드 | 오라클평가 | $/세션 |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for n in names:
        rs = by[n]
        effs = [(r["final"] or {}).get("eff_bbb") for r in rs]
        effs = [e for e in effs if e is not None]
        if not effs:
            print(f"| {LABEL.get(n, n)} | {len(rs)} | — (최종 후보 없음) |||||||")
            continue
        avg, med = sum(effs) / len(effs), st.median(effs)
        ok = sum(1 for r in rs if r["approved"])
        print(f"| {LABEL.get(n, n)} | {len(rs)} | {avg:.3f} | {med:.3f} | "
              f"{min(effs):.3f}–{max(effs):.3f} | {100 * avg / orc:.1f}% | {ok}/{len(rs)} | "
              f"{sum(r['rounds_used'] for r in rs) / len(rs):.1f} | "
              f"{sum(r['n_oracle_evals'] for r in rs) / len(rs):.1f} | "
              f"{sum(r['cost_usd'] for r in rs) / len(rs):.3f} |")

    print("\n### 도구 호출 빈도 (세션당 평균)\n")
    tools = ["evaluate_candidates", "design_candidate", "optimize_developability"]
    print("| 설정 | " + " | ".join(t.replace("_candidates", "").replace("optimize_", "") for t in tools)
          + " | 오류 |")
    print("|---|" + "---|" * (len(tools) + 1))
    for n in names:
        rs = by[n]
        cells = []
        for t in tools:
            cells.append(f"{sum((r.get('tool_calls') or {}).get(t, 0) for r in rs) / len(rs):.1f}")
        err = sum(1 for r in rs if r.get("error"))
        print(f"| {LABEL.get(n, n)} | " + " | ".join(cells) + f" | {err}/{len(rs)} |")

    print("\n### 최종 선택 서열 분포\n")
    for n in names:
        picks = {}
        for r in by[n]:
            f = r["final"] or {}
            k = f"{f.get('linker') or '직접융합'}+{f.get('shuttle', '?')}"
            picks[k] = picks.get(k, 0) + 1
        top = sorted(picks.items(), key=lambda kv: -kv[1])
        print(f"- **{LABEL.get(n, n)}**: " + ", ".join(f"`{k}`×{v}" for k, v in top))

    tot = sum(r["cost_usd"] for r in recs)
    tin = sum(r["tok_in"] for r in recs)
    tout = sum(r["tok_out"] for r in recs)
    print(f"\n총 {len(recs)}세션 · 입력 {tin:,} tok · 출력 {tout:,} tok · **${tot:.3f}** "
          f"(≈{tot * 1460:,.0f}원)")


if __name__ == "__main__":
    main()
