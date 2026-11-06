"""도구 애블레이션 — 에이전트 도구를 하나씩 빼며 성능·효율을 실측한다.

목적: "이 멀티도구 구조가 맞고 효율적인가"를 반증 가능한 형태로 검증.
     각 도구를 끈 세션과 전체 세션의 (a) 최종 유효점수(eff_bbb = BBB − λ·off_target),
     (b) 심사 승인, (c) 소비 라운드·오라클 평가 수, (d) 토큰 비용을 비교한다.

설계 원칙:
  · 운영 코드는 건드리지 않는다 — 도구 차단·프롬프트 조정은 이 스크립트에서 monkeypatch.
  · 차단된 도구를 프롬프트가 계속 지시하면 불공정하므로, 해당 지시를 제거하고
    "이 세션에서 쓸 수 없는 도구" 제약을 명시한다(없는 도구 호출로 라운드를 낭비하지 않게).
  · 비용 상한을 넘기면 자동 중단하고, 결과는 JSONL로 증분 저장한다.

실행:  /usr/bin/python3 scripts/ablation_tools.py --reps 2 --budget-usd 2.0
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import TOXICITY_THRESHOLD, get_settings          # noqa: E402
from core import optimizer_agent_gemini as oag                    # noqa: E402
from core.optimizer_agent_gemini import OFF_TARGET_PENALTY, get_gemini_agent  # noqa: E402

# gemini-3-flash-preview 공개 단가(USD / 1M tokens). thinking 토큰은 출력으로 과금.
PRICE_IN_PER_TOK = 0.50 / 1_000_000
PRICE_OUT_PER_TOK = 3.00 / 1_000_000

TOOLS = ["evaluate_candidates", "design_candidate", "optimize_developability"]

# analyze_structure 는 2026-11 애블레이션(§11) 결과 **도구에서 제거**됐다 —
# 자기 고유 축(셔틀 노출도)에서조차 퇴행이 없었다. 그래서 -struct 팔도 함께 사라졌다.
CONFIGS = {
    "full":      [],                      # 전체 도구 (대조군)
    "-design":   ["design_candidate"],    # 가설 주도 설계 제거
    "-dev":      ["optimize_developability"],
    "-eval":     ["evaluate_candidates"], # 라이브러리 선택 제거 (설계만으로)
    "lib-only":  ["design_candidate", "optimize_developability"],
}


def eff_bbb(row: dict) -> float:
    return row.get("bbb", 0.0) - OFF_TARGET_PENALTY * row.get("sel_off", 0.0)


# ---------------------------------------------------------------------------
# oracle 기준선 — 라이브러리 전수 조합을 에이전트와 '동일한 지표 정의'로 채점
# ---------------------------------------------------------------------------
def oracle_baseline(agent, cargo: str, cache: Path) -> dict:
    if cache.exists():
        return json.loads(cache.read_text())
    from core.config import LINKER_LIBRARY, SHUTTLES
    combos = [{"label": f"{ln}+{sn}", "linker": LINKER_LIBRARY[ln]["seq"], "shuttle": SHUTTLES[sn]["seq"]}
              for ln in LINKER_LIBRARY for sn in SHUTTLES]
    rows = []
    for i in range(0, len(combos), 24):          # evaluate 배치 상한 24
        _, chunk = agent._evaluate(cargo, {"constructs": combos[i:i + 24]})
        rows.extend(chunk)
        print(f"    oracle {len(rows)}/{len(combos)} 채점…", flush=True)
    valid = [r for r in rows if not r["toxic"]]
    best_eff = max(valid, key=eff_bbb) if valid else None
    best_raw = max(valid, key=lambda r: r.get("bbb", 0.0)) if valid else None
    out = {
        "cargo": cargo, "n_combos": len(rows), "n_valid": len(valid),
        "oracle_eff": eff_bbb(best_eff) if best_eff else 0.0,
        "oracle_eff_label": (best_eff or {}).get("label"),
        "oracle_raw_bbb": (best_raw or {}).get("bbb", 0.0),
        "oracle_raw_label": (best_raw or {}).get("label"),
    }
    cache.write_text(json.dumps(out, ensure_ascii=False, indent=2))
    return out


# ---------------------------------------------------------------------------
# 프롬프트 패치 — 차단 도구 지시 제거 + 제약 명시
# ---------------------------------------------------------------------------
def make_prompt_patch(orig_prompt, dis: set[str]):
    def prompt_patch(c, thr, mr):
        s = orig_prompt(c, thr, mr)
        if "design_candidate" in dis:                      # "[필수] design_candidate" 조항 제거
            i = s.find("**[필수] 유망한 후보를")
            if i != -1:
                j = s.find("\n\n", i)
                s = s[:i] + (s[j + 2:] if j != -1 else "")
        for name in dis:                                   # 도구별 사용법 불릿 제거
            i = s.find(f"- {name}:")
            if i != -1:
                cands = [x for x in (s.find("\n- ", i + 1), s.find("\n\n", i)) if x != -1]
                end = min(cands) if cands else len(s)
                s = s[:i] + s[end + 1:]
        if dis:
            s += ("\n\n**[이번 세션 제약 — 중요]** 다음 도구는 이 세션에서 제공되지 않는다: "
                  + ", ".join(sorted(dis)) +
                  ". 위 설명에 그 도구에 관한 지시가 남아 있어도 무시하고, **제공된 도구만으로** "
                  "목적을 최대한 달성하라(없는 도구를 호출하려 시도하지 마라).")
        return s
    return prompt_patch


# ---------------------------------------------------------------------------
# 한 세션 실행 (도구 집합을 지정)
# ---------------------------------------------------------------------------
def run_session(cargo: str, disabled: list[str], max_rounds: int, rep: int) -> dict:
    settings = get_settings()
    agent = get_gemini_agent(settings, max_rounds=max_rounds)
    if agent is None:
        raise SystemExit("Gemini 에이전트 비활성 — GEMINI_API_KEY / 로컬 엔진 확인")

    dis = set(disabled)

    # ① 도구 선언에서 차단 도구 제거
    orig_decls = agent._function_decls
    agent._function_decls = lambda types: [d for d in orig_decls(types) if d.name not in dis]

    # ② 프롬프트: 차단 도구에 대한 지시 제거 + 제약 명시 (모순 방지)
    orig_prompt = oag._gemini_system_prompt
    oag._gemini_system_prompt = make_prompt_patch(orig_prompt, dis)

    # ③ 도구 호출 수 집계
    calls: dict[str, int] = {}
    orig_dispatch = agent._dispatch

    def dispatch(c, name, args):
        calls[name] = calls.get(name, 0) + 1
        return orig_dispatch(c, name, args)

    agent._dispatch = dispatch

    # ④ 토큰·LLM 호출 집계
    usage = {"llm_calls": 0, "designer_turns": 0, "tok_in": 0, "tok_out": 0}
    orig_gen = agent._gen_retry

    def gen(client, model, contents, config, tries=6):
        resp = orig_gen(client, model, contents, config, tries)
        usage["llm_calls"] += 1
        if getattr(config, "tools", None):
            usage["designer_turns"] += 1
        um = getattr(resp, "usage_metadata", None)
        if um is not None:
            usage["tok_in"] += getattr(um, "prompt_token_count", 0) or 0
            usage["tok_out"] += ((getattr(um, "candidates_token_count", 0) or 0)
                                 + (getattr(um, "thoughts_token_count", 0) or 0))
        return resp

    agent._gen_retry = gen

    rec = {
        "config": "+".join(sorted(dis)) or "full", "disabled": sorted(dis), "rep": rep,
        "cargo": cargo, "max_rounds": max_rounds,
        "n_oracle_evals": 0, "rounds_used": 0, "finished": False, "approved": None,
        "final": None, "error": None, "struct_fail": 0,
        # best_seen = 그 세션이 **직접 채점해서 본** 비독성 후보의 최고 유효점수.
        # final < best_seen 이면 "이미 본 더 좋은 후보를 버리고 제출"한 퇴행(regression)이다.
        "best_seen": None, "best_seen_label": None, "best_from_design": None,
    }
    t0 = time.time()
    try:
        for ev in agent.run(cargo):
            if ev.kind == "evaluation":
                rows = ev.data.get("rows", [])
                rec["n_oracle_evals"] += len(rows)
                is_design = ev.text.startswith("🎨")      # design_candidate 유래 채점
                for r in rows:
                    if r.get("toxic"):
                        continue
                    e = round(eff_bbb(r), 4)
                    if rec["best_seen"] is None or e > rec["best_seen"]:
                        rec["best_seen"] = e
                        rec["best_seen_label"] = r.get("label")
                        rec["best_from_design"] = is_design
            elif ev.kind == "progress":
                rec["rounds_used"] = max(rec["rounds_used"], ev.data.get("round", 0))
            elif ev.kind == "structure":
                if "실패" in ev.text or "오류" in ev.text:
                    rec["struct_fail"] += 1
            elif ev.kind == "choice":
                rec["finished"] = True
            elif ev.kind == "critique":
                rec["approved"] = bool(ev.data.get("approve"))
            elif ev.kind == "optimum":
                d = ev.data
                rec["final"] = {
                    "label": d.get("label"), "linker": d.get("linker"), "shuttle": d.get("shuttle"),
                    "bbb": d.get("bbb"), "sel_off": d.get("sel_off"), "eff_bbb": round(eff_bbb(d), 4),
                    "tox": d.get("tox"), "toxic": d.get("toxic"),
                    "instability": d.get("instability"), "mechanism": d.get("mechanism"),
                    "dev_risk": d.get("dev_risk"), "sel_level": d.get("sel_level"),
                    "n_dev_liab": len(d.get("dev_liabilities") or []),
                }
            elif ev.kind == "error":
                rec["error"] = ev.text
    except Exception as exc:  # noqa: BLE001
        rec["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        oag._gemini_system_prompt = orig_prompt       # 전역 복구

    rec["tool_calls"] = calls
    # 평가 도구가 '설계 통로'로 쓰이려다 거부당한 횟수 — 역할 분리의 작동 증거.
    rec["gate_rejections"] = getattr(agent, "gate_rejections", 0)
    rec.update(usage)
    rec["cost_usd"] = round(usage["tok_in"] * PRICE_IN_PER_TOK
                            + usage["tok_out"] * PRICE_OUT_PER_TOK, 5)
    rec["wall_s"] = round(time.time() - t0, 1)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cargo", default="GSNKGAIIGLM")
    ap.add_argument("--configs", default="full,-design,-dev,-eval,lib-only")
    ap.add_argument("--reps", type=int, default=2)
    ap.add_argument("--max-rounds", type=int, default=10)
    ap.add_argument("--budget-usd", type=float, default=2.0)
    ap.add_argument("--out", default="/tmp/ablation.jsonl")
    ap.add_argument("--oracle-cache", default="/tmp/ablation_oracle.json")
    ap.add_argument("--dry-run", action="store_true",
                    help="API 호출 없이 도구 차단·프롬프트 패치가 올바른지만 점검")
    args = ap.parse_args()

    names = [c.strip() for c in args.configs.split(",") if c.strip()]
    bad = [n for n in names if n not in CONFIGS]
    if bad:
        raise SystemExit(f"알 수 없는 설정: {bad} (가능: {list(CONFIGS)})")

    if args.dry_run:
        _dry_run(names, args.cargo, args.max_rounds)
        return

    out = Path(args.out)
    settings = get_settings()
    probe = get_gemini_agent(settings, max_rounds=args.max_rounds)
    if probe is None:
        raise SystemExit("Gemini 에이전트 비활성 — GEMINI_API_KEY / 로컬 엔진 확인")

    print(f"화물 {args.cargo} · 모델 {settings.gemini_model} · 상한 ${args.budget_usd}")
    print("── oracle(라이브러리 전수) 기준선 ──", flush=True)
    orc = oracle_baseline(probe, args.cargo, Path(args.oracle_cache))
    print(f"  전수 {orc['n_combos']}조합 / 비독성 {orc['n_valid']}개")
    print(f"  oracle 유효점수 최적 : {orc['oracle_eff']:.4f}  ({orc['oracle_eff_label']})")
    print(f"  oracle raw BBB 최적  : {orc['oracle_raw_bbb']:.4f}  ({orc['oracle_raw_label']})\n")

    spent = 0.0
    done = []
    for rep in range(1, args.reps + 1):
        for name in names:
            if spent > args.budget_usd:
                print(f"\n⛔ 비용 상한 ${args.budget_usd} 도달(${spent:.3f}) — 중단")
                _summary(done, orc, spent)
                return
            print(f"▶ [{name}] rep{rep} 실행…", flush=True)
            rec = _run_with_network_retry(args.cargo, CONFIGS[name], args.max_rounds, rep)
            rec["oracle_eff"] = orc["oracle_eff"]
            spent += rec["cost_usd"]
            done.append(rec)
            with out.open("a") as f:
                f.write(json.dumps(rec, ensure_ascii=False) + "\n")
            fin = rec["final"] or {}
            reg = ""
            if fin.get("eff_bbb") is not None and rec.get("best_seen") is not None:
                d = rec["best_seen"] - fin["eff_bbb"]
                if d > 1e-6:
                    reg = f" ⬇퇴행 -{d:.3f}(본최고 {rec['best_seen']:.3f})"
            print(f"  eff={fin.get('eff_bbb')}{reg} BBB={fin.get('bbb')} 승인={rec['approved']} "
                  f"라운드={rec['rounds_used']} 평가={rec['n_oracle_evals']} "
                  f"도구={rec['tool_calls']} 게이트거부={rec['gate_rejections']} "
                  f"${rec['cost_usd']:.4f} {rec['wall_s']}s"
                  + (f"  ⚠️{rec['error'][:60]}" if rec["error"] else ""), flush=True)
    _summary(done, orc, spent)


# 네트워크 끊김으로 세션이 죽는 걸 '실험 결과'로 기록하면 측정이 통째로 망가진다.
# 실측: 두 번의 실행에서 각각 12/18·8/15 세션이 DNS 끊김("nodename nor servname")으로
# 라운드 0에서 즉사했다. 이건 도구 애블레이션과 아무 관계 없는 환경 장애다.
# 그래서 **연결 계열 오류만** 골라 DNS 복구를 기다린 뒤 세션을 처음부터 다시 돌린다.
# (모델·오라클 오류는 재시도하지 않는다 — 그건 진짜 결과다.)
_NET_ERRORS = ("ConnectError", "ReadError", "ConnectTimeout", "nodename nor servname",
               "Can't assign requested address", "Temporary failure in name resolution")


def _net_alive(host="generativelanguage.googleapis.com") -> bool:
    import socket
    try:
        socket.gethostbyname(host)
        return True
    except OSError:
        return False


def _run_with_network_retry(cargo, dis, max_rounds, rep, tries=4):
    """연결 장애로 즉사한 세션은 네트워크 복구를 기다려 다시 돌린다."""
    rec = None
    for attempt in range(1, tries + 1):
        rec = run_session(cargo, dis, max_rounds, rep)
        err = rec.get("error") or ""
        if rec.get("final") or not any(k in err for k in _NET_ERRORS):
            return rec
        if attempt == tries:
            break
        wait = 20 * attempt
        print(f"  ⚠️ 연결 장애로 세션 즉사 — {wait}초 뒤 재시도 ({attempt}/{tries - 1})", flush=True)
        waited = 0
        while waited < wait or not _net_alive():
            time.sleep(5)
            waited += 5
            if waited > 300:            # 5분 넘게 안 돌아오면 포기하고 기록
                break
    return rec


def _dry_run(names, cargo, max_rounds):
    """API 호출 없이: 설정별 도구 선언 목록 + 프롬프트의 차단 도구 언급 잔존 여부를 점검."""
    from google.genai import types
    from core.optimizer_agent_gemini import GeminiOptimizationAgent
    settings = get_settings()
    bare = GeminiOptimizationAgent.__new__(GeminiOptimizationAgent)   # 엔진 초기화 없이 선언만
    bare.settings = settings
    all_decls = [d.name for d in bare._function_decls(types)]
    print(f"전체 도구 선언: {all_decls}\n")
    for name in names:
        dis = set(CONFIGS[name])
        decls = [d for d in all_decls if d not in dis]
        prompt = make_prompt_patch(oag._gemini_system_prompt, dis)(cargo, TOXICITY_THRESHOLD, max_rounds)
        leaks = {t: prompt.count(t) for t in dis if prompt.count(t)}
        mand = "**[필수] 유망한 후보를" in prompt
        print(f"[{name}]")
        print(f"  제공 도구 : {decls}")
        print(f"  프롬프트  : {len(prompt)}자 · 차단도구 언급 잔존 {leaks or '없음'}"
              f" · design_candidate 필수조항 {'있음' if mand else '없음'}")
        tail = prompt[-220:].replace("\n", " ⏎ ") if dis else "(제약문 없음)"
        print(f"  말미      : …{tail}\n")


def _summary(done, orc, spent):
    if not done:
        return
    print("\n" + "=" * 78)
    print(f"{'config':10s} {'n':>2s} {'eff(평균)':>10s} {'%oracle':>8s} {'승인':>5s} "
          f"{'라운드':>6s} {'평가':>5s} {'$':>8s}")
    print("-" * 78)
    by: dict[str, list] = {}
    for r in done:
        by.setdefault(r["config"], []).append(r)
    for name, rs in by.items():
        effs = [(r["final"] or {}).get("eff_bbb") for r in rs]
        effs = [e for e in effs if e is not None]
        avg = sum(effs) / len(effs) if effs else 0.0
        ok = sum(1 for r in rs if r["approved"])
        print(f"{name:10s} {len(rs):2d} {avg:10.4f} "
              f"{(100 * avg / orc['oracle_eff'] if orc['oracle_eff'] else 0):7.1f}% "
              f"{ok:d}/{len(rs):d} "
              f"{sum(r['rounds_used'] for r in rs) / len(rs):6.1f} "
              f"{sum(r['n_oracle_evals'] for r in rs) / len(rs):5.1f} "
              f"{sum(r['cost_usd'] for r in rs):8.4f}")
    print("-" * 78)
    print(f"총 비용 ${spent:.4f}  (oracle 유효점수 {orc['oracle_eff']:.4f})")


if __name__ == "__main__":
    main()
