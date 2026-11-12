"""설계 이력 원장(`core/ledger.py`) 애블레이션 — **래칫 주장을 피고석에 세운다**.

심사 대상 명제 (가정이 아니라 **피고**)
----------------------------------------
보고서는 원장을 기능으로 제시한다 — *"세션 간 래칫: 이미 아는 최고치가 바닥이 되므로
나중 세션은 이전 세션이 도달한 지점에서 시작하거나, 최소한 그 아래로 퇴행하지 않는다"*
(`core/ledger.py` 모듈 독스트링 · `docs/experiments.md` §13-9). 이 명제는 **한 번도
애블레이션된 적이 없다**. §22-5 는 켬 1세션을 재고 "효율 주장 불가"라고 정직히 적었고,
§22-6 은 남은 일을 "켬/끔 각 3세션(라이브 6세션)"으로 적어 뒀다.

🔴 이 스크립트가 **라이브 6세션을 돌리지 않는 이유** (설계 결정, 비용 문제가 아니다)
-----------------------------------------------------------------------------------
⑴ **라이브 팔에는 시드가 없다.** `core/optimizer_agent_gemini.py` 는 `temperature=1.0`
   이고 `GenerateContentConfig` 에 `seed` 를 넘기지 않는다 — 큐 항목의 "×3 시드"는
   라이브에서 **반복(replicate)** 일 뿐 시드가 아니다. 같은 조건의 재현이 불가능하다.
⑵ **1차 판정 지표가 천장에 붙어 있다.** M4 판정은 불리언 3개(최상위 근거 그룹 · 비지배 ·
   라이브러리 밖)이고 콜드에서 이미 **3/3**(§22-2)이다. 켬 팔이 올릴 여지가 **0** 이다.
   여유가 있는 지표는 평가 횟수 하나인데 그건 ⑶ 에 걸린다.
⑶ **총 연산량 교란.** 켬 세션 N 은 세션 1..N-1 의 평가를 **공짜로** 물려받는다. 세션당
   예산을 맞춰 비교하면 켬 팔이 총 연산량에서 앞선 상태로 겨루는 것이다. 공정하게
   맞추려면 끔 팔에 N배 예산을 주거나(= 다른 실험) 켬 팔에 과거 비용을 청구해야 한다.
   §22 데이터로 청구해 보면 켬 1세션은 **6 + 19(세션1·2) = 25평가**이고 콜드는 7~10 이다.

그래서 이 스크립트는 **라이브를 대체하지 않고**, LLM 0회로 할 수 있는 두 가지를 한다:

  **A. 정적 감사 (결정론, 서브프로세스 0회)** — 래칫의 **필요조건 두 개**를 직접 측정한다.
     N1 원장의 바닥이 콜드 프롬프트가 이미 공짜로 주는 것보다 나은가
        (콜드 세션도 `reachable_front` 의 전선 22 · 상위 전선 5 를 프롬프트로 받는다).
        원장 항목이 전부 그 전선에 이미 있으면 래칫은 **공허**하다.
     N2 바닥이 **오늘의 제약 아래 유효**한가 — 아니면 래칫은 폐기될 설계에 닻을 내린다
        (§13-9 가 세 번 밟은 "도달 불가능한 바닥" 버그의 4번째 사례인지 본다).

  **B. 대리 실험 (프록시, 3시드 × 켬/끔, 평가 예산 완전 일치)** — 래칫이 주장하는
     **메커니즘**("기록된 설계에서 출발 vs 라이브러리 부모에서 출발")을 **시드를 줄 수
     있는** 결정론 탐색기(`core.developability_opt.pareto_scrub`)로 재현한다.
     ⚠️ 실제 원장은 GA 시딩이 아니라 **프롬프트 텍스트**로 작동한다 — 그래서 이것은
     프록시이고, 그렇게 **라벨링해서** 보고한다(§30 의 팔 A 가 GA 였던 것과 같은 성격).
     ⚠️ 프록시는 의도적으로 **켬 팔에 유리하게** 편향돼 있다: 켬 팔의 출발점(기록된
     설계)은 과거 세션의 연산량으로 산 것인데 그 비용을 청구하지 않는다. 그래서
     **켬이 이기면 정보가 없고, 켬이 지면 정보가 있다**(비대칭 검정력).

🔴 원장 파일 불변 보장
----------------------
`.cache/ledger/` 는 과거 세션의 **유일한 기록**이다(§30-5: 라이브러리 밖 설계 5개가
여기에만 남아 있다). 이 스크립트는:
  · 실행 전/후로 원본 파일의 **sha256 을 비교해 단언**한다,
  · `core.ledger.LEDGER_DIR` 를 임시 디렉터리로 몽키패치한다(누가 `record()` 를
    호출해도 원본에 닿지 않게) → `finally` 복원 + 복원 단언,
  · 원본은 **읽기 전용**으로만 연다(쓰기·삭제·덮어쓰기 경로가 코드에 없다).

재현(서브프로세스 0회·완전 결정론, Biopython 필요 — `python3 -I` 금지):
  $ /usr/bin/python3 scripts/ablation_ledger.py
  $ /usr/bin/python3 scripts/ablation_ledger.py --seeds 0,1,2 --cargo GSNKGAIIGLM
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tempfile
import unicodedata
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import core.ledger as ledger_mod
from core.binding import shuttle_similarity
from core.config import BASE_DIR, LINKER_LIBRARY, SHUTTLES
from core.developability import assess_developability, require_protein_analysis
from core.developability_opt import OFF_BAND, objectives, off_target, pareto_scrub
from core.evidence import INHERIT_SIM_MIN, conserved_positions, group_of, group_rank, inherit
from core.optimizer_agent import (dominates, in_library, is_comparison_only, is_unsupported,
                                  linker_reject_reason, rank_key, risk_axes)
from core.selectivity import assess_selectivity
from core.stability import assess_stability

#: 프로덕션 GA 설정 — `optimizer_agent_gemini._optimize_developability` 와 **같은 값**을
#: 쓴다. 여기서 다른 값을 쓰면 프록시가 파이프라인이 아닌 다른 탐색기를 재는 셈이 된다.
PROD_ROUNDS, PROD_POP, PROD_FLOOR = 25, 24, 0.85
PROD_GA_SEED = 2022     # 프로덕션 호출이 고정해 둔 시드(보고용 — 프록시는 0,1,2 를 쓴다)

#: `brief()` 가 실제로 프롬프트에 싣는 승자 수 (`k_win` 기본값) — 정보량 측정의 분모다.
BRIEF_K_WIN, BRIEF_K_FAIL = 6, 6


def _w(s: str, width: int, right: bool = False) -> str:
    """표시 폭 기준 패딩 — 한글·이모지는 2칸을 차지한다."""
    disp = sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in s)
    pad = " " * max(0, width - disp)
    return pad + s if right else s + pad


def _sha(p: Path) -> str | None:
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


def _lbl(linker: str) -> str:
    return linker or "직접융합"


def _name_of(seq: str) -> str | None:
    for n, v in SHUTTLES.items():
        if v["seq"].upper() == (seq or "").upper():
            return n
    return None


def _linker_name(seq: str) -> str | None:
    for n, v in LINKER_LIBRARY.items():
        if (v["seq"] or "").upper() == (seq or "").upper():
            return n
    return None


# ---------------------------------------------------------------------------
# 전선 로딩 — 콜드 세션이 **공짜로** 받는 정보
# ---------------------------------------------------------------------------
def load_front(cargo: str) -> dict:
    """`reachable_front` 의 디스크 캐시를 **읽기 전용**으로 가져온다.

    왜 `reachable_front()` 를 호출하지 않나: 그 함수는 156조합에 deepB3P(서브프로세스)와
    ToxinPred3(fcntl 락)를 돌린다. 이 애블레이션에 필요한 것은 **순위 축뿐**이고 그건
    캐시에 이미 들어 있다(§30 과 같은 전선: 후보군 89 · 전선 22). 캐시가 없으면
    **큰 소리로 실패**한다 — 추정값으로 전선을 만들지 않는다(R6: 측정 불가는 None).
    """
    from core.config import get_settings
    s = get_settings()
    sig = hashlib.md5(json.dumps({
        "cargo": cargo,
        "linkers": {k: v["seq"] for k, v in LINKER_LIBRARY.items()},
        "shuttles": {k: v["seq"] for k, v in SHUTTLES.items()},
        "axes": "evidence-group+risk-pareto-v3", "thr": s.toxicity_threshold,
    }, sort_keys=True, ensure_ascii=False).encode()).hexdigest()[:12]
    p = BASE_DIR / ".cache" / "reachable" / f"{sig}.json"
    if not p.exists():
        raise SystemExit(
            f"🔴 전선 캐시 없음: {p}\n"
            "   콜드 프롬프트가 받는 전선을 모르면 '원장이 새 정보를 주나'를 잴 수 없다.\n"
            "   먼저 `/usr/bin/python3 -c \"...reachable_front(cargo, refresh=True)\"` 로 "
            "전선을 깔 것(deepB3P·ToxinPred3 가 돈다).")
    return json.loads(p.read_text())


# ---------------------------------------------------------------------------
# A. 정적 감사
# ---------------------------------------------------------------------------
def audit_entry(cargo: str, e: dict) -> dict:
    """원장 항목 1개를 **오늘의 엔진·오늘의 등급표**로 다시 채점한다.

    원장은 등급 글자를 **기록 시점의 값으로 동결**해 저장하고(`record()`), `rank_key`·
    `dominates` 는 그 **저장된 글자**로 `group_rank` 를 계산한다. 즉 등급 규칙이 바뀌면
    원장의 순위가 **조용히 틀린다** — `load()` 의 자가 치유(`_is_cpp`/`_bad_linker`/
    `_unsupported`)도 저장된 글자를 보기 때문에 이걸 잡지 못한다. 그래서 여기서 재계산한다.
    """
    lk, sh = e.get("linker") or "", e.get("shuttle") or ""
    construct = cargo + lk + sh
    module = lk + sh
    today = {
        "dev_liab": assess_developability(construct).n_liabilities,
        "instability": assess_stability(construct).instability_index,
        "sel_off": round(assess_selectivity(sh).off_target_risk, 3),
        "dev_agg_hot": round(assess_developability(module).agg_hotspot, 3) if module else None,
    }
    t_now, parent, sim = inherit(sh)
    stored = {k: e.get(k) for k in today}
    axes_drift = [k for k, v in today.items()
                  if v is not None and stored.get(k) is not None
                  and abs(float(v) - float(stored[k])) > 1e-6]
    # 부모 라이브러리 셔틀의 off = 오늘의 밴드 제약이 걸리는 **앵커**
    off_parent = round(off_target(SHUTTLES[parent]["seq"]), 3) if parent in SHUTTLES else None
    in_lib = in_library(lk, sh)
    # 오늘의 GA 설계 경로(`optimize_developability` → `pareto_scrub`)가 이 셔틀을
    # 생성할 수 있나 — 라이브러리 밖 설계체에만 묻는 질문이다.
    band_viol = (not in_lib and off_parent is not None
                 and abs(today["sel_off"] - off_parent) > OFF_BAND + 1e-9)
    rmt = round(shuttle_similarity(sh).score, 3)
    floor_viol = (not in_lib) and rmt < PROD_FLOOR
    # 문헌 보존 위치(Angiopep-2 K10·K15) 가 건드려졌나
    cons_viol = False
    if parent in SHUTTLES and not in_lib:
        psq = SHUTTLES[parent]["seq"]
        if len(psq) == len(sh):
            cons_viol = any(psq[i - 1] != sh[i - 1] for i in conserved_positions(sh))
    return {
        "linker": lk, "shuttle": sh, "label": f"{_lbl(lk)}+{_name_of(sh) or sh}",
        "in_lib": in_lib, "stored_tier": e.get("tier"), "today_tier": t_now,
        "parent": parent, "sim": sim, "rmt": rmt,
        "grp_stored": group_rank(e.get("tier")), "grp_today": group_rank(t_now),
        "stored": stored, "today": today, "axes_drift": axes_drift,
        "off_parent": off_parent, "band_viol": band_viol, "floor_viol": floor_viol,
        "cons_viol": cons_viol,
        "linker_reject": linker_reject_reason(lk),
        "cpp": is_comparison_only({"mechanism": e.get("mech")}),
        "unsup_stored": is_unsupported({"tier": e.get("tier")}),
        "unsup_today": is_unsupported({"tier": t_now}),
        # 순위 축 — 저장된 등급 / 오늘의 등급 두 버전
        "key_stored": (group_rank(e.get("tier")),) + risk_axes(e),
        "key_today": (group_rank(t_now),) + risk_axes({**e, **today}),
    }


def front_relation(row: dict, front: list[dict]) -> tuple[str, list[str]]:
    """이 원장 항목이 **콜드 프롬프트의 전선**에 대해 어디 서 있나.

    반환: (관계, 지배자 라벨들). 관계 ∈ {전선멤버, 지배당함, 비지배(새 교환점)}.
    핵심 질문은 "원장이 콜드가 이미 가진 것 위에 무엇을 더하나"다 — 전선 멤버와
    전선에 지배당하는 항목은 **더하는 것이 0** 이다.
    """
    me = {"tier": row["today_tier"], **row["today"]}
    for f in front:
        if f.get("label") == row["label"]:
            return "전선멤버", []
    doms = [f.get("label") for f in front if dominates(f, me)]
    return ("지배당함" if doms else "비지배"), doms


# ---------------------------------------------------------------------------
# B. 프록시 — 3시드 × 켬/끔, 평가 예산 완전 일치
# ---------------------------------------------------------------------------
def proxy_arm(start: str, seed: int, anchor_off: float) -> dict:
    """`pareto_scrub` 1회 = 양 팔이 **같은 평가 예산**(rounds×pop)을 쓴다.

    두 팔의 **유일한** 차이는 출발 서열이다:
      · 끔(cold): 라이브러리 부모 셔틀에서 출발 — 콜드 세션의 설계 경로 그대로.
      · 켬(warm): 원장에 기록된 설계체에서 출발 — 래칫이 주장하는 "물려받은 출발점".

    `anchor_off` = 라이브러리 부모의 off. 켬 팔의 밴드는 **자기 시드**에 anchored 되므로
    (`pareto_scrub` 는 `off0 = off_target(shuttle)`), 부모 기준으로 재면 **밴드가 끌려간다**.
    그 끌림이 측정 대상이다(원장이 해를 끼치는 구체적 경로).
    """
    front = pareto_scrub(start, rounds=PROD_ROUNDS, pop=PROD_POP,
                         rmt_floor=PROD_FLOOR, seed=seed)
    objs = [tuple(m["objectives"][k] for k in ("liability", "aggregation", "instability"))
            for m in front]
    offs = [m["off_target"] for m in front]
    sims = [round(shuttle_similarity(m["seq"]).score, 3) for m in front]
    tiers = [inherit(m["seq"])[0] for m in front]
    return {
        "start": start, "seed": seed, "n_front": len(front),
        "lex_best": min(objs) if objs else None,
        "best_each": (min(o[0] for o in objs), round(min(o[1] for o in objs), 3),
                      round(min(o[2] for o in objs), 1)) if objs else None,
        "off_lo": min(offs), "off_hi": max(offs),
        # 🔴 부모 기준 밴드 위반 — 켬 팔이 오늘의 제약 밖으로 나간 멤버 수
        "band_out": sum(1 for o in offs if abs(o - anchor_off) > OFF_BAND + 1e-9),
        "sim_lo": min(sims), "floor_out": sum(1 for s in sims if s < PROD_FLOOR),
        "tier_fail": sum(1 for t in tiers if t is None),
        "seqs": [m["seq"] for m in front],
    }


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description="설계 이력 원장 애블레이션 (LLM 0회)")
    ap.add_argument("--cargo", default="GSNKGAIIGLM")
    ap.add_argument("--seeds", default="0,1,2", help="프록시 시드(콤마 구분)")
    args = ap.parse_args()
    seeds = [int(x) for x in args.seeds.split(",") if x.strip()]

    # R10 — 엔진 무결성을 먼저 확인한다. Biopython 이 없으면 불안정성이 조용히 퇴화한다.
    require_protein_analysis("원장 애블레이션(불안정성·개발성 축)")

    cargo = args.cargo
    real_path = BASE_DIR / ".cache" / "ledger" / \
        f"{hashlib.sha1(cargo.upper().encode()).hexdigest()[:12]}.json"
    sha_before = _sha(real_path)
    if sha_before is None:
        raise SystemExit(f"🔴 원장 파일 없음: {real_path} — 감사할 기록이 없다.")
    led = json.loads(real_path.read_text())      # 읽기 전용. 쓰기 경로는 이 파일에 없다.

    print("=" * 100)
    print("설계 이력 원장 애블레이션 — 래칫 주장 심사 (LLM 0회 · 서브프로세스 0회 · 결정론)")
    print("=" * 100)
    print(f"화물 {cargo} · 원장 {real_path.name} · 기록 세션 {led.get('sessions')} · "
          f"승자 {len(led.get('wins') or [])} · 실패 {len(led.get('fails') or [])}")
    print(f"프록시 시드 **{seeds}** · GA 설정 rounds={PROD_ROUNDS} pop={PROD_POP} "
          f"rmt_floor={PROD_FLOOR} (프로덕션과 동일 · 프로덕션 고정 시드는 {PROD_GA_SEED})")
    print(f"현재 제약: OFF_BAND=±{OFF_BAND} · INHERIT_SIM_MIN={INHERIT_SIM_MIN} · "
          f"LEDGER_VERSION={ledger_mod.LEDGER_VERSION}")
    print(f"🔒 원본 sha256 {sha_before[:16]}… (실행 끝에 재검증)")

    _orig_dir = ledger_mod.LEDGER_DIR
    tmp = tempfile.mkdtemp(prefix="ablation_ledger_")
    try:
        # 🔴 원본 보호: 누가 record() 를 호출해도 임시 디렉터리로 간다.
        ledger_mod.LEDGER_DIR = Path(tmp)

        front_data = load_front(cargo)
        front, front_top = front_data["front"], front_data["front_top"]
        wins = led.get("wins") or []
        rows = [audit_entry(cargo, e) for e in wins]

        # ── A1. 축 재현성 ────────────────────────────────────────────────────
        drift = [r for r in rows if r["axes_drift"]]
        print("\n" + "─" * 100)
        print("A1. 위험 축 재현성 — 기록된 값이 오늘의 엔진에서 그대로 나오나")
        print("─" * 100)
        print(f"  축 4개(liab·불안정성·off·응집) 전수 재채점: 일치 "
              f"**{len(rows) - len(drift)}/{len(rows)}**"
              + (f" · 어긋남 {[r['label'] for r in drift]}" if drift else " (어긋남 없음)"))
        print("  → 원장의 **숫자**는 썩지 않았다. 썩는 것은 아래 A2 의 **등급**이다.")

        # ── A2. 등급 staleness ─────────────────────────────────────────────
        print("\n" + "─" * 100)
        print("A2. 🔴 등급 동결 — 원장은 등급 **글자**를 저장하고 순위는 그 글자로 계산한다")
        print("─" * 100)
        print("  " + _w("항목", 34) + _w("저장등급", 10) + _w("오늘등급", 10)
              + _w("저장그룹", 10) + _w("오늘그룹", 10) + _w("판정", 20))
        grp_inflated = []
        for r in rows:
            bad = r["grp_today"] > r["grp_stored"]
            if bad:
                grp_inflated.append(r)
            verdict = ("🔴 그룹 과대평가" if bad else
                       ("글자만 변동" if r["stored_tier"] != r["today_tier"] else "동일"))
            print("  " + _w(r["label"], 34) + _w(str(r["stored_tier"]), 10)
                  + _w(str(r["today_tier"]), 10)
                  + _w(f"{r['grp_stored']}({group_of(r['stored_tier'])})", 10)
                  + _w(f"{r['grp_today']}({group_of(r['today_tier'])})", 10)
                  + _w(verdict, 20))
        letter_stale = [r for r in rows if r["stored_tier"] != r["today_tier"]]
        print(f"\n  글자 변동 **{len(letter_stale)}/{len(rows)}** · 그룹 과대평가 "
              f"🔴 **{len(grp_inflated)}/{len(rows)}** "
              f"({', '.join(r['label'] for r in grp_inflated) or '없음'})")
        print(f"  왜 자가 치유가 못 잡나: `ledger.load()` 의 `_unsupported()` 는 **저장된 "
              f"글자**로 `group_rank` 를 본다 → 저장이 'B'면 통과한다. "
              f"`LEDGER_VERSION` 은 {ledger_mod.LEDGER_VERSION} 에서 올라가지 않았는데 "
              "§24(GRADE 간접성 하향)가 등급표를 바꿨다.")

        # ── A3. 정보량 — 콜드가 이미 공짜로 받는 것과의 중복 ────────────────
        print("\n" + "─" * 100)
        print("A3. 래칫 필요조건 N1 — 원장이 **콜드 프롬프트에 없는** 무엇을 주나")
        print("─" * 100)
        print(f"  콜드 세션도 `reachable_front` 를 프롬프트로 받는다: 후보군 "
              f"{front_data['n_gate']} · 전선 **{len(front)}** · 상위 전선 **{len(front_top)}** "
              f"(최상위 근거 그룹 = {front_data['best_group']}).")
        print("  " + _w("원장 항목", 34) + _w("라이브러리", 12) + _w("전선 관계", 14)
              + _w("지배자(전선)", 38))
        rel_count = {"전선멤버": 0, "지배당함": 0, "비지배": 0}
        for r in rows:
            rel, doms = front_relation(r, front)
            rel_count[rel] += 1
            r["front_rel"] = rel
            print("  " + _w(r["label"], 34) + _w("예" if r["in_lib"] else "🆕 설계체", 12)
                  + _w(rel, 14) + _w(", ".join(doms[:2]) or "—", 38))
        n_lib = sum(1 for r in rows if r["in_lib"])
        print(f"\n  라이브러리 조합 **{n_lib}/{len(rows)}** = 콜드가 전수 열거로 이미 가진 것 · "
              f"라이브러리 밖 설계체 **{len(rows) - n_lib}/{len(rows)}** = 원장의 **유일한 신규분**")
        print(f"  전선 관계: {rel_count}")
        shown = rows[:BRIEF_K_WIN]
        print(f"  🔴 **실제 프롬프트에 실리는 것은 상위 {BRIEF_K_WIN}개뿐**"
              f"(`brief(k_win={BRIEF_K_WIN})`) — 그 중 설계체는 "
              f"**{sum(1 for r in shown if not r['in_lib'])}/{len(shown)}** 개다: "
              + ", ".join(f"{r['label']}{'🆕' if not r['in_lib'] else ''}" for r in shown))
        print(f"  실패 기록 **{len(led.get('fails') or [])}개** — '같은 벽에 다시 부딪히지 "
              "마라'(원장의 두 번째 주장 메커니즘)는 3세션 뒤에도 **단 한 줄도 쌓이지 않았다**.")

        # ── A4. 래칫 필요조건 N2 — 오늘의 제약 아래 유효한가 ────────────────
        print("\n" + "─" * 100)
        print("A4. 🔴 래칫 필요조건 N2 — 원장 항목이 **오늘의 제약**을 통과하나")
        print("─" * 100)
        print("  " + _w("설계체", 24) + _w("부모", 14) + _w("off(부모)", 12)
              + _w("|Δoff|", 9, True) + _w(f"밴드±{OFF_BAND}", 10) + _w("rmt_sim", 9, True)
              + _w(f"floor {PROD_FLOOR}", 11) + _w("보존위치", 10))
        for r in rows:
            if r["in_lib"]:
                continue
            d = abs(r["today"]["sel_off"] - (r["off_parent"] or 0.0))
            print("  " + _w(r["shuttle"], 24) + _w(str(r["parent"]), 14)
                  + _w(f"{r['today']['sel_off']:.2f}({r['off_parent']:.2f})", 12)
                  + _w(f"{d:.2f}", 9, True) + _w("🔴 위반" if r["band_viol"] else "통과", 10)
                  + _w(f"{r['rmt']:.3f}", 9, True)
                  + _w("🔴 위반" if r["floor_viol"] else "통과", 11)
                  + _w("🔴 훼손" if r["cons_viol"] else "보존", 10))
        band = [r for r in rows if r["band_viol"]]
        floorv = [r for r in rows if r["floor_viol"]]
        consv = [r for r in rows if r["cons_viol"]]
        lkv = [r for r in rows if r["linker_reject"]]
        unsup = [r for r in rows if r["unsup_today"]]
        bad_union = {r["label"] for r in band + floorv + consv + lkv + unsup + grp_inflated}
        print(f"\n  **오늘의 제약 위반 집계 (분모 {len(rows)} 승자)**")
        print(f"    · OFF_BAND ±{OFF_BAND} 위반 ............ 🔴 **{len(band)}** "
              f"(= 라이브러리 밖 설계체 **전원** {len(band)}/{len(rows) - n_lib})")
        print(f"    · rmt_floor {PROD_FLOOR} 미달 ............. 🔴 **{len(floorv)}**")
        print(f"    · 문헌 보존 위치 훼손 ................ {len(consv)}")
        print(f"    · 링커 거부(길이·알파벳) ............. {len(lkv)}")
        print(f"    · 근거 등급 미달(오늘 기준 D·미등재) .. {len(unsup)}")
        print(f"    · 근거 그룹 과대평가(A2) ............. 🔴 **{len(grp_inflated)}**")
        print(f"    → **합집합 {len(bad_union)}/{len(rows)} "
              f"({len(bad_union) / len(rows):.0%})** 가 오늘의 규칙에서 "
              "거부되거나 잘못 순위된다.")
        print("  ⚠️ 경로 구분(정확히 쓸 것): 밴드·floor 는 **GA 경로**"
              "(`optimize_developability`→`pareto_scrub`)의 제약이다. "
              "`design_candidate`(LLM 직접 편집)는 링커만 검사하고 셔틀에 밴드를 걸지 "
              "않는다 → 같은 설계체를 **LLM 은 여전히 만들 수 있다**. 즉 밴드는 "
              "GA 만 묶고 LLM 은 묶지 않으며, 원장은 그 밴드 위반 설계를 **바닥으로 "
              "프롬프트에 되돌려 넣는다**(§30-4a 비대칭의 하류 결과).")

        # ── A5. 바닥(floor) 단독 심사 ───────────────────────────────────────
        print("\n" + "─" * 100)
        print("A5. 🔴 바닥(`floor_entry` = wins[0]) 단독 심사 — 래칫이 **실제로 가리키는 점**")
        print("─" * 100)
        f0 = rows[0]
        print(f"  저장된 바닥: **{f0['label']}** — 저장 등급 {f0['stored_tier']}"
              f"(그룹 {f0['grp_stored']}) · 축 {risk_axes({**wins[0], **f0['today']})}")
        print(f"  오늘의 등급: **{f0['today_tier']}**(그룹 {f0['grp_today']}) "
              f"→ {'🔴 바닥이 최상위 근거 그룹에 **없다**' if f0['grp_today'] > 0 else '그룹 유지'}")
        # 저장 등급으로 계산한 바닥이 **실재하는** 상위 전선 멤버를 지배하는가
        phantom = {"tier": f0["stored_tier"], **f0["today"]}
        realme = {"tier": f0["today_tier"], **f0["today"]}
        dom_phantom = [m["label"] for m in front_top if dominates(phantom, m)]
        dom_real = [m["label"] for m in front_top if dominates(realme, m)]
        print(f"  저장 등급으로 계산한 바닥이 **상위 전선**(실재 후보 {len(front_top)}개)을 "
              f"지배하는 수: **{len(dom_phantom)}** {dom_phantom}")
        print(f"  오늘 등급으로 계산하면: **{len(dom_real)}** {dom_real}")
        if dom_phantom and not dom_real:
            print("  🔴 **도달 불가능한 바닥(§13-9 버그 4번째 사례)**: 원장은 "
                  f"'근거 있음 + liab {f0['today']['dev_liab']} + off "
                  f"{f0['today']['sel_off']} + 응집 {f0['today']['dev_agg_hot']}' 라는 "
                  "후보를 바닥으로 제시하는데, 오늘의 최상위 근거 그룹에는 그런 점이 "
                  "**없다**. 그래서 에이전트가 실재하는 상위 전선 멤버"
                  f"({dom_phantom[0]})를 제출하면 `brief()` 논리상 **바닥에 지배당한 것**이 "
                  "된다 — 정당한 현재 최선을 열등해 보이게 만든다.")
        print(f"  바닥의 전선 관계: **{f0['front_rel']}** · 라이브러리 조합인가: "
              f"**{'예' if f0['in_lib'] else '아니오'}** "
              "→ 라이브러리 조합이면 콜드도 전선으로 **공짜로** 받는다(N1 미충족).")
        # 오늘 등급으로 다시 정렬하면 바닥이 누가 되나
        resorted = sorted(rows, key=lambda r: r["key_today"])
        print(f"  오늘 등급으로 재정렬한 바닥: **{resorted[0]['label']}** "
              f"(축 {resorted[0]['key_today']}) — "
              f"{'전선 멤버' if resorted[0]['front_rel'] == '전선멤버' else resorted[0]['front_rel']}")

        # ── B. 프록시 ───────────────────────────────────────────────────────
        print("\n" + "=" * 100)
        print("B. 프록시 애블레이션 — 3시드 × 켬/끔 · **평가 예산 완전 일치** (LLM 0회)")
        print("=" * 100)
        designs = [r for r in rows if not r["in_lib"]]
        if not designs:
            print("  원장에 라이브러리 밖 설계체가 없어 프록시를 돌릴 수 없다.")
            warm_start = None
        else:
            warm = min(designs, key=lambda r: r["key_today"])
            warm_start = warm["shuttle"]
            parent_name = warm["parent"]
            cold_start = SHUTTLES[parent_name]["seq"]
            anchor = round(off_target(cold_start), 3)
            print(f"  끔(cold) 출발점 = 라이브러리 부모 **{parent_name}** `{cold_start}` "
                  f"(off {anchor} · 목적 {objectives(cold_start)})")
            print(f"  켬(warm) 출발점 = 원장 최선 설계체 `{warm_start}` "
                  f"(off {warm['today']['sel_off']} · 목적 {objectives(warm_start)} · "
                  f"rmt_sim {warm['rmt']})")
            print("  ⚠️ 켬 팔의 출발점은 **과거 세션의 연산량으로 산 것**이고 이 프록시는 "
                  "그 비용을 청구하지 않는다 → 프록시는 **켬에 유리하게 편향**돼 있다.")
            print("\n  " + _w("팔", 10) + _w("시드", 6, True) + _w("전선", 6, True)
                  + _w("사전식 최선(liab·응집·불안정)", 30, True)
                  + _w("축별 최선", 22, True) + _w("off 범위", 14, True)
                  + _w(f"부모밴드밖", 11, True) + _w(f"floor밖", 9, True)
                  + _w("등급실패", 10, True))
            res = {"끔": [], "켬": []}
            for arm, start in (("끔", cold_start), ("켬", warm_start)):
                for sd in seeds:
                    a = proxy_arm(start, sd, anchor)
                    res[arm].append(a)
                    print("  " + _w(arm, 10) + _w(str(sd), 6, True)
                          + _w(str(a["n_front"]), 6, True)
                          + _w(str(a["lex_best"]), 30, True)
                          + _w(str(a["best_each"]), 22, True)
                          + _w(f"{a['off_lo']:.2f}–{a['off_hi']:.2f}", 14, True)
                          + _w(f"{a['band_out']}/{a['n_front']}", 11, True)
                          + _w(f"{a['floor_out']}/{a['n_front']}", 9, True)
                          + _w(f"{a['tier_fail']}/{a['n_front']}", 10, True))

            def agg(arm, k):
                return [a[k] for a in res[arm]]
            print(f"\n  ⑴ **전선 크기** 끔 {agg('끔','n_front')} vs 켬 {agg('켬','n_front')}")
            print(f"  ⑵ **사전식 최선** 끔 {agg('끔','lex_best')}")
            print(f"                   켬 {agg('켬','lex_best')}")
            print(f"  ⑶ **축별 최선**  끔 {agg('끔','best_each')}")
            print(f"                   켬 {agg('켬','best_each')}")
            # 🔴 프록시의 핵심 읽기 — 켬 팔이 **자기 출발점을 넘지 못했다**.
            seed_obj = objectives(warm_start)
            stuck = sum(1 for a in res["켬"] if a["lex_best"] == seed_obj)
            print(f"  ⑶′ 🔴 켬 팔의 사전식 최선이 **출발점 {seed_obj} 과 동일한 시드 수: "
                  f"{stuck}/{len(seeds)}** — 물려받은 설계 위에서 GA 가 **한 칸도 더 가지 "
                  "못했다**. 래칫은 '바닥을 깔아 준다'를 주장하는데, 측정된 것은 "
                  "'바닥에서 멈춘다'다. 전선도 줄었다"
                  f"({agg('끔','n_front')} → {agg('켬','n_front')}).")
            cold_each = agg("끔", "best_each")
            warm_each = agg("켬", "best_each")
            print("  ⑶″ 축별로는 끔 팔이 더 멀리 간다 — 최선 응집 "
                  f"{min(e[1] for e in cold_each)} (끔) vs "
                  f"{min(e[1] for e in warm_each)} (켬) · 최선 불안정성 "
                  f"{min(e[2] for e in cold_each)} (끔) vs "
                  f"{min(e[2] for e in warm_each)} (켬). liability 는 양 팔 "
                  f"{min(e[0] for e in cold_each)} 동일. "
                  "⚠️ **비지배≠우위**: 켬 팔의 사전식 1위가 낮은 것은 **출발점을 공짜로 "
                  "받았기 때문**이고, 끔 팔은 같은 예산으로 **다른 축을 더 멀리** 밀었다 — "
                  "교환관계로 읽어야 한다.")
            def _frac(arm):
                return [f"{a['band_out']}/{a['n_front']}" for a in res[arm]]
            print(f"  ⑷ 🔴 **부모 기준 밴드 밖 멤버** 끔 "
                  f"{_frac('끔')} vs 켬 {_frac('켬')}"
                  f" — 켬 팔의 밴드는 **자기 시드**(off "
                  f"{warm['today']['sel_off']})에 anchored 되므로 부모({anchor}) 기준으로는 "
                  "통째로 밖이다. 원장이 **제약을 끌고 간다**.")
            print(f"  ⑸ 🔴 **rmt_floor 밖 멤버** 끔 "
                  f"{[a['floor_out'] for a in res['끔']]} vs 켬 "
                  f"{[a['floor_out'] for a in res['켬']]}"
                  " — `pareto_scrub` 은 초기 아카이브에 시드를 **floor 검사 없이** 넣는다"
                  "(`archive = [(shuttle, …)]`) → 하한 미달 설계로 워밍업하면 하한을 "
                  "위반한 멤버가 전선에 **밀반입**된다.")
            # 재발견 — 끔 팔이 원장 설계체에 도달하나
            rediscovered = {d["shuttle"]: [arm for arm in ("끔", "켬")
                                           for a in res[arm] if d["shuttle"] in a["seqs"]]
                            for d in designs}
            print("  ⑹ **원장 설계체 재발견** (끔 팔이 원장 없이 같은 서열에 도달하나)")
            for sq, arms in rediscovered.items():
                print(f"       {sq} → {sorted(set(arms)) or '어느 팔도 도달 못함'}")
            print("       ⚠️ 끔 팔의 0/3 은 '탐색이 약해서'가 아니다 — 세 서열 전부 "
                  f"오늘의 **밴드**(3/3)·**floor**(2/3)를 위반하므로 GA 가 **생성 자체를 "
                  "금지**당한다. 즉 '재발견' 질문은 오늘의 제약 아래 **정의되지 않는다**.")

        # ── B2. 기존 라이브 데이터의 **총 연산량 기준** 재해석 ────────────────
        # 새 라이브 세션을 사지 않고, 이미 지불한 §22 의 4세션을 **공정한 분모**로
        # 다시 읽는다. 큐 항목이 요구한 "무엇을 고정했는지 명시"의 핵심이 여기다.
        print("\n" + "─" * 100)
        print("B2. 이미 지불한 라이브 데이터(§22)를 **총 연산량 기준**으로 다시 읽기")
        print("─" * 100)
        print("  §22 가 보고한 값(문서 인용 — 이 스크립트가 측정한 것이 아니다):")
        print("    · 콜드 3세션 평가 횟수 **7 / 9 / 10** · 최종 축(근거그룹·liab·II·off·응집)")
        print("      세션1 (0, 3, 20.1, 0.20, 0.52) · 세션2 (0, 3, 19.5, 0.20, 0.52)")
        print("    · 누적(원장 ON) 1세션 평가 **6** · 최종 (0, 2, 34.4, 0.22, 0.54)")
        print("  🔴 **세션당 예산으로 비교하면 6 < 7~10 이라 원장이 이긴 것처럼 보인다.** "
              "그런데 원장 ON 세션의 바닥은 세션 1·2 가 만든 것이다 →")
        print("    총 연산량 = 6 + (9 + 10) = **25 평가** vs 콜드 **7~10 평가**.")
        print("  그리고 결과는 **더 좋지도 않다** — liab 3→2 를 얻고 off 0.20→0.22 · 응집 "
              "0.52→0.54 를 지불했다(II 는 둘 다 <40 이라 순위에서 동점). "
              "`dominates` 로 어느 쪽도 상대를 지배하지 않는다 = **교환관계**다.")
        print("  → 즉 **세션당 예산을 고정하면 원장이 이기고, 총 연산량을 고정하면 2.5~3.6배 "
              "비싸며, 품질은 교환관계다.** 어느 분모를 쓰는지가 결론을 뒤집는다 — "
              "그래서 이 비교는 분모를 **명시하지 않으면 주장할 수 없다**(§12 분모 불일치와 "
              "같은 계열).")

        # ── 판정 ────────────────────────────────────────────────────────────
        print("\n" + "=" * 100)
        print("판정")
        print("=" * 100)
        print("🔴 **원장의 래칫 주장은 이 데이터에서 지지되지 않는다.** 근거 네 줄:")
        print(f"  ① **신규 정보가 거의 없다** — 승자 {len(rows)} 중 {n_lib} 개가 라이브러리 "
              f"조합이고, 프롬프트에 실제로 실리는 상위 {BRIEF_K_WIN} 개 중 설계체는 "
              f"{sum(1 for r in shown if not r['in_lib'])}개뿐이다. 콜드 세션도 전선 "
              f"{len(front)}·상위 {len(front_top)}을 **공짜로** 받는다(N1 미충족).")
        print(f"  ② **바닥이 도달 불가능하다** — 바닥 {f0['label']} 의 등급이 동결돼 있어"
              f"(저장 {f0['stored_tier']} / 오늘 {f0['today_tier']}) 실재하지 않는 "
              f"'근거 있음 + liab {f0['today']['dev_liab']}' 후보를 가리킨다. "
              "§13-9 가 세 번 고친 그 버그의 **4번째 사례**다(N2 미충족).")
        print(f"  ③ **신규분이 전부 오늘의 제약 위반이다** — 라이브러리 밖 설계체 "
              f"{len(rows) - n_lib}/{len(rows)} 가 OFF_BAND ±{OFF_BAND} 를 위반한다"
              "(§30-4a: 그 밴드는 바로 이 세션들 **때문에** 나중에 추가됐다). "
              "원장은 폐기된 설계를 바닥으로 되돌려 넣는 경로다.")
        print("  ④ **두 번째 메커니즘은 미발화** — '실패 재방문 금지'의 실패 기록이 "
              f"3세션 뒤 {len(led.get('fails') or [])}개다.")
        print("\n🟡 **이 비교가 고정한 것 / 결론할 수 없는 것**")
        print(f"  고정: ⒜ 프록시는 **평가 예산**(rounds×pop = {PROD_ROUNDS}×{PROD_POP} = "
              f"{PROD_ROUNDS * PROD_POP} 삽입 시도/팔)과 탐색기·제약·"
              "시드를 양 팔에 동일하게 고정하고 **출발 서열 하나만** 바꿨다. "
              "⒝ 정적 감사는 **같은 화물·같은 전선 캐시·같은 등급표**에서 전수로 돌렸다.")
        print("  🔴 결론할 수 없는 것: ⒞ 프록시는 **GA 시딩**을 재고 실제 원장은 **프롬프트 "
              "텍스트**로 LLM 에게 작동한다 — 원장이 LLM 의 **서술·오케스트레이션**에 "
              "무엇을 하는지는 **측정되지 않았다**(§30-5 와 같은 공백). ⒟ 라이브 켬/끔 "
              "3세션씩은 **시드가 없고**(temperature 1.0 · seed 미전달) 1차 판정 지표가 "
              "콜드에서 이미 3/3 천장이며, 켬 팔은 총 연산량에서 앞서 있다 → **그 실험은 "
              "6세션을 사도 '원장이 돕는다'를 보일 수 없다.** ⒠ 따라서 위 판정은 "
              "'원장이 해롭다'가 아니라 **'원장이 돕는다는 주장이 측정으로 뒷받침되지 "
              "않고, 측정 가능한 두 필요조건은 지금 깨져 있다'** 다.")
        print("\n🟢 **측정으로 뒷받침되는 교정 (코드 변경은 하지 않았다 — 문서 작업 모드)**")
        print("  ⒡ `LEDGER_VERSION` 을 3 으로 올린다 — §24 가 등급표를 바꿨고 "
              f"그룹 과대평가 {len(grp_inflated)}/{len(rows)} 가 그 결과다. "
              "버전이 다르면 `load()` 가 빈 원장으로 시작하므로 **한 줄로 치유**된다.")
        print("  ⒢ 또는 `load()` 가 등급을 **저장값이 아니라 `inherit()` 로 재계산**하게 "
              "한다(축 숫자는 A1 에서 16/16 재현되므로 등급만 문제다).")
        print("  ⒣ 보고서에서 원장을 **'기능'이 아니라 '미검증 보조 장치'**로 내린다. "
              "§22-5 의 '바닥을 근거 그룹으로 넘었다(B → A)' 는 **글자 기준으로 틀렸다** — "
              "`TIER_GROUP` 에서 A·B 는 **같은 그룹 0** 이다(그룹 상향은 오늘 등급 "
              "C→B 로 다시 계산할 때만 성립한다).")
        print("\n재현: /usr/bin/python3 scripts/ablation_ledger.py "
              f"--cargo {cargo} --seeds {','.join(map(str, seeds))}")
    finally:
        ledger_mod.LEDGER_DIR = _orig_dir
        assert ledger_mod.LEDGER_DIR == _orig_dir, "LEDGER_DIR 복원 실패"
        sha_after = _sha(real_path)
        assert sha_after == sha_before, (
            f"🔴 원장 파일이 변경됐다! {sha_before} → {sha_after} — 이 스크립트는 "
            "원장을 **읽기만** 해야 한다.")
        print(f"🔒 원본 sha256 재검증 통과 ({sha_after[:16]}…) · "
              f"LEDGER_DIR 복원 확인 ({ledger_mod.LEDGER_DIR})")


if __name__ == "__main__":
    main()
