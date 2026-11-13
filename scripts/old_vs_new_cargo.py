#!/usr/bin/env python3
"""화물 1종을 고정하고 **옛 코드 vs 지금 코드**를 구성체 수준에서 비교한다.

🔴 두 단계로 쓴다 — 측정은 **각 트리에서** 돌리고, 비교는 현재 트리에서 한 번만 한다.

    OLD=/path/to/worktree@d671c60
    cd $OLD   && PYTHONPATH=. /usr/bin/python3 <this> --measure --out /tmp/old.json
    cd $NEW   && PYTHONPATH=. /usr/bin/python3 <this> --measure --out /tmp/new.json
    cd $NEW   && PYTHONPATH=. /usr/bin/python3 <this> --compare /tmp/old.json /tmp/new.json

왜 두 트리인가 — §30 의 단위 불일치를 피하는 방법
--------------------------------------------------------------------------
옛 결정 규칙(`eff = bbb − 0.25 × sel_off`)의 두 입력은 지금도 계산되므로 현재 트리에서
재현할 수 있어 보인다. 하지만 `assess_selectivity`·`assess_developability` 의 **정의가
바뀌었을 수 있다**(liability 는 실제로 §38-1 에서 바뀌었다). 그래서 이 스크립트는
두 트리에서 **같은 원시 입력**을 각각 측정하고, 비교 단계에서 **입력이 동일한지 먼저
검증**한다. 동일하면 차이는 **결정 규칙 하나**로 귀속되고, 다르면 수치 비교를 포기하고
**판정만** 비교한다.

🔴 AUC 의 의미가 두 트리에서 다르다 — 이것이 이 측정의 핵심 결과다
--------------------------------------------------------------------------
옛 코드는 CPP 를 **후보로 두고 점수로** 가렸다(하드 게이트 = 독성 + `bbb ≥ 0.15`).
지금 코드는 `is_candidate` 가 CPP 와 등급 D 를 **라벨로 배제**한다. 따라서 지금 코드의
AUC(RMT>CPP) 는 **구조적으로 1.0** 이고 "판별을 잘한다"가 아니라 **"판별하지 않는다"** 다.
🚫 두 트리의 AUC 를 같은 축에 올려 "개선됐다"로 쓰지 말 것.

유료 LLM 0회. 결정론(시드 없음 — 전수 계산).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# 🔴 **현재 작업 디렉터리를 맨 앞에 둔다.** 이 스크립트는 옛 워크트리에서도 실행되는데,
# 스크립트 위치(새 트리)를 sys.path 앞에 넣으면 옛 트리에서 돌려도 **새 트리 core 를
# 임포트해** 비교가 조용히 무의미해진다(실제로 한 번 당했다). cwd 우선 + 자기 위치는 뒤로.
import os
sys.path.insert(0, os.getcwd())
sys.path.append(str(Path(__file__).resolve().parent.parent))

import Bio  # noqa: F401  🔴 없으면 불안정성·off 가 경고 없이 대체값으로 떨어진다

from core.config import (DEFAULT_CARGO, LINKER_LIBRARY, SHUTTLES, bbb_scoring_seq,
                         get_settings)
from core.predictors import DeepB3PLocalPredictor
from core.toxicity import ToxinPred3LocalPredictor
from core.selectivity import assess_selectivity
from core.developability import assess_developability
from core.stability import assess_stability
from core.binding import shuttle_similarity

OLD_OFF_PENALTY = 0.25   # 옛 config.OFF_TARGET_PENALTY (현재 트리에선 폐지 상수)
OLD_BBB_FLOOR = 0.15     # 옛 config.MIN_BBB_SIGNAL (현재 트리에선 폐지)


def _tree_id() -> dict:
    """이 작업 트리가 옛 구조인지 새 구조인지 — 상수·모듈 존재로 판별한다."""
    import core.config as C
    info = {"has_min_bbb_signal": hasattr(C, "MIN_BBB_SIGNAL"),
            "has_off_target_penalty": hasattr(C, "OFF_TARGET_PENALTY"),
            "model_max_len": getattr(C, "MODEL_MAX_LEN", None),
            "junction_flank": getattr(C, "JUNCTION_FLANK", None)}
    try:
        from core.evidence import SHUTTLE_EVIDENCE
        info["has_evidence"] = True
        info["n_evidence"] = len(SHUTTLE_EVIDENCE)
    except Exception:
        info["has_evidence"] = False
    info["structure"] = "new" if info["has_evidence"] else "old"
    info["core_path"] = str(Path(C.__file__).resolve())
    info["cwd"] = os.getcwd()
    return info


def measure(cargo: str) -> dict:
    tree = _tree_id()
    combos = [(ln, sn) for ln in LINKER_LIBRARY for sn in SHUTTLES]
    seqs_full, seqs_junc, seqs_mod = [], [], []
    for ln, sn in combos:
        lk, sh = LINKER_LIBRARY[ln]["seq"], SHUTTLES[sn]["seq"]
        seqs_full.append((cargo + lk + sh).upper())
        seqs_junc.append(bbb_scoring_seq(cargo, lk, sh))
        seqs_mod.append((lk + sh).upper())          # 설계 모듈(링커+셔틀)

    uniq_sh = []
    for _, sn in combos:
        if SHUTTLES[sn]["seq"] not in uniq_sh:
            uniq_sh.append(SHUTTLES[sn]["seq"])

    dp = DeepB3PLocalPredictor()
    preds = dp.predict_many(seqs_junc + uniq_sh)
    bbb = [p.bbb_permeability for p in preds[:len(combos)]]
    sh_bbb = {s: preds[len(combos) + i].bbb_permeability for i, s in enumerate(uniq_sh)}
    thr = float(get_settings().toxicity_threshold)
    tox = [t.risk for t in ToxinPred3LocalPredictor(thr).predict_many(seqs_full)]

    # 셔틀 단독 축(셔틀마다 1회만)
    sel_cache, bind_cache = {}, {}
    for _, sn in combos:
        s = SHUTTLES[sn]["seq"]
        if s not in sel_cache:
            sel_cache[s] = assess_selectivity(s)
            bind_cache[s] = shuttle_similarity(s)

    rows = []
    for i, (ln, sn) in enumerate(combos):
        lk, sh = LINKER_LIBRARY[ln]["seq"], SHUTTLES[sn]["seq"]
        selr, bind = sel_cache[sh], bind_cache[sh]
        dev_full = assess_developability(seqs_full[i])
        dev_mod = assess_developability(seqs_mod[i]) if seqs_mod[i] else None
        stab = assess_stability(seqs_full[i])
        row = {
            "label": f"{ln}+{sn}", "linker_name": ln, "shuttle": sn,
            "bbb": round(bbb[i], 6), "shuttle_bbb": round(sh_bbb[sh], 6),
            "tox": round(tox[i], 6),
            "sel_off": round(float(selr.off_target_risk), 6),
            "selectivity": round(float(selr.selectivity), 6),
            "mechanism": getattr(bind, "mechanism", None),
            "instability": round(float(stab.instability_index), 6),
            "dev_liab": int(dev_full.n_liabilities),
            "net_charge": round(float(dev_full.net_charge), 6),
            "dev_agg_hot_full": _g(dev_full, "agg_hotspot"),
            "dev_agg_hot": _g(dev_mod, "agg_hotspot") if dev_mod else _g(dev_full, "agg_hotspot"),
        }
        if tree["has_evidence"]:
            from core.evidence import SHUTTLE_EVIDENCE
            ev = SHUTTLE_EVIDENCE.get(sn)
            row["tier"] = ev.tier if ev else None
        rows.append(row)
    tree["toxicity_threshold"] = thr
    return {"tree": tree, "cargo": cargo, "n": len(rows), "rows": rows}


def _g(obj, name):
    v = getattr(obj, name, None)
    return round(float(v), 6) if isinstance(v, (int, float)) else None


def auc(pos, neg) -> float:
    """순위 기반 AUC = P(pos > neg) (Mann-Whitney). bench_selectivity.py 와 같은 정의."""
    c = t = 0
    for p in pos:
        for n in neg:
            t += 1
            c += 1 if p > n else (0.5 if p == n else 0)
    return c / t if t else 0.5


# 🔴 기전 분류의 출처는 **셔틀 `note` 의 RMT/CPP 토큰**이다(core/config.py 주석 참고).
#    `shuttle_similarity(...).mechanism` 은 이 용도로 신뢰할 수 없었다(전수 RMT 로 나왔다).
#    두 트리에서 CPP 집합이 동일함을 확인했다: {Penetratin, SynB1, TAT} · RMT 9종.
def _cpp_names() -> set:
    from core.config import SHUTTLES
    return {n for n in SHUTTLES if "CPP" in SHUTTLES[n]["note"]}


_CPP = None


def _is_cpp(r) -> bool:
    global _CPP
    if _CPP is None:
        _CPP = _cpp_names()
    return r.get("shuttle") in _CPP


def compare(p_old: str, p_new: str) -> None:
    old, new = (json.loads(Path(p).read_text()) for p in (p_old, p_new))
    O = {r["label"]: r for r in old["rows"]}
    N = {r["label"]: r for r in new["rows"]}
    print("=" * 78)
    print(f"화물: {old['cargo']}   (새 트리: {new['cargo']})")
    print(f"옛 트리: {old['tree']}")
    print(f"새 트리: {new['tree']}")
    assert old["cargo"] == new["cargo"], "화물이 다르면 비교 불가"
    print(f"조합 수: 옛 {old['n']} · 새 {new['n']}")

    # ── 1. 공유 원시 입력이 동일한가 (§30 회피의 전제) ─────────────────────
    print("\n" + "=" * 78)
    print("1. 공유 원시 입력 동일성 검증 — 다르면 수치 비교를 포기한다")
    print("=" * 78)
    shared = ["bbb", "shuttle_bbb", "tox", "sel_off", "selectivity",
              "instability", "dev_liab", "net_charge", "dev_agg_hot_full"]
    identical = {}
    for f in shared:
        diffs = [k for k in O if k in N and O[k].get(f) != N[k].get(f)]
        identical[f] = not diffs
        mark = "✅ 동일" if not diffs else f"🔴 {len(diffs)}/{len(O)} 다름"
        extra = ""
        if diffs:
            k = diffs[0]
            extra = f"   예: {k}  옛 {O[k].get(f)} → 새 {N[k].get(f)}"
        print(f"  {f:<20} {mark}{extra}")

    # ── 2. 하드 게이트 ────────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("2. 하드 게이트 — 몇 개가 최종 후보로 남나")
    print("=" * 78)
    TOX_THR = old["tree"].get("toxicity_threshold") or 0.38
    assert TOX_THR == (new["tree"].get("toxicity_threshold") or 0.38), "독성 임계가 다르다"
    o_tox = [r for r in old["rows"] if r["tox"] <= TOX_THR]
    o_gate = [r for r in o_tox if r["bbb"] >= OLD_BBB_FLOOR]
    n_tox = [r for r in new["rows"] if r["tox"] <= TOX_THR]
    n_gate = [r for r in n_tox if not _is_cpp(r) and r.get("tier") in ("A", "B", "C")]
    print(f"  옛: 독성 통과 {len(o_tox)}/{old['n']} → BBB ≥ {OLD_BBB_FLOOR} 통과 "
          f"**{len(o_gate)}/{old['n']}**")
    print(f"  새: 독성 통과 {len(n_tox)}/{new['n']} → CPP·등급D 라벨 배제 후 "
          f"**{len(n_gate)}/{new['n']}**")
    print(f"\n  🔴 게이트의 **종류가 다르다**: 옛 = 점수 하한(bbb) · 새 = 라벨 배제(CPP·등급 D).")

    # ── 3. AUC(RMT > CPP) ─────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("3. AUC(RMT > CPP) — 구성체 수준, 같은 화물")
    print("=" * 78)
    rmt_o = [r for r in old["rows"] if not _is_cpp(r)]
    cpp_o = [r for r in old["rows"] if _is_cpp(r)]
    print(f"  평가셋: RMT 구성체 {len(rmt_o)} · CPP 구성체 {len(cpp_o)}")
    print(f"\n  {'축':<34}{'RMT평균':>10}{'CPP평균':>10}{'AUC':>8}")
    axes = [
        ("옛 결정 규칙  eff = bbb − 0.25·off", lambda r: r["bbb"] - OLD_OFF_PENALTY * r["sel_off"]),
        ("옛 1차 신호   bbb(융합 deepB3P)", lambda r: r["bbb"]),
        ("선택성 축     selectivity", lambda r: r["selectivity"]),
        ("off 축        −sel_off", lambda r: -r["sel_off"]),
    ]
    if not cpp_o:
        print("  🔴 CPP 구성체가 0개다 — 라벨 추출 실패. 비교를 중단한다.")
        return
    for name, f in axes:
        p, q = [f(r) for r in rmt_o], [f(r) for r in cpp_o]
        print(f"  {name:<34}{sum(p)/len(p):>10.3f}{sum(q)/len(q):>10.3f}{auc(p, q):>8.3f}")
    print("\n  🔴 새 코드는 이 표에 올릴 수 없다 — CPP 를 라벨로 배제하므로 점수 비교가 없다.")
    print("     새 코드의 AUC(RMT>CPP) 는 정의상 1.000 이고 **판별력의 증거가 아니다**.")

    # ── 4. 1위가 바뀌는가 ─────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("4. 1위 — 두 규칙이 같은 구성체를 고르나")
    print("=" * 78)
    o_best = max(o_gate, key=lambda r: r["bbb"] - OLD_OFF_PENALTY * r["sel_off"]) if o_gate else None
    if o_best:
        eff = o_best["bbb"] - OLD_OFF_PENALTY * o_best["sel_off"]
        print(f"  옛 1위: **{o_best['label']}**  (eff {eff:.3f} = bbb {o_best['bbb']:.3f} "
              f"− 0.25×off {o_best['sel_off']:.3f})")
        nb = N.get(o_best["label"])
        if nb:
            cand = (nb["tox"] <= TOX_THR) and not _is_cpp(nb) and nb.get("tier") in ("A", "B", "C")
            print(f"          새 체계에서: 등급 {nb.get('tier')} · 기전 {nb.get('mechanism')} "
                  f"· 최종 후보 자격 {'있음' if cand else '**없음**'}")
    from core.optimizer_agent import rank_key, dominates
    for r in new["rows"]:
        r["toxic"] = r["tox"] > TOX_THR
    n_best = min(n_gate, key=rank_key) if n_gate else None
    if n_best:
        print(f"  새 1위: **{n_best['label']}**  (rank_key {rank_key(n_best)})")
        ob = O.get(n_best["label"])
        if ob:
            eff = ob["bbb"] - OLD_OFF_PENALTY * ob["sel_off"]
            rank = 1 + sum(1 for r in o_gate
                           if (r["bbb"] - OLD_OFF_PENALTY * r["sel_off"]) > eff)
            print(f"          옛 체계에서: eff {eff:.3f} · 게이트 통과 {len(o_gate)}개 중 "
                  f"**{rank}위** · bbb {ob['bbb']:.3f} "
                  f"({'하한 통과' if ob['bbb'] >= OLD_BBB_FLOOR else '**하한 미달**'})")

    # ── 5. 전선 크기 ──────────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("5. 비지배 전선 크기 (새 체계 기준)")
    print("=" * 78)
    front = [q for q in n_gate if not any(dominates(r, q) for r in n_gate if r is not q)]
    print(f"  새: 후보 {len(n_gate)} 중 비지배 **{len(front)}개**")
    print("  옛 체계에는 위험 축 Pareto 가 없다(단일 eff 스칼라) → 전선 개념 자체가 없다.")
    print("\n" + "=" * 78)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--measure", action="store_true")
    ap.add_argument("--compare", nargs=2, metavar=("OLD", "NEW"))
    ap.add_argument("--out")
    ap.add_argument("--cargo", default=DEFAULT_CARGO)
    a = ap.parse_args()
    if a.measure:
        res = measure(a.cargo)
        Path(a.out).write_text(json.dumps(res, ensure_ascii=False, indent=1))
        print(f"[{res['tree']['structure']}] {res['n']}개 측정 → {a.out}")
    elif a.compare:
        compare(*a.compare)
    else:
        ap.error("--measure 또는 --compare 중 하나")
