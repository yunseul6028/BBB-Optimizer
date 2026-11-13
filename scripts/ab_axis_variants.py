#!/usr/bin/env python3
"""대기열 ⑴·⑶ 의 A/B — 순위 축을 바꾸면 **전선이 움직이는가**를 라이브러리 전수로 잰다.

    PYTHONPATH=. /usr/bin/python3 scripts/ab_axis_variants.py

왜 이 측정인가 (§38-4 의 교훈)
---------------------------------------------------------------------------
§38-4 는 `dev_liab` 의 정의를 바꿨을 때 **축 값은 움직였는데 전선 23개 집합이 비트
단위로 동일**했다. 그래서 "축을 고쳤다"는 것만으로 결과가 바뀐다고 말할 수 없다.
이 스크립트는 대기열 두 항목에 대해 같은 A/B 를 돌린다.

재는 변형 (전부 **순위 축만** 건드린다 — 게이트·근거등급·BBB 는 불변)
---------------------------------------------------------------------------
  base : 프로덕션 그대로
  v1   : ⑴ `dev_liab` 의 **응집 핫스팟 플래그 제거**(§43-3 중복 계수)
         근거: `dev_liab` 는 **전체 구성체**로 재므로 화물(GSNKGAIIGLM)의 자체 핫스팟
         0.80 하나로 셔틀과 무관하게 플래그가 켜진다. 같은 양을 응집 축(`dev_agg_hot`,
         설계 모듈)이 이미 세고 있으므로 **순위 축 4개 중 2개가 같은 양을 두 번 센다.**
  v3   : ⑶ AGGRESCAN **프롤린 규칙**(§41-2) — 원문 HS 정의는 "창 안에 프롤린이 있으면
         HS 가 아니다"(응집 차단기)다. 프롤린을 포함한 창을 핫스팟 최대값 후보에서
         제외한다. `dev_agg_hot`(모듈)과 `dev_liab` 의 핫스팟 플래그(구성체) **양쪽**에
         적용한다 — 한쪽만 하면 범위 불일치가 새로 생긴다.
  v13  : 둘 다

🔴 v3 의 경계 사례 — 프롤린 없는 창이 **하나도 없는** 서열
---------------------------------------------------------------------------
T7(`HAIYPRH`)·THR(`THRPPMWSPVWP`) 는 직접융합일 때 모듈의 모든 5aa 창에 P 가 있다.
원문 정의대로면 **HS 가 없다**. 그런데 그걸 축 0.0 으로 쓰면 "응집 위험 최저"가 되어
해당 셔틀이 공짜로 유리해진다 — 내가 고른 관례 하나가 결과를 만든다. 그래서 **두 관례를
모두 돌리고 둘이 갈리는지 본다**:
  floor    : HS 없음 → 축 0.0
  fallback : HS 없음 → 현행대로 전체 창 최대값(= 축 불변)
둘이 같은 답을 주면 관례는 결론에 무관하다. 갈리면 그 사실 자체가 결과다.

🚫 이 스크립트는 코드를 바꾸지 않는다. 축 값을 **바깥에서 다시 계산**해 행에 넣고,
   심사 함수(`dominates`·`rank_key`·`is_candidate`)는 프로덕션 것을 그대로 쓴다.
유료 LLM 0회. 시드 없음(전수).
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.getcwd())
sys.path.append(str(Path(__file__).resolve().parent.parent))

import Bio  # noqa: F401  🔴 없으면 불안정성·off 가 경고 없이 대체값으로 떨어진다(§은 silent degradation)

from core.config import (DEFAULT_CARGO, LINKER_LIBRARY, SHUTTLES, bbb_scoring_seq,
                         get_settings)
from core.developability import (_AGG_WINDOW, _AGGRESCAN_A3V, _A3V_MAX, _A3V_MIN,
                                 _clamp, assess_developability)
from core.evidence import SHUTTLE_EVIDENCE, group_of, group_rank
from core.binding import shuttle_similarity
from core.delivery import assess_delivery
from core.optimizer_agent import (dominates, is_candidate, is_comparison_only,
                                  is_unsupported, rank_key, risk_axes)
from core.predictors import DeepB3PLocalPredictor
from core.selectivity import assess_selectivity
from core.stability import assess_stability
from core.toxicity import ToxinPred3LocalPredictor

HOT_FLAG = "응집 핫스팟"
VARIANTS = ["base", "v1", "v3", "v13"]


def hotspot(seq: str, proline_rule: bool, convention: str = "floor"):
    """정규화 핫스팟 값. `proline_rule` 이면 P 를 포함한 창을 후보에서 제외한다.

    반환 (hot, hot_raw, n_windows, n_eligible). 적격 창이 0개면 convention 에 따른다.
    ⚠️ 프로덕션과 **같은 산술**이어야 한다 — 누적합이 아니라 창 평균을 직접 계산해
    부동소수 누적 오차까지 피한다(값 하나 차이로 지배 관계가 뒤집힐 수 있다).
    """
    vals = [_AGGRESCAN_A3V.get(c, 0.0) for c in seq]
    w = min(_AGG_WINDOW, len(vals))
    wins = [(i, sum(vals[i:i + w]) / w) for i in range(len(vals) - w + 1)]
    elig = [m for i, m in wins if not (proline_rule and "P" in seq[i:i + w])]
    if elig:
        hot_raw = max(elig)
    elif convention == "floor":
        hot_raw = _A3V_MIN          # HS 없음 → 척도 하한 (정규화 0.0)
    else:
        hot_raw = max(m for _, m in wins)   # fallback: 현행 그대로
    hot = round(_clamp((hot_raw - _A3V_MIN) / (_A3V_MAX - _A3V_MIN)), 2)
    return hot, hot_raw, len(wins), len(elig)


def measure(cargo: str, convention: str) -> dict:
    combos = [(ln, sn) for ln in LINKER_LIBRARY for sn in SHUTTLES]
    full = [(cargo + LINKER_LIBRARY[ln]["seq"] + SHUTTLES[sn]["seq"]).upper()
            for ln, sn in combos]
    junc = [bbb_scoring_seq(cargo, LINKER_LIBRARY[ln]["seq"], SHUTTLES[sn]["seq"])
            for ln, sn in combos]
    mod = [(LINKER_LIBRARY[ln]["seq"] + SHUTTLES[sn]["seq"]).upper() for ln, sn in combos]

    uniq_sh = []
    for _, sn in combos:
        if SHUTTLES[sn]["seq"] not in uniq_sh:
            uniq_sh.append(SHUTTLES[sn]["seq"])
    _pr = DeepB3PLocalPredictor().predict_many(junc + uniq_sh)
    bbb = [p.bbb_permeability for p in _pr[:len(combos)]]
    sh_bbb = {s: _pr[len(combos) + i].bbb_permeability for i, s in enumerate(uniq_sh)}
    thr = float(get_settings().toxicity_threshold)
    tox = [t.risk for t in ToxinPred3LocalPredictor(thr).predict_many(full)]

    sel_c, bind_c = {}, {}
    for _, sn in combos:
        s = SHUTTLES[sn]["seq"]
        if s not in sel_c:
            sel_c[s], bind_c[s] = assess_selectivity(s), shuttle_similarity(s)

    rows, noelig = [], {"mod": 0, "full": 0}
    for i, (ln, sn) in enumerate(combos):
        sh = SHUTTLES[sn]["seq"]
        dev_f = assess_developability(full[i])
        dev_m = assess_developability(mod[i]) if mod[i] else dev_f
        stab = assess_stability(full[i])
        ev = SHUTTLE_EVIDENCE.get(sn)

        n_hot = sum(1 for x in dev_f.liabilities if HOT_FLAG in x)   # 구성체 핫스팟 플래그
        liab_nohot = int(dev_f.n_liabilities) - n_hot

        # v3: 프롤린 규칙을 **두 범위에 똑같이** 적용
        hm, _, _, em = hotspot(mod[i] or full[i], True, convention)
        hf, hf_raw, _, ef = hotspot(full[i], True, convention)
        if em == 0:
            noelig["mod"] += 1
        if ef == 0:
            noelig["full"] += 1
        # 프롤린 규칙 하의 구성체 핫스팟 플래그 (프로덕션과 같은 임계 0.70)
        n_hot_v3 = 1 if hf >= 0.70 else 0

        base = {
            "label": f"{ln}+{sn}", "linker": ln, "shuttle": sn,
            "bbb": round(bbb[i], 6), "tox": round(tox[i], 6), "toxic": tox[i] > thr,
            # 🔴 `shuttle_similarity().mechanism` 은 `'세포투과(CPP)'` 라서
            #    `is_comparison_only` 의 `startswith("CPP")` 를 **통과하지 못한다**
            #    (old_vs_new_cargo.py:151 이 기록한 함정 — 한 번 더 당했다).
            #    프로덕션과 **같은 경로**인 `assess_delivery().mechanism`('CPP')을 쓴다.
            "mechanism": assess_delivery(
                sh_bbb[sh], bbb[i], target=bind_c[sh].target,
                mechanism=bind_c[sh].mechanism).mechanism,
            "tier": ev.tier if ev else None,
            "instability": round(float(stab.instability_index), 6),
            "sel_off": round(float(sel_c[sh].off_target_risk), 6),
            # 축 원재료
            "_liab_prod": int(dev_f.n_liabilities), "_liab_nohot": liab_nohot,
            "_n_hot": n_hot, "_n_hot_v3": n_hot_v3,
            "_agg_prod": float(dev_m.agg_hotspot), "_agg_v3": float(hm),
            "_liab_labels": [x for x in dev_f.liabilities],
        }
        for v in VARIANTS:
            nh = base["_n_hot_v3"] if v in ("v3", "v13") else base["_n_hot"]
            liab = base["_liab_nohot"] + (0 if v in ("v1", "v13") else nh)
            agg = base["_agg_v3"] if v in ("v3", "v13") else base["_agg_prod"]
            base[f"dev_liab__{v}"] = liab
            base[f"dev_agg_hot__{v}"] = agg
        rows.append(base)
    return {"cargo": cargo, "convention": convention, "thr": thr,
            "no_eligible_window": noelig, "rows": rows}


def judge(rows: list, variant: str) -> dict:
    """프로덕션 심사 경로를 그대로 재현한다 — 변형 축 값만 바꿔 끼운다."""
    rs = []
    for r in rows:
        q = dict(r)
        q["dev_liab"] = r[f"dev_liab__{variant}"]
        q["dev_agg_hot"] = r[f"dev_agg_hot__{variant}"]
        rs.append(q)
    gate_all = [r for r in rs if not r["toxic"]]
    gate = [r for r in gate_all if is_candidate(r)]
    front = sorted([r for r in gate
                    if not any(dominates(q, r) for q in gate if q is not r)], key=rank_key)
    best_g = min((group_rank(r.get("tier")) for r in gate), default=None)
    front_top = [r for r in front if group_rank(r.get("tier")) == best_g] if gate else []
    groups = {}
    for r in gate:
        groups[group_of(r.get("tier"))] = groups.get(group_of(r.get("tier")), 0) + 1
    return {
        "chain": [len(rs), len(gate_all), len(gate)],
        "front_n": len(front), "front": [r["label"] for r in front],
        "front_top_n": len(front_top), "front_top": [r["label"] for r in front_top],
        "top1": front[0]["label"] if front else None,
        "top1_axes": list(risk_axes(front[0])) if front else None,
        "groups": dict(sorted(groups.items())),
        "cpp_in_gate": sum(1 for r in gate if is_comparison_only(r)),
        "front_liab_dist": _dist(r["dev_liab"] for r in front),
        "front_agg_dist": _dist(round(r["dev_agg_hot"], 2) for r in front),
        "gate_rank": [r["label"] for r in sorted(gate, key=rank_key)],
    }


def _dist(it) -> dict:
    d = {}
    for x in it:
        d[str(x)] = d.get(str(x), 0) + 1
    return dict(sorted(d.items()))


def main() -> None:
    out = {}
    for conv in ("floor", "fallback"):
        m = measure(DEFAULT_CARGO, conv)
        res = {v: judge(m["rows"], v) for v in VARIANTS}
        out[conv] = {"meta": {k: m[k] for k in ("cargo", "convention", "thr",
                                                "no_eligible_window")},
                     "judge": res, "rows": m["rows"]}
    p = Path(os.environ.get("AB_OUT", "/tmp/ab_axis_variants.json"))
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1))
    print(f"저장: {p}")

    for conv in ("floor", "fallback"):
        d = out[conv]
        print(f"\n{'='*74}\n관례 = {conv}  (적격 창 0개: 모듈 {d['meta']['no_eligible_window']['mod']}"
              f" · 구성체 {d['meta']['no_eligible_window']['full']})\n{'='*74}")
        b = d["judge"]["base"]
        print(f"{'':6s} {'체인':>14s} {'전선':>5s} {'front_top':>10s} {'전선집합':>8s} "
              f"{'1위':>26s}  그룹분포")
        for v in VARIANTS:
            j = d["judge"][v]
            same = "동일" if j["front"] == b["front"] else f"Δ{len(set(j['front'])^set(b['front']))}"
            print(f"{v:6s} {'→'.join(map(str, j['chain'])):>14s} {j['front_n']:>5d} "
                  f"{j['front_top_n']:>10d} {same:>8s} {str(j['top1']):>26s}  {j['groups']}")
        for v in VARIANTS:
            j = d["judge"][v]
            print(f"  {v:5s} 전선 liab 분포 {j['front_liab_dist']} · 응집 {j['front_agg_dist']}")


if __name__ == "__main__":
    main()
