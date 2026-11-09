"""자율 설계 에이전트 — 공용 도구 백엔드 베이스.

에이전트의 **실제 생물학 계산 도구**(evaluate/structure)를 담는 베이스 클래스.
브레인(LLM) 루프는 `optimizer_agent_gemini.py`의 `GeminiOptimizationAgent`가 이 클래스를
상속해 구현한다. (Claude 브레인 버전은 `with-claude` 브랜치에 보존 — main은 Gemini 전용.)

  - evaluate  : deepB3P(BBB) + ToxinPred3(독성) + 안정성 + 수용체유사 + 개발성(용해도·응집 포함) + 선택성 배치 채점
  - structure : ESMFold 폴딩 → 셔틀 구조 노출도 (느림 ~10초, 최종 후보에만)

셔틀은 **de-novo 생성하지 않는다**(검증 리간드만). 라이브러리 밖 탐색은 잔기 편집(design_candidate)으로
수행하며, 개발성 다목적 최적화는 optimize_developability가 담당한다.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .config import (
    LINKER_LIBRARY,
    SHUTTLES,
    Settings,
    bbb_scoring_seq,
)
from .predictors import get_predictor
from .toxicity import get_toxicity_predictor

MAX_CANDIDATES_PER_STEP = 24

# ── 도구 역할 분리: evaluate는 '선택', design은 '생성' ─────────────────────
# `evaluate_candidates`는 원래 **라이브러리 조합 채점(선택)** 도구인데, 임의 서열을 그대로
# 받아주고 있었다. 그래서 에이전트가 **평가 도구를 설계 통로로 전용**했고, 애블레이션에서
# design_candidate를 차단한 조건조차 라이브러리 밖 서열(`G+HAIYPRH`)을 최종 선택했다
# → 설계 도구의 기여를 **측정 자체가 불가능**했다(실측: evaluate만 조건이 후보군 oracle의 112%).
# 그래서 evaluate는 **라이브러리 멤버만** 받고, 라이브러리 밖은 design_candidate로 보낸다.
# (design_candidate는 가설·예측을 요구하므로 '설계'가 기록에 남는다.)
LIB_LINKERS = frozenset((v["seq"] or "").upper() for v in LINKER_LIBRARY.values())
LIB_SHUTTLES = frozenset(v["seq"].upper() for v in SHUTTLES.values())


def in_library(linker: str, shuttle: str) -> bool:
    """이 조합이 라이브러리 전수 공간 안에 있는가(= evaluate로 채점 가능한가)."""
    return (linker or "").upper() in LIB_LINKERS and (shuttle or "").upper() in LIB_SHUTTLES

# ── 셔틀 역할: CPP는 '비교군', RMT는 '후보군' ──────────────────────────────
# CPP(TAT·SynB1·Penetratin 등 양이온성 막투과)는 **비특이 흡수 비교군(baseline)**으로만 쓴다.
# 라이브러리에는 남겨 두되(비교 근거가 되므로) **최종 결과·순위·러닝 베스트에서는 제외**한다.
# 근거: deepB3P가 양이온성 CPP를 고평가하는 편향이 있어(실측 CPP 평균 BBB 0.357 vs RMT 0.058),
# 효능 점수만 좇으면 비특이 후보로 수렴한다. 선택성·off-target 실측도 CPP가 일관되게 나쁘다.
def is_comparison_only(row: dict) -> bool:
    """이 후보가 **비교군(CPP)**인가 — True면 최종 결과로 올리지 않는다."""
    return str(row.get("mechanism") or "").upper().startswith("CPP")


def is_candidate(row: dict) -> bool:
    """최종 후보가 될 수 있는가 — 비독성 + 비교군(CPP) 아님."""
    return (not row.get("toxic")) and not is_comparison_only(row)


@dataclass
class AgentEvent:
    kind: str  # reasoning|text|plan|evaluation|structure|generation|progress|choice|critique|final|optimum|error
    text: str = ""
    data: dict = field(default_factory=dict)


class OptimizationAgent:
    """도구 백엔드(공용). LLM 루프는 서브클래스(GeminiOptimizationAgent)가 구현."""

    # max_rounds 기본값은 팩토리 `get_gemini_agent`(=10, 안전 상한)와 맞춘다.
    # 근거는 optimizer_agent_gemini.py 의 get_gemini_agent 주석 참고.
    def __init__(self, settings: Settings, max_rounds: int = 10):
        self.settings = settings
        self.max_rounds = max_rounds
        self.predictor = get_predictor(settings)
        self.tox_predictor = get_toxicity_predictor(settings)

    # ---- 도구 백엔드 (실제 생물학 계산) ------------------------------------
    def _evaluate(self, cargo, tool_input, library_only: bool = False):
        """후보들을 채점한다.

        library_only=True 면 **라이브러리 조합만** 받는다(= evaluate_candidates 도구).
        라이브러리 밖 서열은 채점하지 않고 design_candidate로 되돌려보낸다 — 평가 도구가
        설계 통로로 쓰이는 것을 막아 설계 도구의 기여를 측정 가능하게 만든다.
        design_candidate·최종 후보 채점 경로는 library_only=False 로 호출해 제약을 받지 않는다.
        """
        items = (tool_input.get("assemblies") or tool_input.get("constructs") or [])[:MAX_CANDIDATES_PER_STEP]
        clean, rejected = [], []
        for c in items:
            lk = "".join(ch for ch in (c.get("prefix", c.get("linker", "")) or "").upper() if ch.isalpha())
            sh = "".join(ch for ch in (c.get("suffix", c.get("shuttle", "")) or "").upper() if ch.isalpha())
            if library_only and not in_library(lk, sh):
                rejected.append(f"{c.get('label', '?')}({lk or '직접융합'}+{sh})")
                continue
            clean.append({"label": c.get("label", "?"), "linker": lk, "shuttle": sh,
                          "sequence": cargo + lk + sh})
        reject_note = ""
        if rejected:
            # 게이트 발동 횟수 — 역할 분리가 실제로 작동하는지 보는 계측값
            # (에이전트가 평가 도구로 설계를 시도한 횟수. 애블레이션이 이 값을 읽는다.)
            self.gate_rejections = getattr(self, "gate_rejections", 0) + len(rejected)
            reject_note = (
                f"\n\n⚠️ **라이브러리 밖이라 채점하지 않은 조합 {len(rejected)}개**: "
                + ", ".join(rejected[:8]) + ("…" if len(rejected) > 8 else "")
                + "\n`evaluate_candidates`는 **라이브러리 전수 조합(선택)** 전용 도구다. "
                  "직접 편집·설계한 서열은 **`design_candidate`로 제출**하라 "
                  "(가설·예측을 함께 적어야 한다). 라이브러리 링커·셔틀 목록은 시스템 프롬프트 참고.")
        if not clean:
            return (("채점한 후보가 없다." + reject_note).strip(), [])
        from .binding import shuttle_similarity
        from .stability import assess_stability
        from .developability import assess_developability
        from .selectivity import assess_selectivity
        from .delivery import assess_delivery
        # 융합 연결부위 deepB3P + 셔틀 단독 deepB3P를 **한 배치**로 추론(전달 축 분해용).
        junction_seqs = [bbb_scoring_seq(cargo, c["linker"], c["shuttle"]) for c in clean]
        uniq_shuttles: list[str] = []
        for c in clean:
            if c["shuttle"] and c["shuttle"] not in uniq_shuttles:
                uniq_shuttles.append(c["shuttle"])
        preds = self.predictor.predict_many(junction_seqs + uniq_shuttles)
        bbb = preds[: len(clean)]
        shuttle_bbb = {s: preds[len(clean) + i].bbb_permeability for i, s in enumerate(uniq_shuttles)}
        tox = self.tox_predictor.predict_many([c["sequence"] for c in clean])
        thr = self.settings.toxicity_threshold
        rows, lines = [], ["label | primary | shuttleBBB | mech | preserv | penalty | index | match | status | assembly"]
        for c, p, t in zip(clean, bbb, tox):
            toxic = t.risk > thr
            bind = shuttle_similarity(c["shuttle"])
            stab = assess_stability(c["sequence"])
            dev = assess_developability(c["sequence"])
            selr = assess_selectivity(c["shuttle"])   # off-target은 셔틀이 주도
            # 전달 축 분해 — deepB3P를 '셔틀 내재 × 융합 보존 × 메커니즘 타당도'로 나눈다.
            si = shuttle_bbb.get(c["shuttle"], 0.0)
            dax = assess_delivery(si, p.bbb_permeability, target=bind.target, mechanism=bind.mechanism)
            rows.append({"label": c["label"], "linker": c["linker"], "shuttle": c["shuttle"],
                         "sequence": c["sequence"], "bbb": round(p.bbb_permeability, 4),
                         "tox": round(t.risk, 4), "toxic": toxic,
                         "bind_ref": bind.best_ref, "bind_score": bind.score,
                         "instability": stab.instability_index, "stable": stab.stable,
                         "dev_risk": dev.risk_level, "dev_liab": dev.n_liabilities,
                         "dev_charge": dev.net_charge, "dev_agg": dev.agg_score,
                         "dev_liabilities": dev.liabilities,
                         "sel_off": selr.off_target_risk, "selectivity": selr.selectivity,
                         "sel_level": selr.risk_level, "sel_mech": selr.mechanism,
                         "sel_drivers": selr.drivers,
                         "sol_score": dev.sol_score, "sol_level": dev.sol_level,
                         "shuttle_bbb": round(si, 4), "mechanism": dax.mechanism,
                         "is_rmt": dax.is_rmt, "preservation": dax.preservation,
                         "avidity": dax.avidity, "deepb3p_valid": dax.deepb3p_validity,
                         "delivery_basis": dax.basis})
            _pres = f"{dax.preservation:.2f}" if dax.preservation is not None else "n/a"
            _mech = dax.mechanism + ("(비교군·최종선택 불가)" if dax.mechanism == "CPP" else "")
            lines.append(f"{c['label']} | {p.bbb_permeability:.3f} | {si:.3f} | {_mech} | "
                         f"{_pres} | {t.risk:.3f} | {stab.instability_index} | {bind.score:.2f} | "
                         f"{'FAIL(penalty)' if toxic else 'ok'} | {c['sequence']}")
        return ("\n".join(lines)
                + f"\n(primary: deepB3P 융합점수 0-1 — **RMT 셔틀엔 약한 프록시**(수용체 결합이 병목) | "
                  "shuttleBBB: 셔틀 단독 deepB3P | mech: RMT(수용체매개)/CPP(막투과) | "
                  f"preserv: 융합 보존도(융합/셔틀) | penalty: ≤{thr:.2f} else FAIL | "
                  "index: lower better, <40 good | match: 수용체유사도 0-1 | "
                  "**CPP는 비특이 흡수 비교군(baseline)이라 최종 선택 대상이 아니다 — "
                  "RMT 후보 중에서 고를 것** | "
                  "→ RMT는 primary 절대값보다 mech·preserv로 판단)" + reject_note), rows

    def _structure(self, cargo, tool_input):
        from .structure import analyze_construct
        lk = "".join(ch for ch in (tool_input.get("prefix", tool_input.get("linker", "")) or "").upper() if ch.isalpha())
        sh = "".join(ch for ch in (tool_input.get("suffix", tool_input.get("shuttle", "")) or "").upper() if ch.isalpha())
        sr = analyze_construct(cargo, lk, sh)
        if sr.error:
            return f"inspect failed: {sr.error}", {"error": sr.error}
        text = (f"exposure={sr.shuttle_exposure} (>=0.35 = exposed), "
                f"confidence(suffix/all)={sr.shuttle_plddt}/{sr.mean_plddt}. "
                f"verdict: {'exposed' if sr.exposed else 'buried/low-confidence'}")
        return text, {"linker": lk, "shuttle": sh, "exposure": sr.shuttle_exposure,
                      "shuttle_plddt": sr.shuttle_plddt, "mean_plddt": sr.mean_plddt,
                      "verdict": sr.verdict, "exposed": sr.exposed}


    # ---- 도달 가능 집합 (심사자에게 줄 기준점) ---------------------------
    def reachable_front(self, cargo: str, refresh: bool = False) -> dict:
        """라이브러리 전수 조합을 **에이전트와 동일한 지표 정의**로 채점해 '도달 가능 집합'과
        그 Pareto 전선을 반환한다. LLM을 쓰지 않는 로컬 계산이며 디스크에 캐시한다.

        왜 필요한가(실측 근거): 심사(Critic)가 **절대적 이상**과 비교하면 승인 가능 영역이
        비어버린다. 화물 GSNKGAIIGLM 실측에서 심사 기준을 AND로 묶으면 130조합 중 0개가
        남았다 — BBB↔off_target +0.615, BBB↔불안정성 +0.496, off↔불안정성 +0.670으로
        축들이 서로 양의 상관이라 동시 만족이 구조적으로 불가능하기 때문이다. 그래서 심사자에게
        **"이 문제에서 실제로 달성 가능한 최선"**을 함께 넘겨, '이상 대비'가 아니라
        **'달성 가능 최선 대비'**로 판정하게 한다.

        전선 축: BBB↑ · off_target↓ · 불안정성↓ (하드 게이트 통과 후보 중 비지배).
        하드 게이트는 **독성 하나뿐**이다. BBB 하한(MIN_BBB_SIGNAL)은 2026-11-10 폐지했다 —
        근거 등급 A 셔틀(Angiopep-2)을 13/13 전멸시키고 등급 D(Tet1)를 통과시키는
        양이온성 필터였다(사유: `core/config.py` 의 폐지 주석). 따라서 전선은 **비독성 전체**
        에서 만들어지고, near-zero 설계체는 순서 특이성 게이트가 막는다.
        """
        import hashlib
        import json

        from .config import BASE_DIR, LINKER_LIBRARY, OFF_TARGET_PENALTY, SHUTTLES

        sig = hashlib.md5(json.dumps({
            "cargo": cargo,
            "linkers": {k: v["seq"] for k, v in LINKER_LIBRARY.items()},
            "shuttles": {k: v["seq"] for k, v in SHUTTLES.items()},
            # 캐시 키: 하한 폐지로 게이트 구성이 바뀌었으므로 "floor" 를 버전 태그로 둔다
            # (옛 캐시를 재사용하면 폐지 전 전선이 섞인다).
            "floor": "no-bbb-floor-v2", "thr": self.settings.toxicity_threshold,
        }, sort_keys=True, ensure_ascii=False).encode()).hexdigest()[:12]
        cache = BASE_DIR / ".cache" / "reachable" / f"{sig}.json"
        if cache.exists() and not refresh:
            try:
                return json.loads(cache.read_text())
            except Exception:  # noqa: BLE001 - 캐시 손상 시 재계산
                pass

        combos = [{"label": f"{ln}+{sn}", "linker": LINKER_LIBRARY[ln]["seq"],
                   "shuttle": SHUTTLES[sn]["seq"]}
                  for ln in LINKER_LIBRARY for sn in SHUTTLES]
        rows = []
        for i in range(0, len(combos), MAX_CANDIDATES_PER_STEP):
            _, chunk = self._evaluate(cargo, {"constructs": combos[i:i + MAX_CANDIDATES_PER_STEP]})
            rows.extend(chunk)

        nontoxic = [r for r in rows if not r.get("toxic")]
        # BBB 하한 폐지 → 하드 게이트는 독성뿐이므로 게이트 통과 = 비독성 전체.
        gate_all = list(nontoxic)
        # 전선은 **후보군(CPP 제외)** 안에서만 만든다 — CPP는 최종 선택 대상이 아니므로
        # 전선에 끼면 "비교군에 지배당했다"는 잘못된 기각 사유를 만들 수 있다.
        gate = [r for r in gate_all if is_candidate(r)]
        gate_cpp = sorted([r for r in gate_all if is_comparison_only(r)],
                          key=lambda r: -(r.get("bbb") or 0.0))

        # 불안정성(Guruprasad instability index)은 문헌에서 **40 기준 이분형**으로 해석한다
        # (<40 안정 / ≥40 불안정). 45.7과 52.2를 다른 값으로 취급하는 건 지표 과독해이고,
        # 그렇게 두면 **효능이 형편없는데 불안정성만 조금 낮은 후보가 '비지배'로 전선에 끼어**
        # 하한을 무너뜨린다(실측: 유효점수 0.098·II 21.7짜리가 비지배로 승인된 사례).
        # 그래서 전선 계산에서는 0(<40) / 1(≥40) 두 값만 쓴다.
        _II_STABLE = 40.0

        def axes(r):
            ii = float(r.get("instability") or 999.0)
            return (r.get("bbb") or 0.0, r.get("sel_off") or 0.0,
                    0 if ii < _II_STABLE else 1)

        def is_dominated(r, pool):
            b1, o1, i1 = axes(r)
            for q in pool:
                if q is r:
                    continue
                b2, o2, i2 = axes(q)
                if b2 >= b1 and o2 <= o1 and i2 <= i1 and (b2 > b1 or o2 < o1 or i2 < i1):
                    return True
            return False

        front = [r for r in gate if not is_dominated(r, gate)]

        def brief(r):
            return {"label": r.get("label"), "bbb": round(r.get("bbb") or 0.0, 3),
                    "off": round(r.get("sel_off") or 0.0, 2),
                    "instability": r.get("instability"), "mechanism": r.get("mechanism"),
                    "dev_risk": r.get("dev_risk"), "shuttle": r.get("shuttle")}

        def is_rmt(r):
            return str(r.get("mechanism") or "").upper().startswith("RMT")

        # CPP(TAT·SynB1·Penetratin 등)는 **비특이 흡수 비교군**이고 RMT가 후보군이다.
        # 게이트를 통과한 RMT가 하나라도 있으면 최종 후보가 CPP일 이유가 없으므로,
        # 심사가 그 대안을 지목해 REVISE할 수 있도록 목록을 함께 넘긴다.
        gate_rmt = sorted([r for r in gate if is_rmt(r)], key=lambda r: -(r.get("bbb") or 0.0))

        # 전선의 **유효점수 범위** — 심사의 '달성 가능 최선 대비' 판정 기준점이 된다.
        # front_eff_min 아래로 떨어진 후보는 비지배여도(사소한 축에서만 앞서도) 결정축에서
        # 전선 밖으로 이탈한 것이다. 임의 상수가 아니라 **전선 자체에서 도출된 값**이다.
        f_effs = [(r.get("bbb") or 0.0) - OFF_TARGET_PENALTY * (r.get("sel_off") or 0.0)
                  for r in front]
        out = {
            "cargo": cargo, "bbb_floor": None,  # 폐지됨 — 심사 텍스트에서 "독성만"으로 표시
            "front_eff_min": round(min(f_effs), 4) if f_effs else None,
            "front_eff_max": round(max(f_effs), 4) if f_effs else None,
            "n_total": len(rows), "n_nontoxic": len(nontoxic),
            "n_gate": len(gate), "n_gate_all": len(gate_all),
            "n_gate_rmt": len(gate_rmt),
            "rmt_options": [brief(r) for r in gate_rmt[:6]],
            # 비교군(CPP) — 참고용 상한. 최종 후보로 제시하지 않는다.
            "cpp_baseline": [brief(r) for r in gate_cpp[:3]],
            "front": [brief(r) for r in sorted(front, key=lambda r: -(r.get("bbb") or 0.0))],
            "best": {
                "bbb_max": round(max((r.get("bbb") or 0.0) for r in gate), 3) if gate else None,
                "off_min": round(min((r.get("sel_off") or 0.0) for r in gate), 2) if gate else None,
                "instability_min": min((r.get("instability") or 999) for r in gate) if gate else None,
            },
        }
        try:
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps(out, ensure_ascii=False, indent=1))
        except Exception:  # noqa: BLE001 - 캐시 실패는 기능에 영향 없음
            pass
        return out
