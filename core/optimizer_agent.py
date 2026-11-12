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
from .evidence import inherit as evidence_inherit
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

AA20 = frozenset("ACDEFGHIKLMNPQRSTVWY")

# 생성 링커의 **길이 상한** = 라이브러리 최장 링커와 같게 둔다(A(EAAAK)3A = 17aa).
# 이것이 유일한 하드 제약이다.
#
# 🔴 **근거를 2026-11-12 에 바꿨다 (§44).** 옛 근거는 **"독성 게이트의 길이 희석"**이었다 —
# *"링커를 늘리면 구성체 독성이 **단조 감소**한다(GGS 0.27 → (GGGGS)3 0.20). 길이로 희석하면
# 독성 셔틀이 게이트를 우회할 수 있다."* **대응 널로 심사했더니 그 근거가 무너졌다:**
#   · "단조"가 **거짓**이다 — 26점 계열 **56개 중 단조 비증가 0개**. 인용된 근거는
#     (화물 1, 셔틀 1, **2점**)이고 2점으로는 단조성을 잴 수 없다. 세 점만 찍어도 깨진다
#     (L=0 **0.245** → GGS(3aa) **0.270** → (GGGGS)3(15aa) 0.205 — 첫 걸음이 **상승**).
#   · 방향은 길이가 아니라 **조성**이 정한다. GS·polyG 는 ρ≈−0.7 로 내려가지만
#     **polyK 는 ρ=+0.890 으로 올라간다**(ToxinPred3 Model 1 은 AAC+DPC 조성 모델이라
#     링커를 붙이면 구성체 조성이 **링커 쪽으로 평균화**된다 — 길이는 혼합 가중치다).
#   · 🔴 **상한을 정하는 어휘가 희석되지 않는다** — `A(EAAAK)3A` 의 EAAAK 계열은
#     ρ 중위 **+0.220**, 음수 4/14. **상한의 *길이*는 EAAAK 에서 오고 근거는 GS 에서 왔다.**
#   · 상한이 통과를 막는 쌍 **0/14**. 겨냥한 우회는 **상한 안에서 전부 완결**된다
#     (Aβ+B6: L=13 에서 0.445 탈락 → 0.115 통과). 상한을 25 로 풀어도 판정 변화 **0/14**
#     (희석이 L≈12~15 에서 포화). 그리고 **17 에서 꺾이는 양이 하나도 없다.**
#   · 반대 방향은 전혀 막지 못한다 — polyK 링커는 **L=3** 에서 통과 구성체를 탈락시키고,
#     **등급 A 셔틀 Angiopep-2 조차 0.265 → 0.440 으로 올린다**(8/56 계열).
#
# ✅ **남은 근거는 하나뿐이고, 그것으로 충분하다 — 대응 널의 범위 선택.**
# 라이브러리와 **같은 길이 운동장**을 쓰면 생성 링커를 라이브러리 링커와 **같은 조건에서**
# 겨루게 할 수 있다(§34 가 그 운동장 위에서 돌아간다). 길이가 늘면 축이 나빠지므로
# (§34-8: `dev_liab` ρ +0.576 · 응집 +0.368, 음수 셀 0/32) **널이 길이에 따라 쉬워지고**,
# 길이를 풀면 "생성이 이겼다"가 길이 조성의 함수가 된다. 이것은 **독성과 무관한 논거**다.
# 🔴 따라서 이 상한을 **"독성 희석"으로 설명하면 안 된다**(§44-7).
# 상한의 값 자체도 측정값이 아니다 — 아래 `max(...)` 가 `A(EAAAK)3A` 의 길이를 집는 것이고,
# 라이브러리에 더 긴 링커를 넣으면 상한이 따라 움직인다.
#
# ⚠️ 다른 축들은 길이 해킹에 이미 면역이거나 자기 방어가 있다(실측 근거):
#   · off-target  — 링커에 **완전 불변**(T7 13/13 모두 0.140 · Angiopep-2 13/13 모두 0.090). 셔틀 주도.
#   · 응집 핫스팟 — 국소 최대라 저응집 잔기를 덧붙여 낮출 수 없다(폭 0.59~0.60).
#   · liability   — **개수**라서 길이로 희석되지 않고, 늘리면 오히려 증가한다
#                   (WGGGGS 가 T7 에서 liab 1→2 = W 산화 부위를 링커가 추가한 실측).
#   · 불안정성    — 강직 반복을 늘리면 계속 좋아지지만(EAAAK 22.9 → (EAAAK)2 20.6 →
#                   A(EAAAK)3A 19.6) **이진 축**이라 임계 40 아래로 들어간 뒤의 이득이 0 이다.
#                   이진화가 여기서 길이 해킹 방어로 작동한다(연속화하면 이 방어가 사라진다).
LINKER_MAX_LEN = max((len(v["seq"] or "") for v in LINKER_LIBRARY.values()), default=17)


def linker_reject_reason(linker: str) -> str | None:
    """생성 링커가 거부되는 이유. 통과면 None.

    라이브러리 **밖 링커를 허용**한다(2026-11-11). 옛 차단의 근거는 §13 의 BBB 접합부
    아티팩트("아무 잔기나 하나 끼우면 융합 BBB 0.526 → 평균 0.743")였고, BBB 를 순위에서
    뺀 뒤(M6·M8) 그 동기가 사라졌다 — 점수를 올릴 수 없는 축을 긁을 이유가 없다.
    남은 제약은 **길이와 알파벳**뿐이다(위 LINKER_MAX_LEN 주석의 실측 근거 참고).
    """
    lk = (linker or "").upper()
    if not lk:
        return None                                  # 직접융합(빈 문자열)은 항상 허용
    bad = sorted(set(lk) - AA20)
    if bad:
        return f"표준 20 아미노산이 아닌 문자 포함: {''.join(bad)}"
    if len(lk) > LINKER_MAX_LEN:
        return (f"길이 {len(lk)}aa > 상한 {LINKER_MAX_LEN}aa. 상한은 라이브러리 최장 링커와 "
                f"같게 뒀다 — 링커를 늘리면 구성체 독성이 희석되어(실측 0.27→0.20) "
                f"**독성 하드 게이트를 우회**할 수 있다")
    return None


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


# ── 문헌 근거 그룹: 양(陽)의 신호의 **유일한** 출처 ──────────────────────────
# 왜 이것뿐인가(실측): 신뢰하는 나머지 축(독성·응집·liability·불안정성·off-target)은 전부
# **"나쁘지 않음"을 재는 위험 스크린**이다. 그래서 불활성 쓰레기가 모든 축에서 이긴다 —
#   poly-Glu : off 0.00 · 선택성 1.00 · 응집 0.12 · liability 0
#   Angiopep-2: off 0.09 · 선택성 0.91 · 응집 0.44 · liability 2
# 불안정성으로도 못 가른다(308 vs 71.5 — 둘 다 >40 이라 어떤 임계값도 등급 A 를 같이 죽인다).
# 유일하게 "좋음"을 재던 축은 BBB 였는데 M6(셔틀 간)·M8(링커 간)이 둘 다 반증했다.
# 따라서 merit 는 **문헌 인용**에서만 온다 — 서열 편집으로 올릴 수 없는 양이다.
def is_unsupported(row: dict) -> bool:
    """문헌 근거가 **없어서** 최종 후보가 될 수 없는가 (등급 D = 근거 없음, 또는 미등재)."""
    from .evidence import group_rank
    return group_rank(row.get("tier")) > 1   # 0 근거 있음 · 1 주의 · 2 근거 없음 · 3 미등재


def is_candidate(row: dict) -> bool:
    """최종 후보가 될 수 있는가 — 독성 게이트 통과 + 비교군(CPP) 아님 + 문헌 근거 있음.

    🔴 `not row["toxic"]` 은 **"비독성"이 아니다**(§47). 그 플래그는
    `ToxinPred3 ML Score > 0.38` 하나이고, 외부 라벨 16개 프로브에서 특이도는 **6/10**
    이었다(위양성: 옥시토신 0.620 · 바소프레신 0.645 · 서브스탠스 P 1.000 ·
    브라디키닌 1.000 — 임상 약과 사람의 내인성 신경펩타이드, 즉 우리 후보가 속한 클래스).
    독소 조각 검출도 같은 길이 무작위 대조 초과 **+9.2pp(Fisher 단측 p = 0.243)** 로
    유의하지 않았다. 그러므로 이 항이 재는 것은 독성 판정이 아니라 **조성 필터** 통과다.

    등급 **D(근거 없음)·미등재**는 CPP 와 **같은 취급**이다: 라이브러리·전선에는 남겨
    비교 근거로 쓰되 최종 후보로는 올리지 않는다. Tet1(등급 D)이 M2·M4 에서 우리 효능 축의
    1위였다는 것이 바로 이 장치가 필요한 이유다 — BBB 통과 단독 근거가 없는 셔틀이
    점수로는 최고였다.
    """
    return (not row.get("toxic")) and not is_comparison_only(row) and not is_unsupported(row)


# 불안정성(Guruprasad instability index)의 문헌 해석 기준 — **40 기준 이분형**.
# 45.7과 52.2를 다른 값으로 취급하는 건 지표 과독해다(연속값으로 두면 효능이 형편없는데
# 불안정성만 조금 낮은 후보가 '비지배'로 전선에 끼어 하한을 무너뜨린다 — 실측 사례 있음).
II_STABLE = 40.0


def risk_axes(r: dict) -> tuple:
    """위험 스크린 축 (전부 **작을수록 좋음**). 축들 사이에 우선순위는 **없다**.

    ⚠️ 응집 축은 **설계 모듈(링커+셔틀)의 AGGRESCAN 핫스팟**이다. 두 번 좁혔고 둘 다 실측 근거가 있다:
      ① 전역 평균 → 국소 최대(핫스팟). 전역 평균은 **길이로 희석**된다 — M10 실측에서
         후보군 셔틀 8종 전부 링커 길이↔응집 순위상관이 음수(평균 -0.697)였고, 그 결과
         전선이 **가장 긴 링커로 전원 수렴**했다(M8 의 BBB 길이 아티팩트와 구조가 같다).
      ② 전체 구성체 → **설계 모듈**. 화물은 고정이고 이 화물(GSNKGAIIGLM)은 `AIIGLM`
         때문에 자체 핫스팟이 0.80 이다. 전체 구성체로 재면 12종 셔틀 **전부 0.80** 으로
         묶여 축이 죽는다(실측). 셔틀 단독은 0.26~0.80 으로 잘 갈린다.
         → `bbb_scoring_seq`/`JUNCTION_FLANK=0` 과 **같은 범위 선택**이다: 우리가 바꿀 수
           있는 모듈을 잰다. 화물 자체의 응집 위험은 **보고**하고 한계로 명시한다.
    핫스팟은 국소 최대이므로 저응집 잔기를 덧붙여 낮출 수 없다 — 길이 해킹에 면역이다.
    """
    ii = float(r.get("instability") or 999.0)
    return (int(r.get("dev_liab") or 0),            # 개발성 liability 개수(전체 구성체)
            0 if ii < II_STABLE else 1,             # 불안정성(이진, 전체 구성체)
            float(r.get("sel_off") or 0.0),         # off-target 위험(셔틀 주도)
            float(r.get("dev_agg_hot") or 0.0))     # 응집 핫스팟(설계 모듈)


def dominates(q: dict, r: dict) -> bool:
    """q 가 r 을 **지배**하나 — 근거 그룹 + 위험 축 전부에서 같거나 낫고 하나는 엄격히 낫다.

    ⚠️ BBB 점수는 축에 **들어가지 않는다**(M6·M8). 여기 있는 축은 전부 '나쁘지 않음'을
    재므로, 근거 그룹이 없으면 불활성 쓰레기가 전부를 지배한다 — 그래서 그룹이 첫 축이다.
    """
    from .evidence import group_rank
    a = (group_rank(q.get("tier")),) + risk_axes(q)
    b = (group_rank(r.get("tier")),) + risk_axes(r)
    return all(x <= y for x, y in zip(a, b)) and a != b


def rank_key(r: dict) -> tuple:
    """표시·대체후보 선택용 **결정론 정렬 키** (작을수록 좋음).

    ⚠️ 위험 축 사이의 순서는 **표시 관례이고 주장이 아니다** — liability 가 off-target 보다
    중요하다는 근거는 없다. 그래서 탈락·게이트 판정에는 이 키를 쓰지 말고 `dominates` 를
    쓴다. 이 키는 "어떤 하나를 결정론적으로 골라야 할 때"(러닝 베스트·대체 제시)만 쓴다.
    """
    from .evidence import group_rank
    return (group_rank(r.get("tier")),) + risk_axes(r)


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
        rows, lines = [], ["label | 근거 | liab | index | off | agg | mech | BBB(참고) | penalty | status | assembly"]
        for c, p, t in zip(clean, bbb, tox):
            # 🔴 `>` 와 `>=` 가 갈리는 자리다(§47-4). ToxinPred3 ML Score 는 **0.005 간격으로
            # 양자화**돼 있고 임계 0.38 은 **하필 그 눈금 위**에 있다 → 점수가 정확히 0.380 인
            # 서열은 `>` 에서 통과하고 `>=` 에서 탈락한다. 실측 영향: 라이브러리 전수 312셀 중
            # **3셀**, 외부 독소 1종(Bombolitin-I, 0.380 → 현재 "비독성"으로 통과).
            # `off_target` 의 `+1e-9`(§42-4b)와 **같은 종류의 경계 정착**이다.
            # 연산자 결정은 동결 해제 후 — 바꾸면 점수 수치가 변해 M1~M10 재측정이 필요하다.
            toxic = t.risk > thr
            bind = shuttle_similarity(c["shuttle"])
            stab = assess_stability(c["sequence"])
            dev = assess_developability(c["sequence"])
            # 응집 순위 축 — **설계 모듈(링커+셔틀)**의 핫스팟. 이유는 `risk_axes` 주석.
            _mod = (c["linker"] or "") + (c["shuttle"] or "")
            dev_mod = assess_developability(_mod) if _mod else dev
            selr = assess_selectivity(c["shuttle"])   # off-target은 셔틀이 주도
            # 문헌 근거 등급 — **merit 의 유일한 출처**. 라이브러리 밖 설계체는
            # 가장 닮은 라이브러리 셔틀에서 **상속**받는다(유사도 하한은 M9 실측).
            _tier, _tparent, _tsim = evidence_inherit(c["shuttle"])
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
                         # 순위 축(모듈 핫스팟) + 보고치(전체 구성체 핫스팟)
                         "dev_agg_hot": dev_mod.agg_hotspot,
                         "dev_agg_hot_full": dev.agg_hotspot,
                         "dev_liabilities": dev.liabilities,
                         "sel_off": selr.off_target_risk, "selectivity": selr.selectivity,
                         "sel_level": selr.risk_level, "sel_mech": selr.mechanism,
                         "sel_drivers": selr.drivers,
                         "tier": _tier, "tier_parent": _tparent, "tier_sim": _tsim,
                         "sol_score": dev.sol_score, "sol_level": dev.sol_level,
                         "shuttle_bbb": round(si, 4), "mechanism": dax.mechanism,
                         "is_rmt": dax.is_rmt, "preservation": dax.preservation,
                         "avidity": dax.avidity, "deepb3p_valid": dax.deepb3p_validity,
                         "delivery_basis": dax.basis})
            _pres = f"{dax.preservation:.2f}" if dax.preservation is not None else "n/a"
            _mech = dax.mechanism + ("(비교군·최종선택 불가)" if dax.mechanism == "CPP" else "")
            _ev = (f"{_tier}" + ("" if _tsim >= 1.0 else f"←{_tparent}·{_tsim:.2f}")) \
                if _tier else (f"미등재(최근접 {_tparent}·{_tsim:.2f})" if _tparent else "미등재")
            lines.append(f"{c['label']} | {_ev} | {dev.n_liabilities} | {stab.instability_index} | "
                         f"{selr.off_target_risk:.2f} | {dev_mod.agg_hotspot:.2f} | {_mech} | "
                         f"{p.bbb_permeability:.3f} | {t.risk:.3f} | "
                         f"{'FAIL(penalty)' if toxic else 'ok'} | {c['sequence']}")
        return ("\n".join(lines)
                + "\n(**근거**: 문헌 전달 근거 등급 — 순위의 유일한 양(陽)의 신호. "
                  "A·B=근거 있음 / C=주의(한계 명시 필수) / D=근거 없음(최종 후보 불가) / "
                  "미등재=라이브러리 밖이고 상속 실패(최종 후보 불가). "
                  "`A←Angiopep-2·0.95` = 그 셔틀에서 상속받음(유사도). | "
                  "liab: 개발성 liability 개수 ↓ | index: 불안정성 <40 good | "
                  "off: off-target 위험 ↓ | agg: **설계 모듈(링커+셔틀)의 응집 핫스팟** ↓ "
                  "— 전역 평균이 아니다(길이로 희석되므로). 화물 자체 핫스팟은 별도 보고 | "
                  "mech: RMT/CPP | "
                  f"penalty: 독성 ≤{thr:.2f} else FAIL)\n"
                + "⚠️ **BBB(참고)는 순위에 쓰지 않는다** — deepB3P 는 셔틀 간(M6)·링커 간(M8) "
                  "순위를 매길 수 없음이 실측으로 확인됐다(M8: 라이브러리 링커가 같은 길이 "
                  "무작위 링커를 이기는 비율 0/108). **보고치로만** 쓰고, '점수가 높다'를 "
                  "선택 근거로 삼지 마라.\n"
                + "**CPP 는 비특이 흡수 비교군이라 최종 선택 대상이 아니다. 등급 D·미등재도 "
                  "같은 취급이다** — 후보는 근거 등급 A·B·C 의 RMT 중에서 고를 것."
                + reject_note), rows

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
        남았다. 그래서 심사자에게 **"이 문제에서 실제로 달성 가능한 최선"**을 함께 넘겨,
        '이상 대비'가 아니라 **'달성 가능 최선 대비'**로 판정하게 한다.

        ── 2026-11-10 개편: 결정축에서 BBB 를 빼고 **문헌 근거 등급**을 1순위로 ──
        옛 전선 축은 `BBB↑ · off↓ · 불안정성↓` 이고 유효점수 `eff = BBB − 0.25·off` 로
        하한을 만들었다. 둘 다 폐지했다 — BBB 는 셔틀 간(M6)·링커 간(M8) 순위를 매길 수
        없음이 실측으로 확인됐기 때문이다(M8: 라이브러리 링커가 같은 길이 무작위 링커를
        이기는 비율 **0/108**, 링커 길이↔BBB 순위상관 평균 **-0.495** = '링커 선택'이
        사실은 '짧게 하기'였다).

        새 구조는 **사전식**이다:
          ① 하드 게이트 — 독성뿐.
          ② 후보 자격   — CPP(비특이 흡수 비교군) 제외 + **등급 D·미등재 제외**.
          ③ 전선        — 근거 그룹 + 위험 축(liability·불안정성·off·응집) Pareto 비지배.
          ④ 최종 후보는 **달성 가능한 최상위 근거 그룹**에서 나와야 한다(`front_top`).
        위험 축끼리는 우선순위를 두지 않는다(`dominates` 주석 참고).
        """
        import hashlib
        import json

        from .config import BASE_DIR, LINKER_LIBRARY, SHUTTLES
        from .evidence import GROUP_ORDER, group_of, group_rank

        sig = hashlib.md5(json.dumps({
            "cargo": cargo,
            "linkers": {k: v["seq"] for k, v in LINKER_LIBRARY.items()},
            "shuttles": {k: v["seq"] for k, v in SHUTTLES.items()},
            # 캐시 키 버전 태그 — 전선의 **축 정의가 바뀌면** 반드시 올린다.
            # (옛 캐시를 재사용하면 BBB 기반 전선이 섞인다.)
            "axes": "evidence-group+risk-pareto-v3", "thr": self.settings.toxicity_threshold,
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
        # 하드 게이트는 독성 게이트 하나뿐 → 게이트 통과 = 그 게이트를 통과한 전체.
        # 🔴 "비독성 전체"로 읽으면 안 된다(§47) — 게이트는 `ML Score > 0.38` 이라는
        # **조성 필터**이고, 외부 라벨 특이도 6/10 · 독소 조각 검출 8/24(무작위 29/120,
        # 초과 +9.2pp, p = 0.243)다. 아래 `nontoxic` 이라는 이름도 같은 이유로 과하다
        # (이름 변경은 동결 해제 후 — `docs/experiments.md` §37-4 대기열).
        gate_all = list(nontoxic)
        # 후보군 = 비교군(CPP) 아님 + 문헌 근거 있음(등급 A·B·C). 전선은 이 안에서만 만든다 —
        # 비교군·근거 없음이 전선에 끼면 "지배당했다"는 잘못된 기각 사유를 만든다.
        gate = [r for r in gate_all if is_candidate(r)]
        gate_cpp = sorted([r for r in gate_all if is_comparison_only(r)], key=rank_key)
        # 등급 D·미등재 — CPP 와 같은 취급(비교 근거로만). Tet1 이 여기 온다.
        gate_unsup = sorted([r for r in gate_all
                             if is_unsupported(r) and not is_comparison_only(r)], key=rank_key)

        front = sorted([r for r in gate
                        if not any(dominates(q, r) for q in gate if q is not r)], key=rank_key)

        # ④ 최종 후보가 나와야 하는 집합 — **달성 가능한 최상위 근거 그룹**의 전선.
        # 왜 그룹을 먼저 고정하나: 근거 등급은 우리 점수가 아니라 인용에서 오는 유일한 양의
        # 신호다. '주의(C)' 후보가 위험 축이 좋다는 이유로 '근거 있음(A·B)'을 제치면,
        # 위험 축은 전부 "나쁘지 않음"을 재므로 **불활성 쪽이 이기는 구조**로 되돌아간다.
        best_g = min((group_rank(r.get("tier")) for r in gate), default=None)
        front_top = [r for r in front if group_rank(r.get("tier")) == best_g] if gate else []

        def brief(r):
            return {"label": r.get("label"), "tier": r.get("tier"),
                    "group": group_of(r.get("tier")),
                    # ⚠️ 키 이름을 행(row)과 **같게** 유지한다 — 심사의 지배 게이트가
                    # 전선 요약을 `dominates()` 에 그대로 넣기 때문이다(어댑터 금지).
                    "dev_liab": r.get("dev_liab"), "instability": r.get("instability"),
                    "sel_off": round(r.get("sel_off") or 0.0, 2),
                    "dev_agg_hot": r.get("dev_agg_hot"), "mechanism": r.get("mechanism"),
                    "dev_risk": r.get("dev_risk"), "shuttle": r.get("shuttle"),
                    # BBB 는 **보고치**로만 실어 보낸다(순위·게이트에 쓰지 않는다).
                    "bbb_report": round(r.get("bbb") or 0.0, 3)}

        def is_rmt(r):
            return str(r.get("mechanism") or "").upper().startswith("RMT")

        gate_rmt = sorted([r for r in gate if is_rmt(r)], key=rank_key)

        # 근거 그룹별 수 — 심사자가 "더 나은 근거의 대안이 있나"를 바로 볼 수 있게.
        # ⚠️ 분모는 **비교군(CPP)을 뺀 게이트 통과 전체**다(= 독성 게이트 통과.
        # "비독성"이 아니다 — §47). `gate` 로 세면 `gate` 가 이미 등급 D 를
        # 제외하므로 '근거 없음'이 **항상 0** 으로 표시된다(= 실제로 있는 걸 없다고 보고하는 셈).
        _pool = [r for r in gate_all if not is_comparison_only(r)]
        by_group = {g: len([r for r in _pool if group_of(r.get("tier")) == g])
                    for g in GROUP_ORDER}
        by_group["미등재"] = len([r for r in _pool if r.get("tier") is None])

        out = {
            "cargo": cargo,
            "n_total": len(rows), "n_nontoxic": len(nontoxic),
            "n_gate": len(gate), "n_gate_all": len(gate_all),
            "n_gate_rmt": len(gate_rmt),
            "by_group": by_group,
            "best_group": (GROUP_ORDER[best_g]
                           if best_g is not None and best_g < len(GROUP_ORDER) else None),
            "rmt_options": [brief(r) for r in gate_rmt[:6]],
            # 비교군(CPP) — 참고용 상한. 최종 후보로 제시하지 않는다.
            "cpp_baseline": [brief(r) for r in gate_cpp[:3]],
            # 등급 D·미등재 — 참고용. 최종 후보로 제시하지 않는다.
            "unsupported": [brief(r) for r in gate_unsup[:3]],
            "front": [brief(r) for r in front],
            "front_top": [brief(r) for r in front_top],
            "best": {
                "liab_min": min((r.get("dev_liab") or 0) for r in gate) if gate else None,
                "off_min": round(min((r.get("sel_off") or 0.0) for r in gate), 2) if gate else None,
                "instability_min": min((r.get("instability") or 999) for r in gate) if gate else None,
                "agg_min": min((r.get("dev_agg_hot") or 0.0) for r in gate) if gate else None,
            },
        }
        try:
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps(out, ensure_ascii=False, indent=1))
        except Exception:  # noqa: BLE001 - 캐시 실패는 기능에 영향 없음
            pass
        return out
