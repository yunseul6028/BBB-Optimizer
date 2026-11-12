"""자율 설계 에이전트 — Gemini 브레인 (experiment/gemini 브랜치).

Claude 버전(`optimizer_agent.py`)과 **동일한 도구 백엔드**(_evaluate/_structure/_generate)를
재사용하되, 브레인을 Google Gemini(function calling)로 바꾼 실험 어댑터.

핵심 실험 질문: **Anthropic 안전 분류기가 오탐 거부한 "BBB 통과 구조 최적화"를 Gemini는
그대로(자연스러운 생물학 프레이밍) 받아주는가?** → 여기서는 중립화 없이 **정직한 BBB 프레이밍**
을 그대로 쓴다(도구명·프롬프트 모두 생물학 용어). 거부되면 Claude처럼 중립화가 필요하다는 결론.

이벤트(AgentEvent) 종류·데이터 키는 Claude 버전과 동일 → 앱 UI 그대로 재사용.
⚠️ GEMINI_API_KEY(또는 GOOGLE_API_KEY)가 있을 때만 동작.
"""

from __future__ import annotations

import time
from typing import Generator

from .config import LINKER_LIBRARY, MODEL_MAX_LEN, SHUTTLES, Settings
from .optimizer_agent import (AgentEvent, OptimizationAgent,
                              dominates, is_candidate, is_comparison_only, is_unsupported,
                              rank_key, risk_axes)
from .evidence import INHERIT_SIM_MIN
from .evidence import critic_rule as evidence_critic_rule
from .evidence import group_of as evidence_group
from .evidence import group_rank as evidence_group_rank
from .evidence import inherit as evidence_inherit
from .evidence import prompt_brief as evidence_brief
from .evidence import tier_of as evidence_tier

# 일시적 오류(모델 과부하·분당 레이트리밋·응답 지연) — 백오프 재시도 대상
_TRANSIENT = ("503", "UNAVAILABLE", "429", "RESOURCE_EXHAUSTED", "overloaded", "high demand",
              "504", "DEADLINE_EXCEEDED", "timeout", "timed out")

# Gemini 요청당 최대 대기(ms). 응답이 이 시간 내 안 오면 예외 → 재시도/에러 표시.
# (미설정 시 stall된 연결에서 UI가 무한 "실행 중"에 걸린다.)
_REQUEST_TIMEOUT_MS = 120_000

# ⚰️ 유효점수 `eff = BBB − λ·off_target` 은 2026-11-10 **폐지**했다(λ=OFF_TARGET_PENALTY).
#    그 식의 전제는 "BBB 점수로 후보를 줄 세울 수 있다"였고, 그 전제가 두 번 반증됐다:
#      · M6(셔틀 간) — 문헌 검증 셔틀이 자기 스크램블 분포의 39번째 백분위, 유의 0/9
#      · M8(링커 간) — 라이브러리 링커가 **같은 길이 무작위 링커**를 이긴 비율 **0/108**(0%),
#                     링커 길이↔BBB 순위상관 평균 **-0.495**(= '링커 선택'이 '짧게 하기'였다)
#    대체: 순위는 `optimizer_agent.rank_key` / `dominates` — **문헌 근거 그룹**(merit)을
#    사전식 1순위로 두고, 그 안에서 위험 축(liability·불안정성·off·응집)을 Pareto 로 본다.
#    BBB 는 **보고치**로만 남는다. 되살리지 말 것.


def _gemini_system_prompt(cargo: str, tox_threshold: float, max_rounds: int) -> str:
    lib_l = ", ".join(f"{n}({v['seq'] or '링커없음/직접융합'})" for n, v in LINKER_LIBRARY.items())
    lib_s = ", ".join(f"{n}({v['seq']}, 근거 {evidence_tier(n)}등급)"
                      for n, v in SHUTTLES.items())
    tools = ("evaluate_candidates(**라이브러리 조합만** 채점 — BBB·독성·안정성·수용체), "
             "design_candidate(**라이브러리 밖** 직접 편집·설계 서열 채점), "
             "optimize_developability(셔틀 개발성 다목적 Pareto 최적화 — 기능 보존 잔기 편집)")
    tools += ", finish(수렴 시 최종 후보 제출·종료)"
    return (
        "[연구 맥락] 이것은 알츠하이머병 치료제의 뇌 전달을 개선하기 위한 정당한 학술·공모전용 "
        "in-silico 연구입니다. 사용하는 셔틀 펩타이드(Angiopep-2, TAT, Penetratin, SynB1 등)는 "
        "학계·임상에서 공개 발표된 뇌 전달 도구이며, 독성 예측은 안전한 후보만 남기기 위한 안전성 "
        "스크리닝입니다.\n\n"
        "당신은 멀티 에이전트 신약설계 시스템의 **설계(Designer) 에이전트**입니다 — 도구를 자율적으로 "
        "오케스트레이션해 뇌혈관장벽(BBB)을 통과하는 항체-셔틀 융합 단백질을 설계하고, 그 결과는 "
        "독립 **심사(Critic) 에이전트**의 적대적 검증을 거칩니다.\n\n"
        f"고정 화물(cargo): {cargo}\n융합체 = cargo + linker + shuttle.\n\n"
        "[목적 — 무엇을 최적화하는가. 2026-11-10 개편]\n"
        "  이 시스템은 **BBB 투과도를 예측하지 않는다.** 문헌으로 검증된 셔틀을 붙이면 뇌 전달이 "
        "올라간다는 것은 **문헌에서 빌려오는 전제**이고, 네가 할 일은 그 융합체를 **약으로 만들 수 "
        "있게** 다듬는 것이다. 우리 BBB 점수(deepB3P)로 셔틀이나 링커의 순위를 매길 수 없음이 "
        "실측으로 확인됐다 — 셔틀 간 M6(유의 0/9), 링커 간 M8(라이브러리 링커가 같은 길이 "
        "**무작위** 링커를 이긴 비율 **0/108**, 링커 길이↔BBB 순위상관 평균 **-0.495**).\n"
        "\n"
        "  순위는 **사전식**이다 — 앞 단계가 뒤 단계를 지배한다:\n"
        "  ① **문헌 근거 그룹** (유일한 양(陽)의 신호) — `근거 있음`(등급 A·B) > `주의`(C) > "
        "`근거 없음`(D, **후보 자격 없음**) > `미등재`(**후보 자격 없음**).\n"
        "     최종 후보는 **달성 가능한 최상위 그룹**에서 나와야 한다. 근거가 더 단단한 대안이 "
        "존재하는데 아래 그룹을 고르면 심사가 자동 반려한다.\n"
        "  ② **위험 스크린** (같은 그룹 안에서, 서로 우선순위 **없음** — Pareto):\n"
        f"     · 독성(ToxinPred3) ≤ {tox_threshold:.2f} — 이건 **하드 게이트**다(탈락)\n"
        "     · 개발성 liability 개수 ↓ · 불안정성 지수 <40 · off-target 위험 ↓ · 응집 ↓\n"
        "     이 축들은 전부 **'나쁘지 않음'을 재는 스크린**이다. 좋음을 재지 않는다 — 그래서 "
        "불활성 쓰레기(poly-Glu 등)가 전 축에서 이긴다. 그게 ①이 먼저 오는 이유다.\n"
        "\n"
        "  **[셔틀 역할] CPP(TAT·Penetratin·SynB1 등 양이온성 막투과)는 비특이 흡수 "
        "**비교군(baseline)**이다. 비교로 채점하는 건 좋지만 **최종 후보로 제출하지 마라.** "
        "등급 **D**(Tet1)와 **미등재** 설계체도 같은 취급이다 — 제출하면 자동 반려된다.\n"
        "\n"
        "  **[설계체의 자리 — 중요]** 라이브러리 밖 설계체에는 인용이 없다. 그래서 "
        "**가장 닮은 라이브러리 셔틀의 등급을 상속**받는다 — 서열 유사도가 "
        f"**{INHERIT_SIM_MIN:.2f} 이상**일 때만(실측 간극: 정당한 개발성 변이체 최소 0.830 vs "
        "스크램블 0.529·무작위 0.136). 즉 **설계의 정당한 역할은 효능을 올리는 것이 아니라, "
        "부모의 근거를 유지한 채 liability·응집·불안정성을 줄여 '만들 수 있게' 만드는 것**이다. "
        "조성을 쌓거나 길이를 늘리는 방향은 **아무 효과가 없다** — 그 점수를 순위에 쓰지 않는다.\n\n"
        "[BBB 점수 — 보고치다. 순위에 쓰지 마라]\n"
        "  evaluate_candidates 는 deepB3P 를 분해해 보여준다: 융합 연결부위 점수·셔틀 단독·"
        "mech(RMT/CPP)·preserv(융합 보존도). **이 값들은 최종 보고서의 '참고 수치'이고 선택 "
        "근거가 아니다.** 반증 내역:\n"
        "   · M1 — 양성 vs 무작위 AUC 0.807 인데 양성 vs **스크램블**(조성 동일)은 0.618 "
        "→ 점수의 상당 부분이 **아미노산 조성**에서 온다(= 양이온 잔기를 쌓으면 오른다).\n"
        "   · M6 — 융합 수준에서 문헌 검증 셔틀이 자기 스크램블 분포의 39번째 백분위, 유의 0/9, "
        "문헌 음성 대조군(Angiopep-7) **역전**.\n"
        "   · M8 — 링커조차 못 고른다. 라이브러리 링커가 같은 길이 **무작위** 링커를 이긴 비율 "
        "**0/108**(우연이면 5%), 경험 p 중위 ~0.86 = 설계된 링커가 무작위보다 **나쁘다**. "
        "그리고 길이↔점수 순위상관 -0.495 — 6/9 셔틀에서 '최적 링커'가 **직접융합(링커 없음)**"
        "이었다. 즉 옛 시스템이 '링커 없음'을 권한 것은 **아티팩트**였다.\n"
        "  👉 그래서 **'BBB 점수가 낮다/높다'를 선택 근거로 쓰면 심사가 반박한다.** "
        "낮다고 버리지도, 높다고 고르지도 마라. 최종 보고서에는 보고치로 적고 "
        "'이 셔틀을 고른 근거(등급·인용)'를 따로 밝혀라.\n"
        "  링커는 무엇으로 고르고 **무엇을 설계하나**(실측 2026-11-11, 라이브러리 13종 × 셔틀):\n"
        "   · **off-target·응집 핫스팟은 링커에 반응하지 않는다** — T7 에서 off 가 13/13 모두 "
        "0.140, Angiopep-2 는 13/13 모두 0.090 으로 **완전히 동일**했고, 응집 핫스팟도 폭이 "
        "0.59~0.60 뿐이었다(국소 최대라 잔기를 덧붙여 낮출 수 없다). **링커로 이 두 축을 "
        "고치려 하지 마라 — 움직이지 않는다.**\n"
        "   · **liability 는 링커가 실제로 움직인다** — `WGGGGS` 는 T7 에서 liab 1→2, "
        "Angiopep-2 에서 3→4 로 **올렸다**(W 산화 부위 추가). 화물과 셔틀은 고정이므로 "
        "**접합부**(화물 끝⊕링커⊕셔틀 시작)에서 NG/NS·DG/DP·산화 M/W 모티프가 생기지 않게 "
        "하는 것은 **링커만 할 수 있는 일**이다. 여기가 링커 설계의 정당한 표적이다.\n"
        "   · **불안정성은 조성이 방향을 정한다** — 유연 GS 계열은 길수록 나빠지고"
        "(GGS 30.0 → (GGGGS)3 52.7), 강직 EAAAK 계열은 길수록 좋아진다"
        "(EAAAK 22.9 → A(EAAAK)3A 19.6). 단 이 축은 **이진**(임계 40)이라 40 아래로 "
        "들어간 뒤에는 더 내려도 이득이 0 이다 — 길이를 늘려 숫자만 깎지 마라.\n"
        "   · 그리고 **문헌·표준 관행**(GGGGS = 표준 유연 / EAAAK = 강직 α-헬릭스)과 "
        "**구조**(링커의 실제 기능은 점수를 올리는 게 아니라 **셔틀이 수용체에 닿게 하는 것**).\n\n"
        f"사용 가능한 도구: {tools}\n"
        "권장 워크플로우(자율 판단): evaluate_candidates로 라이브러리 조합(**링커 유무 — 직접융합"
        "(링커 빈칸)도 포함**)을 폭넓게 스크리닝 → "
        "유망 방향으로 **design_candidate로 잔기 수준 편집**(변이·트리밍·연장·하이브리드)해 라이브러리 "
        "밖까지 정밀 설계·재평가 → **충분히 수렴했다고 "
        "판단하면 스스로 finish**를 호출해 최종 융합체 1개(라벨·링커·셔틀·전체 근거)를 제출하라.\n"
        "  **[도구 역할 분리 — 반드시 지킬 것] `evaluate_candidates`는 위 라이브러리 링커·셔틀 "
        "**원본 서열 조합만** 채점한다(선택 도구). 한 잔기라도 편집한 서열을 넣으면 **채점을 거부**하고 "
        "되돌려보낸다. 직접 설계한 서열은 반드시 `design_candidate`로 제출하라 — 가설·예측을 적어야 "
        "한다. 이는 '고르기'와 '만들기'를 기록상 분리해 설계의 기여를 측정하기 위한 것이다.**\n"
        "- design_candidate: 라이브러리에 얽매이지 말고 셔틀/링커의 특정 잔기를 자유롭게 편집한 서열을 "
        "제안·채점할 수 있다(진짜 설계). 왜 그렇게 편집했는지 edit_note에 남겨라.\n"
        "- optimize_developability: 최종 셔틀을 확정하기 전에 **개발성(liability·응집·불안정성)을 신뢰 "
        "오라클로 다목적 Pareto 최적화**하고 싶을 때. 기능(rmt_sim)을 보존하는 변이체 전선을 받아, "
        "개발성이 개선된 셔틀 하나를 골라 design_candidate/finish에 쓰라. 이는 **BBB 효능이 아니라 "
        "'만들 수 있게(제조·안정성)' 다듬는 축** — BBB는 검증 라이브러리 셔틀 선택이 담당한다.\n"
        "- finish: 고정 스텝을 다 쓸 필요 없다. 목적을 충분히 만족했다고 판단하면 언제든 finish로 최종 "
        "후보를 제출하라. 제출하면 **독립 심사(비평) 에이전트가 적대적으로 검증**한다 — 심사가 REVISE(개선 "
        "요구)를 내면 그 지적을 반영해 개선하고 다시 finish하라. APPROVE(승인)가 나오면 확정된다.\n\n"
        "[설계 방법 — 가설 주도 루프 (중요)] design_candidate로 편집할 때마다 **무작위 변이가 아니라 "
        "근거 있는 가설**로 움직여라:\n"
        "  ① 관찰: 직전 평가에서 무엇이 문제인가(어느 축이 약한가 — off_target·BBB·preserv·안정성 등).\n"
        "  ② 가설(hypothesis): 그 문제를 **어떤 편집으로 왜** 고칠 수 있는가 (예: '응집 핫스팟이 "
        "방향족 FFY 연속에서 오니 F3→D로 극성 잔기를 넣으면 핫스팟↓, rmt_sim은 보존').\n"
        "  🔴 **off_target 을 낮추려고 양전하(K·R)를 지우는 편집은 금지다.** 그 축은 대응 널 심사에서 "
        "*악화 검출은 신뢰할 수 있지만 개선 신호는 신뢰할 수 없다*고 판정됐다(무작위 편집의 46%가 "
        "위험 감소 없이 off를 내리고, CPP의 전하만 지우면 TAT가 0.95→0.00으로 '안전'해진다). "
        "게다가 Angiopep-2의 K10·K15는 문헌이 기능 필수라고 지목한 라이신이다(Angiopep-7=K10R/K15R은 "
        "뇌 전달 ~9.5배 낮다). off는 **올리지 않는 것**이 목표이고, 내리는 것은 성과가 아니다.\n"
        "  ③ 예측(prediction): 그 편집이 **어느 축을 어떻게** 바꿀지 (방향·대략 수치).\n"
        "  ④ 오라클 채점(자동) → ⑤ **반성**: 예측이 맞았는지 **실제 점수와 대조**하라. 틀렸으면 어느 "
        "축에서 왜 빗나갔는지 말하고 **가설을 수정해 다음 설계**로 이어가라.\n"
        "**오라클이 네 가설을 반증하면 솔직히 인정하고 방향을 바꿔라 — 점수를 우기거나 지어내지 마라.** "
        "너의 가치는 '좋은 가설을 세우고 신뢰 오라클로 검증·교정하는 추론'에 있다(무작위 탐색이 아님).\n"
        "**[필수] 유망한 후보를 찾으면 finish 전에 반드시 `design_candidate`로 가설 주도 편집을 "
        "최소 1회 시도해 best를 개선하라** — 라이브러리 조합만으로 끝내지 마라(네 손으로 가설을 세워 설계하라).\n\n"
        f"{evidence_brief()}\n\n"
        f"참고 라이브러리 — 링커: {lib_l}\n  셔틀: {lib_s}\n"
        f"주의: construct가 {MODEL_MAX_LEN}aa 초과면 BBB는 연결부위(링커+셔틀)로 계산된다. deepB3P는 "
        f"짧은 펩타이드 학습이라 절대값보다 상대 비교가 신뢰됨. 최대 {max_rounds}행동. 한국어로 간결히."
    )


class GeminiOptimizationAgent(OptimizationAgent):
    """Claude 백엔드 재사용 + Gemini function-calling 루프."""

    def __init__(self, settings, max_rounds: int = 10, use_ledger: bool = True,
                 ledger_write: bool = True):
        super().__init__(settings, max_rounds)
        # 설계 이력 원장(core/ledger.py) 사용 여부.
        # ⚠️ **벤치마크·애블레이션은 반드시 False**로 고정한다. 원장을 켠 세션은 이전
        # 세션의 탐색을 물려받으므로 콜드 스타트 세션과 **비교 불가**다('18평가로 134.8%'
        # 같은 효율 수치의 의미가 달라진다). 섞으면 §12 분모 불일치와 같은 거짓 비교가 된다.
        self.use_ledger = use_ledger
        # 읽기는 하되 **쓰지 않는** 모드. 켬/끔 A/B 실험 전용이다 —
        # 쓰기를 켜두면 ON 팔의 2·3세션이 **서로 다른(자라난) 원장**을 보게 되어
        # 같은 조건의 반복이 아니게 된다. 두 팔이 **같은 바닥**을 보게 동결한다.
        self.ledger_write = ledger_write

    def _function_decls(self, types):
        T = types.Type
        S = types.Schema
        decls = [
            types.FunctionDeclaration(
                name="evaluate_candidates",
                description=("후보 융합체(cargo+linker+shuttle)들의 BBB 투과 점수(deepB3P), 독성"
                             "(ToxinPred3), 융합체 안정성(불안정성 지수), 수용체 유사도를 배치로 채점."),
                parameters=S(type=T.OBJECT, properties={
                    "rationale": S(type=T.STRING, description="이 후보들을 낸 근거"),
                    "constructs": S(type=T.ARRAY, items=S(type=T.OBJECT, properties={
                        "label": S(type=T.STRING),
                        "linker": S(type=T.STRING),
                        "shuttle": S(type=T.STRING),
                    }, required=["label", "linker", "shuttle"])),
                }, required=["rationale", "constructs"]),
            ),
            types.FunctionDeclaration(
                name="design_candidate",
                description=("**직접 편집·설계한 셔틀/링커**를 채점한다. 편집 전 "
                             "**가설(무엇을 왜 바꾸나)·예측(어느 축이 어떻게 변할까)**을 반드시 제시하라 "
                             "— 무작위 변이가 아니라 근거 있는 설계. 셔틀은 잔기 변이·트리밍·연장·"
                             "하이브리드 등 자유 설계이고, **링커도 자유 설계**가 가능하다"
                             "(표준 20 아미노산·최대 17aa). ⚠️ **링커로 무엇을 노릴 수 있는지는 "
                             "축마다 다르다**(실측): off-target 과 응집 핫스팟은 링커에 거의 반응하지 "
                             "않는다(라이브러리 13종에서 off 는 13/13 동일). 링커가 실제로 움직이는 "
                             "축은 **접합부 liability**(화물 끝⊕링커⊕셔틀 시작에서 NG/NS/DG/DP·산화 "
                             "M/W 같은 모티프가 생기거나 사라진다 — 화물·셔틀은 고정이므로 이걸 "
                             "통제할 수 있는 모듈은 링커뿐이다)와 **불안정성**이다. "
                             "evaluate와 동일 지표 반환."),
                parameters=S(type=T.OBJECT, properties={
                    "label": S(type=T.STRING),
                    "linker": S(type=T.STRING, description="자유 설계 가능 — 직접융합은 빈 문자열. "
                               "표준 20 아미노산·최대 17aa(라이브러리 최장과 동일)."),
                    "shuttle": S(type=T.STRING, description="자유 설계 가능(근거 등급은 셔틀에서 상속된다)"),
                    "hypothesis": S(type=T.STRING, description="무엇을 왜 편집하는가 — 직전 평가의 어느 "
                                    "문제(축)를 어떤 편집으로 고치려는지"),
                    "prediction": S(type=T.STRING, description="그 편집이 어느 축을 어떻게 바꿀지 예측(방향·대략치)"),
                    "edit_note": S(type=T.STRING, description="구체적 편집 내용(어느 잔기를 어떻게)"),
                }, required=["label", "linker", "shuttle", "hypothesis", "prediction"]),
            ),
        ]
        # analyze_structure(ESMFold 구조 노출도)는 **에이전트 도구에서 제거**했다.
        # 근거(애블레이션 n=3, docs/experiments.md §11): 이 도구를 뺀 조건이 **자기 고유 축인
        # 셔틀 노출도에서조차 퇴행하지 않았다**(0.689 vs 대조군 0.693). 유효점수는 오히려 높았다
        # (0.577 vs 0.469). 즉 이 도구는 에이전트의 판단을 개선하지 않으면서 호출당 ~10초와
        # 공개 API 장애 위험만 추가하고 있었다. 점수 공식을 바꿔 쓸모를 만드는 대신 도구를 뺀다.
        # 구조 근거 자체는 버리지 않고 **최종 후보 1건에만 best-effort로 접어 보고**한다
        # (아래 run()의 구조 증거 단계) — 에이전트의 의사결정에는 들어가지 않는다.
        decls.append(types.FunctionDeclaration(
            name="optimize_developability",
            description=("한 셔틀 서열을 **순위 위험 4축**(liability·응집·불안정성·**off-target**)에서 "
                         "신뢰 가능한 규칙/척도 오라클로 **다목적(Pareto) 최적화**한다. 검증 기능을 "
                         "보존(rmt_sim 하한) + **문헌 근거 등급 상속 보존**하는 보존적 잔기 치환으로, "
                         "서로 상충하는 네 축을 동시에 줄인 변이체들의 트레이드오프 전선을 반환하고 "
                         "**그 전선을 실제로 채점해 러닝 베스트·래칫에 반영한다**. "
                         "목적함수가 심사의 축과 **같다** — 네가 직접 편집할 때 자주 잃는 off-target 까지 "
                         "함께 본다. **BBB 효능이 아니라 '만들 수 있게(제조·안정성·선택성)' 다듬는 용도** — "
                         "최종 후보 셔틀을 확정하기 전에 호출하라."),
            parameters=S(type=T.OBJECT, properties={
                "shuttle": S(type=T.STRING, description="개발성을 다듬을 셔틀 서열"),
                "linker": S(type=T.STRING, description=(
                    "함께 쓸 링커 서열(선택). liability·불안정성은 구성체 전체, 응집은 "
                    "모듈에서 재므로 링커가 채점 결과를 바꾼다. 생략하면 직접융합으로 채점한다.")),
            }, required=["shuttle"]),
        ))
        decls.append(types.FunctionDeclaration(
            name="finish",
            description=("충분히 수렴했다고 판단하면 호출해 최종 후보를 제출하고 종료한다. 고정 스텝을 "
                         "다 쓸 필요 없다."),
            parameters=S(type=T.OBJECT, properties={
                "chosen_label": S(type=T.STRING),
                "chosen_linker": S(type=T.STRING),
                "chosen_shuttle": S(type=T.STRING),
                "final_report": S(type=T.STRING, description="최종 융합체·4지표·선정 근거"),
            }, required=["chosen_linker", "chosen_shuttle", "final_report"]),
        ))
        return decls

    def _dispatch(self, cargo, name, args):
        """도구 실행 → (LLM에 돌려줄 텍스트, 방출할 AgentEvent 또는 None)."""
        if name == "evaluate_candidates":
            # library_only=True — 평가 도구는 **선택** 전용. 라이브러리 밖은 design_candidate로.
            text, rows = self._evaluate(cargo, args, library_only=True)
            return text, AgentEvent("evaluation", text=args.get("rationale", ""), data={"rows": rows})
        if name == "design_candidate":
            # ── 링커 제약 (2026-11-11: 라이브러리 밖 링커 **허용**) ─────────────
            # 옛 차단의 근거는 §13 의 BBB 접합부 아티팩트였다("화물+T7 에 아무 잔기나 하나
            # 끼우면 융합 BBB 0.526 → 평균 0.743, 20종 중 18종 상승"). 그때는 BBB 가 순위에
            # 있었으니 열어두면 에이전트가 셔틀 설계 대신 접합부를 긁었다(승자 전원 `V+...`).
            # **BBB 를 순위에서 뺀 뒤(M6·M8) 그 동기가 사라졌다** — 점수를 올릴 수 없는 축을
            # 긁을 이유가 없다. 그리고 차단은 애초에 우회 가능했다: `링커 Y + 셔틀 HAIYPRH` 와
            # `직접융합 + 셔틀 YHAIYPRH` 는 **같은 분자**이고, 실제로 에이전트가 Y 를 셔틀 앞으로
            # 옮겨 그대로 통과했다. 링커/셔틀 구분은 장부상의 경계일 뿐 화학적 경계가 아니다.
            #
            # 여는 이유(실측): 링커가 정당하게 기여할 수 있는 자리가 **실재한다** —
            # liability 는 **개수**라서 링커가 바꿀 수 있고(WGGGGS 가 T7 에서 liab 1→2 로
            # W 산화 부위를 추가한 실측), 화물과 셔틀은 고정이므로 **접합부 liability 를
            # 통제할 수 있는 모듈은 링커뿐**이다. 반면 라이브러리 13종은 그 축에서 거의
            # 균일해(T7 에서 서로 다른 위험축 조합 4개/13링커) 선택만으로는 손댈 수 없었다.
            #
            # 남은 제약은 **길이·알파벳뿐**이고 근거는 독성 게이트의 길이 희석이다
            # (`optimizer_agent.linker_reject_reason` 주석에 축별 실측 근거를 적어뒀다).
            # 장부 정규화 목적은 유지된다: 설계는 여전히 `design_candidate` 로만 들어오고
            # (가설·예측이 기록에 남는다), 근거 등급 상속 판정은 **셔틀 필드**를 보므로
            # 링커를 열어도 그 통제는 그대로 걸린다.
            from .optimizer_agent import linker_reject_reason
            lk_raw = (args.get("linker", "") or "")
            lk = "".join(ch for ch in lk_raw.upper() if ch.isalpha())
            _why = linker_reject_reason(lk)
            if _why is not None:
                from .optimizer_agent import LINKER_MAX_LEN
                self.gate_rejections = getattr(self, "gate_rejections", 0) + 1
                msg = (f"⚠️ 링커 `{lk}`를 **채점하지 않았다** — {_why}.\n"
                       f"링커는 자유 설계가 가능하지만 **표준 20 아미노산 · 최대 "
                       f"{LINKER_MAX_LEN}aa** 안에서만 가능하다. 링커로 노릴 수 있는 축은 "
                       f"**접합부 liability**(모티프를 만들지 않는 조성)와 **불안정성**이고, "
                       f"off-target 과 응집 핫스팟은 **라이브러리 링커 어휘에서만** 거의 불변이고, 그것도 **내리는 방향으로만** 참이다(§34-5: 자유 생성에서 무작위 링커의 53~56%가 핫스팟을 **올린다**; 23,040 관측 중 내린 사례는 **0건**). 즉 링커로 응집을 **낮출 수는 없고 높일 수는 있다** — 이 축에서 할 일은 개선이 아니라 **악화 회피**다.")
                return msg, AgentEvent("text", text=f"🚫 링커 `{lk}` 거부 — {_why}")
            text, rows = self._evaluate(cargo, {"constructs": [{
                "label": args.get("label", "설계"), "linker": lk,
                "shuttle": args.get("shuttle", "")}]})
            hyp = args.get("hypothesis", "").strip()
            pred = args.get("prediction", "").strip()
            note = args.get("edit_note", "").strip()
            detail = " · ".join(x for x in (
                f"가설: {hyp}" if hyp else "",
                f"예측: {pred}" if pred else "",
                f"편집: {note}" if note else "") if x)
            # ── 근거 상속 즉시 피드백 ──
            # 설계체가 **최종 후보 자격이 있나**(= 부모 셔틀의 근거를 물려받았나)를 그 자리에서
            # 알려준다. 심사까지 끌고 가면 라운드를 낭비한다. (등급은 이미 `_evaluate` 가
            # 행에 실어 두므로 여기서는 설명만 붙인다.)
            text += self._inherit_note(args.get("shuttle", ""))
            text += self._junction_note(cargo, lk, args.get("shuttle", ""))
            reflect = ("\n\n[반성 요구] 위 실제 점수를 네 **예측과 대조**하라: 예측이 맞았나? "
                       "틀렸다면 어느 축에서 왜 빗나갔나? 그 교훈으로 다음 가설을 세워 설계를 "
                       "이어가거나, 충분히 수렴했으면 finish하라. (점수를 우기지 말고 오라클 결과를 따르라.)")
            return text + reflect, AgentEvent("evaluation",
                                    text="🎨 정밀 설계" + (f" — {detail}" if detail else ""),
                                    data={"rows": rows})
        if name == "optimize_developability":
            return self._optimize_developability(args, cargo)
        return f"unknown tool: {name}", None

    # ── 문헌 근거 등급 상속 — "이 설계체가 어느 셔틀의 후손인가" ──────────────
    # 이 장치가 옛 **순서 특이성(스크램블 델타)** + **접합부 길이 z** 게이트를 승계한다.
    #
    # 왜 교체했나: 두 게이트가 묻던 것은 "이 BBB 점수가 서열 순서에서 왔나 / 조성·길이에서
    # 왔나"였다. 즉 둘 다 **BBB 점수를 순위 근거로 쓴다는 전제** 위에 서 있었다. 그 전제가
    # 두 번 반증됐다 — M6(셔틀 간 순위 무의미) · M8(링커 간도 무의미: 라이브러리 링커가 같은
    # 길이 무작위 링커를 이긴 비율 **0/108**). 전제가 사라지면 질문도 사라진다. 점수를 아예
    # 순위에 쓰지 않으므로 "조성으로 점수를 올리는" 해킹 경로 자체가 닫힌다.
    #
    # 남는 진짜 질문은 하나다: **설계체가 문헌 근거를 물려받을 자격이 있나.**
    # 등급은 우리가 가진 유일한 양(陽)의 신호이고 설계체에는 인용이 없다. 그대로 두면 모든
    # 설계체가 후보 자격을 잃어 **개발성 liability 스크러빙처럼 정당한 기여까지** 막힌다.
    # 그래서 '가장 닮은 라이브러리 셔틀의 등급을 상속'시키되, 유사도 하한을 둔다 — 하한이
    # 없으면 등급이 다시 해킹 가능한 축이 된다.
    # 임계값 0.70 은 측정값이다(M9, docs/experiments.md §18): 정당한 개발성 변이체 55개의
    # 부모 유사도 최소 **0.830** vs 스크램블 0.529 · homopolymer 0.429 · 무작위 0.136 ·
    # 타 셔틀 간 0.088 → 간극 [0.529, 0.830] 안쪽. 최근접이 부모가 아닌 양성 **0/55**.
    def _junction_note(self, cargo: str, linker: str, shuttle: str) -> str:
        """접합부를 걸친 liability 모티프를 **결정론적으로** 알려준다.

        왜 프롬프트 규칙이 아니라 검출인가: 링커가 고칠 수 있는 liability 는 정확히
        **접합부를 걸친 것들**인데(화물·셔틀 내부 모티프는 링커로 못 없앤다), 그게 있는지는
        화물 끝 잔기와 셔틀 첫 잔기에 달려 있어 **화물마다 다르다**. LLM 이 매번 올바르게
        계산하길 기대하는 대신 검출해서 알려준다('검출+복구 한 쌍').

        실측 근거와 `GGS` 가 실패하는 이유는 `developability.junction_liabilities` 참고.
        """
        try:
            from .developability import junction_liabilities
            js = junction_liabilities(cargo, linker, shuttle)
        except Exception:  # noqa: BLE001 — 보조 표시다. 실패가 본 평가를 막아선 안 된다.
            return ""
        if not js:
            return ("\n\n[접합부 검사] ✅ 접합부를 걸친 liability 모티프 **없음** — 이 조합에서는 "
                    "링커를 바꿔 liability 를 줄일 여지가 없다(화물·셔틀 **내부** 모티프는 링커로 "
                    "못 없앤다). 링커로 liability 를 고치려 하지 말 것.")
        items = " · ".join(f"**{j['label']} `{j['motif']}`**@{j['pos']}({j['spans']})" for j in js)
        c = "".join(ch for ch in (cargo or "").upper() if ch.isalpha())
        sh = "".join(ch for ch in (shuttle or "").upper() if ch.isalpha())
        return ("\n\n[접합부 검사] ⚠️ **접합부에서 생긴 liability "
                f"{len(js)}개**: {items}\n"
                f"  이것은 화물·셔틀 내부 문제가 아니라 **붙이는 방식** 때문에 생긴 것이므로 "
                f"**링커를 바꿔 없앨 수 있다** — 지금 이 시스템에서 링커가 고칠 수 있는 유일한 "
                f"liability 다.\n"
                f"  화물 끝 `...{c[-3:] if len(c) >= 3 else c}` ⊕ 셔틀 시작 "
                f"`{sh[:3] if len(sh) >= 3 else sh}...` 이므로, **링커의 첫/끝 잔기**가 이 둘과 "
                f"모티프를 만들지 않게 골라라. ⚠️ 링커를 '끼우는 것'만으로는 안 된다 — 실측으로 "
                f"`GGS`(G 로 시작)는 N 으로 끝나는 화물과 다시 `NG` 를 만들어 **실패**했고, "
                f"`EAAAK`·`AAAAA`(E/A 로 시작)는 제거에 **성공**했다.")

    def _inherit_note(self, shuttle: str) -> str:
        """design_candidate 응답에 붙일 근거 상속 한 단락. 실패해도 예외 없음."""
        try:
            tier, parent, sim = evidence_inherit(shuttle)
        except Exception:  # noqa: BLE001 — 보조 표시다. 실패가 본 평가를 막아선 안 된다.
            return ""
        if sim >= 1.0 and tier:
            return (f"\n\n[문헌 근거] ✅ 라이브러리 셔틀 **{parent}** — 등급 **{tier}**"
                    f"({evidence_group(tier)}).")
        if tier:
            return (f"\n\n[문헌 근거] ✅ **{parent}** 의 후손으로 인정 — 유사도 {sim:.2f} ≥ "
                    f"{INHERIT_SIM_MIN:.2f} → 등급 **{tier}**({evidence_group(tier)}) 상속. "
                    f"부모의 결합 모티프를 보존했다는 뜻이다(기능 보존 제약).")
        near = f"최근접 {parent}·유사도 {sim:.2f}" if parent else "닮은 셔틀 없음"
        return (f"\n\n[문헌 근거] ⚠️ **미등재** — {near} < {INHERIT_SIM_MIN:.2f} 이라 등급을 "
                f"상속받지 못했다. 이 서열은 **최종 후보가 될 수 없다**(심사가 자동 반려). "
                f"문헌 근거는 이 프로젝트에서 '좋다'를 말할 수 있는 **유일한** 근거이고, "
                f"설계체에는 인용이 없으므로 부모 셔틀의 모티프를 보존해야만 물려받는다. "
                f"→ 라이브러리 셔틀에서 **보존적으로** 편집하라(liability 치환 등). "
                f"조성을 쌓거나 길이를 늘리는 방향은 효과가 없다 — BBB 점수는 순위에 쓰지 않는다.")

    def _inherit_gate(self, choice) -> tuple[bool, str]:
        """심사용 근거 상속 게이트 → (미승인시키나, 게이트 문구).

        **하드 게이트다.** 측정이 아니라 자격 판정이므로 '실패 시 통과'가 아니다 —
        등급을 확인할 수 없으면 후보로 올리지 않는다.
        """
        tier = choice.get("tier")
        if tier:
            return False, ""
        parent, sim = choice.get("tier_parent"), choice.get("tier_sim") or 0.0
        near = f"가장 닮은 라이브러리 셔틀은 **{parent}**(유사도 {sim:.2f})" if parent else \
               "닮은 라이브러리 셔틀이 없다"
        return True, (
            f"\n\n[시스템 게이트] 최종 후보의 셔틀이 **문헌 근거 미등재**다 — {near}로 "
            f"상속 하한 {INHERIT_SIM_MIN:.2f} 에 미달 → 자동 미승인(REVISE).\n"
            f"이 프로젝트에서 '이 셔틀이 좋다'를 말할 수 있는 근거는 **문헌 인용뿐**이다"
            f"(우리 BBB 점수는 셔틀 간·링커 간 순위를 매길 수 없음이 M6·M8 로 확인됐다). "
            f"설계체는 부모 셔틀의 결합 모티프를 보존했을 때만 그 근거를 물려받는다. "
            f"하한 {INHERIT_SIM_MIN:.2f} 은 실측 간극에서 왔다(정당한 개발성 변이체 최소 0.830 "
            f"vs 스크램블 0.529·homopolymer 0.429·무작위 0.136).\n"
            f"요구: 라이브러리 셔틀에서 **보존적 치환**(liability 제거 등)으로 다시 설계하거나, "
            f"라이브러리 셔틀을 그대로 쓰라.")

    def _optimize_developability(self, args, cargo=""):
        """셔틀의 개발성(liability·응집·불안정성)을 다목적 Pareto로 최적화 — 기능 보존.
        BBB가 아니라 '제조·안정성' 다듬기. 신뢰 오라클만 사용(규칙·AGGRESCAN·Guruprasad)."""
        from .developability_opt import OFF_BAND, objectives, pareto_scrub
        sh = "".join(ch for ch in (args.get("shuttle", "") or "").upper() if ch.isalpha())
        if len(sh) < 5:
            return ("셔틀 서열이 너무 짧아 개발성 최적화 불가.",
                    AgentEvent("text", text="🧪 개발성 최적화: 셔틀 서열 부족"))
        try:
            front = pareto_scrub(sh, rounds=25, pop=24, rmt_floor=0.85, seed=2022)
        except ValueError as exc:
            # 🔴 `pareto_scrub` 이 **의도적으로** 던지는 설명 있는 실패(§32·§33-4) —
            #    시드 자신이 기능 보존 하한을 위반하는 경우(참조 명부 누락 또는
            #    하한 미달 서열 워밍업). 옛 코드는 `ValueError` 라는 **타입 이름만**
            #    돌려줘서, 조용한 무작동을 고쳐 놓고도 **이유가 다시 가려졌다.**
            #    메시지를 그대로 통과시켜 에이전트가 원인을 읽고 대응하게 한다.
            return (f"개발성 최적화 불가 — {exc}",
                    AgentEvent("text", text=f"🧪 개발성 최적화 불가: 기능 보존 하한 위반 (셔틀 `{sh}`)"))
        except Exception as exc:  # noqa: BLE001
            return (f"개발성 최적화 실패: {type(exc).__name__}: {exc}",
                    AgentEvent("text", text="🧪 개발성 최적화 실패"))
        seed = objectives(sh)
        front, n_dropped, parent = self._filter_front_by_inheritance(sh, front)
        if not front:
            msg = (f"셔틀 {sh}의 개발성 전선이 **전부 근거 상속 하한에 걸려 비었다**"
                   f"(부모 유사도 < {INHERIT_SIM_MIN:.2f}). 변이체가 부모 셔틀의 모티프를 "
                   "너무 많이 잃었다는 뜻이다 — 원본을 그대로 쓰거나 다른 셔틀을 고려하라.")
            return msg, AgentEvent("text", text=f"🧪 개발성 최적화: 근거 상속 미달로 전선 없음 (셔틀 `{sh}`)")

        def _fmt(o):
            return (f"liability {o['liability']}·응집 {o['aggregation']}"
                    f"·불안정 {o['instability']}")
        rows = [f"- `{m['seq']}` — {_fmt(m['objectives'])} · off {m['off_target']}"
                f"(제약 {m['off_seed']}±{OFF_BAND}) · rmt_sim {m['rmt_sim']} · "
                f"**근거 {m['tier']}**(부모 유사도 {m['tier_sim']:.2f}) · 편집 {m['n_edits']}"
                for m in front[:6]]
        body = "\n".join(rows)
        drop = (f" · 근거 상속 미달 **{n_dropped}개 제외**" if n_dropped else "")
        pn = f"`{parent}`" if parent else "(라이브러리 밖 시드)"

        # ── 전선을 **실제로 채점한다** (2026-11-11, `docs/experiments.md` §25) ──────
        # 옛 구현은 전선을 **텍스트로만** 보여주고 끝냈다. 그래서 그 후보들은
        # `seen_rows`·러닝 베스트·원장에 **한 번도 들어가지 않았고**, 래칫 게이트도
        # 그것들을 몰랐다. 실측(M4 콜드 3세션 원장): 설계체 5건이 전부 에이전트의
        # `design_candidate` 산출물(`F2Y·F3R·E18Q`, off 0.20~0.23)이고 **스크럽 전선
        # 모양의 후보는 한 건도 없다** — 도구 출력이 결정에 닿지 않았다는 증거다.
        # 이제 채점해 `evaluation` 이벤트로 내보내면 세 가지가 같이 해결된다:
        #   ⑴ 러닝 베스트가 스크럽 후보를 본다 → 래칫이 **퇴행을 잡을 수 있다**
        #   ⑵ 원장에 남아 **다음 세션이 물려받는다**
        #   ⑶ UI 랭킹·심사가 같은 행 형식으로 함께 본다
        # 링커는 에이전트가 쓰려는 것을 받고(없으면 직접융합) — liability·불안정성은
        # **구성체 전체**, 응집은 **모듈**에서 재므로 링커가 결과를 바꾼다.
        lk = "".join(ch for ch in (args.get("linker", "") or "").upper() if ch.isalpha())
        scored: list[dict] = []
        try:
            _, scored = self._evaluate(cargo, {"constructs": [
                {"label": f"스크럽{i + 1}", "linker": lk, "shuttle": m["seq"]}
                for i, m in enumerate(front[:6])]})
        except Exception:  # noqa: BLE001
            scored = []   # 채점이 실패해도 전선 제시는 유지한다(보고는 끊기지 않게)

        # 순위 축에서 **동일한 점은 하나만** 남긴다 — 옵티마이저는 불안정성을 연속으로
        # 보지만 순위는 임계 40 이진이라, 전선에 "순위 관점에서 같은" 후보가 여럿 생긴다.
        uniq, seen_ax = [], set()
        for r in scored:
            ax = risk_axes(r)
            if ax in seen_ax:
                continue
            seen_ax.add(ax)
            uniq.append(r)
        dedup = (f" · 순위 축 중복 {len(scored) - len(uniq)}개 접음"
                 if len(uniq) < len(scored) else "")
        score_note = ""
        if uniq:
            best_line = min(uniq, key=rank_key)
            score_note = (
                f"\n\n**채점 결과**(순위 축 = 근거 그룹 → liab·불안정성(이진)·off·응집){dedup}: "
                + " / ".join(f"`{r['shuttle']}` {self._summary(r)}" for r in uniq[:4])
                + f"\n→ 이 중 순위 최상은 `{best_line['shuttle']}` 이다. "
                  "**채점됐으므로 러닝 베스트·원장·래칫이 이 후보들을 본다** — "
                  "이보다 나쁜 것을 제출하면 래칫에 걸린다.")

        md = (f"**🧪 개발성 다목적(Pareto) 최적화** — 셔틀 `{sh}` {pn}\n"
              f"시드: liability {seed[0]}·응집 {seed[1]}·불안정 {seed[2]}·**off {seed[3]}**. "
              f"기능 보존(rmt_sim≥0.85) + **근거 상속 보존(부모 유사도 ≥ "
              f"{INHERIT_SIM_MIN:.2f})**하며 **순위 위험 4축**(liability·응집·불안정성·"
              f"off-target) 동시 최소화 → 비지배 전선 {len(front)}개{drop}:\n{body}"
              f"{score_note}")
        llm = (f"셔틀 {sh}의 개발성 Pareto 전선 {len(front)}개{drop}. "
               f"제약: rmt_sim≥0.85 **및 부모({pn}) 유사도 ≥ {INHERIT_SIM_MIN:.2f}"
               f"(= 문헌 근거 등급을 물려받을 수 있는 선)**. 트레이드오프 변이체:\n{body}\n"
               f"{score_note}\n"
               "하나를 골라 design_candidate/finish의 셔틀로 쓰면 개발성(제조·안정성)이 개선되고, "
               "**부모의 근거 등급은 그대로 유지된다**. 이것이 이 시스템에서 설계가 기여하는 "
               "정당한 자리다 — 효능을 주장하는 것이 아니라 **같은 근거를 더 만들 수 있게** 만드는 것. "
               "🔴 **off-target 을 함께 보라** — 양이온 잔기(R·K)를 더하면 off 가 올라가고 "
               "그것은 순위 축이다. 네가 직접 편집할 때 가장 자주 잃는 축이 이것이다. "
               "⚠️ BBB 점수가 변하는 것은 신경 쓰지 마라(순위에 쓰지 않는다).")
        if uniq:
            return llm, AgentEvent("evaluation", text=md, data={"rows": uniq})
        return llm, AgentEvent("text", text=md)

    # ⚰️ DEV_BBB_FLOOR_RATIO(=0.80, 시드 융합 BBB 대비) 는 2026-11-10 **폐지**했다.
    # 그 하한의 목적은 "개발성으로 다듬다가 **효능을 팔아버리는 것**"을 막는 것이었고
    # (실측 §11: T7→HAVYPNH 가 유사도는 지키면서 융합 BBB 를 0.526→0.104 로 떨어뜨렸다),
    # 그 '효능'의 측정자가 BBB 였다. M6·M8 이 그 측정자를 무효화했으므로 하한도 함께 폐지한다.
    # 대체: **근거 등급 상속 하한**(유사도 ≥ INHERIT_SIM_MIN). 지키려던 것이 원래
    # "부모 셔틀의 기능을 잃지 않았나"였으므로, 그걸 BBB 대신 **모티프 보존**으로 직접 잰다.
    # (M9 실측: pareto_scrub 변이체 55개 전부 부모 유사도 ≥0.830 이고 최근접이 부모였다
    #  → 정당한 스크러빙은 이 하한에 걸리지 않는다. 즉 상방을 막지 않는다.)

    def _filter_front_by_inheritance(self, seed_seq, front):
        """개발성 전선에서 **부모 셔틀의 근거를 물려받지 못한 변이체**를 걷어낸다.

        BBB 를 부르지 않으므로 deepB3P subprocess 호출이 **0회**다(옛 구현은 전선 전체를
        배치 채점했다). 반환: (살아남은 전선[등급·유사도 포함], 제외 개수, 부모 셔틀 이름)
        """
        from .evidence import name_of_seq, similarity
        parent = name_of_seq(seed_seq)
        kept = []
        for m in front:
            tier, par, sim = evidence_inherit(m["seq"])
            if tier is None:
                continue
            # 부모가 바뀐 변이체도 버린다 — 다른 셔틀의 근거를 물려받는 건 상속이 아니다.
            if parent is not None and par != parent:
                continue
            kept.append({**m, "tier": tier, "tier_sim": sim})
        return kept, len(front) - len(kept), parent

    @staticmethod
    def _summary(r) -> str:
        """후보 한 줄 요약 — **근거 등급 + 위험 축**. BBB 는 괄호에 참고로만 적는다."""
        t = r.get("tier")
        ev = f"근거 {t}" if t else "근거 미등재"
        if (r.get("tier_sim") or 1.0) < 1.0 and r.get("tier_parent"):
            ev += f"←{r['tier_parent']}·{r['tier_sim']:.2f}"
        return (f"{ev} · liab {r.get('dev_liab', '?')} · 불안정성 {r.get('instability', '?')} · "
                f"off {r.get('sel_off', 0):.2f} · 응집핫스팟 {r.get('dev_agg_hot', '?')} "
                f"(BBB 참고 {r.get('bbb', 0):.3f})")

    def _score_choice(self, cargo, fa):
        """finish가 제출한 최종 후보를 채점해 row(dict)로 반환. 실패 시 None."""
        lk = "".join(ch for ch in (fa.get("chosen_linker", "") or "").upper() if ch.isalpha())
        sh = "".join(ch for ch in (fa.get("chosen_shuttle", "") or "").upper() if ch.isalpha())
        if not sh:
            return None
        _, rows = self._evaluate(cargo, {"constructs": [{
            "label": fa.get("chosen_label", "최종 선택"), "linker": lk, "shuttle": sh}]})
        if not rows:
            return None
        row = rows[0]
        row["agent_pick"] = True
        return row

    # ---- 심사(비평) 에이전트 — 독립 페르소나·독립 컨텍스트 ----------------
    @staticmethod
    def _front_text(front: dict) -> str:
        """도달 가능 집합(Pareto 전선)을 심사자에게 줄 텍스트로 만든다."""
        b = front.get("best") or {}
        bg = front.get("best_group")
        grp = front.get("by_group") or {}
        lines = [
            f"[도달 가능 집합 — 라이브러리 전수 {front.get('n_total', '?')}조합을 **동일 지표**로 채점한 결과]",
            f"  비독성 {front.get('n_nontoxic', '?')}개 → 하드 게이트(**독성뿐**) 통과 "
            f"{front.get('n_gate_all', '?')}개 → 그중 **후보군 {front.get('n_gate', '?')}개** "
            f"(CPP 비교군·근거 등급 D·미등재 제외. 엄밀 RMT {front.get('n_gate_rmt', 0)}개)",
            "  후보군의 **문헌 근거 그룹 분포**: "
            + " · ".join(f"{g} {n}개" for g, n in grp.items() if n) or "—",
            f"  → **달성 가능한 최상위 근거 그룹 = `{bg}`**. 최종 후보는 여기서 나와야 합니다.",
            f"  위험 축별 **달성 한계**: liability 최소 {b.get('liab_min')} · "
            f"불안정성 최소 {b.get('instability_min')} · off_target 최소 {b.get('off_min')} · "
            f"응집 최소 {b.get('agg_min')}",
        ]
        ft = front.get("front_top") or []
        if ft:
            lines.append(f"  **최상위 근거 그룹(`{bg}`)의 Pareto 전선** {len(ft)}개 "
                         "— 최종 후보가 여기 있거나, 이들에게 지배당하지 않아야 합니다:")
            for r in ft[:8]:
                lines.append(
                    f"    · {r['label']}: 근거 **{r['tier']}** · liab {r['dev_liab']} · 불안정성 "
                    f"{r['instability']} · off {r['sel_off']} · 응집 {r['dev_agg_hot']} · {r['mechanism']} "
                    f"(BBB 참고 {r['bbb_report']})")
        lower = [r for r in (front.get("front") or []) if r.get("group") != bg]
        if lower:
            lines.append("  [아래 근거 그룹의 비지배 후보 — 참고. 위 그룹에 후보가 있으면 "
                         "이쪽을 최종안으로 쓸 수 없습니다]")
            for r in lower[:4]:
                lines.append(f"    · {r['label']}: 근거 **{r['tier']}**({r['group']}) · "
                             f"liab {r['dev_liab']} · off {r['sel_off']}")
        unsup = front.get("unsupported") or []
        if unsup:
            lines.append("  [근거 등급 D(= BBB 통과 단독 근거 없음) — **최종 선택 대상 아님**]")
            for r in unsup:
                lines.append(f"    · {r['label']}: 근거 **{r['tier']}** (BBB 참고 "
                             f"{r['bbb_report']}) — 점수가 높아도 근거가 없습니다")
        cpp = front.get("cpp_baseline") or []
        if cpp:
            lines.append("  [비교군(CPP) 참고치 — **최종 선택 대상 아님**. 비특이 흡수 상한]")
            for r in cpp:
                lines.append(f"    · {r['label']}: off {r['sel_off']} · 불안정성 "
                             f"{r['instability']} (BBB 참고 {r['bbb_report']}, 비교군)")
            lines.append("    ※ CPP·등급 D 의 **BBB 점수가 높은 것이 바로 그 점수를 순위에 "
                         "쓸 수 없는 이유**입니다(M1 조성 편향). 그걸 이유로 후보군을 깎지 마십시오.")
        lines.append("  → **이 집합 밖의 값을 요구하지 마십시오.** 여기 적힌 **후보군** 수치가 "
                     "이 문제에서 실제로 달성 가능한 최선입니다.")
        return "\n".join(lines)

    def _critic_prompt(self, front_text: str):
        thr = self.settings.toxicity_threshold
        return (
            "당신은 신약 후보를 **적대적으로 검증하는 안전·엄밀성 심사 에이전트**입니다. 설계 "
            "에이전트가 제안한 최종 BBB 융합체 후보를 의심하고 반박하세요.\n\n"
            "[판정 원칙 — 가장 중요]\n"
            "**절대적 이상과 비교하지 마십시오.** 아래 '도달 가능 집합'이 이 문제에서 실제로 달성 "
            "가능한 최선입니다. 판정 기준은 '이상적인가'가 아니라 **'달성 가능한 최선 대비 합리적인 "
            "선택인가'**입니다.\n\n"
            "[이 시스템이 무엇을 주장하고 무엇을 주장하지 않는가 — 반드시 숙지]\n"
            "  이 시스템은 **BBB 투과도를 예측하지 않습니다.** 문헌 검증 셔틀을 붙이면 뇌 전달이 "
            "올라간다는 것은 **문헌에서 빌려온 전제**이고, 최적화 대상은 그 융합체의 "
            "**개발성·안전성**입니다. 우리 BBB 점수(deepB3P)는 순위를 매기는 데 쓸 수 없습니다:\n"
            "   · M1 — 양성 vs 무작위 AUC 0.807 vs 양성 vs **스크램블** 0.618 → 점수의 상당 부분이 "
            "**아미노산 조성**에서 온다.\n"
            "   · M6 — 융합 수준에서 문헌 검증 셔틀이 자기 스크램블 분포의 39번째 백분위, 유의 0/9, "
            "문헌 음성 대조군 **역전**.\n"
            "   · M8 — 링커 선택도 불가. 라이브러리 링커가 같은 길이 **무작위** 링커를 이긴 비율 "
            "**0/108**(우연이면 5%), 길이↔점수 순위상관 **-0.495**.\n"
            "  👉 따라서 **설계자가 'BBB 점수가 높다/낮다'를 선택 근거로 들면 반박하십시오.** "
            "그리고 **당신도 BBB 점수가 낮다는 이유로 REVISE 하지 마십시오** — 탈락 사유가 아닙니다.\n\n"
            f"{front_text}\n\n"
            "[하드 게이트 — 이 둘만 탈락 사유. 위반하면 예외 없이 REVISE]\n"
            f"  ① **독성**: ToxinPred3 > {thr:.2f} → 탈락.\n"
            "  ② **문헌 근거 등급 보유**: 최종 후보의 셔틀이 라이브러리 셔틀이거나, 설계체라면 "
            f"가장 닮은 라이브러리 셔틀과 유사도 **{INHERIT_SIM_MIN:.2f} 이상**으로 그 등급을 "
            "**상속**받아야 합니다. 미등재(상속 실패)·등급 **D**(BBB 통과 단독 근거 없음)는 탈락.\n"
            "    왜: 문헌 인용은 이 시스템에서 '이 후보가 좋다'를 말할 수 있는 **유일한 근거**입니다. "
            "신뢰하는 나머지 축(독성·응집·liability·불안정성·off-target)은 전부 **'나쁘지 않음'만** "
            "재므로, 근거를 요구하지 않으면 **불활성 쓰레기가 모든 축에서 이깁니다** "
            "(실측: poly-Glu 는 off 0.00·선택성 1.00·응집 0.12·liability 0 으로 Angiopep-2 를 "
            "전 축에서 이깁니다).\n\n"
            "[셔틀 역할 — CPP·등급 D 는 비교군이다]\n"
            "  **RMT**(수용체 매개)가 후보군이고, **CPP**(TAT·SynB1·Penetratin)는 비특이 흡수 "
            "**비교군**입니다. 최종 후보가 CPP라면: 게이트를 통과한 RMT가 **하나라도 있으면 REVISE**"
            "하고 그 대안을 지목하십시오.\n\n"
            f"{evidence_critic_rule()}"
            "[트레이드오프 축 — 탈락 사유가 아니라 '한계로 명시할' 사유]\n"
            "  불안정성 · 개발성 liability/양이온 · off-target · 응집 · 설계의 실질.\n"
            "  이 축들에서는 제출 후보가 전선에서 **지배당하는지만** 보십시오.\n"
            "   · **지배당함** = 근거 그룹이 같거나 더 좋고 위험 축 전부가 같거나 더 나은 "
            "**달성 가능** 후보가 실재 → REVISE, 그 후보를 지목.\n"
            "   · **비지배** = 트레이드오프 상의 합리적 지점 → 한계를 명시하고 조건부 승인.\n"
            "   · **근거 하향** = 비지배여도 **달성 가능한 최상위 근거 그룹보다 낮은 그룹**을 "
            "골랐으면 REVISE. 위험 축이 조금 좋다는 이유로 근거가 약한 셔틀을 고를 수는 없습니다 "
            "— 위험 축은 '좋음'을 재지 않기 때문입니다.\n"
            "   · 등급 **C(주의)** 후보를 승인할 때는 **왜 주의인지**(문헌의 고효율 형태가 이 서열과 "
            "다름·정량값 미확보 등)를 한계로 **반드시 명시**하게 하십시오.\n\n"
            "[판정 — 마지막 줄에 정확히 하나]\n"
            "  `VERDICT: APPROVE`              하드 게이트 통과 + 비지배 + 최상위 근거 그룹 + 명시할 한계 없음\n"
            "  `VERDICT: APPROVE_WITH_LIMITS`  하드 게이트 통과 + 비지배이나 명시할 한계가 있음 "
            "(정상적인 승인 형태)\n"
            "  `VERDICT: REVISE`               하드 게이트 위반, 전선에서 지배당함, 근거 하향, "
            "또는 RMT 대안이 있는데 CPP를 고름\n"
            "**애매하면 REVISE가 아니라 `APPROVE_WITH_LIMITS`로 하고 한계를 적으십시오.** 단 하드 게이트 "
            "위반·근거 하향·'RMT 대안이 있는데 CPP 선택'은 예외 없이 REVISE입니다.\n"
            "3~5문장으로 핵심 평가(무엇이 한계이고 그래도 왜 합리적인지, 또는 어느 후보에 지배당하는지)를 "
            "쓰고 **마지막 줄에 VERDICT**를 쓰십시오."
        )

    # ⚰️ RATCHET_MARGIN(=0.02) 은 2026-11-10 **폐지**했다. 그것은 유효점수(BBB 기반)
    # **차이**를 재던 여유폭이고, 잡음 있는 연속 점수를 비교하니 필요했다. 래칫을 **지배 관계**
    # (근거 그룹 + 위험 축 Pareto)로 바꾸면서 사라졌다 — 축들이 전부 결정론적 측정이라
    # 잡음 여유폭이 필요 없고, 축 사이의 임의 가중도 함께 없어진다.

    def _critic_review(self, client, model, types, cargo, choice, best_seen=None):
        """최종 후보를 별도 컨텍스트로 적대적 검증.

        반환 `(approve: bool, critique_text: str, tier: str)`.
        tier ∈ {APPROVE, APPROVE_WITH_LIMITS, REVISE} — **계층 판정**이다.
        하드 게이트(**독성 · 문헌 근거 등급 보유**)만 탈락 사유이고, 나머지 축은
        '도달 가능 집합(Pareto 전선) 대비 지배당하는가'로만 본다. 절대적 이상과 비교하면 승인 가능 영역이 비기 때문이다
        (실측: 심사 기준 AND → 130조합 중 0개. `reachable_front` 주석 참고).
        """
        front = self.reachable_front(cargo)
        dev_list = choice.get("dev_liabilities") or []
        m = (f"화물={cargo}, 링커={choice.get('linker') or '—'}, 셔틀={choice.get('shuttle')}, "
             f"전체서열={choice.get('sequence')}\n"
             f"지표 — BBB점수={choice['bbb']:.3f}(0~1), 독성={choice['tox']:.3f}, "
             f"불안정성지수={choice['instability']}, "
             f"수용체유사도({choice.get('bind_ref')})={choice.get('bind_score', 0):.2f}, "
             f"독성판정={'탈락' if choice['toxic'] else '통과'}\n"
             f"전달 분해 — 메커니즘={choice.get('mechanism', '?')}, 셔틀단독BBB={choice.get('shuttle_bbb', '?')}, "
             f"융합보존={choice.get('preservation', 'n/a')}, deepB3P 타당도=\"{choice.get('deepb3p_valid', '?')}\"\n"
             f"개발성 — 위험={choice.get('dev_risk', '?')}, 순전하={choice.get('dev_charge', '?')}, "
             f"응집 핫스팟(설계모듈)={choice.get('dev_agg_hot', '?')}"
             f"/전체구성체={choice.get('dev_agg_hot_full', '?')}, "
             f"liability({len(dev_list)}): "
             + ("; ".join(dev_list) if dev_list else "없음") + "\n"
             f"선택성/off-target — 위험={choice.get('sel_level', '?')}"
             f"(선택성 {choice.get('selectivity', '?')}), {choice.get('sel_mech', '?')}\n"
             f"용해도 — {choice.get('sol_level', '?')}(점수 {choice.get('sol_score', '?')})")
        cfg = types.GenerateContentConfig(
            system_instruction=self._critic_prompt(self._front_text(front)), temperature=1.0,
            thinking_config=types.ThinkingConfig(include_thoughts=False))
        contents = [types.Content(role="user", parts=[types.Part.from_text(
            text="설계 에이전트가 제출한 최종 후보다. 도달 가능 집합 대비로 검증하고 VERDICT를 "
                 "내려라.\n" + m)])]
        resp = self._gen_retry(client, model, contents, cfg)
        cand = (resp.candidates or [None])[0]
        parts = (cand.content.parts if cand and cand.content else None) or []
        text = "\n\n".join(p.text for p in parts
                           if getattr(p, "text", None) and not getattr(p, "thought", False))
        tail = (text or "").upper().rsplit("VERDICT:", 1)[-1]
        if "APPROVE_WITH_LIMITS" in tail or "APPROVE WITH LIMITS" in tail:
            tier = "APPROVE_WITH_LIMITS"
        elif "APPROVE" in tail:
            tier = "APPROVE"
        else:
            tier = "REVISE"

        # ── 결정론적 게이트 (LLM 판정을 덮어쓴다) ─────────────────────────
        # 하드 게이트는 둘뿐: ① 독성 ② **문헌 근거 등급 보유**(상속 포함).
        # 옛 ②(순서 특이성 측정 가능성)·전선 이탈(유효점수 하한)은 2026-11-10 폐지했다 —
        # 둘 다 "BBB 점수로 순위를 매길 수 있다"는 전제에 서 있었고 M6·M8 이 반증했다.
        # 그 역할은 ②가 승계한다: 미등재(상속 실패) 설계체 = homopolymer·조성 쌓기·무작위가
        # 전부 여기서 걸린다(M9 실측 간극 [0.529, 0.830]).
        from .evidence import group_rank as _grank
        if choice.get("toxic"):
            text = (text or "") + (
                f"\n\n[시스템 게이트] 독성 {choice.get('tox', 0.0):.3f} > 임계값 "
                f"{self.settings.toxicity_threshold:.2f} → 안전성 하드 게이트 위반으로 자동 미승인(REVISE).")
            tier = "REVISE"
        elif (ig := self._inherit_gate(choice))[0]:
            # ── 근거 등급 보유 게이트 — '인용 없는 서열' 금지 ──
            # 이 게이트가 옛 BBB 하한·순서 특이성·접합부 길이 통제의 역할을 **한꺼번에** 맡는다.
            # 근거: 신뢰하는 나머지 축은 전부 '나쁘지 않음'만 재므로(위험 스크린), 근거를
            # 요구하지 않으면 불활성 쓰레기가 전 축에서 이긴다 — 실측 poly-Glu(off 0.00·
            # 선택성 1.00·응집 0.12·liability 0)가 Angiopep-2(0.09·0.91·0.44·2)를 전부 이긴다.
            text = (text or "") + ig[1]
            tier = "REVISE"
        elif is_unsupported(choice):
            # 등급 D = **BBB 통과 단독 근거 없음**(Tet1). 상속은 됐지만 물려받은 근거가 없다.
            # 실측 근거: M2·M4 에서 우리 효능 축의 후보군 1위가 `직접융합+Tet1`(BBB 0.568)이었다
            # — 점수 1위가 라이브러리에서 BBB 근거가 가장 약한 셔틀이었다.
            alts = ", ".join(f"{r['label']}(근거 {r['tier']})"
                             for r in (front.get("front_top") or [])[:3])
            text = (text or "") + (
                f"\n\n[시스템 게이트] 최종 후보의 셔틀이 **근거 등급 "
                f"{choice.get('tier')}(= BBB 통과 단독 근거 없음)**이다 → 자동 미승인(REVISE). "
                f"BBB 점수가 높게 나오는 것은 선택 근거가 아니다(M1 조성 편향·M6·M8). "
                f"달성 가능한 대안: {alts}.")
            tier = "REVISE"
        elif (str(choice.get("mechanism") or "").upper().startswith("CPP")
              and (front.get("n_gate_rmt") or 0) > 0):
            # CPP는 비특이 흡수 **비교군**이다. 게이트를 통과한 RMT 후보가 실제로 존재하면
            # 비교군을 최종안으로 확정할 이유가 없다 → 대안을 지목해 REVISE.
            alts = ", ".join(f"{r['label']}(근거 {r['tier']}·liab {r['dev_liab']})"
                             for r in (front.get("rmt_options") or [])[:3])
            text = (text or "") + (
                f"\n\n[시스템 게이트] 최종 후보가 **CPP(비특이 흡수 비교군)**인데, 하드 게이트를 통과한 "
                f"**RMT 후보가 {front.get('n_gate_rmt')}개 실재**한다 → 비교군을 최종안으로 확정할 수 "
                f"없으므로 자동 미승인(REVISE). 달성 가능한 RMT 대안: {alts}.")
            tier = "REVISE"
        elif (front.get("best_group") is not None
              and _grank(choice.get("tier")) > _grank((front.get("front_top") or [{}])[0].get("tier"))):
            # ── 근거 하향 게이트 — 옛 '전선 이탈'의 승계자 ──
            # 비지배여도 **달성 가능한 최상위 근거 그룹보다 낮은 그룹**을 골랐으면 승인하지 않는다.
            # 위험 축(liability·off·불안정성·응집)이 조금 좋다는 이유로 근거가 약한 셔틀을 고를 수
            # 없다 — 그 축들은 '좋음'을 재지 않기 때문이다. 하한은 임의 상수가 아니라
            # **전선 자체에서 도출**된 값이다(달성 가능한 최상위 그룹).
            alts = ", ".join(f"{r['label']}(근거 {r['tier']}·liab {r['dev_liab']}·off {r['sel_off']})"
                             for r in (front.get("front_top") or [])[:3])
            text = (text or "") + (
                f"\n\n[시스템 게이트] 제출 후보의 문헌 근거 그룹이 "
                f"**`{evidence_group(choice.get('tier'))}`**(등급 {choice.get('tier')})인데, "
                f"이 문제에서 달성 가능한 최상위 그룹은 **`{front.get('best_group')}`**이다 → "
                f"**근거 하향**이므로 자동 미승인(REVISE). 달성 가능한 대안: {alts}.")
            tier = "REVISE"
        elif (dom := next((r for r in (front.get("front_top") or [])
                           if dominates(r, choice)), None)) is not None:
            # ── 지배 게이트 — 모든 축에서 더 나은 '달성 가능' 후보가 실재 ──
            # 근거 그룹 + 위험 축 전부에서 같거나 낫고 하나는 엄격히 나은 후보가 있으면
            # 트레이드오프가 아니라 **그냥 못한 선택**이다. BBB 는 축에 들어가지 않는다.
            text = (text or "") + (
                f"\n\n[시스템 게이트] 제출 후보가 전선의 **`{dom['label']}`**"
                f"(근거 {dom['tier']}·liab {dom['dev_liab']}·불안정성 {dom['instability']}·"
                f"off {dom['sel_off']}·응집 {dom['dev_agg_hot']})에게 **모든 축에서 지배당한다** → "
                f"트레이드오프가 아니라 열등한 선택이므로 자동 미승인(REVISE). 그 후보를 쓰거나, "
                f"그것이 지배하지 못하는 설계를 제출하라.")
            tier = "REVISE"
        elif best_seen is not None and dominates(best_seen, choice):
            # ── 래칫(ratchet) 게이트 — '이미 본 더 좋은 후보를 버리고 제출' 금지 ──
            # 애블레이션 실측(docs/experiments.md §11): 16세션 중 **7세션이 퇴행**했고
            # 편집 채널이 있는 조건에서만 일어났다(라이브러리 선택만 하는 조건은 퇴행 0건).
            # 원인은 **되돌아갈 장치가 없는 것**이다. 탐색은 자유롭게 하되 **제출은 이미 본
            # 최고치에 지배당하지 못하게** 막는다 — 설계의 상방은 살리고 하방만 자른다.
            # 2026-11-10: 비교를 유효점수 차이(여유폭 0.02)에서 **지배 관계**로 바꿨다.
            # 축들이 전부 결정론적 측정이라 잡음 여유폭이 필요 없고, 축 간 임의 가중도 사라진다.
            text = (text or "") + (
                f"\n\n[시스템 게이트] 제출 후보가 **이 세션이 이미 채점해 본** 후보 "
                f"`{best_seen.get('linker') or '직접융합'}+{best_seen.get('shuttle', '?')}`"
                f"({self._summary(best_seen)})에게 **모든 축에서 지배당한다** → 더 좋은 걸 "
                f"**이미 손에 쥐고도** 더 나쁜 것을 확정하려 했으므로 자동 미승인(REVISE). "
                f"그 후보로 되돌아가거나, 그보다 나은 설계를 찾아 제출하라.")
            tier = "REVISE"

        return tier != "REVISE", text, tier

    def _recommend_on_failure(self, client, model, types, cargo, choice, critique_text):
        """최종 미승인(전 후보 반려)일 때, 실패 진단 + 시스템 내 구체적 대안을 권고한다.

        화물(치료 목표)은 고정이므로 화물 교체가 아니라 셔틀·링커·잔기 편집·개발성 최적화 관점의
        '다음 수'를 제시한다. 실패 경로에서만 1회 호출되어 비용 부담이 작다. 반환: 권고안 텍스트.
        """
        dev_list = choice.get("dev_liabilities") or []
        m = (f"화물={cargo}, 링커={choice.get('linker') or '—'}, 셔틀={choice.get('shuttle')}\n"
             f"BBB={choice['bbb']:.3f}, 독성={choice['tox']:.3f}, 불안정성={choice['instability']}, "
             f"선택성={choice.get('selectivity', '?')}({choice.get('sel_level', '?')}), "
             f"개발성위험={choice.get('dev_risk', '?')}, 용해도={choice.get('sol_level', '?')}, "
             f"liability={len(dev_list)}개")
        sys = (
            "너는 BBB 융합체 설계의 자문(advisory) 에이전트다. 최종 후보가 심사에서 **최종 미승인**됐다. "
            "화물(치료 목표)은 고정이므로 **화물 교체는 권하지 마라.** 실패 원인을 데이터로 진단하고, "
            "**우리 시스템 안에서 실행 가능한 구체적 다음 수**를 우선순위로 제시하라.\n"
            "선택 가능한 대안(근거를 들어 골라라):\n"
            "  1) design_candidate로 **라이브러리 셔틀의 잔기 편집**(검증 리간드에서 출발한 보존적 변이·"
            "트리밍) — 셔틀은 de-novo 생성하지 말고 검증된 셔틀을 개량한다.\n"
            "  2) 링커 교체 — 셔틀이 구조상 파묻힐(occluded) 땐 (GGGGS)3 등 긴 유연 링커로 노출 개선, "
            "반대로 너무 길어 불안정하면 직접융합/짧은 링커로.\n"
            "  3) 라이브러리 내 다른 셔틀 — 선택성(RMT형) 높은 쪽으로.\n"
            "  4) 실패가 화물 자체 특성(독성·불안정·과소수성) 때문이면 그 사실을 명시하고, 셔틀로는 "
            "해결 불가임을 정직히 밝힌 뒤 화물 측 개량 또는 다른 전달 모달리티(향후)를 권고.\n"
            "출력: 한국어로 '### 🔧 권고안(대안)' 제목 아래 ①실패 진단 1~2문장 + ②우선순위 대안 2~3개"
            "(각 한 줄, 근거 포함). **새 서열을 지어내지 마라**(생성은 도구의 몫).")
        cfg = types.GenerateContentConfig(
            system_instruction=sys, temperature=1.0,
            thinking_config=types.ThinkingConfig(include_thoughts=False))
        contents = [types.Content(role="user", parts=[types.Part.from_text(
            text=f"최종 미승인된 후보:\n{m}\n\n심사 반박:\n{critique_text}\n\n권고안을 작성하라.")])]
        resp = self._gen_retry(client, model, contents, cfg)
        cand = (resp.candidates or [None])[0]
        parts = (cand.content.parts if cand and cand.content else None) or []
        return "\n\n".join(p.text for p in parts
                           if getattr(p, "text", None) and not getattr(p, "thought", False))

    def _gen_retry(self, client, model, contents, config, tries=6):
        """일시적 오류(503/429)는 백오프 재시도. 그 외는 즉시 raise."""
        last = None
        for attempt in range(tries):
            try:
                return client.models.generate_content(model=model, contents=contents, config=config)
            except Exception as exc:  # noqa: BLE001
                last = exc
                if attempt < tries - 1 and any(t in str(exc) for t in _TRANSIENT):
                    time.sleep(2 * (attempt + 1))  # 2·4·6초
                    continue
                raise
        raise last  # pragma: no cover

    def run(self, cargo: str) -> Generator[AgentEvent, None, None]:
        cargo = cargo.strip().upper()
        key = self.settings.gemini_api_key
        if not key:
            yield AgentEvent("error", "GEMINI_API_KEY(또는 GOOGLE_API_KEY) 미설정.")
            return
        try:
            from google import genai
            from google.genai import types
        except ImportError:
            yield AgentEvent("error", "google-genai 미설치: pip install google-genai")
            return

        client = genai.Client(
            api_key=key,
            http_options=types.HttpOptions(timeout=_REQUEST_TIMEOUT_MS),
        )
        model = self.settings.gemini_model
        tools = types.Tool(function_declarations=self._function_decls(types))
        system = _gemini_system_prompt(cargo, self.settings.toxicity_threshold,
                                       self.max_rounds)

        def _cfg(with_tools=True):
            kw = dict(system_instruction=system, temperature=1.0,
                      thinking_config=types.ThinkingConfig(include_thoughts=True))
            if with_tools:
                kw["tools"] = [tools]
            return types.GenerateContentConfig(**kw)

        task = (f"화물 `{cargo}`에 최적인 링커·셔틀을 설계하라. 도구를 자율적으로 사용해 "
                "BBB·독성·구조·안정성을 모두 만족하는 최종 융합체를 찾아라.")
        # 설계 이력 요약 — 지난 세션들이 실측한 상위 후보와 실패를 **바닥**으로 깔아준다.
        # 세션 내 래칫이 '한 세션 안에서 퇴행 금지'라면, 이건 그걸 세션 경계 밖으로 확장한 것이다.
        ledger_brief = ""
        if self.use_ledger:
            try:
                from .ledger import brief as _ledger_brief
                ledger_brief = _ledger_brief(cargo)
            except Exception:  # noqa: BLE001
                ledger_brief = ""   # 보조 장치다. 실패해도 세션은 그대로 돌아간다
        if ledger_brief:
            yield AgentEvent("text", text="📒 설계 이력을 불러왔습니다 — 지난 세션들의 실측 "
                                          "상위 후보·실패 기록을 출발점으로 씁니다.")
        contents = [types.Content(role="user", parts=[types.Part.from_text(
            text=task + ledger_brief)])]
        seen_rows = []     # 이번 세션이 채점한 모든 행 — 원장 기록용
        best = None        # 후보군(CPP 제외) 러닝 베스트
        best_cpp = None    # 비교군 최고치 — 보고용이며 최종 선택이 되지 않는다

        def _blocked(cand, resp):
            fr = str(getattr(cand, "finish_reason", "") or "")
            if any(x in fr.upper() for x in ("SAFETY", "BLOCK", "PROHIBITED", "RECITATION")):
                return fr
            fb = getattr(resp, "prompt_feedback", None)
            if fb and getattr(fb, "block_reason", None):
                return f"prompt blocked: {fb.block_reason}"
            return None

        def _visible(parts):
            return "\n\n".join(p.text for p in parts
                               if getattr(p, "text", None) and not getattr(p, "thought", False))

        def _user(text):
            return types.Content(role="user", parts=[types.Part.from_text(text=text)])

        try:
            # ── 0. 계획(PLAN): 도구를 쓰기 전에 전략 선언 ──
            contents.append(_user("도구를 호출하기 전에, 최적화 전략을 2~3문장으로 밝혀라 — 어떤 "
                                  "링커·셔틀 방향을 우선 탐색하고 왜인지. 아직 도구는 호출하지 마라."))
            presp = self._gen_retry(client, model, contents, _cfg(False))
            pcand = (presp.candidates or [None])[0]
            blk = _blocked(pcand, presp) if pcand is not None else "빈 응답"
            if blk:
                yield AgentEvent("error", f"Gemini 차단/거부: {blk}")
                return
            pparts = (pcand.content.parts if pcand and pcand.content else None) or []
            for p in pparts:
                if getattr(p, "text", None) and getattr(p, "thought", False):
                    yield AgentEvent("reasoning", p.text)
            contents.append(pcand.content)
            yield AgentEvent("plan", _visible(pparts))
            contents.append(_user("이제 전략에 따라 도구를 자율적으로 사용해 실행하라."))

            # ── 1. 설계자 실행 + 자기종료(finish) → 심사(Critic) → REVISE 시 재설계 ──
            choice = final_report = critique_text = None
            approve, tier = True, "APPROVE"
            critiques, MAX_CRITIQUES, turn, grace = 0, 1, 0, 0
            while turn < self.max_rounds + grace:
                turn += 1
                resp = self._gen_retry(client, model, contents, _cfg(True))
                cand = (resp.candidates or [None])[0]
                blk = _blocked(cand, resp) if cand is not None else "빈 응답"
                if blk:
                    yield AgentEvent("error", f"Gemini 차단/거부: {blk}")
                    return
                parts = (cand.content.parts if cand and cand.content else None) or []
                contents.append(cand.content)

                fcalls = []
                for p in parts:
                    if getattr(p, "function_call", None):
                        fcalls.append(p.function_call)
                    elif getattr(p, "text", None):
                        yield AgentEvent("reasoning" if getattr(p, "thought", False) else "text", p.text)

                if not fcalls:  # 도구 없이 결론 텍스트 → 종료로 간주
                    final_report = _visible(parts)
                    break

                fr_parts, finish_fc = [], None
                for fc in fcalls:
                    if fc.name == "finish":
                        finish_fc = fc
                        fr_parts.append(types.Part.from_function_response(
                            name="finish", response={"result": "acknowledged; self-evaluating"}))
                        continue
                    text, ev = self._dispatch(cargo, fc.name, dict(fc.args or {}))
                    if ev is not None:
                        yield ev
                        if ev.kind == "evaluation":
                            seen_rows.extend(ev.data["rows"])
                            for r in ev.data["rows"]:
                                # 러닝 베스트는 **후보군(CPP 비교군 제외)**에서만 고른다.
                                # CPP는 비교 근거로만 쓰고 best_cpp에 따로 보관한다.
                                # 순서 특이성 탈락 후보는 러닝 베스트가 될 수 없다 —
                                # 그래야 최종 교체(복구)가 **통과하는** 후보로 되돌아간다.
                                # `is_candidate` 가 독성·CPP·**근거 미등재/등급 D** 를 모두
                                # 걸러낸다. 순위는 rank_key(근거 그룹 → 위험 축).
                                if is_candidate(r) and (best is None
                                                        or rank_key(r) < rank_key(best)):
                                    best = r
                                elif is_comparison_only(r) and not r.get("toxic") and (
                                        best_cpp is None or rank_key(r) < rank_key(best_cpp)):
                                    best_cpp = r
                    fr_parts.append(types.Part.from_function_response(
                        name=fc.name, response={"result": text}))

                if finish_fc is None:
                    # 예산이 거의 소진됐는데 아직 finish를 안 했으면 → 다음 턴에 finish 유도
                    # (제출해야 독립 심사를 받으므로, 마지막에라도 최적 후보를 내게 한다).
                    if (self.max_rounds + grace) - turn <= 1:
                        fr_parts.append(types.Part.from_text(
                            text="⚠️ 행동 예산이 거의 소진됐다. 다음 호출에서 반드시 `finish`를 불러 "
                                 "지금까지의 최적 후보를 제출하라 — 제출해야 독립 심사(Critic)를 받고 확정된다."))
                    contents.append(types.Content(role="user", parts=fr_parts))
                    if best:
                        yield AgentEvent("progress", data={"round": turn, "best_bbb": best["bbb"]})
                    continue

                # finish → 최종 후보 채점 → 독립 심사(Critic) → APPROVE/REVISE
                fa = dict(finish_fc.args or {})
                final_report = fa.get("final_report", "") or final_report
                ch = self._score_choice(cargo, fa)
                if ch is not None:
                    choice = ch
                if choice is not None:
                    # best(후보군 러닝 베스트)를 함께 넘겨 **퇴행 제출**을 막는다(래칫 게이트).
                    approve, critique_text, tier = self._critic_review(
                        client, model, types, cargo, choice, best_seen=best)
                else:
                    approve, critique_text, tier = True, "후보 채점 실패로 심사를 생략합니다.", "APPROVE"
                yield AgentEvent("progress", data={"round": turn,
                                                   "best_bbb": (choice or best or {}).get("bbb", 0)})
                if (not approve) and critiques < MAX_CRITIQUES:
                    critiques += 1
                    grace += 1  # 재설계용 여분 턴 부여 → 마지막 스텝에서 REVISE가 나와도 개선 1회 보장
                    contents.append(types.Content(role="user", parts=fr_parts + [types.Part.from_text(
                        text="독립 **심사 에이전트가 REVISE(개선 요구)**를 냈다. 심사 반박:\n"
                             f"{critique_text}\n이 지적을 반영해 design_candidate 등으로 후보를 개선하고, "
                             "끝나면 다시 finish를 호출하라.")]))
                    continue
                contents.append(types.Content(role="user", parts=fr_parts))  # 함수응답 기록 후 종료
                break

            # 강제 종료 대비: 최종 보고서·심사 확보
            if final_report is None:
                contents.append(_user("행동 예산을 모두 사용했습니다. 지금까지의 최적 융합체(전체 서열, "
                                      "네 지표)와 근거로 최종 보고서를 쓰세요. 도구 호출 금지."))
                resp = self._gen_retry(client, model, contents, _cfg(False))
                cand = (resp.candidates or [None])[0]
                parts = (cand.content.parts if cand and cand.content else None) or []
                contents.append(cand.content)
                final_report = _visible(parts)
            # ── 비교군(CPP) 제출은 최종 결과로 올리지 않는다 ──────────────────
            # 심사는 이미 자동 미승인(REVISE)을 냈고 재설계 기회도 소진된 상태다. 이때 CPP를
            # 포디움 1위로 올리면 "비교군이 결과로 나온다" — 그래서 후보군 러닝 베스트로 교체한다.
            # (후보군을 하나도 못 찾은 경우에만 결과 없음으로 둔다.)
            swapped_out = False          # 제출 후보를 러닝 베스트로 바꿨는가(재심사 필요)
            if choice is not None and is_comparison_only(choice):
                swapped = best is not None
                swapped_out = swapped
                final_report = (final_report or "") + (
                    "\n\n> **[시스템] 비교군 대체** — 제출된 최종 후보 "
                    f"`{choice.get('shuttle', '?')}`는 비특이 흡수 **비교군(CPP)**이므로 최종 결과로 "
                    "확정하지 않습니다." + (
                        f" 후보군(RMT) 중 최고 후보 `{best.get('shuttle', '?')}`"
                        f"({self._summary(best)})로 대체해 제시합니다."
                        if swapped else
                        " 후보군(RMT)에서 유효한 후보를 찾지 못해 최종 후보 없이 종료합니다."))
                choice = None
            # ── 근거 없는 제출은 러닝 베스트로 되돌린다 (게이트의 '복구' 절반) ──
            # 게이트는 위반을 **검출**해 REVISE를 내지만, 재설계 기회(MAX_CRITIQUES)를 소진하면
            # 그 후보가 그대로 최종이 된다 → 승인 0건으로 끝난다(실측: 순서 게이트 때 2세션 모두,
            # 래칫 때 4세션 중 3세션). 검출만으로는 못 고친다. 그래서 CPP 대체와 **같은 방식으로**
            # 이미 채점해 본 최고 후보로 바꿔 제시한다. 탐색 기록(best)은 오라클 실측이다.
            elif choice is not None and best is not None and (
                    choice.get("tier") is None or is_unsupported(choice)):
                why = ("문헌 근거 **미등재**(가장 닮은 라이브러리 셔틀과 유사도 "
                       f"{(choice.get('tier_sim') or 0.0):.2f} < {INHERIT_SIM_MIN:.2f})"
                       if choice.get("tier") is None else
                       f"근거 등급 **{choice.get('tier')}**(= BBB 통과 단독 근거 없음)")
                final_report = (final_report or "") + (
                    f"\n\n> **[시스템] 근거 미달 대체** — 제출 후보 "
                    f"`{(choice.get('linker') or '직접융합')}+{choice.get('shuttle', '?')}`는 "
                    f"{why}이므로 최종 결과로 확정하지 않습니다. 이 시스템에서 '좋다'를 말할 수 "
                    f"있는 근거는 문헌 인용뿐이고(우리 BBB 점수는 M6·M8 로 순위 근거가 될 수 "
                    f"없음이 확인됐습니다), 근거를 요구하지 않으면 불활성 서열이 모든 위험 축에서 "
                    f"이깁니다. 근거를 보유한 "
                    f"`{(best.get('linker') or '직접융합')}+{best.get('shuttle', '?')}`"
                    f"({self._summary(best)})로 대체해 제시합니다.")
                choice = None
                swapped_out = True
            elif choice is not None and best is not None and dominates(best, choice):
                final_report = (final_report or "") + (
                    "\n\n> **[시스템] 퇴행 대체** — 제출된 최종 후보 "
                    f"`{choice.get('linker') or '직접융합'}+{choice.get('shuttle', '?')}`"
                    f"({self._summary(choice)})가 **이 세션이 이미 채점해 본** "
                    f"`{best.get('linker') or '직접융합'}+{best.get('shuttle', '?')}`"
                    f"({self._summary(best)})에게 **모든 축에서 지배당하므로** 후자로 대체해 "
                    "제시합니다. 설계 탐색 중 더 나쁜 후보로 확정되는 것을 막는 장치입니다.")
                choice = None          # → final_pick = best
                swapped_out = True
            final_pick = choice or best
            # 교체가 일어났으면 **교체된 후보로 다시 심사**한다 — 이전 판정은 버린 후보의 것이다.
            # (교체가 없고 이미 판정이 있으면 재호출하지 않는다 — 불필요한 LLM 비용.)
            if final_pick is not None and (critique_text is None or swapped_out):
                approve, critique_text, tier = self._critic_review(
                    client, model, types, cargo, final_pick, best_seen=best)

            # 최종 미승인(막다른 경우) → 실패 진단 + 대안 권고를 최종 보고서에 덧붙인다(실패 때만 1회).
            if final_pick is not None and not approve:
                try:
                    rec = self._recommend_on_failure(client, model, types, cargo,
                                                     final_pick, critique_text or "")
                    if rec:
                        final_report = (final_report or "") + "\n\n" + rec
                except Exception:  # noqa: BLE001 - 권고는 부가기능, 실패해도 결과는 낸다
                    pass

            # ── 설계 이력 원장에 이번 세션을 기록 ──
            # 다음 세션이 이 결과를 **바닥**으로 깔고 시작한다(세션 간 래칫). 실패 사유도
            # 같이 남겨 같은 벽에 다시 부딪히는 라운드를 아낀다.
            if self.use_ledger and self.ledger_write:
                try:
                    from .ledger import failure_notes as _ledger_notes
                    from .ledger import record as _ledger_record
                    _ledger_record(cargo, seen_rows, final_pick, _ledger_notes(seen_rows))
                except Exception:  # noqa: BLE001
                    pass    # 기록 실패가 결과를 무효화하진 않는다

            # 비교군 최고치는 '참고선'으로만 보고서에 남긴다 — 결과 카드에는 올리지 않는다.
            if best_cpp is not None:
                final_report = (final_report or "") + (
                    f"\n\n> **비교군 참고선(CPP, 최종 선택 대상 아님)** — "
                    f"`{best_cpp.get('shuttle', '?')}` — {self._summary(best_cpp)}. "
                    "비특이 흡수 상한을 보여주는 대조군이며, 기전(비수용체)과 선택성 때문에 "
                    "후보군에서 제외됩니다. **BBB 점수가 높게 나오는 것이 바로 그 점수를 "
                    "순위에 쓸 수 없는 이유입니다**(M1 조성 편향).")

            # ── 결과 방출 (한 번씩): 설계자 결론 + 심사 판정 ──
            # 심사 최종 판정을 최종 후보에 각인 → UI가 "승인/미승인"을 정직하게 표시한다.
            # finish로 고른 choice가 있으면 그것을, 없으면(예산 소진) best를 최종 선택으로
            # 승격한다 → 심사 판정이 항상 포디움 1위 + 하단 심사 카드에 노출된다.
            pick = choice if choice is not None else final_pick
            if pick is not None:
                pick["critic_approved"] = approve
                pick["critic_tier"] = tier
                pick["agent_pick"] = True
                yield AgentEvent("choice", data=pick)
            yield AgentEvent("final", final_report or "")
            if final_pick is not None:
                # 심사(적대적 검증) 카드는 최종 후보가 있으면 항상 노출한다 —
                # 모델이 빈 텍스트를 돌려줘도 판정 문구를 채워 카드가 사라지지 않게.
                yield AgentEvent(
                    "critique",
                    critique_text or ("하드 게이트(독성·BBB 하한)를 통과하고 도달 가능 전선에서 "
                                      "지배당하지 않습니다. 최종 선택을 승인합니다."
                                      if approve else "근거가 불충분합니다. 개선이 필요합니다."),
                    data={"approve": approve, "tier": tier})
                final_pick["critic_approved"] = approve
                final_pick["critic_tier"] = tier
                # ── 구조 증거(ESMFold) — **에이전트 판단 밖**, 최종 후보 1건 보고용 ──
                # analyze_structure 를 도구에서 뺀 뒤에도 구조 근거 자체는 버리지 않는다.
                # 다만 (a) 의사결정에 넣지 않고 (b) 1회만 접고 (c) 실패하면 조용히 건너뛴다.
                # 공개 API 장애가 결과를 막지 못하게 하려는 것이다.
                try:
                    stext, sdata = self._structure(
                        cargo, {"linker": final_pick.get("linker", ""),
                                "shuttle": final_pick.get("shuttle", "")})
                    if not (sdata or {}).get("error"):
                        final_pick["structure"] = sdata
                        yield AgentEvent("structure", text=stext, data=sdata)
                except Exception:  # noqa: BLE001 — 구조는 참고 근거일 뿐, 실패해도 결과를 막지 않는다
                    pass
                yield AgentEvent("optimum", data=final_pick)
        except Exception as exc:  # noqa: BLE001
            yield AgentEvent("error", f"Gemini 호출 오류: {type(exc).__name__}: {exc}")


# max_rounds 기본값 = 10 — 에이전트의 실제 종료는 `finish`(자기 수렴 선언)이고, 이 값은 그걸 못 할 때를
# 대비한 **안전 상한(budget ceiling)**이다. 근거:
#   · 수렴 바닥 6 — M4 효율 실측상 6라운드면 oracle 최적의 96.8% 도달(benchmark_agent.m4_agent_efficiency).
#   · 여유분 +4 — 가설 주도 설계(finish 전 1회 필수)·구조/개발성 도구 라운드 + CRITIC→REVISE 재설계
#     1~2회를 덮기 위함(grace +1은 마지막 스텝 REVISE만 커버).
#   · 상한 성격 — 보통 5~6라운드에 finish로 먼저 끝나 10에는 잘 안 닿음 → 비용 폭주 없는 백스톱.
# (주의: M4 벤치는 "타이트 예산에서도 효율적"임을 보이려 의도적으로 6을 쓴다 — 운영 기본값 10과 역할이 다름.)
def get_gemini_agent(settings: Settings, max_rounds: int = 10,
                     use_ledger: bool = True,
                     ledger_write: bool = True) -> GeminiOptimizationAgent | None:
    """use_ledger: 설계 이력 원장 사용. **벤치마크·애블레이션은 False** (비교 불가해짐).

    ledger_write: False 면 읽기만 한다 — 켬/끔 A/B 에서 두 팔이 **같은 바닥**을 보게 동결.
    """
    return (GeminiOptimizationAgent(settings, max_rounds, use_ledger, ledger_write)
            if settings.use_gemini_agent else None)
