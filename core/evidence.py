"""셔틀별 **문헌 전달 근거** — 순위의 근거를 우리 점수에서 인용으로 옮긴다.

왜 이 모듈이 있는가
-------------------
M6(`docs/experiments.md` §14)에서 **융합 수준 deepB3P 점수에는 서열 순서 신호가
측정되지 않는다**는 것이 확인됐다(경험 p 중위 0.392, 유의 셔틀 0/9, 문헌 음성 대조군
Angiopep-7 역전). 따라서 **"어느 셔틀이 더 좋은가"를 우리 점수로 말할 수 없다.**

그 주장의 정당한 출처는 문헌이다. 검증된 셔틀을 붙이면 뇌 전달이 올라간다는 것은
실험으로 확립된 사실이고, 우리는 그것을 **전제로 빌려 쓴다**. 이 모듈이 그 차용증이다.

⚠️ `fold` 를 왜 float 가 아니라 문자열로 두는가 (의도된 설계)
-------------------------------------------------------------
문헌의 "n배"는 **셔틀끼리 비교할 수 없다.** 전부 다르기 때문이다:
    · 종(랫/마우스) · 화물(저분자/단백질/mRNA/AAV) · 운반체(접합체/리포솜/LNP/나노입자)
    · 측정법(in situ perfusion Kin / %ID/g / 형질도입 효율 / 치료 효능)
    · **비교군**(유리 약물 대비냐, 무표적 운반체 대비냐 — 이게 가장 크게 다르다)
float 로 두면 누군가(사람이든 LLM이든) 반드시 정렬한다. 이 프로젝트에서 "측정 실패를
0 으로 인코딩"해 세 번 데인 것과 같은 함정의 반대쪽 얼굴이다 — **비교 불가능한 수치를
정렬 가능한 타입으로 인코딩하는 것.** 그래서 표시용 문자열로 못 박고, 순위에는
`tier`(근거 등급)만 쓴다.

근거 등급(tier)
---------------
    A — 1차 문헌의 정량 뇌 전달 증가 + 수용체 기전 확인 + 임상/전임상 진행
    B — 1차 문헌의 정량값 또는 포화성(=수용체 매개) 근거 있음. 비교군·측정계가 제한적
    C — 뇌 전달 보고는 있으나 **문헌의 고효율 형태가 라이브러리 서열과 다르거나**
        (D형·retro-enantio·분지형) 정량값을 확인하지 못함
    D — **BBB 통과 단독 근거 없음.** 다른 표적(뉴런 등)용 리간드이거나,
        문헌에서 BBB 투과 펩타이드와 **병용**을 전제로 쓰인다

⚠️ 등급은 **두 단계**로 결정된다(2026-11-11 추가) — ① 근거의 단단함 ② **간접성 하향**.
②는 GRADE 의 indirectness 를 가져온 것으로, **근거의 전달 포맷이 우리 구성체와
다르면 등급을 내린다.** 등급은 셔틀 속성으로 붙는데 순위가 매기는 것은 구성체라서,
그 사이에 검증되지 않은 외삽이 있었다. 라이브러리 RMT 9종 중 **유전자 융합 근거를
가진 것은 ApoE 하나**이고, 등급 A 를 만들던 Angiopep-2 의 `clinical` 은 **소분자
접합체**(ANG1005)에서 온 것이었다. 상세는 `derive_tier` 독스트링.

조사 시점: 2026-11-10 (간접성 축 2026-11-11). 전부 2차 확인(리뷰·초록·서지) 수준이며 원문 표/그림의 수치를
직접 대조한 것은 아니다. 수치를 논문에 쓰려면 1차 출처를 열어 확인해야 한다.
"""

from __future__ import annotations

from dataclasses import dataclass

from .config import SHUTTLES

# 등급 글자(보고용). ⚠️ **순위 계산에 직접 쓰지 말 것** — 아래 GROUP_ORDER 를 쓴다.
TIER_ORDER = ("A", "B", "C", "D")

# ── 순위에 실제로 쓰는 **거친 3단계** ───────────────────────────────────────
# 왜 글자 4단계를 그대로 쓰지 않는가: A/B/C/D 매핑은 1인·1세션 큐레이션이고,
# 특히 **A/B 경계가 가장 무르다.** 그런데 후보군(RMT) 9종 중 A 는 1종뿐이라
# 그 경계 하나가 "Angiopep-2 단독 1위" vs "여러 종 동급 경쟁"을 혼자 결정해 버린다.
# 즉 **가장 약한 판단이 결과를 단독으로 정한다.** 그래서 주장의 해상도를 증거의
# 해상도에 맞춘다 — 글자는 표시하되, 순위는 아래 3단계로만 가른다.
# (이것은 이 프로젝트가 세 번 데인 함정 — "비교 불가능한 값을 정렬 가능한 타입으로
#  인코딩" — 의 같은 얼굴이다. MIN_BBB_SIGNAL 폐지·`fold: str` 과 같은 이유다.)
GROUP_ORDER = ("근거 있음", "주의", "근거 없음")
#: 등급 글자 → 순위용 그룹
TIER_GROUP = {"A": "근거 있음", "B": "근거 있음", "C": "주의", "D": "근거 없음"}
#: 등급이 아예 없는 서열(라이브러리 밖이고 상속도 못 한 설계체)
GROUP_UNLISTED = "미등재"

# 설계체가 라이브러리 셔틀의 **등급을 상속**받기 위한 서열 유사도 하한.
# 임의 상수가 아니라 측정값이다 (M9, docs/experiments.md §18):
#   · 양성(pareto_scrub 개발성 변이체 55개) 부모 유사도 **최소 0.830**
#   · 음성 최대 — 스크램블 0.529 · 같은 길이 무작위 0.136 · homopolymer 0.429
#   · 타 라이브러리 셔틀 간 최대 **0.088** (교차 오상속 위험 사실상 없음)
#   · 최근접 라이브러리 셔틀이 부모가 아닌 양성: **0/55**
# → 간극 [0.529, 0.830] 안쪽. 음성 쪽으로 여유를 더 두어 0.70.
INHERIT_SIM_MIN = 0.70


@dataclass(frozen=True)
class Evidence:
    """한 셔틀의 문헌 전달 근거."""

    tier: str            # "A" | "B" | "C" | "D" — TIER_ORDER 참조
    fold: str | None     # 표시용 문자열. None = 정량값 미확보. ⚠️ 절대 float 로 바꾸지 말 것
    context: str         # 그 수치가 '무엇 대 무엇'인지 — 종·화물·운반체·측정법·비교군
    ref: str             # 1차 출처(확인 가능한 수준까지)
    caveat: str = ""     # 정직하게 같이 읽어야 하는 한계


@dataclass(frozen=True)
class Criteria:
    """등급을 **도출**하는 5개 감사 항목. 각 항목은 예/아니오이고 출처로 확인 가능하다.

    왜 이렇게 두는가: 등급을 **단언**하면 "내가 그렇게 봤다"가 근거가 된다. 그러면
    이견이 취향 다툼이 되고, 1인 큐레이션의 오류가 드러나지 않는다. 5개 불리언에서
    **계산**하면 이견이 "②가 예라고 본다"는 **확인 가능한 다툼**이 된다.

    판정 기준(엄격하게 적용할 것):
      quant      1차 문헌의 **정량 뇌 전달값**을 확인했나. 비교군이 명시돼야 하고,
                 효능 판독값(유전자 침묵률 등)·2차 인용만 있는 값은 **아니오**다.
      mech       수용체 매개 **기전이 실험으로 확인**됐나 — 차단 항체·경쟁 저해·
                 포화성·수용체 KO 중 하나. 단순히 "표적이 X다"는 **아니오**다.
      clinical   임상 단계에 진입했나(등록번호로 확인 가능).
      form_match 문헌에서 **효율이 입증된 형태**가 라이브러리 서열과 같은가.
                 D형·retro-enantio·분지형이 효능 형태인데 우리가 선형 L형이면 아니오.
      bbb_alone  **BBB 통과 단독 근거**가 있나. BBB를 우회하는 투여 경로(ICV 등)나
                 다른 셔틀과의 병용이 전제면 아니오.

    그리고 **간접성(indirectness)** — 위 5개와 성격이 다르다. 5개는 "근거가 얼마나
    단단한가"를 재고, 간접성은 "**그 근거가 우리 구성체에 대한 것인가**"를 잰다.
    GRADE 의 indirectness 영역을 그대로 가져왔다 (아래 `derive_tier` 주석 참조).
    """

    quant: bool
    mech: bool
    clinical: bool
    form_match: bool
    bbb_alone: bool
    #: 전달 포맷 간접성 단계. 우리 구성체는 **단일 폴리펩타이드 유전자 융합체**
    #: (화물 + 링커 + 셔틀, 셔틀 1벌)이다. 근거가 그것과 얼마나 먼 포맷에서 나왔나:
    #:   0 = 직접   — 유전자 융합 단백질(셔틀 1벌)
    #:   1 = 1단계  — 공유 접합체(소분자·약물) 또는 화물 없는 펩타이드 단독 수송
    #:   2 = 2단계  — **다가 나노캐리어**(나노입자·리포솜·LNP·폴리머솜·AAV·덴드리머)
    #:                또는 BBB 우회 투여 경로
    indirectness: int = 0
    why: str = ""        # 각 항목을 그렇게 판정한 근거 한 줄


def _base_tier(c: Criteria) -> str:
    """간접성 하향 **이전**의 등급 — 근거 자체가 얼마나 단단한가만 본다."""
    if not c.bbb_alone:
        return "D"                      # BBB 통과 단독 근거 없음
    if not c.form_match:
        return "C"                      # 문헌의 효능 형태가 우리 서열과 다름
    if c.quant and c.mech and c.clinical:
        return "A"
    if c.quant or c.mech:
        return "B"
    return "C"                          # 정량값·기전 모두 미확보


def derive_tier(c: Criteria) -> str:
    """감사 항목 → 등급 글자. **이 함수가 등급의 정의다**(위 표는 그 결과일 뿐).

    두 단계다: ① 근거의 단단함(`_base_tier`) → ② **간접성 하향**.

    왜 ②가 필요한가 (GRADE indirectness)
    -------------------------------------
    등급은 **셔틀의 속성**으로 붙는데, 순위가 매기는 것은 **구성체**(화물+링커+셔틀)다.
    그 사이에 아무도 검증하지 않은 외삽이 있었다 — 근거의 **전달 포맷**이 우리가 만드는
    것과 다르다. 실제로 라이브러리 RMT 9종의 근거 포맷을 열어 보면:

        · 유전자 융합 단백질 … ApoE(159-167)2 (IDUA 융합) — **1종뿐**
        · 공유 접합체/펩타이드 단독 … Angiopep-2(ANG1005=파클리탁셀 접합) · Leptin30 · SynB1
        · 다가 나노캐리어 … T7(DGL-PEG 벡터) · B6(PEG-PLA 입자) · RVG29(LNP) ·
                            CDX(리포솜) · THR(금나노입자·AAV) · Tet1(폴리머솜·ICV)

    즉 **최상위 등급 A 를 만드는 `clinical` 이 소분자 접합체에서 왔고**, 포맷이 일치하는
    근거를 가진 셔틀은 등급 B 하나였다. 나노캐리어 근거는 **표면에 셔틀을 여러 벌 제시한
    avidity** 에 의존하는데 우리 구성체는 셔틀이 **1벌**이라 기전 수준에서 다르다.

    이 분야도 이미 알고 있다 — BBB 셔틀 리뷰는 "**화물이 셔틀 효율에 영향을 주고 그 역도
    성립한다**"고 명시하고(Chem Soc Rev 2016), 항체에 셔틀을 유전자 융합한 최근 연구는
    **펩타이드 복제 수를 늘려도 BBB 수송이 개선되지 않고** 수율·안정성만 떨어졌다고
    보고한다(Mol Pharmaceutics 2026, cyclic SLS–trastuzumab). 포맷 간 이전은 방향도
    크기도 보장되지 않는다.

    GRADE 를 그대로 따른다
    ----------------------
    GRADE 는 목표 질문(target PICO)과 근거의 PICO 가 다를 때를 **간접성**으로 보고
    확실성 등급을 내린다. **몇 단계 내릴지에 대한 수치 규칙은 없다** — 판단이고, 대신
    그 판단을 **명시·감사 가능하게** 남기도록 요구한다. 그래서 여기서도 단계를 셔틀별로
    선언하고 `why` 에 이유를 적는다. (GRADE 는 "직접 근거를 찾다가 없는 것"과 "애초에
    간접 근거를 가져다 쓰는 것"을 구별해 **후자**를 하향 대상으로 본다 — 우리는 후자다.)

    왜 C 에서 멈추나 (floor)
    ------------------------
    D 는 확실성 단계가 아니라 **범위 실패**("BBB 통과 단독 근거가 없다")라서 의미가
    다르다. 간접성으로 D 까지 내리면 두 뜻이 섞인다. 그래서 하향은 **A→B→C 안에서만**
    하고 C 를 바닥으로 둔다. C 로 내려가면 순위 그룹이 「근거 있음」→「주의」로 바뀌므로
    (TIER_GROUP) **최종 후보 자격을 잃는 것과 같은 효과**다 — 충분히 센 결과다.
    """
    base = _base_tier(c)
    if base == "D" or c.indirectness <= 0:
        return base
    # A → B → C, C 가 바닥
    i = TIER_ORDER.index(base)
    return TIER_ORDER[min(i + c.indirectness, TIER_ORDER.index("C"))]


#: 셔틀 이름 → 근거. 키는 `config.SHUTTLES` 와 정확히 일치한다(아래 임포트 시 검증).
SHUTTLE_EVIDENCE: dict[str, Evidence] = {
    "Angiopep-2": Evidence(
        tier="B",
        fold="~9.5× / 86×",
        context="**~9.5×** = Angiopep-2 의 뇌 전달이 음성 대조군 **Angiopep-7**(K10R·K15R) 대비 "
                "우세한 배수 — 같은 실험 안의 **점 변이 쌍**이라 라이브러리에서 가장 해석이 깨끗한 "
                "수치다(Demeule 2008 JPET). **86×** = ANG1005(파클리탁셀 접합체)의 뇌 유입계수 Kin "
                "vs **유리 파클리탁셀**, 랫 in situ perfusion — 비교군이 유리 약물이므로 "
                "'무표적 운반체 대비'가 아니다.",
        ref="Demeule et al. 2008 J Pharmacol Exp Ther (DOI 10.1124/jpet.107.131318; "
            "Angiopep 계열·Angiopep-7 대조, ~9.5×); Demeule et al. 2008 J Neurochem "
            "(LRP1 매개·투과를 트랜스페린/락토페린/아비딘보다 높게 보고); "
            "ANG1005 Kin 86× (리뷰 인용, Oncotarget 2016). ANG1005는 임상 2상까지 진행.",
        caveat="라이브러리 중 근거가 가장 단단하다. 점 변이 음성 대조군 Angiopep-7(K10R·K15R)이 "
               "존재하며, 우리 효능 축은 그 쌍을 **역전**시킨다(M6 §14-2).",
    ),
    "ApoE(159-167)2": Evidence(
        tier="B",
        fold="3.6× / ~20×",
        context="IDUA(리소좀 효소) 융합체의 세포 pulse-chase 수송 3.6×; 간에서 발현시킨 융합체의 "
                "**뇌 IDUA 양**이 셔틀 없는 발현 대비 ~20×(마우스).",
        ref="Wang et al. 2013 (ApoE(159-167) 단량체 LAVYQAGAR, 18aa 탠덤 이량체, MPS I 맥락). "
            "세포 BBB 모델에서 이량체 > 단량체. Teixidó 그룹 2020: rAAV8 뇌 전달 증가 펩타이드 중 하나.",
        caveat="🔴 **이름 충돌 주의.** 일부 문헌은 'ApoE(159-167)2'를 **(LRKLRKRLL)₂**로 적는다 — "
               "성숙 단백질 141-149를 신호펩타이드 18aa 포함 전구체 번호로 센 것이다. 우리 서열은 "
               "**LAVYQAGAR 이량체**(성숙체 159-167, Wang 2013)로 전혀 다른 분자다. "
               "(LRKLRKRLL)₂는 강한 양이온성이라 deepB3P 아티팩트를 크게 받는다 — 혼동하면 "
               "근거와 점수가 둘 다 어긋난다.",
    ),
    "Leptin30": Evidence(
        tier="C",
        fold=None,
        context="렙틴 단편 스크리닝에서 **뇌:혈장 비가 최고**이고 전장 렙틴과 동등. 흡수가 "
                "**포화성**(=수용체 매개)이고, 방사능이 모세혈관보다 **뇌 실질**에 분포. "
                "61-90 구간을 절단하면 뇌 축적이 유의하게 감소(=서열 의존).",
        ref="렙틴 유래 뇌 흡수 단편 스크린(leptin30 = 잔기 61-90); "
            "Liu et al. 2010 Biomaterials 31:5246 (leptin30-PEG-DGL 뇌 표적 유전자 전달).",
        caveat="절대 정량값이 아니라 **상대 비교**만 확인됐다. 30aa로 라이브러리에서 가장 길어 "
               "개발성·합성 부담이 크다.",
    ),
    "T7": Evidence(
        tier="C",
        fold="2.17×",
        context="DGL-PEG-T7 벡터의 **in vivo 유전자 침묵 효율** vs T7 없는 벡터. "
                "⚠️ 뇌 흡수량이 아니라 **효능 판독값**이다.",
        ref="Kuang et al. 2013. 기전: TfR1의 트랜스페린과 **다른 에피토프**에 결합 → "
            "내인성 Tf와 비경쟁(라이브러리 note에 이미 반영).",
        caveat="뇌 축적량을 무표적 대조와 직접 비교한 수치를 확보하지 못했다. "
               "7aa로 가장 짧아 개발성은 유리하다.",
    ),
    "B6": Evidence(
        tier="C",
        fold=None,
        context="B6-PEG-PLA 나노입자가 무표적 입자보다 뇌 축적이 '훨씬 높음'(형광, 정성). "
                "AD 모델에서 **저용량으로도** 학습장애·콜린성 손상·해마 신경 손실 개선.",
        ref="Liu et al. 2013 Bioconjug Chem 24:997 (NAP 전달). hTfR 세포외도메인 "
            "9-mer 파지디스플레이 유래. B6 vs T12 vs T7 빈크리스틴 리포솜 head-to-head 존재"
            "(PubMed 34933080) — **수치 미확인, 원문 확인 필요**.",
        caveat="정량 비교값을 확보하지 못했다. head-to-head 논문이 있으므로 "
               "라이브러리 내 TfR 셔틀(T7·B6·THR) 서열 간 근거 비교는 그 논문으로 메울 수 있다.",
    ),
    "RVG29": Evidence(
        tier="C",
        fold="~70×",
        context="RVG29-LNP의 **뇌 루시페라제 mRNA 형질도입** vs 무표적 LNP, 정맥 투여, "
                "건강 성체 마우스. (2025 리뷰의 2차 인용 — 원 논문 미확인.)",
        ref="Tang·Han·Mitchell 2025 리뷰(Drug Deliv Transl Res)가 인용; "
            "Han et al. Nano Letters(펩타이드 기능화 LNP 뇌 mRNA 전달). "
            "siRNA 계열은 Kumar 그룹 원 논문 계보.",
        caveat="🔴 **효능으로 이어지지 않는다는 보고가 있다.** PLGA 나노입자 연구에서 겉보기 "
               "뇌 전달은 증가했지만 **직교소 뇌종양 모델에서 치료 효능 이득 없음**이었고, "
               "RVG29 리간드의 **수용액 불안정성** 때문에 증가가 지속되지 않았다. "
               "면역원성·단백질분해 민감성도 보고됐다. 29aa에 Cys 1개 → 이량화 위험.",
    ),
    "CDX": Evidence(
        tier="C",
        fold=None,
        context="nAChR 결합 IC50 84.5 nM. **D-retro-inverso 형(DCDX)**이 L-모체보다 in vitro BBB "
                "단층 transcytosis가 유의하게 높고, DCDX-PEG 리포솜의 in vitro transcytosis "
                "~7.9%(α-bungarotoxin으로 차단 = 수용체 특이).",
        ref="Wei et al. 2015 Angew Chem Int Ed(DCDX); candoxin(코브라 독) 유래 16-mer.",
        caveat="🔴 **문헌에서 작동하는 형태는 D형(DCDX)이다.** 라이브러리의 L형 CDX는 "
               "리소좀·혈청 안정성이 낮은 쪽이고, L형 단독의 정량 뇌 전달값을 확보하지 못했다.",
    ),
    "THR": Evidence(
        tier="C",
        fold="2.6×",
        context="셔틀 2개를 단 **분지형 retro-enantio(THRre)**가 BBB 세포 모델에서 모델 단백질 수송 2.6×. "
                "선형 THR 자체는 금나노입자를 in vitro·in vivo 로 BBB 통과시켰고 AAV8 뇌 전달도 "
                "유의하게 증가시켰다(비리온 직접 결합 기전).",
        ref="Teixidó/Giralt 그룹(THR·THRre·분지형 BBB 셔틀); BBB 셔틀 펩타이드가 전신 투여 후 "
            "AAV 뇌 전달을 높인다는 2018 PMC5999585 계열 연구.",
        caveat="🔴 **선형 L-THR은 혈청 프로테아제 안정성이 낮다**(라이브러리 note에 이미 반영). "
               "문헌의 고효율 형태는 **retro-enantio·분지형**이고 우리 라이브러리는 선형 L형이다. "
               "2.6×도 그 변형체의 세포 모델 수치다.",
    ),
    "Tet1": Evidence(
        tier="D",
        fold=None,
        context="GT1b 갱글리오사이드·스핑고미엘린 결합 **뉴런 표적** 리간드. 뇌실내(ICV) 투여 시 "
                "형질도입 뉴런 58% vs 무표적 벡터 13% — 단 ICV는 **BBB를 우회하는 경로**다. "
                "Tf/Tet1 이중리간드 폴리머솜의 BBB PS 0.28 mL/h/g·0.08 %ID/g는 **복합 구성체** "
                "값이라 Tet1 단독에 귀속할 수 없다.",
        ref="파상풍독소 중쇄 C단편 모사 12-mer 파지디스플레이; Tf/Tet1 폴리머솜 BBB 수송 연구; "
            "B6/Tet1·TGN-Tet1(TPL) 이중리간드 연구.",
        caveat="🔴 **BBB 셔틀로서 단독 근거가 없다.** Tet1-LNP가 뇌가 아니라 **말초 신경**으로 "
               "주로 간다는 보고가 있고, 문헌에서 Tet1은 거의 항상 BBB 투과 펩타이드(B6·TGN)와 "
               "**짝지어** 쓰인다 — Tet1이 BBB를 통과시키는 역할이 아니라, 통과 후 **뉴런을 고르는** "
               "역할이기 때문이다.\n"
               "⚠️ **벤치마크 함의**: M2·M4 의 후보군 oracle 최적이 `직접융합+Tet1`(융합 BBB 0.568)"
               "이다 — 즉 **우리 효능 축의 1위가 라이브러리에서 BBB 근거가 가장 약한 셔틀**이다. "
               "M6 와 완전히 독립된 경로(통계가 아니라 문헌)에서 같은 결론에 도달한다.",
    ),
    "TAT": Evidence(
        tier="C",
        fold=None,
        context="AuNP-TAT 뇌 흡수 0.02–0.06 %ID/g(낮음). TAT-ANG 융합체가 추가 뇌 추출 분석에서 "
                "BBB 통과에 실패했다는 보고도 있다. 비특이 막투과이므로 뇌 **선택성**이 없다.",
        ref="HIV-1 TAT 유래 13-mer CPP. 뇌 전달 맥락의 정량 perfusion 쌍(TAT 단독 vs 접합체, "
            "같은 실험 내)을 확보하지 못했다.",
        caveat="**CPP = 비교군.** 에이전트의 최종 선택 대상이 아니다(비특이 흡수 기저선). "
               "deepB3P가 양이온성을 고평가하므로 점수만 보면 상위에 온다 — 그게 이 축을 "
               "최적화 목표로 쓸 수 없는 이유다.",
    ),
    "Penetratin": Evidence(
        tier="C",
        fold="~6× / ~20×",
        context="**D-penetratin** 접합 독소루비신의 뇌 흡수 ~6×, 모세혈관 고갈법으로 본 "
                "**뇌 실질** 도달량 ~20×(랫 in situ perfusion). 저자들은 P-gp 유출 회피로 해석.",
        ref="Rousselle et al. 2000 Mol Pharmacol 57:679.",
        caveat="그 수치는 **D형**이고 라이브러리는 L형 Penetratin이다. **CPP = 비교군.** "
               "M6에서 W48F 점 변이 쌍(야생형 > W48F)을 셔틀 단독으로는 구분 못 했다"
               "(0.985 vs 0.986).",
    ),
    "SynB1": Evidence(
        tier="C",
        fold="~6× / ~20× / ~30×",
        context="SynB1 접합 독소루비신 뇌 흡수 ~6×, 모세혈관 고갈 기준 실질 ~20×(랫 in situ "
                "perfusion); 마우스 perfusion 연구에서는 ~30×. 두 수치는 종·설계가 달라 "
                "직접 비교 불가.",
        ref="Rousselle et al. 2000 Mol Pharmacol 57:679(랫, ~6×/~20×); "
            "Rousselle et al. 2001 J Pharmacol Exp Ther(마우스, L-SynB1·L-SynB3·D-SynB3 각 ~30×).",
        caveat="증가분의 상당 부분이 **P-gp 유출 회피**로 설명된다 — 즉 '뇌로 가는 능력'과 "
               "'퍼내지지 않는 능력'이 섞여 있다. **CPP = 비교군.** M2 전체 oracle "
               "(직접융합+SynB1, 0.960)이 바로 이 비교군이다.",
    ),
}


#: 셔틀 이름 → 등급 **감사 항목**. `SHUTTLE_EVIDENCE[n].tier` 는 여기서 도출되며,
#: 둘이 어긋나면 임포트가 실패한다(아래 검증). 등급을 바꾸고 싶으면 **여기를** 고칠 것.
SHUTTLE_CRITERIA: dict[str, Criteria] = {
    "Angiopep-2": Criteria(
        quant=True, mech=True, clinical=True, form_match=True, bbb_alone=True,
        indirectness=1,
        why="정량: Demeule 2008 JPET 의 Angiopep-7 점변이쌍 대비 ~9.5×(같은 실험 내 대조). "
            "기전: LRP1 매개 확인. 임상: ANG1005 2상. 형태: 선형 L 19aa 그대로."
        "간접성 1: 임상(ANG1005)은 **파클리탁셀 접합체**이고 ~9.5× 는 화물 없는 **펩타이드 단독** 수송이다. 단일 복제 유전자 융합 근거는 없다.",
    ),
    "ApoE(159-167)2": Criteria(
        quant=True, mech=False, clinical=False, form_match=True, bbb_alone=True,
        indirectness=0,
        why="정량: IDUA 융합 3.6×(세포 pulse-chase)·뇌 IDUA ~20×(마우스). "
            "기전: 이량체>단량체는 보였으나 수용체 차단·경쟁 실험을 확인하지 못했다. "
            "형태: LAVYQAGAR 탠덤 이량체 = 라이브러리 서열(Wang 2013)."
        "간접성 0: IDUA **유전자 융합 단백질** = 우리 구성체와 같은 포맷.",
    ),
    "Leptin30": Criteria(
        quant=False, mech=True, clinical=False, form_match=True, bbb_alone=True,
        indirectness=1,
        why="정량: 절대값 없음 — 뇌:혈장 비의 **상대 비교**만. "
            "기전: 흡수가 **포화성**이고 61-90 절단 시 뇌 축적 유의 감소(서열 의존) → 수용체 매개."
        "간접성 1: 뇌:혈장 비·포화성은 **방사표지 펩타이드 단독** 수치이고, 유전자 전달은 PEG-PLL 덴드리그라프트(다가)였다.",
    ),
    "T7": Criteria(
        quant=False, mech=True, clinical=False, form_match=True, bbb_alone=True,
        indirectness=2,
        why="정량: 2.17× 는 **유전자 침묵 효능 판독값**이고 뇌 흡수량이 아니다 → 아니오. "
            "기전: TfR1 의 Tf 와 **다른 에피토프**에 결합(내인성 Tf 비경쟁) 확인."
        "간접성 2: 2.17× 는 **DGL-PEG-T7 다가 나노벡터**의 효능 판독값이다. 나노입자 표면 다가 제시의 avidity 에 의존 → 셔틀 1벌 융합으로 이전 미보장.",
    ),
    "B6": Criteria(
        quant=False, mech=False, clinical=False, form_match=True, bbb_alone=True,
        indirectness=2,
        why="정량: 무표적 입자 대비 '훨씬 높음'(형광 **정성**)뿐. 기전: hTfR 파지디스플레이 "
            "유래라 표적은 알지만 차단·경쟁 실험 미확인. head-to-head 원문(PubMed 34933080) "
            "**미확보** — 이 두 칸이 채워지면 B 로 올라간다."
        "간접성 2: **B6-PEG-PLA 나노입자** 뇌 축적.",
    ),
    "RVG29": Criteria(
        quant=False, mech=False, clinical=False, form_match=True, bbb_alone=True,
        indirectness=2,
        why="정량: ~70× 는 2025 리뷰의 **2차 인용**이고 원 논문을 확인하지 못했다 → 아니오. "
            "기전: nAChR 표적이라고만 되어 있고 차단 실험 미확인. "
            "⚠️ 별도로 '겉보기 전달은 늘었으나 치료 효능 이득 없음' 보고가 있다."
        "간접성 2: **RVG29-LNP** mRNA 형질도입.",
    ),
    "CDX": Criteria(
        quant=False, mech=True, clinical=False, form_match=False, bbb_alone=True,
        indirectness=2,
        why="기전: α-bungarotoxin 으로 **차단**됨 = 수용체 특이 확인. "
            "형태: 문헌에서 작동하는 형태는 **D-retro-inverso(DCDX)** 이고 라이브러리는 선형 L형 "
            "→ form_match 아니오이므로 기전이 확인돼도 C."
        "간접성 2: **DCDX-PEG 리포솜** transcytosis (게다가 D형).",
    ),
    "THR": Criteria(
        quant=False, mech=False, clinical=False, form_match=False, bbb_alone=True,
        indirectness=2,
        why="정량: 2.6× 는 **분지형 retro-enantio(THRre)** 의 세포 모델 수치로 우리 형태가 아니다. "
            "형태: 고효율 형태가 retro-enantio·분지형 → 아니오. 선형 L-THR 은 혈청 안정성도 낮다."
        "간접성 2: 분지형 THRre 세포 모델 / 선형 THR 은 **금나노입자·AAV** 비리온 결합.",
    ),
    "Tet1": Criteria(
        quant=False, mech=False, clinical=False, form_match=True, bbb_alone=False,
        indirectness=2,
        why="bbb_alone **아니오**: 형질도입 58% 수치는 **ICV(뇌실내)** 투여로 BBB 를 우회한 경로이고, "
            "문헌에서 Tet1 은 거의 항상 BBB 투과 펩타이드(B6·TGN)와 **병용**된다. Tet1 의 역할은 "
            "BBB 통과가 아니라 통과 후 **뉴런 선택**이다 → 등급 D."
        "간접성 2: **폴리머솜** 복합 구성체 + ICV 우회 경로. (이미 bbb_alone 아니오로 D)",
    ),
    "TAT": Criteria(
        quant=False, mech=False, clinical=False, form_match=True, bbb_alone=True,
        indirectness=2,
        why="정량: AuNP-TAT 0.02-0.06 %ID/g 로 낮고 같은 실험 내 정량 쌍을 확보하지 못했다. "
            "기전: **비특이 막투과**라 수용체 기전이 애초에 없다. (CPP = 비교군)"
        "간접성 2: **AuNP-TAT** 입자. (CPP = 비교군)",
    ),
    "Penetratin": Criteria(
        quant=False, mech=False, clinical=False, form_match=False, bbb_alone=True,
        indirectness=1,
        why="형태: ~6×/~20× 는 **D-penetratin** 수치이고 라이브러리는 L형 → 아니오. (CPP = 비교군)"
        "간접성 1: **독소루비신 공유 접합체**. (CPP = 비교군)",
    ),
    "SynB1": Criteria(
        quant=True, mech=False, clinical=False, form_match=True, bbb_alone=True,
        indirectness=1,
        why="정량: **L-SynB1** 자체의 랫 ~6×/실질 ~20×(Rousselle 2000)·마우스 ~30×(2001) → 예. "
            "기전: 증가분이 **P-gp 유출 회피**로 설명돼 수용체 매개가 아니다. "
            "⚠️ 등급은 '근거의 단단함'이지 '선택해도 된다'가 아니다 — **CPP 는 비교군**이라 "
            "최종 후보 자격은 기전 분리(is_comparison_only)가 따로 막는다."
        "간접성 1: **독소루비신 공유 접합체** — 단백질 융합 근거 아님. (CPP = 비교군)",
    ),
}


# --- 등급은 **도출값**이다: 선언(Evidence.tier)과 감사 항목이 어긋나면 임포트 실패 ---
_cmiss = set(SHUTTLES) - set(SHUTTLE_CRITERIA)
_cextra = set(SHUTTLE_CRITERIA) - set(SHUTTLES)
if _cmiss or _cextra:
    raise RuntimeError(
        f"SHUTTLE_CRITERIA 가 SHUTTLES 와 어긋난다 — 누락 {sorted(_cmiss)} / 잉여 {sorted(_cextra)}")
_drift = {n: (SHUTTLE_EVIDENCE[n].tier, derive_tier(SHUTTLE_CRITERIA[n]))
          for n in SHUTTLES
          if SHUTTLE_EVIDENCE[n].tier != derive_tier(SHUTTLE_CRITERIA[n])}
if _drift:
    raise RuntimeError(
        "등급 글자가 감사 항목에서 도출되지 않는다(선언, 도출): "
        + ", ".join(f"{n}{v}" for n, v in sorted(_drift.items()))
        + " — 등급은 SHUTTLE_CRITERIA 를 고쳐서 바꿀 것이고, 글자만 손으로 바꾸지 말 것.")
del _cmiss, _cextra, _drift


# --- 드리프트 방지: 두 딕셔너리는 전부 정적 리터럴이므로 임포트 시 검사해도 안전하다 ---
_missing = set(SHUTTLES) - set(SHUTTLE_EVIDENCE)
_extra = set(SHUTTLE_EVIDENCE) - set(SHUTTLES)
if _missing or _extra:  # pragma: no cover - 편집 드리프트 즉시 적발
    raise RuntimeError(
        f"core.evidence 와 config.SHUTTLES 키 불일치 — 근거 누락 {sorted(_missing)} / "
        f"라이브러리에 없는 항목 {sorted(_extra)}"
    )
_bad = {n: e.tier for n, e in SHUTTLE_EVIDENCE.items() if e.tier not in TIER_ORDER}
if _bad:  # pragma: no cover
    raise RuntimeError(f"core.evidence 알 수 없는 근거 등급: {_bad}")
del _missing, _extra, _bad


def evidence_for(shuttle_name: str) -> Evidence | None:
    """셔틀 **이름**으로 근거를 찾는다(서열이 아니라 이름 — 행의 `shuttle` 키가 이름이다)."""
    return SHUTTLE_EVIDENCE.get(shuttle_name)


def tier_of(shuttle_name: str) -> str | None:
    e = SHUTTLE_EVIDENCE.get(shuttle_name)
    return e.tier if e else None


def tier_rank(shuttle_name: str) -> int:
    """정렬용 순위(작을수록 근거가 강함). 미등재는 맨 뒤."""
    t = tier_of(shuttle_name)
    return TIER_ORDER.index(t) if t in TIER_ORDER else len(TIER_ORDER)


def prompt_brief(max_chars: int = 2600) -> str:
    """에이전트 프롬프트에 넣을 근거 요약.

    목적은 "어느 셔틀이 더 좋은가"의 근거를 **우리 점수에서 인용으로** 옮기는 것이다.
    등급만 쓰고 배수는 넣지 않는다 — 배수는 셔틀끼리 비교할 수 없기 때문이다.
    """
    lines = [
        "[문헌 전달 근거 등급 — 순위를 정당화할 때 쓸 유일한 근거]",
        "우리 효능 점수(deepB3P)는 **융합 수준에서 셔틀 간 순위를 매길 수 없다**(검증 M6: "
        "문헌 검증 셔틀이 자기 스크램블 분포의 39번째 백분위, 유의 0/9, 음성 대조군 역전). "
        "따라서 **'이 셔틀이 더 좋다'는 주장은 점수가 아니라 아래 근거 등급으로** 하십시오.",
        "  A = 정량 뇌 전달 + 기전 확인 + 임상 진행 · B = 정량값 또는 포화성 근거 있음",
        "  C = 문헌의 고효율 형태가 이 서열과 다름(D형·분지형) 또는 정량값 미확보",
        "  D = **BBB 통과 단독 근거 없음** (다른 표적용 리간드이거나 병용 전제)",
        "  ⚠️ 여기에 **간접성 하향**이 더해진다 — 근거의 전달 포맷이 우리 구성체"
        "(단일 폴리펩타이드 유전자 융합, 셔틀 1벌)와 멀면 등급을 내린다. "
        "다가 나노캐리어(나노입자·리포솜·LNP·AAV) 근거는 −2단계, 공유 접합체·펩타이드 "
        "단독 근거는 −1단계다. 그래서 **현재 등급 A 는 없다** — 유전자 융합 포맷의 "
        "근거를 가진 셔틀은 ApoE(159-167)2 하나이고, Angiopep-2 의 임상 근거는 "
        "파클리탁셀 접합체(ANG1005)다. 포맷이 다른 근거를 '임상 통과'로 내세우지 마십시오.",
        "",
    ]
    for tier in TIER_ORDER:
        names = [n for n in SHUTTLES if SHUTTLE_EVIDENCE[n].tier == tier]
        if names:
            lines.append(f"  · 등급 {tier}: {', '.join(names)}")
    lines += [
        "",
        "주의할 개별 항목:",
        "  · **Tet1 = 등급 D**: 뉴런 표적(GT1b) 리간드이고 BBB 통과 근거가 없다. 점수가 높게 "
        "나와도 그것만으로 최종안으로 고르지 말고, 고른다면 한계로 **명시**하십시오.",
        "  · **CDX·THR = 등급 C**: 문헌의 고효율 형태는 D형/retro-enantio/분지형이고 "
        "라이브러리는 선형 L형이다. THR은 혈청 안정성도 낮다.",
        "  · **RVG29**: 뇌 전달 증가 보고는 크지만 **치료 효능 이득이 없었다는 보고**와 "
        "리간드 수용액 불안정성이 함께 있다.",
        "  · **TAT·Penetratin·SynB1 = CPP 비교군**: 선택 대상이 아니다.",
    ]
    out = "\n".join(lines)
    return out if len(out) <= max_chars else out[:max_chars] + " …"


def report() -> str:
    """사람이 읽는 전체 표(문서·검토용)."""
    rows = sorted(SHUTTLES, key=lambda n: (tier_rank(n), n))
    out = ["| 셔틀 | 등급 | 보고된 증가 | 측정 맥락(비교군) | 출처 |",
           "|---|---|---|---|---|"]
    for n in rows:
        e = SHUTTLE_EVIDENCE[n]
        ctx = e.context.replace("\n", " ")
        out.append(f"| {n} | **{e.tier}** | {e.fold or '—'} | {ctx} | {e.ref} |")
    out.append("")
    out.append("⚠️ 배수는 **셔틀끼리 비교할 수 없다** — 종·화물·운반체·측정법·비교군이 모두 다르다.")
    return "\n".join(out)


if __name__ == "__main__":  # pragma: no cover
    print(report())
    print()
    print(prompt_brief())


def names_by_tier(tier: str) -> list[str]:
    """해당 근거 등급의 셔틀 이름들(라이브러리 선언 순서)."""
    return [n for n in SHUTTLES if SHUTTLE_EVIDENCE[n].tier == tier]


def critic_rule() -> str:
    """심사(Critic) 프롬프트용 압축 규칙 — 설계자용 `prompt_brief` 보다 짧게."""
    d = ", ".join(names_by_tier("D")) or "없음"
    c = ", ".join(names_by_tier("C")) or "없음"
    ab = ", ".join(names_by_tier("A") + names_by_tier("B")) or "없음"
    return (
        "[문헌 근거 등급 — 셔틀 순위 주장의 유일한 근거]\n"
        "  융합 BBB 점수는 **셔틀 간 순위를 매길 수 없다**(검증 M6: 문헌 검증 셔틀이 자기 스크램블 "
        "분포의 39번째 백분위, 유의 0/9, 음성 대조군 역전). 따라서 '이 셔틀이 더 좋다'는 주장을 "
        "**점수로 하면 반박하십시오.** 근거는 아래 등급뿐입니다.\n"
        f"   · 근거 있음(A·B): {ab}\n"
        f"   · 주의(C — 문헌의 고효율 형태가 이 서열과 다름·정량값 미확보): {c}\n"
        f"   · **근거 없음(D — BBB 통과 단독 근거 없음)**: {d}\n"
        "  최종 후보가 **등급 D**라면: 게이트를 통과한 A·B 대안이 하나라도 있으면 그쪽을 지목하고, "
        "없을 때만 조건부 승인하되 **'BBB 셔틀로서 문헌 근거가 없다'를 한계로 반드시 명시**하게 "
        "하십시오. 등급은 탈락 사유가 아니라 **명시 사유**입니다(하드 게이트는 독성·효능 신호 둘뿐).\n\n"
    )


# ── 서열 → 등급 (행의 `shuttle` 키는 **서열**이라 이름 조회로는 닿지 않는다) ──────
#: 라이브러리 셔틀 서열 → 이름. 동일 서열이 둘이면 먼저 선언된 쪽.
_SEQ_TO_NAME: dict[str, str] = {}
for _n, _v in SHUTTLES.items():
    _SEQ_TO_NAME.setdefault(_v["seq"].upper(), _n)
del _n, _v


def _clean(seq: str) -> str:
    return "".join(ch for ch in (seq or "").upper() if ch.isalpha())


def name_of_seq(seq: str) -> str | None:
    """라이브러리 셔틀 **서열**로 이름을 찾는다(정확 일치만)."""
    return _SEQ_TO_NAME.get(_clean(seq))


def similarity(seq: str, ref: str) -> float:
    """정규화 서열 유사도 0~1. `binding.shuttle_similarity` 와 **같은 정의**(BLOSUM62 전역정렬)."""
    a, b = _clean(seq), _clean(ref)
    if not a or not b:
        return 0.0
    if a == b:
        return 1.0
    try:
        from .binding import _aligner
        aln = _aligner()
        selfs = aln.score(b, b)
        return max(0.0, min(1.0, aln.score(a, b) / selfs)) if selfs > 0 else 0.0
    except Exception:  # noqa: BLE001 — 정렬기 없으면 상속을 주지 않는다(보수적)
        return 0.0


def inherit(seq: str) -> tuple[str | None, str | None, float]:
    """이 셔틀 서열의 **근거 등급**을 구한다 → (등급, 부모 셔틀 이름, 유사도).

    · 라이브러리 서열이면 자기 등급 그대로 (유사도 1.0).
    · 라이브러리 밖(설계체)이면 **가장 닮은 라이브러리 셔틀의 등급을 상속**한다.
      단 유사도가 `INHERIT_SIM_MIN` 이상일 때만 — 아니면 (None, 최근접, 유사도).

    왜 상속이 필요한가: 등급은 **문헌에서 빌려온 유일한 양(陽)의 신호**이고 설계체에는
    인용이 없다. 그러면 모든 설계체가 등급 0 이 되어 **개발성 liability 스크러빙처럼
    정당한 기여까지 후보 자격을 잃는다**(그 편집은 효능을 주장하지 않고 '만들 수 있게'
    다듬는 것이므로 부모의 근거를 그대로 물려받는 것이 맞다).

    왜 유사도 하한이 필요한가: 하한이 없으면 아무 서열이나 '가장 닮은' 셔틀의 등급을
    물려받아 **등급이 다시 해킹 가능한 축이 된다.** 하한은 측정값이다 — M9 실측에서
    정당한 개발성 변이체(55개)의 부모 유사도 최소 0.830 vs 스크램블 최대 0.529 ·
    homopolymer 0.429 · 무작위 0.136 · 타 셔틀 간 0.088 로 간극이 비어 있다.
    """
    s = _clean(seq)
    if not s:
        return None, None, 0.0
    n = _SEQ_TO_NAME.get(s)
    if n is not None:
        return SHUTTLE_EVIDENCE[n].tier, n, 1.0
    best_n, best_s = None, 0.0
    for nm, v in SHUTTLES.items():
        sc = similarity(s, v["seq"])
        if sc > best_s:
            best_n, best_s = nm, sc
    if best_n is None:
        return None, None, 0.0
    if best_s >= INHERIT_SIM_MIN:
        if not _pins_intact(s, best_n):
            # 🔴 문헌이 기능 필수라 말한 자리가 깨졌다 → **상속 거부** (§31).
            return None, best_n, round(best_s, 3)
        return SHUTTLE_EVIDENCE[best_n].tier, best_n, round(best_s, 3)
    return None, best_n, round(best_s, 3)


def _pins_intact(seq: str, parent_name: str) -> bool:
    """부모의 **문헌 고정 위치**가 이 서열에 그대로 남아 있나 (§31).

    🔴 왜 유사도만으로는 안 되나 — 측정으로 닫힌 질문이다.
    옛 `inherit()` 은 유사도 하나만 봤다. 그래서:

      · **Angiopep-7**(`K10R/K15R`, 문헌상 뇌 전달 ~9.5배 낮은 **알려진 기능 상실체**)이
        유사도 **0.9434** 로 원본과 **같은 등급 B** 를 받았다.
      · K10·K15 를 파괴하는 변이 **399개를 전수 열거하니 399/399(100%)가 등급을 받았다**
        (유사도 0.849~0.972 — 임계 0.70 근처에도 가지 않는다).
      · **임계값으로는 고칠 수 없다.** 핀 파괴 변이체의 유사도 중위(0.887)가 정당한
        개발성 변이체의 중위(0.858)**보다 높다** — 분포가 겹치는 게 아니라 **역전**이다.
        Angiopep-7(0.943)을 막으려면 임계를 0.943 초과로 올려야 하고 그러면 정당 변이체
        25/25 가 **전부** 탈락한다.

    원인은 명확하다 — **K→R 은 BLOSUM62 +2 의 보존적 치환**이라 유사도가 거의 안 떨어진다.
    `CONSERVED_POSITIONS` 독스트링이 이미 적어 둔 결론인데 **상속 경로에 반영되지 않았다.**
    그래서 유사도를 손대는 대신 **위치를 직접 확인**한다.

    ⚠️ **길이가 바뀌면 거부한다.** 길이가 다르면 위치 대응이 깨져 핀을 **검증할 수 없다**.
       옛 동작은 `conserved_positions()` 가 `()` 를 돌려주는 바람에 가드가 사라졌고,
       단일 삽입 400개 전수에서 **400/400 이 등급을 유지하면서 가드 0%** 였다 —
       **삽입 1개로 우회**할 수 있었다. 최상위 근거 그룹에서 '검증 불가'는 '통과'가 아니다.

    ⚠️ **실효 범위는 1/12 이다.** `CONSERVED_POSITIONS` 에 문헌 위치가 있는 셔틀은
       **Angiopep-2 하나**뿐이다. 나머지 11종에 대해 이 검사는 **아무것도 걸러내지 못한다**.
       구조는 고쳤지만 데이터가 없다 — 과대 주장 금지(§31-6).
    ⚠️ 문헌이 시험한 것은 **이중 변이**(K10R/K15R)다. 한 자리만 바꾼 경우의 직접 근거는
       없지만, 동결 결정은 보수적으로 **자리마다** 적용한다(섞어 따지지 않는다).
    """
    pins = CONSERVED_POSITIONS.get(parent_name) or ()
    if not pins:
        return True                      # 문헌 위치가 없는 셔틀 → 검사할 것이 없다
    parent = _clean(SHUTTLES[parent_name]["seq"])
    if len(seq) != len(parent):
        return False                     # 위치 대응 불가 = 검증 불가 → 거부
    return all(1 <= p <= len(seq) and seq[p - 1] == parent[p - 1] for p in pins)


# ── 순위용 거친 3단계 ────────────────────────────────────────────────────────
def group_of(tier: str | None) -> str:
    """등급 글자 → 순위용 그룹. None(미등재)도 받는다."""
    return TIER_GROUP.get(tier or "", GROUP_UNLISTED)


def group_rank(tier: str | None) -> int:
    """정렬용 그룹 순위(작을수록 근거가 강함). 미등재는 맨 뒤.

    ⚠️ **등급 글자(A/B/C/D)로 직접 정렬하지 말고 이 함수를 쓸 것.** A/B 경계는
    1인 큐레이션에서 가장 무른 판단이고, 후보군 9종 중 A 가 1종뿐이라 그 경계가
    결과를 단독으로 결정해 버린다(GROUP_ORDER 주석 참고).
    """
    g = group_of(tier)
    return GROUP_ORDER.index(g) if g in GROUP_ORDER else len(GROUP_ORDER)


def rubric_report() -> str:
    """감사용 표 — 등급이 **어떤 항목에서 도출됐는지** 사람이 검토할 수 있게."""
    cols = [("quant", "정량"), ("mech", "기전"), ("clinical", "임상"),
            ("form_match", "형태일치"), ("bbb_alone", "BBB단독")]
    out = ["| 셔틀 | " + " | ".join(c[1] for c in cols)
           + " | 기본등급 | 간접성 | → 최종등급 | 순위 그룹 |",
           "|---|" + "---|" * (len(cols) + 4)]
    for n in sorted(SHUTTLES, key=lambda k: (tier_rank(k), k)):
        c = SHUTTLE_CRITERIA[n]
        cells = ["✅" if getattr(c, k) else "—" for k, _ in cols]
        base, t = _base_tier(c), derive_tier(c)
        ind = "직접" if c.indirectness == 0 else f"−{c.indirectness}단계"
        out.append(f"| {n} | " + " | ".join(cells)
                   + f" | {base} | {ind} | **{t}** | {group_of(t)} |")
    out += [
        "",
        "판정 규칙(`derive_tier`) — **두 단계**:",
        "① 기본등급: BBB단독 아니오 → **D** · 형태일치 아니오 → **C** · "
        "정량·기전·임상 전부 → **A** · 정량 또는 기전 → **B** · 그 외 → **C**",
        "② **간접성 하향**(GRADE indirectness): 근거의 **전달 포맷**이 우리 구성체"
        "(단일 폴리펩타이드 유전자 융합, 셔틀 1벌)와 멀수록 단계를 내린다 — "
        "0 직접(유전자 융합) / −1 공유 접합체·펩타이드 단독 / −2 다가 나노캐리어·BBB 우회. "
        "**C 가 바닥**이다(D 는 확실성 단계가 아니라 범위 실패라서 섞지 않는다).",
        "",
        "⚠️ 이 표는 **1인·1세션 큐레이션**이고 외부 검증을 받지 않았다. 각 칸은 출처로 "
        "확인 가능하도록 `SHUTTLE_CRITERIA[*].why` 에 판정 근거를 적어 두었다 — 이견이 "
        "있으면 **칸을 걸고** 다투면 된다(등급 글자를 손으로 바꾸면 임포트가 실패한다).",
    ]
    return "\n".join(out)


# ---------------------------------------------------------------------------
# 문헌이 **기능 필수**라고 지목한 잔기 위치 (1-based, 셔틀 서열 기준)
# ---------------------------------------------------------------------------
#: 왜 이게 필요한가 (2026-11-11, `docs/experiments.md` §25)
#:
#: `pareto_scrub` 의 기능 보존 제약은 `rmt_sim ≥ 0.85` — BLOSUM62 **전역 정렬 평균**이다.
#: 평균이라 **치명적 치환 2~3개가 동일 잔기 16개 뒤에 희석된다.** M10 에서 배운 병과
#: 같은 병이다("길이로 정규화된 지표는 희석된다 → 국소 핫스팟으로 교체"). 실측:
#:
#:   · 목적함수 3축일 때  Angiopep-2 변이체 165개 중 **20개(12.1%)** 가 K10/K15 를 바꿨다
#:   · off-target 을 목적에 넣자 254개 중 **47개(18.5%)**, 그중 **전하 역전 27건**
#:
#: 악화 메커니즘이 분명하다 — `off = (전하밀도 ∨ 소수성) × (1 − 0.6·rmt_sim)` 이므로
#: off 를 최소화하면 **양전하를 지우는 것이 보상**되고, 하필 K10·K15 가 라이신이다.
#: 즉 **옵티마이저가 문헌이 기능이라고 말한 잔기를 정확히 지우는 법을 배웠다.**
#:
#: 그리고 **BLOSUM 으로는 못 막는다**: 알려진 기능 상실 돌연변이가 `K10R/K15R`
#: (Angiopep-7)이고 K→R 은 BLOSUM62 **+2 = 보존적 치환**이다. 즉 *우리 옵티마이저가
#: 하도록 설계된 종류의 편집*이 그 기능 상실을 만든다. 유사도 축으로는 구별 불가이고
#: (README: "선택성·BLOSUM 등 모든 서열 축이 Angiopep-2 와 Angiopep-7 을 구별하지 못했다"),
#: ECHA read-across 가 화학 유사성만으로 풀지 않고 **기전 유사성**을 함께 요구하는 것과
#: 같은 이유다. 그래서 제약의 출처를 **유사도에서 문헌으로** 옮긴다.
#:
#: ⚠️ **문헌이 집어 준 곳만 적는다.** 비어 있는 셔틀은 "보존할 필요가 없다"가 아니라
#: **"어디가 중요한지 우리가 모른다"**는 뜻이다 — 이 간극 자체가 한계로 보고돼야 한다.
CONSERVED_POSITIONS: dict[str, tuple[int, ...]] = {
    # Angiopep-2 `TFFYGGSRGKRNNFKTEEY` — K10·K15.
    # 근거: **Angiopep-7 = K10R/K15R** 이고 문헌상 뇌 전달이 ~9.5배 낮다(Demeule 2008 JPET,
    # ⚠️ **2차 출처·원문 미대조**이며 같은 사실이 '뇌 실질 축적 ~85% 감소'(≈6.7배)로도 적혀
    #    있다. 두 수치가 문서 안에 공존한다 — 동결 결정은 **방향**(K10R/K15R 이 뚜렷히 열등)에
    #    의존하고 배수에 의존하지 않지만, 배수를 인용할 때는 원문을 대조해야 한다.
    # 같은 실험 내 점 변이 쌍 = 라이브러리에서 가장 깨끗한 음성 대조). 두 라이신을
    # **같은 양전하 잔기(R)로만 바꿔도** 기능이 무너진다면, 그 자리는 건드리지 않는다.
    "Angiopep-2": (10, 15),
}


def conserved_positions(shuttle_seq: str) -> tuple[int, ...]:
    """이 셔틀(또는 그 변이체)에서 **건드리면 안 되는 1-based 위치**.

    라이브러리 셔틀과 길이가 같을 때만 위치가 1:1 로 대응하므로, 길이가 다르면
    보수적으로 **빈 튜플**을 돌려준다(틀린 위치를 동결하는 것이 더 나쁘다).
    """
    nm = name_of_seq(shuttle_seq)
    if nm is None:
        nm, _, _sim = None, None, 0.0
        best, bname = 0.0, None
        for cand, ev in SHUTTLE_EVIDENCE.items():
            if cand not in CONSERVED_POSITIONS:
                continue
            sq = SHUTTLES[cand]["seq"] if isinstance(SHUTTLES[cand], dict) else SHUTTLES[cand]
            if len(sq) != len(shuttle_seq):
                continue
            sc = similarity(shuttle_seq, sq)
            if sc > best:
                best, bname = sc, cand
        nm = bname if best >= INHERIT_SIM_MIN else None
    if nm is None:
        return ()
    sq = SHUTTLES[nm]["seq"] if isinstance(SHUTTLES[nm], dict) else SHUTTLES[nm]
    if len(sq) != len(shuttle_seq):
        return ()
    return CONSERVED_POSITIONS.get(nm, ())
