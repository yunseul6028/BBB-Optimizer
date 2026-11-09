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

조사 시점: 2026-11-10. 전부 2차 확인(리뷰·초록·서지) 수준이며 원문 표/그림의 수치를
직접 대조한 것은 아니다. 수치를 논문에 쓰려면 1차 출처를 열어 확인해야 한다.
"""

from __future__ import annotations

from dataclasses import dataclass

from .config import SHUTTLES

# 등급이 높은 순. 순위 정당화에 쓸 수 있는 **유일한** 정렬 키다.
TIER_ORDER = ("A", "B", "C", "D")


@dataclass(frozen=True)
class Evidence:
    """한 셔틀의 문헌 전달 근거."""

    tier: str            # "A" | "B" | "C" | "D" — TIER_ORDER 참조
    fold: str | None     # 표시용 문자열. None = 정량값 미확보. ⚠️ 절대 float 로 바꾸지 말 것
    context: str         # 그 수치가 '무엇 대 무엇'인지 — 종·화물·운반체·측정법·비교군
    ref: str             # 1차 출처(확인 가능한 수준까지)
    caveat: str = ""     # 정직하게 같이 읽어야 하는 한계


#: 셔틀 이름 → 근거. 키는 `config.SHUTTLES` 와 정확히 일치한다(아래 임포트 시 검증).
SHUTTLE_EVIDENCE: dict[str, Evidence] = {
    "Angiopep-2": Evidence(
        tier="A",
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
        tier="B",
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
        tier="B",
        fold="2.17×",
        context="DGL-PEG-T7 벡터의 **in vivo 유전자 침묵 효율** vs T7 없는 벡터. "
                "⚠️ 뇌 흡수량이 아니라 **효능 판독값**이다.",
        ref="Kuang et al. 2013. 기전: TfR1의 트랜스페린과 **다른 에피토프**에 결합 → "
            "내인성 Tf와 비경쟁(라이브러리 note에 이미 반영).",
        caveat="뇌 축적량을 무표적 대조와 직접 비교한 수치를 확보하지 못했다. "
               "7aa로 가장 짧아 개발성은 유리하다.",
    ),
    "B6": Evidence(
        tier="B",
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
        tier="B",
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
