#!/usr/bin/env python3
"""수용체 수준 off-target 참조 — **외부 발현 데이터**로 서열 축을 대조한다 (§50).

왜 이 스크립트가 있는가
-----------------------
`core/selectivity.py` 의 `off_target` 축은 off-target 위험을 **서열에서만** 답한다:

    raw = 1 - (1 - cationicity) * (1 - lipophilicity)
    off = clamp(raw * (1 - 0.6 * rmt_sim))

이 식에는 *"RMT 면 선택적이다"* 라는 **전제**가 박혀 있다(`rmt_sim` 할인항). 그런데
수용체 매개 수송(RMT) 셔틀의 off-target 위험을 정하는 것은 셔틀의 전하가 아니라
**그 수용체가 어디에 발현하는가**다. 라이브러리 12종 중 **3종(T7·B6·THR)이 모두
TfR1(TFRC)** 을 노리는데, TFRC 는 적혈구 전구세포·골수·간·태반에 널리 발현하고
**TfR1 결합 자체가 급성 reticulocyte 감소의 알려진 원인**이다(아래 `LIABILITIES`).
즉 생물학이 축과 **반대**를 말할 수 있다.

이 스크립트는 **외부 대조군**이다. 축을 바꾸지 않는다 — `core/` 는 한 줄도 건드리지
않으며, 여기서 나온 결론은 §50 의 **권고**로만 남는다.

설계 원칙 — `core/evidence.py` 를 그대로 베낀다
-----------------------------------------------
`core/evidence.py` 가 문헌 근거를 다루는 방식(불리언 감사 항목 → 규칙으로 등급 도출 →
선언값과 어긋나면 **임포트 실패**)을 그대로 쓴다. 이유도 같다: 등급을 **단언**하면
"내가 그렇게 봤다"가 근거가 되고 이견이 취향 다툼이 된다. **계산**하면 이견이
"이 칸이 틀렸다"는 **확인 가능한 다툼**이 된다.

그래서 세 겹이다:
  ① `HpaRecord` — 외부 출처에서 **그대로 적어 온 수치**(조회일·출처 ID 포함)
  ② `criteria()` — ①에서 **계산되는** 5개 불리언 (손으로 쓰지 않는다)
  ③ `derive_tier()` — ②에서 **계산되는** 거친 3단계. `Receptor.tier` 선언값과
     어긋나면 임포트가 실패한다.

⚠️ 단계(tier)를 **연속 점수로 만들지 않는 이유**는 `core/evidence.py` 의 `fold: str`
   과 같다 — 비교 가능성이 보증되지 않은 값을 정렬 가능한 타입으로 인코딩하면
   누군가(사람이든 LLM이든) **반드시 정렬한다**. nTPM 은 조직 간 비교가 어느 정도
   되지만 "간 121.9 vs 골수 98.6" 이 off-target 위험의 비율이라는 보증은 **없다**.
   그래서 순위에는 **3단계**만 쓴다.

⚠️ **배수(fold-change)는 이 파일에 하나도 없다.** 문헌에서 가져온 말초 독성 근거는
   `Liability.direction` 에 **방향만** 적는다(`docs/protocol.md` §5 보고 규약).
   §40 에서 확인 불가한 배수(~9.5×)로 한 번 데였다.

출처 (전부 이 세션에서 실제로 조회함)
--------------------------------------
  HPA   Human Protein Atlas **v26.0**(릴리스 2026-11-11 · Ensembl 109),
        `https://www.proteinatlas.org/api/search_download.php` 공개 API,
        컬럼 `rnats/rnatd/rnatsm/rnascsm/rnabcsm/prts/prtd/relih` + `t_RNA_<조직>`.
        조회일 **2026-11-12**. 조직 분모는 HPA consensus tissue **47종** 중
        CNS 10종(amygdala·basal ganglia·cerebellum·cerebral cortex·choroid plexus·
        hippocampal formation·hypothalamus·midbrain·spinal cord·retina)을 뺀
        **37종**이다. (pituitary gland 는 BBB 밖이라 비-CNS 로 센다.)
  PM    PubMed/NCBI E-utilities 로 **초록 본문을 읽은** 1차 문헌. PMID 명시.

🔴 **읽지 못한 것은 적지 않았다.** 본문 표·그림을 대조한 논문은 **0편**이고, 아래
   `Liability` 의 `verified=True` 는 *"초록 본문을 이 세션에서 읽었다"* 까지만
   뜻한다. 정량값은 그래서 하나도 싣지 않는다.

실행
----
    /usr/bin/python3 scripts/receptor_expression.py
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import SHUTTLES            # noqa: E402  (읽기만 한다)
from core.selectivity import assess_selectivity  # noqa: E402  (실제 축 함수를 호출)

# ---------------------------------------------------------------------------
# 출처 명부
# ---------------------------------------------------------------------------
SOURCES: dict[str, str] = {
    "HPA26": "Human Protein Atlas v26.0 (릴리스 2026-11-11 · Ensembl 109) — "
             "consensus tissue RNA (nTPM) + Tissue Atlas 항체 기반 단백질 수준 주석. "
             "공개 API `proteinatlas.org/api/search_download.php`, 조회일 2026-11-12. "
             "HPA consensus 는 HPA 자체 RNA-seq 과 GTEx 를 통합한 데이터셋이다 "
             "(즉 GTEx 를 따로 조회하지 않고 HPA consensus 를 썼다).",
    "PM23636093": "Couch JA et al. 2013 Sci Transl Med 5(183):183ra57 "
                  "(DOI 10.1126/scitranslmed.3005338, PMID 23636093) — "
                  "\"Addressing safety liabilities of TfR bispecific antibodies that "
                  "cross the blood-brain barrier\". 초록 본문 조회(2026-11-12).",
    "PM41287279": "Smith BA et al. 2025 MAbs 17(1):2592422 (PMID 41287279) — "
                  "anti-TfR1 항체의 친화도·발렌시를 체계적으로 바꿔 뇌 흡수와 "
                  "**reticulocyte depletion** 을 함께 측정. 초록 본문 조회(2026-11-12).",
    "PM24411731": "Niewoehner J et al. 2014 Neuron 81(1):49-60 (PMID 24411731) — "
                  "monovalent TfR 결합(Brain Shuttle). 초록 본문 조회(2026-11-12).",
    "PM18626063": "Lillis AP et al. 2008 Physiol Rev 88(3):887-918 (PMID 18626063) — "
                  "\"LRP1 ... is a large endocytic receptor that is **widely expressed "
                  "in several tissues**\"; 조직 특이 KO 가 혈관·CNS·**대식세포**·"
                  "**지방세포**에서의 기여를 보인다. 초록 본문 조회(2026-11-12).",
}

#: HPA consensus tissue 중 CNS 로 셈하지 않는 조직 분모. 이 값이 바뀌면 모든 불리언이
#: 바뀌므로 **수치와 같은 줄에 조건을 적는다**(`docs/protocol.md` §5).
N_NONCNS_TISSUES = 37

#: 「광범위」판정의 조직 수 하한. 비-CNS 37조직의 **절반 이상**(≥20)에서 nTPM ≥ 10.
#: ⚠️ 이 두 숫자(20·10)는 **판단이고 외부 근거가 없다.** HPA 자신의 "Detected"
#:    하한이 nTPM ≥ 1 이라 그보다 보수적으로 잡았을 뿐이다. 밴드 폭 0.05 와 같은
#:    종류의 자유 파라미터이므로 §50 의 한계에 명시한다.
BROAD_MIN_TISSUES = 20
#: 「말초 핫스팟」판정의 nTPM 하한 — 비-CNS 조직 중 하나라도 이 값 이상.
HOTSPOT_NTPM = 30.0
#: 「조직 수」를 셀 때의 nTPM 하한.
PRESENT_NTPM = 10.0

# ---------------------------------------------------------------------------
# 거친 3단계 (+ 미결정)
# ---------------------------------------------------------------------------
#: 넓을수록 앞. 순위·비교에는 **이 문자열만** 쓴다(연속 점수를 만들지 않는다).
TIER_ORDER = ("광범위", "편중", "좁음")
TIER_BROAD, TIER_FOCAL, TIER_NARROW = TIER_ORDER
#: 수용체 라벨 자체를 조회할 수 없는 경우. 「좁음」이 **아니다** — "모른다"다.
TIER_UNKNOWN = "미결정"


@dataclass(frozen=True)
class HpaRecord:
    """HPA 에서 **그대로 적어 온** 외부 수치. 손으로 해석하지 않는다."""

    ensembl: str
    uniprot: str
    rna_specificity: str          # HPA "RNA tissue specificity"
    rna_distribution: str         # HPA "RNA tissue distribution"
    protein_specificity: str      # HPA "Protein tissue specificity"(항체 IH)
    protein_distribution: str     # HPA "Protein tissue distribution"
    reliability_ih: str | None    # HPA "Reliability (IH)" — 항체 신뢰도
    n_noncns_ge10: int            # 비-CNS 37조직 중 nTPM ≥ 10 인 조직 수
    n_noncns_ge30: int            # 같은 분모에서 nTPM ≥ 30 인 조직 수
    top_noncns: tuple             # (조직, nTPM) — 비-CNS 최대
    top_cns: tuple                # (조직, nTPM) — CNS 최대 (비교용)
    haematopoietic: tuple = ()    # HPA 가 집어 준 혈액·골수 세포형 (세포형, 값, 단위)
    source: str = "HPA26"


@dataclass(frozen=True)
class Liability:
    """말초 위험의 **1차 문헌 근거**. 🔴 수치를 적지 않는다 — 방향만.

    ⚠️ `kind` 를 분리하는 이유: *"뇌 전달 목적으로 이 수용체를 때렸을 때 말초 부작용이
    관찰됐다"*(`kind="부작용"`)와 *"이 수용체가 넓게 발현한다"*(`kind="발현"`)는 **전혀
    다른 강도의 근거**다. 한 칸에 섞으면 후자가 전자로 읽힌다 — 라이브러리에서 전자가
    있는 수용체는 **TFRC 하나뿐**이고, 그 사실이 흐려지면 이 절의 요점이 사라진다.
    """

    kind: str           # "부작용" = 전달 맥락의 실제 AE 보고 / "발현" = 넓은 발현 근거만
    direction: str      # 방향 한 줄. "n배"·"n%" 금지
    verified: bool      # 이 세션에서 출처 본문(초록)을 실제로 읽었나
    refs: tuple         # SOURCES 키들
    unverified_note: str = ""   # 읽지 못한 것이 있으면 무엇인지


@dataclass(frozen=True)
class Receptor:
    """한 수용체의 조직 발현 프로파일과 **도출된** off-target 단계."""

    gene: str | None            # HGNC 기호. None = 단백질 수용체가 아니거나 미지정
    label: str                  # 사람이 읽는 이름
    basis: str                  # "hpa" | "mechanism" | "unresolvable"
    tier: str                   # 선언값 — derive_tier 와 어긋나면 임포트 실패
    hpa: HpaRecord | None = None
    liability: Liability | None = None
    why: str = ""               # 이 판정을 그렇게 본 근거 한 줄
    candidates: tuple = ()      # basis="unresolvable" 일 때 후보 유전자들


# ---------------------------------------------------------------------------
# ② 불리언 — ①에서 **계산된다**. 손으로 선언하지 않는다.
# ---------------------------------------------------------------------------
CRITERIA_LABELS = (
    ("detected_all", "전조직검출"),
    ("broad_rna", f"광범위RNA(≥{BROAD_MIN_TISSUES}/{N_NONCNS_TISSUES}조직)"),
    ("peripheral_hotspot", f"말초핫스팟(≥{HOTSPOT_NTPM:.0f}nTPM)"),
    ("haematopoietic", "혈액·골수세포"),
    ("protein_outside_cns", "단백질검출(IH)"),
)


def criteria(h: HpaRecord) -> dict:
    """HPA 기록 → 5개 감사 불리언. **이 함수가 불리언의 정의다.**"""
    return {
        "detected_all": h.rna_distribution == "Detected in all",
        "broad_rna": h.n_noncns_ge10 >= BROAD_MIN_TISSUES,
        "peripheral_hotspot": h.n_noncns_ge30 >= 1,
        "haematopoietic": bool(h.haematopoietic),
        "protein_outside_cns": h.protein_distribution != "Not detected",
    }


def derive_tier(r: Receptor) -> str:
    """수용체 → 거친 3단계. **이 함수가 단계의 정의다**(위 표는 그 결과일 뿐).

    규칙 (2개 불리언만 쓴다):
      · `basis="mechanism"`(수용체가 없는 CPP) → **광범위**.
        수용체가 없다는 것은 조직 제한이 **없다**는 뜻이므로 가장 넓은 칸이다.
      · `basis="unresolvable"`(수용체 라벨을 유전자로 환원할 수 없음) → **미결정**.
        🔴 「좁음」으로 내려보내지 않는다 — 조회 불가는 안전이 아니다(래칫 R6:
        "측정 불가는 `None` 이다").
      · `broad_rna` ∧ `peripheral_hotspot` → **광범위**
      · `peripheral_hotspot` 만 → **편중** (말초 한두 장기에 몰려 있다)
      · 둘 다 아니면 → **좁음**

    ⚠️ 나머지 세 불리언(`detected_all`·`haematopoietic`·`protein_outside_cns`)은
       **단계에 넣지 않는다.** 보고만 한다. 이유: `haematopoietic` 과 `liability`
       의 **거짓은 "없다"가 아니라 "우리가 못 찾았다"** 이고, 그것을 단계에 넣으면
       조사 부족이 안전으로 인코딩된다. 같은 함정의 반대쪽 얼굴이 §19-4 의
       폐지된 BBB 하한이다.
    """
    if r.basis == "mechanism":
        return TIER_BROAD
    if r.basis == "unresolvable" or r.hpa is None:
        return TIER_UNKNOWN
    c = criteria(r.hpa)
    if c["broad_rna"] and c["peripheral_hotspot"]:
        return TIER_BROAD
    if c["peripheral_hotspot"]:
        return TIER_FOCAL
    return TIER_NARROW


def tier_rank(tier: str) -> int:
    """넓은 쪽이 작다. 미결정은 맨 뒤(= 비교에서 제외하기 위한 자리)."""
    return TIER_ORDER.index(tier) if tier in TIER_ORDER else len(TIER_ORDER)


# ---------------------------------------------------------------------------
# ① 외부 데이터 — 전부 2026-11-12 조회값
# ---------------------------------------------------------------------------
RECEPTORS: dict[str, Receptor] = {
    "TFRC": Receptor(
        gene="TFRC", label="TfR1 (transferrin receptor 1, CD71)", basis="hpa",
        tier=TIER_BROAD,
        hpa=HpaRecord(
            ensembl="ENSG00000072274", uniprot="P02786",
            rna_specificity="Low tissue specificity", rna_distribution="Detected in all",
            protein_specificity="Tissue enriched", protein_distribution="Detected in all",
            reliability_ih="Enhanced",
            n_noncns_ge10=35, n_noncns_ge30=13,
            top_noncns=("bone marrow", 98.6), top_cns=("retina", 27.3),
            haematopoietic=(("Erythrocyte progenitors", 976.6, "단일세포 nCPM"),
                            ("Syncytiotrophoblasts", 728.5, "단일세포 nCPM"),
                            ("non-classical monocyte", 9.2, "혈액세포 nTPM")),
        ),
        liability=Liability(
            kind="부작용",
            direction="뇌 전달 목적의 anti-TfR1 항체 투여 뒤 **급성 임상 증상과 순환 "
                      "reticulocyte 감소**가 관찰되며, Fc 효과기 기능 제거가 증상을 "
                      "완화하고 감소를 **부분적으로만** 되살린다(잔여분은 보체 매개). "
                      "TfR1 친화도를 **낮추면** 뇌 흡수와 전신 노출·안전성이 함께 개선된다. "
                      "후속 연구는 친화도·발렌시·에피토프를 바꿔 뇌 흡수와 reticulocyte "
                      "감소를 함께 측정한다.",
            verified=True, refs=("PM23636093", "PM41287279", "PM24411731"),
            unverified_note="🔴 초록만 읽었다. 감소 폭·용량·시점 등 **정량값은 본문 "
                            "표/그림 대조가 필요해 하나도 싣지 않았다**(방향만).",
        ),
        why="HPA: 비-CNS 37조직 중 **35조직이 nTPM ≥ 10**(≥30 은 13조직), 최대는 "
            "**골수 98.6 nTPM** 으로 CNS 최대(retina 27.3)보다 높다. 단일세포에서 "
            "**적혈구 전구세포 976.6 nCPM**. 단백질 수준도 'Detected in all'(IH 신뢰도 "
            "Enhanced). → 광범위 + 조혈 + 문헌 말초 부작용까지 **세 칸 전부 켜진다**.",
    ),
    "LRP1": Receptor(
        gene="LRP1", label="LRP1 (LDL receptor-related protein 1)", basis="hpa",
        tier=TIER_BROAD,
        hpa=HpaRecord(
            ensembl="ENSG00000123384", uniprot="Q07954",
            rna_specificity="Low tissue specificity", rna_distribution="Detected in all",
            protein_specificity="Low tissue specificity", protein_distribution="Detected in all",
            reliability_ih="Supported",
            n_noncns_ge10=28, n_noncns_ge30=8,
            top_noncns=("adipose tissue", 102.4), top_cns=("cerebellum", 34.1),
            haematopoietic=(("classical monocyte", 4.3, "혈액세포 nTPM"),
                            ("myeloid DC", 1.6, "혈액세포 nTPM")),
        ),
        liability=Liability(
            kind="발현",
            direction="LRP1 은 **여러 조직에 널리 발현**하는 대형 엔도사이토시스 수용체이고 "
                      "30종 이상의 리간드를 인식한다. 조직 특이 KO 는 혈관·CNS·**대식세포**·"
                      "**지방세포**에서의 기여를 보인다. ⚠️ 이것은 '수용체가 넓다'는 근거이고 "
                      "**Angiopep-2 셔틀의 말초 부작용 보고가 아니다** — 그런 보고는 찾지 못했다.",
            verified=True, refs=("PM18626063",),
            unverified_note="🔴 Angiopep-2/ANG1005 의 **말초 off-target 부작용 1차 보고를 "
                            "확보하지 못했다.** 이 칸의 거짓(부작용 미확인)은 **안전이 아니라 "
                            "미조사**다.",
        ),
        why="HPA: 28/37 조직이 nTPM ≥ 10, 최대 **지방조직 102.4 nTPM** 이 CNS 최대"
            "(cerebellum 34.1)의 3배 수준. 단백질도 'Low tissue specificity · Detected in all'. "
            "→ 광범위.",
    ),
    "LDLR": Receptor(
        gene="LDLR", label="LDLR (low density lipoprotein receptor)", basis="hpa",
        tier=TIER_BROAD,
        hpa=HpaRecord(
            ensembl="ENSG00000130164", uniprot="P01130",
            rna_specificity="Tissue enhanced", rna_distribution="Detected in all",
            protein_specificity="Group enriched", protein_distribution="Detected in some",
            reliability_ih="Approved",
            n_noncns_ge10=25, n_noncns_ge30=8,
            top_noncns=("adrenal gland", 112.1), top_cns=("choroid plexus", 32.8),
            haematopoietic=(("basophil", 46.2, "혈액세포 nTPM"),
                            ("granulocytes", 46.2, "혈액계통 nTPM")),
        ),
        liability=None,
        why="HPA: 25/37 조직이 nTPM ≥ 10, 최대 **부신 112.1 nTPM**(간 56.4 · 폐 51.7). "
            "RNA 는 'Tissue enhanced' 지만 **그 enhanced 조직이 부신**이라 CNS 선택성이 "
            "아니다. → 광범위. ⚠️ 말초 부작용 1차 근거는 **미조사**.",
    ),
    "LEPR": Receptor(
        gene="LEPR", label="LepR (leptin receptor, ObR)", basis="hpa",
        tier=TIER_FOCAL,
        hpa=HpaRecord(
            ensembl="ENSG00000116678", uniprot="P48357",
            rna_specificity="Tissue enriched", rna_distribution="Detected in all",
            protein_specificity="Tissue enriched", protein_distribution="Detected in single",
            reliability_ih="Approved",
            n_noncns_ge10=12, n_noncns_ge30=2,
            top_noncns=("liver", 121.9), top_cns=("choroid plexus", 16.2),
            haematopoietic=(("Bone marrow MSCs", 2714.9, "단일세포 nCPM"),
                            ("basophil", 5.0, "혈액세포 nTPM")),
        ),
        liability=None,
        why="HPA: 12/37 조직만 nTPM ≥ 10 이라 **TFRC·LRP1·LDLR 보다 좁다**. 그러나 "
            "**간 121.9 nTPM** 은 이 표에서 두 번째로 큰 말초 단일 핫스팟이고, 단백질 "
            "수준은 'Tissue enriched / Detected in **single**' = **간 하나**다. "
            "단일세포에서 **골수 MSC 2714.9 nCPM**(이 표 최대값). → 편중.",
    ),
    "CHRNA7": Receptor(
        gene="CHRNA7", label="α7 nAChR (nicotinic acetylcholine receptor α7)", basis="hpa",
        tier=TIER_NARROW,
        hpa=HpaRecord(
            ensembl="ENSG00000175344", uniprot="P36544",
            rna_specificity="Low tissue specificity", rna_distribution="Detected in many",
            protein_specificity="Not detected", protein_distribution="Not detected",
            reliability_ih="Enhanced",
            n_noncns_ge10=0, n_noncns_ge30=0,
            top_noncns=("small intestine", 9.1), top_cns=("choroid plexus", 3.7),
            haematopoietic=(),
        ),
        liability=Liability(
            kind="발현",
            direction="α7nAChR 은 **백혈구에 발현해 콜린성 항염 경로**를 매개한다고 "
                      "보고된다(대식세포 NF-κB 억제·호중구 이동 억제). ⚠️ 이는 HPA 의 "
                      "혈액세포 RNA 'Not detected' 와 **방향이 어긋난다** — 출처 간 "
                      "불일치이고 한쪽을 고르지 않는다(`docs/protocol.md` 출처 표기 규칙).",
            verified=False, refs=(),
            unverified_note="🔴 **1차 출처를 확보하지 못했다.** 이 세션에서 읽은 것은 "
                            "bioRxiv 프리프린트의 배경 서술이라 `verified=False` 로 둔다. "
                            "경로의 1차 문헌(Wang 2003 Nature 계열)을 열지 못했다.",
        ),
        why="HPA: 비-CNS 37조직에서 **nTPM ≥ 10 인 조직이 0개**(최대 소장 9.1), CNS 최대도 "
            "3.7 로 낮다. 단백질 수준은 **'Not detected'**(IH 신뢰도 Enhanced). "
            "→ 좁음. 🔴 **단 이것은 '뇌에 있다'는 근거가 아니다** — 같은 표가 **CNS 에서도** "
            "낮다고 말한다. 좁음은 '표적이 확인됐다'가 아니라 **'말초 발현 신호가 작다'** 뿐이다.",
    ),
    #: RVG29 의 라벨이 "nAChR(아형 미지정)" 이라 **유전자로 환원되지 않는다.**
    #: 아형을 고르는 순간 단계가 바뀐다 — 그래서 단계를 주지 않고 후보를 병기한다.
    "NACHR_UNSPEC": Receptor(
        gene=None, label="nAChR — **아형 미지정**(라이브러리 라벨 그대로)",
        basis="unresolvable", tier=TIER_UNKNOWN,
        candidates=(
            ("CHRNA7", "ENSG00000175344", "비-CNS ≥10nTPM **0/37** · 단백질 Not detected → 좁음"),
            ("CHRNA4", "ENSG00000101204", "비-CNS ≥10nTPM **2/37** 이지만 **간 82.6 nTPM** "
                                          "(단일세포 Hepatocytes 70.5 nCPM) → **편중**"),
            ("CHRNB2", "ENSG00000160716", "비-CNS ≥10nTPM **0/37**(최대 뇌하수체 3.4) → 좁음"),
        ),
        why="🔴 **아형에 따라 단계가 갈린다**: CHRNA7·CHRNB2 는 좁음인데 **CHRNA4 는 간 "
            "82.6 nTPM 으로 편중**이다. `config.SHUTTLES` 의 RVG29 라벨이 "
            "'nAChR(니코틴성 수용체)' 로만 적혀 있어 **어느 쪽인지 고를 근거가 우리에게 없다**. "
            "문헌도 RVG29 의 수용체를 nAChR 로 쓰면서 아형을 확정하지 않는다. "
            "→ 미결정. (후보 전부 조회했고 수치는 위에 적었다.)",
    ),
    #: Tet1 의 표적은 **단백질이 아니다** — HPA 에 유전자 페이지가 없다.
    "GT1B": Receptor(
        gene=None, label="GT1b 갱글리오사이드 (당지질 — **단백질 수용체가 아니다**)",
        basis="unresolvable", tier=TIER_UNKNOWN,
        candidates=(
            ("ST3GAL2", "ENSG00000157350", "GT1b **생합성 효소**. 비-CNS ≥10nTPM **1/37** · "
                                           "단백질 Tissue enriched(lymphoid tissue) — "
                                           "🔴 **효소 발현은 당지질 분포가 아니다**"),
            ("B3GALT4", "(HPA 조회 실패 — API 가 이 Ensembl ID 로 0행 반환)",
             "🔴 **조회 못 했다.** 대체하지 않고 공백으로 남긴다"),
        ),
        why="🔴 **발현 데이터로 답할 수 없는 표적이다.** GT1b 는 당지질이라 HPA/GTEx 의 "
            "유전자 발현 축에 **존재하지 않는다**. 생합성 효소(ST3GAL2 등)를 대리로 쓰면 "
            "'효소 RNA = 당지질 표면 밀도' 라는 **검증되지 않은 외삽**을 하게 된다 — "
            "이 프로젝트가 세 번 데인 '측정 불가를 측정값으로 인코딩' 과 같은 함정이다. "
            "→ 미결정. (별도로 `core/evidence.py` 가 Tet1 을 **근거 등급 D** 로 두고 "
            "Tet1-LNP 가 뇌가 아니라 **말초 신경**으로 간다는 보고를 적어 둔다.)",
    ),
    #: CPP 는 수용체가 없다 = 조직 제한이 없다. 이 칸만 `basis="mechanism"` 이다.
    "NO_RECEPTOR": Receptor(
        gene=None, label="수용체 없음 — 막 직접투과(CPP)", basis="mechanism",
        tier=TIER_BROAD,
        why="수용체를 쓰지 않으므로 **조직 제한의 출처 자체가 없다.** 발현 데이터로 좁힐 "
            "대상이 없다는 뜻이고, 그래서 가장 넓은 칸이다. (이 칸은 HPA 를 쓰지 않는다 — "
            "`basis='mechanism'`.)",
    ),
}


# ---------------------------------------------------------------------------
# 셔틀 → 수용체 (라이브러리 라벨에서 옮긴 것이고, **우리 측정이 아니다**)
# ---------------------------------------------------------------------------
#: 🔴 이 매핑의 출처는 `core/config.py SHUTTLES[*]["target"]` 의 **문헌 라벨**이다.
#:    우리가 결합을 측정한 것이 아니다. 라벨이 틀리면 §50 전체가 틀린다.
SHUTTLE_RECEPTORS: dict[str, tuple] = {
    "Angiopep-2":     ("LRP1",),
    "ApoE(159-167)2": ("LDLR", "LRP1"),       # 라벨이 "LDLR/LRP1"
    "Leptin30":       ("LEPR",),
    "T7":             ("TFRC",),
    "B6":             ("TFRC",),
    "RVG29":          ("NACHR_UNSPEC",),
    "CDX":            ("CHRNA7",),
    "THR":            ("TFRC",),
    "Tet1":           ("GT1B",),
    "TAT":            ("NO_RECEPTOR",),
    "Penetratin":     ("NO_RECEPTOR",),
    "SynB1":          ("NO_RECEPTOR",),
}


# ---------------------------------------------------------------------------
# 드리프트 방지 — `core/evidence.py` 와 같은 임포트 시 검증
# ---------------------------------------------------------------------------
_miss = set(SHUTTLES) - set(SHUTTLE_RECEPTORS)
_extra = set(SHUTTLE_RECEPTORS) - set(SHUTTLES)
if _miss or _extra:  # pragma: no cover
    raise RuntimeError(
        f"SHUTTLE_RECEPTORS 가 config.SHUTTLES 와 어긋난다 — 누락 {sorted(_miss)} / "
        f"잉여 {sorted(_extra)}")

_unknown = {s: [k for k in ks if k not in RECEPTORS] for s, ks in SHUTTLE_RECEPTORS.items()}
_unknown = {s: v for s, v in _unknown.items() if v}
if _unknown:  # pragma: no cover
    raise RuntimeError(f"RECEPTORS 에 없는 수용체 키를 참조한다: {_unknown}")

_bad_tier = {k: r.tier for k, r in RECEPTORS.items()
             if r.tier not in TIER_ORDER and r.tier != TIER_UNKNOWN}
if _bad_tier:  # pragma: no cover
    raise RuntimeError(f"알 수 없는 단계: {_bad_tier}")

_bad_src = {k: r.hpa.source for k, r in RECEPTORS.items()
            if r.hpa is not None and r.hpa.source not in SOURCES}
if _bad_src:  # pragma: no cover
    raise RuntimeError(f"SOURCES 에 없는 출처 ID: {_bad_src}")

_bad_ref = {k: [x for x in r.liability.refs if x not in SOURCES]
            for k, r in RECEPTORS.items() if r.liability is not None}
_bad_ref = {k: v for k, v in _bad_ref.items() if v}
if _bad_ref:  # pragma: no cover
    raise RuntimeError(f"SOURCES 에 없는 인용 ID: {_bad_ref}")

# 🔴 **단계는 도출값이다**: 선언(`Receptor.tier`)과 규칙(`derive_tier`)이 어긋나면 임포트 실패.
#    단계를 바꾸고 싶으면 HPA 기록이나 규칙을 고칠 것이고, 글자만 손으로 바꾸지 말 것.
_drift = {k: (r.tier, derive_tier(r)) for k, r in RECEPTORS.items()
          if r.tier != derive_tier(r)}
if _drift:  # pragma: no cover
    raise RuntimeError(
        "단계가 규칙에서 도출되지 않는다(선언, 도출): "
        + ", ".join(f"{k}{v}" for k, v in sorted(_drift.items())))
del _miss, _extra, _unknown, _bad_tier, _bad_src, _bad_ref, _drift


# ---------------------------------------------------------------------------
# 대조
# ---------------------------------------------------------------------------
#: 서열 축의 판정 경계. `core/selectivity.py` 의 `level` 과 **같은 임계**를 쓴다
#: (높음 ≥ 0.6 · 보통 ≥ 0.33 · 그 외 낮음). 여기서 다시 정의하지 않고 그 함수의
#: `risk_level` 문자열을 그대로 읽는다 — 임계가 바뀌면 같이 바뀌어야 하기 때문이다.
AXIS_SELECTIVE = "낮음"
AXIS_RISKY = "높음"


@dataclass
class Row:
    shuttle: str
    seq: str
    off: float
    axis_level: str
    receptor_keys: tuple
    receptor_label: str
    tier: str
    flag: str
    note: str = ""
    criteria: dict = field(default_factory=dict)


def shuttle_tier(name: str) -> tuple:
    """셔틀의 수용체 단계 — 여러 수용체를 노리면 **가장 넓은 쪽**을 쓴다(보수적).

    미결정만 있으면 미결정이다. 미결정과 실제 단계가 섞이면 실제 단계를 쓰되
    그 사실을 `note` 에 남긴다.
    """
    keys = SHUTTLE_RECEPTORS[name]
    tiers = [derive_tier(RECEPTORS[k]) for k in keys]
    known = [t for t in tiers if t in TIER_ORDER]
    if not known:
        return TIER_UNKNOWN, keys
    return min(known, key=tier_rank), keys


def compare() -> list:
    """12종 전수 대조 — 서열 축(실제 함수 호출) vs 수용체 단계."""
    rows = []
    for name, v in SHUTTLES.items():
        seq = v["seq"]
        r = assess_selectivity(seq)
        tier, keys = shuttle_tier(name)
        labels = " + ".join(RECEPTORS[k].label for k in keys)

        flag, note = "—", ""
        if tier == TIER_UNKNOWN:
            flag = "⚠️ 미결정"
            note = "수용체 라벨을 발현 데이터로 환원할 수 없다"
        elif r.risk_level == AXIS_SELECTIVE and tier == TIER_BROAD:
            flag = "🔴 축 낙관"
            note = (f"축은 off {r.off_target_risk:.2f}(낮음)인데 수용체는 광범위")
        elif r.risk_level == AXIS_RISKY and tier == TIER_NARROW:
            flag = "🔴 축 비관"
            note = f"축은 off {r.off_target_risk:.2f}(높음)인데 수용체는 좁음"
        elif r.risk_level == AXIS_RISKY and tier == TIER_BROAD:
            flag = "🟢 일치"
            note = "축도 수용체도 넓다"
        elif r.risk_level == AXIS_SELECTIVE and tier == TIER_NARROW:
            flag = "🟢 일치"
            note = "축도 수용체도 좁다"
        else:
            flag = "◻ 중간"
            note = f"축 {r.risk_level} · 수용체 {tier}"

        crit = {}
        for k in keys:
            rec = RECEPTORS[k]
            if rec.hpa is not None:
                crit[k] = criteria(rec.hpa)

        rows.append(Row(shuttle=name, seq=seq, off=r.off_target_risk,
                        axis_level=r.risk_level, receptor_keys=keys,
                        receptor_label=labels, tier=tier, flag=flag, note=note,
                        criteria=crit))
    return rows


# ---------------------------------------------------------------------------
# 보고
# ---------------------------------------------------------------------------
def receptor_report() -> str:
    """수용체별 외부 데이터 표 — 수치와 **측정 조건을 같은 줄에** 쓴다."""
    out = [
        f"| 수용체 | 유전자 / Ensembl | HPA RNA 조직특이성 | 비-CNS {N_NONCNS_TISSUES}조직 중 "
        f"nTPM≥{PRESENT_NTPM:.0f} | ≥{HOTSPOT_NTPM:.0f} | 비-CNS 최대 (조직, nTPM) | "
        "CNS 최대 (조직, nTPM) | 단백질(IH) | → 단계 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for k in sorted(RECEPTORS, key=lambda x: (tier_rank(derive_tier(RECEPTORS[x])), x)):
        r = RECEPTORS[k]
        if r.hpa is None:
            out.append(f"| {r.label} | — | — | — | — | — | — | — | **{derive_tier(r)}** |")
            continue
        h = r.hpa
        out.append(
            f"| {r.label} | {r.gene} / {h.ensembl} | {h.rna_specificity} | "
            f"**{h.n_noncns_ge10}**/{N_NONCNS_TISSUES} | {h.n_noncns_ge30} | "
            f"{h.top_noncns[0]} {h.top_noncns[1]} | {h.top_cns[0]} {h.top_cns[1]} | "
            f"{h.protein_distribution} ({h.protein_specificity}, 신뢰도 {h.reliability_ih}) | "
            f"**{derive_tier(r)}** |")
    out += [
        "",
        f"측정 조건: HPA v26.0 consensus tissue nTPM · 조회 2026-11-12 · 분모 = "
        f"consensus 47조직 − CNS 10조직 = **{N_NONCNS_TISSUES}조직** · "
        f"단계 규칙 = `broad_rna(≥{BROAD_MIN_TISSUES}조직) ∧ hotspot(≥{HOTSPOT_NTPM:.0f}nTPM)` "
        f"→ 광범위 / `hotspot` 만 → 편중 / 그 외 → 좁음.",
    ]
    return "\n".join(out)


def rubric_report() -> str:
    """감사용 — 단계가 **어떤 불리언에서 도출됐는지**."""
    out = ["| 수용체 | " + " | ".join(lbl for _, lbl in CRITERIA_LABELS)
           + " | → 단계 | 문헌 근거(종류) |",
           "|---|" + "---|" * (len(CRITERIA_LABELS) + 2)]
    for k in sorted(RECEPTORS, key=lambda x: (tier_rank(derive_tier(RECEPTORS[x])), x)):
        r = RECEPTORS[k]
        if r.hpa is None:
            cells = ["—"] * len(CRITERIA_LABELS)
        else:
            c = criteria(r.hpa)
            cells = ["✅" if c[key] else "—" for key, _ in CRITERIA_LABELS]
        if r.liability is None:
            lia = "🔴 미조사"
        elif r.liability.verified:
            lia = f"✅ {r.liability.kind} · 방향 확인(초록)"
        else:
            lia = f"⚠️ {r.liability.kind} · 1차 출처 미확보"
        out.append(f"| {r.label} | " + " | ".join(cells)
                   + f" | **{derive_tier(r)}** | {lia} |")
    out += [
        "",
        "🔴 **마지막 두 칸(`혈액·골수세포`·`문헌 근거`)의 거짓은 '없다'가 아니라 "
        "'우리가 못 찾았다'다.** 그래서 **단계 계산에 들어가지 않는다** — 조사 부족이 "
        "안전으로 인코딩되면 폐지된 `MIN_BBB_SIGNAL`(§19-4)과 같은 사고가 된다.",
        "⚠️ 종류 구분: **부작용** = 뇌 전달 목적으로 그 수용체를 때렸을 때 말초 AE 가 "
        "관찰됐다는 보고 / **발현** = 넓게 발현한다는 근거만. 라이브러리에서 **부작용** "
        "근거를 가진 수용체는 **TFRC 하나뿐**이다.",
    ]
    return "\n".join(out)


def comparison_report(rows: list) -> str:
    out = ["| 셔틀 | 서열 축 off | 축 판정 | 표적 수용체(문헌 라벨) | 수용체 단계 | 불일치 |",
           "|---|---|---|---|---|---|"]
    for r in rows:
        out.append(f"| {r.shuttle} | {r.off:.2f} | {r.axis_level} | {r.receptor_label} | "
                   f"**{r.tier}** | {r.flag} |")
    return "\n".join(out)


def headline(rows: list) -> str:
    opt = [r for r in rows if r.flag == "🔴 축 낙관"]
    pes = [r for r in rows if r.flag == "🔴 축 비관"]
    unk = [r for r in rows if r.flag == "⚠️ 미결정"]
    agr = [r for r in rows if r.flag == "🟢 일치"]
    mid = [r for r in rows if r.flag == "◻ 중간"]

    # 같은 수용체를 노리는 셔틀들 사이의 축 분산 — **배수로 쓰지 않는다**(§5 보고 규약)
    groups = {}
    for r in rows:
        if r.tier == TIER_UNKNOWN:
            continue
        groups.setdefault(r.receptor_keys, []).append(r)
    spread = []
    for keys, g in sorted(groups.items()):
        if len(g) < 2:
            continue
        lo, hi = min(x.off for x in g), max(x.off for x in g)
        spread.append(f"{RECEPTORS[keys[0]].label}: {len(g)}종 — "
                      f"off {lo:.2f} ~ {hi:.2f} ("
                      + ", ".join(f"{x.shuttle} {x.off:.2f}" for x in sorted(g, key=lambda y: y.off))
                      + ")")

    lines = [
        f"셔틀 **{len(rows)}종** 중 서열 축이 '낮음'으로 부르는데 표적 수용체는 "
        f"**광범위**(비-CNS {N_NONCNS_TISSUES}조직 중 ≥{BROAD_MIN_TISSUES}조직에서 "
        f"nTPM ≥ {PRESENT_NTPM:.0f})인 것이 **{len(opt)}종**: "
        + (", ".join(r.shuttle for r in opt) or "없음") + ".",
        f"반대 방향(축 '높음' · 수용체 좁음) 불일치는 **{len(pes)}종**"
        + (f"({', '.join(r.shuttle for r in pes)})" if pes else "") + ".",
        f"방향이 일치하는 것 **{len(agr)}종**"
        + (f"({', '.join(r.shuttle for r in agr)})" if agr else "")
        + f" · 중간 **{len(mid)}종**"
        + (f"({', '.join(r.shuttle for r in mid)})" if mid else "")
        + f" · 수용체 라벨을 발현 데이터로 환원할 수 없는 것 **{len(unk)}종**"
        + (f"({', '.join(r.shuttle for r in unk)})" if unk else "") + ".",
        "",
        "같은 수용체를 공유하는 셔틀 사이에서 서열 축이 벌리는 폭 "
        "(🔴 **배수로 쓰지 말 것** — 교정되지 않은 축이다):",
    ]
    lines += ["  · " + s for s in spread]
    return "\n".join(lines)


def axis_decomposition() -> str:
    """서열 축이 **실제 라이브러리에서** 무엇을 재고 있는지 분해한다.

    축의 식은 `raw = 1 − (1 − 전하밀도항)(1 − 친유성항)` 이라 "전하 **또는** 친유성"
    으로 읽히지만, 친유성항은 `GRAVY > 0` 일 때만 켜진다. 라이브러리 12종에서 실제로
    몇 개가 켜지는지 세어 본다 — 이것은 수용체 데이터와 무관한 **코드 내부 사실**이라
    여기서 재확인만 한다(계산은 전부 `core/selectivity.py` 의 공개 함수에서 온다).
    """
    from core.selectivity import _clamp  # noqa: PLC0415 — 축 내부 상수를 재구현하지 않기 위해
    rows, off_lip = [], 0
    for name, v in SHUTTLES.items():
        r = assess_selectivity(v["seq"])
        cat = _clamp(r.charge_density / 0.5) if r.charge_density > 0 else 0.0
        lip = _clamp(r.gravy / 0.8) if r.gravy > 0 else 0.0
        charge_only = round(_clamp(cat * (1 - 0.6 * r.rmt_similarity)), 2)
        same = charge_only == r.off_target_risk
        off_lip += int(lip == 0.0)
        rows.append((name, r.off_target_risk, r.charge_density, r.gravy, cat, lip,
                     charge_only, same))
    out = ["| 셔틀 | off | 전하밀도 | GRAVY | 전하항 | 친유성항 | 전하항만으로 계산한 off | 동일? |",
           "|---|---|---|---|---|---|---|---|"]
    for n, off, dens, gravy, cat, lip, co, same in rows:
        out.append(f"| {n} | {off:.2f} | {dens:+.2f} | {gravy:+.2f} | {cat:.3f} | "
                   f"{lip:.3f} | {co:.2f} | {'✅' if same else '🔴 다름'} |")
    n_same = sum(1 for *_, s in rows if s)
    out += [
        "",
        f"🔴 **친유성항이 0 인 셔틀 = {off_lip}/{len(rows)}** (GRAVY ≤ 0 이면 항이 꺼진다). "
        f"그래서 off 값이 전하항만으로 계산한 값과 **정확히 같은 셔틀이 {n_same}/{len(rows)}** 이고, "
        "다른 것은 ApoE(159-167)2 하나다(GRAVY +0.41). "
        "즉 *\"전하 **또는** 친유성\"* 이라는 서술은 식에서는 맞지만 **이 라이브러리에서는 "
        "사실상 전하 하나**다. 측정 조건: 라이브러리 12종 셔틀 단독 서열 · "
        "`core/selectivity.assess_selectivity` 직접 호출 · Biopython GRAVY 정상 동작 확인.",
    ]
    return "\n".join(out)


def limitations() -> str:
    return "\n".join([
        "이 대조가 **보지 못하는 것**:",
        "ⓐ **RNA ≠ 단백질 ≠ 혈액에서 접근 가능**. HPA consensus 는 벌크 RNA 이고, "
        "단백질 칸은 항체 IH 주석이다. 혈관 내강 쪽에서 혈액 중 분자가 실제로 닿는 "
        "수용체 밀도는 **어느 칸에도 없다**.",
        "ⓑ **뇌 혈관내피는 어떤 벌크 조직 시료에서도 극소 분획**이다. 'cerebral cortex "
        "nTPM' 은 거의 전부 뉴런·글리아 신호이므로, 이 표의 CNS 열은 **BBB 발현의 대리가 "
        "아니다**. 즉 이 대조는 분자가 **말초 어디에 끌려갈 수 있나**만 말하고 "
        "**뇌에 얼마나 가나**는 말하지 않는다.",
        "ⓒ **셔틀→수용체 매핑은 우리가 물려받은 문헌 라벨**이고 측정이 아니다"
        "(`core/config.py SHUTTLES[*]['target']`). 라벨이 틀리면 이 절 전체가 틀린다. "
        "실제로 12종 중 **2종**(RVG29 의 아형 미지정 · Tet1 의 비단백질 표적)은 "
        "라벨만으로 조회가 불가능했다.",
        "ⓓ **발현 ≠ 흡수**. 수용체가 있어도 트랜스사이토시스가 일어나지 않을 수 있고, "
        "없어도 비특이 흡수는 일어난다 — 이 표는 **기회의 지도**이고 유량(flux)이 아니다.",
        "ⓔ **단계 임계 2개(≥20조직 · ≥30nTPM)에 외부 근거가 없다.** 우리가 고른 값이고, "
        "바꾸면 LDLR(25조직)·LEPR(12조직)의 칸이 먼저 움직인다.",
        "ⓕ **정량 대조 0편.** 모든 문헌 근거는 초록까지만 읽었고, 그래서 배수·백분율을 "
        "하나도 싣지 않았다(§40 의 교훈).",
        "ⓖ **종 차이 미반영.** 라이브러리 셔틀의 전달 근거는 대부분 마우스·랫이고 "
        "이 표는 **사람** 발현이다.",
    ])


def main() -> int:
    rows = compare()
    print("=" * 78)
    print("§50 수용체 수준 off-target 외부 대조 — 서열 축 vs HPA 발현 데이터")
    print("=" * 78)
    print()
    print("[1] 수용체별 외부 발현 데이터")
    print(receptor_report())
    print()
    print("[2] 단계 도출 감사표")
    print(rubric_report())
    print()
    print("[3] 12종 전수 대조")
    print(comparison_report(rows))
    print()
    print("[4] 서열 축 분해 — 이 라이브러리에서 축이 실제로 재는 것")
    print(axis_decomposition())
    print()
    print("[5] 한 줄 결론")
    print(headline(rows))
    print()
    print("[6] " + limitations())
    print()
    print("[7] 출처")
    for k, v in SOURCES.items():
        print(f"  · {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
