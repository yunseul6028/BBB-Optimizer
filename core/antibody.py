"""항체/나노바디 셔틀 파이프라인 (Tier-3 모달리티) — 골격.

펩타이드 셔틀과 달리 항체 셔틀(예: 항-TfR Fab/scFv/VHH)은 **3D CDR 구조 + 친화도(KD)**로
수용체에 결합하며, BBB 통과 여부는 **서열 조성으로 예측 불가**하다(그래서 deepB3P를 못 쓴다).
본 모듈은 그 평가 파이프라인의 **골격**을 제공한다 — 서열 DB(실제 항-TfR 항체·나노바디 서열)는
**추후 주입**하고, 지금은 표현·계산 로직·플러그 지점만 갖춘다.

파이프라인 단계:
  1) 표현       : AntibodyShuttle (포맷/표적/사슬/CDR/친화도/결합가/출처)
  2) 서열 기반 축(가능): 개발성(항체 liability)·용해도·전하 — 기존 엔진 재사용
                        (⚠️ 엔진은 펩타이드 보정이라 항체엔 근사; 항체 특화 TAP는 향후)
  3) 친화도 축(Tier-3 핵심): TfR **친화도 sweet-spot** 모델 (KD 입력 시) — 너무 세면 뇌 쪽에서
                        방출 안 돼 갇히고, 너무 약하면 결합 실패. 문헌상 중간~약 친화가 뇌흡수에 유리.
  4) 구조       : IgFold/ABodyBuilder(단일 항체) 또는 Fab–TfR **co-folding**(AF-Multimer/Boltz) 훅
                        — 외부·미장착이면 우아하게 skip.
  5) 결합가     : monovalent(1) vs bivalent(2). TfR 브레인셔틀은 **monovalent가 수송에 유리**한
                        것이 잘 알려짐(문헌) → bivalent는 감점.
  6) 집계       : 계산된 축 + Tier-3/DB로 보류된 축을 정직하게 구분해 반환.

⚠️ 모든 수치는 서열/휴리스틱 프록시이며 **실측(SPR/BLI·뇌흡수) 대체가 아니다.** 항체 서열·KD 출처는
   반드시 실재 확인(연구윤리) 후 주입한다.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

# 대표 RMT 항체 표적(참고) — 실제 브레인셔틀들이 노리는 수용체.
ANTIBODY_TARGETS = ["TfR", "IGF1R", "CD98hc"]

# 항체 도메인은 보통 ≥~110aa(VHH)~450aa(Fab), 펩타이드 화물/셔틀은 ≤~50aa. 길이가 주 신호.
ANTIBODY_MIN_LEN = 60
_AB_N_TERM = ("EVQL", "QVQL", "DVQL", "QVKL", "EVKL", "QITL", "AVQL",  # VH 관용 N말단
              "DIQM", "DIVM", "EIVL", "DIVL", "QSVL", "QAVL", "SYEL", "NFML")  # VL 관용


def detect_modality(seq: str) -> tuple[str, str]:
    """서열을 보고 'peptide'(펩타이드 화물) vs 'antibody'(항체/나노바디)로 자동 분기.

    길이가 주 신호(항체 도메인 ≥~110aa, 펩타이드 ≤~50aa), 항체 관용 N말단 모티프로 보강.
    반환 (modality, 근거 문자열).
    """
    s = "".join(ch for ch in (seq or "").upper() if ch.isalpha())
    n = len(s)
    if not s:
        return "peptide", "빈 입력"
    motif = s[:4] in _AB_N_TERM
    if n >= ANTIBODY_MIN_LEN or (n >= 45 and motif):
        return "antibody", f"길이 {n}aa{'·항체 N말단 모티프' if motif else ''} → 항체/나노바디로 판단"
    return "peptide", f"길이 {n}aa → 펩타이드 화물로 판단"


# ── CDR3(파라토프) 근사 추출 — dependency-free ──────────────────────────────
_H_FR4 = re.compile(r"WG.G")   # 중쇄 FR4 시작(WGQG/WGAG/WGKG…) → CDR-H3 뒤 앵커
_L_FR4 = re.compile(r"FG.G")   # 경쇄 FR4 시작(FGGG/FGQG…) → CDR-L3 뒤 앵커


def _cdr3(seq: str | None, fr4_re) -> str:
    """FR3 말단 보존 Cys ~ FR4 앵커 사이 = CDR3 루프(근사). 못 찾으면 ''."""
    s = "".join(ch for ch in (seq or "").upper() if ch.isalpha())
    if not s:
        return ""
    matches = list(fr4_re.finditer(s))
    if not matches:
        return ""
    # FR4(WG.G/FG.G)는 V도메인 C말단 근처 → **마지막 앵커**를 취해 CDR 내부의 spurious 모티프
    # (예: H3 안의 'WGLG', VHH FR3의 'WGGG')를 회피한다.
    m = matches[-1]
    cys = s.rfind("C", 0, m.start())
    if cys < 0 or (m.start() - cys) > 40:   # Cys~앵커가 비정상적으로 멀면 신뢰 불가
        return ""
    return s[cys + 1 : m.start()]


def extract_cdrs(vh: str | None, vl: str | None = None) -> dict:
    """항체 서열에서 결합을 주도하는 CDR3(H3·L3)를 dependency-free 근사 추출.

    CDR3(특히 H3)는 파라토프 결합의 지배적 결정인자다. 보존 Cys(FR3)~WG.G/FG.G(FR4) 앵커로 잡는다.
    ⚠️ Kabat 넘버링(ANARCI) 없이 모티프 근사 — H1/H2/L1/L2는 신뢰도가 낮아 잡지 않는다.
    반환 {h3, l3, paratope, method}. 못 잡으면 paratope=전체 Fv(fallback).
    """
    h3 = _cdr3(vh, _H_FR4)
    l3 = _cdr3(vl, _L_FR4) if vl else ""
    para = h3 + l3
    if para:
        method = "CDR3 근사(H3" + ("+L3" if l3 else "") + ")"
    else:
        para = "".join(ch for ch in ((vh or "") + (vl or "")).upper() if ch.isalpha())
        method = "CDR3 추출 실패 → 전체 Fv(fallback)"
    return {"h3": h3, "l3": l3, "paratope": para, "method": method}


# 항체/나노바디 셔틀 라이브러리. 스키마:
#   {name: {"target": str, "fmt": "VHH|scFv|Fab|IgG", "vh": str|None, "vl": str|None,
#           "kd_nM": float|None, "valency": 1|2, "source": str}}
# ⚠️ **서열(vh/vl)·KD는 지어내지 않는다.** 아래 프리셋은 문헌으로 공개된 '설계 사실'
#   (표적·포맷·결합가)만 담고, 서열/친화도는 None → 검증된 값을 확인 후 채운다.
ANTIBODY_SHUTTLES: dict[str, dict] = {
    # --- 서열 보유 (사용자 제공 · 표적/PDB 주석은 원문 대조 권장) ---
    "8D3 (항-마우스 TfR1 · scFv)": {
        "target": "Mouse TfR1", "fmt": "scFv", "valency": 1, "kd_nM": None,
        "vh": ("EVQLQQSGAELVKPGASVKLSCTASGFNIKDTYIHWVKQRPEQGLEWIGRIDPANGNTKY"
               "DPKFQGKATITADTSSNTAYLQLSSLTSEDTAVYYCARGLYGYFDVWGAGTTVTVSS"),
        "vl": ("DIQMTQSPASLSASVGETVTITCRASENIYSYLAWYQQKQGKSPQLLVYNAKTLAEGVPS"
               "RFSGSGSGTQFSLKINSLQPEDFGSYYCQHHYGTPFTFGGGTKLEIK"),
        "source": "사용자 제공 · PDB 6D04. 표적/구조 원문 대조 권장.",
    },
    "Ri18b (항-마우스 TfR1 · VHH)": {
        "target": "Mouse TfR1", "fmt": "VHH", "valency": 1, "kd_nM": None,
        "vh": ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAIRWSGGSTYY"
               "ADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGRGSYWYFDYWGQGTQVTVSS"),
        "vl": None,
        "source": "사용자 제공 · 저-nM 중간 친화(정성). 원문 대조 권장.",
    },
    "FC5 (BBB 투과 VHH)": {
        "target": "CEACAM/시알산 (원문 확인)", "fmt": "VHH", "valency": 1, "kd_nM": None,
        "vh": ("EVQLQASGGGLVQAGGSLRLSCAASGFKITHYTMGWFRQAPGKEREFVAGITWGGGSTYY"
               "ADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYVCAAGSTSTATMGGYWGQGTQVTVSS"),
        "vl": None,
        "source": "사용자 제공 · PDB 3U85. 표적 주석 원문 대조 권장.",
    },
    "FC44 (BBB 투과 VHH)": {
        "target": "BBB 수용체 (원문 확인)", "fmt": "VHH", "valency": 1, "kd_nM": None,
        "vh": ("EVQLQASGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSSINNGGGSTYY"
               "ADSVKGRFTISRDNAKNTLYLQMNSLKPEDTAVYYCAKDRLVTAYYFDYWGQGTQVTVSS"),
        "vl": None,
        "source": "사용자 제공 · 고투과 VHH. 표적 주석 원문 대조 권장.",
    },
    "HIRMAb 83-14 (항-인간 인슐린수용체)": {
        "target": "Human INSR", "fmt": "scFv", "valency": 1, "kd_nM": None,
        "vh": ("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVAVISYDGSNKYY"
               "ADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARWGLGYYFDYWGQGTLVTVSS"),
        "vl": ("DIQMTQSPSSLSASVGDRVTITCRASQGISSWLAWYQQKPEKAPKSLIYAASSLQSGVPS"
               "RFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYPPTFGGGTKVEIK"),
        "source": "사용자 제공 · Pardridge HIRMAb 83-14 계열(IGF/INSR RMT). 원문 대조 권장.",
    },
    # --- 서열 미공개, 실측 KD 보유 (sweet-spot 계산용) ---
    "BBB00533 (항-TfR VHH · KD 207nM)": {
        "target": "TfR", "fmt": "VHH", "vh": None, "vl": None, "kd_nM": 207.0,
        "valency": 1,
        "source": "PMC10300862 — 인간 TfR KD≈207 nM, monovalent. 서열 원문 미공개.",
    },
    "BBB00515 (항-TfR VHH · KD 1184nM)": {
        "target": "TfR", "fmt": "VHH", "vh": None, "vl": None, "kd_nM": 1184.0,
        "valency": 1,
        "source": "PMC10300862 — 인간 TfR KD≈1184 nM(저친화), monovalent. 서열 원문 미공개.",
    },
    "Trontinemab (monovalent 항-TfR)": {
        "target": "TfR", "fmt": "Fab", "vh": None, "vl": None, "kd_nM": None,
        "valency": 1,
        "source": "Roche brainshuttle — monovalent 설계(문헌). 서열·KD는 검증 후 입력.",
    },
}


@dataclass
class AntibodyShuttle:
    """항체 셔틀 표현. 서열(vh/vl)·kd는 DB 주입 전엔 None일 수 있다."""
    name: str
    target: str = "TfR"
    fmt: str = "VHH"                 # VHH(나노바디)/scFv/Fab
    vh: str | None = None            # VH 서열
    vl: str | None = None            # VL 서열(VHH면 None)
    cdrs: list | None = None         # [(chain, cdr_no, seq)] 선택
    kd_nM: float | None = None       # 실측/예측 친화도 (nM)
    valency: int = 1                 # 1=monovalent(권장), 2=bivalent
    source: str = ""                 # 출처(논문/DB) — 실재 확인 필수

    @property
    def sequence(self) -> str:
        """서열 기반 축 계산용 결합 서열(VH+VL). 없으면 빈 문자열."""
        return ((self.vh or "") + (self.vl or "")).upper()


@dataclass
class AntibodyAssessment:
    name: str
    target: str
    fmt: str
    valency: int
    # 서열 기반(가능 시)
    developability: dict | None = None
    solubility: dict | None = None
    net_charge: float | None = None
    # CDR(파라토프) 기반 — 결합 루프에 집중한 서열 축
    cdr_h3: str = ""
    cdr_l3: str = ""
    paratope: str = ""
    cdr_method: str = ""
    nonspecific: dict | None = None      # 파라토프 비특이 흡착 프록시(전하·친유성)
    deepb3p_cdr: float | None = None     # (선택) CDR CPP-유사도 = 비특이흡착 참고(⚠️BBB 투과도 아님)
    # 친화도 기반(Tier-3)
    sweetspot: float | None = None       # 0~1 (KD 있을 때만)
    kd_nM: float | None = None
    # 구조(훅)
    structure: dict | None = None
    # 결합가
    valency_note: str = ""
    # 서열 Tier 초기 랭킹(triage — 결합·구조 미포함)
    initial_rank: float = 0.0
    # 정직한 축 구분
    computable: list = field(default_factory=list)   # 실제 계산한 축
    pending: list = field(default_factory=list)       # Tier-3/DB 필요로 보류
    notes: str = ""


# ── 3) 친화도 sweet-spot 모델 ────────────────────────────────────────────────
def transcytosis_sweetspot(kd_nM: float | None,
                           optimum_nM: float = 100.0,
                           width_decades: float = 1.0) -> float | None:
    """TfR 친화도 sweet-spot의 **정성 프록시**(0~1). KD 없으면 None.

    문헌 관찰: TfR 친화도가 **너무 높으면(저 KD)** 뇌 쪽에서 방출되지 못해 갇히고, **너무 낮으면**
    결합 자체가 안 돼 뇌흡수가 준다 → 중간~약 친화가 유리. 로그 KD에 대한 종형(bell) 근사이며,
    optimum_nM·width는 **튜닝 파라미터**(실측으로 보정해야 함, 지금은 기본 휴리스틱).
    """
    if kd_nM is None or kd_nM <= 0:
        return None
    x = math.log10(kd_nM) - math.log10(optimum_nM)
    return round(math.exp(-(x * x) / (2 * width_decades * width_decades)), 3)


# ── 4) 구조 예측 훅 (외부·미장착) ────────────────────────────────────────────
def predict_antibody_structure(ab: AntibodyShuttle) -> dict:
    """항체 구조 예측 플러그 지점. IgFold/ABodyBuilder(단일) 또는 Fab–TfR co-folding
    (AlphaFold-Multimer/Boltz)을 연동하면 이 함수만 구현하면 된다. 미장착이면 available=False."""
    return {
        "available": False,
        "note": ("항체 구조 예측 미장착 — 단일 항체는 IgFold/ABodyBuilder, 복합체 결합 타당성은 "
                 "Fab–{0} co-folding(AF-Multimer/Boltz) 연동 필요(Tier-3).").format(ab.target),
    }


# ── 5) 결합가 로직 ───────────────────────────────────────────────────────────
def _valency_note(valency: int, target: str) -> str:
    if valency >= 2:
        return (f"bivalent({valency}가) — {target} 브레인셔틀은 **monovalent가 수송에 유리**한 것이 "
                "문헌상 알려져 있어 감점 요인(수용체 가교→분해경로 유도 위험).")
    return "monovalent(1가) — 수용체 매개 수송에 유리한 권장 형식."


# ── 6) 파이프라인 집계 ───────────────────────────────────────────────────────
def assess_antibody_shuttle(ab: AntibodyShuttle, bbb_scorer=None) -> AntibodyAssessment:
    """항체 셔틀을 서열 Tier 축으로 평가하고, 보류 축을 정직하게 표시해 반환.

    bbb_scorer: (선택) seq→0~1 콜러블(deepB3P 등). 주어지면 **CDR(파라토프)의 CPP-유사 성향**을
                참고 축으로 계산 — ⚠️ 이는 '비특이 흡착 프록시'이지 BBB 투과도 예측이 아니다
                (높을수록 끈적 = off-target 우려). 결합/전달은 친화도+구조(Tier-3)가 결정한다.
    """
    res = AntibodyAssessment(name=ab.name, target=ab.target, fmt=ab.fmt, valency=ab.valency,
                             kd_nM=ab.kd_nM)

    # 1) CDR3(파라토프) 근사 추출 — 특이성 축은 결합 루프에 집중해 본다.
    cd = extract_cdrs(ab.vh, ab.vl)
    res.cdr_h3, res.cdr_l3 = cd["h3"], cd["l3"]
    res.paratope, res.cdr_method = cd["paratope"], cd["method"]

    # 2) 개발성·용해도 — Fv 전체(liability는 어디서든 의미, 용해도는 분자 물성). 엔진은 펩타이드 보정 근사.
    seq = ab.sequence
    if seq:
        from .developability import assess_developability
        from .solubility import assess_solubility
        dev = assess_developability(seq)
        sol = assess_solubility(seq)
        if not dev.error:
            res.developability = {"risk": dev.risk_level, "n_liab": dev.n_liabilities,
                                  "liabilities": dev.liabilities, "agg": dev.agg_score}
            res.net_charge = dev.net_charge
            res.computable.append("개발성(Fv)")
        if not sol.error:
            res.solubility = {"score": sol.score, "level": sol.level}
            res.computable.append("용해도")
        res.notes += "⚠️ 서열 기반 축은 펩타이드 보정 엔진의 근사(항체 특화 TAP 향후). "
    else:
        res.pending.append("개발성·용해도(서열 DB 주입 후 계산)")

    # 3) 비특이 흡착 프록시 — **파라토프(CDR)의 전하/친유성**. 결합 루프가 끈적하면 off-target·빠른 청소.
    if res.paratope:
        from .selectivity import assess_selectivity
        sel = assess_selectivity(res.paratope)
        res.nonspecific = {"off_risk": sel.off_target_risk, "specificity": sel.selectivity,
                           "level": sel.risk_level, "drivers": sel.drivers}
        res.computable.append("비특이 흡착(파라토프)")

    # 4) (선택) CDR deepB3P = **CPP-유사 비특이 흡착 프록시** — BBB 투과도 아님(높을수록 끈적).
    if bbb_scorer and res.paratope:
        try:
            res.deepb3p_cdr = round(float(bbb_scorer(res.paratope)), 3)
            res.computable.append("CDR CPP-유사도(참고·비특이흡착)")
        except Exception:  # noqa: BLE001
            pass

    # 5) 친화도 sweet-spot — KD 있을 때만.
    ss = transcytosis_sweetspot(ab.kd_nM)
    if ss is not None:
        res.sweetspot = ss
        res.computable.append("친화도 sweet-spot")
    else:
        res.pending.append("친화도 sweet-spot(KD 실측/예측 필요, Tier-3)")

    # 6) 구조 — 훅.
    res.structure = predict_antibody_structure(ab)
    if not res.structure.get("available"):
        res.pending.append("구조/결합 타당성(IgFold·AF-Multimer 연동, Tier-3)")

    # 7) 결합가.
    res.valency_note = _valency_note(ab.valency, ab.target)
    res.computable.append("결합가")

    # 8) 서열 Tier 초기 랭킹(triage) — 결합 친화도·구조 미포함(Tier-3에서 보강).
    res.initial_rank = _initial_rank(res)

    # BBB 통과 절대값은 서열로 예측 불가 — 항상 보류(Tier-3 핵심: 친화도+구조).
    res.pending.append("BBB 통과/뇌흡수 절대값(서열 예측 불가 — 친화도+구조 기반 Tier-3)")
    return res


def _initial_rank(res: AntibodyAssessment) -> float:
    """서열 Tier 초기 랭킹(0~1, 높을수록 유망). 결합·구조 미포함 — triage 순서용.

    특이성(비특이 흡착 낮음) 0.35 · 개발성 0.25 · 용해도 0.20 · 결합가 0.10 · avidity 0.10.
    (deepB3P-CDR은 '참고' 축이라 점수에 넣지 않는다.)
    """
    spec = res.nonspecific["specificity"] if res.nonspecific else 0.5
    devmap = {"낮음": 1.0, "보통": 0.6, "높음": 0.25}
    dev = devmap.get(res.developability["risk"], 0.5) if res.developability else 0.5
    sol = res.solubility["score"] if res.solubility else 0.5
    val = 1.0 if res.valency < 2 else 0.6
    avid = res.sweetspot if res.sweetspot is not None else 0.5
    score = 0.35 * spec + 0.25 * dev + 0.20 * sol + 0.10 * val + 0.10 * avid
    return round(max(0.0, min(1.0, score)), 3)


def registry_shuttles() -> list[AntibodyShuttle]:
    """ANTIBODY_SHUTTLES(설정) → AntibodyShuttle 목록. DB 주입 전엔 빈 목록."""
    out = []
    for name, v in ANTIBODY_SHUTTLES.items():
        out.append(AntibodyShuttle(
            name=name, target=v.get("target", "TfR"), fmt=v.get("fmt", "VHH"),
            vh=v.get("vh"), vl=v.get("vl"), kd_nM=v.get("kd_nM"),
            valency=int(v.get("valency", 1)), source=v.get("source", "")))
    return out


def rank_antibody_shuttles(shuttles: list | None = None, bbb_scorer=None) -> list:
    """항체 셔틀들을 **서열 Tier 초기 랭킹**으로 정렬(결합·구조 미포함 triage).

    서열(vh) 있는 것만 평가한다. bbb_scorer를 주면 CDR CPP-유사도(비특이흡착 참고)도 채운다.
    ⚠️ 이 랭킹은 BBB 투과 확률이 아니라 '만들기 좋고 비특이 흡착 낮은' 후보의 우선순위다.
    실제 뇌 전달은 친화도(sweet-spot)+구조(Tier-3)가 확정한다.
    """
    abs_list = shuttles if shuttles is not None else registry_shuttles()
    assessed = [assess_antibody_shuttle(a, bbb_scorer=bbb_scorer) for a in abs_list if a.sequence]
    return sorted(assessed, key=lambda r: r.initial_rank, reverse=True)
