"""설정 및 상수.

로컬 deepB3P(BBB) / ToxinPred3(독성)가 감지되면 자동으로 실측 모드로 전환.
링커·셔틀은 '라이브러리'로 정의하고, 에이전트가 전수 조합(링커 × 셔틀)을 탐색한다.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

# 프로젝트 루트 (core/config.py → core → 루트)
BASE_DIR = Path(__file__).resolve().parent.parent


def _load_dotenv(path: Path) -> None:
    """의존성 없이 .env를 os.environ에 로드한다.

    - 이미 설정된 실제 환경변수는 유지(setdefault).
    - 중복 키는 **먼저 나온 non-empty 값**이 이긴다(빈 줄이 실제 값을 덮지 않도록).
    - 인라인 주석(#)·양쪽 따옴표 제거.
    """
    if not path.exists():
        return
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, val = line.split("=", 1)
        key, val = key.strip(), val.strip()
        if val and val[0] not in "\"'" and "#" in val:  # 인라인 주석 제거
            val = val.split("#", 1)[0].strip()
        if len(val) >= 2 and val[0] in "\"'" and val[-1] == val[0]:  # 따옴표 제거
            val = val[1:-1]
        if not val:  # 빈 값은 무시 → 중복 빈 줄이 실제 키를 덮어쓰지 않음
            continue
        os.environ.setdefault(key, val)


_load_dotenv(BASE_DIR / ".env")

# 로컬 deepB3P (펩타이드 BBB 예측) 자원 경로
DEEPB3P_REPO = BASE_DIR / "vendor" / "deepB3P"
DEEPB3P_PYTHON = BASE_DIR / "vendor" / ".venv-deepb3p" / "bin" / "python"
DEEPB3P_RUNNER = DEEPB3P_REPO / "_run_predict.py"
DEEPB3P_CKPT = DEEPB3P_REPO / "model" / "checkpoint" / "512d_16ff_32k_1nl_2h_0.0001lr_0.1p"

# 로컬 ToxinPred3 (펩타이드 독성 예측) 자원 경로
TOXINPRED3_REPO = BASE_DIR / "vendor" / "toxinpred3"
TOXINPRED3_PYTHON = BASE_DIR / "vendor" / ".venv-toxinpred3" / "bin" / "python"
TOXINPRED3_SCRIPT = TOXINPRED3_REPO / "toxinpred3.py"
TOXINPRED3_MODEL_PKL = TOXINPRED3_REPO / "model" / "toxinpred3.0_model.pkl"
TOXINPRED3_MODEL_ARG = 1   # 1: ML-only(AAC+DPC, 외부도구 불필요), 2: Hybrid(MERCI 필요)

# ---------------------------------------------------------------------------
# 도메인 상수
# ---------------------------------------------------------------------------
VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
TOXICITY_THRESHOLD = 0.38          # 부작용 탈락 기준선 (ToxinPred3 기본 임계값에 맞춤)
# ⚰️ MIN_BBB_SIGNAL (= 0.15) 는 2026-11-10 **폐지**했다. 하지 말아야 할 일을 하고 있었다:
#    이름은 "최소 효능 신호"였지만 실제로 한 일은 **양이온성 필터**다. 실측(화물 GSNKGAIIGLM,
#    링커 13종 전조합)에서
#      · 문헌 근거 등급 **A**(ANG1005 임상 2상) Angiopep-2 → 최고 융합 BBB 0.102, **13/13 전멸**
#      · 등급 B ApoE(159-167)2 0.136 · Leptin30 0.006 → 역시 전멸
#      · 등급 **D**(BBB 통과 근거 없음) Tet1 0.568 → **통과**
#      · 상위 3위가 전부 CPP 비교군(SynB1 0.960 · TAT 0.929 · Penetratin 0.877)
#    즉 **근거가 가장 단단한 셔틀을 설계 공간에서 제거하고 근거 없는 셔틀을 통과시켰다.**
#    Angiopep-2 가 낮은 이유는 투과를 못 해서가 아니라 C말단이 산성(…KTEEY)이라서다.
#    근본 원인: M6(docs/experiments.md §14)에서 **융합 수준 deepB3P 에 서열 신호가 측정되지
#    않는다**가 확인됐으므로, 그 점수의 절대 임계값으로 후보를 탈락시킬 근거가 없다.
#    👉 대체: ⑴ 셔틀 선택 근거는 `core/evidence.py` 의 **문헌 근거 등급**(M7),
#            ⑵ near-zero 쓰레기 설계체는 **순서 특이성 측정 가능성** 하드 게이트가 막는다
#               (라이브러리 밖 셔틀이 homopolymer·반복서열이면 스크램블이 0개 = 측정 불가 → REVISE).
#    되살리지 말 것. 되살린다면 먼저 등급 A 셔틀이 통과하는지 확인할 것.

# 유효점수 λ — eff_bbb = BBB − λ·off_target. **후보군 내부**의 선택성 저울질용이다.
# (비교군 CPP 배제는 λ가 담당하지 않는다: 실측상 λ=0.25의 최대 영향력 0.25 < CPP↔RMT
#  BBB 격차 0.434 → 어떤 off 값으로도 역전 불가. 그래서 배제는 역할 분리가 맡는다.
#  optimizer_agent.is_comparison_only 참고.)
OFF_TARGET_PENALTY = 0.25
MODEL_MAX_LEN = 50                 # deepB3P seq_len — 초과분은 잘림(truncate)

# --- 링커 라이브러리 (융합체의 링커 축) ---
LINKER_LIBRARY = {
    "직접융합":    {"seq": "",               "kind": "none",        "note": "링커 없이 화물-셔틀 직접 연결(0aa, direct fusion)"},
    "GGS":        {"seq": "GGS",             "kind": "flexible",    "note": "짧은 유연 링커(3aa)"},
    "GGGGS":      {"seq": "GGGGS",           "kind": "flexible",    "note": "고전 GS 링커(5aa)"},
    "GGGGGS":     {"seq": "GGGGGS",          "kind": "flexible",    "note": "표준 유연 링커(6aa)"},
    "(GGGGS)2":   {"seq": "GGGGSGGGGS",      "kind": "flexible",    "note": "긴 유연 링커(10aa)"},
    "(GGGGS)3":   {"seq": "GGGGSGGGGSGGGGS", "kind": "flexible",    "note": "매우 긴 유연 링커(15aa)"},
    "YGGGGS":     {"seq": "YGGGGS",          "kind": "aromatic",    "note": "N-말단 타이로신(방향족)"},
    "WGGGGS":     {"seq": "WGGGGS",          "kind": "aromatic",    "note": "N-말단 트립토판(방향족)"},
    "FGGGGS":     {"seq": "FGGGGS",          "kind": "aromatic",    "note": "N-말단 페닐알라닌(방향족)"},
    "VGGGGS":     {"seq": "VGGGGS",          "kind": "hydrophobic", "note": "N-말단 발린(소수성)"},
    "EAAAK":      {"seq": "EAAAK",           "kind": "rigid",       "note": "α-헬릭스 강직 링커(5aa)"},
    "(EAAAK)2":   {"seq": "EAAAKEAAAK",      "kind": "rigid",       "note": "강직 링커(10aa)"},
    "A(EAAAK)3A": {"seq": "AEAAAKEAAAKEAAAKA", "kind": "rigid",     "note": "강직 헬릭스 링커(17aa)"},
}

# --- 셔틀 라이브러리 (BBB 투과 수송 펩타이드, 전수 탐색) ---
# ⚠️ `note` 문자열은 **기전 분류에 쓰인다**(`"CPP" in note` / `"RMT" in note` —
#    benchmark.py·bench_selectivity.py·directional_pairs.py·scramble_delta.py). 편집할 때
#    CPP/RMT 토큰을 지우거나 RMT 항목에 "CPP" 글자를 넣지 말 것.
# 📚 셔틀별 **문헌 전달 근거 등급·인용·측정 맥락**은 `core/evidence.py` 가 단일 출처다.
#    등급을 여기 복제하지 말 것(드리프트). 융합 BBB 점수로는 셔틀 간 순위를 매길 수 없다(M6).
SHUTTLES = {
    "Angiopep-2":     {"seq": "TFFYGGSRGKRNNFKTEEY", "target": "LRP1",       "note": "LRP1 타겟 대표 BBB 셔틀(19aa, RMT)"},
    "ApoE(159-167)2": {"seq": "LAVYQAGARLAVYQAGAR",  "target": "LDLR/LRP1",  "note": "ApoE 유래 RMT 셔틀·탠덤 이량체(18aa, 단량체 LAVYQAGAR = 성숙체 159-167). 🔴 이름 충돌: 일부 문헌이 같은 이름으로 (LRKLRKRLL)2(=성숙체 141-149, 전구체 번호 159-167)를 지칭한다 — 강양이온성 다른 분자이므로 혼동 금지"},
    "Leptin30":       {"seq": "YQQILTSMPSRNVIQISNDLENLRDLLHVL", "target": "렙틴수용체(LepR)", "note": "렙틴 61-90 유래 RMT 셔틀(30aa)"},
    "T7":             {"seq": "HAIYPRH",             "target": "TfR(트랜스페린 수용체)",   "note": "TfR 표적 셔틀·파지디스플레이(7aa, RMT) — Tf와 다른 결합부위로 내인성 Tf 비경쟁"},
    "B6":             {"seq": "GHKAKGPRK",           "target": "TfR1(트랜스페린 수용체)",  "note": "hTfR 세포외도메인 파지디스플레이 유래 monovalent 셔틀(9aa, RMT)"},
    "RVG29":          {"seq": "YTIWMPENPRPGTPCDIFTNSRGKRASNG", "target": "nAChR(니코틴성 수용체)", "note": "광견병 바이러스 당단백 유래 nAChR 셔틀(29aa, RMT)"},
    "CDX":            {"seq": "FKESWREARGTRIERG",    "target": "α7 nAChR(니코틴성 수용체)", "note": "Candoxin 유래 α7 nAChR 셔틀(16aa, RMT) — sub-μM 결합(IC50 84.5nM)·낮은 off-target. ⚠️ 근거 등급 C: 문헌에서 뇌 전달이 입증된 형태는 D-retro-inverso(DCDX)이고 이 라이브러리 서열은 덜 안정한 선형 L형"},
    "THR":            {"seq": "THRPPMWSPVWP",        "target": "TfR1(트랜스페린 수용체)",  "note": "TfR1 표적 셔틀·파지디스플레이 12-mer(12aa, RMT) — 금나노입자·AAV8·siRNA 뇌전달 보고. 주의: 혈청 안정성이 낮다는 보고(불안정성 지수도 높음) → 한계로 명시할 것"},
    "Tet1":           {"seq": "HLNILSTLWKYR",        "target": "GT1b 갱글리오사이드 수용체(뉴런)", "note": "파상풍독소 모사 뉴런 표적 셔틀(12aa, RMT) — 🔴 근거 등급 D: GT1b/스핑고미엘린 결합 뉴런 표적 리간드이고 BBB 통과 **단독 근거가 없다**. Tet1-LNP가 뇌가 아니라 말초 신경으로 주로 간다는 보고가 있고, 문헌에서는 거의 항상 BBB 투과 펩타이드(B6·TGN)와 병용된다. 비교·탐색용으로만 포함 — 자세한 인용은 core/evidence.py"},
    "TAT":            {"seq": "GRKKRRQRRRPPQ",       "target": "세포투과",    "note": "HIV-1 TAT 세포투과 펩타이드(13aa, CPP)"},
    "Penetratin":     {"seq": "RQIKIWFQNRRMKWKK",    "target": "세포투과",    "note": "Antennapedia 세포투과 펩타이드(16aa, CPP)"},
    "SynB1":          {"seq": "RGGRLSYSRRRFSTSTGR",  "target": "세포투과",    "note": "SynB 벡터 계열(18aa, CPP)"},
}

# 예시 화물(cargo) 펩타이드 — Aβ(25-35) 신경독성 단편(11aa, 알츠하이머 맥락). 데모용.
DEFAULT_CARGO = "GSNKGAIIGLM"


# 연결부위 계산 시 남길 화물 꼬리 길이(aa). 0 = 화물 완전 제외(링커+셔틀 = BBB 모듈).
# deepB3P는 친수성 화물 잔기에 민감해 flank>0이면 셔틀 신호가 희석되므로 0 권장.
JUNCTION_FLANK = 0


def junction_window(sequence: str, max_len: int = MODEL_MAX_LEN) -> str:
    """원시 서열 안전망: 최대 길이 초과 시 C말단(링커+셔틀) 윈도우를 취한다."""
    seq = sequence.upper()
    return seq if len(seq) <= max_len else seq[-max_len:]


def bbb_scoring_seq(cargo: str, linker: str, shuttle: str,
                    max_len: int = MODEL_MAX_LEN, flank: int = JUNCTION_FLANK) -> str:
    """BBB 계산에 넣을 서열. construct(=cargo+linker+shuttle)가 max_len 이하면 전체 그대로,
    초과하면 **연결부위(링커+셔틀=BBB 모듈, + 화물 꼬리 flank aa)** 만 남긴다. flank=0이면
    화물을 완전히 빼 셔틀 신호를 깨끗이 본다. (독성은 전체 서열로 별도 계산.)"""
    full = (cargo + linker + shuttle).upper()
    if len(full) <= max_len:
        return full
    tail = cargo[-flank:].upper() if flank > 0 else ""
    return (tail + linker + shuttle).upper()[-max_len:]


# ---------------------------------------------------------------------------
# 런타임 설정
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Settings:
    deepb3_api_url: str | None = None
    deepb3_api_key: str | None = None
    gemini_api_key: str | None = None
    gemini_model: str = "gemini-3-flash-preview"
    toxicity_threshold: float = TOXICITY_THRESHOLD

    @property
    def use_deepb3p_local(self) -> bool:
        return (
            DEEPB3P_PYTHON.exists()
            and DEEPB3P_RUNNER.exists()
            and (DEEPB3P_CKPT / "deepb3p_1.pth").exists()
        )

    @property
    def use_toxinpred3_local(self) -> bool:
        return (
            TOXINPRED3_PYTHON.exists()
            and TOXINPRED3_SCRIPT.exists()
            and TOXINPRED3_MODEL_PKL.exists()
        )

    @property
    def use_real_predictor(self) -> bool:
        return bool(self.deepb3_api_url and self.deepb3_api_key)

    @property
    def use_gemini_agent(self) -> bool:
        """Gemini 브레인 자율 에이전트 사용 가능 여부."""
        return bool(self.gemini_api_key) and self.use_deepb3p_local and self.use_toxinpred3_local


def get_settings() -> Settings:
    return Settings(
        deepb3_api_url=os.getenv("DEEPB3_API_URL"),
        deepb3_api_key=os.getenv("DEEPB3_API_KEY"),
        gemini_api_key=os.getenv("GEMINI_API_KEY") or os.getenv("GOOGLE_API_KEY"),
        gemini_model=os.getenv("GEMINI_MODEL", "gemini-3-flash-preview"),
    )
