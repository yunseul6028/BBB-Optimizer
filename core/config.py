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

# 로컬 deepB3P (펩타이드 BBB 예측) 자원 경로
DEEPB3P_REPO = BASE_DIR / "vendor" / "deepB3P"
DEEPB3P_PYTHON = BASE_DIR / "vendor" / ".venv-deepb3p" / "bin" / "python"
DEEPB3P_RUNNER = DEEPB3P_REPO / "_run_predict.py"
DEEPB3P_CKPT = DEEPB3P_REPO / "model" / "checkpoint" / "512d_16ff_32k_1nl_2h_0.0001lr_0.1p"

# 로컬 FBGAN 생성 최적화 (잠재공간 진화 루프) 자원 경로
FBGAN_RUNNER = DEEPB3P_REPO / "_run_fbgan.py"
FBGAN_GEN_CKPT = DEEPB3P_REPO / "fbgan" / "checkpoint" / "G_weights_1000.pth"

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
MODEL_MAX_LEN = 50                 # deepB3P seq_len — 초과분은 잘림(truncate)
TOP_N = 3                          # 상위 몇 개 조합을 최종 추천으로 보여줄지

# --- 링커 라이브러리 (융합체의 링커 축, 전수 탐색) ---
STANDARD_LINKER_NAME = "GGGGGS"
LINKER_LIBRARY = {
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
SHUTTLES = {
    "Angiopep-2": {"seq": "TFFYGGSRGKRNNFKTEEY", "target": "LRP1",       "note": "LRP1 타겟 대표 BBB 셔틀(19aa)"},
    "TAT":        {"seq": "GRKKRRQRRRPPQ",       "target": "세포투과",    "note": "HIV-1 TAT 세포투과 펩타이드(13aa)"},
    "Penetratin": {"seq": "RQIKIWFQNRRMKWKK",    "target": "세포투과",    "note": "Antennapedia 세포투과 펩타이드(16aa)"},
    "SynB1":      {"seq": "RGGRLSYSRRRFSTSTGR",  "target": "세포투과",    "note": "SynB 벡터 계열(18aa)"},
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
    llm_api_key: str | None = None
    llm_model: str = "claude-opus-4-8"
    toxicity_threshold: float = TOXICITY_THRESHOLD
    top_n: int = TOP_N

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
    def use_fbgan_local(self) -> bool:
        """로컬 FBGAN 생성 최적화 사용 가능 여부 (생성기 가중치 + deepB3P + ToxinPred3)."""
        return (
            self.use_deepb3p_local
            and self.use_toxinpred3_local
            and FBGAN_RUNNER.exists()
            and FBGAN_GEN_CKPT.exists()
        )

    @property
    def use_real_predictor(self) -> bool:
        return bool(self.deepb3_api_url and self.deepb3_api_key)

    @property
    def use_llm_strategy(self) -> bool:
        return bool(self.llm_api_key)

    @property
    def use_llm_agent(self) -> bool:
        """자율 가설 에이전트(LLM 뇌 + 로컬 엔진 검증) 사용 가능 여부."""
        return bool(self.llm_api_key) and self.use_deepb3p_local and self.use_toxinpred3_local


def get_settings() -> Settings:
    return Settings(
        deepb3_api_url=os.getenv("DEEPB3_API_URL"),
        deepb3_api_key=os.getenv("DEEPB3_API_KEY"),
        llm_api_key=os.getenv("LLM_API_KEY") or os.getenv("ANTHROPIC_API_KEY"),
        llm_model=os.getenv("LLM_MODEL", "claude-opus-4-8"),
    )
