#!/usr/bin/env python3
"""도구 4개가 **부르면 실제로 돌아가는가** — LLM 0회, 유료 호출 0회.

왜 이 스크립트가 있나 (§49)
---------------------------
2026-11-12, `optimize_developability` 가 **호출될 때마다 100% 터지고 있었다**.
원인은 한 줄이었다 — §26 이 off-target 을 목적 벡터에서 빼 **밴드 제약**으로 옮기면서
`objectives()` 가 4-튜플에서 3-튜플이 됐는데, 그 벡터를 출력하는 문자열이
`seed[3]` 을 **off 로 계속 읽고 있었다**(`IndexError: tuple index out of range`).

그 실패는 **조용했다**:

* 예외가 `run()` 의 `except` 로 올라가 화면에는 **한 줄 요약**만 떴고,
* 요약은 접힌 status 안에 들어가 **페이지 본문에서 보이지도 않았으며**,
* 세션은 "최적화 완료"로 끝나고 **포디움 3개가 정상적으로 그려졌다**.
  사라진 것은 결과가 아니라 **탐색**이었다.

그리고 18시간 동안 아무도 몰랐다. 이 도구를 실제로 밟는 경로가 **유료 LLM 세션뿐**
이었기 때문이다. 그래서 **LLM 없이 도구를 직접 밟는** 이 스모크 테스트를 만든다.
`_dispatch` 는 모델이 부르는 것과 **같은 입구**이므로, 여기서 통과하면
"도구가 죽어 있는" 종류의 실패는 유료 세션 전에 잡힌다.

🔴 이것은 **정확성 테스트가 아니다.** 값이 맞는지는 보지 않는다 — 터지지 않고,
빈 문자열이 아니고, 핵심 키가 들어 있는지만 본다(= 조용한 죽음만 잡는다).
"""
from __future__ import annotations

import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import LINKER_LIBRARY, SHUTTLES, get_settings  # noqa: E402
from core.optimizer_agent_gemini import GeminiOptimizationAgent  # noqa: E402

CARGO = "GSNKGAIIGLM"
SHUTTLE = SHUTTLES["Angiopep-2"]["seq"]
LINKER = LINKER_LIBRARY["GGGGS"]["seq"] if "GGGGS" in LINKER_LIBRARY else "GGGGS"

CASES = [
    ("evaluate_candidates", {
        "rationale": "스모크 — 라이브러리 조합 1개",
        "constructs": [{"label": "smoke", "linker": LINKER, "shuttle": SHUTTLE}],
    }, ["liab", "|"]),
    ("design_candidate", {
        "label": "smoke-design", "linker": LINKER, "shuttle": SHUTTLE,
        "hypothesis": "스모크 — 편집 없이 원본 제출",
        "prediction": "지표 변화 없음",
    }, ["근거", "반성"]),
    ("optimize_developability", {"shuttle": SHUTTLE, "linker": LINKER},
     ["전선", "rmt_sim"]),
]


def main() -> int:
    agent = GeminiOptimizationAgent(get_settings(), max_rounds=6, use_ledger=False)
    print("=" * 72)
    print("도구 스모크 — LLM 0회 (`_dispatch` 를 모델과 같은 입구로 직접 부른다)")
    print("=" * 72)
    bad = 0
    for name, args, must in CASES:
        try:
            text, _ev = agent._dispatch(CARGO, name, dict(args))
        except Exception:  # noqa: BLE001
            print(f"🔴 {name}: 예외")
            traceback.print_exc()
            bad += 1
            continue
        if not (text or "").strip():
            print(f"🔴 {name}: 빈 응답")
            bad += 1
            continue
        missing = [m for m in must if m not in text]
        if missing:
            print(f"🔴 {name}: 응답에 {missing} 없음 — 내용이 비었을 수 있다")
            print("   ↳ " + text[:200].replace("\n", " "))
            bad += 1
            continue
        print(f"✅ {name}: {len(text)}자 · 핵심 키 {must} 확인")
    print("-" * 72)
    print("모두 통과" if not bad else f"🔴 실패 {bad}/{len(CASES)}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
