#!/usr/bin/env python3
"""LLM·사용자가 읽는 **모든 문자열**에 폐기·반증된 주장이 남아 있는지 전수 점검한다 (§67).

왜 이 스크립트가 존재하나
--------------------------------------------------------------------------
폐기된 수치가 프롬프트에 살아남는 사고가 **다섯 번** 났다
(§48-3 · §60 · §61 · §67 에서 둘). 네 번은 "눈에 띈 한 곳"을 고치는 식으로
처리했고, 그래서 다섯 번째가 또 났다. 고칠 대상은 문자열이 아니라 **점검 절차**다.

두 겹으로 본다 — 어느 한 겹만으로는 빠진다.
  ⑴ **AST 전수** — `ast.walk` 로 모든 문자열 리터럴을 꺼낸다. 함수 안인지 독스트링인지
     도구 설명인지 예외 메시지인지 가리지 않는다. §67 에서 나온 29곳 중 **절반이
     여기서만** 걸렸다(도구 응답 문자열 · `ValueError` 메시지).
  ⑵ **실제 입력 덤프** — 에이전트를 실제로 띄워 설계자 시스템 프롬프트 · 심사자 프롬프트 ·
     도구 선언 설명을 그대로 찍는다. f-string 으로 조립된 **최종 모습**은 여기서만 보인다.

그리고 §67-2 의 교훈: **LLM 의 출력을 읽는 쪽(측정 코드)도 같은 오염에 노출된다.**
`GATE_MARKERS` 가 심사자 프롬프트의 어휘와 겹치면 LLM 산문이 마커로 오인된다 → ⑶ 에서 본다.

재현: /usr/bin/python3 scripts/audit_prompt_strings.py
"""
from __future__ import annotations

import ast
import io
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

def _targets() -> list:
    """점검 대상 — **전 파일이다.**

    🔴 2026-11-16 정정 — 처음 쓴 이 스크립트는 파일 10개를 **손으로 적은 목록**으로
    돌면서 스스로 "전수"라고 출력했다. 실제로는 59개 중 10개였고, 빠진 것 중에
    **`app.py`(사용자가 보는 화면 전부)** 가 있었다. 거기서 9곳이 나왔다 —
    데모 심사문이 *"BBB·독성·**선택성** 근거가 일관되고"* 로 승인을 찍고 있었고,
    후보 카드가 *"**선택성** 높음"* 을 그대로 보여주고 있었다.
    → **목록을 손으로 적지 않는다.** 디렉터리를 훑는다.
    """
    out = []
    for rel in ("app.py", "benchmark.py", "benchmark_agent.py"):
        if os.path.exists(os.path.join(ROOT, rel)):
            out.append(rel)
    for d in ("core", "scripts"):
        dp = os.path.join(ROOT, d)
        if not os.path.isdir(dp):
            continue
        for f in sorted(os.listdir(dp)):
            if f.endswith(".py") and f != "__init__.py":
                out.append(f"{d}/{f}")
    # 이 스크립트 자신은 뺀다 — `BAD` 목록 자체가 전부 적중한다(자기 참조 21건).
    return [f for f in out if not f.endswith("audit_prompt_strings.py")]

# 폐기·반증된 주장을 담은 말들. **늘려야 하는 목록이다** — 새 측정이 무언가를 깰 때마다
# 그 말을 여기 추가하는 것이 §67 절차의 핵심이고, 그래야 다음에 자동으로 걸린다.
BAD = [
    "비독성", "독성 게이트", "독성뿐",                 # §47-6 조성 필터
    "기능 보존", "기능(rmt_sim)", "검증 기능", "결합 모티프를 보존",   # §55 상속 게이트
    "선택성", "선택적",                                # §50 off 축은 전하 하나
    "BBB 하한",                                        # 2026-11-10 폐지
    "무작위보다 나쁘", "상대 비교가 신뢰",             # §53-5 · §57 deepB3P
    "0/108", "5/5 지배",                               # §53 폐기
    "AGGRESCAN 핫스팟",                                # §41 이름 불일치
    "6세션", "3종 중 2종", "0/3", "4/9",               # §64 · §65 · §67 중간값
    "BBB 를 예측", "BBB를 예측",                       # §56 M6
]

# 같은 문맥에 이 표시가 있으면 **정정문**이므로 통과시킨다.
OK_MARK = ["🔴", "⚰️", "§", "아니다", "없다", "말 것", "쓰지 마", "않는다", "폐기", "정정", "옛 ", "금지"]

# 🔴 숫자 비(比)는 **부분 문자열로 찾으면 안 된다** — `"0/3"` 이 `"0/37"`(§50 수용체
#    발현 대조)에 걸려 위양성 2건을 냈다. 뒤에 숫자가 붙으면 다른 수다.
_NUMERIC = {"0/3", "4/9", "0/108", "6세션", "3종 중 2종"}


def _real_hit(text: str, bad: str, at: int) -> bool:
    if bad not in _NUMERIC:
        return True
    tail = text[at + len(bad): at + len(bad) + 1]
    return not tail.isdigit()


def _string_nodes(tree) -> list:
    """문자열을 `(텍스트, 줄번호)` 로 꺼낸다 — **f-string 은 통째로 이어 붙인다.**

    🔴 2026-11-16 정정 — 처음에는 `ast.Constant` 만 훑었다. 그런데 파이썬은 f-string 을
    `JoinedStr` 로 **조각내서** 담는다: `f"BBB 하한({x}) 통과 ⚰️ 폐지됨"` 은
    `"BBB 하한("` · `") 통과 ⚰️ 폐지됨"` 두 조각이 되고, 조각 하나만 보면 **같은 줄에 있는
    정정 표시를 못 본다.** 실제로 `old_structure_angiopep.py` 가 묘비(⚰️)를 달았는데도
    계속 적중으로 떴다. f-string 은 이어 붙여서 **한 문장으로** 판단한다.
    """
    out, inner = [], set()
    for node in ast.walk(tree):
        if isinstance(node, ast.JoinedStr):
            parts = []
            for v in node.values:
                if isinstance(v, ast.Constant) and isinstance(v.value, str):
                    parts.append(v.value)
                    inner.add(id(v))        # 조각은 따로 세지 않는다
            if parts:
                out.append(("{}".join(parts), node.lineno))
    for node in ast.walk(tree):
        if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                and id(node) not in inner):
            out.append((node.value, node.lineno))
    return out


def scan_ast() -> list:
    hits = []
    for rel in _targets():
        path = os.path.join(ROOT, rel)
        if not os.path.exists(path):
            print(f"  ⚠️ 없음: {rel}")
            continue
        src = io.open(path, encoding="utf-8").read()
        for node in _string_nodes(ast.parse(src)):
            v, lineno = node
            for bad in BAD:
                i = v.find(bad)
                if i < 0:
                    continue
                if not _real_hit(v, bad, i):
                    continue
                ctx = v[max(0, i - 120): i + 180]
                if any(o in ctx for o in OK_MARK):
                    continue
                hits.append((rel, lineno, bad, ctx.replace("\n", " | ")[:220]))
    return hits


def scan_markers() -> list:
    """⑶ 측정 도구 오염 — 게이트 마커가 심사자 프롬프트의 **어휘와 겹치면** LLM 이 따라 쓴다."""
    try:
        from benchmark_agent import GATE_MARKERS
    except Exception as exc:                    # pragma: no cover
        print(f"  ⚠️ GATE_MARKERS 를 못 읽었다: {exc}")
        return []
    src = io.open(os.path.join(ROOT, "core/optimizer_agent_gemini.py"), encoding="utf-8").read()
    i = src.find("def _critic_prompt")
    j = src.find("\n    def ", i + 10)
    prompt_src = src[i:j] if i >= 0 else ""
    return [(nm, mk) for nm, mk in GATE_MARKERS if mk in prompt_src]


def main() -> int:
    print("=" * 78)
    print("§67 프롬프트 문자열 전수 점검")
    print("=" * 78)

    tgt = _targets()
    print(f"\n[⑴ AST 전수] 파일 {len(tgt)}개 × 폐기어 {len(BAD)}종")
    hits = scan_ast()
    if hits:
        print(f"  🔴 적중 {len(hits)}건 — 하나씩 보고 **정정문인지 미정정인지** 판단하라:")
        for rel, ln, bad, ctx in hits:
            print(f"    {rel}:{ln}  [{bad}]\n       {ctx}")
    else:
        print("  ✅ 미정정 0건")

    print("\n[⑶ 측정 도구 오염] GATE_MARKERS × 심사자 시스템 프롬프트")
    bad_mk = scan_markers()
    if bad_mk:
        print("  🔴 심사자 프롬프트가 **가르치는 말**을 마커로 쓰고 있다 — LLM 산문에 걸린다:")
        for nm, mk in bad_mk:
            print(f"    {nm:12s} ← {mk!r}  → 게이트 **고유 머리말**로 바꿔라")
    else:
        print("  ✅ 오염 0건 (전부 고유 머리말)")

    print("\n[⑵ 실제 입력 덤프] — 이 스크립트가 못 보는 겹")
    print("  f-string 으로 조립된 **최종 프롬프트**는 에이전트를 띄워야 보인다.")
    print("  benchmark_agent.py --dry-run-m5 로 게이트 표를, §67-7 의 덤프로 프롬프트를 확인한다.")

    print("\n" + "=" * 78)
    print("판정: ⑴ 적중은 **자동 실패가 아니다** — 정정문도 같이 걸린다. 눈으로 가른다.")
    print("      ⑶ 적중은 **무조건 버그다** — 측정값이 오염된다.")
    print("=" * 78)
    return 1 if bad_mk else 0


if __name__ == "__main__":
    raise SystemExit(main())
