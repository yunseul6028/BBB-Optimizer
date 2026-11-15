#!/usr/bin/env python3
# AUDIT-EXEMPT: 이 파일은 금지어 목록 자체를 담는다 (점검 대상에서 제외)
"""금지어 목록을 **`docs/수치원장.md` Ⅳ 표에서 자동 추출**한다 (§67-11).

왜 이 모듈이 있나 — §67 의 마지막 뿌리
--------------------------------------------------------------------------
§67 에서 같은 실수를 **네 번** 했다: 파일 목록(59개 중 10개) · 폐기어 목록
(`선택성 축` 누락) · 훑는 깊이(`docs/figs/`) · 자기 제외(파일 이름 하나).
§65 의 GA 대조 목록까지 치면 다섯 번이고 전부 모양이 같다 —
**내가 손으로 적어 둔 것이 세상보다 늦게 자란다.**

앞의 넷은 `os.walk` 와 역할 기반 제외로 닫았다. 남은 하나가 **폐기어 목록**이었다.
그 목록은 이미 다른 곳에 **정식으로** 존재한다: `docs/수치원장.md` **Ⅳ. 🚫 쓰면 안 되는
숫자**. 새 측정이 무언가를 깰 때 그 표에 한 행을 적는 것이 이 프로젝트의 규율이다.
→ 그러면 **점검기가 그 표를 읽으면 된다.** 표가 자라면 점검도 같이 자란다.

추출 규칙
--------------------------------------------------------------------------
Ⅳ 표 1열에서 ⑴ `"..."` 안의 문구 ⑵ `**...**` 안의 문구 ⑶ `` `...` `` 안의 문구를
꺼낸다. 🔴 🆕 같은 표시와 꼬리 괄호는 떼고, **2글자 미만·30글자 초과는 버린다**
(너무 짧으면 아무 데나 걸리고, 너무 길면 어디에도 안 걸려 조용히 무용지물이 된다).

여기에 **`_NAMES`(손으로 적은 보충)** 를 더한다. 이건 원장 표가 "문장"으로 적는 것들의
**이름 형태**다 — 예컨대 표에는 *"게이트 통과 = 비독성"* 이라고 적히지만 코드에서
찾아야 하는 것은 `비독성` 네 글자다. 🔴 **`_NAMES` 는 자라지 않는 목록이라는 것을
알고 쓴다** — 새 폐기는 `_NAMES` 가 아니라 **원장 Ⅳ 표**에 적는다.
"""
from __future__ import annotations

import io
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEDGER = os.path.join(ROOT, "docs", "수치원장.md")

# 🔴 손으로 적은 보충 — **원장 표가 문장으로 적는 것의 "이름" 형태**만 둔다.
#    새로 무언가가 깨지면 여기가 아니라 **원장 Ⅳ 표**에 적는다.
_NAMES = [
    "비독성", "독성 게이트", "독성뿐",            # §47-6 → 조성 필터
    "기능 보존", "기능(rmt_sim)", "검증 기능", "결합 모티프를 보존",   # §55 → 부모 유사도 하한
    "선택성 축", "선택성 AUC", "수용체 선택적", "선택성↑",             # §50 → 전하 기반 off-target
    "상대 비교가 신뢰", "무작위보다 나쁘",       # §53-5 · §57
    "BBB 하한",                                   # 2026-11-10 폐지
]

# 숫자 비(比)는 **부분 문자열로 찾으면 안 된다** — `0/3` 이 `0/37` 에 걸렸다(§67-8).
_NUMERIC_RE = re.compile(r"^[\d]+\s*/\s*[\d]+$|^\d+세션$")

_STRIP = ("🔴", "🆕", "🚫", "⚰️", "⚠️", "**", "“", "”", '"', "`", " ")


def _clean(t: str) -> str:
    t = t.strip()
    for ch in ("🔴", "🆕", "🚫", "⚰️", "⚠️"):
        t = t.replace(ch, "")
    t = re.sub(r"\([^)]*\)\s*$", "", t)          # 꼬리 괄호 (Pareto) 등
    return t.strip().strip("“”\"`").strip()


def from_ledger() -> list:
    """원장 Ⅳ 표 1열에서 금지 문구를 꺼낸다. 원장이 없으면 빈 목록(조용히 실패하지 않게 호출부에서 센다)."""
    if not os.path.exists(LEDGER):
        return []
    lines = io.open(LEDGER, encoding="utf-8").read().split("\n")
    try:
        i = next(k for k, l in enumerate(lines) if l.startswith("## Ⅳ"))
    except StopIteration:
        return []
    j = next((k for k, l in enumerate(lines[i + 1:], i + 1) if l.startswith("## ")), len(lines))

    out = []
    for line in lines[i:j]:
        if not line.startswith("|") or re.match(r"^\|[\s\-|]+\|$", line):
            continue
        col1 = line.split("|")[1]
        if _clean(col1) in ("숫자", ""):
            continue
        # ⑴ 따옴표 ⑵ 굵게 ⑶ 코드 — 셋 다 훑는다
        for pat in (r"[“\"]([^”\"]{2,60})[”\"]", r"\*\*([^*]{2,60})\*\*", r"`([^`]{2,60})`"):
            for m in re.findall(pat, col1):
                t = _clean(m)
                if not (2 <= len(t) <= 30):
                    continue
                # 🔴 잡음 거르기 — 실측으로 둘이 나왔다:
                #  ⑴ `front_top` 같은 **코드 식별자**. 금지 문구가 아니라 그 문구가
                #     언급한 변수다. 그대로 두면 정상 코드 전체가 적중한다.
                #  ⑵ `§56 의 쌍 표 "2/6` 처럼 **§ 참조가 섞인 조각**. 따옴표가
                #     짝이 안 맞아 잘린 것이고, 문구가 아니라 출처다.
                if re.fullmatch(r"[a-z_][a-z0-9_]*", t):
                    continue
                if "§" in t:
                    continue
                out.append(t)
    # 중복 제거(순서 유지)
    seen, uniq = set(), []
    for t in out:
        if t not in seen:
            seen.add(t)
            uniq.append(t)
    return uniq


def load() -> tuple:
    """`(금지어 목록, 숫자비 집합, 출처별 개수)` 를 돌려준다."""
    led = from_ledger()
    terms, seen = [], set()
    for t in led + _NAMES:
        if t not in seen:
            seen.add(t)
            terms.append(t)
    numeric = {t for t in terms if _NUMERIC_RE.match(t)}
    return terms, numeric, {"원장 Ⅳ": len(led), "손 보충": len(_NAMES), "합계(중복 제거)": len(terms)}


if __name__ == "__main__":
    terms, numeric, src = load()
    print(f"출처: {src}")
    print(f"숫자 비(경계 검사 대상) {len(numeric)}개: {sorted(numeric)}")
    print("\n전체 목록:")
    for t in terms:
        print(f"  {'[숫자]' if t in numeric else '      '} {t}")
