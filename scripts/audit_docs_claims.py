#!/usr/bin/env python3
"""문서(`*.md`)에 **폐기·반증된 주장**이 현재형으로 남아 있는지 점검한다 (§67-9).

`scripts/audit_prompt_strings.py` 의 자매편이다. 그쪽은 **코드 안의 문자열**을 보고,
이쪽은 **사람이 읽는 산문**을 본다. 둘 다 필요한 이유는 §67 에서 다섯 번 증명됐다 —
고친 자리가 아니라 **안 본 자리**에서 계속 나왔다.

문서는 코드와 성격이 다르다
--------------------------------------------------------------------------
`docs/experiments.md` 는 **append-only 실험 로그**다. 거기에는 폐기된 주장이 **일부러**
남아 있다 — 무엇이 언제 왜 깨졌는지가 기록의 본체이기 때문이다. 그래서 로그는 **따로**
세고, 판정 기준도 다르다:

  · **로그(`experiments.md`)**     → 폐기 수치가 **정정 표시 없이** 떠 있으면 문제.
  · **현재형 문서**(README · 기획서 · 원장 · 흐름 · 논문개요 등)
                                   → 폐기 수치가 **있다는 것 자체**가 검토 대상.

판정은 자동이 아니다. 적중은 **"여기 사람이 봐야 한다"** 는 뜻이고, 정정문·역사 서술이면
통과다. 이 스크립트는 **눈이 가야 할 자리를 좁혀 줄 뿐**이다.

재현: /usr/bin/python3 scripts/audit_docs_claims.py
"""
from __future__ import annotations

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 로그 — 폐기 주장이 **기록으로서** 남아 있는 것이 정상인 파일
LOGS = {"docs/experiments.md", "docs/protocol.md"}

# 폐기·반증된 주장. 🔴 **늘려야 하는 목록이다** — 새 측정이 무언가를 깰 때마다 추가한다.
# (코드 쪽 `audit_prompt_strings.BAD` 와 **같이** 늘려야 한다.)
BAD = [
    # 이름 오칭
    "비독성", "독성 게이트", "기능 보존", "결합 모티프를 보존", "AGGRESCAN 핫스팟",
    # 반증된 주장
    "무작위보다 나쁘", "상대 비교가 신뢰", "BBB 하한", "BBB 를 예측", "BBB를 예측",
    "수용체 선택적", "선택성↑", "선택성 축", "선택성 AUC",
    # 🔴 이 세 줄이 §67 의 교훈이다 — **깨질 때마다 여기 추가한다.** 추가하면 로그([B])가
    #    시끄러워지는데, 그건 정상이다. [A](현재형 문서)가 조용한지가 판정 기준이다.
    # 폐기·갱신된 수치 (뒤에 숫자가 붙으면 다른 수 → 경계 검사)
    "0/108", "5/5 지배", "6/6 재현", "3종 중 2종", "13/15", "0/3", "4/9",
    # 금지 어법
    "전선을 넘었다", "GA 를 이긴다", "GA를 이긴다", "~9.5배", "독성 단조",
]
_NUMERIC = {"0/108", "5/5 지배", "6/6 재현", "3종 중 2종", "13/15", "0/3", "4/9"}

# 같은 줄(또는 바로 앞뒤 줄)에 이게 있으면 정정문·역사 서술로 보고 통과시킨다.
OK_MARK = ["🔴", "⚰️", "🚫", "§", "아니다", "없다", "말 것", "쓰지 마", "쓰지 않", "않는다",
           "폐기", "정정", "옛 ", "금지", "깨졌", "철회", "대체", "한때", "이전", "미반영"]


def _targets() -> list:
    """🔴 목록을 **손으로 적지 않는다**(§67-8 의 교훈). 디렉터리를 훑는다."""
    out = []
    for f in sorted(os.listdir(ROOT)):
        if f.endswith(".md"):
            out.append(f)
    d = os.path.join(ROOT, "docs")
    if os.path.isdir(d):
        out += [f"docs/{f}" for f in sorted(os.listdir(d)) if f.endswith(".md")]
    return out


def _real(line: str, bad: str, at: int) -> bool:
    """숫자 비(比)는 부분 문자열로 찾으면 안 된다 — `0/3` 이 `0/37` 에 걸렸다(§67-8)."""
    if bad not in _NUMERIC:
        return True
    tail = line[at + len(bad): at + len(bad) + 1]
    return not tail.isdigit()


def scan(rel: str) -> list:
    path = os.path.join(ROOT, rel)
    lines = open(path, encoding="utf-8").read().split("\n")
    hits = []
    for i, line in enumerate(lines):
        # 문맥 = 앞뒤 **3줄**까지. 🔴 ±1 줄로는 좁았다 — 인용 블록(`>`) 아래에 붙인
        # 정정 단락이 2~3줄 떨어져 있어 위양성이 났다(EVALUATION.md 의 "5/5 지배").
        ctx = "\n".join(lines[max(0, i - 3): i + 4])
        for bad in BAD:
            at = line.find(bad)
            if at < 0 or not _real(line, bad, at):
                continue
            if any(o in ctx for o in OK_MARK):
                continue
            hits.append((i + 1, bad, line.strip()[:170]))
    return hits


def main() -> int:
    print("=" * 78)
    print("§67-9 문서 주장 점검 — 폐기·반증된 주장이 **정정 없이** 남아 있나")
    print("=" * 78)

    tgt = _targets()
    cur = [f for f in tgt if f not in LOGS]
    logs = [f for f in tgt if f in LOGS]
    total = 0

    print(f"\n[A] 현재형 문서 {len(cur)}개 — 폐기 주장이 **있다는 것 자체**가 검토 대상")
    for rel in cur:
        h = scan(rel)
        total += len(h)
        if not h:
            continue
        print(f"\n  ── {rel}  ({len(h)}건)")
        for ln, bad, txt in h:
            print(f"     :{ln}  [{bad}]  {txt}")
    if total == 0:
        print("  ✅ 0건")

    print(f"\n[B] 실험 로그 {len(logs)}개 — 폐기 주장이 **기록으로** 남는 것이 정상")
    print("    (정정 표시 없이 뜬 줄만 센다. 적중이 0 일 필요는 없다.)")
    for rel in logs:
        h = scan(rel)
        print(f"  ── {rel}: {len(h)}건" + ("" if h else "  ✅"))
        for ln, bad, txt in h[:25]:
            print(f"     :{ln}  [{bad}]  {txt}")
        if len(h) > 25:
            print(f"     … 외 {len(h) - 25}건")

    print("\n" + "=" * 78)
    print("판정: 적중은 **자동 실패가 아니다** — '사람이 봐야 할 자리'를 좁혀 줄 뿐이다.")
    print("      🔴 `BAD` 는 늘려야 하는 목록이고, 코드 쪽 audit_prompt_strings.BAD 와")
    print("         **같이** 늘려야 한다.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
