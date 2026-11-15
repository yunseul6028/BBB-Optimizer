#!/usr/bin/env python3
# AUDIT-EXEMPT: 이 파일은 금지어 목록 자체를 담는다 (점검 대상에서 제외)
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
sys.path.insert(0, os.path.join(ROOT, "scripts"))

# 로그 — 폐기 주장이 **기록으로서** 남아 있는 것이 정상인 파일
LOGS = {"docs/experiments.md", "docs/protocol.md"}

# 🔴 **폐기어 목록은 더 이상 여기 손으로 적지 않는다**(§67-11) — 코드 점검기와 **같은**
#    모듈을 쓴다. 전에는 두 파일에 따로 적혀 있어서 "같이 늘려야 한다"고 주석만 달아
#    뒀는데, 그 주석이 지켜진 적이 없다(`선택성 축` 이 코드 쪽에만 있었다).
from _forbidden import load as _load_forbidden        # noqa: E402

BAD, _NUMERIC, _SRC = _load_forbidden()

# 같은 줄(또는 바로 앞뒤 줄)에 이게 있으면 정정문·역사 서술로 보고 통과시킨다.
OK_MARK = ["🔴", "⚰️", "🚫", "§", "아니다", "없다", "말 것", "쓰지 마", "쓰지 않", "않는다",
           "폐기", "정정", "옛 ", "금지", "깨졌", "철회", "대체", "한때", "이전", "미반영"]


# 🔴 2026-11-16 **세 번째 정정** — 처음엔 `ROOT/*.md` + `docs/*.md` 만 봤다.
#    `os.listdir` 는 **한 겹만** 본다 → `docs/figs/*.md`(표 3개)와 `.html` 그림이
#    통째로 빠져 있었다. 손으로 적은 목록을 디렉터리 훑기로 바꿨더니, 이번엔
#    **훑는 깊이**가 손으로 정해져 있었던 것이다. 같은 병의 세 번째 얼굴이다.
#    → `os.walk` 로 **전부** 훑고, 확장자도 늘린다.
_EXT = (".md", ".html", ".txt")
_SKIP_DIRS = {".git", ".cache", "node_modules", "__pycache__", ".venv", "venv",
              "deepb3p", "toxinpred3", "List_CitedPapers"}


def _targets() -> list:
    """🔴 목록도 **깊이도** 손으로 적지 않는다(§67-8 · §67-10 의 교훈)."""
    out = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames
                       if d not in _SKIP_DIRS and not d.startswith(".")]
        for f in sorted(filenames):
            if f.endswith(_EXT):
                rel = os.path.relpath(os.path.join(dirpath, f), ROOT)
                out.append(rel)
    return sorted(out)


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

    print(f"\n[폐기어] {len(BAD)}종  출처 {_SRC}")
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
