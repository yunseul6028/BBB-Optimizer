#!/usr/bin/env python3
# AUDIT-EXEMPT: 이 파일은 금지어 목록을 다룬다 (문자열 점검 대상에서 제외)
"""수치원장 **Ⅳ(🚫 쓰면 안 되는 숫자) ↔ Ⅴ(절별 인용 원장)** 교차 점검 (§73-3).

왜 이 스크립트가 존재하나
--------------------------------------------------------------------------
§73 에서 **Ⅴ 절별 인용 원장만** 낡아 있었다. `0/108`(§53-5 폐기) · `12 : 0`(§49-2
폐기) 이 Ⅳ 에서는 🚫 인데 Ⅴ 에서는 **그냥 인용 가능한 행**으로 남아 있었다.
Ⅴ 는 "논문 n절에 어떤 숫자를 쓸지" 를 적은 표 — **글 쓸 때 실제로 펴 놓고 보는 표**다.
그대로 쓰면 폐기한 숫자가 원고에 되돌아온다.

🔴 그런데 Ⅴ 를 봐 주는 자동 점검이 **하나도 없었다.** `audit_docs_claims.py` 는
수치원장을 일부러 제외한다 — 폐기 수치를 **등재하는** 문서라 전수 적중하기 때문이다.
맞는 결정이지만, 그래서 Ⅴ 가 사각지대였다. → **같은 문서 안에서 Ⅳ 와 Ⅴ 를 맞대 본다.**

이 점검이 §38-1 의 교훈을 따르는 방식: 경고 주석은 재발을 막지 못했다(§38-1 · §68-4 ·
§69-3 · §72-2 · §73 — 다섯 번). **묻는 쪽을 코드로 만든다.**

무엇을 보나
--------------------------------------------------------------------------
⑴ 🔴 **문구** — Ⅳ 1열의 금지 문구가 Ⅴ 의 한 행에 **폐기 표시 없이** 나타나면 실패.
⑵ 🔴 **수치** — Ⅳ 1열에서 꺼낸 수(비 `12 : 0` · 비율 `98%` · 소수 `0.96`)가 Ⅴ 의
   한 행에 **폐기 표시 없이** 나타나면 실패. §73 의 네 칸 중 셋이 여기서 걸린다.
⑶ 🟡 **대체 수치의 집** — Ⅳ 가 "대신 이걸 써라"고 제시한 수가 Ⅴ 에 **등재돼 있나**.
   실패가 아니라 안내다. 🔴 §73 의 진짜 구멍이 여기였다 — `4/108` 과 `57/96` 은
   §53-5 에서 측정되고 §72 에서 재확인됐는데 **Ⅴ 에 자리가 없었다.** 폐기만 표시하고
   대체를 등재하지 않으면, 글 쓸 때 펼 표에는 **쓸 수 있는 수가 비어 있다.**

⑷ 🔴 **정의 세대** — Ⅴ 에서 `liab`(liability)를 인용하는 행은 ①/②/③ 중
   어느 정의인지 밝혀야 한다(Ⅲ-2 · 정의가 **두 번** 바뀌었다).

판정 단위는 **행이 아니라 칸**이다. 🔴 처음엔 행으로 했는데 **네 칸을 놓쳤다** —
`| 0/108 vs 41/96 | … | 41/96 은 재현 불가 |` 처럼 **다른 수치**에 붙은 정정 표시가
행 전체를 통과시켰다. §67-11 과 같은 모양이다(규칙을 넓게 잡아 조용히 샜다).

재현: /usr/bin/python3 scripts/audit_ledger_crossref.py
"""
from __future__ import annotations

import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 🔴 경로를 인자로 받는다 — **고장 검사(negative test)를 할 수 있어야** 한다.
#    §49 의 교훈: 아무것도 못 잡는 하네스는 조용히 ✅ 를 찍는다. 정정 **이전**
#    판본에 돌려서 🔴 가 나오는지 확인하는 것이 이 스크립트의 자기 증명이다:
#      git show HEAD~1:docs/수치원장.md > /tmp/old.md
#      /usr/bin/python3 scripts/audit_ledger_crossref.py /tmp/old.md   # → 🔴 가 나와야 한다
LEDGER = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.join(ROOT, "docs", "수치원장.md")

import _forbidden                                     # noqa: E402
# 🔴 **한 파일만 본다** — 금지 문구도 검사 대상 파일의 Ⅳ 에서 꺼낸다.
#    (현행 Ⅳ + 옛 Ⅴ 를 섞으면 "어느 판본을 점검한 건가"가 흐려진다.)
_forbidden.LEDGER = LEDGER
_iv_phrases = _forbidden.from_ledger

# 같은 **행** 에 이 표시가 있으면 폐기를 알고 적은 것이므로 통과.
OK_MARK = (
    "🚫", "⚰️", "폐기", "금지", "쓰지 말", "말 것", "아니다", "않는다", "없다",
    "재현 불가", "재현되지", "미기록", "시드 의존", "정정", "옛 정의", "옛 ", "보류",
    "충돌", "미측정", "조건", "밝힐 것", "확인 필수", "대체",
)

# 🔴 7절은 **설계상 면제**다 — 제목이 "측정 도구가 틀렸던 사건들"이고 모든 행이
#    `옛값 → 새값` 정정 화살표다. 즉 그 절에서 폐기 수치의 등장은 **내용 자체**다.
#    (면제를 절 제목으로 선언한다 — 줄 번호로 적으면 문서가 자랄 때 어긋난다.)
EXEMPT_SUBSECTION = "7절"

_NUM_RE = re.compile(r"\d+\s*:\s*\d+|\d+\s*/\s*\d+|\d+(?:\.\d+)?\s*%|0\.\d{2,3}")


def _sections() -> tuple:
    """`(Ⅳ 행 목록, Ⅴ 행 목록[(절, 행)])`."""
    lines = io.open(LEDGER, encoding="utf-8").read().split("\n")

    def span(head: str) -> tuple:
        i = next(k for k, l in enumerate(lines) if l.startswith(head))
        j = next((k for k, l in enumerate(lines[i + 1:], i + 1) if l.startswith("## ")), len(lines))
        return i, j

    i4, j4 = span("## Ⅳ")
    i5, j5 = span("## Ⅴ")

    def rows(a: int, b: int) -> list:
        out, sub = [], ""
        for l in lines[a:b]:
            if l.startswith("### "):
                sub = l[4:].split("—")[0].strip()
            if not l.startswith("|") or re.match(r"^\|[\s\-|]+\|$", l):
                continue
            if l.split("|")[1].strip() in ("숫자", "공백", "자리"):
                continue
            out.append((sub, l))
        return out

    return rows(i4, j4), rows(i5, j5)


def _iv_numbers(iv_rows: list) -> list:
    """Ⅳ 1열(= 금지되는 양)에서만 수를 꺼낸다. 2열·3열에는 **대체 수치**가 섞여 있다."""
    out, seen = [], set()
    for _, line in iv_rows:
        for m in _NUM_RE.findall(line.split("|")[1]):
            t = re.sub(r"\s+", " ", m).strip()
            if t not in seen:
                seen.add(t)
                out.append(t)
    return out


def _iv_replacements(iv_rows: list) -> list:
    """Ⅳ 3열(대신)에서 꺼낸 **쓸 수 있는** 수 — Ⅴ 에 집이 있나 보려고."""
    out, seen = [], set()
    for _, line in iv_rows:
        cols = line.split("|")
        if len(cols) < 4:
            continue
        for m in _NUM_RE.findall(cols[3]):
            t = re.sub(r"\s+", " ", m).strip()
            if t in seen:
                continue
            seen.add(t)
            out.append(t)
    return out


def _hit(row: str, term: str) -> bool:
    """행에 term 이 있나. 숫자 비는 **경계**를 본다(`0/3` 이 `0/37` 에 걸린 사고, §67-8)."""
    if not re.fullmatch(r"\d+\s*[:/]\s*\d+", term):
        return term in row
    pat = re.escape(term).replace(r"\ ", r"\s*")
    return re.search(r"(?<!\d)" + pat + r"(?!\d)", row) is not None


def main() -> int:
    iv, v = _sections()
    phrases = [p for p in _iv_phrases() if len(p) >= 4]
    numbers = _iv_numbers(iv)

    print("=" * 78)
    print("수치원장 Ⅳ ↔ Ⅴ 교차 점검 (§73-3)")
    print("=" * 78)
    print(f"Ⅳ 행 {len(iv)}개 → 금지 문구 {len(phrases)}종 · 금지 수치 {len(numbers)}종")
    print(f"Ⅴ 행 {len(v)}개  (면제: {EXEMPT_SUBSECTION} — 제목이 '측정 도구가 틀렸던 사건들')")

    fails = []
    for kind, terms in (("문구", phrases), ("수치", numbers)):
        print(f"\n[{'⑴' if kind == '문구' else '⑵'} 금지 {kind}가 Ⅴ 에 표시 없이 있나]")
        bad = []
        for sub, row in v:
            if sub == EXEMPT_SUBSECTION:
                continue
            for cell in row.split("|"):
                if any(m in cell for m in OK_MARK):
                    continue                     # 이 **칸**이 스스로 폐기를 밝혔다
                for t in terms:
                    if _hit(cell, t):
                        bad.append((sub, t, row))
        if bad:
            for sub, t, row in bad:
                print(f"  🔴 {sub} · [{t}]")
                print(f"     {row.strip()[:160]}")
            fails.extend(bad)
        else:
            print(f"  ✅ 0건")

    print("\n[⑶ 🟡 Ⅳ 가 '대신 쓰라'고 한 수치가 Ⅴ 에 등재됐나]")
    vtext = "\n".join(r for _, r in v)
    missing = [t for t in _iv_replacements(iv) if not _hit(vtext, t)]
    if missing:
        print(f"  🟡 Ⅴ 에 없는 대체 수치 {len(missing)}종 — 실패가 아니다(절에 넣을 자리가")
        print(f"     없을 수도 있다). 🔴 다만 **폐기만 표시하고 대체를 등재하지 않은 상태**가")
        print(f"     §73 의 구멍이었다. 하나씩 '쓸 자리가 없나'를 눈으로 가른다:")
        for t in missing:
            print(f"       · {t}")
    else:
        print("  ✅ 전부 Ⅴ 에 있다")

    print("\n[⑷ `liab` 인용 행이 정의 세대(①/②/③)를 밝혔나]")
    gen_bad = []
    for sub, row in v:
        if sub == EXEMPT_SUBSECTION or "liab" not in row.lower():
            continue
        if not any(m in row for m in ("①", "②", "③", "정의")):
            gen_bad.append((sub, "정의 세대 미표기", row))
    if gen_bad:
        for sub, t, row in gen_bad:
            print(f"  🔴 {sub} · {t}")
            print(f"     {row.strip()[:160]}")
        fails.extend(gen_bad)
    else:
        print("  ✅ 0건")

    print("\n" + "=" * 78)
    if fails:
        print(f"판정: 🔴 **{len(fails)}건** — Ⅴ 는 글 쓸 때 펴는 표다. 폐기 표시 없는 행은")
        print("      그대로 원고에 들어간다(§73). 해당 행에 🚫·폐기 사유·대체 수치를 적는다.")
    else:
        print("판정: ✅ Ⅳ 의 금지 수치가 Ⅴ 에 표시 없이 남아 있지 않다.")
        print("      🔴 ⑶ 의 🟡 는 자동으로 못 가른다 — 눈으로 본다.")
    print("=" * 78)
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
