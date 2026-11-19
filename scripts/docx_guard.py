#!/usr/bin/env python3
"""원고 docx 를 **사람과 내가 번갈아 고칠 때** 충돌을 막는 검사기.

사용자가 워드에서 직접 고치고, 내가 각주·그림·표를 붙인다. 두 쪽이 같은 파일을
건드리므로 다음 세 가지가 조용히 깨질 수 있다:

  ⑴ 워드가 파일을 연 채로 내가 쓰면 → 사용자의 다음 저장이 내 변경을 **덮는다**
  ⑵ 내가 예전에 풀어둔 사본으로 작업하면 → 사용자의 수정을 **잃는다**
  ⑶ 내가 덮어쓰면 되돌릴 수 없다

이 스크립트는 ⑴을 **검출**하고, ⑶에 대비해 자동 백업을 남기며, ⑵를 막기 위해
파일 지문을 기록해 "내가 마지막으로 본 것과 같은가"를 답한다.

  check    지금 써도 되는지 (워드 잠금 · 지문 대조)
  baseline 현재 상태를 기준으로 기록
  backup   타임스탬프 사본 보관
"""
from __future__ import annotations
import hashlib, json, os, re, subprocess, sys, zipfile
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STATE = os.path.join(ROOT, ".cache", "docx_guard.json")     # .cache 는 gitignore 대상
BACKUP = os.path.join(ROOT, ".cache", "docx_backups")
T = re.compile(r"<w:t(?: [^>]*)?>(.*?)</w:t>", re.S)


def fingerprint(path: str) -> dict:
    """내용 지문 — 바이트 해시는 워드가 저장만 해도 바뀌므로 본문도 같이 센다."""
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    with zipfile.ZipFile(path) as z:
        doc = z.read("word/document.xml").decode("utf-8")
        try:
            fn = z.read("word/footnotes.xml").decode("utf-8")
        except KeyError:
            fn = ""
    text = "".join(T.findall(doc))
    return {
        "sha256": h,
        "bytes": os.path.getsize(path),
        "paras": doc.count("<w:p ") + doc.count("<w:p>"),
        "footnote_refs": doc.count("<w:footnoteReference"),
        "footnotes": len(re.findall(r'<w:footnote w:id="\d+"', fn)),
        "drawings": doc.count("<w:drawing>"),
        "tables": doc.count("<w:tbl>"),
        "text_sha": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "text_len": len(text),
    }


def word_has_open(path: str) -> bool:
    """워드(또는 LibreOffice)가 이 파일을 연 채인가. `~$` 잠금 파일도 함께 본다."""
    d, n = os.path.split(os.path.abspath(path))
    # 워드의 잠금 파일명은 앞 2글자를 떼고 `~$` 를 붙인 형태다
    if os.path.exists(os.path.join(d, "~$" + n[2:])) or os.path.exists(os.path.join(d, "~$" + n)):
        return True
    try:
        out = subprocess.run(["lsof", "--", os.path.abspath(path)],
                             capture_output=True, text=True, timeout=20).stdout
    except Exception:
        return False
    return any(l for l in out.splitlines()[1:] if l.strip())


def load() -> dict:
    return json.load(open(STATE, encoding="utf-8")) if os.path.exists(STATE) else {}


def save(st: dict) -> None:
    os.makedirs(os.path.dirname(STATE), exist_ok=True)
    json.dump(st, open(STATE, "w", encoding="utf-8"), ensure_ascii=False, indent=1)


def do_backup(path: str) -> str:
    os.makedirs(BACKUP, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    dst = os.path.join(BACKUP, f"{stamp}__{os.path.basename(path)}")
    with open(path, "rb") as a, open(dst, "wb") as b:
        b.write(a.read())
    return dst


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__); return 2
    cmd, path = sys.argv[1], sys.argv[2]
    if not os.path.exists(path):
        print(f"🔴 파일이 없다: {path}"); return 1
    key = os.path.relpath(os.path.abspath(path), ROOT)
    st, fp = load(), fingerprint(path)

    if cmd == "baseline":
        st[key] = fp | {"at": datetime.now().isoformat(timespec="seconds")}
        save(st)
        print(f"✅ 기준 기록 — {key}")
        print(f"   본문 {fp['text_len']}자 · 각주 {fp['footnote_refs']} · "
              f"그림 {fp['drawings']} · 표 {fp['tables']}")
        return 0

    if cmd == "backup":
        print("✅ 백업:", do_backup(path)); return 0

    if cmd != "check":
        print(f"🔴 모르는 명령: {cmd}"); return 2

    bad = False
    if word_has_open(path):
        print("🔴 워드가 이 파일을 **열어 두고 있다** — 지금 쓰면 다음 저장에 덮인다.")
        print("   사용자가 닫을 때까지 쓰지 말 것.")
        bad = True
    else:
        print("✅ 워드 잠금 없음")

    base = st.get(key)
    if not base:
        print("⚠️  기준 기록이 없다 — 먼저 `baseline` 을 돌릴 것.")
        return 1 if bad else 0

    if fp["text_sha"] == base["text_sha"]:
        print(f"✅ 본문이 내가 마지막으로 본 것과 같다 (기준 {base['at']})")
    else:
        bad = True
        print(f"🔴 본문이 **바뀌었다** (기준 {base['at']})")
        print(f"   글자수 {base['text_len']} → {fp['text_len']} ({fp['text_len']-base['text_len']:+d})")
        print("   🔴 내 작업 사본을 버리고 **이 파일을 다시 읽어서** 시작할 것.")
    for k, label in (("footnote_refs", "각주"), ("drawings", "그림"), ("tables", "표")):
        if fp[k] != base[k]:
            print(f"   ⚠️  {label} {base[k]} → {fp[k]}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
