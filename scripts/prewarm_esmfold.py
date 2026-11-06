#!/usr/bin/env python3
"""ESMFold 캐시 프리워밍 — **최종 선택 가능한 후보만** 미리 접어 둔다.

왜 (좁은) 프리워밍인가:
    ESMFold 공개 API(api.esmatlas.com)는 자주 502/503/504를 뱉거나 통째로 내려간다.
    다만 `fold_esmfold`에 이미 **재시도 3회·회당 25초·총 60초 상한·백오프**와 sha256 디스크
    캐시가 있고, 구조 축은 **보조 지표**라 실패해도 결과는 나온다. 따라서 프리워밍은
    "기능상 필수"가 아니라 **시연 중 구조 카드가 비지 않게 하는 보험**이다.

    그래서 전수(링커×셔틀 전 조합)를 접지 않는다. 대신 **하드 게이트(독성·BBB 하한)를 통과한
    후보군**만 접는다 — CPP(비교군)는 최종 선택 대상이 아니므로 제외된다. 보통 5~15개라
    1~2분이면 끝나고, 에이전트가 무엇을 고르든 캐시 히트가 된다.

    ⚠️ 2026-11 이후 역할이 더 줄었다: 애블레이션 결과 `analyze_structure`를 **에이전트 도구에서
    제거**했으므로(docs/experiments.md §11) 한 세션의 폴딩은 **최종 후보 1건**뿐이다. 즉 이
    스크립트는 이제 "시연 직전 1건을 미리 데워 두는" 용도이고, 안 돌려도 결과는 그대로 나온다.

무엇을 하나:
    `OptimizationAgent.reachable_front(cargo)`로 게이트 통과 후보군을 구해(로컬 계산·캐시)
    그 조합만 접는다. `--all`을 주면 예전처럼 전 조합을 접는다(권장하지 않음).

사용법:
    python scripts/prewarm_esmfold.py                 # DEFAULT_CARGO, 후보군만
    python scripts/prewarm_esmfold.py GSNKGAIIGLM     # 특정 화물
    python scripts/prewarm_esmfold.py GSNKGAIIGLM 1.5 # + 호출 간 간격(초, 서버 예의)
    python scripts/prewarm_esmfold.py GSNKGAIIGLM 1.0 --all   # 전 조합(구버전 동작)

서버가 죽어 있으면 실패한 조합은 건너뛰고 계속 진행한다(부분 성공도 캐시에 남는다).
서버가 회복된 뒤 다시 실행하면 남은 것만 채운다(캐시된 건 재요청 안 함).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

# scripts/ 에서 실행해도 프로젝트 루트의 core 패키지를 import 할 수 있도록.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.config import DEFAULT_CARGO, LINKER_LIBRARY, SHUTTLES  # noqa: E402
from core.structure import CACHE_DIR, fold_esmfold  # noqa: E402


def main() -> int:
    _args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cargo = (_args[0] if _args else DEFAULT_CARGO).upper()
    _pos = [a for a in sys.argv[1:] if not a.startswith("--")]
    gap = float(_pos[1]) if len(_pos) > 1 else 1.0  # 호출 간 간격(초)

    want_all = "--all" in sys.argv
    if want_all:
        combos = [(ln, lv["seq"], sn, sv["seq"])
                  for ln, lv in LINKER_LIBRARY.items()
                  for sn, sv in SHUTTLES.items()]
        scope = f"전 조합 {len(combos)}개(--all)"
    else:
        # 하드 게이트를 통과한 **후보군**만 — 최종 선택이 될 수 있는 조합이 곧 구조 검증 대상이다.
        from core.config import get_settings
        from core.optimizer_agent_gemini import get_gemini_agent
        from core.optimizer_agent import OptimizationAgent
        settings = get_settings()
        agent = get_gemini_agent(settings) or OptimizationAgent(settings)
        front = agent.reachable_front(cargo)
        seqs = {r["label"]: r["shuttle"] for r in (front.get("rmt_options") or [])}
        for r in (front.get("front") or []):
            seqs.setdefault(r["label"], r["shuttle"])
        combos = []
        for label, sseq in seqs.items():
            ln = label.split("+", 1)[0]
            lseq = LINKER_LIBRARY.get(ln, {}).get("seq", "")
            sn = label.split("+", 1)[1] if "+" in label else label
            combos.append((ln, lseq, sn, sseq))
        scope = (f"게이트 통과 후보군 {len(combos)}개"
                 f" (전수 {front.get('n_total', '?')}조합 중 · CPP 비교군 제외)")

    total = len(combos)
    if total == 0:
        print(f"화물 {cargo}: 게이트를 통과한 후보가 없어 프리워밍할 조합이 없습니다.")
        return 0
    print(f"화물 {cargo} ({len(cargo)}aa) · {scope} 프리워밍 → {CACHE_DIR}")

    ok = cached = fail = 0
    for i, (ln, lseq, sn, sseq) in enumerate(combos, 1):
        seq = cargo + lseq + sseq
        # 이미 캐시돼 있으면 fold_esmfold가 즉시 반환하므로 네트워크·간격을 낭비하지 않는다.
        import hashlib
        key = CACHE_DIR / (hashlib.sha256(
            "".join(c for c in seq.upper() if c.isalpha()).encode()).hexdigest() + ".pdb")
        was_cached = key.exists()

        pdb, err = fold_esmfold(seq)
        tag = f"[{i:>3}/{total}] {sn} · {ln} ({len(seq)}aa)"
        if err:
            fail += 1
            print(f"{tag} — 실패: {err.splitlines()[0][:70]}")
        elif was_cached:
            cached += 1
            print(f"{tag} — 캐시됨(건너뜀)")
        else:
            ok += 1
            print(f"{tag} — 접음·캐시 저장")
            if gap > 0 and i < total:
                time.sleep(gap)  # 새로 접은 경우에만 서버에 예의상 간격

    print(f"\n완료 — 새로 {ok} · 기존캐시 {cached} · 실패 {fail} / 총 {total}")
    if fail:
        print("실패분은 ESMFold 서버가 회복된 뒤 이 스크립트를 다시 실행하면 채워집니다.")
    return 1 if fail and not (ok or cached) else 0


if __name__ == "__main__":
    raise SystemExit(main())
