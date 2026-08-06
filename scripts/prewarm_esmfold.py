#!/usr/bin/env python3
"""ESMFold 캐시 프리워밍 — 데모 화물의 모든 링커×셔틀 조합을 미리 접어 둔다.

왜 필요한가:
    ESMFold 공개 API(api.esmatlas.com)는 자주 502/503/504를 뱉거나 통째로 내려간다.
    발표/시연 중에 구조 축이 실패하면 UX가 나빠진다(결과 자체는 보조 지표라 안 막히지만).
    구조 예측기는 sha256(전체 서열) 키로 `.cache/esmfold`에 결과를 캐싱하므로, **미리**
    유력 조합을 접어 캐시에 넣어두면 그 화물 데모는 네트워크 상태와 무관하게 항상 구조가 뜬다.

무엇을 하나:
    주어진 화물(기본: config.DEFAULT_CARGO)에 대해 라이브러리의 모든 링커×셔틀 조합
    (기본 13×6=78개)을 접어 캐시에 저장한다. 이미 캐시된 건 즉시 건너뛴다(fold_esmfold 내장).
    에이전트는 상위 1~3개만 구조 검증하지만, 전 조합을 데워두면 무엇을 고르든 캐시 히트다.

사용법:
    python scripts/prewarm_esmfold.py                # DEFAULT_CARGO
    python scripts/prewarm_esmfold.py GSNKGAIIGLM    # 특정 화물
    python scripts/prewarm_esmfold.py GSNKGAIIGLM 1.5 # + 호출 간 간격(초, 서버 예의)

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
    cargo = (sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CARGO).upper()
    gap = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0  # 호출 간 간격(초)

    combos = [(ln, lv["seq"], sn, sv["seq"])
              for ln, lv in LINKER_LIBRARY.items()
              for sn, sv in SHUTTLES.items()]
    total = len(combos)
    print(f"화물 {cargo} ({len(cargo)}aa) · 조합 {total}개 프리워밍 → {CACHE_DIR}")

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
