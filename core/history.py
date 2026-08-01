"""자율 설계 에이전트 실행 기록의 디스크 영구 저장.

세션 메모리(`st.session_state`)는 브라우저 새로고침·서버 재시작 시 사라진다. 그래서 실행
레코드(화물·이벤트·라운드 등)를 JSON으로 남겨, 앱을 다시 열었을 때 **지난 실행을 그대로
복원**할 수 있게 한다.

UI(app.py)는 `save_run(rec)` 한 번, `load_runs()` 한 번만 부르면 되고(각 1줄), 직렬화·정리·
경로 관리는 전부 여기서 담당한다. 저장 위치는 `.cache/agent_runs/`(gitignore됨).

레코드 포맷(app.py의 `_rec`와 동일):
    {"cargo": str, "events": [AgentEvent], "rounds": int, "brain": str, "label": str}
"""

from __future__ import annotations

import json
import time

from .config import BASE_DIR
from .optimizer_agent import AgentEvent

HISTORY_DIR = BASE_DIR / ".cache" / "agent_runs"
KEEP = 3  # 최근 몇 개 실행을 남길지 (app.py의 '최근 3개' 정책과 일치)


def _event_to_dict(ev: AgentEvent) -> dict:
    return {"kind": ev.kind, "text": ev.text, "data": ev.data}


def _event_from_dict(d: dict) -> AgentEvent:
    return AgentEvent(kind=d.get("kind", ""), text=d.get("text", ""),
                      data=d.get("data") or {})


def _prune() -> None:
    """파일명이 저장 시각(ms)이라 정렬하면 시간순 → 오래된 것부터 삭제, 최근 KEEP개만 유지."""
    files = sorted(HISTORY_DIR.glob("*.json"))
    for f in files[:-KEEP]:
        try:
            f.unlink()
        except OSError:
            pass


def save_run(rec: dict) -> None:
    """실행 레코드를 JSON으로 저장하고 최근 KEEP개만 남긴다.

    기록 저장 실패가 에이전트 실행 자체를 막지 않도록 모든 예외를 조용히 삼킨다.
    """
    try:
        HISTORY_DIR.mkdir(parents=True, exist_ok=True)
        payload = {
            "cargo": rec.get("cargo", ""),
            "rounds": rec.get("rounds", 0),
            "brain": rec.get("brain", ""),
            "label": rec.get("label", ""),
            "events": [_event_to_dict(e) for e in rec.get("events", [])],
        }
        # 파일명 = 저장 시각(ns) → 사전순 정렬이 곧 시간순 정렬(고정 19자리).
        # 같은 시각에 여러 번 저장돼도 겹치지 않도록 이미 있으면 뒤에 카운터를 붙인다.
        stamp = time.time_ns()
        path = HISTORY_DIR / f"{stamp:019d}.json"
        bump = 0
        while path.exists():
            bump += 1
            path = HISTORY_DIR / f"{stamp:019d}_{bump}.json"
        path.write_text(json.dumps(payload, ensure_ascii=False, default=str))
        _prune()
    except Exception:  # noqa: BLE001 - 기록은 부가기능, 실패해도 실행은 계속
        pass


def load_runs(limit: int = KEEP) -> list[dict]:
    """저장된 실행 레코드를 **오래된→최신** 순으로 반환한다(app.py의 append 순서와 동일).

    events는 AgentEvent 객체로 복원되어, 새로 실행한 기록과 똑같이 렌더된다.
    파일이 없거나 깨졌으면 조용히 건너뛴다.
    """
    if not HISTORY_DIR.exists():
        return []
    out: list[dict] = []
    for f in sorted(HISTORY_DIR.glob("*.json"))[-limit:]:
        try:
            d = json.loads(f.read_text())
            d["events"] = [_event_from_dict(e) for e in d.get("events", [])]
            out.append(d)
        except Exception:  # noqa: BLE001 - 깨진 기록은 무시
            continue
    return out
