"""설계 이력 원장(ledger) — 세션 간 누적 학습.

## 왜 필요한가

실측(docs/experiments.md §11.2): 같은 화물·같은 설정으로 4세션을 돌리면 유효점수가
**0.408~0.635** 로 흩어진다(평균 0.498 = 후보군 최적의 101%). 1세션이 134.8% 를 찍고
다른 세션이 83% 를 찍는다. 원인은 **매 실행이 콜드 스타트**라는 것이다. 모든 세션이
T7·Tet1 을 처음부터 다시 발견하고, 나쁜 궤적은 나쁜 수를 **반복**한다. 영속되는 건
캐시(`reachable`, `esmfold`)뿐이고 **설계 지식은 한 줄도 안 남았다**.

원장은 (조합, 점수, 가설, 결과)를 화물별로 남기고 다음 세션 시작에 요약해 넣는다.
효과는 두 가지다.

1. **세션 간 래칫** — 이미 아는 최고치가 바닥이 된다. 세션 내 래칫이 '한 세션 안에서
   퇴행 금지'였다면, 원장은 그걸 세션 경계 밖으로 확장한다.
2. **실패 재방문 금지** — "독성으로 탈락", "접합부에서 효능 붕괴", "조성 쌓기로 미승인"
   같은 결과를 알고 시작하므로 같은 벽에 다시 부딪히는 라운드를 아낀다.

## 다양성 붕괴를 어떻게 막나

승자만 보여주면 에이전트가 그걸 그대로 제출하고 끝낸다(탐색 0). 그래서 요약에
**실패와 그 이유를 같이** 싣고, 프롬프트에서 **"이미 아는 것을 다시 제출하는 것은
가치가 없다 — 이 바닥을 넘어서라"**로 못박는다. 원장은 목표가 아니라 **하한**이다.

## 벤치마크에서는 반드시 끈다

원장을 켠 세션은 콜드 스타트 세션과 **비교 불가**다(이전 세션의 탐색을 물려받았으므로
'18평가로 134.8%' 같은 효율 수치의 의미가 달라진다). 그래서 `benchmark_agent.py` 와
애블레이션은 `use_ledger=False` 로 고정하고, 원장 효과는 **별도 실험**으로 켬/끔을
나란히 재야 한다. 이걸 섞으면 §12 의 분모 불일치와 같은 종류의 거짓 비교가 된다.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

from .config import BASE_DIR

LEDGER_DIR = BASE_DIR / ".cache" / "ledger"

# 원장 1개가 무한정 커지지 않게 상한을 둔다. 넘으면 유효점수 하위부터 버리되
# **실패 기록은 따로 보존**한다(실패는 점수가 낮아서 먼저 잘리는데, 정작 재방문
# 방지에 쓸모가 있는 쪽이다).
MAX_WINS = 60
MAX_FAILS = 40


def _path(cargo: str) -> Path:
    """화물별 파일. 점수는 화물에 의존하므로 **화물이 다르면 다른 원장**이다."""
    h = hashlib.sha1((cargo or "").upper().encode()).hexdigest()[:12]
    return LEDGER_DIR / f"{h}.json"


def load(cargo: str) -> dict:
    """원장 읽기. 파일이 없거나 깨졌으면 **빈 원장**을 돌려준다(실행을 막지 않는다)."""
    p = _path(cargo)
    if not p.exists():
        return {"cargo": cargo, "sessions": 0, "wins": [], "fails": []}
    try:
        d = json.loads(p.read_text())
        if not isinstance(d, dict):
            raise ValueError("원장 형식 오류")
        d.setdefault("cargo", cargo)
        d.setdefault("sessions", 0)
        d.setdefault("wins", [])
        d.setdefault("fails", [])
        # 과거 버전이 비교군(CPP)이나 라이브러리 밖 링커를 승자로 기록했을 수 있다
        # → 읽을 때 걸러낸다(자가 치유). 둘 다 **도달 불가능한 바닥**을 만든다.
        d["wins"] = [e for e in d["wins"]
                     if not _is_cpp(e) and not _bad_linker(e) and not _length_hack(e)]
        return d
    except Exception:  # noqa: BLE001
        # 손상된 원장은 조용히 무시한다. 보조 장치가 본 실행을 죽여선 안 된다.
        return {"cargo": cargo, "sessions": 0, "wins": [], "fails": []}


def _key(e: dict) -> tuple[str, str]:
    return ((e.get("linker") or "").upper(), (e.get("shuttle") or "").upper())


def _bad_linker(e: dict) -> bool:
    """라이브러리 밖 링커인가 — 승자 목록에서 빼야 한다.

    실측(docs/experiments.md §13): 화물+T7 에 **아무 잔기나 하나** 끼우면 융합 BBB 가
    0.526 → 평균 0.743(20종 중 18종 상승)으로 뛴다. deepB3P 의 접합부 민감도이고, 발린
    하나가 TfR 결합을 바꿀 메커니즘은 없다. 원장이 이런 조합을 '승자'로 들고 있으면
    다음 세션에 **도달 불가능하고 무의미한 바닥**을 깔아준다(실제로 바닥이 0.918까지
    올라갔다 — 전원 `V+...`). 읽을 때 걸러내 자가 치유한다.
    """
    from .optimizer_agent import LIB_LINKERS
    lk = (e.get("linker") or "").upper()
    return bool(lk) and lk not in LIB_LINKERS


def _length_hack(e: dict) -> bool:
    """부모 셔틀에 길이만 덧붙인 조합인가 — 승자 목록에서 빼야 한다.

    실측(docs/experiments.md §13): 부모 T7 에 **무작위 1잔기**를 덧붙인 20종의 융합 BBB 가
    평균 0.743(부모 0.526), z≥2 는 0개였다. 즉 상승분은 '어느 잔기를 골랐나'가 아니라
    '길이를 늘렸다'에서 온다. 에이전트가 실제로 `YHAIYPRH`(eff 0.835)를 승자로 올렸고,
    그건 `링커 Y + T7` 과 **같은 분자**다.

    `iz`(접합부 길이 통제 z)가 기록돼 있으면 그 값으로 판정한다. 옛 기록엔 없으므로
    **부모가 붙은 채 z 가 없으면 보수적으로 버린다** — 통과했다는 증거가 없는 것이다.
    """
    from .order_specificity import Z_MIN, find_parent
    from .optimizer_agent import LIB_SHUTTLES
    found = find_parent((e.get("shuttle") or "").upper(), LIB_SHUTTLES)
    if found is None or (found[1] + found[2]) <= 0:
        return False
    z = e.get("iz")
    return z is None or z < Z_MIN


def _is_cpp(e: dict) -> bool:
    """비교군(CPP)인가 — 원장의 **승자 목록에서 반드시 빼야 한다**.

    실측으로 밟은 버그: 첫 end-to-end 검증에서 1위가 `직접융합+TAT`(eff 0.692)로
    올라왔다. 에이전트는 CPP 를 최종 선택할 수 없으므로(비특이 흡수 비교군) 이 바닥은
    **도달 불가능**하다 — "0.692 를 넘어라"고 시켜놓고 넘을 수 있는 길을 막아둔 셈이다.
    docs/experiments.md §12 의 분모 불일치와 **같은 종류의 오류**다. 기준점은 언제나
    에이전트가 실제로 선택할 수 있는 집합에서 나와야 한다.
    """
    return str(e.get("mech") or e.get("mechanism") or "").upper().startswith("CPP")


def record(cargo: str, rows: list[dict], final: dict | None = None,
           notes: dict | None = None) -> dict:
    """세션 결과를 원장에 합친다. 같은 (링커, 셔틀)은 **유효점수 높은 쪽만** 남긴다.

    `rows` 는 이번 세션이 채점한 행들(evaluate/design 모두), `final` 은 최종 선택,
    `notes` 는 (링커,셔틀)→사유 문자열(예: "순서 특이성 미달", "접합부 효능 붕괴").
    """
    d = load(cargo)
    notes = notes or {}
    wins = {_key(e): e for e in d["wins"]}
    fails = {_key(e): e for e in d["fails"]}

    for r in rows or []:
        eff = (r.get("bbb") or 0.0) - 0.25 * (r.get("sel_off") or 0.0)
        e = {"linker": r.get("linker") or "", "shuttle": r.get("shuttle") or "",
             "eff": round(eff, 4), "bbb": round(r.get("bbb") or 0.0, 4),
             "off": round(r.get("sel_off") or 0.0, 3),
             "inst": r.get("instability"), "toxic": bool(r.get("toxic")),
             "mech": r.get("mechanism") or r.get("sel_mech"),
             "delta": r.get("order_delta"), "iz": r.get("insert_z"),
             "ts": int(time.time())}
        k = _key(e)
        reason = notes.get(k) or notes.get((k[0], k[1]))
        if e["toxic"] or (e["bbb"] < 0.15) or reason:
            e["why"] = reason or ("독성 탈락" if e["toxic"] else "효능 하한 미달")
            prev = fails.get(k)
            if prev is None or (prev.get("eff") or 0) < e["eff"]:
                fails[k] = e
            wins.pop(k, None)
            continue
        if _is_cpp(e) or _bad_linker(e) or _length_hack(e):
            continue    # 비교군·라이브러리 밖 링커·길이 해킹은 선택 대상이 아니다
        if k in fails:
            continue    # 한 번 실패로 기록된 조합은 승자 목록으로 올리지 않는다
        prev = wins.get(k)
        if prev is None or (prev.get("eff") or 0) < e["eff"]:
            wins[k] = e

    d["wins"] = sorted(wins.values(), key=lambda e: -(e.get("eff") or 0))[:MAX_WINS]
    d["fails"] = sorted(fails.values(), key=lambda e: -(e.get("eff") or 0))[:MAX_FAILS]
    d["sessions"] = int(d.get("sessions") or 0) + 1
    if final is not None:
        d["last_final"] = {"linker": final.get("linker") or "",
                           "shuttle": final.get("shuttle") or "",
                           "eff": round((final.get("bbb") or 0.0)
                                        - 0.25 * (final.get("sel_off") or 0.0), 4)}
    try:
        LEDGER_DIR.mkdir(parents=True, exist_ok=True)
        _path(cargo).write_text(json.dumps(d, ensure_ascii=False, indent=1))
    except Exception:  # noqa: BLE001
        pass            # 기록 실패가 결과를 무효화하진 않는다
    return d


def floor_eff(cargo: str) -> float | None:
    """이 화물에서 지금까지 확인된 최고 유효점수 = 이번 세션이 넘어야 할 **바닥**."""
    d = load(cargo)
    return d["wins"][0]["eff"] if d["wins"] else None


def brief(cargo: str, k_win: int = 6, k_fail: int = 6) -> str:
    """세션 시작 프롬프트에 넣을 요약. 원장이 비면 빈 문자열."""
    d = load(cargo)
    if not d["wins"] and not d["fails"]:
        return ""
    out = [f"\n\n[설계 이력 — 이 화물로 이미 {d['sessions']}세션을 돌렸다]",
           "아래는 **지난 세션들이 실제로 채점해 본** 결과다(네 추측이 아니라 오라클 실측값).",
           "비교군(CPP)은 선택 대상이 아니므로 **목록에서 제외**했다 — 바닥은 네가 실제로 "
           "선택할 수 있는 후보군에서만 나온다."]
    if d["wins"]:
        out.append(f"\n확인된 상위 후보 (유효점수 = BBB − 0.25·off_target):")
        for e in d["wins"][:k_win]:
            dl = "" if e.get("delta") is None else f"·순서델타 {e['delta']:+.3f}"
            out.append(f"  - {e['linker'] or '직접융합'}+{e['shuttle']} → eff {e['eff']:.3f} "
                       f"(BBB {e['bbb']:.3f}·off {e['off']:.2f}{dl})")
        out.append(f"\n⚠️ **바닥은 {d['wins'][0]['eff']:.3f} 다.** 이걸 다시 찾아 제출하는 것은 "
                   "가치가 0이다 — 이미 아는 값이다. 네 할 일은 **이 바닥을 넘는 것**이고, "
                   "라이브러리 조합으로는 더 못 넘기면 `design_candidate`로 설계해서 넘어라. "
                   "넘지 못했으면 바닥 후보를 그대로 제출하되 '넘지 못했다'고 정직히 쓰라.")
    if d["fails"]:
        out.append("\n이미 실패로 확인된 것들 (같은 벽에 다시 부딪히지 마라):")
        for e in d["fails"][:k_fail]:
            out.append(f"  - {e['linker'] or '직접융합'}+{e['shuttle']} → {e.get('why', '탈락')}"
                       f" (BBB {e['bbb']:.3f})")
    out.append("\n이 이력은 **참고 자료**다. 점수를 우기는 근거로 쓰지 말고, 탐색의 "
               "출발점으로만 쓰라 — 최종 판정은 이번 세션의 실측과 독립 심사가 한다.")
    return "\n".join(out)
