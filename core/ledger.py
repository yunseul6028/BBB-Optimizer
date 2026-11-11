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

# 원장 스키마/목적함수 버전. **목적함수가 바뀌면 반드시 올린다** — 옛 기록은 폐지된
# 축(유효점수 = BBB − 0.25·off)으로 매긴 순위라서, 그대로 '바닥'으로 깔면 새 목적함수를
# 오염시킨다(reachable 캐시 키를 올리는 것과 같은 이유). 버전이 다르면 **빈 원장**으로 시작한다.
LEDGER_VERSION = 2   # v2: 문헌 근거 그룹 + 위험 축 Pareto (2026-11-10)


def _path(cargo: str) -> Path:
    """화물별 파일. 점수는 화물에 의존하므로 **화물이 다르면 다른 원장**이다."""
    h = hashlib.sha1((cargo or "").upper().encode()).hexdigest()[:12]
    return LEDGER_DIR / f"{h}.json"


def load(cargo: str) -> dict:
    """원장 읽기. 파일이 없거나 깨졌으면 **빈 원장**을 돌려준다(실행을 막지 않는다)."""
    p = _path(cargo)
    empty = {"cargo": cargo, "version": LEDGER_VERSION, "sessions": 0, "wins": [], "fails": []}
    if not p.exists():
        return empty
    try:
        d = json.loads(p.read_text())
        if not isinstance(d, dict):
            raise ValueError("원장 형식 오류")
        if int(d.get("version") or 1) != LEDGER_VERSION:
            # 옛 목적함수로 매긴 바닥은 새 목적함수에서 **의미가 없다**(BBB 기반 eff).
            # 조용히 변환하지 않고 버린다 — 잘못된 바닥은 없는 바닥보다 나쁘다.
            return empty
        d.setdefault("cargo", cargo)
        d.setdefault("sessions", 0)
        d.setdefault("wins", [])
        d.setdefault("fails", [])
        # 과거 버전이 비교군(CPP)·라이브러리 밖 링커·근거 미달을 승자로 기록했을 수 있다
        # → 읽을 때 걸러낸다(자가 치유). 모두 **도달 불가능한 바닥**을 만든다.
        d["wins"] = [e for e in d["wins"]
                     if not _is_cpp(e) and not _bad_linker(e) and not _unsupported(e)]
        return d
    except Exception:  # noqa: BLE001
        # 손상된 원장은 조용히 무시한다. 보조 장치가 본 실행을 죽여선 안 된다.
        return empty


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


def _unsupported(e: dict) -> bool:
    """문헌 근거가 없어 **승자 목록에 올릴 수 없는** 항목인가.

    이 함수가 옛 `_length_hack`(부모 셔틀에 길이만 덧붙인 조합)을 **승계한다**.
    옛 판정의 근거는 "길이를 늘리면 BBB 가 오른다"(실측 §13: 부모 T7 에 무작위 1잔기를
    붙인 20종의 융합 BBB 평균 0.743 vs 부모 0.526)였고, 그 BBB 를 순위에 쓰지 않게 된
    지금(M6·M8)은 **길이를 늘릴 동기 자체가 없다.** 남는 질문은 하나다 — 이 조합이
    문헌 근거를 보유했나(등급 A·B·C). 미등재·등급 D 는 바닥이 될 수 없다.

    ⚠️ 옛 기록은 `tier` 가 없으므로 여기서 걸러진다 — 의도된 동작이다(LEDGER_VERSION 참고).
    """
    from .evidence import group_rank
    return group_rank(e.get("tier")) > 1


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

    from .optimizer_agent import rank_key
    for r in rows or []:
        # ⚠️ 키 이름을 행(row)과 **같게** 유지한다 — `rank_key` 를 그대로 쓰기 위해서다.
        e = {"linker": r.get("linker") or "", "shuttle": r.get("shuttle") or "",
             "tier": r.get("tier"), "tier_parent": r.get("tier_parent"),
             "tier_sim": r.get("tier_sim"),
             "dev_liab": r.get("dev_liab"), "instability": r.get("instability"),
             "sel_off": round(r.get("sel_off") or 0.0, 3),
             "dev_agg_hot": r.get("dev_agg_hot"),
             "toxic": bool(r.get("toxic")),
             "mech": r.get("mechanism") or r.get("sel_mech"),
             # BBB 는 **보고치**로만 남긴다(순위·바닥 계산에 쓰지 않는다).
             "bbb_report": round(r.get("bbb") or 0.0, 4),
             "ts": int(time.time())}
        k = _key(e)
        reason = notes.get(k) or notes.get((k[0], k[1]))
        # ⚠️ 옛 구현에는 `e["bbb"] < 0.15` 가 실패 조건으로 남아 있었다 — 폐지한
        # MIN_BBB_SIGNAL 을 원장이 **되살려 쓰고 있었다**(등급 A 셔틀을 전멸시키던 그 값이다).
        # 제거했다. 실패 사유는 독성과 게이트 사유(notes)뿐이다.
        if e["toxic"] or reason:
            e["why"] = reason or "독성 탈락"
            prev = fails.get(k)
            if prev is None or rank_key(e) < rank_key(prev):
                fails[k] = e
            wins.pop(k, None)
            continue
        if _is_cpp(e) or _bad_linker(e) or _unsupported(e):
            continue    # 비교군·라이브러리 밖 링커·근거 미달은 선택 대상이 아니다
        if k in fails:
            continue    # 한 번 실패로 기록된 조합은 승자 목록으로 올리지 않는다
        prev = wins.get(k)
        if prev is None or rank_key(e) < rank_key(prev):
            wins[k] = e

    d["version"] = LEDGER_VERSION
    d["wins"] = sorted(wins.values(), key=rank_key)[:MAX_WINS]
    d["fails"] = sorted(fails.values(), key=rank_key)[:MAX_FAILS]
    d["sessions"] = int(d.get("sessions") or 0) + 1
    if final is not None:
        d["last_final"] = {"linker": final.get("linker") or "",
                           "shuttle": final.get("shuttle") or "",
                           "tier": final.get("tier"),
                           "dev_liab": final.get("dev_liab"),
                           "sel_off": final.get("sel_off")}
    try:
        LEDGER_DIR.mkdir(parents=True, exist_ok=True)
        _path(cargo).write_text(json.dumps(d, ensure_ascii=False, indent=1))
    except Exception:  # noqa: BLE001
        pass            # 기록 실패가 결과를 무효화하진 않는다
    return d


def failure_notes(rows: list[dict]) -> dict:
    """채점한 행들에서 **(링커,셔틀) → 실패 사유** 표를 만든다.

    `record()` 의 `notes` 인자로 그대로 넣는다. 라이브러리 셔틀은 제외한다 — 사유를
    남길 대상은 에이전트가 **새로 만든** 후보다(라이브러리 조합은 전선 계산이 이미 판정한다).

    세션 종료 경로(`optimizer_agent_gemini`)와 원장 시딩(`benchmark_agent --seed-ledger`)이
    **같은 함수**를 쓰도록 여기 둔다. 두 곳에 복제하면 사유 문구가 갈라지고, 그러면
    원장이 세션마다 다른 기준으로 실패를 기록한다(§19-4 의 상수 복제와 같은 계열의 함정).
    """
    from .evidence import INHERIT_SIM_MIN
    from .optimizer_agent import LIB_SHUTTLES, is_unsupported

    notes: dict = {}
    for r in rows or []:
        sh = (r.get("shuttle") or "").upper()
        if not sh or sh in LIB_SHUTTLES:
            continue
        k = ((r.get("linker") or "").upper(), sh)
        if r.get("tier") is None:
            par, sim = r.get("tier_parent"), r.get("tier_sim") or 0.0
            notes[k] = (f"문헌 근거 미등재 — 최근접 {par or '없음'}·유사도 "
                        f"{sim:.2f} < {INHERIT_SIM_MIN:.2f} (후보 자격 없음)")
        elif is_unsupported(r):
            notes[k] = f"근거 등급 {r.get('tier')} = 근거 없음 (후보 자격 없음)"
    return notes


def floor_entry(cargo: str) -> dict | None:
    """이 화물에서 지금까지 확인된 **최선 후보** = 이번 세션이 넘어야 할 바닥.

    스칼라가 아니라 항목을 돌려준다 — 새 목적함수는 단일 점수가 아니라
    (근거 그룹 → 위험 축 Pareto)이므로, 넘었는지는 `optimizer_agent.dominates` 로 본다.
    """
    d = load(cargo)
    return d["wins"][0] if d["wins"] else None


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
        out.append("\n확인된 상위 후보 (근거 등급 → liability·불안정성·off·응집 순):")
        for e in d["wins"][:k_win]:
            inh = ("" if (e.get("tier_sim") or 1.0) >= 1.0 or not e.get("tier_parent")
                   else f"←{e['tier_parent']}·{e['tier_sim']:.2f}")
            out.append(f"  - {e['linker'] or '직접융합'}+{e['shuttle']} → 근거 "
                       f"**{e.get('tier')}**{inh} · liab {e.get('dev_liab')} · 불안정성 "
                       f"{e.get('instability')} · off {e.get('sel_off')} · 응집 "
                       f"{e.get('dev_agg_hot')} (BBB 참고 {e.get('bbb_report')})")
        t = d["wins"][0]
        out.append(f"\n⚠️ **바닥은 `{t['linker'] or '직접융합'}+{t['shuttle']}`"
                   f"(근거 {t.get('tier')}·liab {t.get('dev_liab')}·off {t.get('sel_off')}) 다.** "
                   "이걸 다시 찾아 제출하는 것은 가치가 0이다 — 이미 아는 값이다. 네 할 일은 "
                   "**이 후보에게 지배당하지 않는 것**이고(근거 등급이 같거나 높고 위험 축에서 "
                   "하나라도 더 나아야 한다), 라이브러리 조합으로 더 못 넘기면 "
                   "`design_candidate`·`optimize_developability`로 **부모의 근거를 유지한 채 "
                   "liability 를 줄여** 넘어라. 넘지 못했으면 바닥 후보를 그대로 제출하되 "
                   "'넘지 못했다'고 정직히 쓰라.\n"
                   "⚠️ **BBB 참고값을 올리는 방향으로 움직이지 마라** — 그 점수는 순위에 쓰지 "
                   "않는다(M6·M8).")
    if d["fails"]:
        out.append("\n이미 실패로 확인된 것들 (같은 벽에 다시 부딪히지 마라):")
        for e in d["fails"][:k_fail]:
            out.append(f"  - {e['linker'] or '직접융합'}+{e['shuttle']} → {e.get('why', '탈락')}")
    out.append("\n이 이력은 **참고 자료**다. 점수를 우기는 근거로 쓰지 말고, 탐색의 "
               "출발점으로만 쓰라 — 최종 판정은 이번 세션의 실측과 독립 심사가 한다.")
    return "\n".join(out)
