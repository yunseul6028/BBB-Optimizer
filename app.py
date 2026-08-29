"""
BBB-Optimize AI Agent — Streamlit UI (thin layer)
=================================================
화물(cargo) 서열을 받아, **자율 설계 에이전트**가 링커·셔틀을 조립·잔기 편집·서열 진화 도구로
오케스트레이션하고 deepB3P(BBB)·ToxinPred3(독성)·선택성 등 8축으로 평가해 최종 융합체를 설계한다.
UI는 얇은 표시·스트리밍 계층이며 계산은 core/*.py.

실행:
    pip install -r requirements.txt
    streamlit run app.py
"""

from pathlib import Path

import altair as alt
import pandas as pd
import streamlit as st

from core import get_settings
from core.config import (
    DEFAULT_CARGO,
    LINKER_LIBRARY,
    MODEL_MAX_LEN,
    SHUTTLES,
    VALID_AMINO_ACIDS,
)
from core.optimizer_agent_gemini import get_gemini_agent, OFF_TARGET_PENALTY
from core.schemas import Verdict

st.set_page_config(page_title="BBB-Optimize AI Agent", page_icon="🧬", layout="wide")


def _inject_theme():
    css = Path(__file__).parent / "assets" / "theme.css"
    if css.exists():
        st.html(f"<style>{css.read_text(encoding='utf-8')}</style>")


_inject_theme()
settings = get_settings()

st.title("BBB-Optimize AI Agent")
st.markdown("화물 펩타이드를 입력하면 AI가 BBB를 통과하는 **항체–셔틀 융합 단백질**을 설계합니다.")


def _dot(ok):
    return "연동" if ok else "대기"


def _clean_cargo(raw):
    """화물 입력 정리: FASTA 헤더(>) 줄 제거 + 모든 공백·줄바꿈·탭 제거 + 대문자화.
    (복붙 시 딸려오는 줄바꿈·띄어쓰기를 자동으로 걷어낸다.)"""
    lines = [ln for ln in (raw or "").splitlines() if not ln.strip().startswith(">")]
    return "".join("".join(lines).split()).upper()


def _cargo_error(cargo):
    """화물 서열이 유효하면 None, 아니면 어떤 문자가 문제인지 구체적 메시지를 돌려준다."""
    if not cargo:
        return "⚠️ 화물 펩타이드 서열을 입력해 주세요."
    bad = sorted(set(cargo) - VALID_AMINO_ACIDS)
    if bad:
        shown = ", ".join(("공백/줄바꿈" if not c.strip() else f"`{c}`") for c in bad)
        return (f"⚠️ 표준 20종 아미노산이 아닌 문자가 있습니다: {shown}. "
                "1글자 코드(ACDEFGHIKLMNPQRSTVWY)만 사용하세요 — 항체 서열의 **X(불명 잔기)**·"
                "공백·복붙 특수문자 등을 확인하세요.")
    return None


st.caption(
    f"엔진 — deepB3P {_dot(settings.use_deepb3p_local)} · "
    f"ToxinPred3 {_dot(settings.use_toxinpred3_local)} · ESMFold 연동(API) · "
    f"Gemini {_dot(settings.use_gemini_agent)}"
    f"{'' if settings.use_gemini_agent else ' · 키 필요'}"
)
with st.expander("동작 방식 · 지표 안내"):
    st.markdown(
        "- **설계 에이전트**가 평가·잔기 편집·구조 도구로 후보를 만들고, **심사 에이전트**가 적대적으로 검증합니다.\n"
        "- **BBB 점수**는 deepB3P 예측 확률(0–100)로 후보 비교용이며, 실제 투과율은 실험(Papp/logBB)이 필요합니다.\n"
        "- `융합체 = 화물 + 링커 + 셔틀`"
    )

n_combos = len(LINKER_LIBRARY) * len(SHUTTLES)

st.divider()

# --- 입력 (중앙 정렬 히어로 · 라이브러리 전체 자동 탐색) --------------------
_hl, _hero, _hr = st.columns([1, 2, 1])
with _hero:
    cargo_input = st.text_input(
        "화물(cargo) 서열",
        value=DEFAULT_CARGO,
        help="뇌로 전달할 화물 서열을 넣으면 자율 설계 에이전트가 링커·셔틀을 붙여 BBB 투과 융합체를 "
             "설계합니다. 항체·나노바디도 화물로 처리됩니다. (표준 20종 1글자 코드)",
    )
    _seq = _clean_cargo(cargo_input)

    # 실행 플래그·기본값
    agent_run = False
    agent_rounds = 10

    # === 자율 설계 에이전트 (입력은 항상 화물) ===
    if settings.use_gemini_agent:
        agent_rounds = st.slider("최대 스텝 수", 4, 12, 10,
                                 help="에이전트가 도구를 호출하는 최대 횟수. 낮추면 빠르고 저렴하지만, "
                                      "너무 낮으면 finish 전에 예산이 소진돼 심사가 축약될 수 있습니다.")
        st.caption(f"모델 — Gemini ({settings.gemini_model})")
    else:
        st.caption("데모 모드 — API 키 없이 예시 결과로 화면을 미리봅니다.")
    agent_run = st.button("자율 설계 에이전트 실행", type="primary", width='stretch')
    st.caption("에이전트가 BBB·독성·구조·수용체·잔기 편집·서열 진화 도구를 스스로 오케스트레이션해 "
               "최종 융합체를 설계합니다.")
    with st.expander(f"라이브러리 구성 (링커 {len(LINKER_LIBRARY)} · 셔틀 {len(SHUTTLES)})"):
        st.markdown("**링커**")
        st.dataframe(
            {"링커": list(LINKER_LIBRARY), "서열": [v["seq"] for v in LINKER_LIBRARY.values()],
             "종류": [v["kind"] for v in LINKER_LIBRARY.values()],
             "설명": [v["note"] for v in LINKER_LIBRARY.values()]},
            width='stretch', hide_index=True,
        )
        st.markdown("**셔틀**")
        st.dataframe(
            {"셔틀": list(SHUTTLES), "서열": [v["seq"] for v in SHUTTLES.values()],
             "타겟": [v["target"] for v in SHUTTLES.values()],
             "설명": [v["note"] for v in SHUTTLES.values()]},
            width='stretch', hide_index=True,
        )

    _cargo_preview = _seq
    if _cargo_preview:
        longest = max(len(v["seq"]) for v in SHUTTLES.values())
        longest_l = max(len(v["seq"]) for v in LINKER_LIBRARY.values())
        est = len(_cargo_preview) + longest_l + longest
        if est > MODEL_MAX_LEN:
            st.info(
                f"일부 조합은 {MODEL_MAX_LEN}aa를 넘습니다(≈{est}aa). 이런 조합은 BBB를 "
                f"연결부위(링커+셔틀)만으로 계산해(화물 벌크 제외) 셔틀 신호가 희석되지 않습니다. "
                f"독성은 전체 서열로 계산합니다."
            )

# 상태 색 언어(theme.css와 동일): 초록=좋음 · 주황=애매 · 빨강=위험 · 보라=더 알아야함
_VERDICT_PILL = {
    Verdict.ACCEPTED:       ("채택(베스트)",   "#e7f6ed", "#0b6b34"),  # 초록
    Verdict.SUBOPTIMAL:     ("후순위",          "#fdf1de", "#8a5406"),  # 주황
    Verdict.REJECTED_TOXIC: ("독성 탈락",       "#fbe9e9", "#a32020"),  # 빨강
    Verdict.REFERENCE:      ("화물 단독(기준)", "#efeafc", "#4a2fa0"),  # 보라
}
# 판정 컬럼 셀 배경 (라벨 텍스트 → (배경, 글자))
_VERDICT_CELL = {label: (bg, fg) for label, bg, fg in _VERDICT_PILL.values()}
_PASS_CELL = {"통과": ("#e7f6ed", "#0b6b34"), "독성": ("#fbe9e9", "#a32020")}


def _pill(text: str, bg: str, fg: str) -> str:
    """색 알약(pill) 배지 HTML. st.markdown(..., unsafe_allow_html=True)로 렌더."""
    return (f'<span style="display:inline-block;font-size:12px;font-weight:600;'
            f'padding:3px 10px;border-radius:9999px;background:{bg};color:{fg};'
            f'vertical-align:middle;">{text}</span>')


def _verdict_pill(verdict) -> str:
    text, bg, fg = _VERDICT_PILL.get(verdict, ("-", "#eee", "#333"))
    return _pill(text, bg, fg)


def _approve_pill(approve: bool) -> str:
    return (_pill("승인 · APPROVE", "#e7f6ed", "#0b6b34") if approve
            else _pill("개선 요구 · REVISE", "#fbe9e9", "#a32020"))


def _style_verdict(data: dict, cellmap: dict, col: str = "판정"):
    """dict → 판정 컬럼을 상태색으로 배경 칠한 pandas Styler."""
    df = pd.DataFrame(data)

    def _c(v):
        pair = cellmap.get(v)
        return f"background-color:{pair[0]};color:{pair[1]};font-weight:600;" if pair else ""

    return df.style.map(_c, subset=[col])


class _Ev:
    """실제 에이전트 이벤트와 동일한 형태(kind/text/data)의 경량 이벤트 객체 (데모용)."""
    __slots__ = ("kind", "text", "data")

    def __init__(self, kind, text="", data=None):
        self.kind = kind
        self.text = text
        self.data = data or {}


def _demo_agent_events(cargo):
    """API 키 없이 결과 화면을 미리보기 위한 예시(더미) 이벤트 스트림.
    실제 Gemini 에이전트가 내보내는 이벤트와 같은 kind/text/data 형태로 흉내낸다."""
    cargo = cargo or "GHRPYD"
    sh, lk = list(SHUTTLES), list(LINKER_LIBRARY)

    def _mk(shuttle_name, linker_name, bbb, tox, **extra):
        s, l = SHUTTLES[shuttle_name], LINKER_LIBRARY[linker_name]
        row = {
            "label": f"{shuttle_name} · {linker_name}",
            "linker": l["seq"], "shuttle": s["seq"],
            "shuttle_name": shuttle_name, "linker_name": linker_name,
            "sequence": cargo + l["seq"] + s["seq"],
            "bbb": bbb, "tox": tox, "toxic": tox > settings.toxicity_threshold,
            "bind_ref": s.get("target", "TfR"), "bind_score": extra.get("bind_score", 0.72),
            "instability": extra.get("instability", 34), "stable": extra.get("stable", True),
            "dev_risk": extra.get("dev_risk", "낮음"), "dev_liab": extra.get("dev_liab", 1),
            "dev_charge": extra.get("dev_charge", "+2"),
            "dev_liabilities": extra.get("dev_liabilities", []),
            "sol_level": extra.get("sol_level", "양호"), "sol_score": extra.get("sol_score", 0.71),
            "selectivity": extra.get("selectivity", "높음"),
            "sel_level": extra.get("sel_level", "낮음"),
            "sel_mech": extra.get("sel_mech", f"{s.get('target', 'TfR')} 특이 결합"),
        }
        return row

    r_win = _mk(sh[0], lk[0], 0.94, 0.06, bind_score=0.81, instability=31)
    r2 = _mk(sh[1 % len(sh)], lk[1 % len(lk)], 0.88, 0.12, bind_score=0.69,
             instability=38, dev_risk="중간", selectivity="중간", sel_level="중간")
    r3 = _mk(sh[2 % len(sh)], lk[2 % len(lk)], 0.83, 0.09, bind_score=0.60, instability=41,
             dev_risk="중간")
    r_tox = _mk(sh[3 % len(sh)], lk[0], 0.71, 0.63, bind_score=0.55, stable=False,
                dev_risk="높음", dev_liab=3, sel_level="높음", selectivity="낮음")
    rows = [r_win, r2, r3, r_tox]

    yield _Ev("plan", text=(
        "**전략** — 라이브러리의 셔틀·링커 조합을 전수 평가해 BBB 투과·독성·안정성·"
        "수용체 결합을 함께 만족하는 융합체를 찾습니다. 상위 후보는 구조(ESMFold)로 재검증합니다."))
    yield _Ev("reasoning", text=(
        "화물이 짧아 셔틀의 표면 노출이 결합에 결정적입니다. 유연 링커(GS 계열)를 우선 검토합니다."))
    yield _Ev("evaluation", text="1차 후보 평가", data={"rows": rows})
    for i, b in enumerate([0.72, 0.85, 0.90, 0.94]):
        yield _Ev("progress", data={"best_bbb": b, "round": i + 1})
    yield _Ev("structure", text="상위 후보에서 셔틀이 표면에 충분히 노출됩니다.",
              data={"exposed": True})
    win = dict(r_win)
    win["agent_pick"] = True
    yield _Ev("choice", data=win)
    yield _Ev("final", text=(
        f"**최종 보고서** — `{r_win['shuttle_name']} · {r_win['linker_name']}` 조합이 "
        f"BBB {r_win['bbb']*100:.0f}점, 독성 {r_win['tox']*100:.0f}%로 최적입니다. "
        "안정성·선택성 모두 안전 범위이며, 실제 적용 전 합성·In Vitro 검증이 필요합니다."))
    yield _Ev("critique", text=(
        "BBB·독성·선택성 근거가 일관되고 임계값을 만족합니다. 최종 선택을 승인합니다."),
        data={"approve": True})
    yield _Ev("reflection", text=(
        "추가로 검증 셔틀의 잔기 편집(design_candidate)이나 링커·셔틀 co-evolution으로 "
        "라이브러리 시드에서 더 나은 후보를 탐색할 수 있습니다."))


def _emit_agent_event(ev):
    """최적화 에이전트 이벤트 한 개 렌더 (plan/reflection/final/optimum/progress는 요약에서도 처리)."""
    if ev.kind == "plan":
        with st.container(border=True):
            st.markdown("**계획 — 에이전트 전략**")
            st.markdown(ev.text)
    elif ev.kind == "reflection":
        with st.container(border=True):
            st.markdown("**자기평가 — 에이전트 반성**")
            st.markdown(ev.text)
    elif ev.kind == "critique":
        with st.container(border=True):
            st.markdown(
                f"**심사 에이전트 — 적대적 검증** &nbsp; {_approve_pill(ev.data.get('approve'))}",
                unsafe_allow_html=True,
            )
            st.markdown(ev.text)
    elif ev.kind == "reasoning":
        with st.container(border=True):
            st.caption(ev.text)
    elif ev.kind == "text":
        st.markdown(ev.text)
    elif ev.kind == "evaluation":
        st.markdown(f"**도구 — 후보 평가 (BBB·독성)**: {ev.text or '(개선 후보)'}")
        rows = ev.data.get("rows", [])
        if rows:
            def _sh(s):
                return (s[:14] + "…") if len(s) > 15 else (s or "—")
            st.dataframe(
                _style_verdict({
                    "라벨": [r["label"] for r in rows],
                    "링커": [r["linker"] or "—" for r in rows],
                    "셔틀": [_sh(r["shuttle"]) for r in rows],
                    "BBB점": [f"{r['bbb']*100:.0f}" for r in rows],
                    "독성": [f"{r['tox']*100:.0f}%" for r in rows],
                    "수용체": [f"{r.get('bind_ref','?')}·{r.get('bind_score',0):.2f}" for r in rows],
                    "안정성(II)": [f"{r.get('instability','?')}·"
                                   f"{'안정' if r.get('stable') else '불안정'}" for r in rows],
                    "개발성": [f"{r.get('dev_risk','?')}·L{r.get('dev_liab','?')}·q{r.get('dev_charge','?')}"
                               for r in rows],
                    "판정": ["독성" if r["toxic"] else "통과" for r in rows],
                }, _PASS_CELL),
                width='stretch', hide_index=True,
            )
    elif ev.kind == "structure":
        d = ev.data
        icon = "노출" if d.get("exposed") else "가림 우려"
        st.markdown(f"**도구 — 구조 검증 (ESMFold)** · {icon} — {ev.text}")
    elif ev.kind == "error":
        st.error(ev.text)


def _eff_bbb(r):
    """에이전트 결정축과 동일한 **유효점수** = BBB − λ·off_target_risk.
    raw BBB가 최고여도 비특이(off-target) CPP면 감점돼 밀린다(reward-hacking 회피).
    표시축(랭킹)을 이 결정축과 일치시켜 '왜 1위를 안 골랐나' 불일치를 없앤다."""
    return r.get("bbb", 0.0) - OFF_TARGET_PENALTY * r.get("sel_off", 0.0)


def _extract_podium(events, n=3):
    """상위 후보 n개. 에이전트 결정축과 동일한 **유효점수(eff = BBB − λ·off_target)** 순으로
    정렬한다 — raw BBB 최고여도 off-target CPP면 밀린다. 에이전트가 finish로 고른 최종
    선택(choice)이 있으면 1위로 고정하고, 나머지는 유효점수 순으로 채운다(중복 제거)."""
    seen = {}
    for e in events:
        if e.kind == "evaluation":
            for r in e.data.get("rows", []):
                if r.get("toxic"):
                    continue
                seq = r["sequence"]
                if seq not in seen or _eff_bbb(r) > _eff_bbb(seen[seq]):
                    seen[seq] = r
    ranked = sorted(seen.values(), key=_eff_bbb, reverse=True)
    choice = next((e.data for e in events if e.kind == "choice"), None)
    if choice:
        cseq = choice.get("sequence")
        rest = [r for r in ranked if r.get("sequence") != cseq]
        return [choice] + rest[:n - 1]
    return ranked[:n]

def _render_agent_summary(events, cargo):
    """최적화 궤적 + 상위 3 후보 카드(+온디맨드 분석) + 보고서."""
    progress = [e.data for e in events if e.kind == "progress"]
    podium = _extract_podium(events)
    plan = next((e.text for e in events if e.kind == "plan"), None)
    final = next((e.text for e in events if e.kind == "final"), None)
    reflection = next((e.text for e in events if e.kind == "reflection"), None)
    critique = next((e for e in events if e.kind == "critique"), None)
    if plan:
        with st.container(border=True):
            st.markdown("### 에이전트 전략")
            st.markdown(plan)
    if len(progress) >= 2:
        st.markdown("##### 최적화 궤적")
        _traj = pd.DataFrame({
            "라운드": list(range(1, len(progress) + 1)),
            "BBB 투과 점수": [round(p["best_bbb"] * 100) for p in progress],
        })
        _chart = (
            alt.Chart(_traj)
            .mark_line(point=True, color="#30405c")
            .encode(
                x=alt.X("라운드:O", title="라운드"),
                y=alt.Y("BBB 투과 점수:Q", title="BBB 투과 점수",
                        scale=alt.Scale(domain=[0, 100])),  # 0~100 고정, 음수 제거
            )
        )
        st.altair_chart(_chart, use_container_width=True)
    if podium:
        st.markdown("##### 상위 후보")
        st.caption("에이전트가 고른 최종 선택(심사 판정 포함)이 1위이고, 나머지는 **유효점수"
                   "(BBB − off-target 감점) 순** 차순위입니다 — raw BBB 최고여도 비특이(off-target) "
                   "CPP는 밀립니다(에이전트 결정축과 동일).")
        cols = st.columns(len(podium))
        for i, cand in enumerate(podium):
            with cols[i]:
                approved = cand.get("critic_approved")
                is_pick = cand.get("agent_pick") or approved is not None
                with st.container(border=True, key=f"best-card-{i}"):
                    if is_pick and approved is False:
                        # 1위(최종 선택)로는 제시하되, 심사 미승인임을 명시한다.
                        st.markdown("### 에이전트 최종 선택")
                        st.caption("심사 에이전트가 개선요구(REVISE) — 확정 아님. 스텝 수를 늘려 "
                                   "재실행하면 개선안을 냅니다.")
                        _v = Verdict.SUBOPTIMAL
                    elif is_pick:
                        st.markdown("### 에이전트 최종 선택")
                        _v = Verdict.ACCEPTED
                    else:
                        st.markdown(f"### {i + 1}위")
                        _v = (Verdict.REJECTED_TOXIC if cand.get("toxic")
                              else Verdict.SUBOPTIMAL)
                    st.markdown(_verdict_pill(_v), unsafe_allow_html=True)
                    st.code(cand["sequence"], language="text")
                    m1, m2 = st.columns(2)
                    m1.metric("BBB 투과 점수", f"{cand['bbb']*100:.0f}")
                    m2.metric("독성", f"{cand['tox']*100:.0f}%", delta="안전", delta_color="off")
                    if cand.get("mechanism"):
                        _mech = cand["mechanism"]
                        _pres = cand.get("preservation")
                        _pres_s = f"{_pres*100:.0f}" if isinstance(_pres, (int, float)) else "n/a"
                        _avid = cand.get("avidity")
                        st.caption(
                            f"전달 분해 — 셔틀 내재 **{cand.get('shuttle_bbb', 0)*100:.0f}** × 융합보존 "
                            f"**{_pres_s}** · 메커니즘 **{_mech}**"
                            + (f" · avidity {_avid:.2f}" if isinstance(_avid, (int, float)) else "")
                            + (f"  \n:gray[{'RMT — 실제 전달은 수용체 결합·avidity가 결정, deepB3P는 참고치' if cand.get('is_rmt') else 'CPP — deepB3P 비교적 유효하나 비특이(off-target) 주의' if _mech == 'CPP' else 'deepB3P 절대값 신뢰 낮음 — 상대 비교'}]"))
                    st.caption(
                        f"불안정성 {cand.get('instability', '?')}"
                        f"({'안정' if cand.get('stable') else '불안정'}) · "
                        f"수용체 {cand.get('bind_ref', '?')}·{cand.get('bind_score', 0):.2f}")
                    _devs = cand.get("dev_liabilities") or []
                    st.caption(
                        f"개발성 위험 **{cand.get('dev_risk', '?')}** · 순전하 {cand.get('dev_charge', '?')} · "
                        f"liability {cand.get('dev_liab', 0)}개 · 용해도 **{cand.get('sol_level', '?')}**"
                        f"({cand.get('sol_score', '?')})"
                        + (f" ({', '.join(_devs[:2])}…)" if _devs else ""))
                    st.caption(
                        f"선택성 **{cand.get('selectivity', '?')}** · off-target 위험 "
                        f"**{cand.get('sel_level', '?')}** · {cand.get('sel_mech', '?')}")
                    st.caption(f"링커 `{cand['linker'] or '—'}` · 셔틀 `{cand['shuttle'] or '—'}`")
    else:
        st.warning("독성 임계값을 통과한 후보가 없습니다. 스텝 수를 늘리거나 화물을 바꿔 다시 시도해 보세요.")
    if final:
        with st.container(border=True):
            st.markdown("### 설계 에이전트 — 최종 보고서")
            st.markdown(final)
    if critique is not None:
        with st.container(border=True):
            st.markdown(
                f"### 심사 에이전트 — 적대적 검증  &nbsp; {_approve_pill(critique.data.get('approve'))}",
                unsafe_allow_html=True,
            )
            st.markdown(critique.text)
    if reflection:
        with st.container(border=True):
            st.markdown("### 에이전트 자기평가")
            st.markdown(reflection)

# --- 실행 -------------------------------------------------------------------
# --- 자율 가설 에이전트 실행(버튼): 라이브 스트리밍 + session_state 저장 ---
if agent_run:
    cargo = _clean_cargo(cargo_input)
    _err = _cargo_error(cargo)
    if _err:
        st.error(_err)
        st.stop()
    if settings.use_gemini_agent:
        event_source = get_gemini_agent(settings, max_rounds=agent_rounds).run(cargo)
        brain = f"Gemini ({settings.gemini_model})"
    else:
        event_source = _demo_agent_events(cargo)   # API 키 없음 → 예시 데이터로 화면 미리보기
        brain = "데모 모드 (예시 데이터)"
    st.session_state["agent_analysis"] = {}  # 새 실행 → 이전 온디맨드 분석 초기화
    st.divider()
    st.subheader("자율 설계 에이전트")
    st.caption(
        f"화물 `{cargo}` · **{brain}**가 BBB·독성·구조·잔기 편집·서열 진화 도구를 자율 오케스트레이션해 "
        f"최종 융합체를 탐색 (최대 {agent_rounds}스텝)"
    )
    events = []
    with st.status("에이전트가 후보를 제안하고 평가 중...", expanded=True) as status:
        for ev in event_source:
            events.append(ev)
            _emit_agent_event(ev)  # 라이브 스트리밍(완료 시 접힘 → 과정 기록)
            if ev.kind == "error":
                status.update(label="중단됨", state="error")
        status.update(label="최적화 완료 · 위 status를 펼치면 과정 기록", state="complete",
                      expanded=False)
    _render_agent_summary(events, cargo)
    if settings.use_gemini_agent:
        st.caption("최적화 판단은 deepB3P·ToxinPred3 예측 기반입니다. 실제 합성·검증 필요.")
    else:
        st.caption("※ 예시(더미) 데이터입니다 — 실제 예측·설계 결과가 아니라 화면 미리보기용입니다.")
    _rec = {"cargo": cargo, "events": events, "rounds": agent_rounds, "brain": brain}
    _opt = (next((e.data for e in events if e.kind == "choice"), None)
            or next((e.data for e in events if e.kind == "optimum"), None))
    _sm = (f"{_opt.get('sequence', '')[:18]}… BBB {_opt.get('bbb', 0)*100:.0f}"
           if _opt else "결과 없음")
    _hist = st.session_state.setdefault("agent_runs", [])
    _rec["label"] = f"#{len(_hist)+1} · 화물 {cargo[:9]} → {_sm}"
    _hist.append(_rec)
    del _hist[:-3]                       # 최근 3개만 유지
    st.session_state["agent"] = _rec
    st.session_state["view"] = "agent"
    st.stop()

_view = st.session_state.get("view")
_agent = st.session_state.get("agent")
if _view == "agent" and _agent:
    # 이전 최적화 실행 결과 다시 렌더 (rerun 유지)
    st.divider()
    st.subheader("자율 설계 에이전트")
    _runs = st.session_state.get("agent_runs", [])
    if len(_runs) > 1:                    # 실행 기록 선택 (최근 3개, 최신 먼저)
        _rev = _runs[::-1]
        _labels = [r.get("label", "실행") for r in _rev]
        _pick = st.selectbox("실행 기록 (최근 3개)", _labels, index=0)
        _agent = _rev[_labels.index(_pick)]
    st.caption(f"화물 `{_agent['cargo']}` · {_agent.get('brain', 'LLM')} 다중 도구 자율 설계 "
               f"(최대 {_agent['rounds']}스텝)")
    _render_agent_summary(_agent["events"], _agent["cargo"])
    with st.expander("중간 과정 기록 — 스텝별 추론·평가·검증 다시 보기"):
        for ev in _agent["events"]:
            if ev.kind not in ("plan", "reflection", "critique"):  # 위 요약에 이미 표시
                _emit_agent_event(ev)
    st.caption("최적화 판단은 deepB3P·ToxinPred3 예측 기반입니다. 실제 합성·검증 필요.")
else:
    if settings.use_gemini_agent:
        with st.container(key="yy-more-input1"):
            st.info("화물 서열을 입력하고 ‘자율 설계 에이전트 실행’ 버튼을 누르세요.")
    else:
        with st.container(key="yy-more-input2"):
            st.info(
                "화물 서열을 입력하세요. 자율 설계 에이전트를 쓰려면 Gemini API 키가 "
                "필요합니다(위 안내 참고). 데모 모드로도 예시 결과를 미리볼 수 있습니다."
            )
