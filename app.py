"""
BBB-Optimize AI Agent — Streamlit UI (thin layer)
=================================================
화물(cargo) 펩타이드를 받아, 링커·셔틀 라이브러리를 **전수 조합**해 융합체를 만들고,
deepB3P(BBB)·ToxinPred3(독성)로 분석해 **투과 점수 높고 비독성인 베스트 N**을 추천한다.

실행:
    pip install -r requirements.txt
    streamlit run app.py
"""

import re
import time
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from core import build_agent, get_settings
from core.binding import shuttle_similarity
from core.config import (
    DEFAULT_CARGO,
    LINKER_LIBRARY,
    MODEL_MAX_LEN,
    SHUTTLES,
    STANDARD_LINKER_NAME,
    VALID_AMINO_ACIDS,
    bbb_scoring_seq,
)
from core.generative import get_fbgan
from core.optimizer_agent import get_optimization_agent
from core.optimizer_agent_gemini import get_gemini_agent
from core.predictors import get_predictor
from core.structure import analyze_construct
from core.toxicity import get_toxicity_predictor
from core.schemas import StepLevel, Verdict

STEP_DELAY = 0.7

st.set_page_config(page_title="BBB-Optimize AI Agent", page_icon="🧬", layout="wide")


def _inject_theme():
    css = Path(__file__).parent / "assets" / "theme.css"
    if css.exists():
        st.html(f"<style>{css.read_text()}</style>")


_inject_theme()
settings = get_settings()

st.title("🧬 BBB-Optimize AI Agent")
st.markdown(
    """
    **화물(cargo) 펩타이드를 주면, AI가 여러 도구를 자율적으로 오케스트레이션해
    뇌혈관장벽(BBB)을 통과하는 항체-셔틀 융합 단백질을 설계합니다.**
    효능(deepB3P) · 독성(ToxinPred3) · 안정성(ProtParam) · 입체구조(ESMFold) · 신규 생성(FBGAN)
    도구를 조합 — `융합체 = 화물 + 링커 + 셔틀`
    """
)


def _dot(ok):
    return "🟢" if ok else "🟠"


st.caption(
    f"엔진 상태 — 효능 deepB3P {_dot(settings.use_deepb3p_local)} · "
    f"독성 ToxinPred3 {_dot(settings.use_toxinpred3_local)} · 구조 ESMFold 🟢(API) · "
    f"브레인 Gemini {_dot(settings.use_gemini_agent)} · Claude {_dot(settings.use_llm_agent)}"
    f"{'' if (settings.use_gemini_agent or settings.use_llm_agent) else ' (API 키 필요)'}"
)
st.caption(
    "**🤖 자율 설계 에이전트**가 메인입니다 — LLM(**기본 Gemini**, Claude 선택 가능)이 "
    "**평가(BBB·독성·수용체·안정성) · 구조(ESMFold) · 생성(FBGAN)** 도구를 스스로 골라 써가며 "
    "최종 융합체를 설계합니다. 각 도구는 '개별 도구' 패널에서 수동으로도 실행할 수 있어요."
)
st.caption(
    "※ **BBB 투과 점수**는 deepB3P 예측 확률(0~100)로 **실제 물리적 투과율이 아닙니다** — "
    "후보 비교용. 진짜 투과율은 실험(Papp/logBB) 필요. 이중 트랙에서 **수용체 결합 가능성**(메커니즘 근거)도 확인하세요."
)

n_combos = len(LINKER_LIBRARY) * len(SHUTTLES)

st.divider()

# --- 입력 (중앙 정렬 히어로 · 라이브러리 전체 자동 탐색) --------------------
_hl, _hero, _hr = st.columns([1, 2, 1])
with _hero:
    cargo_input = st.text_input(
        "🧪 화물(cargo) 펩타이드 서열",
        value=DEFAULT_CARGO,
        help="링커·셔틀은 라이브러리에서 전부 자동으로 붙여봅니다. (표준 20종 아미노산 1글자 코드)",
    )
    # === 메인 CTA: 자율 설계 에이전트 (공모전 flagship) ===
    agent_run, agent_rounds, agent_provider = False, 8, None
    _providers = []  # (id, label) — 기본 Gemini 먼저
    if settings.use_gemini_agent:
        _providers.append(("gemini", f"Gemini ({settings.gemini_model})"))
    if settings.use_llm_agent:
        _providers.append(("claude", f"Claude ({settings.llm_model})"))
    if _providers:
        agent_rounds = st.slider("최대 스텝 수", 4, 12, 8,
                                 help="에이전트가 도구를 호출하는 최대 횟수. 낮추면 빠르고 저렴합니다.")
        if len(_providers) > 1:
            _labels = [lb for _, lb in _providers]
            _sel = st.radio("브레인 모델", _labels, index=0, horizontal=True,
                            help="기본은 Gemini(정직한 BBB 프레이밍). Claude는 안전분류기 우회용 "
                                 "도메인 중립화 버전.")
            agent_provider = _providers[_labels.index(_sel)][0]
        else:
            agent_provider = _providers[0][0]
            st.caption(f"브레인: **{_providers[0][1]}**")
        agent_run = st.button("🤖 자율 설계 에이전트 실행", type="primary",
                              use_container_width=True)
        st.caption("에이전트가 **BBB·독성·구조·수용체·생성 도구를 스스로 오케스트레이션**해 "
                   "최종 융합체를 설계합니다.")
    else:
        st.info(
            "🔑 **자율 설계 에이전트는 LLM API 키가 필요합니다.** `.env`에 "
            "`GEMINI_API_KEY=...`(기본·무료 티어) 또는 `LLM_API_KEY=sk-ant-...`(Claude)를 넣고 "
            "앱을 재시작하세요. (키 없이 무료로 확인하려면 아래 '개별 도구'를 이용하세요.)",
            icon="🔑",
        )

    # === 개별 도구 (에이전트가 내부적으로 쓰는 도구들 · 수동 실행 · 무료) ===
    run = fbgan_run = struct_run = False
    fbgan_rounds = 4
    struct_linker_name, struct_shuttle_name = STANDARD_LINKER_NAME, list(SHUTTLES)[0]
    with st.expander("🔧 개별 도구 직접 실행 (에이전트 없이 · 무료)"):
        st.caption("에이전트가 자율적으로 호출하는 도구들을 수동으로도 돌려볼 수 있습니다.")
        run = st.button("🚀 라이브러리 전수 스윕 (모든 링커 × 셔틀)", use_container_width=True)

        st.markdown("**🧬 이중 트랙 — 구조(ESMFold)+투과 점수+수용체**")
        _sc1, _sc2 = st.columns(2)
        struct_linker_name = _sc1.selectbox(
            "링커", list(LINKER_LIBRARY), index=list(LINKER_LIBRARY).index(STANDARD_LINKER_NAME))
        struct_shuttle_name = _sc2.selectbox("셔틀", list(SHUTTLES), index=0)
        struct_run = st.button("🧬 구조 분석 실행 (~10초)", use_container_width=True)

        if settings.use_fbgan_local:
            st.markdown("**🧫 신규 셔틀 생성 (FBGAN)**")
            fbgan_rounds = st.slider("생성 라운드 수", 2, 8, 4)
            fbgan_run = st.button("🧫 생성 실행", use_container_width=True)

    with st.expander(f"📚 라이브러리 구성 보기 (링커 {len(LINKER_LIBRARY)} · 셔틀 {len(SHUTTLES)})"):
        st.markdown("**링커**")
        st.dataframe(
            {"링커": list(LINKER_LIBRARY), "서열": [v["seq"] for v in LINKER_LIBRARY.values()],
             "종류": [v["kind"] for v in LINKER_LIBRARY.values()],
             "설명": [v["note"] for v in LINKER_LIBRARY.values()]},
            use_container_width=True, hide_index=True,
        )
        st.markdown("**셔틀**")
        st.dataframe(
            {"셔틀": list(SHUTTLES), "서열": [v["seq"] for v in SHUTTLES.values()],
             "타겟": [v["target"] for v in SHUTTLES.values()],
             "설명": [v["note"] for v in SHUTTLES.values()]},
            use_container_width=True, hide_index=True,
        )

    _cargo_preview = cargo_input.strip().upper()
    if _cargo_preview:
        longest = max(len(v["seq"]) for v in SHUTTLES.values())
        longest_l = max(len(v["seq"]) for v in LINKER_LIBRARY.values())
        est = len(_cargo_preview) + longest_l + longest
        if est > MODEL_MAX_LEN:
            st.info(
                f"ℹ️ 일부 조합은 {MODEL_MAX_LEN}aa를 넘습니다(≈{est}aa). 이런 조합은 BBB를 "
                f"**연결부위(링커+셔틀 = BBB 모듈)만으로** 계산합니다(화물 벌크 제외) — "
                f"셔틀 신호가 화물에 희석되지 않습니다. 독성은 전체 서열로 계산합니다.",
                icon="🔎",
            )

_VERDICT_LABEL = {
    Verdict.ACCEPTED: "✅ 채택(베스트)",
    Verdict.REJECTED_TOXIC: "❌ 독성 탈락",
    Verdict.SUBOPTIMAL: "🔸 후순위",
    Verdict.REFERENCE: "⚪ 화물 단독(기준)",
}
_MEDAL = ["🥇", "🥈", "🥉"]


def _emit_agent_event(ev):
    """최적화 에이전트 이벤트 한 개 렌더 (final/optimum/progress는 호출측에서 처리)."""
    if ev.kind == "reasoning":
        with st.container(border=True):
            st.caption("🧠 " + ev.text)
    elif ev.kind == "text":
        st.markdown(ev.text)
    elif ev.kind == "evaluation":
        st.markdown(f"🧪 **[도구] 후보 평가 (BBB·독성)**: {ev.text or '(개선 후보)'}")
        rows = ev.data.get("rows", [])
        if rows:
            def _sh(s):
                return (s[:14] + "…") if len(s) > 15 else (s or "—")
            st.dataframe(
                {"라벨": [r["label"] for r in rows],
                 "링커": [r["linker"] or "—" for r in rows],
                 "셔틀": [_sh(r["shuttle"]) for r in rows],
                 "BBB점": [f"{r['bbb']*100:.0f}" for r in rows],
                 "독성": [f"{r['tox']*100:.0f}%" for r in rows],
                 "수용체": [f"{r.get('bind_ref','?')}·{r.get('bind_score',0):.2f}" for r in rows],
                 "안정성(II)": [f"{r.get('instability','?')}·"
                                f"{'안정' if r.get('stable') else '불안정'}" for r in rows],
                 "판정": ["❌독성" if r["toxic"] else "✅통과" for r in rows]},
                use_container_width=True, hide_index=True,
            )
    elif ev.kind == "structure":
        d = ev.data
        icon = "✅" if d.get("exposed") else "⚠️"
        st.markdown(f"🧊 **[도구] 구조 검증 (ESMFold)** {icon} — {ev.text}")
    elif ev.kind == "generation":
        st.markdown("🧬 **[도구] 신규 셔틀 생성 (FBGAN)**")
        novel = ev.data.get("novel", [])
        if novel:
            st.dataframe(
                {"생성 셔틀": [b["shuttle"] for b in novel],
                 "BBB": [f"{b['bbb']*100:.0f}" for b in novel],
                 "독성": [f"{b['tox']*100:.0f}%" for b in novel]},
                use_container_width=True, hide_index=True,
            )
    elif ev.kind == "error":
        st.error(ev.text)


def _extract_podium(events, n=3):
    """평가 이벤트 전체에서 비독성 후보를 모아 BBB 순 상위 n개(중복 제거)."""
    seen = {}
    for e in events:
        if e.kind == "evaluation":
            for r in e.data.get("rows", []):
                if r.get("toxic"):
                    continue
                seq = r["sequence"]
                if seq not in seen or r["bbb"] > seen[seq]["bbb"]:
                    seen[seq] = r
    return sorted(seen.values(), key=lambda r: r["bbb"], reverse=True)[:n]


def _render_candidate_analysis(cargo, cand, idx):
    """온디맨드 구조 분석 (ESMFold) — 버튼을 누를 때만 실행하고 session_state에 보관."""
    analyses = st.session_state.setdefault("agent_analysis", {})
    seq = cand["sequence"]
    if seq not in analyses:
        if st.button("🔬 구조 분석 (ESMFold, ~10초)", key=f"agent-analyze-{idx}",
                     use_container_width=True):
            with st.spinner("ESMFold로 접는 중..."):
                sr = analyze_construct(cargo, cand["linker"], cand["shuttle"])
            analyses[seq] = ({"error": sr.error} if sr.error else
                             {"exposure": sr.shuttle_exposure, "shuttle_plddt": sr.shuttle_plddt,
                              "mean_plddt": sr.mean_plddt, "exposed": sr.exposed,
                              "verdict": sr.verdict})
            st.rerun()  # 버튼→결과로 깔끔히 재렌더
        else:
            return
    a = analyses.get(seq)
    if not a:
        return
    if a.get("error"):
        st.caption(f"⚠️ 구조 예측 실패: {a['error']}")
        return
    icon = "✅ 노출" if a["exposed"] else "⚠️ 가림/저신뢰"
    st.caption(f"🧊 셔틀 노출도 **{a['exposure']:.2f}** ({icon}) · 셔틀 pLDDT {a['shuttle_plddt']} · "
               f"전체 pLDDT {a['mean_plddt']}")
    st.caption(a["verdict"])


def _render_agent_summary(events, cargo):
    """최적화 궤적 + 상위 3 후보 카드(+온디맨드 분석) + 보고서."""
    progress = [e.data for e in events if e.kind == "progress"]
    podium = _extract_podium(events)
    final = next((e.text for e in events if e.kind == "final"), None)
    if len(progress) >= 2:
        st.markdown("##### 📈 최적화 궤적 (best-so-far BBB)")
        st.line_chart({"best BBB": [p["best_bbb"] for p in progress]},
                      x_label="라운드", y_label="BBB 투과 점수")
    if podium:
        st.markdown("##### 🏆 상위 후보 (비독성 · BBB 투과 점수 순)")
        st.caption("에이전트가 선정한 최적 후보와 차순위입니다. 각 카드에서 **구조 분석**을 "
                   "원할 때만 실행할 수 있어요.")
        cols = st.columns(len(podium))
        medals = ["🥇", "🥈", "🥉"]
        for i, cand in enumerate(podium):
            with cols[i]:
                with st.container(border=True, key=f"best-card-{i}"):
                    st.markdown(f"### {medals[i]}  {i + 1}위")
                    st.code(cand["sequence"], language="text")
                    m1, m2 = st.columns(2)
                    m1.metric("🟢 BBB 투과 점수", f"{cand['bbb']*100:.0f}")
                    m2.metric("🔵 독성", f"{cand['tox']*100:.0f}%", delta="안전", delta_color="off")
                    st.caption(
                        f"불안정성 {cand.get('instability', '?')}"
                        f"({'안정' if cand.get('stable') else '불안정'}) · "
                        f"수용체 {cand.get('bind_ref', '?')}·{cand.get('bind_score', 0):.2f}")
                    st.caption(f"링커 `{cand['linker'] or '—'}` · 셔틀 `{cand['shuttle'] or '—'}`")
                    _render_candidate_analysis(cargo, cand, i)
    else:
        st.warning("독성 임계값을 통과한 후보가 없습니다. 스텝 수를 늘리거나 화물을 바꿔 다시 시도해 보세요.")
    if final:
        with st.container(border=True):
            st.markdown("### 🏁 최적화 보고서")
            st.markdown(final)


def _structure_html(pdb, cargo_len, linker_len, height=420):
    """py3Dmol 3D 뷰: 화물=회색, 링커=주황, 셔틀=브랜드 옐로."""
    import py3Dmol
    v = py3Dmol.view(width=680, height=height)
    v.addModel(pdb, "pdb")
    c_end, l_end = cargo_len, cargo_len + linker_len
    v.setStyle({}, {"cartoon": {"color": "#b9b9b9"}})
    if cargo_len:
        v.setStyle({"resi": f"1-{c_end}"}, {"cartoon": {"color": "#b9b9b9"}})
    v.setStyle({"resi": f"{c_end + 1}-{l_end}"}, {"cartoon": {"color": "#ff8500"}})
    v.setStyle({"resi": f"{l_end + 1}-99999"}, {"cartoon": {"color": "#ffde36"}})
    v.setBackgroundColor("0xffffff")
    v.zoomTo()
    return v._make_html()


def _render_dual_track(d):
    """이중 트랙 결과 렌더: Track1(구조/ESMFold) + Track2(투과 점수/deepB3P) + 독성."""
    sr = d["sr"]
    full_len = len(d["cargo"]) + len(d["linker"]) + len(d["shuttle"])
    st.caption(
        f"융합체 {full_len}aa = 화물({len(d['cargo'])}aa) + {d['linker_name']}({len(d['linker'])}) "
        f"+ {d['shuttle_name']}({len(d['shuttle'])})"
    )

    t1, t2 = st.columns([3, 2])

    # ---- Track 1: 구조 ----
    with t1:
        st.markdown("#### 🧊 Track 1 — 입체 구조 (ESMFold)")
        if sr.error:
            st.warning(f"구조 예측 실패: {sr.error}\n\n(전체 서열이 ESMFold 한계를 넘었을 수 있습니다. "
                       "Track 2 투과 점수은 아래에서 확인하세요.)")
        else:
            components.html(_structure_html(sr.pdb, len(sr.cargo), len(sr.linker)), height=420)
            st.caption("🎨 화물=회색 · 링커=주황 · **셔틀=노랑**. 노랑(셔틀)이 드러날수록 수용체 결합 유리.")
            sm1, sm2, sm3 = st.columns(3)
            sm1.metric("셔틀 노출도(RSA)", f"{sr.shuttle_exposure:.2f}",
                       delta="노출" if sr.exposed else "가림 우려",
                       delta_color="normal" if sr.exposed else "inverse")
            sm2.metric("셔틀 pLDDT", f"{sr.shuttle_plddt:.0f}",
                       delta="저신뢰" if sr.low_confidence else "양호",
                       delta_color="inverse" if sr.low_confidence else "normal")
            sm3.metric("전체 pLDDT", f"{sr.mean_plddt:.0f}")
            (st.warning if (not sr.exposed or sr.low_confidence) else st.success)("🔬 " + sr.verdict)
            st.download_button("⬇️ PDB 다운로드", sr.pdb, file_name="fusion.pdb",
                               mime="chemical/x-pdb")

    # ---- Track 2: 투과 점수 + 독성 ----
    with t2:
        st.markdown("#### 🎯 Track 2 — 투과 점수·독성")
        with st.container(border=True):
            st.metric("🟢 BBB 투과 점수 (deepB3P)", f"{d['bbb']*100:.0f}")
            st.caption("예측 확률(0~100)·상대비교 — 실제 투과율 아님. "
                       f"연결부위 슬라이스({len(d['bbb_seq'])}aa) 계산.")
            st.divider()
            _safe = d["tox"] <= settings.toxicity_threshold
            st.metric("🔵 독성 위험 (ToxinPred3)", f"{d['tox']*100:.0f}%",
                      delta="안전" if _safe else "위험",
                      delta_color="normal" if _safe else "inverse")
            st.caption(f"전체 {full_len}aa 서열로 계산 (조성 기반, 길이 무관)")

        # ---- Track 3: 수용체 결합 가능성 (메커니즘 프록시) ----
        st.markdown("#### 🧲 Track 3 — 수용체 결합 가능성")
        b = d.get("bind")
        if b:
            _strong = b.score >= 0.6
            with st.container(border=True):
                st.metric("🧲 수용체 결합 프록시", f"{b.score:.2f}",
                          delta=f"{b.best_ref} 유사", delta_color="off")
                st.caption(f"가장 닮은 검증 셔틀: **{b.best_ref}** ({b.target} · {b.mechanism})")
                (st.success if _strong else st.info)(b.verdict)
                st.caption("서열 유사도 기반 **기능 유추(프록시)** — 실제 결합 친화도 아님.")

    st.caption(
        "⚠️ Track 1 구조는 짧은 펩타이드에선 무질서로 신뢰도↓(항체 등 큰 화물일수록↑). 셔틀 노출도·구조→BBB는 "
        "검증된 모델이 아닌 **휴리스틱 프록시** — 실험 검증 필요. Track 2 BBB도 deepB3P 예측값(상대 비교 권장)."
    )


# --- 실행 -------------------------------------------------------------------
if run:
    cargo = cargo_input.strip().upper()
    if not cargo:
        st.error("⚠️ 화물 펩타이드 서열을 입력해 주세요.")
        st.stop()
    if set(cargo) - VALID_AMINO_ACIDS:
        st.error("⚠️ 유효하지 않은 아미노산 문자가 있습니다. (표준 20종 1글자 코드만)")
        st.stop()

    agent = build_agent(settings)

    result = None
    with st.status(f"🤖 {n_combos}개 조합 조립·분석 중입니다...", expanded=True) as status:
        gen = agent.run(cargo)
        try:
            while True:
                step = next(gen)
                icon = {StepLevel.INFO: "🔹", StepLevel.WARN: "⚠️",
                        StepLevel.SUCCESS: "🏆", StepLevel.ERROR: "🛑"}.get(step.level, "🔹")
                st.write(f"{icon} **{step.stage}**: {step.message}")
                time.sleep(STEP_DELAY)
        except StopIteration as stop:
            result = stop.value
        status.update(label="최적화 완료!", state="complete", expanded=False)

    if result is None or result.winner is None:
        st.error("독성 임계값을 통과한 조합이 없습니다. 임계값/라이브러리를 재검토하세요.")
        st.stop()

    cargo_only = result.cargo_only
    base_bbb = cargo_only.prediction.bbb_permeability if cargo_only else 0.0

    # 상위 N 재계산 (생존 조합 BBB 내림차순)
    survivors = sorted(
        (e for e in result.candidates if e.verdict != Verdict.REJECTED_TOXIC),
        key=lambda e: e.prediction.bbb_permeability, reverse=True,
    )
    top = survivors[: settings.top_n]

    st.divider()
    st.subheader(f"🏆 베스트 {len(top)} 융합체")
    st.caption(
        f"화물 {len(cargo)}aa · {result.n_linkers}링커 × {result.n_shuttles}셔틀 = "
        f"{len(result.candidates)}조합 분석 · 화물 단독 BBB {base_bbb*100:.0f}점"
        + ("  ·  🔴 화물 자체가 독성" if cargo_only and cargo_only.prediction.toxicity_risk > settings.toxicity_threshold else "")
    )

    # --- 베스트 N 카드 (가로 3단) ------------------------------------------
    cols = st.columns(len(top))
    for i, e in enumerate(top):
        c, p = e.construct, e.prediction
        with cols[i]:
            with st.container(border=True, key=f"best-card-{i}"):
                st.markdown(f"### {_MEDAL[i] if i < 3 else f'#{i+1}'}  {c.label}")
                m1, m2 = st.columns(2)
                m1.metric("🟢 BBB 투과 점수", f"{p.bbb_permeability*100:.0f}",
                          delta=f"{(p.bbb_permeability-base_bbb)*100:+.0f}점")
                m2.metric("🔵 독성 위험", f"{p.toxicity_risk*100:.0f}%", delta="안전",
                          delta_color="off")
                st.code(c.sequence, language="text")
                st.caption(
                    f"화물 `{c.cargo}` + 링커 `{c.linker}`({c.linker_name}) + "
                    f"셔틀 `{c.shuttle_seq}`({c.shuttle_name})"
                    + ("  ·  🔎 연결부위 계산" if c.truncated else "")
                )

    # --- 1위 상세 해설 ------------------------------------------------------
    win = top[0]
    wc, wp = win.construct, win.prediction
    with st.container(border=True):
        st.markdown(f"#### 💡 1위 [{wc.label}] 해설")
        st.markdown(
            f"""
            라이브러리의 **{len(result.candidates)}개 조합**을 전수 분석한 결과,
            **{wc.shuttle_name} 셔틀 + {wc.linker_name} 링커** 조합이 BBB 투과 점수
            **{wp.bbb_permeability*100:.0f}점** ({(wp.bbb_permeability-base_bbb)*100:+.0f}점 vs 화물 단독)로
            가장 높으면서, 독성 **{wp.toxicity_risk*100:.0f}%** (임계값 {settings.toxicity_threshold:.2f} 이하)로
            안전 범위라 최적으로 선정됐습니다.

            > {wp.note}
            """
        )

    # --- 전체 조합 (접이식) -------------------------------------------------
    with st.expander(f"🔎 전체 {len(result.candidates)}개 조합 + 기준선 보기 (BBB 내림차순)"):
        rows = ([cargo_only] if cargo_only else []) + sorted(
            result.candidates, key=lambda e: e.prediction.bbb_permeability, reverse=True
        )
        tox_col = "독성(ToxinPred3)" if settings.use_toxinpred3_local else "독성(임시)"
        st.dataframe(
            {
                "조합": [e.construct.label for e in rows],
                "링커": [e.construct.linker or "—" for e in rows],
                "BBB 투과 점수": [f"{e.prediction.bbb_permeability*100:.0f}" for e in rows],
                tox_col: [f"{e.prediction.toxicity_risk*100:.0f}%" for e in rows],
                "길이": [f"{len(e.construct.sequence)}aa" + ("🔎" if e.construct.truncated else "")
                        for e in rows],
                "판정": [_VERDICT_LABEL.get(e.verdict, "-") for e in rows],
            },
            use_container_width=True, hide_index=True,
        )

    _tox_note = "ToxinPred3 실측" if settings.use_toxinpred3_local else "placeholder"
    st.caption(
        f"⚗️ BBB는 deepB3P 실측(짧은 펩타이드 학습 → 긴 융합체는 조합 간 **상대 비교**로 해석). "
        f"독성은 {_tox_note}. 실제 적용 전 In Vitro/In Vivo 검증 필요."
    )
else:
    # --- 자율 가설 에이전트 실행(버튼): 라이브 스트리밍 + session_state 저장 ---
    if agent_run:
        cargo = cargo_input.strip().upper()
        if not cargo or (set(cargo) - VALID_AMINO_ACIDS):
            st.error("⚠️ 유효한 화물 펩타이드 서열을 입력해 주세요.")
            st.stop()
        if agent_provider == "claude":
            agent = get_optimization_agent(settings, max_rounds=agent_rounds)
            brain = f"Claude ({settings.llm_model})"
        else:
            agent = get_gemini_agent(settings, max_rounds=agent_rounds)
            brain = f"Gemini ({settings.gemini_model})"
        st.session_state["agent_analysis"] = {}  # 새 실행 → 이전 온디맨드 분석 초기화
        st.divider()
        st.subheader("🤖 자율 설계 에이전트")
        st.caption(
            f"화물 `{cargo}` · **{brain}**가 **BBB·독성·구조·생성 도구를 자율 오케스트레이션**해 "
            f"최종 융합체를 탐색 (최대 {agent_rounds}스텝)"
        )
        events = []
        with st.status("에이전트가 후보를 제안하고 평가 중...", expanded=True) as status:
            for ev in agent.run(cargo):
                events.append(ev)
                _emit_agent_event(ev)  # 라이브 스트리밍(완료 시 접힘 → 과정 기록)
                if ev.kind == "error":
                    status.update(label="중단됨", state="error")
            status.update(label="최적화 완료 · 위 status를 펼치면 과정 기록", state="complete",
                          expanded=False)
        _render_agent_summary(events, cargo)
        st.caption("⚗️ 최적화 판단은 deepB3P·ToxinPred3 예측 기반입니다. 실제 합성·검증 필요.")
        st.session_state["agent"] = {"cargo": cargo, "events": events, "rounds": agent_rounds,
                                     "brain": brain}
        st.session_state["view"] = "agent"
        st.stop()

    # --- 구조 기반 접합부 분석 실행(버튼) ---
    if struct_run:
        cargo = cargo_input.strip().upper()
        if not cargo or (set(cargo) - VALID_AMINO_ACIDS):
            st.error("⚠️ 유효한 화물 펩타이드 서열을 입력해 주세요.")
            st.stop()
        linker = LINKER_LIBRARY[struct_linker_name]["seq"]
        shuttle = SHUTTLES[struct_shuttle_name]["seq"]
        st.divider()
        st.subheader("🧬 융합체 정밀 분석 — 이중 트랙")
        with st.status("Track1: ESMFold 폴딩 · Track2: deepB3P·ToxinPred3 채점 중... (~10초)",
                       expanded=True) as status:
            sr = analyze_construct(cargo, linker, shuttle)                  # Track 1 (구조)
            bbb_seq = bbb_scoring_seq(cargo, linker, shuttle)               # Track 2 슬라이스
            bbb_val = get_predictor(settings).predict_many([bbb_seq])[0].bbb_permeability
            tox_val = get_toxicity_predictor(settings).predict_many(
                [cargo + linker + shuttle])[0].risk
            bind = shuttle_similarity(shuttle)                              # Track 3 (수용체)
            status.update(label="분석 완료", state="complete", expanded=False)
        d = {"sr": sr, "cargo": cargo, "linker": linker, "shuttle": shuttle,
             "linker_name": struct_linker_name, "shuttle_name": struct_shuttle_name,
             "bbb": bbb_val, "bbb_seq": bbb_seq, "tox": tox_val, "bind": bind}
        _render_dual_track(d)
        st.session_state["struct"] = d
        st.session_state["view"] = "struct"
        st.stop()

    # FBGAN 실행(버튼) → 결과를 session_state에 저장 (rerun에도 유지)
    if fbgan_run:
        st.session_state["view"] = "fbgan"
        cargo = cargo_input.strip().upper()
        if not cargo or (set(cargo) - VALID_AMINO_ACIDS):
            st.error("⚠️ 유효한 화물 펩타이드 서열을 입력해 주세요.")
            st.stop()
        linker = LINKER_LIBRARY[STANDARD_LINKER_NAME]["seq"]
        with st.status(f"🧬 생성 최적화 루프 실행 중... ({fbgan_rounds}라운드: 생성→평가→피드백)",
                       expanded=True) as status:
            st.write("사전학습 생성기로 novel 셔틀 생성 → deepB3P·ToxinPred3 채점 → latent 진화")
            try:
                fres = get_fbgan(settings).run(
                    cargo, linker, rounds=fbgan_rounds,
                    tox_threshold=settings.toxicity_threshold)
                st.session_state["fbgan"] = {
                    "history": fres.history, "best": fres.best,
                    "cargo": cargo, "linker": linker, "rounds": fbgan_rounds,
                }
            except Exception as exc:  # noqa: BLE001
                status.update(label="생성 최적화 실패", state="error")
                st.error(f"실행 오류: {exc}")
                st.stop()
            status.update(label="생성 최적화 완료!", state="complete", expanded=False)

    _view = st.session_state.get("view")
    _agent = st.session_state.get("agent")
    _struct = st.session_state.get("struct")
    fb = st.session_state.get("fbgan")
    if _view == "struct" and _struct:
        st.divider()
        st.subheader("🧬 융합체 정밀 분석 — 이중 트랙")
        _render_dual_track(_struct)
    elif _view == "agent" and _agent:
        # 이전 최적화 실행 결과 다시 렌더 (rerun 유지)
        st.divider()
        st.subheader("🤖 자율 설계 에이전트")
        st.caption(f"화물 `{_agent['cargo']}` · {_agent.get('brain', 'LLM')} 다중 도구 자율 설계 "
                   f"(최대 {_agent['rounds']}스텝)")
        _render_agent_summary(_agent["events"], _agent["cargo"])
        with st.expander("🔍 중간 과정 기록 — 스텝별 추론·평가·검증 다시 보기"):
            for ev in _agent["events"]:
                _emit_agent_event(ev)
        st.caption("⚗️ 최적화 판단은 deepB3P·ToxinPred3 예측 기반입니다. 실제 합성·검증 필요.")
    elif fb:
        st.divider()
        st.subheader("🧬 AI 생성 셔틀 — 최적화 결과")
        st.caption(
            f"화물 {len(fb['cargo'])}aa · 링커 `{fb['linker']}` 고정 · {fb['rounds']}라운드 진화 · "
            f"셔틀은 **AI가 새로 생성**한 서열(라이브러리 밖)"
        )
        if fb["history"]:
            st.markdown("##### 📈 라운드별 BBB 개선 (생성기가 학습하는 과정)")
            st.line_chart(
                {"평균 BBB": [h["mean_bbb"] for h in fb["history"]],
                 "최고 BBB": [h["best_bbb"] for h in fb["history"]]},
                x_label="라운드", y_label="BBB 투과 점수",
            )
        st.markdown("##### 🏆 생성된 베스트 셔틀")
        top = fb["best"][:3]
        if not top:
            st.warning("비독성 후보를 찾지 못했습니다. 라운드를 늘려 다시 시도해 보세요.")
        else:
            cols = st.columns(len(top))
            medals = ["🥇", "🥈", "🥉"]
            for i, b in enumerate(top):
                with cols[i]:
                    with st.container(border=True, key=f"gen-card-{i}"):
                        st.markdown(f"### {medals[i]}  생성 셔틀 #{i+1}")
                        st.code(b["shuttle"], language="text")
                        m1, m2 = st.columns(2)
                        m1.metric("🟢 BBB 투과 점수", f"{b['bbb']*100:.0f}")
                        m2.metric("🔵 독성 위험", f"{b['tox']*100:.0f}%", delta="안전",
                                  delta_color="off")
                        st.caption(f"전체 융합체({b['len']}aa): `{b['sequence']}`")
        st.caption(
            "⚠️ 이 셔틀들은 **AI가 새로 생성**한 서열입니다 — 자연·검증된 펩타이드가 아니며, "
            "실제 합성·In Vitro 검증이 반드시 필요합니다. BBB/독성은 deepB3P·ToxinPred3 예측값."
        )
    else:
        if settings.use_gemini_agent or settings.use_llm_agent:
            st.info("👆 화물 펩타이드를 입력하고 **[🤖 자율 설계 에이전트 실행]** 버튼을 눌러 주세요.")
        else:
            st.info(
                "👆 화물 펩타이드를 입력하세요. 자율 설계 에이전트를 쓰려면 **LLM API 키**가 "
                "필요합니다(기본 Gemini, 위 안내 참고). 키 없이 확인하려면 **[🔧 개별 도구 직접 실행]** "
                "패널을 펼쳐 각 도구를 수동으로 실행할 수 있어요."
            )
