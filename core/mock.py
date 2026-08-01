"""디자인 미리보기용 목업 실행 — 엔진·LLM·API 없이 결과 화면을 렌더한다.

프론트(디자인) 작업자가 백엔드(deepB3P·ToxinPred3·Gemini 키) 설치 없이도
포디움 카드·심사 박스·과정 기록 등 **결과 화면 전체를 보고 디자인**할 수 있게 한다.
(값은 실제가 아닌 예시 목업 — `pip install streamlit`만으로 동작.)
"""

from __future__ import annotations

from .optimizer_agent import AgentEvent


def _row(label, linker, shuttle, cargo, **kw):
    base = dict(
        label=label, linker=linker, shuttle=shuttle, sequence=cargo + linker + shuttle,
        bbb=0.42, tox=0.24, toxic=False, bind_ref="Angiopep-2", bind_score=1.0,
        instability=45.3, stable=False,
        dev_risk="보통", dev_liab=2, dev_charge=2.0, dev_agg=0.3,
        dev_liabilities=["프로테아제 절단부위 KR@10", "탈아마이드화 NN@12"],
        sel_off=0.09, selectivity=0.91, sel_level="낮음",
        sel_mech="RMT형(Angiopep-2 유사·수용체 선택적)",
        sel_drivers=["RMT 수용체 표적 → 선택성↑"],
        sol_score=0.78, sol_level="높음")
    base.update(kw)
    return base


def mock_agent_run(cargo: str = "GSNKGAIIGLM") -> dict:
    """세션 저장 포맷과 동일한 목업 실행 레코드."""
    r1 = _row("Angio_EAAAK2", "EAAAKEAAAK", "TFFYGGSRGKRNNFKTEEY", cargo)
    r2 = _row("SynB1_Direct", "", "RGGRLSYSRRRFSTSTGR", cargo, bbb=0.96, tox=0.175,
              bind_ref="SynB1", bind_score=1.0, instability=53.9,
              selectivity=0.34, sel_level="높음", sel_mech="CPP형(SynB1 유사·비특이 막투과)",
              sel_drivers=["양전하 밀도 0.33 — 비특이 정전 흡수 위험"], dev_risk="낮음", dev_liab=1)
    r3 = _row("ApoE_G4S", "GGGGS", "LAVYQAGARLAVYQAGAR", cargo, bbb=0.06, tox=0.175,
              bind_ref="ApoE(159-167)2", bind_score=1.0, selectivity=0.75,
              sel_mech="RMT형(ApoE(159-167)2 유사·수용체 선택적)",
              dev_risk="낮음", sol_level="보통", sol_score=0.6)
    choice = dict(r1)
    choice["agent_pick"] = True

    events = [
        AgentEvent("plan", "먼저 검증된 RMT 셔틀(Angiopep-2, ApoE)을 유연·강직 링커와 조합해 폭넓게 "
                   "스크리닝하고, 유망 후보를 design_candidate로 잔기 수준 정밀 최적화하겠습니다."),
        AgentEvent("evaluation", "RMT·CPP 셔틀 초기 스크리닝(링커 유무 포함)",
                   data={"rows": [r1, r2, r3]}),
        AgentEvent("progress", data={"round": 1, "best_bbb": 0.96}),
        AgentEvent("structure", "셔틀 노출도=0.69 (>=0.35 = exposed), 셔틀 pLDDT 55/전체 62. 판정: exposed",
                   data={"exposed": True}),
        AgentEvent("generation", "novel 셔틀 생성",
                   data={"novel": [{"shuttle": "RYRVK", "bbb": 0.77, "tox": 0.19},
                                   {"shuttle": "VYRKKR", "bbb": 0.77, "tox": 0.25}]}),
        AgentEvent("progress", data={"round": 2, "best_bbb": 0.96}),
        AgentEvent("choice", data=choice),
        AgentEvent("final",
                   "## 결론\n\n**최종 융합체:** `GSNKGAIIGLMEAAAKEAAAKTFFYGGSRGKRNNFKTEEY`\n\n"
                   "- BBB 투과 점수 **0.42** · 독성 **0.24** · 선택성 **0.91(RMT)** · 용해도 **높음**\n\n"
                   "Angiopep-2는 LRP1 수용체 선택적이라 off-target 위험이 낮아 최종 후보로 선정합니다."),
        AgentEvent("critique",
                   "제안 후보는 RMT 수용체 선택성(0.91)이 높아 off-target 위험이 낮고 독성·용해도도 "
                   "양호합니다. 다만 불안정성 지수 45.3으로 <40 기준을 살짝 초과하니 링커 길이 조정 "
                   "여지가 있습니다.\n\nVERDICT: APPROVE", data={"approve": True}),
        AgentEvent("optimum", data=choice),
    ]
    return {"cargo": cargo, "events": events, "rounds": 8,
            "brain": "목업(디자인 미리보기)", "label": "🎨 디자인 미리보기 (목업)"}
