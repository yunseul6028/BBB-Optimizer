"""수용체/수송 메커니즘 프록시 — 검증된 BBB 셔틀과의 유사도.

deepB3P 분류 확률보다 **메커니즘에 근거한** 지표를 제공한다. 진짜 투과의 병목은
셔틀이 수용체에 결합하거나(RMT, 예: Angiopep-2 → LRP1) 막을 투과(CPP)하는 것.
후보 셔틀이 **검증된 BBB 셔틀과 얼마나 닮았는지**를 서열 정렬(BLOSUM62)로 재서,
"어떤 알려진 수송 메커니즘에 해당할 가능성"을 근사한다.

⚠️ 실제 결합 친화도·투과계수가 아니다. 서열 유사도 기반 **기능 유추(프록시)** —
   Angiopep과 전혀 다른 새 결합자는 낮게 나올 수 있다. 실험 검증 필요.
"""

from __future__ import annotations

from dataclasses import dataclass

# 검증된 BBB 셔틀 참조 세트 (메커니즘·타깃 라벨 포함)
#
# 🔴 **이 명부는 `config.SHUTTLES` 와 어긋나면 안 된다** (§32). THR·Tet1 이 여기 없어서
#    `rmt_similarity` 가 **0.0** 이 됐고, 그 결과 ⑴ off 가 2.5배 과대평가되어 THR 이 전선에서
#    배제됐고 ⑵ `rmt_floor = 0.85` 를 **시드 자신이 위반**해 `pareto_scrub` 이 두 셔틀에 대해
#    **조용히 아무것도 하지 않았다**(전선 = 시드 1개). 길이 불일치는 원인이 아니라 증폭기였다 —
#    T7(7aa)·B6(9aa)도 짧지만 **자기 자신이 명부에 있어** 1.0 을 받는다. 즉 축이 재고 있던 것은
#    유사도가 아니라 **명부 등재 여부**였다. 아래 `assert_roster_covers_library()` 로 고정한다.
#
# ⚠️ **Tet1 은 의도적으로 넣지 않는다.** `config.SHUTTLES` 의 Tet1 주석과 `evidence.py` 가
#    적듯 Tet1 은 **근거 등급 D** — GT1b 결합 뉴런 표적 리간드이고 **BBB 통과 단독 근거가 없다**
#    (문헌에서 거의 항상 B6·TGN 과 병용). 이 명부의 이름이 "**검증된** 셔틀 참조 세트"이므로
#    등급 D 서열을 넣으면 명부 자체의 근거가 무너지고, RMT/CPP 분류의 기준점이 오염된다.
#    그 대가로 Tet1 은 계속 rmt_sim = 0.0 이고 off 가 과대평가된다 — 이는 **버그가 아니라
#    받아들인 비용**이다. 다만 그로 인한 `pareto_scrub` 무작동은 조용해선 안 된다(§32-4 처방 ⑶).
REFERENCE_SHUTTLES = {
    "Angiopep-2": {"seq": "TFFYGGSRGKRNNFKTEEY", "mech": "수용체매개 수송(RMT)", "target": "LRP1"},
    "Angiopep-1": {"seq": "TFFYGGCRGKRNNFKTEEY", "mech": "수용체매개 수송(RMT)", "target": "LRP1"},
    "ApoE(159-167)2": {"seq": "LAVYQAGARLAVYQAGAR", "mech": "수용체매개 수송(RMT)", "target": "LDLR/LRP1"},
    "Leptin30": {"seq": "YQQILTSMPSRNVIQISNDLENLRDLLHVL", "mech": "수용체매개 수송(RMT)", "target": "렙틴수용체(LepR)"},
    "T7":       {"seq": "HAIYPRH",             "mech": "수용체매개 수송(RMT)", "target": "TfR(트랜스페린 수용체)"},
    "B6":       {"seq": "GHKAKGPRK",           "mech": "수용체매개 수송(RMT)", "target": "TfR1(트랜스페린 수용체)"},
    "RVG29":    {"seq": "YTIWMPENPRPGTPCDIFTNSRGKRASNG", "mech": "수용체매개 수송(RMT)", "target": "nAChR(니코틴성 수용체)"},
    "CDX":      {"seq": "FKESWREARGTRIERG",    "mech": "수용체매개 수송(RMT)", "target": "α7 nAChR(니코틴성 수용체)"},
    "THR":      {"seq": "THRPPMWSPVWP",         "mech": "수용체매개 수송(RMT)", "target": "TfR1(트랜스페린 수용체)"},
    "TAT":        {"seq": "GRKKRRQRRRPPQ",       "mech": "세포투과(CPP)",       "target": "막 직접투과"},
    "Penetratin": {"seq": "RQIKIWFQNRRMKWKK",    "mech": "세포투과(CPP)",       "target": "막 직접투과"},
    "SynB1":      {"seq": "RGGRLSYSRRRFSTSTGR",  "mech": "세포투과(CPP)",       "target": "막 직접투과"},
}

_ALIGNER = None


def _aligner():
    global _ALIGNER
    if _ALIGNER is None:
        from Bio.Align import PairwiseAligner, substitution_matrices
        a = PairwiseAligner()
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score = -10.0
        a.extend_gap_score = -0.5
        a.mode = "global"
        _ALIGNER = a
    return _ALIGNER


@dataclass
class BindingResult:
    score: float = 0.0        # 0~1, 검증 셔틀과의 최대 정규화 유사도
    best_ref: str = ""        # 가장 닮은 참조 셔틀
    mechanism: str = ""       # 그 셔틀의 수송 메커니즘
    target: str = ""          # 타깃(LRP1 등)
    verdict: str = ""
    per_ref: list = None      # [{name, score}]


def shuttle_similarity(shuttle: str) -> BindingResult:
    seq = "".join(ch for ch in (shuttle or "").upper() if ch.isalpha())
    if not seq:
        return BindingResult(verdict="빈 셔틀")
    try:
        aln = _aligner()
    except Exception as exc:  # noqa: BLE001
        return BindingResult(verdict=f"정렬기 오류: {exc}")

    per, best = [], None
    for name, r in REFERENCE_SHUTTLES.items():
        ref = r["seq"]
        try:
            raw = aln.score(seq, ref)
            selfs = aln.score(ref, ref)
        except Exception:  # noqa: BLE001
            continue
        norm = max(0.0, min(1.0, raw / selfs)) if selfs > 0 else 0.0
        per.append({"name": name, "score": round(norm, 3)})
        if best is None or norm > best.score:
            best = BindingResult(score=round(norm, 3), best_ref=name,
                                 mechanism=r["mech"], target=r["target"])
    if best is None:
        return BindingResult(verdict="계산 실패")
    best.per_ref = sorted(per, key=lambda x: x["score"], reverse=True)

    if best.score >= 0.6:
        best.verdict = (f"'{best.best_ref}'와 강한 유사 → {best.target} {best.mechanism} "
                        "메커니즘 가능성 높음")
    elif best.score >= 0.35:
        best.verdict = f"'{best.best_ref}'와 부분 유사 → {best.mechanism} 가능성 보통"
    else:
        # 🔴 **거짓 귀속 차단** (§32-3). 위 루프는 `best is None or norm > best.score` 로
        #    갱신하므로 **첫 참조가 항상 best 를 선점**한다. 점수가 전부 0.0 이어도
        #    `best_ref='Angiopep-2'` · `target='LRP1'` · `mech='RMT'` 가 붙어서,
        #    **어떤 쓰레기 서열도 "Angiopep-2 유사 · LRP1 · RMT" 로 보고됐다**
        #    (`shuttle_similarity("PPPPPPPPPPPP")` 로 재현된다). `verdict` 문자열만
        #    "메커니즘 불확실"로 올바르게 말하고 **행에 실리는 라벨은 틀렸다.**
        #    `delivery.classify_mechanism` 은 **점수를 보지 않고 문자열만** 보므로
        #    `is_rmt=True` 까지 전파됐다(`selectivity.best_is_rmt` 는 score ≥ 0.4 로 막혀 있었다).
        #    이 분기에서는 verdict 가 참조를 언급하지 않으므로 라벨을 비우는 것이 정합적이다.
        best.best_ref, best.mechanism, best.target = "", "", ""
        best.verdict = "알려진 BBB 셔틀과 약한 유사 → 메커니즘 불확실(신규 가능성)"
    return best


def assert_roster_covers_library() -> None:
    """참조 명부가 라이브러리 셔틀을 덮는지 단언 — 드리프트 조기 검출 (§32).

    명부 누락은 **조용히** off 를 과대평가하고 `pareto_scrub` 을 무력화한다. 그래서
    "조용한 퇴화는 큰 소리 실패로"(§7 계열) 원칙을 명부에도 적용한다.
    **의도된 예외는 이유와 함께 여기 적는다** — 목록에 없는 누락만 실패시킨다.
    """
    from .config import SHUTTLES
    intentional = {
        # 근거 등급 D — BBB 통과 단독 근거 없음. 위 REFERENCE_SHUTTLES 주석 참조.
        "Tet1": "근거 등급 D(GT1b 뉴런 표적, BBB 단독 근거 없음) — 검증 명부에 넣지 않는다",
    }
    missing = {n for n in SHUTTLES if n not in REFERENCE_SHUTTLES} - set(intentional)
    if missing:
        raise RuntimeError(
            f"참조 명부 누락: {sorted(missing)} — `core/binding.REFERENCE_SHUTTLES` 에 없다. "
            "그대로 두면 rmt_similarity 가 0.0 이 되어 off 가 과대평가되고 "
            "rmt_floor 제약이 시드 자신을 거부해 개발성 최적화가 조용히 무작동한다(§32). "
            "명부에 추가하거나, 제외가 의도라면 `assert_roster_covers_library` 의 "
            "`intentional` 에 **이유와 함께** 등록하라."
        )
