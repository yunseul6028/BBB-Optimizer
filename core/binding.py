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
REFERENCE_SHUTTLES = {
    "Angiopep-2": {"seq": "TFFYGGSRGKRNNFKTEEY", "mech": "수용체매개 수송(RMT)", "target": "LRP1"},
    "Angiopep-1": {"seq": "TFFYGGCRGKRNNFKTEEY", "mech": "수용체매개 수송(RMT)", "target": "LRP1"},
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
        best.verdict = "알려진 BBB 셔틀과 약한 유사 → 메커니즘 불확실(신규 가능성)"
    return best
