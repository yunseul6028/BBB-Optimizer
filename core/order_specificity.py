"""순서 특이성(스크램블 델타) — deepB3P 점수가 '순서'에서 왔나 '조성'에서 왔나.

## 왜 필요한가

M1 실측: 양성 vs 무작위 AUC 0.807 인데 양성 vs **스크램블**(조성 동일·순서 파괴)은
0.618 뿐이다. 즉 deepB3P 점수의 상당 부분이 **아미노산 조성**에서 온다. 그래서
아르기닌을 쌓으면 점수가 오른다 — 보상 해킹의 메커니즘이 이것이다.

대응: 서열의 **자기 스크램블을 대조군**으로 같이 채점하고

    delta = bbb(서열) - 평균 bbb(스크램블 n개)

을 본다. 조성 쌓기는 섞어도 점수가 그대로라 delta ≈ 0 이 된다.

## 실측 교정 결과 (scripts/scramble_delta.py, 화물 GSNKGAIIGLM, n=20)

**셔틀 단독** — 분리 명확:

| 그룹 | 델타 중위 | 범위 |
|---|---|---|
| 후보군 RMT 9종 | +0.289 | −0.020 ~ +0.788 |
| 비교군 CPP 3종 | +0.064 | +0.001 ~ +0.756 |
| 조성 쌓기 프로브 | +0.000 | +0.000 ~ +0.018 |

**융합체** — 분리 실패(RMT 중위 +0.033, 9종 중 4종 음수, CPP 중위가 RMT보다 높음).
접합부에서 효능이 무너지는 셔틀은 원본 점수가 0 근처라 델타가 양수일 수 없다.
→ **융합 수준 게이트는 폐기**했다. 셔틀 단독으로만 본다.

## 왜 라이브러리 셔틀에는 게이트를 걸지 않나

프로브를 잡는 어떤 마진(≥ +0.02)도 **THR(+0.014)을 같이 잡는다**. THR 의 1.000 점도
조성이라는 뜻이지만, THR 은 문헌 보고된 TfR1 파지디스플레이 12-mer 다. 즉 그 서열의
근거는 deepB3P 가 아니라 **문헌**에 있다. 모델의 조성 성분을 못 믿는다는 이유로
published 셔틀을 떨어뜨리는 건 앞뒤가 안 맞는다.

그래서 선을 이렇게 긋는다:
- **라이브러리 셔틀** = 근거가 외부(문헌)에 있음 → 게이트 면제, 지표로만 보고
- **설계체**(`design_candidate` 산출물) = 근거가 모델뿐 → **하드 게이트**(마진 `DELTA_MIN`)

## 측정 불가를 0 으로 때우지 말 것

동일잔기 반복서열(`RRRRRRRRR`)은 섞어도 자기 자신이라 유효 스크램블이 **0개**다.
이때 평균을 0.0 으로 기본값 처리하면 `delta = 1.000 - 0 = +1.000` 이라는 **가짜 최고점**이
나온다(실제로 한 번 밟았다 — poly-Arg 가 전체 1위로 찍혔다). 측정 불가는 **None**이다.
"""
from __future__ import annotations

import random
import statistics

# 설계체 하드 게이트 마진. 교정 근거: 조성 쌓기 프로브 최대 +0.018,
# 정당한 RMT 셔틀 8/9 가 +0.014 이상(중위 +0.289). 프로브 상한 위로 여유를 둔 값.
DELTA_MIN = 0.05

# 설계체당 스크램블 개수. 20개로 교정했으나 루프 안에서는 예측 비용이 붙으므로 12개로
# 줄였다 — 평균 추정에는 충분하고(프로브들의 sd ≤ 0.006), 호출당 추가 예측 12회다.
N_SCRAMBLE = 12


def scrambles(seq: str, n: int = N_SCRAMBLE, seed: int = 0) -> list[str]:
    """조성 보존·순서 파괴 변이를 최대 n개. **원본과 같은 것은 제외**.

    동일잔기 반복서열이면 빈 리스트를 돌려준다(측정 불가 — 위 모듈 docstring 참고).
    """
    rng = random.Random(f"{seq}:{seed}")
    out: set[str] = set()
    tries = 0
    while len(out) < n and tries < n * 50:
        tries += 1
        chars = list(seq)
        rng.shuffle(chars)
        s = "".join(chars)
        if s != seq:
            out.add(s)
    return list(out)


def order_delta(seq: str, predictor, n: int = N_SCRAMBLE, seed: int = 0):
    """(delta, real_bbb, scramble_mean, n_used). 측정 불가면 delta=None.

    `seq` 는 **셔틀 단독** 서열이다(융합체가 아니다 — 융합 수준은 분리가 안 된다).
    예측 실패 시에도 delta=None 을 돌려준다. 절대 0.0 으로 때우지 않는다.
    """
    seq = "".join(ch for ch in (seq or "").upper() if ch.isalpha())
    if len(seq) < 4:
        return None, None, None, 0
    scr = scrambles(seq, n, seed)
    if not scr:
        # 동일잔기 반복서열 → 스크램블 대조가 존재하지 않는다.
        try:
            real = predictor.predict_many([seq])[0].bbb_permeability
        except Exception:  # noqa: BLE001
            real = None
        return None, real, None, 0
    try:
        preds = [p.bbb_permeability for p in predictor.predict_many([seq] + scr)]
    except Exception:  # noqa: BLE001
        return None, None, None, 0
    real, sv = preds[0], preds[1:]
    mean = statistics.fmean(sv)
    return real - mean, real, mean, len(sv)


def verdict(delta, is_library: bool) -> tuple[bool, str]:
    """(통과 여부, 사람이 읽을 한 줄). 라이브러리 셔틀은 항상 통과(근거가 문헌에 있음)."""
    if is_library:
        d = "측정 불가" if delta is None else f"{delta:+.3f}"
        return True, f"순서 특이성 델타 {d} (라이브러리 셔틀 — 근거가 문헌에 있어 게이트 면제)"
    if delta is None:
        # 설계체인데 스크램블 대조를 만들 수 없다 = 동일잔기 반복서열.
        # 그 자체가 조성 쌓기의 극단형이므로 통과시키지 않는다.
        return False, ("순서 특이성 **측정 불가** — 동일잔기 반복서열이라 스크램블 대조가 "
                       "존재하지 않는다. 이는 조성 쌓기의 극단형이다.")
    if delta < DELTA_MIN:
        return False, (f"순서 특이성 델타 {delta:+.3f} < 요구 {DELTA_MIN:+.3f} — 이 설계체의 "
                       f"BBB 점수는 **서열 순서가 아니라 아미노산 조성**에서 왔다. 같은 조성을 "
                       f"무작위로 섞어도 점수가 거의 같다는 뜻이다(조성 쌓기).")
    return True, f"순서 특이성 델타 {delta:+.3f} ≥ {DELTA_MIN:+.3f} (순서가 점수를 만든다)"


# ──────────────────────────────────────────────────────────────────────────────
# 길이 대응 널(matched null) — '접합부에 잔기를 끼워 점수를 올리는' 해킹 차단
# ──────────────────────────────────────────────────────────────────────────────
# 위 스크램블 델타는 **순서**의 대응 널이다. 그런데 에이전트는 순서를 건드리지 않고
# **길이**로 같은 짓을 할 수 있다. 실측(docs/experiments.md §13):
#
#   화물 GSNKGAIIGLM + T7(HAIYPRH) 융합 BBB = 0.526
#   앞에 **무작위 1잔기**를 끼운 20종: 평균 0.743 · sd 0.116 · 범위 0.471~0.933
#   → 끼우는 것 자체로 평균 **+0.217 이 공짜로** 오른다
#
# 어느 잔기를 고르느냐는 그 노이즈 밴드 안의 운이다: 에이전트가 고른 Y 는 z=+1.14,
# 가장 높은 R 도 z=+1.64 로 **20종 중 z≥2 는 0개**다. 발린이나 타이로신 하나가 TfR
# 결합이나 transcytosis 를 바꿀 메커니즘도 없다. 즉 이 축에는 설계 신호가 없다.
#
# 중요: 링커 필드를 막는 것으로는 **못 막는다**. `링커 Y + 셔틀 HAIYPRH` 와
# `직접융합 + 셔틀 YHAIYPRH` 는 **같은 분자**이고 점수도 똑같이 0.875 다(실측).
# 링커/셔틀 구분은 장부상의 경계일 뿐 화학적 경계가 아니라서, 필드를 막으면
# 에이전트가 옆 필드로 옮긴다(실제로 그렇게 우회했다). 그래서 통제는 **구성체(construct)
# 수준**에서, 부모 셔틀 대비 **같은 길이만큼 무작위로 늘린 널**과 비교해 걸어야 한다.
Z_MIN = 2.0          # 널 대비 요구 z. 2.0 이면 1잔기 삽입 클래스 전체가 막힌다(실측)
N_NULL = 20          # 널 표본 수

_AA = "ACDEFGHIKLMNPQRSTVWY"


def find_parent(shuttle: str, library: set[str] | frozenset[str]):
    """설계 셔틀이 '라이브러리 셔틀 + 양끝 덧붙임'이면 (부모, 앞덧길이, 뒤덧길이).

    없으면 None. 가장 긴 부모를 고른다(가장 보수적인 해석 = 덧붙임이 가장 적은 쪽).
    """
    sh = (shuttle or "").upper()
    best = None
    for par in library:
        if not par or par == sh or par not in sh:
            continue
        i = sh.index(par)
        cand = (par, i, len(sh) - i - len(par))
        if best is None or len(par) > len(best[0]):
            best = cand
    return best


def insertion_z(cargo: str, shuttle: str, predictor, library,
                n: int = N_NULL, seed: int = 0):
    """부모 대비 **길이만 늘린** 설계인지 검사 → (z, 부모, 덧붙인 길이, 널평균, 실제값).

    부모가 없거나(완전 신규 설계) 길이가 안 늘었으면(순수 치환) z=None 을 돌려준다
    — 그 경우는 이 통제의 대상이 아니다(순서 델타가 담당한다).
    실패 시에도 None. 절대 0 으로 때우지 않는다.
    """
    import random as _random
    import statistics as _st
    from .config import bbb_scoring_seq

    sh = "".join(ch for ch in (shuttle or "").upper() if ch.isalpha())
    found = find_parent(sh, library)
    if found is None:
        return None, None, 0, None, None
    parent, npre, npost = found
    k = npre + npost
    if k <= 0:
        return None, parent, 0, None, None

    rng = _random.Random(f"{cargo}:{sh}:{seed}")
    nulls = []
    for _ in range(n):
        pre = "".join(rng.choice(_AA) for _ in range(npre))
        post = "".join(rng.choice(_AA) for _ in range(npost))
        nulls.append(pre + parent + post)
    try:
        preds = [p.bbb_permeability for p in predictor.predict_many(
            [bbb_scoring_seq(cargo, "", s) for s in [sh] + nulls])]
    except Exception:  # noqa: BLE001
        return None, parent, k, None, None
    real, nv = preds[0], preds[1:]
    mean = _st.fmean(nv)
    sd = _st.pstdev(nv)
    if sd < 1e-6:
        # 널이 전혀 흩어지지 않으면 z 가 정의되지 않는다 → 측정 불가.
        return None, parent, k, mean, real
    return (real - mean) / sd, parent, k, mean, real


def insertion_verdict(z, parent, k, mean, real) -> tuple[bool, str]:
    """(통과 여부, 설명). z 가 None 이면 이 통제의 대상이 아니므로 통과."""
    if z is None:
        return True, ""
    if z >= Z_MIN:
        return True, (f"접합부 길이 통제 통과 — 부모 `{parent}` 대비 {k}잔기 추가, "
                      f"같은 길이 무작위 널 평균 {mean:.3f} 대비 z {z:+.2f} ≥ {Z_MIN}")
    return False, (f"**접합부 길이 통제 미달** — 부모 `{parent}` 에 {k}잔기를 덧붙여 "
                   f"BBB {real:.3f} 를 얻었지만, 같은 길이를 **무작위로** 덧붙인 널의 평균이 "
                   f"{mean:.3f} 다(z {z:+.2f} < {Z_MIN}). 즉 이 상승은 어떤 잔기를 골랐기 "
                   f"때문이 아니라 **길이를 늘렸기 때문**이다. deepB3P 는 접합부 길이에 "
                   f"민감하지만, 잔기 하나가 수용체 결합이나 transcytosis 를 바꿀 메커니즘은 "
                   f"없다. 길이를 늘리는 대신 셔틀 **내부 잔기를 치환**하라.")
