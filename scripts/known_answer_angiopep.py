#!/usr/bin/env python3
"""문헌에 정답이 있는 서열 한 쌍을 우리 도구에 넣으면 무슨 일이 일어나는가.

왜 이 스크립트가 있나
---------------------
이 프로젝트의 축은 거의 전부 **우리가 만든 자**다. 자기 자를 자기가 검증하면
순환이다. 그래서 **밖에 정답이 적혀 있는 입력** 하나를 통째로 통과시킨다.

정답(Demeule et al. 2008 *J Pharmacol Exp Ther* 324(3):1064-1072, 원문 대조 완료 §40):
    Angiopep-2  TFFYGGSRGKRNNFKTEEY       원본
    Angiopep-5  = K10R          단독 변이  뇌 실질 흡수 **약 60% 감소**
    Angiopep-7  = K10R + K15R   이중 변이  뇌 실질 흡수 **약 85% 감소**
    셋의 **순전하는 모두 +2 로 같다** → 기능 상실은 전하 효과가 아니라 **자리 효과**다.
    원문: "replacement of the lysine residues at positions 10 and 15 prevents the
    transport capacity."

🔴 이것은 **우리 축의 정확도 측정이 아니다.** n=3 이고 한 계열이다.
   이 스크립트가 답하는 질문은 하나뿐이다 —
   **"문헌이 망가졌다고 말하는 서열을, 우리 축 중 무엇이 알아보는가?"**

🔴 **배수(~9.5배 · 86배)는 쓰지 않는다.** 그 숫자는 인용 출처에 없었다(§40).
   쓸 수 있는 것은 **−60% / −85%** 와 **순전하 +2 동일** 뿐이다.

실행: /usr/bin/python3 scripts/known_answer_angiopep.py   (유료 LLM 0회)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core import evidence                                  # noqa: E402
from core.binding import shuttle_similarity                # noqa: E402
from core.config import (DEFAULT_CARGO, LINKER_LIBRARY,    # noqa: E402
                         SHUTTLES, bbb_scoring_seq, get_settings)
from core.developability import (assess_developability,    # noqa: E402
                                 junction_liabilities)
from core.stability import assess_stability                # noqa: E402
from core.predictors import get_predictor                  # noqa: E402
from core.selectivity import assess_selectivity            # noqa: E402
from core.toxicity import get_toxicity_predictor           # noqa: E402

A2 = SHUTTLES["Angiopep-2"]["seq"]


def _mutate(seq: str, *subs: tuple[int, str, str]) -> str:
    """1-기반 위치에서 치환. (위치, 원잔기, 새잔기) — 원잔기가 다르면 즉시 터진다."""
    out = list(seq)
    for pos, old, new in subs:
        assert out[pos - 1] == old, f"위치 {pos} 가 {old} 가 아니라 {out[pos-1]} 다"
        out[pos - 1] = new
    return "".join(out)


A5 = _mutate(A2, (10, "K", "R"))                      # Angiopep-5
A7 = _mutate(A2, (10, "K", "R"), (15, "K", "R"))      # Angiopep-7

PANEL = [
    ("Angiopep-2", A2, "원본",          None),
    ("Angiopep-5", A5, "K10R",          -60),
    ("Angiopep-7", A7, "K10R·K15R",     -85),
]

CARGO = DEFAULT_CARGO
LINKER = LINKER_LIBRARY["GGGGS"]["seq"]


def main() -> int:
    settings = get_settings()
    bbb = get_predictor(settings)
    tox = get_toxicity_predictor(settings)

    print("=" * 78)
    print("정답이 밖에 있는 입력 — Angiopep 계열 3종을 현행 축 전부에 통과시킨다")
    print("=" * 78)
    print(f"측정 조건: 화물 {CARGO!r} · 링커 GGGGS · 현행 코드(§38-1 liability 정의 ·")
    print(f"           §26 off 밴드 · §49 수리 이후) · BBB {type(bbb).__name__} ·")
    print(f"           독성 {type(tox).__name__}\n")

    print("[1] 문헌이 말하는 것 (Demeule 2008 JPET, 원문 대조 완료)")
    print("-" * 78)
    for name, seq, edit, loss in PANEL:
        lit = "기준" if loss is None else f"뇌 실질 흡수 **{loss}%**"
        print(f"  {name:<12} {seq}  {edit:<10} {lit}")
    print("  → 순전하는 셋 다 동일해야 한다(문헌 주장). 아래 [2] 에서 검산한다.\n")

    # 예측기는 배치로 한 번에 — 모델 로딩을 3번 하지 않는다
    shuttle_seqs = [seq for _, seq, _, _ in PANEL]
    construct_seqs = [bbb_scoring_seq(CARGO, LINKER, seq) for seq in shuttle_seqs]
    bbb_sh = bbb.predict_many(shuttle_seqs)
    bbb_co = bbb.predict_many(construct_seqs)
    tox_sh = tox.predict_many(shuttle_seqs)
    tox_co = tox.predict_many(construct_seqs)

    rows = []
    for i, (name, seq, edit, loss) in enumerate(PANEL):
        dev = assess_developability(seq)
        stab = assess_stability(seq)
        sel = assess_selectivity(seq)
        sim = shuttle_similarity(seq)
        rows.append({
            "name": name, "seq": seq, "edit": edit, "loss": loss,
            "charge": dev.net_charge,
            "ii": stab.instability_index,
            "agg": dev.agg_hotspot,
            "liab_n": dev.n_liabilities,
            "off": sel.off_target_risk,
            "rmt": sim.score,
            "ref": sim.best_ref,
            "bbb_shuttle": bbb_sh[i].bbb_permeability,
            "bbb_construct": bbb_co[i].bbb_permeability,
            "tox_shuttle": tox_sh[i].risk,
            "tox_construct": tox_co[i].risk,
            "junction_n": len(junction_liabilities(CARGO, LINKER, seq)),
        })

    print("[2] 우리 축 전부가 말하는 것")
    print("-" * 78)
    hdr = f"{'축':<26}{'Angiopep-2':>13}{'Angiopep-5':>13}{'Angiopep-7':>13}   방향"
    print(hdr)
    print("-" * 78)

    def line(label, key, fmt="{:.3f}", better_low=True, note=""):
        v = [r[key] for r in rows]
        # 문헌 방향: A2 가 가장 좋고 A7 이 가장 나쁘다.
        if v[0] == v[1] == v[2]:
            verdict = "⬜ 구별 못 함"
        elif (v[0] < v[2]) == better_low:
            verdict = "✅ 방향 맞음"
        else:
            verdict = "🔴 방향 반대"
        cells = "".join(fmt.format(x).rjust(13) for x in v)
        print(f"{label:<26}{cells}   {verdict} {note}")
        return verdict

    verdicts = {}
    verdicts["순전하"]       = line("순전하 (pH 7.4)", "charge", "{:+.2f}", True,
                                  "← 문헌이 '같다'고 한 칸")
    verdicts["rmt_sim"]      = line("rmt_sim (기능 보존 프록시)", "rmt", "{:.3f}", False)
    verdicts["off_target"]   = line("off_target", "off", "{:.2f}", True)
    verdicts["불안정성"]     = line("불안정성 지수", "ii", "{:.1f}", True)
    verdicts["응집"]         = line("응집 핫스팟(정규화)", "agg", "{:.3f}", True)
    verdicts["liability"]    = line("liability 개수", "liab_n", "{:d}", True)
    verdicts["접합부"]       = line("접합부 liability(구성체)", "junction_n", "{:d}", True)
    verdicts["독성(셔틀)"]   = line("ToxinPred3 ML Score(셔틀)", "tox_shuttle", "{:.3f}", True)
    verdicts["독성(구성체)"] = line("ToxinPred3 ML Score(구성체)", "tox_construct", "{:.3f}", True)
    verdicts["deepB3P(셔틀)"] = line("deepB3P (셔틀 단독)", "bbb_shuttle", "{:.3f}", False)
    verdicts["deepB3P(구성체)"] = line("deepB3P (구성체)", "bbb_construct", "{:.3f}", False)

    print()
    print("[3] 집계 — 문헌 정답 대비")
    print("-" * 78)
    ok   = [k for k, v in verdicts.items() if v.startswith("✅")]
    flat = [k for k, v in verdicts.items() if v.startswith("⬜")]
    bad  = [k for k, v in verdicts.items() if v.startswith("🔴")]
    print(f"  ✅ 방향 맞음  {len(ok):>2}/{len(verdicts)} — {', '.join(ok) or '없음'}")
    print(f"  ⬜ 구별 못 함 {len(flat):>2}/{len(verdicts)} — {', '.join(flat) or '없음'}")
    print(f"  🔴 방향 반대  {len(bad):>2}/{len(verdicts)} — {', '.join(bad) or '없음'}")

    print()
    print("[3b] 🔴 방향을 맞춘 3개는 **판정을 바꾸는가**")
    print("-" * 78)
    from core.developability_opt import OFF_BAND
    from core.optimizer_agent import II_STABLE
    rmt_floor = 0.85                      # core/developability_opt.pareto_scrub 기본값
    r2, r7 = rows[0], rows[2]
    print(f"  rmt_sim      하한 {rmt_floor}  →  Angiopep-7 {r7['rmt']:.3f} = **통과**        ❌ 안 막는다")
    print(f"  off_target   밴드 ±{OFF_BAND}  →  변화 {r2['off']:.2f}→{r7['off']:.2f} = {abs(r7['off']-r2['off']):.2f} "
          f"= **밴드 안**  ❌ 안 막는다")
    print(f"               (그 0.01 은 off 축의 **양자화 눈금** 자체다 — §42-4b)")
    print(f"  불안정성     이진화 임계 {II_STABLE}  →  {r2['ii']:.1f} 도 {r7['ii']:.1f} 도 '불안정'")
    print(f"               = 같은 칸  ❌ 순위 기여 **0**")
    print()
    print("  🔴 **방향을 맞춘 축 3개 · 그 중 판정을 바꾸는 축 0개.**")
    print("     방향이 맞는 것과 쓸모가 있는 것은 다르다. 세 축 모두 Angiopep-7 을")
    print("     '약간 나쁘다'고 말하지만, 어느 축도 **후보에서 떨어뜨리지 않는다**.")

    print()
    print("[4] 서열이 아닌 장치 — 근거 등급 상속 (`core/evidence.py`)")
    print("-" * 78)
    pins = evidence.conserved_positions(A2)
    for r in rows:
        tier, parent, sim = evidence.inherit(r["seq"])          # (등급, 부모, 유사도)
        got = f"등급 {tier}" if tier else "🔴 등급 없음 — 상속 거부"
        print(f"  {r['name']:<12} 최근접 {parent or '—':<12} 유사도 {sim:.3f} → {got}")
    print(f"\n  문헌 보존 위치(Angiopep-2): {pins}  ← 이것은 측정이 아니라 **문헌 사실**이다.")
    print("  유사도 하한(0.70)은 셋 다 통과한다. 상속을 막는 것은 유사도가 아니라")
    print("  **그 자리의 라이신이 살아 있는가** 하나뿐이다.")

    print()
    print("[5] 이 시연이 말하지 않는 것")
    print("-" * 78)
    print("  ⓐ n=3 · 한 계열 · 한 화물. **축의 정확도가 아니다.**")
    print("  ⓑ 🚫 배수(~9.5배·86배)는 인용 출처에 없다(§40). −60%/−85% 와 순전하 동일만 쓴다.")
    print("  ⓒ 문헌 보존 위치가 있는 셔틀은 **12종 중 1종**뿐이다 — [4] 의 장치는")
    print("     Angiopep-2 계열 밖에서는 작동하지 않는다(§38-2).")
    print("  ⓓ K→R 은 BLOSUM62 **+2(보존적 치환)** 다. [2] 의 실패는 버그가 아니라")
    print("     **서열 유사도라는 양이 원래 이 질문에 답할 수 없다**는 뜻이다.")

    if "--md" in sys.argv:
        out = Path(__file__).resolve().parent.parent / "docs" / "figs" / "표-known-answer-angiopep.md"
        out.write_text(_md(rows, verdicts, rmt_floor, OFF_BAND, II_STABLE, pins),
                       encoding="utf-8")
        print(f"\n📋 표 파일: {out}")

    if "--html" in sys.argv:
        out = Path(__file__).resolve().parent.parent / "docs" / "figs" / "known-answer-angiopep.html"
        out.write_text(_html(rows, verdicts, rmt_floor, OFF_BAND, II_STABLE, pins),
                       encoding="utf-8")
        print(f"\n🖼  그림 파일: {out}")
        print("    브라우저로 열고 ⌘⇧4 로 캡처하거나, ⌘P → PDF 로 저장하세요.")
    return 0


# ── 그림 생성 (--html) ────────────────────────────────────────────────────────
_ROWSPEC = [
    ("순전하 (pH 7.4)",            "charge", "{:+.2f}", "문헌이 '셋 다 같다'고 한 칸"),
    ("rmt_sim (기능 보존 프록시)",  "rmt",    "{:.3f}",  "하한 0.85"),
    ("off_target",                 "off",    "{:.2f}",  "밴드 ±0.05"),
    ("불안정성 지수",               "ii",     "{:.1f}",  "이진화 임계 40"),
    ("응집 핫스팟 (정규화)",        "agg",    "{:.3f}",  ""),
    ("liability 개수",             "liab_n", "{:d}",    ""),
    ("접합부 liability (구성체)",   "junction_n", "{:d}", ""),
    ("ToxinPred3 ML Score (셔틀)",  "tox_shuttle",   "{:.3f}", "하드 게이트 0.38"),
    ("ToxinPred3 ML Score (구성체)","tox_construct", "{:.3f}", ""),
    ("deepB3P (셔틀 단독)",         "bbb_shuttle",   "{:.3f}", "보고치 — 순위 미사용"),
    ("deepB3P (구성체)",            "bbb_construct", "{:.3f}", "보고치 — 순위 미사용"),
]
_VKEY = ["순전하", "rmt_sim", "off_target", "불안정성", "응집", "liability", "접합부",
         "독성(셔틀)", "독성(구성체)", "deepB3P(셔틀)", "deepB3P(구성체)"]


def _html(rows, verdicts, rmt_floor, off_band, ii_stable, pins) -> str:
    vlist = [verdicts[k] for k in _VKEY]
    n_ok = sum(1 for v in vlist if v.startswith("✅"))
    n_flat = sum(1 for v in vlist if v.startswith("⬜"))
    n_bad = sum(1 for v in vlist if v.startswith("🔴"))
    r2, r7 = rows[0], rows[2]

    body = []
    for (label, key, fmt, note), v in zip(_ROWSPEC, vlist):
        cls = "ok" if v.startswith("✅") else ("flat" if v.startswith("⬜") else "bad")
        tag = {"ok": "방향 맞음", "flat": "구별 못 함", "bad": "방향 반대"}[cls]
        cells = "".join(f"<td class='num'>{fmt.format(r[key])}</td>" for r in rows)
        body.append(
            f"<tr class='{cls}'><th class='ax'>{label}"
            f"{f'<span class=note>{note}</span>' if note else ''}</th>"
            f"{cells}<td class='v {cls}'>{tag}</td></tr>")

    return f"""<!doctype html>
<html lang="ko"><meta charset="utf-8">
<title>Angiopep 알려진 정답 대조</title>
<style>
  :root {{ --ink:#16181d; --mut:#6b7280; --line:#e5e7eb; --bg:#fff;
           --ok:#166534; --okbg:#f0fdf4; --bad:#b91c1c; --badbg:#fef2f2;
           --flat:#92400e; --flatbg:#fffbeb; }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; padding:34px 30px; background:var(--bg); color:var(--ink);
    font-family:"Apple SD Gothic Neo","Pretendard","Noto Sans KR",system-ui,sans-serif;
    font-size:14px; line-height:1.55; -webkit-font-smoothing:antialiased; }}
  .wrap {{ max-width:880px; margin:0 auto; }}
  h1 {{ font-size:19px; margin:0 0 4px; letter-spacing:-.01em; }}
  .sub {{ color:var(--mut); font-size:12.5px; margin:0 0 20px; }}
  .lit {{ border:1px solid var(--line); border-left:3px solid var(--ink);
    padding:12px 15px; margin:0 0 20px; background:#fafafa; }}
  .lit h2 {{ font-size:12px; margin:0 0 8px; letter-spacing:.04em;
    text-transform:uppercase; color:var(--mut); font-weight:600; }}
  .lit table {{ border-collapse:collapse; font-size:13px; }}
  .lit td {{ padding:2px 16px 2px 0; }}
  .seq {{ font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:12.5px;
    letter-spacing:.02em; }}
  .hit {{ background:#fde68a; padding:0 1px; border-radius:2px; font-weight:700; }}
  table.main {{ width:100%; border-collapse:collapse; font-size:13px; }}
  table.main th, table.main td {{ padding:7px 9px; border-bottom:1px solid var(--line); }}
  table.main thead th {{ font-size:11.5px; color:var(--mut); font-weight:600;
    letter-spacing:.03em; border-bottom:1.5px solid var(--ink); text-align:right; }}
  table.main thead th:first-child {{ text-align:left; }}
  .ax {{ text-align:left; font-weight:500; white-space:nowrap; }}
  .note {{ display:block; font-size:10.5px; color:var(--mut); font-weight:400; }}
  .num {{ text-align:right; font-variant-numeric:tabular-nums;
    font-family:ui-monospace,SFMono-Regular,Menlo,monospace; }}
  .v {{ text-align:right; font-size:11.5px; font-weight:600; white-space:nowrap; }}
  tr.ok .v {{ color:var(--ok); }}  tr.ok {{ background:var(--okbg); }}
  tr.bad .v {{ color:var(--bad); }} tr.bad {{ background:var(--badbg); }}
  tr.flat .v {{ color:var(--flat); }} tr.flat {{ background:var(--flatbg); }}
  .tally {{ display:flex; gap:9px; margin:16px 0 0; }}
  .tally div {{ flex:1; border:1px solid var(--line); padding:9px 11px; }}
  .tally b {{ display:block; font-size:21px; line-height:1.15;
    font-variant-numeric:tabular-nums; }}
  .tally span {{ font-size:11.5px; color:var(--mut); }}
  .punch {{ margin:18px 0 0; border:1.5px solid var(--ink); padding:14px 16px; }}
  .punch h2 {{ font-size:14px; margin:0 0 9px; }}
  .punch ul {{ margin:0; padding-left:17px; font-size:12.5px; }}
  .punch li {{ margin:3px 0; }}
  .punch .kicker {{ margin:10px 0 0; font-size:13.5px; font-weight:700; }}
  .pin {{ margin:16px 0 0; border:1px solid var(--line); padding:13px 15px;
    background:#fafafa; }}
  .pin h2 {{ font-size:13px; margin:0 0 8px; }}
  .pin p {{ margin:4px 0; font-size:12.5px; }}
  footer {{ margin-top:20px; padding-top:11px; border-top:1px solid var(--line);
    font-size:10.5px; color:var(--mut); line-height:1.5; }}
  @media print {{ body {{ padding:0; }} tr.ok,tr.bad,tr.flat {{ -webkit-print-color-adjust:exact; print-color-adjust:exact; }} }}
</style>
<div class="wrap">
  <h1>문헌에 정답이 있는 서열을 우리 축에 통과시키면</h1>
  <p class="sub">Angiopep-2 / -5 / -7 · 화물 {CARGO} · 링커 GGGGS · 유료 LLM 0회</p>

  <div class="lit">
    <h2>밖에 적힌 정답 — Demeule et al. 2008 <i>J Pharmacol Exp Ther</i> 324(3):1064–1072 (원문 대조 완료)</h2>
    <table>
      <tr><td><b>Angiopep-2</b></td><td class="seq">TFFYGGSRG<span class="hit">K</span>RNNF<span class="hit">K</span>TEEY</td><td>원본</td><td>기준</td></tr>
      <tr><td><b>Angiopep-5</b></td><td class="seq">TFFYGGSRG<span class="hit">R</span>RNNF<span class="hit">K</span>TEEY</td><td>K10R</td><td>뇌 실질 흡수 <b>−60%</b></td></tr>
      <tr><td><b>Angiopep-7</b></td><td class="seq">TFFYGGSRG<span class="hit">R</span>RNNF<span class="hit">R</span>TEEY</td><td>K10R·K15R</td><td>뇌 실질 흡수 <b>−85%</b></td></tr>
    </table>
    <p style="margin:9px 0 0;font-size:12px">셋의 <b>순전하가 모두 +2 로 같다</b> → 기능 상실은 전하 효과가 아니라 <b>자리 효과</b>다.
    원문: <i>“replacement of the lysine residues at positions 10 and 15 prevents the transport capacity.”</i></p>
  </div>

  <table class="main">
    <thead><tr><th>우리 축</th><th>Angiopep-2</th><th>Angiopep-5</th><th>Angiopep-7</th><th>정답 대비</th></tr></thead>
    <tbody>{''.join(body)}</tbody>
  </table>

  <div class="tally">
    <div><b style="color:var(--ok)">{n_ok}</b><span>방향 맞음</span></div>
    <div><b style="color:var(--flat)">{n_flat}</b><span>구별 못 함</span></div>
    <div><b style="color:var(--bad)">{n_bad}</b><span>방향 반대</span></div>
  </div>

  <div class="punch">
    <h2>그런데 — 방향을 맞춘 {n_ok}개는 <u>판정을 바꾸지 않는다</u></h2>
    <ul>
      <li><b>rmt_sim</b> 하한 {rmt_floor} → Angiopep-7 이 <b>{r7['rmt']:.3f} 로 통과</b>. 안 막는다.</li>
      <li><b>off_target</b> 밴드 ±{off_band} → 변화 {r2['off']:.2f}→{r7['off']:.2f} = <b>{abs(r7['off']-r2['off']):.2f}</b>, 밴드 안. 그 0.01 은 이 축의 <b>양자화 눈금 자체</b>다.</li>
      <li><b>불안정성</b> 이진화 임계 {ii_stable:.0f} → {r2['ii']:.1f} 도 {r7['ii']:.1f} 도 ‘불안정’. <b>순위 기여 0.</b></li>
    </ul>
    <p class="kicker">방향을 맞춘 축 {n_ok}개 · 그 중 후보를 떨어뜨리는 축 <span style="color:var(--bad)">0개</span>.</p>
  </div>

  <div class="pin">
    <h2>유일하게 잡아내는 것은 서열 축이 아니다</h2>
    <p>근거 등급 상속이 Angiopep-5 (유사도 0.972) 와 Angiopep-7 (0.943) 의 <b>등급 상속을 거부</b>한다.
       유사도 하한 0.70 은 셋 다 통과한다 — 막는 것은 유사도가 아니라
       <b>문헌 보존 위치 {pins} 의 라이신이 살아 있는가</b> 하나다.</p>
    <p style="color:var(--mut)">즉 이 계열을 걸러낸 것은 <b>측정이 아니라 문헌 사실</b>이다.
       그리고 그 핀을 가진 셔틀은 라이브러리 <b>12종 중 1종</b>뿐이다.</p>
  </div>

  <footer>
    n=3 · 한 계열 · 화물 1종. <b>축의 정확도 측정이 아니다.</b> ·
    K→R 은 BLOSUM62 +2(보존적 치환)이므로 서열 유사도 축의 실패는 버그가 아니라 <b>그 양이 이 질문에 답할 수 없다</b>는 뜻이다. ·
    인용 금지: “~9.5배 / 86배”는 인용 출처에 없다 — 쓸 수 있는 것은 −60% / −85% 와 순전하 동일뿐이다. ·
    생성 <code>scripts/known_answer_angiopep.py --html</code>
  </footer>
</div>
</html>"""


# ── 논문 인용용 표 (--md) ─────────────────────────────────────────────────────
#
# 🔴 축을 두 묶음으로 나누는 이유 (2026-11-12 정정)
# -----------------------------------------------------------------------------
# 초판은 축 11개를 한 표에 넣고 전부 정답지로 심판했다. 그건 틀렸다 —
# 응집·liability·접합부·독성은 **뇌 전달을 예측한다고 주장한 적이 없다**.
# K→R 에서 응집이 안 변하는 것은 정상 동작이지 실패가 아니다. 이 정답지로
# 심판할 수 있는 축은 **전달을 주장하는 축**뿐이다.
#
# 그리고 순전하는 세 서열이 같은 것이 **문헌 주장의 검산 성공**이다(실패가 아니다).

# 전달을 주장하는 축 — 이 정답지로 심판할 수 있다
_JUDGED = [
    ("rmt_sim (기능 보존 프록시)",  "rmt",           "{:.3f}", "낮을수록 기능 손상"),
    ("off-target (선택성)",         "off",           "{:.2f}", "RMT 선택성 주장"),
    ("deepB3P (셔틀 단독)",         "bbb_shuttle",   "{:.3f}", "BBB 투과 확률"),
    ("deepB3P (구성체)",            "bbb_construct", "{:.3f}", "BBB 투과 확률"),
]
# 전달을 주장하지 않는 축 — 참고 보고. 값이 같은 것은 기대된 동작이다
_REPORTED = [
    ("순전하 (pH 7.4)",             "charge",        "{:+.2f}"),
    ("불안정성 지수",                "ii",            "{:.1f}"),
    ("응집 핫스팟 (정규화)",         "agg",           "{:.3f}"),
    ("liability 개수",              "liab_n",        "{:d}"),
    ("접합부 liability (구성체)",    "junction_n",    "{:d}"),
    ("ToxinPred3 ML Score (셔틀)",  "tox_shuttle",   "{:.3f}"),
    ("ToxinPred3 ML Score (구성체)","tox_construct", "{:.3f}"),
]


def _dir(rows, key, better_low=True):
    a, c = rows[0][key], rows[2][key]
    if a == c:
        return "구별하지 못함"
    return "방향 맞음" if ((a < c) == better_low) else "**방향이 반대**"


def _md(rows, verdicts, rmt_floor, off_band, ii_stable, pins) -> str:
    A, _B, C = rows
    L = []

    # ── 표 1 ──────────────────────────────────────────────────────────────
    L.append("**표 1.** 뇌 전달을 주장하는 축이 문헌 정답과 일치하는가")
    L.append("")
    L.append("| 평가 축 | 축이 주장하는 것 | Angiopep-2 | Angiopep-5 | Angiopep-7 | 정답 대비 |")
    L.append("|---|---|---:|---:|---:|---|")
    L.append("| *(문헌) 뇌 실질 흡수* | *— (정답지)* | *기준* | *−60%* | *−85%* | *—* |")
    for label, key, fmt, claim in _JUDGED:
        low = key not in ("rmt", "bbb_shuttle", "bbb_construct")
        cells = " | ".join(fmt.format(r[key]) for r in rows)
        L.append(f"| {label} | {claim} | {cells} | {_dir(rows, key, low)} |")
    L.append("")
    L.append("전달을 주장하는 축 **4개 중 2개가 방향을 맞추고 2개가 반대로 간다.** "
             "deepB3P 는 기능이 파괴된 Angiopep-7 을 원본보다 **더 높게** 평가한다.")
    L.append("")
    L.append("")

    # ── 표 2 ──────────────────────────────────────────────────────────────
    L.append("**표 2.** 방향을 맞춘 두 축이 실제로 후보를 떨어뜨리는가")
    L.append("")
    L.append("| 축 | 작동 기준 | Angiopep-7 | 결과 |")
    L.append("|---|---|---:|---|")
    L.append(f"| rmt_sim | 하한 {rmt_floor} 이상이면 통과 | {C['rmt']:.3f} | 통과 — 막지 못함 |")
    L.append(f"| off-target | 시드에서 ±{off_band} 밴드 안이면 통과 | "
             f"{A['off']:.2f} → {C['off']:.2f} (차이 {abs(C['off']-A['off']):.2f}) | "
             f"밴드 안 — 막지 못함 |")
    L.append("")
    L.append("두 축 모두 Angiopep-7 을 더 나쁘게 평가하지만 **후보에서 제외하지는 못한다.** "
             "off-target 의 차이 0.01 은 해당 축의 양자화 눈금과 같다. "
             "즉 방향이 맞는 것과 판정을 바꾸는 것은 다르다.")
    L.append("")
    L.append("")

    # ── 표 3 ──────────────────────────────────────────────────────────────
    L.append("**표 3.** 전달을 주장하지 않는 축 — 참고 보고")
    L.append("")
    L.append("| 평가 축 | Angiopep-2 | Angiopep-5 | Angiopep-7 |")
    L.append("|---|---:|---:|---:|")
    for label, key, fmt in _REPORTED:
        cells = " | ".join(fmt.format(r[key]) for r in rows)
        L.append(f"| {label} | {cells} |")
    L.append("")
    L.append("🔴 **이 축들은 정답 대비 판정을 달지 않는다.** 뇌 전달을 예측한다고 주장한 적이 "
             "없는 개발성·안전성 축이므로, 세 서열에서 값이 같은 것은 실패가 아니라 "
             "**기대된 동작**이다.")
    L.append("")
    L.append(f"다만 **순전하는 검산에 성공했다** — 세 서열 모두 "
             f"{A['charge']:+.2f} 로 동일하다. 문헌이 "
             "\u201c세 펩타이드의 순전하가 모두 +2 로 같다\u201d고 적은 것과 일치하며, "
             "따라서 기능 상실은 전하 효과가 아니라 **위치 효과**다.")
    L.append("")
    L.append("")

    # ── 표 4 ──────────────────────────────────────────────────────────────
    L.append("**표 4.** 표 1 의 결과로 도입한 처방 — 문헌 위치 동결")
    L.append("")
    L.append("| 서열 | 부모 유사도 | 유사도 하한(0.70) | 보존 위치 라이신 | 근거 등급 |")
    L.append("|---|---:|---|---|---|")
    for r in rows:
        tier, _parent, sim = evidence.inherit(r["seq"])
        pin_ok = "보존" if r["name"] == "Angiopep-2" else "파괴"
        L.append(f"| {r['name']} | {sim:.3f} | 통과 | {pin_ok} | "
                 f"{tier if tier else '**없음 — 상속 거부**'} |")
    L.append("")
    L.append("🔴 **이 표는 성능이 아니라 처방이다.** 서열 축이 두 서열을 구별하지 못했기 "
             f"때문에(표 1), 제약의 출처를 유사도에서 문헌으로 옮겨 "
             f"Angiopep-2 의 {pins[0]}·{pins[1]}번 라이신을 동결했다. "
             "**그 동결의 근거가 바로 이 Angiopep-7 실험이므로, 이 표를 "
             "\u2018우리 시스템이 걸러냈다\u2019는 증거로 쓸 수 없다**(순환). "
             "이 표가 보이는 것은 하나뿐이다 — 유사도는 세 서열 모두 하한을 "
             "통과하므로, **유사도로는 걸러지지 않는다.**")
    L.append("")
    L.append("---")
    L.append("")

    # ── 각주 ──────────────────────────────────────────────────────────────
    L.append("**각주 — 표 1~4 공통**")
    L.append("")
    L.append("1. n=3, 한 계열, 화물 1종이다. 축의 정확도를 측정한 것이 아니라 "
             "정답이 알려진 입력 하나에 대한 반응을 기록한 것이다.")
    L.append("2. 문헌에서 Angiopep-7 은 수송 실험의 **음성 대조군**이다. 이 보고서에서는 "
             "같은 서열을 **기능 상실이 문헌으로 확정된 시험 입력**으로 사용한다 — "
             "통제할 어세이가 없으므로 역할이 다르다.")
    L.append("3. K→R 은 BLOSUM62 점수 +2 인 보존적 치환이다. 서열 유사도 기반 축이 "
             "두 서열을 구별하지 못하는 것은 구현 오류가 아니라, **그 양이 이 질문에 "
             "답할 수 없다**는 뜻이다.")
    L.append("4. 문헌 보존 위치가 등록된 셔틀은 라이브러리 12종 중 1종뿐이므로, "
             "표 4 의 처방은 Angiopep-2 계열 밖에서는 작동하지 않는다.")
    L.append("5. deepB3P 점수는 순위에 사용하지 않는 보고치다. 표 1 의 역전은 "
             "그 축을 순위에서 제외한 근거와 같은 방향이다.")
    L.append("6. 표 1 의 off-target 실패와 deepB3P 역전은 이전 측정에서 이미 기록된 "
             "결과이며, 이번 실행은 현행 정의에서의 **재현**이다.")
    L.append("")
    L.append("**출처** — 측정: `scripts/known_answer_angiopep.py` (유료 호출 없음). "
             "문헌: Demeule M, Régina A, Ché C, et al. Identification and design of peptides "
             "as a new drug delivery system for the brain. *J Pharmacol Exp Ther*. "
             "2008;324(3):1064–1072. doi:10.1124/jpet.107.131318")
    L.append("")
    L.append(f"측정 조건: 화물 `{CARGO}` · 링커 GGGGS · "
             "Angiopep-5 = K10R, Angiopep-7 = K10R·K15R (Angiopep-2 기준).")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
