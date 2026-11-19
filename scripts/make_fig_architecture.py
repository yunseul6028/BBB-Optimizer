#!/usr/bin/env python3
"""[그림] 시스템 아키텍처 — 인쇄용 SVG 생성기.

🔴 화면용(Artifact)과 요구가 다르다:
  · **흑백 인쇄**에서 구분돼야 한다 → 색상(hue)에 의존하지 않고 **명도 + 선 종류**로 가른다.
  · 본문 글꼴과 섞이지 않게 **고정 산세리프**를 쓴다.
  · 벡터여야 축소·확대에 안 깨진다.

출력: docs/figs/arch-pipeline.svg  (+ scripts 가 PNG/PDF 를 따로 만든다)
고치려면 이 스크립트를 고치고 다시 돌린다 — SVG 를 손으로 편집하지 말 것.
"""
from pathlib import Path

W = 1000
FONT = "'Apple SD Gothic Neo','Malgun Gothic','Noto Sans KR',sans-serif"
MONO = "'IBM Plex Mono','Menlo',monospace"

# ── 흑백 안전 팔레트 (명도만으로 구분된다) ──
BAND_FILL, BAND_LINE = "#fafafa", "#b4b4b4"
DET_FILL,  DET_LINE  = "#e9e9e9", "#555555"   # 결정론 = 회색 채움
LLM_FILL,  LLM_LINE  = "#ffffff", "#111111"   # LLM   = 흰 채움 + 굵은 테두리
IO_FILL,   IO_LINE   = "#d5d5d5", "#333333"
REP_LINE             = "#8a8a8a"              # 보고치 = 파선 테두리
INK, SUB = "#111111", "#555555"

out = []
def esc(s): return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def band(x, y, w, h, label):
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" fill="{BAND_FILL}" '
               f'stroke="{BAND_LINE}" stroke-width="1" stroke-dasharray="5 4"/>')
    out.append(f'<text x="{x+14}" y="{y+21}" font-family="{FONT}" font-size="14.5" '
               f'font-weight="700" fill="{INK}">{esc(label)}</text>')

def box(x, y, w, h, lines, kind="det", tag=None):
    fill, line, sw, dash = DET_FILL, DET_LINE, 1.1, ""
    if kind == "llm":  fill, line, sw = LLM_FILL, LLM_LINE, 2.2
    if kind == "io":   fill, line, sw = IO_FILL, IO_LINE, 1.4
    if kind == "rep":  fill, line, dash = "#ffffff", REP_LINE, ' stroke-dasharray="5 3"'
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="3" fill="{fill}" '
               f'stroke="{line}" stroke-width="{sw}"{dash}/>')
    n = len(lines); lh = 17.5
    y0 = y + h/2 - (n-1)*lh/2 + 5
    for i, (t, st) in enumerate(lines):
        fs, fw, col, ff = 14, "400", INK, FONT
        if st == "h":  fs, fw = 14.5, "700"
        if st == "s":  fs, col = 12.5, SUB
        if st == "m":  ff, fs, col = MONO, 12, SUB
        out.append(f'<text x="{x+w/2}" y="{y0+i*lh}" text-anchor="middle" font-family="{ff}" '
                   f'font-size="{fs}" font-weight="{fw}" fill="{col}">{esc(t)}</text>')
    if tag:
        tw = 34
        out.append(f'<rect x="{x+w-tw-8}" y="{y+7}" width="{tw}" height="15" rx="2" '
                   f'fill="#111111"/>')
        out.append(f'<text x="{x+w-tw/2-8}" y="{y+18.5}" text-anchor="middle" '
                   f'font-family="{FONT}" font-size="10" font-weight="700" fill="#ffffff">{tag}</text>')

def arrow(x, y1, y2, label=None, dashed=False):
    d = ' stroke-dasharray="6 4"' if dashed else ""
    out.append(f'<line x1="{x}" y1="{y1}" x2="{x}" y2="{y2-7}" stroke="{INK}" '
               f'stroke-width="1.4"{d} marker-end="url(#a)"/>')
    if label:
        out.append(f'<rect x="{x+7}" y="{(y1+y2)/2-9}" width="{len(label)*7.6+10}" height="17" '
                   f'rx="2" fill="#ffffff" stroke="none"/>')
        out.append(f'<text x="{x+12}" y="{(y1+y2)/2+3.5}" font-family="{FONT}" font-size="12" '
                   f'fill="{SUB}">{esc(label)}</text>')

CX, BX, BW = W/2, 30, W-60
y = 16

# 입력
box(CX-130, y, 260, 40, [("사용자 · 화물 펩타이드 입력", "h")], "io"); y += 40
arrow(CX, y, y+26); y += 26

# ① UI
band(BX, y, BW, 92, "① UI — app.py (Streamlit)")
box(BX+22, y+32, 420, 46, [("입력 = 화물 단일 경로", ""), ("모달리티 분기 없음", "s")])
box(BX+476, y+32, 420, 46, [("실시간 이벤트 스트리밍", ""), ("결과 카드 · 과정 기록", "s")])
ui_bottom = y + 92; y = ui_bottom
arrow(CX, y, y+26); y += 26

# ② 설계
band(BX, y, BW, 100, "② 설계 — 자율 에이전트 · 유일한 LLM 생성 지점")
box(BX+230, y+32, 480, 54,
    [("Designer 루프", "h"), ("PLAN → 도구 호출 → finish · 가설 · 예측 · 편집", "s")], "llm", "LLM")
ag_top, ag_bottom = y+32, y+100; y += 100
arrow(CX, y, y+26, "도구 호출"); y += 26

# ③ 채점 — 엔진 6종을 낱개로 (한 상자에 몰아넣으면 글자가 엉긴다)
band(BX, y, BW, 164, "③ 채점 — 엔진 6종 · 전부 결정론")
COLW, COLG = 288, 16
cols = [BX+22, BX+22+COLW+COLG, BX+22+2*(COLW+COLG)]
eng = [("ToxinPred3", "조성 필터", "det"), ("개발성", "liability · 응집", "det"),
       ("ProtParam", "불안정성", "det"), ("off-target", "전하 · GRAVY", "det"),
       ("BLOSUM62", "부모 유사도", "det"),
       ("deepB3P · BBB 투과", "보고치 — 순위에 들어가지 않는다", "rep")]
for i, (t, sub, kind) in enumerate(eng):
    bx, by = cols[i % 3], y + 32 + (i // 3) * 64
    box(bx, by, COLW, 54, [(t, "h"), (sub, "s")], kind)
y += 164
arrow(CX, y, y+26, "채점 결과"); y += 26

# ④ 심사
band(BX, y, BW, 228, "④ 심사 — LLM 판정을 결정론이 덮는다")
box(BX+230, y+32, 480, 46,
    [("Critic · 별도 컨텍스트", "h"), ("APPROVE / WITH_LIMITS / REVISE", "s")], "llm", "LLM")
arrow(CX, y+78, y+100)
gate_top = y+100
box(BX+120, y+100, 700, 70,
    [("결정론 게이트 7 — 단방향", "h"),
     ("tox · inherit · tier_D · cpp_but_rmt", "m"),
     ("evidence_down · dominated_front · ratchet", "m")])
gate_bottom = y+170
arrow(CX, y+170, y+192, "재설계 기회 소진")
box(BX+230, y+192, 480, 46,
    [("복구 — 제출을 러닝 베스트로 교체", "h"),
     ("CPP · 미등재/등급 D · 지배당함", "s")])
y += 228
arrow(CX, y, y+26); y += 26

# ⑤ 순위
band(BX, y, BW, 140, "⑤ 순위 — 사전식 · 가중합이 아니다")
box(BX+180, y+32, 580, 40, [("1순위 · 문헌 근거 그룹", "h")])
arrow(CX, y+72, y+90)
box(BX+120, y+90, 700, 42,
    [("2순위 · 위험 4축 Pareto 비지배", "h"),
     ("liability · 불안정성(40 이진) · off-target · 응집", "s")])
y += 140
arrow(CX, y, y+26); y += 26
box(CX-150, y, 300, 40, [("최종 후보 + 심사 판정", "h")], "io"); y += 40

# REVISE 되먹임 — 오른쪽 바깥으로 돌아 올라간다
rx = W - 16
out.append(f'<path d="M {BX+820} {gate_top+35} H {rx} V {ag_top+27} H {BX+710+7}" '
           f'fill="none" stroke="{INK}" stroke-width="1.4" stroke-dasharray="6 4" '
           f'marker-end="url(#a)"/>')
# 라벨은 ②와 ③ 띠 **사이의 빈 줄**에 놓는다 — 띠 안에 두면 상자를 덮는다.
lab_y = ag_bottom + 13
out.append(f'<rect x="{rx-128}" y="{lab_y-9}" width="118" height="18" rx="2" fill="#ffffff"/>')
out.append(f'<text x="{rx-14}" y="{lab_y+4}" text-anchor="end" font-family="{FONT}" '
           f'font-size="11.5" font-weight="700" fill="{INK}">REVISE · 재설계 1회</text>')

H = y + 18
svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
       f'viewBox="0 0 {W} {H}">'
       f'<defs><marker id="a" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
       f'markerHeight="6" orient="auto-start-reverse">'
       f'<path d="M 0 0 L 10 5 L 0 10 z" fill="{INK}"/></marker></defs>'
       f'<rect width="{W}" height="{H}" fill="#ffffff"/>' + "".join(out) + "</svg>")
p = Path("docs/figs/arch-pipeline.svg"); p.write_text(svg, encoding="utf-8")
print(f"{p} · {W}×{H}")
