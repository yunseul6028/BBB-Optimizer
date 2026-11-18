#!/usr/bin/env python3
"""§81 — 심사 게이트 발화 기록을 읽어 "재설계 기회를 늘리면 되는 문제인가"에 답한다.

배경: 게이트는 위반을 **검출**해 REVISE 를 낼 뿐이고, 재설계 기회(MAX_CRITIQUES=1)를
소진하면 그 후보가 그대로 최종이 된다. 그래서 **검출 + 복구**를 한 쌍으로 두었다.
"그럼 재설계 기회를 더 주면 되는 것 아닌가?" 를 가르는 유일한 관측이 이것이다:

    재시도가 **같은 게이트**에 다시 걸렸나, **다른 게이트**였나?

  · 같은 게이트 반복  → 예산 문제가 **아니다**(못 고치거나 안 고치는 조건)
  · 매번 다른 게이트  → 수렴 중일 수 있다. 예산이 실제로 도울 수 있다

🔴 구조적으로는 예산을 올려도 구멍이 **옮겨질 뿐 사라지지 않는다** — 어떤 유한 k 에도
`k+1` 번째 실패가 남는다. 이 보고서는 "복구를 뗄 수 있나"가 아니라 **"예산을 올리면
복구 발동 빈도가 줄어드나"** 에만 답한다.

사용: /usr/bin/python3 scripts/critic_gate_report.py
"""
import json, sys, collections
from pathlib import Path

LOG = Path(__file__).resolve().parent.parent / ".cache" / "critic_gates.jsonl"
GATE_KO = {"tox": "① 독성(조성 필터)", "inherit": "② 근거 등급 보유(상속)",
           "tier_D": "③ 등급 D(BBB 통과 단독 근거 없음)", "cpp_but_rmt": "④ CPP인데 RMT 실재",
           "evidence_down": "⑤ 근거 하향", "dominated_front": "⑥ 전선에 지배당함",
           "ratchet": "⑦ 래칫(이미 본 최고에 지배당함)", "APPROVE": "— (게이트 미발화)"}


def main() -> int:
    if not LOG.exists():
        print(f"기록이 아직 없다: {LOG}")
        print("라이브 세션을 한 번도 돌리지 않았다는 뜻이다. §81 이후 실행부터 쌓인다.")
        return 0
    rows = [json.loads(l) for l in LOG.read_text(encoding="utf-8").splitlines() if l.strip()]
    if not rows:
        print("기록 파일이 비어 있다.")
        return 0

    n = len(rows)
    print(f"\n세션 {n}건  ({LOG})\n" + "=" * 74)

    # ── ⒜ 핵심 질문 — 재시도가 같은 칸에 다시 걸렸나 ──
    # 🔴 `gate_log` 길이만으로는 "재설계를 했는가"가 갈리지 않는다. `swapped_out` 이면
    #    복구 교체 **뒤** 재심사가 1건 더 붙기 때문이다. 루프 안 심사 횟수 = 길이 − (교체면 1).
    #    (옛 필터는 '게이트 발화 2회 이상'이라, **한 번 걸리고 재설계 후 통과**한 세션 —
    #     가장 알고 싶은 사례 — 을 통째로 놓쳤다. 2026-11-18 수정.)
    def loop_calls(r):
        return len(r["gate_log"]) - (1 if r.get("swapped_out") else 0)

    fired = [r for r in rows if any(g != "APPROVE" for g in r["gate_log"])]
    retried = [r for r in rows if (r.get("critiques_used") or 0) >= 1 or loop_calls(r) >= 2]
    print(f"\n⒜ 🔴 핵심 — 게이트가 쏜 세션: {len(fired)}/{n} · 그중 **재제출까지 간** 세션: "
          f"{sum(1 for r in retried if loop_calls(r) >= 2)}/{len(fired) or 1}")
    unused = [r for r in fired if loop_calls(r) < 2]
    if unused:
        print(f"   🔴 **재설계 기회를 주고도 다시 제출하지 않은 세션: {len(unused)}**")
        print("      (루프 안 심사가 1회뿐 = REVISE 뒤 `finish` 를 다시 부르지 않고 예산을 태웠다)")
        print("      → 이 경우 **예산을 늘려도 소용없다**. 쓰지 않는 기회를 더 주는 셈이다.")
        for r in unused:
            print(f"        · {r['ts']}  {r['gate_log']}  swap={r.get('swapped_out')}  {r['pick']}")
    resub = [r for r in retried if loop_calls(r) >= 2]
    if not resub:
        print("\n   재제출이 한 번도 없어 **'같은 게이트냐 다른 게이트냐'는 아직 답할 수 없다**.")
        print("   (분모 0 — 어떤 결론도 내지 말 것)")
    else:
        retried = resub
        same = [r for r in retried if len(set(g for g in r["gate_log"] if g != "APPROVE")) == 1]
        print(f"   같은 게이트 반복 : {len(same)}/{len(retried)}  → 예산 문제가 **아니라는** 쪽")
        print(f"   다른 게이트로 이동: {len(retried)-len(same)}/{len(retried)}  → 예산이 도울 수 있다는 쪽")
        for r in retried:
            arrow = " → ".join(GATE_KO.get(g, g).split()[0] for g in r["gate_log"])
            print(f"     · {r['ts']}  {r['pick']:<28} {arrow}")
        if len(retried) < 5:
            print(f"   ⚠️ 분모가 {len(retried)} 이다. 방향만 보고 **비율로 쓰지 말 것**.")

    # ── ⒝ 게이트별 발화 ──
    c = collections.Counter(g for r in rows for g in r["gate_log"] if g != "APPROVE")
    print(f"\n⒝ 게이트별 발화 (총 {sum(c.values())}회)")
    for g, v in c.most_common():
        print(f"   {v:>3}  {GATE_KO.get(g, g)}")
    silent = [g for g in GATE_KO if g != "APPROVE" and g not in c]
    if silent:
        print(f"   🔴 한 번도 안 쏜 게이트 {len(silent)}개: "
              + ", ".join(GATE_KO[g] for g in silent))
        print("      → M5 의 '7/7' 을 **이 기록으로 재주장하지 말 것**(관측되지 않은 칸이 있다).")

    # ── ⒞ 복구 발동 · LLM vs 게이트 불일치(§77) ──
    sw = sum(1 for r in rows if r.get("swapped_out"))
    ap = sum(1 for r in rows if r.get("approved"))
    print(f"\n⒞ 복구(러닝 베스트 교체) 발동: {sw}/{n} · 최종 승인: {ap}/{n}")
    dis = [r for r in rows if r.get("llm_tier") and r["llm_tier"] != r.get("tier")]
    print(f"\n⒟ §77 — LLM 자신의 판정 ≠ 게이트 판정: {len(dis)}/{n}")
    if dis:
        for r in dis[:6]:
            print(f"     · {r['ts']}  LLM={r['llm_tier']:<20} 게이트={r['tier']}  ({r.get('gate_final')})")
        print("   → LLM 심사자가 게이트와 **다른 판단**을 한 사례다. 기여를 말할 수 있는 첫 재료.")
    else:
        print("   → 불일치 0건. 🔴 '게이트가 LLM 과 같은 판단을 한다'가 아니라 "
              "**아직 갈린 사례가 없다**로 읽을 것.")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
