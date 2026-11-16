"""§64-4 분산 가설을 **새 라이브 없이** 셀 수 있는 데까지 센다 (§68-8).

    PYTHONPATH=. /usr/bin/python3 scripts/var_batch_dispersion.py   # 유료 0회 · ~6분

가설: 프롬프트 정정 후 **묶음 간 분산**이 커졌다(0.027 → 0.132).
🔴 교란 후보: 묶음당 세션이 **3개뿐**이라 **이항 잡음만으로도** 묶음 간 분산이 생기고,
   그 크기는 평균이 0.5 에 가까울수록 **커진다**. 평균이 0.867 → 0.417 로 움직였으므로
   분산 증가의 일부(혹은 전부)는 **평균 이동의 산술적 부산물**일 수 있다.
"""
import itertools, random, statistics as st

PRE  = [3, 3, 2, 2, 3]      # TFEY 수 / 묶음(각 3세션)
POST = [0, 1, 3, 1]
N = 3                        # 묶음당 세션 수

def pvar(counts):
    xs = [c / N for c in counts]
    m = sum(xs) / len(xs)
    return sum((x - m) ** 2 for x in xs) / len(xs)   # 모분산(문서와 같은 정의)

vp, vq = pvar(PRE), pvar(POST)
mp = sum(PRE) / (N * len(PRE)); mq = sum(POST) / (N * len(POST))
print(f"관측  변경전 평균 {mp:.4f} · 분산 {vp:.4f}   |   변경후 평균 {mq:.4f} · 분산 {vq:.4f}")
print(f"      분산 차 {vq - vp:+.4f} · 분산 비 {vq / vp:.2f}배\n")

# ── ① 묶음 단위 정확 순열검정 (분산 차) — 라벨만 섞는다
pool = PRE + POST
obs = vq - vp
cnt_one = cnt_two = 0; tot = 0
for idx in itertools.combinations(range(len(pool)), len(PRE)):
    a = [pool[i] for i in idx]
    b = [pool[i] for i in range(len(pool)) if i not in idx]
    d = pvar(b) - pvar(a)
    tot += 1
    if d >= obs - 1e-12: cnt_one += 1
    if abs(d) >= abs(obs) - 1e-12: cnt_two += 1
print(f"① 정확 순열검정(묶음 라벨 교환, 전체 {tot}가지)")
print(f"   단측 p = {cnt_one}/{tot} = {cnt_one/tot:.3f}  ·  양측 p = {cnt_two}/{tot} = {cnt_two/tot:.3f}")
print(f"   🔴 달성 가능한 **최소** p = 1/{tot} = {1/tot:.3f} — 분모가 이미 한계다\n")

# ── ② 이항 잡음만으로 기대되는 묶음 간 분산 (교란 확인)
print("② 교란 검사 — 묶음당 3세션이면 **이항 잡음만으로도** 분산이 생긴다")
for name, m, obsv, k in (("변경 전", mp, vp, len(PRE)), ("변경 후", mq, vq, len(POST))):
    exp_bin = m * (1 - m) / N
    print(f"   {name}: 평균 {m:.3f} → 이항 기대 분산 **{exp_bin:.4f}** vs 관측 {obsv:.4f} "
          f"({'과소' if obsv < exp_bin else '과대'}산포 {obsv/exp_bin:.2f}배)")
print("   ⇒ 평균이 0.5 에 가까워지면 이항 기대 분산 자체가 커진다 — "
      "분산 증가의 일부는 **평균 이동의 산술적 부산물**이다.\n")

# ── ③ 과산포 자체를 검정한다 (각 팔의 평균을 고정한 모수 부트스트랩)
#    귀무: 각 묶음이 같은 p(팔 평균)에서 독립 이항 → 추가 산포 없음
rng = random.Random(0)
B = 200000
print(f"③ 과산포 검정 — 귀무 '묶음 간 추가 변동 없음(순수 이항)', 모수 부트스트랩 {B:,}회")
for name, m, obsv, k in (("변경 전", mp, vp, len(PRE)), ("변경 후", mq, vq, len(POST))):
    ge = 0
    for _ in range(B):
        sim = [sum(1 for _ in range(N) if rng.random() < m) for _ in range(k)]
        if pvar(sim) >= obsv - 1e-12: ge += 1
    print(f"   {name}: p = {ge/B:.3f}  (관측 분산 {obsv:.4f} 이상이 나올 확률)")
print()

# ── ④ 검정력 — 묶음을 몇 개 더 모아야 하나 (관측값을 참값으로 가정한 낙관적 추정)
print("④ 검정력 추정 — 관측값(과산포 비)을 **참값으로 가정한 낙관적** 수치")
def sim_batch(m, extra, rng):
    """베타-이항: 평균 m, 과산포를 extra 배로."""
    if extra <= 1.0:
        p = m
    else:
        # var(p) = m(1-m)*rho 가 되도록 베타 모수를 잡는다
        rho = min(0.99, (extra - 1) / (N - 1))
        a = m * (1 / rho - 1); b = (1 - m) * (1 / rho - 1)
        p = rng.betavariate(max(a, 1e-6), max(b, 1e-6))
    return sum(1 for _ in range(N) if rng.random() < p)

ex_pre = vp / (mp * (1 - mp) / N); ex_post = vq / (mq * (1 - mq) / N)
print(f"   가정: 변경 전 과산포 {ex_pre:.2f}배 · 변경 후 {ex_post:.2f}배 (관측 그대로)")
for k in (5, 8, 12, 20, 30):
    hit = 0; T = 2000
    for _ in range(T):
        A = [sim_batch(mp, ex_pre, rng) for _ in range(k)]
        Bv = [sim_batch(mq, ex_post, rng) for _ in range(k)]
        pool2 = A + Bv; o = pvar(Bv) - pvar(A)
        # 순열 근사(1000회)
        ge = 0; R = 1000
        for _ in range(R):
            rng.shuffle(pool2)
            if pvar(pool2[k:]) - pvar(pool2[:k]) >= o - 1e-12: ge += 1
        if (ge + 1) / (R + 1) <= 0.05: hit += 1
    print(f"   팔당 묶음 {k:>2}개 → 검정력 ≈ **{hit/T:.0%}**  (= 라이브 세션 {2*k*N}회)")
