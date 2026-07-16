# -- coding: utf-8 --
"""FBGAN 생성 최적화 루프 (잠재공간 진화 피드백).

사전학습 생성기(G_weights_1000.pth)를 고정 디코더로 쓰고, latent z 를 진화시켜
'cargo + 링커 + 생성셔틀'의 BBB(deepB3P)를 높이고 독성(ToxinPred3)을 낮춘다.

  라운드: z 샘플 → G 디코딩(=셔틀 후보) → 조립 → BBB·독성 채점
          → fitness = 비독성이면 BBB, 독성이면 페널티
          → 상위 elite z 선택 → 주변 변이 + 소량 랜덤으로 다음 세대

deepB3P는 동일 venv에서 in-process 채점(5-fold 1회 로딩), ToxinPred3는 별 venv subprocess.

Usage:
  python _run_fbgan.py <cargo> <linker> <rounds> <pop> <elite> <tox_thr> <toxpy> <toxrepo> <toxscript> <out.json>
"""
import functools
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import torch

torch.load = functools.partial(torch.load, map_location="cpu")

from fbgan.models import Generator
from model.deepb3p import DeepB3P
from utils.amino_acid import from_amino_acid_to_id, from_id_to_amino_acid, remove_pad_from_str_seq
from utils.config_transformer import Config as B3Config
from utils.utils import SeqDataset
from torch.utils.data import DataLoader

SEQ_LEN, N_CHARS, HIDDEN, MODEL_MAX = 50, 21, 128, 50


# ---------------- deepB3P in-process 채점기 (5-fold 1회 로딩) ----------------
class BBBScorer:
    def __init__(self):
        p = B3Config(d_model=512, d_ff=16, d_k=32, n_layers=1, n_heads=2, lr=0.0001, drop=0.1)
        p.make_dir()
        logger = p.set_logging(str(p.model_file / "deepb3p.log"))
        self.device = p.device
        self.models = []
        for i in range(1, p.kFold + 1):
            m = DeepB3P(p, logger)
            m.model.reset_parameters()
            m.load_model(p.model_file / f"deepb3p_{i}.pth")
            m.model.eval()
            self.models.append(m)

    def _encode(self, seqs):
        arr = []
        for s in seqs:
            # 50aa 초과는 연결부위(C말단) 윈도우로 — 앞(화물)이 아닌 뒤(링커+셔틀) 보존
            s = (s if len(s) <= SEQ_LEN else s[-SEQ_LEN:]).ljust(SEQ_LEN, "0")
            arr.append(np.asarray(from_amino_acid_to_id(s)))
        return np.stack(arr).astype(np.int64)

    def score(self, seqs):
        feats = self._encode(seqs)
        loader = DataLoader(SeqDataset(feats, np.ones(len(seqs), dtype=np.int64)), batch_size=len(seqs))
        probs_sum = np.zeros(len(seqs))
        with torch.no_grad():
            for m in self.models:
                for x, _ in loader:
                    out = m.model(x.to(self.device))
                    probs_sum += out[:, 1].cpu().numpy()
        return probs_sum / len(self.models)


# ---------------- ToxinPred3 subprocess 채점기 ----------------
def score_tox(seqs, toxpy, toxrepo, toxscript):
    with tempfile.TemporaryDirectory() as tmp:
        fa, out = Path(tmp) / "in.fa", Path(tmp) / "out.csv"
        with open(fa, "w") as f:
            for i, s in enumerate(seqs):
                f.write(f">{i}\n{s}\n")
        r = subprocess.run([toxpy, toxscript, "-i", str(fa), "-o", str(out), "-m", "1", "-d", "2"],
                           cwd=toxrepo, capture_output=True, text=True, timeout=300)
        if r.returncode != 0 or not out.exists():
            raise RuntimeError(f"ToxinPred3 실패: {r.stderr[-300:]}")
        import csv
        rows = {row["ID"]: float(row["ML Score"]) for row in csv.DictReader(open(out))}
    return [rows[str(i)] for i in range(len(seqs))]


JUNCTION_FLANK = 0   # 0 = 화물 완전 제외(링커+셔틀 = BBB 모듈)


def bbb_region(cargo, linker, shuttle):
    """BBB 계산용: 전체가 SEQ_LEN 초과면 연결부위(링커+셔틀=BBB 모듈)만 남김."""
    full = cargo + linker + shuttle
    if len(full) <= SEQ_LEN:
        return full
    tail = cargo[-JUNCTION_FLANK:] if JUNCTION_FLANK > 0 else ""
    return (tail + linker + shuttle)[-SEQ_LEN:]


def decode(G, z):
    G.bs = z.shape[0]
    G.eval()
    with torch.no_grad():
        hard = np.argmax(G(torch.from_numpy(z).float()).cpu().numpy(), axis=2)
    seqs = []
    for row in hard:
        s = remove_pad_from_str_seq(list("".join(from_id_to_amino_acid(row))), "0")
        seqs.append(s if s else "A")   # 빈 서열 방지
    return seqs


def main(cargo, linker, rounds, pop, elite, tox_thr, toxpy, toxrepo, toxscript, out_path):
    G = Generator(n_chars=N_CHARS, seq_len=SEQ_LEN, bs=pop, hidden=HIDDEN)
    G.load_state_dict(torch.load("fbgan/checkpoint/G_weights_1000.pth"))
    bbb = BBBScorer()

    rng = np.random.RandomState(2022)
    z = rng.randn(pop, 128)
    sigma = 0.6
    history, best = [], {}

    for rd in range(rounds):
        shuttles = decode(G, z)
        constructs = [cargo + linker + s for s in shuttles]
        bbb_scores = bbb.score([bbb_region(cargo, linker, s) for s in shuttles])
        tox_scores = score_tox(constructs, toxpy, toxrepo, toxscript)

        fitness = np.array([
            b if t <= tox_thr else b - 1.0            # 독성이면 페널티
            for b, t in zip(bbb_scores, tox_scores)
        ])
        order = np.argsort(fitness)[::-1]

        # 전역 베스트 갱신 (비독성 중 최고 BBB)
        for i in order:
            if tox_scores[i] <= tox_thr:
                cand = {"shuttle": shuttles[i], "sequence": constructs[i],
                        "bbb": round(float(bbb_scores[i]), 4), "tox": round(float(tox_scores[i]), 4),
                        "len": len(constructs[i])}
                if not best or cand["bbb"] > best["bbb"]:
                    best = cand
                break

        safe = [b for b, t in zip(bbb_scores, tox_scores) if t <= tox_thr]
        history.append({"round": rd + 1,
                        "mean_bbb": round(float(np.mean(bbb_scores)), 4),
                        "best_bbb": round(float(np.max(bbb_scores)), 4),
                        "n_safe": len(safe)})

        # --- 피드백: 상위 elite z 주변 변이 + 소량 랜덤 재샘플 ---
        elites = z[order[:elite]]
        children = []
        while len(children) < pop - elite - 2:
            parent = elites[rng.randint(elite)]
            children.append(parent + sigma * rng.randn(128))
        z = np.vstack([elites, np.array(children), rng.randn(2, 128)])  # elite 보존 + 탐색용 랜덤2
        sigma *= 0.85   # 점진 수렴

    # 마지막 세대 상위 후보들도 수집
    final_seqs = decode(G, z)
    fc = [cargo + linker + s for s in final_seqs]
    fb = bbb.score([bbb_region(cargo, linker, s) for s in final_seqs])
    ft = score_tox(fc, toxpy, toxrepo, toxscript)
    pool = [{"shuttle": s, "sequence": c, "bbb": round(float(b), 4), "tox": round(float(t), 4), "len": len(c)}
            for s, c, b, t in zip(final_seqs, fc, fb, ft) if t <= tox_thr]
    if best:
        pool.append(best)
    # 중복 제거 + BBB 내림차순 top
    seen, top = set(), []
    for item in sorted(pool, key=lambda x: x["bbb"], reverse=True):
        if item["sequence"] in seen:
            continue
        seen.add(item["sequence"]); top.append(item)

    json.dump({"history": history, "best": top[:5], "cargo": cargo, "linker": linker},
              open(out_path, "w"), ensure_ascii=False)
    print("OK")


if __name__ == "__main__":
    a = sys.argv
    main(a[1], a[2], int(a[3]), int(a[4]), int(a[5]), float(a[6]), a[7], a[8], a[9], a[10])
