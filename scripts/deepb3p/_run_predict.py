# -- coding: utf-8 --
"""CPU 추론 런너 (BBB-Optimizer에서 subprocess로 호출).

deepB3P의 predict_user.predict() 로직을 그대로 따르되:
  - torch.load 를 map_location='cpu' 로 강제 (GPU 저장 체크포인트를 CPU에서 로드)
  - 출력 경로를 인자로 받음 (prob.txt 하드코딩 회피, 동시 호출 충돌 방지)

Usage: python _run_predict.py <input.fasta> <output.csv>
"""
import functools
import sys

import torch

# GPU에서 저장된 5-fold 체크포인트를 CPU-only 환경에서 로드하기 위한 패치.
torch.load = functools.partial(torch.load, map_location="cpu")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from torch.utils.data import DataLoader  # noqa: E402

from model.deepb3p import DeepB3P  # noqa: E402
from utils.amino_acid import fasta_to_numpy  # noqa: E402
from utils.config_transformer import Config  # noqa: E402
from utils.utils import SeqDataset, transfer  # noqa: E402


def get_seqs(path):
    from Bio.SeqIO import parse
    return [str(r.seq) for r in parse(path, "fasta")]


def main(fasta_path: str, out_path: str):
    # deepb3p_1~5.pth 를 학습한 하이퍼파라미터 (체크포인트 폴더명과 일치해야 함)
    params = Config(d_model=512, d_ff=16, d_k=32, n_layers=1, n_heads=2, lr=0.0001, drop=0.1)
    params.make_dir()  # reload=True 라 체크포인트 삭제 분기는 타지 않음
    logger = params.set_logging(str(params.model_file / "deepb3p.log"))

    feats, label = fasta_to_numpy(fasta_path, label=1)
    if feats is None:
        raise SystemExit("no sequences parsed from fasta")

    loader = DataLoader(SeqDataset(feats, label), batch_size=len(feats))

    # 5-fold 앙상블 확률 평균
    prob_df = pd.DataFrame()
    for i in range(1, params.kFold + 1):
        model = DeepB3P(params, logger)
        model.model.reset_parameters()
        model.load_model(params.model_file / f"deepb3p_{i}.pth")
        with torch.no_grad():
            model.model.eval()
            probs = []
            for x, y in loader:
                x = x.to(params.device)
                out = model.model(x)
                probs.extend(out[:, 1].cpu().numpy())
        prob_df[i] = probs

    avg = np.round(prob_df.mean(axis=1), 4)
    res = pd.DataFrame({"peptide": get_seqs(fasta_path), "prob": avg})
    res.to_csv(out_path, index=False)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python _run_predict.py <input.fasta> <output.csv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
    print("OK")
