"""구조 기반 접합부 분석 (ESMFold + 셔틀 용매노출).

서열 모델(deepB3P)이 못 보는 것을 본다: **융합체를 실제로 접었을 때 셔틀이 표면에
노출됐는가(수용체 결합 가능) vs 화물에 가려졌는가(occluded).**

  1) ESMFold 공개 API로 융합체(cargo+linker+shuttle) 구조 예측 → PDB + pLDDT
  2) Biopython Shrake-Rupley SASA로 셔틀 영역의 상대 용매노출(RSA) 계산
  3) 셔틀 노출도 + 구조 신뢰도(pLDDT)로 판정

⚠️ 짧은 펩타이드는 무질서 경향(pLDDT 낮음) → 정적 구조는 참고용. 검증된 구조→BBB
   모델은 없으며 셔틀 노출도는 **기능 프록시(휴리스틱)**. 실험 검증 필요.
"""

from __future__ import annotations

import hashlib
import io
from dataclasses import dataclass, field

from .config import BASE_DIR

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
CACHE_DIR = BASE_DIR / ".cache" / "esmfold"

# 잔기별 최대 용매접근면적 (Tien et al. 2013, theoretical, Å²) — RSA 정규화용
MAX_ASA = {
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLU": 223,
    "GLN": 225, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
    "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285,
    "TYR": 263, "VAL": 174,
}


@dataclass
class StructureResult:
    cargo: str = ""
    linker: str = ""
    shuttle: str = ""
    sequence: str = ""
    pdb: str = ""
    n_residues: int = 0
    mean_plddt: float = 0.0        # 0~100
    shuttle_plddt: float = 0.0
    shuttle_exposure: float = 0.0  # 셔틀 영역 평균 RSA (0~1)
    per_residue: list = field(default_factory=list)  # [{i,resn,region,rsa,plddt}]
    verdict: str = ""
    exposed: bool = False
    low_confidence: bool = False
    error: str = ""


def fold_esmfold(sequence: str, timeout: float = 25.0, tries: int = 3,
                 deadline: float = 60.0) -> tuple[str, str]:
    """서열 → ESMFold PDB. (pdb, error) 반환. sha256 캐싱 + 일시적 서버오류 백오프 재시도.

    ESMFold 공개 API(api.esmatlas.com)는 종종 502/503/**504(게이트웨이 타임아웃)** 나 응답
    지연을 내는데, 대부분 일시적이라 짧게 기다렸다 재시도하면 성공한다. 4xx(클라이언트 오류)는
    재시도 안 함. 다만 서버가 통째로 내려간 경우 재시도는 소용없으므로 **빨리 실패**가 중요하다:
      · per-try timeout을 짧게(기본 25초) 잡고,
      · `deadline`(기본 60초) 총 대기 상한을 두어 서버가 hang해도 UI가 오래 멎지 않게 한다.
    구조 노출도는 **보조 지표**라 실패해도 나머지 결과(BBB·독성·전달 분해)엔 영향이 없다.
    """
    import time
    import requests
    seq = "".join(ch for ch in sequence.upper() if ch.isalpha())
    cache = CACHE_DIR / (hashlib.sha256(seq.encode()).hexdigest() + ".pdb")
    if cache.exists():
        return cache.read_text(), ""

    last = ""
    start = time.monotonic()
    for attempt in range(tries):
        # 남은 예산이 per-try 타임아웃보다 작으면 그만큼만 기다린다(총 deadline 준수).
        remaining = deadline - (time.monotonic() - start)
        if remaining <= 0:
            break
        try:
            r = requests.post(ESMFOLD_URL, data=seq, timeout=min(timeout, remaining))
            if r.status_code in (429, 500, 502, 503, 504):
                last = f"서버 {r.status_code}(일시적)"          # transient → 재시도
            elif not r.ok:
                return "", f"ESMFold API 오류: HTTP {r.status_code}"  # 4xx 등 → 즉시 실패
            else:
                pdb = r.text
                if "ATOM" not in pdb:
                    last = f"응답 이상: {pdb[:120]}"
                else:
                    CACHE_DIR.mkdir(parents=True, exist_ok=True)
                    cache.write_text(pdb)
                    return pdb, ""
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as exc:
            last = type(exc).__name__                          # 네트워크 지연 → 재시도
        except Exception as exc:  # noqa: BLE001
            return "", f"ESMFold API 오류: {type(exc).__name__}: {exc}"
        # 다음 시도 전 짧은 백오프(2초, 4초) — 단, 총 deadline을 넘기지 않을 때만.
        if attempt < tries - 1 and (deadline - (time.monotonic() - start)) > 2:
            time.sleep(2 * (attempt + 1))

    return "", (f"ESMFold 서버 일시 오류({last or '응답 없음'}) — 재시도 후에도 실패했습니다. "
                "잠시 후 다시 시도하세요. (구조 노출도는 보조 지표라 나머지 결과엔 영향 없습니다.)")


def analyze_construct(cargo: str, linker: str, shuttle: str,
                      timeout: float = 25.0) -> StructureResult:
    """융합체를 접고 셔틀 노출도를 분석한다."""
    cargo, linker, shuttle = cargo.upper(), linker.upper(), shuttle.upper()
    seq = cargo + linker + shuttle
    res = StructureResult(cargo=cargo, linker=linker, shuttle=shuttle, sequence=seq)

    pdb, err = fold_esmfold(seq, timeout=timeout)
    if err:
        res.error = err
        return res
    res.pdb = pdb

    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    structure = PDBParser(QUIET=True).get_structure("f", io.StringIO(pdb))
    ShrakeRupley().compute(structure, level="R")
    model = next(structure.get_models())
    residues = [r for r in model.get_residues() if r.id[0] == " "]
    n = len(residues)
    res.n_residues = n

    # 셔틀 = C말단 shuttle_len 잔기 (construct = cargo+linker+shuttle)
    sh_len = min(len(shuttle), n)
    sh_start = n - sh_len

    raw = []
    for r in residues:
        bfacs = [a.get_bfactor() for a in r.get_atoms()]
        raw.append(sum(bfacs) / len(bfacs) if bfacs else 0.0)
    # ESMFold pLDDT가 0~1 스케일이면 0~100으로 정규화
    scale = 100.0 if (raw and max(raw) <= 1.5) else 1.0

    per, plddts = [], []
    for i, r in enumerate(residues):
        maxa = MAX_ASA.get(r.resname)
        rsa = min(1.0, getattr(r, "sasa", 0.0) / maxa) if maxa else 0.0
        plddt = raw[i] * scale
        plddts.append(plddt)
        region = "shuttle" if i >= sh_start else ("linker" if i >= len(cargo) else "cargo")
        per.append({"i": i + 1, "resn": r.resname, "region": region,
                    "rsa": round(rsa, 3), "plddt": round(plddt, 1)})
    res.per_residue = per

    sh = [p for p in per if p["region"] == "shuttle"]
    res.mean_plddt = round(sum(plddts) / n, 1) if n else 0.0
    res.shuttle_exposure = round(sum(p["rsa"] for p in sh) / len(sh), 3) if sh else 0.0
    res.shuttle_plddt = round(sum(p["plddt"] for p in sh) / len(sh), 1) if sh else 0.0

    # 판정 (pLDDT<70 = 저신뢰; AlphaFold/ESMFold 관례)
    res.low_confidence = res.shuttle_plddt < 70.0
    e = res.shuttle_exposure
    res.exposed = e >= 0.35
    if e >= 0.45:
        res.verdict = "셔틀이 잘 노출됨 — 수용체 결합에 유리"
    elif e >= 0.30:
        res.verdict = "셔틀 부분 노출 — 접근성 보통"
    else:
        res.verdict = "셔틀이 파묻힘(occluded) 우려 — 화물이 가릴 수 있음"
    if res.low_confidence:
        res.verdict += " · ⚠️구조 신뢰도 낮음(유연/무질서 가능)"
    return res
