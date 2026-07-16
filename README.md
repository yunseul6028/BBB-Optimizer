# 🧬 BBB-Optimize AI Agent

뇌혈관장벽(BBB)을 통과하는 알츠하이머 치료용 **항체-셔틀 융합 단백질**을 설계하는
**자율 설계 에이전트**. 화물(cargo) 펩타이드를 넣으면, **Claude가 여러 예측 도구를 스스로
골라 써가며** BBB 투과(효능)·독성(안전성)·입체구조·수용체 결합을 모두 만족하는 최종 융합체
(`화물 + 링커 + 셔틀`)로 수렴한다.

## 🤖 메인: 자율 설계 에이전트

Claude(claude-opus-4-8)가 아래 도구들을 **스스로 오케스트레이션**한다 — 스크리닝 → 개선 →
구조 검증 → (필요시) 신규 셔틀 생성 → 최종 융합체 도출. (Anthropic API 키 필요)

| 에이전트 도구 | 하는 일 | 백엔드 엔진 |
|---|---|---|
| `evaluate_candidates` | 후보 융합체의 BBB 투과 점수·독성·수용체 유사도 배치 평가 | deepB3P + ToxinPred3 + BLOSUM62 |
| `analyze_structure` | 최종 후보를 접어 셔틀이 표면에 노출됐는지 검증 (3D) | ESMFold + Biopython |
| `generate_novel_shuttles` | 라이브러리 밖 신규 셔틀을 생성 | FBGAN |

> 각 도구는 앱의 **"🔧 개별 도구 직접 실행"** 패널에서 키 없이 무료로 수동 실행도 가능하다
> (라이브러리 전수 스윕 · 이중 트랙 정밀 분석 · FBGAN 생성).

## 예측 엔진

| 축 | 엔진 | 방식 |
|---|---|---|
| BBB 투과 | [deepB3P](https://github.com/GreatChenLab/deepB3P) | 펩타이드 서열 딥러닝 (로컬) |
| 독성 | [ToxinPred3](https://github.com/raghavagps/toxinpred3) | AAC+DPC 조성 (로컬, GPL-3.0) |
| 생성 | FBGAN (deepB3P 동봉) | 사전학습 생성기 + 잠재공간 진화 |
| 구조 | [ESMFold](https://esmatlas.com/) | 공개 API 폴딩 |
| 추론 | Claude (Anthropic API) | 자율 최적화 에이전트 (선택) |

> ⚠️ deepB3P는 짧은 펩타이드 학습 모델 — 긴 융합체는 조합 간 **상대 비교**로 해석.
> 구조 노출도는 검증된 구조→BBB 모델이 아닌 **휴리스틱 프록시**. 실제 적용 전 In Vitro/In Vivo 검증 필요.

---

## 설치

### 1. 앱 본체
```bash
pip install -r requirements.txt        # streamlit
pip install biopython requests py3Dmol # 구조 분석(Track 1)용
```

### 2. 로컬 예측 엔진 (별도 클론 + 모델 + venv)
제3자 코드·모델 가중치는 이 repo에 포함하지 않는다(라이선스/용량). 아래로 재현:

```bash
mkdir -p vendor && cd vendor

# --- deepB3P (BBB + FBGAN) ---
git clone https://github.com/GreatChenLab/deepB3P.git
python3 -m venv .venv-deepb3p
.venv-deepb3p/bin/pip install -r requirements-deepb3p.txt   # (repo 루트에 있음)
cp ../scripts/deepb3p/_run_predict.py ../scripts/deepb3p/_run_fbgan.py deepB3P/

# --- ToxinPred3 (독성) ---
git clone https://github.com/raghavagps/toxinpred3.git
python3 -m venv .venv-toxinpred3
.venv-toxinpred3/bin/pip install -r requirements-toxinpred3.txt
(cd toxinpred3/model && unzip -o toxinpred3.0_model.pkl.zip)  # 모델 압축 해제
```

`core/config.py`가 `vendor/` 경로를 자동 감지해 로컬 엔진을 켠다.

### 3. LLM 에이전트 (선택)
```bash
pip install anthropic
cp .env.example .env   # .env 에 ANTHROPIC_API_KEY(또는 LLM_API_KEY) 입력
```

## 실행
```bash
streamlit run app.py
```

## 구조
```
app.py                  Streamlit UI (thin layer)
core/
  config.py             설정·상수·경로 자동감지, 연결부위 윈도우
  predictors.py         BBB 예측기 (deepB3P local / mock / API)
  toxicity.py           독성 예측기 (ToxinPred3 local / placeholder)
  agent.py              라이브러리 전수 스윕 오케스트레이션
  generative.py         FBGAN 생성 최적화
  optimizer_agent.py    LLM tool-use 분자 최적화 루프
  structure.py          ESMFold 구조 + 셔틀 노출도
  schemas.py            데이터 구조
scripts/deepb3p/        deepB3P venv에서 도는 러너 (설치 시 vendor/deepB3P/ 로 복사)
assets/theme.css        디자인 시스템(YumYum) 테마
```

## 라이선스 주의
`vendor/`의 ToxinPred3는 **GPL-3.0**, deepB3P는 라이선스 미명시. 이들을 재배포하거나
상업적으로 쓰려면 각 저작자·라이선스를 확인할 것. 이 repo에는 포함하지 않음.
