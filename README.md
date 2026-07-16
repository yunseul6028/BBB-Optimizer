# 🧬 BBB-Optimize AI Agent

뇌혈관장벽(BBB)을 통과하는 알츠하이머 치료용 **항체-셔틀 융합 단백질**을 설계하는
**멀티 에이전트 시스템**. 화물(cargo) 펩타이드를 넣으면, **설계 에이전트가 여러 예측 도구를
스스로 골라 써가며** 후보를 만들고, **심사 에이전트가 이를 적대적으로 검증**해 BBB 투과(효능)·
독성(안전성)·입체구조·수용체 결합·안정성을 모두 만족하는 최종 융합체(`화물 + 링커 + 셔틀`)로 수렴한다.

> **공모전 분야 매핑**: 핵심은 **분야2 (도구 활용 기반의 분자 최적화 루프)** — 예측 도구를 루프로
> 돌려 분자를 최적화. 여기에 **설계–심사 멀티 에이전트** 구조로 **분야4 (융합, 멀티 에이전트
> 시스템)** 를, 후보 가설 생성·검증으로 **분야1** 요소를 함께 아우른다.

## 🤖 메인: 멀티 에이전트 자율 설계 (Designer–Critic)

브레인은 **기본 Gemini**(정직한 BBB 프레이밍) / **Claude 선택 가능**.

- **설계(Designer) 에이전트** — 계획 → 도구로 스크리닝·**정밀 설계**·구조 검증·(필요시) 생성 →
  수렴 시 스스로 `finish`로 최종 후보 제출.
- **심사(Critic) 에이전트** — 독립 페르소나·독립 컨텍스트로 제출된 후보를 **적대적으로 검증**
  (독성 마진·안정성·BBB 신뢰도·구조 근거·개발성 liability). **APPROVE**면 확정, **REVISE**면
  설계자가 지적을 반영해 재설계.

| 에이전트 도구 | 하는 일 | 백엔드 엔진 |
|---|---|---|
| `evaluate_candidates` | 후보 융합체의 BBB·독성·안정성·수용체 유사도 배치 평가 | deepB3P + ToxinPred3 + ProtParam + BLOSUM62 |
| `design_candidate` | 라이브러리 밖 **잔기 수준 편집** 서열을 직접 설계·채점 | 위 엔진 재사용 |
| `analyze_structure` | 후보를 접어 셔틀이 표면에 노출됐는지 검증 | ESMFold + Biopython |
| `generate_novel_shuttles` | 라이브러리 밖 신규 셔틀 생성 | FBGAN |
| `finish` | 수렴 판단 시 최종 후보 제출·종료 (자기 종료) | — |

> 각 도구는 앱의 **"🔧 개별 도구 직접 실행"** 패널에서 키 없이 무료로 수동 실행도 가능하다.

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
