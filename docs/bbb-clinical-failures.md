# BBB 투과 병목과 CNS 임상 실패 — 통계·사례 문헌 정리

> 본 프로젝트(뇌 전달 셔틀 융합체 설계)의 **문제 정의 근거**. 혈액뇌장벽(BBB)이 왜 뇌질환
> 신약개발의 병목인지, 통계와 실제 실패 사례로 정리한다.
>
> ⚠️ **정직한 구분**: "BBB를 못 넘는다"(투과 실패)와 "임상에서 떨어진다"(임상 실패)는 **다른 지표**다.
> 임상 실패는 다요인(효능·안전성·BBB 전달)이며, BBB는 그중 **주요 기여 요인**이지 유일 원인이 아니다.
> 아래 수치는 출처·정의가 제각각이므로 원문 대조 권장.

---

## 1. BBB 투과 병목 — "애초에 뇌에 못 들어간다"
**Pardridge (2005), *NeuroRx* "The blood-brain barrier: bottleneck in brain drug development":**
- **대분자 신경치료제(항체·단백질 등)의 ~100%**가 BBB를 통과하지 못함
- **소분자 약물의 >98%**가 BBB를 통과하지 못함
- 고속 스크리닝(HTS)에서 나온 **후보의 >98%가 BBB 미통과** → in vitro 활성이 좋아도 뇌 내 활성이 없어 프로그램이 종료됨

> ※ 이 98%는 **투과 실패율**이지 임상 실패율이 아님. 단일 저자(Pardridge)의 추정으로, 독립적
> 체계 연구 확인은 제한적. 그래도 "BBB가 1차 관문"이라는 정성적 근거로 널리 인용됨.

**항체 전달 효율**: 전신 투여 항체(~150 kDa)는 **투여량의 약 0.1–0.2%만 뇌 실질에 도달**(RMT 셔틀 기술의 동기).

---

## 2. CNS 임상 소모율(attrition) — 다른 질환보다 훨씬 높음
- **CNS 약물 개발 소모율 ≈ 91%**, 특히 **Phase II 성공률이 매우 낮음** (Leiden Univ. 학위논문)
- CNS 약은 시판까지 **12–16년** (비CNS 10–12년)으로 더 오래 걸리고 소모율도 높음
- **CNS/정신과 질환의 LOA(Phase I→승인 확률) ≈ 6.2%** — 종양학(5.1%) 다음으로 **최하위권** (BIO/Biomedtracker 2006–2015, Amplion 요약). 혈액학이 최고. CNS·종양은 **프로그램 수가 많아** 전체 평균을 끌어내림
- **Phase II가 최대 관문**: 전체 평균 성공확률 ~28% (Biomedtracker 2014–2023)

---

## 3. 알츠하이머 — 99.6% 실패
**Cummings, Morstorf & Zhong (2014), *Alzheimer's Research & Therapy* "few candidates, frequent failures" (PMID 25024750):**
- ClinicalTrials.gov 2002–2012 분석: **Phase 전이 성공률 0.4% → 실패율 99.6%**
- 10년간 **244개 약물만** 시험됨(질병 부담 대비 극소), 2009년 이후 감소
- 성공 사례는 사실상 **1개(메만틴, 2003)** — 질병조절제(disease-modifying)는 전무

> ※ 99.6%는 "Phase 전이 기준 소모율"이지 BBB 단독 실패율이 아님. 단 AD 치료제 다수가 항체·대분자라
> **전달(BBB)이 핵심 장벽** 중 하나로 반복 지적됨.

---

## 4. 실제 실패 사례 — "뇌에 도달하지 못해서"
전달 자체가 실패 원인으로 지목된 대표 사례 (신경영양인자류):

| 약물 | 적응증 | 결과 | 전달 실패 근거 |
|---|---|---|---|
| **GDNF** | 파킨슨병 | Phase 2 1차 지표(UPDRS 25%↑) **실패**; 뇌실 투여도 위약과 차이 없음 | 카테터 주변 확산이 **인간 퍼타멘의 2–9%**에만 도달 — 조직에 충분히 못 퍼짐 |
| **NGF** | 콜린성 뉴런(AD형) | 뇌실(ICV) 투여 **실패** | CSF→뇌실질 **확산 속도 극히 느림(~10⁻⁶ cm/s)** |
| **BDNF** | ALS | 1,000명+ 투여(피하·척수강), **효과 없음** | 피하 투여는 **BBB 투과 불량 + 짧은 혈장 반감기** |

→ 공통 교훈: **분자의 생물학이 아니라 "뇌 표적까지 전달"이 실패 원인.** 그래서 BBB 통과 기술(셔틀)이 필요.

---

## 5. 설계 trade-off — "친유성 올리면 독성 올라간다"
Bayer 계열 분석(119개 약물·108개 후보): BBB 투과를 높이려 **친유성↑·극성↓**하는 전략은 **in vivo 독성 위험이 높은 화학공간**으로 밀어붙여, **투과 이득이 안전성과 상충**함.
→ 본 프로젝트가 "양이온·친유성으로 뚫는 CPP"보다 **수용체 선택적 RMT**를 선호하고, **선택성·개발성 축으로 off-target/독성을 감점**하는 근거.

---

## 6. 본 프로젝트와의 연결
- **문제**: 뇌질환 약은 (1) 애초에 BBB를 못 넘고(>98%), (2) CNS 임상 소모율이 최고 수준(AD 99.6%)이며, (3) 실제로 "뇌에 도달 못 해" 실패한 사례(GDNF·NGF·BDNF)가 많다.
- **해결 방향**: 검증된 **RMT 셔틀**에 화물을 융합해 BBB를 "편승"하되, 친유성·양이온으로 뚫는 CPP의 **off-target·독성 상충**은 선택성·개발성 축으로 피한다.
- **우리 시스템의 자리**: 비싼 임상·실험 전에 **BBB 전달 가능성·선택성·개발성을 in-silico로 정직하게 선별**한다(실험 대체 아님, 우선순위화).

---

## 참고문헌 (원문 대조 권장)
- Pardridge WM. *The blood-brain barrier: bottleneck in brain drug development.* NeuroRx. 2005. — [PMC539316](https://pmc.ncbi.nlm.nih.gov/articles/PMC539316)
- Cummings JL, Morstorf T, Zhong K. *Alzheimer's disease drug-development pipeline: few candidates, frequent failures.* Alzheimer's Res Ther. 2014 (PMID 25024750). — [PMC4095696](https://pmc.ncbi.nlm.nih.gov/articles/4095696) · [FierceBiotech 99.6%](https://www.fiercebiotech.com/r-d/alzheimer-s-r-d-suffers-as-trial-failure-rate-hits-an-astounding-99-6)
- BIO/Biomedtracker/QLS. *Clinical Development Success Rates.* — [go.bio.org](https://go.bio.org/clinicaldevelopment) · [Amplion 요약(CNS 6.2%)](https://www.amplion.com/?p=833)
- 신경영양인자 전달 실패 리뷰. *Front. Aging Neurosci.* 2021. — [fnagi.2021.645583](https://www.frontiersin.org/journals/aging-neuroscience/articles/10.3389/fnagi.2021.645583/full) · [GDNF (Michael J. Fox Foundation)](https://www.michaeljfox.org/news/nature-news-crossing-barrier)
- CNS 소모율·개발기간: Leiden Univ. 학위논문; CNS PK 리뷰(Reichel 2009 등)
