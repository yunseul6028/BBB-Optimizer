# 냠냠픽업 — 디자인 시스템

> Figma 파일 [YumYum-v1](https://www.figma.com/design/dOF4N6r4sQpNs068WXudxP/YumYum-v1)에서 추출한 디자인 토큰 및 컴포넌트 가이드

---

## 1. 브랜드 아이덴티티

### 1-1. 컬러 — 브랜드 / 액센트

| 토큰 | HEX | 용도 |
|---|---|---|
| **`main`** | `#ffde36` | **메인 브랜드 컬러** (YumYum Yellow). CTA, 핵심 액센트, 로고 |
| `main_orange` | `#ff8500` | 서브 브랜드 컬러. 포인트, 강조 |
| `Onboarding_pink` | `#ff91c7` | 온보딩 전용 핑크 |

### 1-2. 옐로 스케일

```
#fff6cf  #fff4c9  #ffe6b8  #ffd78c  #fdbc45  #f5c51d  #ffcd1e  #ffde36
─────── light ────────────────────────────────────────────  brand →
```

### 1-3. 오렌지 / 핑크 / 마젠타

| HEX | 용도 |
|---|---|
| `#ff8500` | 메인 오렌지 (`main_orange`) |
| `#ff8a00` / `#eb7501` | 오렌지 변형 |
| `#ff7752` | 살몬 오렌지 |
| `#ff91c7` | 온보딩 핑크 |
| `#ffb8bf` | 라이트 핑크 |
| `#a40037` / `#ba134b` / `#bc2f5b` | 딥 마젠타 / 와인 |

---

## 2. 시스템 컬러

### 2-1. 뉴트럴 (Neutrals)

| 토큰 | HEX | 사용 빈도 | 용도 |
|---|---|---|---|
| `Neutral/Black` | `#000000` | 최다 | 본문 텍스트, 아이콘 |
| `Neutral/Near-Black` | `#080808` | 높음 | 강조 텍스트 |
| `Neutral/900` | `#222222` | 높음 | 헤더 텍스트 |
| `Neutral/800` | `#2d2d2d` | 중간 | 부제 |
| `Neutral/700` | `#333333` | 중간 | 본문 |
| `Neutral/600` | `#464646` | 중간 | 서브 텍스트 |
| `Neutral/500` | `#6b6b6b` | 중간 | 보조 텍스트 |
| `Neutral/400` | `#777777` | 중간 | 비활성 텍스트, `BG-Gray-30/70` 베이스 |
| `Neutral/350` | `#848484` | 중간 | 플레이스홀더 |
| `Neutral/300` | `#999999` | 중간 | 비활성 |
| `Neutral/200` | `#cbcac8` / `#cecece` | 중간 | 디바이더 |
| `Neutral/150` | `#ededed` | 중간 | 라이트 디바이더 |
| `Neutral/100` | `#f7f7f7` | 중간 | 백그라운드 |
| `Neutral/0` | `#ffffff` | 두번째로 많음 | 카드, 시트, 화이트 베이스 |

### 2-2. 알파 변형 (Variable)

| 변수명 | 정의 |
|---|---|
| `BG-Gray-30` | `#777777` @ 30% 투명도 |
| `BG-Gray-70` | `#777777` @ 70% 투명도 |
| `Yellow-cart` | `#000000` @ 20% 투명도 |

---

## 3. 시맨틱 컬러 (Status / System)

### 3-1. Error / Danger

| HEX | 용도 |
|---|---|
| `#eb2323` | Error 텍스트, 경고 아이콘 |
| `#fcdede` | Error 배경 (라이트) |

### 3-2. Success

| HEX | 용도 |
|---|---|
| `#3cde75` | Success 메인 |
| `#99d162` | Success 강조 |
| `#458a4b` | Success 텍스트 |
| `#c5ff9c` | Success 배경 |

### 3-3. Info / Link

| HEX | 용도 |
|---|---|
| `#2597ba` | 인포 메인 |
| `#005878` / `#006699` | 인포 다크 |
| `#1496c9` | 링크 변형 |
| `#c2e3eb` | 인포 배경 |
| `#2763ff` | 링크 블루 |

### 3-4. Food / Earth Tones (음식 관련 톤)

| HEX | 용도 |
|---|---|
| `#bb7b4d` | 브라운 (커피·베이커리) |
| `#8f4d36` | 다크 브라운 |
| `#996246` | 미디엄 브라운 |
| `#c18161` | 라이트 브라운 |
| `#eb894b` | 오렌지 브라운 |
| `#191403` | 딥 커피 |

---

## 4. 타이포그래피

### 4-1. 폰트 패밀리

| 우선순위 | 폰트 | 용도 |
|---|---|---|
| Primary | **Pretendard** | 한글 본문 / UI 전체 |
| Secondary | **SF Pro** | iOS 시스템 글리프 (시간, 상태바) |
| (제거 권장) | Noto Sans, Inter | 일부 잔존 — 정리 필요 |

### 4-2. 타입 스케일 (제안 토큰)

| 토큰 | Family / Weight | Size | 사용 빈도 | 용도 |
|---|---|---|---|---|
| `caption-xs` | Pretendard Regular | 10 | 145회 | 캡션, 메타 |
| `caption` | Pretendard Regular | 12 | **377회 (최다)** | 보조 텍스트, 라벨 |
| `caption-bold` | Pretendard SemiBold | 12 | 다수 | 강조 라벨 |
| `body-sm` | Pretendard Regular | 14 | 251회 | 본문 small |
| `body-sm-medium` | Pretendard Medium | 14 | 다수 | 본문 small 강조 |
| `body-sm-semibold` | Pretendard SemiBold | 14 | 다수 | 강조 |
| `body` | Pretendard Regular | 16 | 다수 | 본문 |
| `body-medium` | Pretendard Medium | 16 | 다수 | 본문 강조 |
| `body-semibold` | Pretendard SemiBold | 16 | 다수 | 본문 강조 (semi) |
| `body-bold` | Pretendard Bold | 16 | 다수 | 본문 굵게 |
| `heading-sm` | Pretendard SemiBold | 17 / 18 | 다수 | 부제목 |
| `heading` | Pretendard Bold | 18 | 다수 | 섹션 제목 |
| `display` | Pretendard Bold | 28 | 적음 | 큰 헤드라인 |

### 4-3. 사이즈 스텝

`10 / 12 / 14 / 16 / 17 / 18 / 28` (px)

---

## 5. 레이아웃 / 스페이싱

### 5-1. 캔버스

- **모바일 baseline:** `390 × 844` (iPhone 13 클래스 portrait)
- **iOS Safe Area:**
  - 상단 Status Bar: `54px` (`Time` 컴포넌트)
  - 하단 Home Indicator: `21px`

### 5-2. 컨테이너 패턴

| 패턴 | 사이즈 | 용도 |
|---|---|---|
| Full Width | `390` | 풀폭 콘텐츠 |
| Side Padding | `20` (좌우) | 화면 좌우 여백 |
| Onboarding Container | `350 × 463` | 온보딩 콘텐츠 박스 |
| Button Large | `295 × 55` | 메인 CTA |
| Bottom Sheet | `390 × 549` (large), `390 × 260` (small), `390 × 226` (sort), `390 × 205` (terms) | 바텀시트 변형 |
| Popup | `300 × 247` (with content), `300 × 161` (empty) | 팝업 |

---

## 6. 컴포넌트 카탈로그

> 총 **641개**의 로컬 컴포넌트 / 406개 고유 이름

### 6-1. 화면 셸 (Shell)

| 컴포넌트 | 용도 |
|---|---|
| `Basic_screen` / `Base_screen` | 모든 화면의 컨테이너 |
| `Base_screen/Bottom sheet` | 바텀시트 셸 |
| `Base_screen/Home Indicator` | 홈 인디케이터 (21px) |
| `Time` | iOS 상태바 시계 (54px) |
| `Bottom Tab Bar` | 5개 탭 (Home / Search / Favorite / Orders / Profile) |
| `Title Bar` | 페이지 타이틀 |
| `Header style=Default·Search·Typing` | 헤더 변형 |
| `Home header` | 홈 전용 헤더 |
| `Navigation bar` (`320`/`375`/`+375`) | 일반 네비게이션 |

### 6-2. 버튼

| 컴포넌트 | 사이즈 | 용도 |
|---|---|---|
| `Button - Large` | 295 × 55 | 메인 CTA |
| `Button - Small` | — | 보조 액션 |
| `Bottom bnt` | — | 화면 하단 고정 |
| `Under_Bnt` | — | 텍스트 링크 |
| `Start Yumyum btn` | — | 시작 버튼 |
| `Popup_Button` | — | 팝업 내부 |
| `Large/Small button - Icon frame` | — | 아이콘 + 텍스트 |
| `Btn` | — | 범용 |

### 6-3. 입력 (Input)

| 컴포넌트 | 용도 |
|---|---|
| `Input field` | 기본 텍스트 입력 |
| `Input field message icon` | 메시지 아이콘 포함 |
| `Input filed function icon` | 기능 아이콘 포함 |
| `Input filed indicator icon` | 인디케이터 포함 |
| `Login / input` (variant set) | 로그인 입력 |
| `Copoun input` | 쿠폰 입력 |

### 6-4. 선택 컨트롤

| 컴포넌트 | Variants |
|---|---|
| `Check box` | `Checked=Yes/No` |
| `Radio btn` | — |
| `Toggle` | `On` / `Off`, `Selected=Yes/No`, `Show=Yes/No` |
| `Quantity selector` | 수량 선택 (- / +) |
| `Duration selector` | 시간 선택 |
| `Dropdown` | `status=Open/Closed` |

### 6-5. 카드 / 콘텐츠

| 컴포넌트 | 용도 |
|---|---|
| `Restaurant card` | 매장 카드 |
| `Restaurant image` | 매장 이미지 |
| `store card - home` (358 × 195) | 홈 매장 카드 |
| `Home screen / Resauratn card` | 홈 매장 카드 변형 |
| `Favorite restaurant` | 즐겨찾기 카드 |
| `Meal card` / `Meals` / `Meals summary` | 메뉴 카드 |
| `Order card` / `Order flow` / `Order status` | 주문 카드 |
| `Orders summary` | 주문 요약 |
| `Notification card` | 알림 카드 |
| `Banner` | 배너 |
| `List / Categories` | 카테고리 리스트 |
| `Screen tabs` | 스크린 탭 |

### 6-6. 카트 / 쿠폰

| 컴포넌트 | Variants |
|---|---|
| `Add to cart` / `Add to cart section` | — |
| `Cart` | `=Add to / Remove from / Disabled` |
| `Coupon` | `=Added / Not added` |
| `Ordered meal` | 주문된 메뉴 |

### 6-7. 매장 상태 (Restaurant Status)

| Variant | 의미 |
|---|---|
| `Open now` | 영업 중 |
| `Closed` | 영업 종료 |
| `Opnenig time` | 오픈 시간 표시 |
| `Receive orders at` | 주문 접수 시간 |
| `Restaurant call` | 전화 주문 |

### 6-8. 주문 상태 머신

**고객 측 (User App)**
```
pending → recieved → Accepted → Preparing → Done
```

**사장 측 (Restaurant App)**
```
New order → Active
              ├─ Paused
              ├─ Busy
              └─ Closed
```

추가 상태: `Replied`, `started`, `Middle state`

### 6-9. 지도 / 위치

| 컴포넌트 | 용도 |
|---|---|
| `Map_pin-open` / `Map_pin-close` | 매장 핀 (영업/마감) |
| `MiniMap_pin-open` / `MiniMap_pin-close` | 미니맵 |
| `Map_view-Button-On/Off` | 지도/리스트 토글 |
| `Map_view-open` / `Map_view-close` | 지도 뷰 상태 |
| `Map_view-open-Discount` | 할인 강조 지도 |
| `Current location` | 현재 위치 |
| `Location` | `=Default / Detect location` |
| `Mappin` | 일반 핀 |
| `Saved addresses bottom sheet` | 저장 주소 시트 |
| `Address type` | 주소 타입 |

### 6-10. 피드백 / 오버레이

| 컴포넌트 | 용도 |
|---|---|
| `Toast` | 토스트 메시지 |
| `Snackbar content` | 스낵바 |
| `Popup` (`1-line` / `2-line` / `3-pic` / `NO-text`) | 팝업 변형 |
| `Popup_2` | 팝업 변형 |
| `Alert / 현재위치` (140 × 35) | 위치 알림 |

### 6-11. 리뷰 / 평점

| 컴포넌트 | 용도 |
|---|---|
| `Review` / `Review card` / `Review card - client` | 리뷰 |
| `Bussiness owner reply` | 사장님 답글 |
| `Status` (`=Add reply / Replied / Reply form`) | 답글 상태 |
| `Grade` / `Grade_icon` / `Star` | 별점 |

### 6-12. 아바타 / 마스코트

**마스코트 캐릭터** — 브랜드의 핵심 자산
- Yami (냐미) · Yamu (냐무) · Kumo (쿠모) · Pong-ji (빵지) · Hong-ji · Chila (칠라)
- 추가: 무무 / 무지 / 상수 / 흰둥이 / 효희

| 컴포넌트 | 용도 |
|---|---|
| `User avatar` | 일반 사용자 아바타 |
| `User avatar / 사장님` | 사장님 아바타 |
| `Avatar / Image-60` | 60px 아바타 |
| `Charater` | 캐릭터 일러스트 |

### 6-13. 카테고리 칩 (음식 분류)

한식 · 중식 · 일식 · 아시안 · 회 · 고기 · 치킨 · 피자 · 햄버거 · 김밥분식 · 도시락 · 샌드위치 · 샐러드 · 커피음료 · 빵지

### 6-14. 아이콘 시스템

`Icon=` variant — 약 40종

**네비게이션:** Back, Cancel, Down, Up, Right, Left, More
**액션:** Edit, Delete, Plus, Minus, Reload, Reorder, Search, Upload, Download
**상태:** Success, Error, Hint, Hide, View, Mandatory, Maintenance, Repair
**컨텐츠:** Cart, Cart out, Calendar, Date, Time, Won, Flag, Email, Phone number
**개인:** User - Profile, ID, Date of birth, Gender, Lock, Lock - Password, Passport, Language
**위치:** Location, User location
**기타:** Checkout, WhatsApp, Hide, Assign to unit, Change accommodation

### 6-15. 빈/예외 상태

| 컴포넌트 | 용도 |
|---|---|
| `Empty` | 빈 상태 |
| `Empty_icon` | 빈 상태 아이콘 |
| `Warring_icon` | 경고 아이콘 |
| `X-2` | 닫기 |

---

## 7. 인터랙션 패턴

### 7-1. CTA 버튼 위계

1. **메인 CTA** — `Button - Large` (295 × 55, 노란색 `#ffde36` 채움)
2. **하단 고정 CTA** — `Bottom bnt`
3. **보조 액션** — `Button - Small`
4. **텍스트 링크** — `Under_Bnt`

### 7-2. 바텀시트 (Bottom Sheet)

용도별 높이:
- **약관 동의**: 549 (long form), 205 (short)
- **위치 선택**: 260
- **정렬 옵션**: 226

### 7-3. 팝업 (Popup)

- 표준: 300 × 247
- Empty 상태: 300 × 161

### 7-4. 토스트

- 화면 하단에서 슬라이드 업
- 자동 dismiss

---

## 8. 디자인 토큰 정리 권장사항

> 현재 Figma 파일은 로컬 토큰이 매우 적게 정의되어 있음 (paint style 1개, variable 10개). 아래 정리를 권장합니다.

### 8-1. Variable Collection (제안)

```
Color/Brand
  ├─ main          → #ffde36
  ├─ main-orange   → #ff8500
  └─ onboarding-pink → #ff91c7

Color/Neutral
  ├─ 0 (white)     → #ffffff
  ├─ 100           → #f7f7f7
  ├─ 150           → #ededed
  ├─ 200           → #cecece
  ├─ 300           → #999999
  ├─ 400           → #777777
  ├─ 500           → #6b6b6b
  ├─ 600           → #464646
  ├─ 700           → #333333
  ├─ 800           → #2d2d2d
  ├─ 900           → #222222
  └─ 1000 (black)  → #000000

Color/Semantic
  ├─ error         → #eb2323
  ├─ error-bg      → #fcdede
  ├─ success       → #3cde75
  ├─ info          → #2597ba
  └─ link          → #2763ff

Color/Alpha
  ├─ bg-gray-30    → #777777 @ 30%
  ├─ bg-gray-70    → #777777 @ 70%
  └─ yellow-cart   → #000000 @ 20%

Spacing
  ├─ xs            → 4
  ├─ sm            → 8
  ├─ md            → 16
  ├─ lg            → 20
  ├─ xl            → 24
  └─ 2xl           → 32

Radius
  ├─ sm            → 4
  ├─ md            → 8
  ├─ lg            → 12
  └─ full          → 10000

Typography (Text Style)
  ├─ display       → Pretendard Bold 28
  ├─ heading       → Pretendard Bold 18
  ├─ heading-sm    → Pretendard SemiBold 17
  ├─ body          → Pretendard Regular 16
  ├─ body-sm       → Pretendard Regular 14
  └─ caption       → Pretendard Regular 12
```

### 8-2. 정리가 필요한 항목

- [ ] **Text Style 부재** — 폰트 사이즈/웨이트가 raw value로 적용됨. 위 토큰으로 등록 권장
- [ ] **Color Style 1개** — `main` 외에 시멘틱 컬러 미정의. 시스템 컬러를 paint style 또는 variable로 등록 권장
- [ ] **Component naming 정리** — `속성 1=…` 같은 placeholder, `Resauratn` 오타, `Opnenig`/`Bussiness` 오타 등
- [ ] **Property 1=…` 익명 variant** → 의미 있는 property 이름으로 변경

---

## 9. 참고

- 디자인 파일: [YumYum-v1 (Figma)](https://www.figma.com/design/dOF4N6r4sQpNs068WXudxP/YumYum-v1)
- BM Flow 보드: [YumYum BM Flow / 기획 (FigJam)](https://www.figma.com/board/rIh4oJKJMEVByNAh0VFc3m/YumYum-BM-Flow--%EA%B8%B0%ED%9A%8D)
- 서비스 페이지: https://www.yumyum.im
- 화면 명세서: `02-screen-specs.md`
- 개발자 핸드오프: `04-developer-handoff.md`
