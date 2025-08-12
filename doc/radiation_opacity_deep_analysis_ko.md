# AMRVAC의 복사 모듈 및 불투명도에 대한 심층 분석

## 요약

이 문서는 MPI-AMRVAC 코드의 복사 전달 구현 및 불투명도 계산에 대한 포괄적인 분석을 제공합니다. 특히 OPAL 불투명도 테이블과 플럭스 제한 확산(Flux-Limited Diffusion, FLD) 프레임워크 내에서의 통합에 중점을 둡니다. 이 분석은 `mod_fld.t`, `mod_opal_opacity.t`, `mod_cak_opacity.t` 및 관련 모듈을 포함한 소스 코드의 철저한 검토를 기반으로 합니다.

## 목차

1. [서론](#서론)
2. [복사 전달 프레임워크](#복사-전달-프레임워크)
3. [OPAL 불투명도 테이블](#opal-불투명도-테이블)
4. [불투명도 유형과 물리적 과정](#불투명도-유형과-물리적-과정)
5. [구현 세부사항](#구현-세부사항)
6. [복사-물질 결합](#복사-물질-결합)
7. [계산 워크플로우](#계산-워크플로우)
8. [물리적 영역과 응용](#물리적-영역과-응용)
9. [참고문헌 및 추가 자료](#참고문헌-및-추가-자료)

## 서론

AMRVAC(Adaptive Mesh Refinement Versatile Advection Code)는 다양한 복사 전달 영역을 처리하도록 설계된 여러 모듈을 통해 정교한 복사 유체역학 기능을 구현합니다. 이 코드는 세 가지 주요 복사 전달 시나리오를 다룹니다:

1. **광학적으로 두꺼운 전달**: 플럭스 제한 확산(FLD)을 통해
2. **선 구동 바람(Line-driven winds)**: CAK 이론 사용
3. **광학적으로 얇은 냉각**: 표로 작성된 냉각 곡선을 통해

이 문서는 주로 FLD 구현과 불투명도 테이블, 특히 OPAL(Opacity Project At Livermore) 테이블의 사용에 초점을 맞춥니다.

## 복사 전달 프레임워크

### 플럭스 제한 확산(FLD) 모듈

FLD 모듈(`mod_fld.t`)은 광자 평균 자유 경로가 시스템의 특성 스케일보다 훨씬 작은 광학적으로 두꺼운 영역에 적합한 복사 전달을 위한 플럭스 제한 확산 근사를 구현합니다. 이 모듈은 Turner & Stone (2001)을 기반으로 개발되었으며 Moens et al. (2022, A&A 657)에 자세히 설명되어 있습니다.

#### 핵심 구성 요소

1. **복사 에너지 밀도 진화**
   ```fortran
   ∂E_rad/∂t + ∇·F_rad = -cκρ(E_rad - aT^4)
   ```
   여기서:
   - `E_rad`: 복사 에너지 밀도
   - `F_rad`: 복사 플럭스
   - `κ`: 불투명도 (cm²/g)
   - `ρ`: 질량 밀도
   - `T`: 온도
   - `a`: 복사 상수
   - `c`: 빛의 속도

2. **플럭스 제한자 함수**
   복사 플럭스는 다음과 같이 계산됩니다:
   ```fortran
   F_rad = - (c λ(R) / (κ ρ)) ∇E_rad
   ```
   여기서 λ(R)는 다음 사이를 전환하는 플럭스 제한자입니다:
   - 광학적으로 두꺼운 확산: λ = 1/3
   - 자유 스트리밍: λ = 1/R
   - R = |∇E_rad|/(κρE_rad)는 무차원 기울기

### 사용 가능한 플럭스 제한자

AMRVAC는 여러 플럭스 제한자 처방을 구현합니다:
- **Pomraning**: 부드러운 전환을 가진 기본 제한자
- **Diffusion**: 순수 확산 한계 (λ = 1/3)
- **FreeStream**: 전체 R 의존 제한자를 고려

## OPAL 불투명도 테이블

### 테이블 구조와 형식

OPAL 테이블은 온도와 밀도의 함수로 Rosseland 평균 불투명도를 제공합니다. 테이블은 다음과 같이 구성됩니다:

#### 매개변수 공간
- **온도 범위**: log₁₀(T) ∈ [3.75, 8.70]
  - 물리적 범위: ~5,600 K에서 ~5×10⁸ K
  - 표준 테이블에서 70개의 온도 점

- **밀도 매개변수**: log₁₀(R) ∈ [-8.0, 1.0]
  - R = ρ/(T/10⁶ K)³ (CGS 단위)
  - 온도당 19개의 밀도 점

#### 테이블 형식 예시
```
TABLE # 8  X=0.0000 Y=0.9800 Z=0.0200  (울프-레이에 대기)

                    log R
logT  -8.0  -7.5  -7.0  ...  0.5   1.0
3.75  9.999 9.999 -3.233 ... -1.806 -1.346
3.80  9.999 9.999 -3.298 ... -1.687 -1.120
...
```

여기서:
- `X`: 수소 질량 분율
- `Y`: 헬륨 질량 분율
- `Z`: 금속 질량 분율
- 값은 cm²/g 단위의 log₁₀(κ)
- `9.999`: 이러한 조건에 대한 데이터 없음

### R 매개변수의 물리적 해석

밀도 매개변수 R = ρ/(T/10⁶)³는 다음의 결합된 효과를 포착합니다:
1. **흡수체의 수밀도**: ρ에 비례
2. **이온화 상태**: 강한 온도 의존성
3. **복사압 스케일링**: 항성 내부에서의 T³ 관계

이 매개변수화가 특히 유용한 이유:
- 2차원 (ρ,T) 공간을 더 컴팩트한 표현으로 줄임
- 복사압이 지배적인 항성 내부의 물리를 자연스럽게 포착
- 넓은 매개변수 범위에서 부드러운 보간 제공

### AMRVAC에서의 구현

#### 모듈 구조 (`mod_opal_opacity.t`)

```fortran
module mod_opal_opacity
  ! 테이블 차원
  integer, parameter :: iRmin = 2, iRmax = 20
  integer, parameter :: iTmin = 7, iTmax = 76
  
  ! 전역 저장소
  double precision, public :: Kappa_vals(iTmin:iTmax,iRmin:iRmax)
  double precision, public :: logR_list(iRmin:iRmax)
  double precision, public :: logT_list(iTmin:iTmax)
  
  public :: init_opal_table
  public :: set_opal_opacity
end module
```

#### 불투명도 조회 과정

1. **초기화** (한 번 수행):
   ```fortran
   call init_opal_table(tablename)
   ! AMRVAC_DIR/src/rhd/OPAL_tables/에서 테이블 읽기
   ```

2. **런타임 조회**:
   ```fortran
   ! OPAL 관례로 변환
   logR_in = log10(rho/(temp*1d-6)**3)
   logT_in = log10(temp)
   
   ! 이중선형 보간
   call get_kappa(Kappa_vals, logR_list, logT_list, 
                  logR_in, logT_in, logKappa_out)
   
   ! 물리적 단위로 변환
   kappa = 10**logKappa_out
   ```

3. **보간 알고리즘**:
   - (logR, logT) 공간에서 주변 네 개의 격자점 식별
   - log-log 공간에서 이중선형 보간 수행
   - 경계 밖 값에 대한 외삽 처리
   - 누락된 데이터(9.999 값)에 대한 특별 처리

## 불투명도 유형과 물리적 과정

### FLD에서 사용 가능한 불투명도 법칙

AMRVAC의 FLD 모듈은 `fld_opacity_law` 매개변수를 통해 여러 불투명도 처방을 지원합니다:

#### 1. 상수 불투명도 (`'const'`)
```fortran
fld_kappa = fld_kappa0/unit_opacity
```
- 간단하고 계산적으로 효율적
- 테스트 및 이상화된 문제에 유용

#### 2. 톰슨 산란 (`'thomson'`)
```fortran
! 전자 산란(톰슨); κ_es [cm^2 g^-1] = (σ_T/m_p)*(X + 0.5*Y) ≈ 0.2*(1+X)
sigma_thomson = 6.6524587321d-25    ! cm^2 (CGS)
m_p_cgs       = 1.67262192369d-24   ! g (CGS)
kappa_cgs     = (sigma_thomson/m_p_cgs) * (X + 0.5d0*Y)   ! cm^2/g (CGS)
fld_kappa     = kappa_cgs / unit_opacity                  ! 코드 단위
```In doc/radiation_opacity_deep_analysis_ko.md around lines 623 to 625, remove the AI persona authorship line ("작성자: Claude에 의한 기술 분석") and replace it with a neutral authorship or project/team attribution consistent with repository conventions (e.g., "작성자: [프로젝트 팀명]" or omit the authorship line entirely if the repo does not track per-page authors); ensure the surrounding metadata formatting (version and last updated) remains intact and matches the repo's docs style.
- 순수 전자 산란
- 주파수 독립적
- 뜨겁고 완전히 이온화된 플라즈마에서 지배적
  - 완전 이온화된 H/He 혼합물 가정에서 빠른 근사: κ_es ≈ 0.2*(1+X) [cm^2 g^-1]
  - 코드 단위 변환은 `unit_opacity = (unit_length_cm^2 / unit_mass_g)` 정의를 사용하므로, `fld_kappa = κ_cgs / unit_opacity`로 처리됨

#### 3. 크라머스 불투명도 (`'kramers'`)
```fortran
fld_kappa ∝ ρ * T^(-3.5)
```
- 자유-자유(제동복사) 흡수 근사
- 뜨겁고 이온화된 플라즈마에 유효
- 간단한 해석적 형태

#### 4. OPAL 테이블 (`'opal'`)
```fortran
case('opal')
  call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)
  rho0 = w(ix^D,iw_rho)*unit_density
  Temp0 = Temp(ix^D)*unit_temperature
  call set_opal_opacity(rho0,Temp0,kappa)
  fld_kappa(ix^D) = kappa/unit_opacity
```
- 가장 정교한 옵션
- 모든 불투명도 소스 포함
- 온도와 밀도 의존

#### 5. 비등온 (`'non_iso'`)
```fortran
fld_kappa ∝ (ρ/ρ₀) * (T/T₀)^n
```
- 멱법칙 온도 의존성
- 조정 가능한 지수 n (일반적으로 -3.5)

#### 6. FastWind (`'fastwind'`)
- 항성풍 응용에 특화
- 온도 의존 향상 포함

#### 7. 특수/사용자 정의 (`'special'`)
```fortran
call usr_special_opacity(ixI^L, ixO^L, w, x, fld_kappa)
```
- 사용자 정의 불투명도 구현 허용

### CAK 선 불투명도

CAK 모듈(`mod_cak_opacity.t`)은 선 구동 항성풍을 위한 특수 불투명도를 제공합니다:

```fortran
module mod_cak_opacity
  ! CAK 매개변수 (Gayley 1995 표기법)
  double precision :: alpha_vals    ! 멱법칙 지수
  double precision :: Qbar_vals     ! 평균 선 강도
  double precision :: Q0_vals       ! 연속체 정규화
  double precision :: kappae_vals   ! 전자 산란
```

이러한 매개변수는 다음을 특성화합니다:
- **선 분포**: α는 선 강도의 분포를 설명
- **선 힘 승수**: M(t) = k*t^(-α) 여기서 t는 광학 깊이 매개변수
- **연속체 대 선 불투명도**: 비율은 바람 가속 효율성을 결정

### 불투명도에 기여하는 물리적 과정

#### 온도 영역

1. **저온 (T < 10⁴ K)**
   - 분자 흡수 (H₂, H₂O, CO 등)
   - 원자 선 흡수
   - 먼지 입자 (존재하는 경우)

2. **중간 온도 (10⁴ K < T < 10⁶ K)**
   - 속박-자유 전이 (광이온화)
   - 속박-속박 전이 (선 흡수)
   - 부분적으로 이온화된 금속 (특히 철-피크 원소)
   - 10⁵ K 근처의 철 불투명도로 인한 "Z-범프"

3. **고온 (T > 10⁶ K)**
   - 자유-자유 흡수 (제동복사)
   - 전자 산란 (톰슨/콤프턴)
   - 완전히 이온화된 플라즈마

#### 조성 의존성

불투명도는 화학 조성에 강하게 의존합니다:
- **수소 (X)**: 산란을 위한 전자 제공
- **헬륨 (Y)**: 단위 질량당 수소보다 적은 불투명도
- **금속 (Z)**: 특히 중간 온도에서 불투명도를 극적으로 증가

## 복사-물질 결합

### 에너지 교환 메커니즘

복사와 물질 사이의 에너지 상호작용은 다음에 의해 지배됩니다:

```fortran
dE_gas/dt = cκρ(E_rad - aT⁴)
```

이는 AMRVAC가 다양한 방법을 사용하여 해결하는 결합된 시스템으로 이어집니다:

#### 상호작용 방법

1. **즉각적 평형**
   - 즉각적인 열평형 가정
   - E_gas + E_rad = 상수
   - 평형 온도 해결

2. **암시적 시간 적분**
   - 4차 다항식 해결:
   ```fortran
   e_gas⁴ + c₁*e_gas - c₀ = 0
   ```
   여기서:
   - `c₁ = (1 + a₂)/a₁`
   - `c₀ = ((1 + a₂)*e_gas + a₂*E_rad)/a₁`
   - `a₁ = 4κσ_B(γ-1)⁴/ρ³ Δt`
   - `a₂ = cκρΔt`

#### 근 찾기 방법

AMRVAC는 에너지 다항식을 해결하기 위해 세 가지 방법을 구현합니다:

1. **이분법 (Bisection Method)**
   ```fortran
   subroutine Bisection_method(e_gas, E_rad, c0, c1)
     ! 견고하지만 느린 수렴
     ! 브래킷 내에서 근을 찾는 것이 보장됨
   ```

2. **뉴턴-랩슨 방법 (Newton-Raphson Method)**
   ```fortran
   subroutine Newton_method(e_gas, E_rad, c0, c1)
     ! 더 빠른 수렴 (이차)
     ! 초기 추정이 좋지 않으면 실패할 수 있음
   ```

3. **할리 방법 (Halley's Method)**
   ```fortran
   subroutine Halley_method(e_gas, E_rad, c0, c1)
     ! 삼차 수렴
     ! 이 문제에 대해 뉴턴보다 안정적
   ```

### 복사력

물질에 작용하는 복사에 의한 체적 힘 밀도는 다음과 같습니다:

```fortran
f_rad = (κ_F ρ / c) * F_flux
```

- 여기서 `κ_F`는 플럭스 평균 불투명도, `ρ`는 질량 밀도, `F_flux`는 복사 플럭스(`F_flux ≡ F_rad`), `c`는 빛의 속도입니다.

질량당 가속도가 필요한 경우는 다음과 같습니다:

```fortran
a_rad = f_rad / ρ = (κ_F / c) * F_flux
```

아래의 항목들은 `f_rad`(힘 밀도) 또는 이에 상응하는 `a_rad`(질량당 가속도)가 하는 역할을 의미합니다:

이 힘/가속도는:
- 항성풍을 가속
- 복사압 불안정성을 구동
- 복사가 지배하는 영역에서 중력에 대항하여 별을 지지

### 확산 계수

FLD 근사에서의 복사 확산 계수:

```fortran
D = c*λ/(κρ)

이 계수는:
- 복사 에너지 전달 속도를 제어
- 불투명도(κ)와 플럭스 제한자(λ)에 모두 의존
- 광학적으로 두꺼운 영역(R→0): λ→1/3 이므로 D→ c/(3κρ)
- 광학적으로 얇은 영역(R→∞): λ→1/R 이고, 유효 확산 계수 개념은 한계가 있으며 FLD는 |F| ≤ c E의 포화(flux-limiting)로 전환됨
## 계산 워크플로우

### 초기화 단계

1. **매개변수 파일 읽기**
   ```fortran
   namelist /fld_list/ fld_kappa0, fld_opacity_law, 
                       fld_opal_table, fld_fluxlimiter, ...
   ```

2. **불투명도 테이블 초기화** (OPAL 사용 시)
   ```fortran
   if (fld_opacity_law .eq. 'opal') then
     call init_opal_table(fld_opal_table)
   endif
   ```

3. **멀티그리드 솔버 설정** (암시적 확산용)
   ```fortran
   if (fld_diff_scheme .eq. 'mg') then
     use_multigrid = .true.
     phys_implicit_update => Diffuse_E_rad_mg
   endif
   ```

### 런타임 실행

#### 시간 단계별 작업

1. **불투명도 필드 계산**
   ```fortran
   call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
   ```

2. **플럭스 제한자 계산**
   ```fortran
   call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, R)
   ```

3. **복사 에너지 업데이트** (암시적 확산)
   ```fortran
   call Diffuse_E_rad_mg(dt, ixI^L, ixO^L, w, wCT, x)
   ```

4. **에너지 상호작용** (연산자 분할이 아닌 경우)
   ```fortran
   call Energy_interaction(w, wCT, x, ixI^L, ixO^L)
   ```

5. **복사력 적용** (연산자 분할이 아닌 경우)
   ```fortran
   call get_fld_rad_force(qdt, ixI^L, ixO^L, wCT, w, x)
   ```

### 성능 최적화

1. **테이블 캐싱**: OPAL 테이블은 초기화 시 한 번 로드
2. **플럭스 제한자 필터링**: 기울기를 부드럽게 하기 위한 선택적 이동 평균
3. **확산 계수 필터링**: 수치적 노이즈 감소
4. **멀티그리드 가속**: 효율적인 암시적 확산 솔버
5. **연산자 분할**: 경직된 항의 별도 처리

## 물리적 영역과 응용

### 항성 대기

OPAL 불투명도는 다음에 특히 적합합니다:
- **항성 내부**: 복사압이 중요한 곳
- **항성 외피**: 대류에서 복사 전달로의 전환
- **울프-레이에 바람**: 고온, 헬륨이 풍부한 환경

울프-레이에 별의 예시 구성:
```fortran
&fld_list
  fld_opacity_law = 'opal'
  fld_opal_table = 'Y09800'  ! Y=0.98, Z=0.02
  fld_fluxlimiter = 'Pomraning'
  fld_diff_scheme = 'mg'
/
```

### 강착 원반

적절한 불투명도를 가진 FLD는 다음을 모델링할 수 있습니다:
- **원반 수직 구조**: 복사압 지지
- **열적 불안정성**: S-곡선 거동
- **복사 구동 유출**: 내부 원반 영역에서

### 초신성 잔해

고온 충격파의 경우:
```fortran
&fld_list
  fld_opacity_law = 'thomson'  ! 전자 산란이 지배적
  fld_fluxlimiter = 'Pomraning'
  fld_interaction_method = 'Halley'
/
```

### 실험실 천체물리학

제어된 실험을 위한 단순화된 불투명도:
```fortran
&fld_list
  fld_opacity_law = 'const'
  fld_kappa0 = 0.34d0  ! cm²/g
/
```

## 고급 기능

### 선-힘 불투명도

CAK 유형 선 구동 바람의 경우, AMRVAC는 공간적으로 변하는 선-힘 매개변수를 사용할 수 있습니다:

```fortran
if (lineforce_opacities) then
  ! 방향 불투명도를 위한 배열 할당
  allocate(i_opf(ndim))
  do idir = 1,ndim
    i_opf(idir) = var_set_extravar('k'//ind_1, 'k'//ind_1)
  enddo
endif
```

### 사용자 정의 불투명도

사용자 모듈을 통한 사용자 정의 불투명도 구현:

```fortran
subroutine usr_special_opacity(ixI^L, ixO^L, w, x, kappa)
  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
  double precision, intent(out) :: kappa(ixO^S)
  
  ! 사용자 정의 불투명도 계산
  ! 지역 조건, 조성 등에 의존할 수 있음
end subroutine
```

### 특정 물리를 위한 불투명도 수정

1. **불투명도 범프** (예: 철 불투명도 피크):
   ```fortran
   case('bump')
     kappa = kappa0*(1 + n*exp(-log(ρ/ρ0)²/σ²))
   ```

2. **속도 의존 불투명도** (도플러 효과):
   - 빠른 흐름에 중요
   - 유효 광학 깊이 수정

3. **비-LTE 보정**:
   - 열평형으로부터의 이탈
   - 수정된 준위 밀도

## 모범 사례 및 권장 사항

### 불투명도 법칙 선택

1. **항성 응용의 경우**:
   - 정확성을 위해 OPAL 테이블 사용
   - 올바른 조성 (X, Y, Z) 확인
   - 온도/밀도 범위 유효성 확인

2. **고온 플라즈마의 경우** (T > 10⁷ K):
   - 톰슨 산란이 종종 충분
   - 필요한 경우 자유-자유를 위해 크라머스 추가

3. **빠른 탐색의 경우**:
   - 상수 또는 멱법칙 불투명도로 시작
   - 매개변수가 설정되면 테이블로 개선

### 수치적 고려사항

1. **불투명도 하한/상한**:
   ```fortran
   kappa = max(kappa_min, min(kappa_max, kappa))
   ```
   극한 조건에서 수치 문제 방지

2. **부드러운 전환**:
   - 안정성을 위해 플럭스 제한자 필터링 사용
   - 확산 계수 필터링 고려

3. **시간 단계 제약**:
   ```fortran
   dt_rad = min(dt_rad_force, dt_rad_diffusion, dt_energy_exchange)
   ```

### 검증 및 테스트

1. **마샤크 파 테스트**: 차가운 물질로의 복사 확산
2. **평형 확산**: 평형에서 E_rad → aT⁴ 확인
3. **광학 깊이 영역**: τ << 1 및 τ >> 1 한계 테스트
4. **에너지 보존**: 총 (가스 + 복사) 에너지 모니터링

## 일반적인 문제 해결

### 문제 1: OPAL 테이블 경계 초과
**증상**: 외삽 경고 또는 오류
**해결책**: 
- 온도/밀도 범위 확인
- 극단에서 해석적 불투명도로 전환 고려
- 톰슨/크라머스로의 부드러운 전환 구현

### 문제 2: 에너지 상호작용 수렴
**증상**: 근 찾기가 수렴하지 못함
**해결책**:
- `fld_bisect_tol` 증가
- 이분법으로 전환 (가장 견고함)
- 음의 온도/에너지 확인

### 문제 3: 복사 확산 불안정성
**증상**: 복사 에너지의 진동
**해결책**:
- 플럭스 제한자 필터링 활성화
- 복사에 대한 CFL 수 감소
- 불투명도 기울기 확인

### 문제 4: 잘못된 불투명도 단위
**증상**: 비물리적인 복사 전달 속도
**해결책**:
- unit_opacity 일관성 확인
- CGS 대 코드 단위 변환 확인
- 테이블 단위 확인 (OPAL의 경우 cm²/g)

## 결론

AMRVAC의 복사 모듈은 다양한 천체물리학적 영역에서 복사 유체역학을 모델링하기 위한 포괄적인 프레임워크를 제공합니다. OPAL 불투명도 테이블의 통합은 항성 내부에서 실험실 플라즈마까지 복사-물질 상호작용의 정확한 처리를 가능하게 합니다. 주요 강점은 다음과 같습니다:

1. **유연성**: 다양한 물리를 위한 여러 불투명도 법칙
2. **정확성**: 최첨단 OPAL 테이블
3. **효율성**: 최적화된 알고리즘과 솔버
4. **확장성**: 사용자 정의 불투명도 및 수정

불투명도, 복사 전달 및 수치 방법 간의 상호작용을 이해하는 것은 AMRVAC를 사용한 성공적인 복사-유체역학 시뮬레이션에 중요합니다.

## 참고문헌 및 추가 자료

### 주요 참고문헌

1. **AMRVAC FLD 구현**:
   - Moens, N., et al. (2022), "Radiation-hydrodynamics with MPI-AMRVAC: Flux-limited diffusion", A&A, 657, A81

2. **OPAL 불투명도 프로젝트**:
   - Iglesias, C. A. & Rogers, F. J. (1996), "Updated OPAL Opacities", ApJ, 464, 943
   - [OPAL 웹사이트](https://opalopacity.llnl.gov/)

3. **CAK 이론**:
   - Castor, J. I., Abbott, D. C., & Klein, R. I. (1975), "Radiation-driven winds in Of stars", ApJ, 195, 157
   - Gayley, K. G. (1995), "An Improved Line-Force Formula for Stellar Winds", ApJ, 454, 410

4. **FLD 방법**:
   - Turner, N. J. & Stone, J. M. (2001), "A Module for Radiation Hydrodynamic Calculations with ZEUS-2D", ApJS, 135, 95
   - Levermore, C. D. & Pomraning, G. C. (1981), "A flux-limited diffusion theory", ApJ, 248, 321

### 코드 문서

- [AMRVAC 문서](http://amrvac.org/)
- [AMRVAC GitHub 저장소](https://github.com/amrvac/amrvac)

### 관련 주제

1. **복사 유체역학**:
   - Mihalas, D. & Mihalas, B. W. (1984), "Foundations of Radiation Hydrodynamics"
   - Castor, J. I. (2004), "Radiation Hydrodynamics", Cambridge University Press

2. **항성 대기**:
   - Gray, D. F. (2005), "The Observation and Analysis of Stellar Photospheres"
   - Hubeny, I. & Mihalas, D. (2014), "Theory of Stellar Atmospheres"

3. **수치 방법**:
   - Stone, J. M., Mihalas, D., & Norman, M. L. (1992), "ZEUS-2D: A radiation magnetohydrodynamics code", ApJS, 80, 819

---

*문서 버전: 1.0*  
*최종 업데이트: 2024년 12월*  
*작성자: Claude에 의한 기술 분석*