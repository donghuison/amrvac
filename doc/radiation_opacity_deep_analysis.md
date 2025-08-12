# Deep Analysis of Radiation Module and Opacity in AMRVAC

## Executive Summary

This document provides a comprehensive analysis of the radiation transport implementation and opacity calculations in the MPI-AMRVAC code, with particular emphasis on the OPAL opacity tables and their integration within the Flux-Limited Diffusion (FLD) framework. The analysis is based on a thorough examination of the source code, including `mod_fld.t`, `mod_opal_opacity.t`, `mod_cak_opacity.t`, and related modules.

## Table of Contents

1. [Introduction](#introduction)
2. [Radiation Transport Framework](#radiation-transport-framework)
3. [OPAL Opacity Tables](#opal-opacity-tables)
4. [Opacity Types and Physical Processes](#opacity-types-and-physical-processes)
5. [Implementation Details](#implementation-details)
6. [Radiation-Matter Coupling](#radiation-matter-coupling)
7. [Computational Workflow](#computational-workflow)
8. [Physical Regimes and Applications](#physical-regimes-and-applications)
9. [References and Further Reading](#references-and-further-reading)

## Introduction

AMRVAC (Adaptive Mesh Refinement Versatile Advection Code) implements sophisticated radiation hydrodynamics capabilities through multiple modules designed to handle different regimes of radiation transport. The code addresses three primary radiation transport scenarios:

1. **Optically thick transport** via Flux-Limited Diffusion (FLD)
2. **Line-driven winds** using CAK theory
3. **Optically thin cooling** through tabulated cooling curves

This document focuses primarily on the FLD implementation and its use of opacity tables, particularly the OPAL (Opacity Project At Livermore) tables.

## Radiation Transport Framework

### Flux-Limited Diffusion (FLD) Module

The FLD module (`mod_fld.t`) implements the flux-limited diffusion approximation for radiation transport, suitable for optically thick regimes where the photon mean free path is much smaller than the characteristic system scale. The module was developed based on Turner & Stone (2001) and is described in detail in Moens et al. (2022, A&A 657).

#### Key Components

1. **Radiation Energy Density Evolution**
   ```fortran
   ∂E_rad/∂t + ∇·F_rad = -cκρ(E_rad - aT^4)
   ```
   where:
   - `E_rad`: Radiation energy density
   - `F_rad`: Radiation flux
   - `κ`: Opacity (cm²/g)
   - `ρ`: Mass density
   - `T`: Temperature
   - `a`: Radiation constant
   - `c`: Speed of light

2. **Flux Limiter Function**
   The radiation flux is computed as:
   ```fortran
   F_rad = -cλ(R)/(3κρ) ∇E_rad
   ```
   where λ(R) is the flux limiter that transitions between:
   - Optically thick diffusion: λ = 1/3
   - Free streaming: λ = 1/R
   - R = |∇E_rad|/(κρE_rad) is the dimensionless gradient

### Available Flux Limiters

AMRVAC implements several flux limiter prescriptions:
- **Pomraning**: Default limiter with smooth transition
- **Diffusion**: Pure diffusion limit (λ = 1/3)
- **FreeStream**: Considers the full R-dependent limiter

## OPAL Opacity Tables

### Table Structure and Format

OPAL tables provide Rosseland-mean opacities as a function of temperature and density. The tables are structured as follows:

#### Parameter Space
- **Temperature Range**: log₁₀(T) ∈ [3.75, 8.70]
  - Physical range: ~5,600 K to ~5×10⁸ K
  - 70 temperature points in standard tables

- **Density Parameter**: log₁₀(R) ∈ [-8.0, 1.0]
  - R = ρ/(T/10⁶ K)³ in CGS units
  - 19 density points per temperature

#### Table Format Example
```
TABLE # 8  X=0.0000 Y=0.9800 Z=0.0200  (Wolf-Rayet atmosphere)

                    log R
logT  -8.0  -7.5  -7.0  ...  0.5   1.0
3.75  9.999 9.999 -3.233 ... -1.806 -1.346
3.80  9.999 9.999 -3.298 ... -1.687 -1.120
...
```

where:
- `X`: Hydrogen mass fraction
- `Y`: Helium mass fraction
- `Z`: Metal mass fraction
- Values are log₁₀(κ) in cm²/g
- `9.999`: No data available for these conditions

### Physical Interpretation of R Parameter

The density parameter R = ρ/(T/10⁶)³ captures the combined effect of:
1. **Number density of absorbers**: Proportional to ρ
2. **Ionization state**: Strong temperature dependence
3. **Radiation pressure scaling**: T³ relation in stellar interiors

This parameterization is particularly useful because:
- It reduces the 2D (ρ,T) space to a more compact representation
- It naturally captures the physics of stellar interiors where radiation pressure dominates
- It provides smooth interpolation across wide parameter ranges

### Implementation in AMRVAC

#### Module Structure (`mod_opal_opacity.t`)

```fortran
module mod_opal_opacity
  ! Table dimensions
  integer, parameter :: iRmin = 2, iRmax = 20
  integer, parameter :: iTmin = 7, iTmax = 76
  
  ! Global storage
  double precision, public :: Kappa_vals(iTmin:iTmax,iRmin:iRmax)
  double precision, public :: logR_list(iRmin:iRmax)
  double precision, public :: logT_list(iTmin:iTmax)
  
  public :: init_opal_table
  public :: set_opal_opacity
end module
```

#### Opacity Lookup Process

1. **Initialization** (performed once):
   ```fortran
   call init_opal_table(tablename)
   ! Reads table from AMRVAC_DIR/src/rhd/OPAL_tables/
   ```

2. **Runtime Lookup**:
   ```fortran
   ! Convert to OPAL convention
   logR_in = log10(rho/(temp*1d-6)**3)
   logT_in = log10(temp)
   
   ! Bilinear interpolation
   call get_kappa(Kappa_vals, logR_list, logT_list, 
                  logR_in, logT_in, logKappa_out)
   
   ! Convert to physical units
   kappa = 10**logKappa_out
   ```

3. **Interpolation Algorithm**:
   - Identifies four surrounding grid points in (logR, logT) space
   - Performs bilinear interpolation in log-log space
   - Handles extrapolation for out-of-bounds values
   - Special handling for missing data (9.999 values)

## Opacity Types and Physical Processes

### Available Opacity Laws in FLD

AMRVAC's FLD module supports multiple opacity prescriptions through the `fld_opacity_law` parameter:

#### 1. Constant Opacity (`'const'`)
```fortran
fld_kappa = fld_kappa0/unit_opacity
```
- Simple, computationally efficient
- Useful for testing and idealized problems

#### 2. Thomson Scattering (`'thomson'`)
```fortran
sigma_thomson = 6.6524585d-25  ! cm²
fld_kappa = sigma_thomson/m_p * (1+2*Y)/(1+4*Y)
```
- Pure electron scattering
- Frequency-independent
- Dominant in hot, fully ionized plasmas

#### 3. Kramers Opacity (`'kramers'`)
```fortran
fld_kappa ∝ ρ * T^(-3.5)
```
- Approximates free-free (Bremsstrahlung) absorption
- Valid for hot, ionized plasmas
- Simple analytic form

#### 4. OPAL Tables (`'opal'`)
```fortran
case('opal')
  call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)
  rho0 = w(ix^D,iw_rho)*unit_density
  Temp0 = Temp(ix^D)*unit_temperature
  call set_opal_opacity(rho0,Temp0,kappa)
  fld_kappa(ix^D) = kappa/unit_opacity
```
- Most sophisticated option
- Includes all opacity sources
- Temperature and density dependent

#### 5. Non-Isothermal (`'non_iso'`)
```fortran
fld_kappa ∝ (ρ/ρ₀) * (T/T₀)^n
```
- Power-law temperature dependence
- Adjustable exponent n (typically -3.5)

#### 6. FastWind (`'fastwind'`)
- Specific to stellar wind applications
- Includes temperature-dependent enhancements

#### 7. Special/User-Defined (`'special'`)
```fortran
call usr_special_opacity(ixI^L, ixO^L, w, x, fld_kappa)
```
- Allows custom opacity implementations

### CAK Line Opacities

The CAK module (`mod_cak_opacity.t`) provides specialized opacities for line-driven stellar winds:

```fortran
module mod_cak_opacity
  ! CAK parameters (Gayley 1995 notation)
  double precision :: alpha_vals    ! Power-law exponent
  double precision :: Qbar_vals     ! Mean line strength
  double precision :: Q0_vals       ! Continuum normalization
  double precision :: kappae_vals   ! Electron scattering
```

These parameters characterize:
- **Line distribution**: α describes the distribution of line strengths
- **Line force multiplier**: M(t) = k*t^(-α) where t is the optical depth parameter
- **Continuum vs line opacity**: Ratio determines wind acceleration efficiency

### Physical Processes Contributing to Opacity

#### Temperature Regimes

1. **Low Temperature (T < 10⁴ K)**
   - Molecular absorption (H₂, H₂O, CO, etc.)
   - Atomic line absorption
   - Dust grains (if present)

2. **Intermediate Temperature (10⁴ K < T < 10⁶ K)**
   - Bound-free transitions (photoionization)
   - Bound-bound transitions (line absorption)
   - Partially ionized metals (especially iron-peak elements)
   - "Z-bump" around 10⁵ K from iron opacity

3. **High Temperature (T > 10⁶ K)**
   - Free-free absorption (Bremsstrahlung)
   - Electron scattering (Thomson/Compton)
   - Fully ionized plasma

#### Composition Dependence

Opacity strongly depends on chemical composition:
- **Hydrogen (X)**: Provides electrons for scattering
- **Helium (Y)**: Less opacity than hydrogen per unit mass
- **Metals (Z)**: Dramatically increase opacity, especially at intermediate temperatures

## Radiation-Matter Coupling

### Energy Exchange Mechanism

The energy interaction between radiation and matter is governed by:

```fortran
dE_gas/dt = cκρ(E_rad - aT⁴)
```

This leads to a coupled system that AMRVAC solves using various methods:

#### Interaction Methods

1. **Instant Equilibrium**
   - Assumes instant thermal equilibrium
   - E_gas + E_rad = constant
   - Solves for equilibrium temperature

2. **Implicit Time Integration**
   - Solves a 4th-order polynomial:
   ```fortran
   e_gas⁴ + c₁*e_gas - c₀ = 0
   ```
   where:
   - `c₁ = (1 + a₂)/a₁`
   - `c₀ = ((1 + a₂)*e_gas + a₂*E_rad)/a₁`
   - `a₁ = 4κσ_B(γ-1)⁴/ρ³ Δt`
   - `a₂ = cκρΔt`

#### Root-Finding Methods

AMRVAC implements three methods to solve the energy polynomial:

1. **Bisection Method**
   ```fortran
   subroutine Bisection_method(e_gas, E_rad, c0, c1)
     ! Robust but slower convergence
     ! Guaranteed to find root in bracket
   ```

2. **Newton-Raphson Method**
   ```fortran
   subroutine Newton_method(e_gas, E_rad, c0, c1)
     ! Faster convergence (quadratic)
     ! May fail for poor initial guess
   ```

3. **Halley's Method**
   ```fortran
   subroutine Halley_method(e_gas, E_rad, c0, c1)
     ! Cubic convergence
     ! More stable than Newton for this problem
   ```

### Radiation Force

The radiation force on matter is computed as:

```fortran
F_rad = (κρ/c) * F_rad
```

This force:
- Accelerates stellar winds
- Drives radiation pressure instabilities
- Supports stars against gravity in radiation-dominated regimes

### Diffusion Coefficient

The radiation diffusion coefficient in the FLD approximation:

```fortran
D = c*λ/(3κρ)
```

This coefficient:
- Controls the rate of radiation energy transport
- Depends on both opacity and flux limiter
- Transitions from D = c/(3κρ) in optically thick regions to D → c in optically thin regions

## Computational Workflow

### Initialization Phase

1. **Read Parameter File**
   ```fortran
   namelist /fld_list/ fld_kappa0, fld_opacity_law, 
                       fld_opal_table, fld_fluxlimiter, ...
   ```

2. **Initialize Opacity Tables** (if using OPAL)
   ```fortran
   if (fld_opacity_law .eq. 'opal') then
     call init_opal_table(fld_opal_table)
   endif
   ```

3. **Set Up Multigrid Solver** (for implicit diffusion)
   ```fortran
   if (fld_diff_scheme .eq. 'mg') then
     use_multigrid = .true.
     phys_implicit_update => Diffuse_E_rad_mg
   endif
   ```

### Runtime Execution

#### Per Timestep Operations

1. **Compute Opacity Field**
   ```fortran
   call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
   ```

2. **Calculate Flux Limiter**
   ```fortran
   call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, R)
   ```

3. **Update Radiation Energy** (implicit diffusion)
   ```fortran
   call Diffuse_E_rad_mg(dt, ixI^L, ixO^L, w, wCT, x)
   ```

4. **Energy Interaction** (if not operator split)
   ```fortran
   call Energy_interaction(w, wCT, x, ixI^L, ixO^L)
   ```

5. **Apply Radiation Force** (if not operator split)
   ```fortran
   call get_fld_rad_force(qdt, ixI^L, ixO^L, wCT, w, x)
   ```

### Performance Optimizations

1. **Table Caching**: OPAL tables loaded once at initialization
2. **Flux Limiter Filtering**: Optional running average to smooth gradients
3. **Diffusion Coefficient Filtering**: Reduces numerical noise
4. **Multigrid Acceleration**: Efficient implicit diffusion solver
5. **Operator Splitting**: Separate treatment of stiff terms

## Physical Regimes and Applications

### Stellar Atmospheres

OPAL opacities are particularly suited for:
- **Stellar interiors**: Where radiation pressure is significant
- **Stellar envelopes**: Transition from convective to radiative transport
- **Wolf-Rayet winds**: High temperature, helium-rich environments

Example configuration for Wolf-Rayet star:
```fortran
&fld_list
  fld_opacity_law = 'opal'
  fld_opal_table = 'Y09800'  ! Y=0.98, Z=0.02
  fld_fluxlimiter = 'Pomraning'
  fld_diff_scheme = 'mg'
/
```

### Accretion Disks

FLD with appropriate opacities can model:
- **Disk vertical structure**: Radiation pressure support
- **Thermal instabilities**: S-curve behavior
- **Radiation-driven outflows**: From inner disk regions

### Supernova Remnants

For high-temperature shocks:
```fortran
&fld_list
  fld_opacity_law = 'thomson'  ! Electron scattering dominates
  fld_fluxlimiter = 'Pomraning'
  fld_interaction_method = 'Halley'
/
```

### Laboratory Astrophysics

Simplified opacities for controlled experiments:
```fortran
&fld_list
  fld_opacity_law = 'const'
  fld_kappa0 = 0.34d0  ! cm²/g
/
```

## Advanced Features

### Line-Force Opacities

For CAK-type line-driven winds, AMRVAC can use spatially varying line-force parameters:

```fortran
if (lineforce_opacities) then
  ! Allocate arrays for directional opacities
  allocate(i_opf(ndim))
  do idir = 1,ndim
    i_opf(idir) = var_set_extravar('k'//ind_1, 'k'//ind_1)
  enddo
endif
```

### User-Defined Opacities

Custom opacity implementations via user module:

```fortran
subroutine usr_special_opacity(ixI^L, ixO^L, w, x, kappa)
  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
  double precision, intent(out) :: kappa(ixO^S)
  
  ! Custom opacity calculation
  ! Can depend on local conditions, composition, etc.
end subroutine
```

### Opacity Modifications for Specific Physics

1. **Opacity Bumps** (e.g., iron opacity peak):
   ```fortran
   case('bump')
     kappa = kappa0*(1 + n*exp(-log(ρ/ρ0)²/σ²))
   ```

2. **Velocity-Dependent Opacities** (Doppler effects):
   - Important for fast flows
   - Modifies effective optical depth

3. **Non-LTE Corrections**:
   - Departure from thermal equilibrium
   - Modified level populations

## Best Practices and Recommendations

### Choosing Opacity Laws

1. **For Stellar Applications**:
   - Use OPAL tables for accuracy
   - Ensure correct composition (X, Y, Z)
   - Check temperature/density range validity

2. **For High-Temperature Plasmas** (T > 10⁷ K):
   - Thomson scattering often sufficient
   - Add Kramers for free-free if needed

3. **For Rapid Exploration**:
   - Start with constant or power-law opacity
   - Refine with tables once parameters established

### Numerical Considerations

1. **Opacity Floor/Ceiling**:
   ```fortran
   kappa = max(kappa_min, min(kappa_max, kappa))
   ```
   Prevents numerical issues in extreme conditions

2. **Smooth Transitions**:
   - Use flux limiter filtering for stability
   - Consider diffusion coefficient filtering

3. **Time Step Constraints**:
   ```fortran
   dt_rad = min(dt_rad_force, dt_rad_diffusion, dt_energy_exchange)
   ```

### Validation and Testing

1. **Marshak Wave Test**: Radiation diffusion into cold material
2. **Equilibrium Diffusion**: Verify E_rad → aT⁴ in equilibrium
3. **Optical Depth Regimes**: Test τ << 1 and τ >> 1 limits
4. **Energy Conservation**: Monitor total (gas + radiation) energy

## Troubleshooting Common Issues

### Issue 1: OPAL Table Out of Bounds
**Symptom**: Extrapolation warnings or errors
**Solution**: 
- Check temperature/density ranges
- Consider switching to analytical opacity at extremes
- Implement smooth transition to Thomson/Kramers

### Issue 2: Energy Interaction Convergence
**Symptom**: Root-finding fails to converge
**Solution**:
- Increase `fld_bisect_tol`
- Switch to bisection method (most robust)
- Check for negative temperatures/energies

### Issue 3: Radiation Diffusion Instability
**Symptom**: Oscillations in radiation energy
**Solution**:
- Enable flux limiter filtering
- Reduce CFL number for radiation
- Check opacity gradients

### Issue 4: Incorrect Opacity Units
**Symptom**: Unphysical radiation transport rates
**Solution**:
- Verify unit_opacity consistency
- Check CGS vs code units conversion
- Confirm table units (cm²/g for OPAL)

## Conclusion

AMRVAC's radiation module provides a comprehensive framework for modeling radiation hydrodynamics across diverse astrophysical regimes. The integration of OPAL opacity tables enables accurate treatment of radiation-matter interactions from stellar interiors to laboratory plasmas. Key strengths include:

1. **Flexibility**: Multiple opacity laws for different physics
2. **Accuracy**: State-of-the-art OPAL tables
3. **Efficiency**: Optimized algorithms and solvers
4. **Extensibility**: User-defined opacities and modifications

Understanding the interplay between opacity, radiation transport, and numerical methods is crucial for successful radiation-hydrodynamic simulations with AMRVAC.

## References and Further Reading

### Primary References

1. **AMRVAC FLD Implementation**:
   - Moens, N., et al. (2022), "Radiation-hydrodynamics with MPI-AMRVAC: Flux-limited diffusion", A&A, 657, A81

2. **OPAL Opacity Project**:
   - Iglesias, C. A. & Rogers, F. J. (1996), "Updated OPAL Opacities", ApJ, 464, 943
   - [OPAL Website](https://opalopacity.llnl.gov/)

3. **CAK Theory**:
   - Castor, J. I., Abbott, D. C., & Klein, R. I. (1975), "Radiation-driven winds in Of stars", ApJ, 195, 157
   - Gayley, K. G. (1995), "An Improved Line-Force Formula for Stellar Winds", ApJ, 454, 410

4. **FLD Method**:
   - Turner, N. J. & Stone, J. M. (2001), "A Module for Radiation Hydrodynamic Calculations with ZEUS-2D", ApJS, 135, 95
   - Levermore, C. D. & Pomraning, G. C. (1981), "A flux-limited diffusion theory", ApJ, 248, 321

### Code Documentation

- [AMRVAC Documentation](http://amrvac.org/)
- [AMRVAC GitHub Repository](https://github.com/amrvac/amrvac)

### Related Topics

1. **Radiation Hydrodynamics**:
   - Mihalas, D. & Mihalas, B. W. (1984), "Foundations of Radiation Hydrodynamics"
   - Castor, J. I. (2004), "Radiation Hydrodynamics", Cambridge University Press

2. **Stellar Atmospheres**:
   - Gray, D. F. (2005), "The Observation and Analysis of Stellar Photospheres"
   - Hubeny, I. & Mihalas, D. (2014), "Theory of Stellar Atmospheres"

3. **Numerical Methods**:
   - Stone, J. M., Mihalas, D., & Norman, M. L. (1992), "ZEUS-2D: A radiation magnetohydrodynamics code", ApJS, 80, 819

---

*Document Version: 1.0*  
*Last Updated: December 2024*  
*Author: Technical Analysis by Claude*