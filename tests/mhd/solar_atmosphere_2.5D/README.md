# Solar Atmosphere 2.5D MHD Test Case

## 📋 Overview

This test case simulates a 2.5D magnetohydrodynamic (MHD) model of the solar atmosphere with bipolar magnetic field configuration. It serves as a fundamental benchmark for solar physics applications, particularly for studying magnetic reconnection, wave propagation, and energy transport in the solar atmosphere.

### Key Features
- **2.5D MHD simulation** (2D spatial dimensions with 3D vector fields)
- **Multi-layer solar atmosphere** (photosphere, chromosphere, corona)
- **Bipolar arcade magnetic field** configuration
- **Advanced physics modules** (thermal conduction, radiative cooling, gravity)
- **Multiple numerical schemes** validation framework
- **Adaptive Mesh Refinement** (AMR) capability

### Scientific Applications
- Solar flare initiation and evolution
- Coronal heating mechanisms
- Magnetic reconnection studies
- Wave propagation in stratified atmospheres
- Prominence formation and dynamics

## 🌟 Physical Model

### Solar Atmosphere Stratification

The model implements a realistic vertical stratification of the solar atmosphere with three distinct layers:

1. **Photosphere** (bottom layer)
   - Temperature: ~8,000 K
   - High density region
   - Partially ionized plasma

2. **Chromosphere** (middle layer)
   - Temperature: 8,000 - 20,000 K
   - Transition region
   - Strong temperature gradients

3. **Corona** (top layer)
   - Temperature: ~1.5 MK
   - Low density, fully ionized plasma
   - Magnetically dominated region

### Magnetic Field Configuration

The magnetic field follows a bipolar arcade structure:
```
B_x = -B₀ cos(kx·x) exp(-ly·y) cos(θ)
B_y = +B₀ sin(kx·x) exp(-ly·y)
B_z = -B₀ cos(kx·x) exp(-ly·y) sin(θ)
```

Where:
- `B₀ = 20 G` (field strength at photosphere)
- `θ = 60°` (angle to the xy-plane)
- `kx = π/Lx` (horizontal wavenumber)
- `ly = kx·cos(θ)` (vertical decay rate)

### Physical Processes

#### 1. Gravity
Solar gravity with height-dependent acceleration:
```
g(y) = g_sun × (R_sun/(R_sun + y))²
```
- `g_sun = -274 m/s²` (solar surface gravity)
- `R_sun = 696.1 Mm` (solar radius)

#### 2. Thermal Conduction
Anisotropic Spitzer thermal conduction:
- Parallel to magnetic field lines
- Temperature-dependent coefficient: κ ∝ T^(5/2)
- Saturation flux limiting

#### 3. Radiative Cooling
- **Cooling curve**: JCcorona (optically thin losses)
- **Temperature range**: 10^4 - 10^7 K
- **Cooling method**: Exact integration
- **Temperature floor**: 8,000 K

#### 4. Heating
Two heating mechanisms:
- **Background heating**: Exponential decay with height
  - `Q_b = Q₀ exp(-y/H)` where H = 5 Mm
- **Localized heating** (optional): Magnetic field-aligned heating

## 💻 Numerical Implementation

### MHD Equations Solved

The code solves the ideal MHD equations with additional source terms:

```
∂ρ/∂t + ∇·(ρv) = 0
∂(ρv)/∂t + ∇·(ρvv + pI - BB) = ρg
∂e/∂t + ∇·((e+p)v - B(B·v)) = ρg·v + Q_heat - Q_cool + ∇·(κ∇T)
∂B/∂t + ∇×(B×v) = 0
```

### Boundary Conditions

#### Bottom Boundary (y = 0)
- **Type**: Special boundary condition
- **Velocity**: Fixed (v = 0) or driven horizontal flow
- **Density/Pressure**: Hydrostatic stratification
- **Magnetic field**: Fixed or extrapolated

#### Top Boundary (y = 6 Mm)
- **Type**: Special boundary condition
- **Velocity**: Zero-gradient extrapolation
- **Density/Pressure**: Hydrostatic stratification
- **Magnetic field**: Zero-gradient extrapolation

#### Lateral Boundaries (x-direction)
- **Type**: Symmetric/Antisymmetric for different variables
- Preserves arcade structure

### AMR Strategy

- **Maximum refinement level**: 4
- **Refinement criteria**: Error estimator (Löhner)
- **Variables weighted**: ρ (40%), B components (20% each)
- **Special refinement**: Bottom layer fixed at maximum level
- **Block size**: 12×12 cells

## 📁 Test Scenarios

### 1. Relaxation Phase (`relax.par`)
- **Purpose**: Relax initial magnetic field to numerical equilibrium
- **Duration**: 80 time units
- **Key features**: 
  - No driving flows
  - Full physics enabled
  - Constrained transport (CT) for ∇·B = 0

### 2. Convection Simulation (`conv.par`)
- **Purpose**: Study convection-driven dynamics
- **Duration**: 40 time units
- **Key features**:
  - Driven horizontal flows at photosphere
  - Restart from relaxed state
  - UCT-HLL divergence cleaning

### 3. Standard Test (`solar_atm_25D.par`)
- **Purpose**: Regression testing framework
- **Duration**: 10 iterations
- **Key features**:
  - Multiple numerical schemes tested
  - Log file output for validation
  - GLM divergence cleaning

### 4. Full Simulation (`amrvac.par`)
- **Purpose**: Production runs
- **Duration**: 60 time units
- **Key features**:
  - Full physics and AMR
  - VTU output for visualization
  - Auxiliary diagnostic variables

## 🔢 Numerical Schemes Tested

The test framework validates 11 different numerical configurations:

| Scheme | Time Stepper | Flux Scheme | Limiter | Special Features |
|--------|-------------|-------------|---------|------------------|
| `3step_hll_cada_ct` | Three-step RK | HLL | CADA3 | Constrained Transport |
| `2step_tvdlf_mm` | Two-step | TVDLF | MinMod | Basic scheme |
| `3step_hll_cada` | Three-step RK | HLL | CADA3 | Standard |
| `3step_hlld_cada` | Three-step RK | HLLD | CADA3 | Enhanced MHD |
| `4step_hll_mc` | Four-step RK | HLL | MC | High-order time |
| `4step_hllc_ko` | Four-step RK | HLLC | Koren | Contact resolution |
| `rk4_tvdlf_cada` | RK4 | TVDLF | CADA3 | High accuracy |
| `ssprk54_hlld_mp5` | SSP-RK54 | HLLD | MP5 | Very high-order |
| `B0split` | Three-step | HLL | CADA3 | Background field splitting |
| `mhd_internal_e` | Three-step | HLL | CADA3 | Internal energy equation |
| `mhd_hydrodynamic_e` | Three-step | HLL | CADA3 | Hydrodynamic energy |

## 🚀 Usage Instructions

### Prerequisites
```bash
# Set AMRVAC environment
export AMRVAC_DIR=/path/to/amrvac

# Required modules (system-dependent)
module load mpi fortran
```

### Compilation
```bash
# Configure for 2D
$AMRVAC_DIR/setup.pl -d=2

# Compile with optimization
make -j 4

# Or compile with debugging
make ARCH=debug
```

### Running Simulations

#### Basic Run
```bash
# Run relaxation phase
mpirun -np 4 ./amrvac -i relax.par

# Continue with convection
mpirun -np 4 ./amrvac -i conv.par -restart 40

# Run standard test
mpirun -np 4 ./amrvac -i solar_atm_25D.par
```

#### Testing Framework
```bash
# Run all tests
make -f test.make

# Run specific scheme test
make -f test.make solar_atm_25D_3step_hll_cada.log

# Compare with reference
diff solar_atm_25D_3step_hll_cada.log correct_output/
```

### Output Conversion
```bash
# Convert to VTU for ParaView
aiconvert

# Convert specific snapshots
aiconvert 10 20

# Convert with different parameter file
aiconvert conv.par 0 40
```

## ⚙️ Key Parameters

### Physical Parameters
| Parameter | Value | Description | Unit |
|-----------|-------|-------------|------|
| `unit_length` | 10^9 | Length normalization | cm |
| `unit_temperature` | 10^6 | Temperature normalization | K |
| `unit_numberdensity` | 10^9 | Number density normalization | cm^-3 |
| `Busr` | 20.0 | Magnetic field strength | Gauss |
| `usr_grav` | -274.0 | Solar surface gravity | m/s^2 |

### Numerical Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `courantpar` | 0.8 | CFL number |
| `dtmin` | 1e-7 | Minimum timestep |
| `small_pressure` | 1e-8 | Pressure floor |
| `small_density` | 1e-14 | Density floor |
| `refine_threshold` | 0.2 | AMR refinement threshold |
| `derefine_ratio` | 0.15 | AMR derefinement ratio |

### Domain Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `xprobmin1/max1` | -3.0/3.0 | X-domain (Mm) |
| `xprobmin2/max2` | 0.0/6.0 | Y-domain (Mm) |
| `domain_nx1` | 96 | Base resolution in X |
| `domain_nx2` | 96 | Base resolution in Y |
| `block_nx1/2` | 12 | Block size |

## 📊 Diagnostic Variables

The simulation outputs additional diagnostic variables:

1. **Te**: Temperature (K)
2. **Alfv**: Alfvén speed (km/s)
3. **divB**: Divergence of B (numerical error)
4. **beta**: Plasma beta (2p/B²)
5. **bQ**: Background heating rate
6. **rad**: Radiative cooling rate
7. **j1, j2, j3**: Current density components
8. **Tcutoff**: Temperature cutoff indicator (if enabled)

## 🧪 Validation and Testing

### Test Framework
- **Location**: `correct_output/` directory
- **Format**: Log files with conservation quantities
- **Tolerance**: Relative 1e-5, Absolute 1e-8
- **Variables tested**: Mass, momentum, energy, magnetic flux

### Validation Metrics
Each test monitors:
- Conservation of mass, momentum, energy
- Magnetic field divergence
- Maximum values of key variables
- RMS values for stability assessment

### Running Validation
```bash
# Run test and check results
make -f test.make
grep PASSED test_results.log
```

## 🎨 Visualization and Analysis

### Recommended Tools

#### ParaView
```bash
# Open VTU files
paraview solar_atm_*.vtu

# Suggested visualizations:
# - Temperature (color map)
# - Magnetic field lines (streamlines)
# - Velocity vectors
# - Current density (volume rendering)
```

#### Python Analysis
```python
import amrvac_pytools as apt

# Load data
data = apt.load_datfile('solar_atm_0010.dat')

# Access variables
rho = data.get_var('rho')
temp = data.get_var('Te')
bfield = data.get_var('b')

# Plot
import matplotlib.pyplot as plt
plt.pcolormesh(data.x, data.y, temp)
plt.colorbar(label='Temperature (K)')
```

### Key Diagnostics to Monitor
1. **Energy balance**: Heating vs. cooling rates
2. **Magnetic topology**: Field line connectivity
3. **Flow patterns**: Convection cells, jets
4. **Wave activity**: Alfvén, slow/fast modes
5. **Current sheets**: Reconnection sites

## 📚 References

### Method Papers
- Keppens et al. (2012): "Parallel, grid-adaptive approaches for relativistic hydro and magnetohydrodynamics"
- Porth et al. (2014): "MPI-AMRVAC for Solar and Astrophysics"
- Xia et al. (2018): "MPI-AMRVAC 2.0 for Solar and Astrophysical Applications"

### Solar Physics Applications
- Simulation of coronal rain
- Prominence formation studies
- Flare trigger mechanisms
- Coronal heating investigations

## 🤝 Contributing

To modify or extend this test case:

1. **Physics modifications**: Edit `mod_usr.t`
2. **Parameter variations**: Create new `.par` files
3. **Scheme testing**: Add to `SCHEMES` in `test.make`
4. **Validation**: Update `correct_output/` after verification

## 📧 Support

For questions or issues related to this test case:
- AMRVAC website: http://amrvac.org/
- Documentation: http://amrvac.org/md_doc_test_suites.html
- User guide: http://amrvac.org/md_doc_getting_started.html

## 📝 Notes

### Important Considerations
- Initial relaxation phase is crucial for numerical stability
- Bottom boundary refinement prevents spurious oscillations  
- Thermal conduction requires small timesteps (supertimestepping helps)
- Divergence cleaning method affects long-term evolution

### Common Issues and Solutions
- **Negative pressure**: Adjust `small_pressure` or heating rates
- **Divergence B growth**: Switch divergence cleaning method
- **Slow performance**: Reduce refinement levels or domain size
- **Boundary artifacts**: Check special boundary implementation

---
*Last updated: 2024*
*Test case version: MPI-AMRVAC 3.0*