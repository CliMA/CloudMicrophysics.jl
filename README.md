# CloudMicrophysics.jl

CloudMicrophysics.jl provides a library of cloud microphysics and aerosol parameterizations for the [CliMA Earth System Model](https://clima.caltech.edu). It implements bulk microphysics schemes for cloud formation, precipitation, and aerosol processes, designed for high-performance climate simulations.

|||
|-----------------------------:|:-------------------------------------------------|
| **Documentation**            | [![dev][docs-dev-img]][docs-dev-url]             |
| **Docs Build**               | [![docs build][docs-bld-img]][docs-bld-url]      |
| **GHA CI**                   | [![gha ci][gha-ci-img]][gha-ci-url]              |
| **Code Coverage**            | [![codecov][codecov-img]][codecov-url]           |
| **Downloads**                | [![Downloads][dlt-img]][dlt-url]                 |

[docs-bld-img]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/docs.yml/badge.svg
[docs-bld-url]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/docs.yml

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/CloudMicrophysics.jl/dev/

[gha-ci-img]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/CloudMicrophysics.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/CloudMicrophysics.jl

[dlt-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FCloudMicrophysics&query=total_requests&label=Downloads
[dlt-url]: https://juliapkgstats.com/pkg/CloudMicrophysics

## Quick Start

### Installation

```julia
using Pkg
Pkg.add("CloudMicrophysics")
Pkg.add("ClimaParams")
```

### Basic Usage

```julia
import CloudMicrophysics as CM
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Parameters as CMP

# Create microphysics parameters
rain = CMP.Rain(Float64)
vel = CMP.Blk1MVelType(Float64).rain

# Compute rain terminal velocity
ρ = 1.2      # air density [kg/m³]
q_rai = 1e-3 # rain specific content [kg/kg]
v_term = CM1.terminal_velocity(rain, vel, ρ, q_rai)
```

## Key Features

### 🌧️ **Bulk Microphysics Schemes**

- **0-moment scheme**: Simple precipitation removal
- **1-moment scheme**: Marshall-Palmer distributions for rain and snow
- **2-moment scheme**: [Seifert & Beheng (2006)](https://doi.org/10.1007/s00703-005-0112-4) with mass and number concentration
- **P3 scheme**: Predicted particle properties for ice

### 🧊 **Ice Nucleation**

- **Heterogeneous nucleation**: Deposition, immersion freezing ([ABIFM](https://doi.org/10.5194/acp-12-9817-2012))
- **Homogeneous nucleation**: [Koop et al. (2000)](https://doi.org/10.1038/35020537) parameterization
- **INP distributions**: [Frostenberg et al. (2023)](https://doi.org/10.5194/acp-23-10883-2023)

### 💨 **Aerosol Processes**

- **Aerosol activation**: [Abdul-Razzak & Ghan (2000)](https://doi.org/10.1029/1999JD901161) parameterization
- **Aerosol nucleation**: H₂SO₄ and organic nucleation pathways
- **Aerosol model**: Modal distributions with κ-Köhler theory

### ⚡ **High Performance**

- **Type-stable** and **GPU-compatible** (CUDA.jl, AMDGPU.jl)
- **AD-compatible** (ForwardDiff.jl) for differentiable physics
- Optimized for minimal allocations

## Documentation

- **[Getting Started](https://clima.github.io/CloudMicrophysics.jl/dev/GettingStarted/)** - Installation and first steps
- **[API Reference](https://clima.github.io/CloudMicrophysics.jl/dev/API/)** - Detailed function documentation
- **[Microphysics Schemes](https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics1M/)** - Scheme descriptions

## Integration with Climate Models

CloudMicrophysics.jl is used throughout the [CliMA](https://github.com/CliMA) ecosystem:

- [ClimaAtmos](https://github.com/CliMA/ClimaAtmos.jl) - Atmospheric model
- [KinematicDriver](https://github.com/CliMA/KinematicDriver.jl) - 1D/2D kinematic framework
- [Thermodynamics](https://github.com/CliMA/Thermodynamics.jl) - Moist thermodynamics

## Getting Help

For questions, check the [documentation](https://clima.github.io/CloudMicrophysics.jl/dev/) or open an issue on [GitHub](https://github.com/CliMA/CloudMicrophysics.jl).
