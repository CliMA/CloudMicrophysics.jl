# CloudMicrophysics Agent Guide

## Ecosystem Guidelines

Please refer to the shared CliMA agent index for ecosystem-wide rules regarding architecture, performance, code quality, infrastructure, and workflows:

- [docs/dev-guides/AGENTS.md](docs/dev-guides/AGENTS.md) — Shared CliMA agent guidelines.

> Shared guides live at `docs/dev-guides/` and are vendored from the canonical source:
> <https://github.com/CliMA/DeveloperGuides>. Edit shared guides there, not here.

## Before You Act: Agent Autonomy

Before making changes that are externally visible or scientifically consequential (`git push`, version bumps, reproducibility-test edits, CI config changes, public API renames), check [docs/dev-guides/workflow/agent_autonomy.md](docs/dev-guides/workflow/agent_autonomy.md). The boundaries listed there require explicit user approval.

## Repo-Specific Guidelines

### Source layout

| Path | Description |
| ------ | ----------- |
| `src/CloudMicrophysics.jl` | Top-level module; all submodules are included here |
| `src/Microphysics0M.jl` | 0-moment microphysics (threshold removal) |
| `src/Microphysics1M.jl` | 1-moment bulk microphysics |
| `src/Microphysics2M.jl` | 2-moment bulk microphysics |
| `src/MicrophysicsNonEq.jl` | Non-equilibrium condensation/evaporation and cloud terminal velocities |
| `src/BulkMicrophysicsTendencies.jl` | Unified tendency API for bulk schemes |
| `src/P3.jl`, `src/P3_*.jl` | P3 ice microphysics scheme (properties, size distribution, processes, velocity) |
| `src/AerosolActivation.jl` | Aerosol activation (CCN) |
| `src/AerosolModel.jl` | Aerosol size distribution types (`Mode_B`, `Mode_κ`, `AerosolDistribution`) |
| `src/IceNucleation.jl` | Heterogeneous and homogeneous ice nucleation |
| `src/Nucleation.jl` | General nucleation routines |
| `src/Common.jl` | Shared helper functions (terminal velocity, ventilation, etc.) |
| `src/CloudDiagnostics.jl` | Radar reflectivity and effective radius diagnostics |
| `src/DistributionTools.jl` | Generalized gamma and exponential distribution tools |
| `src/ThermodynamicsInterface.jl` | Thin wrapper around Thermodynamics.jl |
| `src/Utilities.jl` | Numerical utilities (`clamp_to_nonneg`, `ϵ_numerics`) |
| `src/PrecipitationSusceptibility.jl` | Precipitation susceptibility diagnostics |
| `src/ArtifactCalling.jl` | Artifact (lookup table) loading |
| `src/show.jl` | Pretty-printing for parameter structs |
| `src/parameters/` | ClimaParams parameter definitions |
| `test/` | Unit and integration tests |
| `docs/` | Documenter.jl documentation |
| `ext/` | Weak-dep extensions (`EmulatorModelsExt` for MLJ/DataFrames) |
| `parcel/` | Parcel model examples |
| `box/` | Box model examples |
| `papers/` | Scripts reproducing published figures |
| `p3_sandbox/` | P3 scheme sandbox scripts |

### Testing

- Run the full test suite with `Pkg.test()` (preferred) or `julia --project=test test/runtests.jl`.
- GPU tests live in `test/gpu_tests.jl` and `test/gpu_clima_core_test.jl`.
- Performance and type-stability tests are in `test/performance_tests.jl` and `test/type_stability_tests.jl`.

## Local norms

- All functions must be type-stable and GPU-safe. Use the `FT` (float type) pattern consistently.
- Parameters are accessed through `ClimaParams`; never hard-code physical constants.
- Run `julia -e 'using JuliaFormatter; format(".")'` before committing code.
- Match existing style: explicit names, narrow imports, comments that explain *why*.
- Follow the ecosystem conventions in [docs/dev-guides/architecture/ecosystem_conventions.md](docs/dev-guides/architecture/ecosystem_conventions.md) for new code and refactor toward them when touching existing code.

## Self-correction

- If the source layout table above is discovered to be stale, update it.
- If the user gives a correction about how work should be done in this repo, add it to `Local norms` or another clearly labeled persistent section in this file so future sessions inherit it.
