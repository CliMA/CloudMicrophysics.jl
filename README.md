# CloudMicrophysics.jl
A package containing a library of cloud microphysics and aerosol parameterizations.
See [our documentation](https://clima.github.io/CloudMicrophysics.jl/dev/) for the list of available schemes.

|||
|---------------------:|:----------------------------------------------|
| **Docs Build**       | [![docs build][docs-bld-img]][docs-bld-url]   |
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url]          |
| **GHA CI**           | [![gha ci][gha-ci-img]][gha-ci-url]           |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |
| **Downloads**        | [![downloads][downloads-img]][downloads-url]  |

[docs-bld-img]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/docs.yml/badge.svg
[docs-bld-url]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/docs.yml

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/CloudMicrophysics.jl/dev/

[gha-ci-img]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/CloudMicrophysics.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/CloudMicrophysics.jl

[downloads-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FCloudMicrophysics&query=total_requests&suffix=%2Ftotal&label=Downloads
[downloads-url]: http://juliapkgstats.com/pkg/CloudMicrophysics

## Installation and running instructions

CloudMicrophysics.jl is a Julia registered package.
It depends on a couple of standard Julia packages as well as
  the [Thermodynamics.jl](https://github.com/CliMA/Thermodynamics.jl) and
  [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl)
  (two Julia packages developed at [CliMA](https://github.com/CliMA)).
See the [Project.toml](https://github.com/CliMA/CloudMicrophysics.jl/blob/main/Project.toml)
  for a full list of CloudMicrophysics.jl's dependencies.

When using the CloudMicrophysics.jl inside your own project,
  the easiest way to obtain the latest stable version of the package
  and it's dependencies is to use the Julia built-in package manager
  (accessed by pressing `]` in the Julia REPL):

```julia
julia>]
pkg> add CloudMicrophysics
pkg> instantiate
```

CloudMicrophysics.jl can be updated to the latest tagged release
  from the package manager.
The package is still under development and changes to API are very possible!

```julia
pkg> update CloudMicrophysics
```

When contributing to the CloudMicrophysics.jl development,
  the easiest way is to clone the [repository](https://github.com/CliMA/CloudMicrophysics.jl)
  and then run it using its project environment.
For example, to get all the needed dependencies and then run all the tests
  you could try:

```julia
julia --project=test

julia>]

pkg> dev .

pkg> instantiate

julia> include("test/runtests.jl")
```

See the [Pkg docs](https://docs.julialang.org/en/v1/stdlib/Pkg/)
  for an overview of basic package manager features.

## Contributing

CloudMicrophysics.jl is being actively developed
  and welcomes contributions and feedback.
There is a variety of projects big and small that are available to take up as
  fun research projects for interested students and other contributors.
Below is a list of possible examples,
  but other suggestions and ideas are always welcome!

- CloudMicrophysics.jl should be tested against a high-resolution model.
  We have chosen [PySDM](https://github.com/atmos-cloud-sim-uj/PySDM)
  as our high-resolution benchmark.
  PySDM is a package for simulating the dynamics of population of particles
  and is based on the [Super-Droplet algorithm](https://doi.org/10.1002/qj.441).
  Possible tasks in this project would include testing the aerosol activation parameterization
  against PySDM in adiabatic parcel setup, or testing the 1-moment
  microphysics parameterization against PySDM in an
  [1-dimensional](https://github.com/CliMA/Kinematic1D.jl) or
  2-dimensional prescribed flow setups.
  This could be extended further into a calibration exercise using the
  [EnsembleKalmanProcesses.jl](https://github.com/CliMA/EnsembleKalmanProcesses.jl) package.
  An example pipeline can be seen in the
  [EKP.jl docs](https://clima.github.io/EnsembleKalmanProcesses.jl/dev/examples/Cloudy_example/)
  where [Cloudy.jl](https://github.com/CliMA/Cloudy.jl) parameters are calibrated.

- The CloudMicrophysics.jl package should be tested against observations.
  We are focusing on the ice-free precipitation first. The tests include
  comparisons against [CFODDs](https://doi.org/10.1175/JAS-D-20-0321.1) and
  against [Stratocumulus LWP(N) patterns](https://doi.org/10.5194/acp-19-10191-2019).

- Adding an aerosol model and coupling it with the aerosol activation
  parameterization.
  [MAM4](https://doi.org/10.5194/gmd-9-505-2016) could be the aerosol model to implement,
  but we are also searching for some simpler alternatives first.
  This is a big project and an opportunity for a more long term contribution.
