# CloudMicrophysics.jl
A package containing a library of cloud microphysics and aerosol parameterizations.

|||
|---------------------:|:----------------------------------------------|
| **Docs Build**       | [![docs build][docs-bld-img]][docs-bld-url]   |
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url]          |
| **GHA CI**           | [![gha ci][gha-ci-img]][gha-ci-url]           |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |

[docs-bld-img]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/docs.yml/badge.svg
[docs-bld-url]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/docs.yml

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/CloudMicrophysics.jl/dev/

[gha-ci-img]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/CloudMicrophysics.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/CloudMicrophysics.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/CloudMicrophysics.jl

## Installation and running instructions

CloudMicrophysics.jl is a Julia registered package.
It depends on a couple of standard Julia packages as well as
  the [Thermodynamics.jl](https://github.com/CliMA/Thermodynamics.jl) and
  [CLIMAParameters.jl](https://github.com/CliMA/CLIMAParameters.jl)
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

The CloudMicrophysics.jl package is being actively developed
  and welcomes contributions and feedback.
There is a variety of projects big and small that are available to take up as
  fun research projects for interested students and other contributors.
Below is a list of possible examples,
  but other suggestions and ideas are always welcome!

- The CloudMicrophysics.jl package should be tested against a high-resolution model.
  We have chosen [PySDM](https://github.com/atmos-cloud-sim-uj/PySDM)
  as our high-resolution benchmark.
  PySDM is a package for simulating the dynamics of population of particles
  and is based on the [Super-Droplet algorithm](https://doi.org/10.1002/qj.441).
  Possible tasks in this project would include testing the aerosol activation parameterization
  against PySDM in an adiabatic parcel setup, or testing the 1-moment
  microphysics parameterization against PySDM in an already implemented
  [1-dimensional](https://github.com/CliMA/Kinematic1D.jl) or
  2-dimensional prescribed flow setup.
  This could be extended further into a calibration exercise using the
  [EnsembleKalmanProcesses.jl](https://github.com/CliMA/EnsembleKalmanProcesses.jl) package.
  An example pipeline can be seen in the
  [EKP.jl docs](https://clima.github.io/EnsembleKalmanProcesses.jl/dev/examples/Cloudy_example/)
  where [Cloudy.jl](https://github.com/CliMA/Cloudy.jl) parameters are calibrated.

- Adding an aerosol model and coupling it with the aerosol activation
  parameterization.
  [MAM4](https://doi.org/10.5194/gmd-9-505-2016) could be the aerosol model to implement.
  This is a big project and an opportunity for a more long term contribution.

- Adding the [P3 scheme](https://doi.org/10.1175/JAS-D-14-0065.1) for ice phase microphysics.

- Adding more accurate fall speed parameterization based on
  [Chen et. al. 2022](https://doi.org/10.1016/j.atmosres.2022.106171).

- Adding precipitation susceptibility tests based on
  [Glassmeier and Lohmann 2016](https://doi.org/10.1175/JAS-D-16-0008.1).
