# Microphysics Flexible

The `MicrophysicsFlexible.jl` module relies on the extension defined in `ext/CloudyExt.jl`, based on a flexible N-moment microphysics scheme built in the external package [Cloudy.jl](https://github.com/CliMA/Cloudy.jl). This option currently handles warm-rain processes including coalescence, condensation/evaporation, and sedimentation (terminal velocity). Unlike typical moment-based schemes which distinguish between categories such as rain and cloud, and which determine rates of conversion between categories (the canonical autoconversion, accretion, and self-collection), this option gives the user the flexibility to define as many or as few moments as they please, with these coalescence-based processes being solved directly without relying on conversion rates. Likewise, the rate of condensation/evaporation is defined through the rate of diffusion of water vapor to/from the surface of droplets defined by the subdistributions which underpin the method. The user has not only the flexibility to specify the number of moments (and therefore the complexity/accuracy) to use, but also the assumed size distributions corresponding to these moments. For instance, one might define a 5-moment implementation using an Exponential mode for smaller cloud droplets, plus a Gamma mode for larger rain droplets. Or, more creatively, perhaps a 12-moment implementation comprised of four Gamma modes.

Options for dynamics and size distributions are under continuous development in the `Cloudy.jl` package, thus only the default and suggested use cases are described in detail here.

## Moments and Sub-Distributions

The prognostic variables of this parameterization are a set of N moments, which can be further divided into P sets of moments, each of which correponds to a subdistribution p. By design these moments begin at order 0 and increase as integers up to the maximum number of parameters for the chosen subdistribution. The first three such default moments have interpretable meanings:
  - ``M_0`` - the number density of droplets [1/m^3]
  - ``M_1`` - the mass density of droplets [kg/m^3]
  - ``M_2`` - proportional to the radar reflectivity [kg^2/m^3]
and can be converted to more canonical definitions of `q_liq` and `q_rai` through numerical integration.

When the user wishes to use more than 2 or 3 total variables to represent the system, these moments must be divided between ``P > 1`` sub-distributions, each of which assumes the form of a particular mathematical distribution, such as an Exponential, Lognormal, or Monodisperse (each of which has two parameters), or a Gamma distribution (which takes 3 parameters). 

## Loading the extension
The package `Cloudy.jl` and its dependencies are not loaded by default when using `CloudMicrophysics.jl`. Rather, one must specify:
```
using CloudMicrophysics
using Cloudy
```
from the Julia REPL. Upon recognizing that `Cloudy.jl` is being loaded, the extension `CloudyExt.jl` will then be loaded and overwrite the function stubs defined in `src/MicrophysicsFlexible.jl`.

## Setting up a system
All the details from the number of moments and type of subdistributions, to the parameterizations of coalescence, condensation, and sedimentation are defined through the `CLSetup` (CLoudySetup) mutable struct. This struct is mutable specifically because certain of its components, such as backend-computed coalescence tendencies, are updated prior to being passed to the timestepper. The components of a `CLSetup` object and their defaults are further described below.

|   component         |   description                              |   default                  |
|---------------------|--------------------------------------------|----------------------------|
| ``pdists``          | Vector of subdistributions corresponding   | ``[Exponential, Gamma]``   |
|                     | to the moments                             |                            |
| ``mom``             | Prognostic mass moments, in the same order |``[1e8 / m^3, 1e-2 kg/m^3,``|
|                     | as the corresponding subdistributions;     |``1e6/m^3, 1e-3 kg/m^3,``   |
|                     | first 2 for Exp, next 3 for Gamma          |``2e-12 kg^2/m^3]``         |
| ``KernelFunc``      | Form of the coalescence kernel function    | ``LongKernelFunction``     |
| ``mass_thresholds`` | Particle size thresholds for coalescence   | ``[10.0, Inf]``            |
|                     | integration                                |                            |
| ``kernel order``    | Polynomial order for the approx. kernel    | ``1``                      |
| ``kernel_limit``    | Size threshold for approx. kernel          | ``500``                    |
| ``vel``             | Power-series coefficients for velocity     | ``[2.0, 1/6]``             |
| ``norms``           | Normalizing number density & mass          | ``[1e6/m^3, 1e-9 kg]``     |

# Examples

## Condensation in a parcel
We demonstrate the capabilities and accuracy of this flexible microphysics method in a simple liquid-phase adiabatic parcel in `parcel/Example_NMoment_condensation.jl`, which is modeled after the analogous `parcel/Example_Liquid_only.jl`. Specifically, three different initial particle distributions with equivalent liquid water mass are specified and compared against the canonical solution of Rogers 1975:
  - `Monodisperse`, with 200/cm$^3$ particles and a fixed size of 8$\mu$m
  - `Gamma` distribution, with the same concentration and mean mass, and a shape parameter $k = 2$
  - A `Mixture` of a smaller `Exponential` and a larger `Gamma` distribution, with 90% of the concentration and 50% of the mass in the smaller mode, and the remainder in the larger mode.

Condensation in `Cloudy.jl` proceeds by largely the same mathematical equation as in `CloudMicrophysics.jl`, the the results are expected to be similar. We demonstrate the evolution of the adiabatic parcel as well as the latter two size distributions.

```@example
include("../../parcel/Example_NMoment_condensation.jl")
```
![](cloudy_parcel.svg)
