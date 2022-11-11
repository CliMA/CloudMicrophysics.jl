# Aerosol Nucleation

Homogeneous aerosol nucleation describes the formation of aerosol particles
  from trace gases.
The process occurs when the collisions and clustering
  of individual trace gas molecules outpace the cluster evaporation rate.
It typically involves sulfuric acid ``H_{2}SO_{4}`` and ammonia ``NH_{3}``
  and is a significant source of ultrafine aerosol particles
  (Aitken mode).

The `Nucleation.jl` module contains parameterization
  of binary ``H_{2}SO_{4}-H_{2}O`` nucleation, as described by [Vehkamaki2002](@cite).
This parameterization approximates theoretical nucleation values with empirical fits.

The theoretical values are given by:
``` math
  J = Z \cdot \rho(1,2) \cdot \exp [\frac{-(W^{*}-W(1,2))}{kT}]
```
where:
 - ``\rho(1, 2)`` is the number concentration of sulfuric acid dihydrate,
 - ``W(1, 2)`` is the formation energy of sulfuric acid dihydrate,
 - ``W^{*}`` is the work needed to form the critical cluster,
 - ``Z`` is a kinetic pre-factor.

TODO - add how the above are approximated by the parameterization.
(at least the example functional form)

TODO - what is J and what are the units of all those variables?

## Example of Aerosol Nucleation from Vehkamaki et. al 2002

The figures below have been reproduced from [Vehkamaki2002](@cite) using the
aerosol nucleation parameterization implemented in src/Nucleation.jl.

```@example
include("VehkamakiPlots.jl")
```
![](Vehk_236_K_55.0_P.svg)
![](Vehk_298_K_38.2_P.svg)
![](Vehk_298_K_52.3_P.svg)
