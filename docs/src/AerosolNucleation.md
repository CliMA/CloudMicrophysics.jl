# Aerosol Nucleation

The `Nucleation.jl` currently contains parameterization
of binary ``H_{2}SO_{4}-H_{2}O`` nucleation, as described by [Vehkamaki2002](@cite).
This parameterization approximates theoretical nucleation values at a fraction of the computational cost.
The theoretical values are given by:
``` math
  J = Z \cdot \rho(1,2) \cdot \exp [\frac{-(W^{*}-W(1,2))}{kT}]
```
where ``ρ(1, 2)`` and ``W(1, 2)`` are, respectively, the number
concentration and the formation energy of a sulfuric acid
dihydrate. ``W^{*}`` is the work needed to form the critical cluster, and Z is a kinetic pre-factor.

## Example of Aerosol Nucleation from Vehkamaki et. al 2002

The figures below have been reproduced from [Vehkamaki2002](@cite) using the
aerosol nucleation parameterization implemented in src/Nucleation.jl.

```@example
include("VehkamakiPlots.jl")
```
![](Vehk_236_K_55.0_P.svg)
![](Vehk_298_K_38.2_P.svg)
![](Vehk_298_K_52.3_P.svg)
