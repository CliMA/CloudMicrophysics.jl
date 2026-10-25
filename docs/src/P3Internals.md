```@meta
CurrentModule = CloudMicrophysics
```

# P3 internals

This page is a reference for the building blocks of the P3 functions on the [API](API.md) page, grouped by what they compute.
Their names, arguments, and behavior can change without notice, so models should call the API functions instead.

## Mass relation

[`P3Scheme.P3Particles`](@ref) evaluates these relations once at construction, from the rime mass fraction and the rime density.
The graupel density follows from the density of the unrimed part, and the size thresholds of the mass regimes follow from the densities.
The aspect ratio of a particle uses the material density of its regime.

```@docs
P3Scheme.get_ρ_d
P3Scheme.get_ρ_g
P3Scheme.get_D_th
P3Scheme.get_D_gr
P3Scheme.get_D_cr
P3Scheme.ϕ_material_density
```

## Slope solve

[`P3Scheme.P3Distribution`](@ref) solves for the slope with `get_distribution_logλ`.
The solve matches `log(L/N)` from `logLdivN`, built from the mass moments of `logmass_gamma_moment` and the gamma moments of `loggamma_moment`, to the ratio of the ice mass and the bounded number.
The shape parameter follows from the slope with `get_μ`, and the intercept from the number with `get_logN₀`.
`logN′ice` evaluates the log of the number concentration at a diameter.

```@docs
P3Scheme.get_distribution_logλ
P3Scheme.logLdivN
P3Scheme.logmass_gamma_moment
P3Scheme.loggamma_moment
P3Scheme.loggamma_inc_moment
P3Scheme.gamma_inc_moment
P3Scheme.get_μ
P3Scheme.get_logN₀
P3Scheme.logN′ice
```

## Liquid-ice collisions

[`P3Scheme.bulk_liquid_ice_collision_sources`](@ref) builds the collision kernel, the local rime density, and the maximum freezing rate for the ice particles of the state.
For each ice diameter, `get_liquid_integrals` integrates the kernel over the liquid size distribution, and `∫liquid_ice_collisions` integrates the result over the ice size distribution and partitions the collected mass between freezing and shedding.

```@docs
P3Scheme.volumetric_collision_rate_integrand
P3Scheme.compute_max_freeze_rate
P3Scheme.compute_local_rime_density
P3Scheme.get_liquid_integrals
P3Scheme.∫liquid_ice_collisions
```

## Quadrature

The P3 integrals split the integration range at the regime thresholds and at the kinks of the integrand, and sum a quadrature over each subinterval.

```@docs
P3Scheme.integrate
P3Scheme.subintervals
P3Scheme.velocity_integral_bounds
P3Scheme.velocity_breakpoints
```
