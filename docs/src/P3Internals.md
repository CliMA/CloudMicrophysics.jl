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
P3Scheme.get_μ
P3Scheme.get_logN₀
P3Scheme.logN′ice
P3Scheme.logλ_bounds
```

## Ventilation

[`P3Scheme.ice_melt`](@ref) reads the ventilation integral of the size distribution from `ice_ventilation_integral`.

```@docs
P3Scheme.ice_ventilation_integral
```

## Liquid-ice collisions

[`P3Scheme.bulk_liquid_ice_collision_sources`](@ref) builds the collision kernel, the local rime density, and the maximum freezing rate for the ice particles of the state.
For each ice diameter, `get_liquid_integrals` integrates the kernel over the liquid size distribution, split at the diameter at which a drop falls as fast as the ice particle.
`∫liquid_ice_collisions` integrates the result over the ice size distribution, split at the wet-growth onset diameters, and partitions the collected mass between freezing and shedding.

```@docs
P3Scheme.collision_cross_section_ice_liquid_coeffs
P3Scheme.collision_cross_section_ice_liquid
P3Scheme.volumetric_collision_rate_integrand
P3Scheme.VolumetricCollisionRate
P3Scheme.compute_max_freeze_rate
P3Scheme.compute_local_rime_density
P3Scheme.CollisionRimeDensity
P3Scheme.rime_density_at
P3Scheme.get_liquid_integrals
P3Scheme.liquid_collision_integrand
P3Scheme.crossover_diameter
P3Scheme.crossing_integral_bounds
P3Scheme.∫liquid_ice_collisions
P3Scheme.wet_growth_onset_diameters
P3Scheme.liquid_ice_collision_quantile
P3Scheme.r_eff_ice
```

## Quadrature

The P3 integrals split the integration range at the regime thresholds and at the kinks of the integrand, and sum a quadrature over each subinterval.

```@docs
P3Scheme.integrate
P3Scheme.subintervals
P3Scheme.velocity_integral_bounds
P3Scheme.insert_breakpoints
P3Scheme.velocity_breakpoints
```

## Lookup tables

[`P3Scheme.P3Tables`](@ref) holds one `UniformTable` per table group, a function with one or more components on a grid uniform along each `TableAxis`, read by multilinear interpolation with `lookup`.
A component that scales as a power of the coordinates is stored as its logarithm.
`tabulate` evaluates a point function at the points of a grid, and `group_table` stores the values of each group in the array type of the values, so a group is converted before the next one is filled.
The rime axes use the coordinate `rime_coordinate`, `integral_logλ_bounds` gives the slopes of the integral tables, and `liquid_mass_bounds` the liquid mean-mass axes of the supported size distributions.

```@docs
P3Scheme.P3Tables
P3Scheme.P3TableGrid
P3Scheme.table_axes
P3Scheme.tabulate
P3Scheme.TableAxis
P3Scheme.axis_point
P3Scheme.axis_points
P3Scheme.grid_point
P3Scheme.GridPointFunction
P3Scheme.evaluate_points
P3Scheme.UniformTable
P3Scheme.inner_axis_table
P3Scheme.component_count
P3Scheme.check_table_size
P3Scheme.component_array
P3Scheme.host_array
P3Scheme.stored_components
P3Scheme.device_array
P3Scheme.lookup
P3Scheme.interpolate
P3Scheme.axis_cell
P3Scheme.uniform_index
P3Scheme.primal_value
P3Scheme.rime_coordinate
P3Scheme.rime_mass_fraction_at
P3Scheme.integral_logλ_bounds
P3Scheme.liquid_mass_bounds
P3Scheme.check_liquid_pdfs
P3Scheme.ice_parameters
P3Scheme.group_table
P3Scheme.table_fingerprint
```

### Inverse slope table

The slope of the size distribution follows from `logLdivN`, which decreases with `logλ`.
`SlopePoint` evaluates `logLdivN` on uniform `logλ` points, and `tabulated_distribution_logλ` inverts the table and refines the result with one Newton step on the closed form.

```@docs
P3Scheme.SlopePoint
P3Scheme.tabulated_distribution_logλ
P3Scheme.bilinear_column
P3Scheme.bisection_steps
```
