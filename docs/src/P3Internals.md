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
`table_groups` gives the point function of each table group, and `tabulate` evaluates it at the points of the axes of the group.
The rime axes use the coordinate `rime_coordinate`, and `integral_logλ_bounds` gives the slopes of the integral tables.

```@docs
P3Scheme.TableAxis
P3Scheme.grid_point
P3Scheme.UniformTable
P3Scheme.lookup
P3Scheme.interpolate
P3Scheme.axis_cell
P3Scheme.uniform_index
P3Scheme.primal_value
P3Scheme.rime_coordinate
P3Scheme.rime_mass_fraction_at
P3Scheme.integral_logλ_bounds
P3Scheme.TABLE_POINTS
P3Scheme.table_groups
P3Scheme.TabulatedP3State
P3Scheme.tabulate
```

### Inverse slope table

The slope of the size distribution follows from `logLdivN`, which decreases with `logλ`.
`slope_point` evaluates `logLdivN` on uniform `logλ` points.
`tabulated_distribution_logλ` interpolates the table bilinearly in the rime coordinate and the rime density, finds the segment of `logλ` that holds the target by a bisection with a fixed number of steps, inverts the target linearly on it, and refines the result with one Newton step on the closed form.

```@docs
P3Scheme.slope_point
P3Scheme.tabulated_distribution_logλ
```

### Ice tables

`ice_point` evaluates the ventilation moment `M_v`, the self-collection rate and the mean fall speeds per unit number, on `logλ`, the rime coordinate, the rime density and the log air density.
`node_self_collection` evaluates the self-collection rate with the node values of each subinterval evaluated once.
`table_state` gives the state at a grid point, `ice_table_key` the coordinates of a state, and `tabulated_ice_integrals` reads the table.

```@docs
P3Scheme.ice_point
P3Scheme.ventilation_moment
P3Scheme.node_self_collection
P3Scheme.table_state
P3Scheme.ice_table_key
P3Scheme.tabulated_ice_integrals
```

### [Cost of the tables](@id P3-cost-of-the-tables)

The build is dominated by the ice tables, whose self-collection integral is a double integral at each of 0.9 million points:

| Group | Points | Build, 16 threads | Build, 1 thread |
|:------|-------:|------------------:|----------------:|
| slope | 32768 | 0.01 s | 0.02 s |
| ice | 912384 | 10 s | 109 s |

These are times in `Float64` on an Apple M4 Max with 12 performance and 4 efficiency cores; in `Float32` the ice tables take 7.7 s on 16 threads and 84 s on one.
On one CPU thread, the rates take

| Function | Quadrature | Tables |
|:---------|-----------:|-------:|
| `P3State`, with the slope | 4.8 μs | 0.60 μs |
| `ice_terminal_velocities` | 4.1 μs | 0.06 μs |
| `ice_self_collection` | 183 μs | 0.06 μs |
| `ice_melt` | 4.6 μs | 0.39 μs |
| `bulk_microphysics_tendencies`, two-moment with P3 ice | 300 μs | 140 μs |

These are `Float32` times with the state built, except for `P3State` and `bulk_microphysics_tendencies`.
The time of `bulk_microphysics_tendencies` is the mean over 400 random states from 235 K to 280 K.
