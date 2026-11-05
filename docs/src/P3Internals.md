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
P3Scheme.max_freeze_rate_scalars
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

## [Lookup tables](@id P3-internals-lookup-tables)

The quantities in the tables and their use are described in [Lookup tables](@ref P3-lookup-tables).
This section describes how the tables are built and read, the conventions that every table follows, and the steps to add a table.

### Building the tables

[`P3Scheme.with_tables`](@ref) builds the tables in five steps:

1. [`P3Scheme.P3TableGrid`](@ref) gives the axes of each table group, a tuple of [`P3Scheme.TableAxis`](@ref), with the bounds derived from the parameters and the numbers of points of [`P3Scheme.TABLE_POINTS`](@ref).
2. [`P3Scheme.table_groups`](@ref) gives the point function of each group, which maps the coordinates of one grid point to the values of the group there by quadrature, and the storage options `logs` and `inner` of the group.
3. [`P3Scheme.tabulate`](@ref) evaluates the point function at every point of the axes of the group, on all threads for an `Array`, and by a broadcast where the array type places its elements for any other array type.
4. [`P3Scheme.UniformTable`](@ref) stores the values in the float type of the parameters.
5. [`P3Scheme.P3Tables`](@ref) holds the tables of all groups and a hash of the parameters they were built from, and `with_tables` returns a copy of the parameters with the tables as their `integration`.

Each group is stored before the next one is evaluated, so the build holds the values of one group at a time.

```@docs
P3Scheme.TABLE_POINTS
P3Scheme.table_groups
P3Scheme.TableAxis
P3Scheme.grid_point
P3Scheme.tabulate
P3Scheme.UniformTable
P3Scheme.integral_logλ_bounds
P3Scheme.liquid_mass_bounds
```

### Reading the tables

A [`P3Scheme.TabulatedP3State`](@ref) is a [`P3Scheme.P3State`](@ref) whose `integration` is a `P3Tables`.
The process functions have a method for this type with the signature of the quadrature method, so the callers do not change.
Each method computes the coordinates of the state, reads the values with [`P3Scheme.lookup`](@ref), and multiplies them by the ice number and by the liquid number or mass of the state:

| Method | Groups |
|:-------|:-------|
| `P3Distribution`, in the construction of a `P3State` | `slope` |
| [`P3Scheme.ice_terminal_velocities`](@ref) | `ice` |
| [`P3Scheme.ice_self_collection`](@ref) | `ice` |
| [`P3Scheme.ice_ventilation_integral`](@ref), in [`P3Scheme.ice_melt`](@ref) | `ice` |
| `∫liquid_ice_collisions`, in [`P3Scheme.bulk_liquid_ice_collision_sources`](@ref) | `cloud_collision`, `rain_collision`, `cloud_rime_density`, `rain_rime_density`, and below freezing the wet-set groups |

[`P3Scheme.ice_table_key`](@ref) gives the coordinates `(logλ, w, ρ_rim, logρₐ)` of a state, which all groups on the ice axes share.

```@docs
P3Scheme.TabulatedP3State
P3Scheme.lookup
P3Scheme.interpolate
P3Scheme.axis_cell
P3Scheme.uniform_index
P3Scheme.primal_value
```

### Conventions

- A point function takes a tuple of the coordinates of one point, in the order of the axes of its group, and returns a number or an `SVector` of components.
- The values are per unit ice number, and per unit liquid number or mass, so that a table depends only on the shapes of the size distributions.
- A component that scales as a power of the coordinates is stored as its logarithm, with `logs`, so that the interpolation is accurate and the value stays positive.
- With `inner = true`, the point function evaluates every point of the first axis at once and returns their values as an `SVector`, so that the points share the expensive part of the integral.
  The rime density tables evaluate the points of `logR` in this way, and the fraction tables those of [`P3Scheme.number_logit`](@ref).
- The coordinates are uniform in the logarithm of a quantity that varies over orders of magnitude, `logλ`, `logρₐ`, `logx̄`, `logR` and `logD`.
- The rime coordinate `w = 1 - √(1 - F_rim)` of [`P3Scheme.rime_coordinate`](@ref) replaces `F_rim`, because the particle properties change fastest near `F_rim = 1`.
- [`P3Scheme.ice_table_key`](@ref) takes the primal value of `logλ`, so that a `ForwardDiff.Dual` state is differentiated through the other coordinates only, as in the quadrature.
- [`P3Scheme.lookup`](@ref) moves a coordinate outside an axis to the nearest end of the axis.
- The tables are in the float type of the parameters, and `Adapt` moves them to a device with the parameters.
- The fingerprint of `P3Tables` hashes the parameters that the tables depend on, and the constructors of the parameters from a parameter dictionary check it, see [`P3Scheme.tables_are_current`](@ref).

```@docs
P3Scheme.rime_coordinate
P3Scheme.rime_mass_fraction_at
```

### Adding a table

A new table group needs

1. a family of numbers of points in [`P3Scheme.TABLE_POINTS`](@ref), and the axes of the group in the constructor of `P3TableGrid`,
2. a point function, and its entry in [`P3Scheme.table_groups`](@ref) with `logs` and `inner`,
3. a method of the process function for [`P3Scheme.TabulatedP3State`](@ref) that reads the group with [`P3Scheme.lookup`](@ref),
4. the parameters of the new integrals in `table_fingerprint`, if the existing ones do not cover them,
5. a test against the quadrature in `test/p3_tables_tests.jl`, and the errors of the new quantities in `test/p3_tables_accuracy.jl`.

`test/p3_tables_accuracy.jl` also measures the errors of a proposed grid, which helps to choose the numbers of points.

### Inverse slope table

The tables use the slope ``λ`` as their coordinate, rather than the mean particle mass ``\bar{x} = L/N``.
The shape parameter ``μ`` changes steeply in a narrow range of the mean mass, so a table on the mean mass would need many points there.
The left panel shows the slope table at the rime density 400 kg/m³: the mean mass decreases monotonically with the slope, so the table has an inverse.
The right panel shows the shape parameter against the mean mass, and the range where its slope against ``\log \bar{x}`` exceeds half of its largest value, in gray.

```@example
include("plots/P3SlopeTable.jl")
nothing # hide
```
![](P3SlopeTable.svg)

[`P3Scheme.slope_point`](@ref) evaluates [`P3Scheme.logLdivN`](@ref) on uniform `logλ` points.
[`P3Scheme.tabulated_distribution_logλ`](@ref) interpolates the table bilinearly in the rime coordinate and the rime density, finds the segment of `logλ` that holds the target by a bisection with a fixed number of steps, inverts the target linearly on it, and refines the result with one Newton step on the closed form.

```@docs
P3Scheme.slope_point
P3Scheme.tabulated_distribution_logλ
```

### Ice tables

[`P3Scheme.ice_point`](@ref) evaluates the ventilation moment `M_v`, the self-collection rate and the mean fall speeds per unit number, on `logλ`, the rime coordinate, the rime density and the log air density.
[`P3Scheme.node_self_collection`](@ref) evaluates the self-collection rate with the node values of each subinterval evaluated once.
[`P3Scheme.table_state`](@ref) gives the state of unit number at a grid point, and [`P3Scheme.tabulated_ice_integrals`](@ref) reads the table.

```@docs
P3Scheme.ice_point
P3Scheme.ventilation_moment
P3Scheme.node_self_collection
P3Scheme.table_state
P3Scheme.ice_table_key
P3Scheme.tabulated_ice_integrals
```

### Liquid-ice tables

[`P3Scheme.collision_point`](@ref) evaluates the collected number and mass of cloud or rain per unit ice number and per unit liquid number or mass, on the four ice axes and the liquid mean mass, without the freeze and shed partition.
A liquid size distribution is its number times a shape that depends only on its mean mass, see [`P3Scheme.liquid_table_scales`](@ref), so the collected rates per unit liquid number or mass depend on the liquid only through `logx̄`.
The local rime density depends on the temperature through the ratio `Rᵢ` of the impact scale to the temperature, so [`P3Scheme.rime_density_point`](@ref) evaluates the mean rime density on the log of that ratio at the mean impact scale `log κ̄` of the collected mass.
[`P3Scheme.tabulated_collection`](@ref) reads the tables for a state, at `logR = log κ̄ - log|T - T_freeze|`.

```@docs
P3Scheme.collision_point
P3Scheme.rime_density_point
P3Scheme.collection_integrals
P3Scheme.particle_collection
P3Scheme.impact_scale
P3Scheme.table_liquid
P3Scheme.liquid_table_scales
P3Scheme.table_log_temperature
P3Scheme.tabulated_collection
```

### Wet-set tables

An ice particle of diameter ``D`` collects the liquid mass ``u(D) = L_c\, \hat{u}_c(D) + L_r\, \hat{u}_r(D)`` per unit time, with ``\hat{u}_l`` the collected mass per unit liquid mass, and it freezes at most
```math
m(D) = c_1 D + c_2 D^{3/2} v(D)^{1/2},
```
the maximum freezing rate ``∂_t\mathcal{M}_\text{max}`` of [Collisions with liquid droplets](@ref P3-collisions).
The coefficients ``c_1`` and ``c_2`` hold the temperature, the air properties and the ventilation coefficients, see [`P3Scheme.max_freeze_rate_scalars`](@ref).
The particles with ``u > m`` form the wet set ``W``, which freezes ``∫_W m N' \mathrm{d}D`` and sheds ``∫_W (u - m) N' \mathrm{d}D``, and the other particles freeze all they collect.

[`P3Scheme.tabulated_wet_set`](@ref) evaluates this partition in four steps:

1. A scan of 25 points uniform in ``\log D``, between the bounds of the collision integrals from the table of [`P3Scheme.collision_bounds_point`](@ref), locates the ends of ``W`` as the zeros of ``\log(u / m)``.
   ``\hat{u}_l`` comes from the per-particle tables of [`P3Scheme.collected_point`](@ref), and ``D^{3/2} v^{1/2}`` from the table of [`P3Scheme.capacity_point`](@ref).
2. The fraction tables of [`P3Scheme.fraction_point`](@ref) and [`P3Scheme.capacity_fraction_point`](@ref) give the fractions of ``∫ \hat{u}_l N' \mathrm{d}D`` and of ``M_v = ∫ D^{3/2} v^{1/2} N' \mathrm{d}D`` below each end of ``W``, on the number logit of the end, see [`P3Scheme.number_logit`](@ref).
   The fraction of the first moment ``∫ D N' \mathrm{d}D`` has a closed form in the incomplete gamma function.
   The frozen mass is the mass collected outside ``W`` plus ``c_1 ∫_W D N' \mathrm{d}D + c_2 ∫_W D^{3/2} v^{1/2} N' \mathrm{d}D``, and the shed mass the mass collected on ``W`` less the frozen mass on ``W``.
   The smaller of the frozen and the shed mass is computed, and the larger is the collected mass less the smaller, which keeps the relative precision of the smaller.
3. The trapezoid rule on the scan points, with the frozen fraction ``\min(m / u, 1)`` at each point, divides the smaller mass between cloud and rain.
4. The rime volume of each species is its frozen mass over its mean local rime density.
   The trapezoid rule gives the ratio of the mean inverse rime density of the frozen mass to that of the collected mass, with the local rime density of each point at the mean impact scale of the mass that a particle collects there, and this ratio scales the rime density table.

```@docs
P3Scheme.collision_bounds_point
P3Scheme.collected_point
P3Scheme.capacity_point
P3Scheme.fraction_point
P3Scheme.capacity_fraction_point
P3Scheme.cumulative_logits
P3Scheme.number_logit
P3Scheme.logistic_fractions
P3Scheme.wet_set_scan_intervals
P3Scheme.tabulated_wet_set
P3Scheme.collision_diameter_bounds
P3Scheme.collision_logξ_lo
```

### [Cost of the tables](@id P3-cost-of-the-tables)

The build is dominated by the liquid-ice collision tables, a double integral at each of 6.2 million points, and by the ice tables, whose self-collection integral is a double integral at each of 0.9 million points:

| Group | Points | Build, 16 threads | Build, 1 thread |
|:------|-------:|------------------:|----------------:|
| slope | 32768 | 0.01 s | 0.02 s |
| ice | 912384 | 10 s | 109 s |
| cloud collisions | 3133440 | 14 s | 151 s |
| rain collisions | 3041280 | 25 s | 273 s |
| rime density | 716800 | 0.9 s | 5.1 s |
| wet set | 2741824 | 4.3 s | 40 s |
| all, with the assembly | | 52 s | 577 s |

These are times in `Float64`, measured in November 2026 on an Apple M4 Max with 12 performance and 4 efficiency cores; in `Float32` the build takes 70 s on 16 threads and 515 s on one.
On one thread of the same machine, the rates take

| Function | Quadrature | Tables |
|:---------|-----------:|-------:|
| `P3State`, with the slope | 4.7 μs | 0.63 μs |
| `ice_terminal_velocities` | 4.2 μs | 0.09 μs |
| `ice_self_collection` | 183 μs | 0.09 μs |
| `ice_melt` | 4.6 μs | 0.42 μs |
| `∫liquid_ice_collisions`, below freezing | 152 μs | 2.4 μs |
| `∫liquid_ice_collisions`, above freezing | 127 μs | 0.26 μs |
| `bulk_microphysics_tendencies`, two-moment with P3 ice | 300 μs | 3.6 μs |

These are `Float32` times with the state built, except for `P3State` and `bulk_microphysics_tendencies`.
The time of `bulk_microphysics_tendencies` is the mean over 400 random states from 235 K to 280 K.
