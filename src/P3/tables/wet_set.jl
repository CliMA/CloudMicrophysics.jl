"""
    CollisionBoundsPoint()

The point function of the bounds of the collision integrals and of the partition over the wet set at a point `(μ,)`.

The value is the `SVector` of the log of the bounds of the scaled diameter `ξ = D λ` of
[`collision_logξ_lo`](@ref) and `collision_logξ_hi`, and of [`partition_logξ_bounds`](@ref).
"""
struct CollisionBoundsPoint end
(::CollisionBoundsPoint)((μ,)) = SA.SVector(collision_logξ_lo(μ), collision_logξ_hi(μ), partition_logξ_bounds(μ)...)

"""
    partition_logξ_bounds(μ)

Return the log of the bounds of the scaled diameter `ξ = D λ` that hold the bulk of the integrals of the freeze and shed partition over the wet set, see [`tabulated_wet_set`](@ref).

The integrands of the partition scale as `D^a` times the size distribution, with `a` from 1 for
the maximum freezing rate of small particles to about 3 for the collected mass of large ones,
so they are gamma distributions in `ξ` of the shapes from `μ + 2` to `μ + 4`. The bounds are the
quantile `1e-4` of the shape `μ + 2` and the quantile `1 - 1e-4` of the shape `μ + 4`.
"""
function partition_logξ_bounds(μ)
    p = 1e-4
    return (log(UT.gamma_inc_inv(μ + 2, p, 1 - p)), log(UT.gamma_inc_inv(μ + 4, 1 - p, p)))
end

"""
    CollectedPoint{P, V, L, G}

The point function of the per-particle collection table of one liquid species at a point `(logD, w, ρ_rim, logρₐ, logx̄)`.

The value is an `SVector` of the mass that one ice particle of diameter `exp(logD)` collects per
unit liquid mass [m³/s] and the impact scale `log κ̄` of the collected mass, the mean of `log κ`
over it (see [`impact_scale`](@ref)), from [`particle_collection`](@ref) with the rule
`generator`. It does not depend on the temperature.

# Construction

    CollectedPoint(params, terminal_velocity, pdf, generator)

# Fields
$(FIELDS)
"""
struct CollectedPoint{P, V, L, G}
    "The P3 parameters, [`CMP.ParametersP3`](@ref)"
    params::P
    "The terminal velocity, [`CMP.Chen2022VelType`](@ref)"
    terminal_velocity::V
    "The size distribution of the liquid species"
    pdf::L
    "The `QuadratureRule` of the integrals"
    generator::G
end
function (f::CollectedPoint)((logD, w, ρ_rim, logρₐ, logx̄))
    particles = P3Particles(f.params, rime_mass_fraction_at(w), ρ_rim)
    liquid = table_liquid(f.terminal_velocity, particles, f.pdf, exp(logρₐ), logx̄)
    r = particle_collection(f.params, particles, liquid, (), f.generator)(exp(logD))
    return SA.SVector(r[2] / liquid.mass, r[3] / r[2])
end

"""
    ParticleRimeDensityPoint{K, P, V, L, G}

The point function of the rime density table of one particle and one liquid species at a point `(logD, w, ρ_rim, logρₐ, logx̄)`.

The value is an `SVector` of the mean local rime density [kg/m³] of the mass that one ice particle
of diameter `exp(logD)` collects, the ratio of its collected mass rate to its rime volume rate of
[`particle_collection`](@ref), at each of the `K` points of `logR`. A point of `logR` stands for
the temperature `T = -exp(log κ̄ - logR)` [°C] at the impact scale `log κ̄` of the mass that the
particle collects.

# Construction

    ParticleRimeDensityPoint(params, terminal_velocity, pdf, generator, logR)

with `logR` the `NTuple` of the points of the `logR` axis.

# Fields
$(FIELDS)
"""
struct ParticleRimeDensityPoint{K, P, V, L, G}
    "The P3 parameters, [`CMP.ParametersP3`](@ref)"
    params::P
    "The terminal velocity, [`CMP.Chen2022VelType`](@ref)"
    terminal_velocity::V
    "The size distribution of the liquid species"
    pdf::L
    "The `QuadratureRule` of the integrals"
    generator::G
    "The points of the `logR` axis"
    logR::NTuple{K, Float64}
end
function (f::ParticleRimeDensityPoint{K})((logD, w, ρ_rim, logρₐ, logx̄)) where {K}
    particles = P3Particles(f.params, rime_mass_fraction_at(w), ρ_rim)
    liquid = table_liquid(f.terminal_velocity, particles, f.pdf, exp(logρₐ), logx̄)
    D = exp(logD)
    rates = particle_collection(f.params, particles, liquid, (), f.generator)(D)
    log_impact = rates[3] / rates[2]
    T°Cs = ntuple(k -> -exp(log_impact - f.logR[k]), Val(K))
    r = particle_collection(f.params, particles, liquid, T°Cs, f.generator)(D)
    return SA.SVector(ntuple(k -> r[2] / r[k + 3], Val(K)))
end

"""
    CapacityPoint{P, V}

The point function of the capacity table, `D^(3/2) v(D)^(1/2)` [m² s^(-1/2)] of one ice particle at a point `(logD, w, ρ_rim, logρₐ)`.

# Construction

    CapacityPoint(params, terminal_velocity)

# Fields
$(FIELDS)
"""
struct CapacityPoint{P, V}
    "The P3 parameters, [`CMP.ParametersP3`](@ref)"
    params::P
    "The terminal velocity, [`CMP.Chen2022VelType`](@ref)"
    terminal_velocity::V
end
function (f::CapacityPoint)((logD, w, ρ_rim, logρₐ))
    particles = P3Particles(f.params, rime_mass_fraction_at(w), ρ_rim)
    D = exp(logD)
    v = ice_particle_terminal_velocity(f.terminal_velocity, particles, exp(logρₐ))(D)
    return D * sqrt(D * v)
end

"""
    FractionPoint{K, P, V, L, G}

The point function of the fraction table of one liquid species, the logits of the fractions of the collected mass below the diameters of the number logits at a point `(logλ, w, ρ_rim, logρₐ, logx̄)`.

The value is the `SVector` of [`cumulative_logits`](@ref) of the collected mass rate of one ice
particle, see [`CollectedPoint`](@ref), at the `K` number logits `logits`, for the state of
[`table_state`](@ref).

# Construction

    FractionPoint(params, terminal_velocity, pdf, generator, logits)

with `logits` the `NTuple` of the points of the `number_logit` axis.

# Fields
$(FIELDS)
"""
struct FractionPoint{K, P, V, L, G}
    "The P3 parameters, [`CMP.ParametersP3`](@ref)"
    params::P
    "The terminal velocity, [`CMP.Chen2022VelType`](@ref)"
    terminal_velocity::V
    "The size distribution of the liquid species"
    pdf::L
    "The `QuadratureRule` of the integrals"
    generator::G
    "The points of the `number_logit` axis"
    logits::NTuple{K, Float64}
end
function (f::FractionPoint)((logλ, w, ρ_rim, logρₐ, logx̄))
    (; particles, dist) = table_state(f.params, f.generator, logλ, w, ρ_rim)
    liquid = table_liquid(f.terminal_velocity, particles, f.pdf, exp(logρₐ), logx̄)
    rates_at = particle_collection(f.params, particles, liquid, (), f.generator)
    return cumulative_logits(D -> rates_at(D)[2], particles, dist, liquid.v_i, f.logits, f.generator)
end

"""
    CapacityFractionPoint{K, P, V, G}

The point function of the capacity fraction table, the logits of the fractions of the ventilation moment `M_v` below the diameters of the number logits at a point `(logλ, w, ρ_rim, logρₐ)`.

The value is the `SVector` of [`cumulative_logits`](@ref) of `D^(3/2) v(D)^(1/2)`, see
[`CapacityPoint`](@ref), at the `K` number logits `logits`, for the state of
[`table_state`](@ref).

# Construction

    CapacityFractionPoint(params, terminal_velocity, generator, logits)

# Fields
$(FIELDS)
"""
struct CapacityFractionPoint{K, P, V, G}
    "The P3 parameters, [`CMP.ParametersP3`](@ref)"
    params::P
    "The terminal velocity, [`CMP.Chen2022VelType`](@ref)"
    terminal_velocity::V
    "The `QuadratureRule` of the integrals"
    generator::G
    "The points of the `number_logit` axis"
    logits::NTuple{K, Float64}
end
function (f::CapacityFractionPoint)((logλ, w, ρ_rim, logρₐ))
    (; particles, dist) = table_state(f.params, f.generator, logλ, w, ρ_rim)
    v_i = ice_particle_terminal_velocity(f.terminal_velocity, particles, exp(logρₐ))
    return cumulative_logits(D -> D * sqrt(D * v_i(D)), particles, dist, v_i, f.logits, f.generator)
end

"""
    group_table(::Val{name}, values, grid, FT)

Return a table of the freeze and shed partition of the liquid-ice collisions, for `name` in

- `:collision_bounds`: the log of the bounds of the scaled diameter `ξ = D λ` of the collision
  integrals on the shape parameter `μ`, the values of [`CollisionBoundsPoint`](@ref),
- `:cloud_collected`, `:rain_collected`: the log of the mass a particle collects per unit liquid
  mass, `u / L` [m³/s], and its impact scale `log κ̄`, on `(logD, w, ρ_rim, logρₐ, logx̄)`, the
  values of [`CollectedPoint`](@ref),
- `:capacity`: the log of `D^(3/2) v(D)^(1/2)` [m² s^(-1/2)] on `(logD, w, ρ_rim, logρₐ)`, the
  values of [`CapacityPoint`](@ref),
- `:cloud_fraction`, `:rain_fraction`: the logit of the fraction of the collected mass that the
  particles below a diameter collect, on `(number_logit, logλ, w, ρ_rim, logρₐ, logx̄)`, the
  values of [`FractionPoint`](@ref),
- `:capacity_fraction`: the logit of the fraction of the ventilation moment `M_v` below a
  diameter, on `(number_logit, logλ, w, ρ_rim, logρₐ)`, the values of
  [`CapacityFractionPoint`](@ref),
- `:cloud_particle_rime_density`, `:rain_particle_rime_density`: the mean local rime density of
  the mass a particle collects [kg/m³] on `(logR, logD, w, ρ_rim, logρₐ, logx̄)`, the values of
  [`ParticleRimeDensityPoint`](@ref).

The axes are those of `grid.wet_set` and `grid.particle_rime_density`, see [`table_axes`](@ref).
An ice particle of diameter `D` freezes the collected liquid mass `u(D) + v(D)`, from cloud and
rain, up to its maximum freezing rate `m(D)` of [`compute_max_freeze_rate`](@ref), and sheds the
rest. The shed mass rate is then an integral over the wet set `W`, where `u + v > m`,

    S = ∫_W (u(D) + v(D) - m(D)) N′(D) dD.

The maximum freezing rate is `m(D) = 2π A (a + b D^(1/2) v(D)^(1/2)) D`, with `A = heat / denom`
of [`max_freeze_rate_scalars`](@ref) and the ventilation coefficients `(a, b)` of
[`CO.ventilation_coefficients`](@ref), so the air properties enter only as factors. Once the ends
of `W` are known, `S` is linear in the integrals of `u`, `v`, `D` and `D^(3/2) v^(1/2)` over `W`.
A diameter enters the fraction tables by the logit of the number fraction of the ice below it,
see [`number_logit`](@ref), which spreads the bulk of every distribution over the axis. The
integral of `D` has a closed form. [`tabulated_wet_set`](@ref) reads the tables.
"""
group_table(::Val{:collision_bounds}, values, grid, ::Type{FT}) where {FT} =
    UniformTable(values, table_axes(grid).collision_bounds, FT)
group_table(::Val{:cloud_collected}, values, grid, ::Type{FT}) where {FT} =
    UniformTable(values, table_axes(grid).cloud_collected, FT; logs = (true, false))
group_table(::Val{:rain_collected}, values, grid, ::Type{FT}) where {FT} =
    UniformTable(values, table_axes(grid).rain_collected, FT; logs = (true, false))
group_table(::Val{:capacity}, values, grid, ::Type{FT}) where {FT} =
    UniformTable(values, table_axes(grid).capacity, FT; logs = (true,))
group_table(::Val{:cloud_fraction}, values, grid, ::Type{FT}) where {FT} =
    inner_axis_table(values, grid.wet_set.number_logit, table_axes(grid).cloud_fraction, FT)
group_table(::Val{:rain_fraction}, values, grid, ::Type{FT}) where {FT} =
    inner_axis_table(values, grid.wet_set.number_logit, table_axes(grid).rain_fraction, FT)
group_table(::Val{:capacity_fraction}, values, grid, ::Type{FT}) where {FT} =
    inner_axis_table(values, grid.wet_set.number_logit, table_axes(grid).capacity_fraction, FT)
group_table(::Val{:cloud_particle_rime_density}, values, grid, ::Type{FT}) where {FT} =
    inner_axis_table(values, grid.particle_rime_density.logR, table_axes(grid).cloud_particle_rime_density, FT)
group_table(::Val{:rain_particle_rime_density}, values, grid, ::Type{FT}) where {FT} =
    inner_axis_table(values, grid.particle_rime_density.logR, table_axes(grid).rain_particle_rime_density, FT)

"""
    cumulative_logits(k, particles, dist, v_i, logits, generator)

Return the `SVector` of the logits of the fractions of `∫ k(D) N′(D) dD` below the diameters at the number logits `logits`.

A number logit gives the diameter below which that fraction of the ice number lies, see
[`number_logit`](@ref).

The integral extends over the bounds of the liquid-ice collisions, see
[`liquid_ice_collision_quantile`](@ref). It is the sum of the integrals of `generator` between
consecutive diameters, the diameters of the number logits, the thresholds of `particles` and the
breakpoints of the terminal velocity `v_i`. The logit of a fraction `F` is `log(F) - log(1 -
F)`, with both terms from sums of intervals, so that neither tail is the difference of two
numbers near the total, and with each term at least `1e-30` of the total.

# Arguments
- `k`: function of the ice diameter `D` [m]
- `particles`: [`P3Particles`](@ref)
- `dist`: [`P3Distribution`](@ref)
- `v_i`: ice terminal velocity, see [`ice_particle_terminal_velocity`](@ref)
- `logits`: `NTuple` of increasing number logits
- `generator`: `QuadratureRule`
"""
function cumulative_logits(k::FK, particles, dist, v_i::V, logits::NTuple{K}, generator) where {FK, V, K}
    n_i = DT.size_distribution(dist)
    p = liquid_ice_collision_quantile(Float64)
    bounds = velocity_integral_bounds(particles, dist, v_i; p)
    D_lo, D_hi = first(bounds), last(bounds)
    breaks = SA.sort(SA.SVector(particles.D_th, particles.D_gr, particles.D_cr, velocity_breakpoints(v_i)...))
    cuts = ntuple(Val(K)) do j
        (F, Q) = logistic_fractions(logits[j])
        return clamp(UT.gamma_inc_inv(dist.μ + 1, F, Q) * exp(-dist.logλ), D_lo, D_hi)
    end
    edges = (D_lo, cuts..., D_hi)
    # The integral between consecutive edges, split at the breakpoints between them
    pieces = SA.SVector(
        ntuple(Val(K + 1)) do i
            a, b = edges[i], edges[i + 1]
            return integrate(D -> n_i(D) * k(D), (a, map(D -> clamp(D, a, b), Tuple(breaks))..., b), generator)
        end,
    )
    smallest = 1e-30 * sum(pieces)
    return SA.SVector(ntuple(Val(K)) do j
        below = sum(ntuple(i -> ifelse(i ≤ j, pieces[i], zero(eltype(pieces))), Val(K + 1)))
        above = sum(ntuple(i -> ifelse(i > j, pieces[i], zero(eltype(pieces))), Val(K + 1)))
        return log(max(below, smallest)) - log(max(above, smallest))
    end)
end

"""
    number_logit(μ, ξ)

Return the logit `log(P) - log(1 - P)` of the number fraction `P` of the ice below the scaled diameter `ξ = D λ`.

`P = P(μ + 1, ξ)` is the regularized lower incomplete gamma function, see [`UT.gamma_inc`](@ref).
The bounds of the liquid-ice collision integrals are at `P = p` and `P = 1 - p`, see
[`liquid_ice_collision_quantile`](@ref).

# Arguments
- `μ`: shape parameter [-]
- `ξ`: scaled diameter [-]
"""
@inline function number_logit(μ, ξ)
    (P, Q) = UT.gamma_inc(μ + 1, ξ)
    return log(max(P, floatmin(P))) - log(max(Q, floatmin(Q)))
end

"""
    logistic_fractions(z)

Return the fraction `F = 1 / (1 + exp(-z))` of the logit `z`, and `1 - F`.

The logit is limited to `±50`, where `F` or `1 - F` is below `2e-22`, so an infinite logit gives a
finite derivative.
"""
@inline function logistic_fractions(z)
    z = clamp(z, -50, 50)
    return (inv(1 + exp(-z)), inv(1 + exp(z)))
end

"""
    fraction_between(z_a, z_b)

Return `F(z_b) - F(z_a)` for the logits `z_a ≤ z_b` of two fractions `F`, see [`logistic_fractions`](@ref).

The difference is formed from the fractions below the upper logit if it is negative, and from
the fractions above the lower one otherwise, so that no tail is the difference of two numbers
near one.
"""
@inline function fraction_between(z_a, z_b)
    (F_a, Q_a) = logistic_fractions(z_a)
    (F_b, Q_b) = logistic_fractions(z_b)
    return max(ifelse(z_b < 0, F_b - F_a, Q_a - Q_b), zero(F_a))
end

"""
    fixed_cells(t, first_axis, x)

Return the cells of the coordinates `x` on the axes of the [`UniformTable`](@ref) `t` from `first_axis` on, see [`axis_cell`](@ref).
"""
@inline fixed_cells(t, first_axis, x) = ntuple(d -> axis_cell(t, first_axis + d - 1, x[d]), Val(length(x)))

"""
    tabulated_wet_set(tables, state, key, cloud, rain, coeffs, logT)

Return the freeze and shed partition of the liquid-ice collisions per unit ice number from the wet-set tables.

As in the quadrature, each particle freezes the fraction `f = min(m / (u + v), 1)` of the mass it
collects from cloud and from rain, with `m` its maximum freezing rate. The wet set `W`, where
`f < 1`, is located as in [`wet_growth_onset_diameters`](@ref), with the ratio of the collected
mass to the maximum freezing rate of one particle from the per-particle tables (see
[`group_table`](@ref)), and its ends follow from primal values (see [`primal_value`](@ref)).

The total frozen and shed masses follow from the fraction tables. The frozen mass is the mass
collected on the dry set plus `c₁ M_1 + c₂ M_v` on `W`, and the shed mass is `U_c + U_v - c₁ M_1 -
c₂ M_v`, with the integrals over `W`. The smaller of the two is used, and the larger is the
collected mass less the smaller, so neither is a small difference of two large numbers. The
division of the smaller between cloud and rain, and the mean local rime density of the frozen
mass of each species, follow from the integrals of the per-particle tables by the rule
`tables.partition_rule` in `log D` on each interval between the ends of `W`, within the bulk of
the integrands, see [`partition_logξ_bounds`](@ref). The rime volume of each species is its
frozen mass divided by that mean density, scaled so that it is the rime volume of the
rime density table where every particle freezes, see [`ParticleRimeDensityPoint`](@ref).

The result is a `NamedTuple` with the frozen mass `cloud_frozen`, `rain_frozen`, the shed mass
`cloud_shed`, `rain_shed`, and the rime volume of the frozen mass `cloud_rime`, `rain_rime`, of
cloud and rain, all per unit ice number.

# Arguments
- `tables`: [`P3Tables`](@ref)
- `state`: [`P3State`](@ref)
- `key`: the ice coordinates of [`ice_table_key`](@ref)
- `cloud`, `rain`: the collected rates of [`tabulated_collection`](@ref)
- `coeffs`: the coefficients `(c₁, c₂)` of the maximum freezing rate `c₁ D + c₂ D^(3/2) v(D)^(1/2)`
  of one particle
- `logT`: the temperature coordinate of [`table_log_temperature`](@ref)
"""
@inline function tabulated_wet_set(tables, state, key, cloud, rain, (c₁, c₂), logT)
    FT = eltype(tables)
    (; groups, partition_rule) = tables
    μ = FT(primal_value(state.dist.μ))
    logλ = FT(primal_value(state.dist.logλ))
    (_, w, ρ_rim, logρₐ) = key
    (_, w_p, ρ_rim_p, logρₐ_p) = map(x -> oftype(μ, primal_value(x)), key)
    t_c, t_r, t_cap = groups.cloud_collected, groups.rain_collected, groups.capacity
    t_ρc, t_ρr = groups.cloud_particle_rime_density, groups.rain_particle_rime_density
    L_c_p, L_r_p, c₁_p, c₂_p = map(x -> oftype(μ, primal_value(x)), (cloud.mass, rain.mass, c₁, c₂))

    # The log ratio of the collected mass to the maximum freezing rate of a particle of diameter `exp(logD)`
    cells_c_p = fixed_cells(t_c, 2, (w_p, ρ_rim_p, logρₐ_p, FT(primal_value(cloud.logx̄))))
    cells_r_p = fixed_cells(t_r, 2, (w_p, ρ_rim_p, logρₐ_p, FT(primal_value(rain.logx̄))))
    cells_cap_p = fixed_cells(t_cap, 2, (w_p, ρ_rim_p, logρₐ_p))
    smallest = floatmin(FT)
    function log_ratio(logD)
        cell_D = axis_cell(t_c, 1, logD)
        u =
            L_c_p * exp(interpolate(t_c, (cell_D, cells_c_p...))[1]) +
            L_r_p * exp(interpolate(t_r, (cell_D, cells_r_p...))[1])
        m = c₁_p * exp(logD) + c₂_p * exp(interpolate(t_cap, (cell_D, cells_cap_p...))[1])
        return log(max(u, smallest)) - log(m)
    end
    (logξ_lo, logξ_hi, logξ_bulk_lo, logξ_bulk_hi) = lookup(groups.collision_bounds, μ)
    logD_lo, logD_hi = logξ_lo - logλ, logξ_hi - logλ
    (; roots, n_found, f_lo) = UT.find_zeros(
        log_ratio, RS.RegulaFalsiMethod, logD_lo, logD_hi, Val(16), Val(2), FixedIterations{FT}(), 3,
    )
    # The three intervals between the bounds and the zeros are wet and dry in turn, and an
    # interval past the last zero is empty
    wet_lo = f_lo > 0
    wet = (wet_lo, !wet_lo, wet_lo)
    found = (n_found ≥ 1, n_found ≥ 2)
    ends = (logD_lo, ifelse(found[1], roots[1], logD_hi), ifelse(found[2], roots[2], logD_hi), logD_hi)

    # The integrals over W and the dry set from the fraction tables, with the logits -Inf and Inf at the bounds
    ξ = map(r -> exp(r + logλ), roots)
    z = map((f, ξ) -> ifelse(f, number_logit(μ, ξ), oftype(μ, Inf)), found, ξ)
    function fractions(t, x)
        F = map((f, z) -> (v = lookup(t, z, x...)[1]; ifelse(f, v, oftype(v, Inf))), found, z)
        zs = (-oftype(F[1], Inf), F..., oftype(F[1], Inf))
        on(set) = sum(ntuple(i -> ifelse(set[i], fraction_between(zs[i], zs[i + 1]), zero(F[1])), Val(3)))
        return (; wet = on(wet), dry = on(map(!, wet)))
    end
    cloud_fraction = fractions(groups.cloud_fraction, (key..., cloud.logx̄))
    rain_fraction = fractions(groups.rain_fraction, (key..., rain.logx̄))
    M_v = fractions(groups.capacity_fraction, key).wet * exp(lookup(groups.ice, key...)[1])
    # The first moment in closed form, (μ + 1) / λ times the regularized incomplete gamma function P(μ + 2, ξ)
    PQ = map(ξ -> UT.gamma_inc(μ + 2, ξ), ξ)
    P = (zero(FT), ifelse(found[1], first(PQ[1]), one(FT)), ifelse(found[2], first(PQ[2]), one(FT)), one(FT))
    Q = (one(FT), ifelse(found[1], last(PQ[1]), zero(FT)), ifelse(found[2], last(PQ[2]), zero(FT)), zero(FT))
    between(i) = max(ifelse(P[i + 1] < one(μ) / 2, P[i + 1] - P[i], Q[i] - Q[i + 1]), zero(μ))
    M_1 = sum(ntuple(i -> ifelse(wet[i], between(i), zero(μ)), Val(3))) * (μ + 1) * exp(-logλ)
    frozen_on_W = c₁ * M_1 + c₂ * M_v
    frozen = cloud_fraction.dry * cloud.M + rain_fraction.dry * rain.M + frozen_on_W
    shed = max(cloud_fraction.wet * cloud.M + rain_fraction.wet * rain.M - frozen_on_W, zero(frozen_on_W))

    # The integrands of the partition per unit log D, up to a constant factor
    cells_c = fixed_cells(t_c, 2, (w, ρ_rim, logρₐ, cloud.logx̄))
    cells_r = fixed_cells(t_r, 2, (w, ρ_rim, logρₐ, rain.logx̄))
    cells_cap = fixed_cells(t_cap, 2, (w, ρ_rim, logρₐ))
    cells_ρc = fixed_cells(t_ρc, 3, (w, ρ_rim, logρₐ, cloud.logx̄))
    cells_ρr = fixed_cells(t_ρr, 3, (w, ρ_rim, logρₐ, rain.logx̄))
    function partition(logD, on_wet_set)
        cell_D = axis_cell(t_c, 1, logD)
        (log_u_c, log_κ_c) = interpolate(t_c, (cell_D, cells_c...))
        (log_u_r, log_κ_r) = interpolate(t_r, (cell_D, cells_r...))
        u_c = cloud.mass * exp(log_u_c)
        u_r = rain.mass * exp(log_u_r)
        m = c₁ * exp(logD) + c₂ * exp(interpolate(t_cap, (cell_D, cells_cap...))[1])
        f = ifelse(on_wet_set, min(UT.guarded_quotient(m, u_c + u_r), one(m)), one(m))
        cell_ρD = axis_cell(t_ρc, 2, logD)
        ρ′_c = interpolate(t_ρc, (axis_cell(t_ρc, 1, log_κ_c - logT), cell_ρD, cells_ρc...))[1]
        ρ′_r = interpolate(t_ρr, (axis_cell(t_ρr, 1, log_κ_r - logT), cell_ρD, cells_ρr...))[1]
        # D N′(D), up to a constant factor, since dD = D d(log D)
        g = exp((μ + 1) * (logD + logλ) - exp(logD + logλ))
        return SA.SVector(
            g * u_c * f, g * u_c * (1 - f), g * u_r * f, g * u_r * (1 - f),
            g * u_c * f / ρ′_c, g * u_c / ρ′_c, g * u_c, g * u_r * f / ρ′_r, g * u_r / ρ′_r, g * u_r,
        )
    end
    lo, hi = logξ_bulk_lo - logλ, logξ_bulk_hi - logλ
    sums = zero(
        SA.SVector{
            10,
            Base.promote_typeof(cloud.mass, rain.mass, c₁, c₂, w, ρ_rim, logρₐ, cloud.logx̄, rain.logx̄, logT),
        },
    )
    for i in 1:3
        a = clamp(ends[i], lo, hi)
        b = clamp(ends[i + 1], a, hi)
        h = (b - a) / 2
        for k in eachindex(partition_rule.nodes)
            sums += h * partition_rule.weights[k] * partition(h * partition_rule.nodes[k] + (a + b) / 2, wet[i])
        end
    end
    (cloud_f, cloud_s, rain_f, rain_s, cloud_fB, cloud_B, cloud_u, rain_fB, rain_B, rain_u) = sums
    # The smaller of the frozen and the shed mass, divided between cloud and rain
    freezes_less = frozen < shed
    smaller = min(frozen, shed, cloud.M + rain.M)
    # A wet set outside the bulk of the integrands sheds in proportion to the mass collected on it
    U_c, U_r = cloud_fraction.wet * cloud.M, rain_fraction.wet * rain.M
    shed_share = ifelse(
        cloud_s + rain_s > 0,
        UT.guarded_quotient(cloud_s, cloud_s + rain_s),
        UT.guarded_quotient(U_c, U_c + U_r),
    )
    cloud_share = ifelse(freezes_less, UT.guarded_quotient(cloud_f, cloud_f + rain_f), shed_share)
    cloud_smaller = clamp(smaller * cloud_share, max(smaller - rain.M, zero(smaller)), min(smaller, cloud.M))
    rain_smaller = smaller - cloud_smaller
    cloud_frozen = ifelse(freezes_less, cloud_smaller, cloud.M - cloud_smaller)
    rain_frozen = ifelse(freezes_less, rain_smaller, rain.M - rain_smaller)
    # The mean inverse local rime density of the frozen mass, relative to that of the collected mass
    relative_inverse_density(fB, f, B, u) = UT.guarded_quotient(UT.guarded_quotient(fB, f), UT.guarded_quotient(B, u))
    return (;
        cloud_frozen,
        rain_frozen,
        cloud_shed = cloud.M - cloud_frozen,
        rain_shed = rain.M - rain_frozen,
        cloud_rime = cloud_frozen / cloud.ρ′ * relative_inverse_density(cloud_fB, cloud_f, cloud_B, cloud_u),
        rain_rime = rain_frozen / rain.ρ′ * relative_inverse_density(rain_fB, rain_f, rain_B, rain_u),
    )
end
