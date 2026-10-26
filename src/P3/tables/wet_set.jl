"""
    collision_bounds_point((μ,))

Return the `SVector` of the logs of the bounds of the scaled diameter `ξ = D λ` of the collision integrals at the shape parameter `μ`, see [`collision_logξ_lo`](@ref).
"""
collision_bounds_point((μ,)) = SA.SVector(collision_logξ_lo(μ), collision_logξ_hi(μ))

"""
    collected_point(params, vel, pdf, generator, (logD, w, ρ_rim, logρₐ, logx̄))

Return the mass of one liquid species that one ice particle of diameter `exp(logD)` collects per unit liquid mass [m³/s], and the impact scale `log κ̄` of that mass.

The rates are those of [`particle_collection`](@ref), and `log κ̄` is the mean of `log κ` over the
collected mass, see [`impact_scale`](@ref). The table stores the log of the mass rate.
"""
function collected_point(params, vel, pdf, generator, (logD, w, ρ_rim, logρₐ, logx̄))
    particles = P3Particles(params, rime_mass_fraction_at(w), ρ_rim)
    liquid = table_liquid(vel, particles, pdf, exp(logρₐ), logx̄)
    r = particle_collection(params, particles, liquid, (), generator)(exp(logD))
    return SA.SVector(r[2] / liquid.mass, r[3] / r[2])
end

"""
    capacity_point(params, vel, (logD, w, ρ_rim, logρₐ))

Return `D^(3/2) v(D)^(1/2)` [m² s^(-1/2)] of one ice particle of diameter `D = exp(logD)`, the integrand of the ventilation moment of [`ventilation_moment`](@ref).
"""
function capacity_point(params, vel, (logD, w, ρ_rim, logρₐ))
    particles = P3Particles(params, rime_mass_fraction_at(w), ρ_rim)
    D = exp(logD)
    v = ice_particle_terminal_velocity(vel, particles, exp(logρₐ))(D)
    return D * sqrt(D * v)
end

"""
    fraction_point(params, vel, pdf, generator, logits, (logλ, w, ρ_rim, logρₐ, logx̄))

Return the `SVector` of the logits of the fractions of the mass of one liquid species that the ice below the diameters of the number logits `logits` collects.

The fractions are those of [`cumulative_logits`](@ref) of the mass rate of one particle of
[`particle_collection`](@ref), for the state of [`table_state`](@ref).
"""
function fraction_point(params, vel, pdf, generator, logits, (logλ, w, ρ_rim, logρₐ, logx̄))
    (; particles, dist) = table_state(params, generator, logλ, w, ρ_rim)
    liquid = table_liquid(vel, particles, pdf, exp(logρₐ), logx̄)
    rates_at = particle_collection(params, particles, liquid, (), generator)
    return cumulative_logits(D -> rates_at(D)[2], particles, dist, liquid.v_i, logits, generator)
end

"""
    capacity_fraction_point(params, vel, generator, logits, (logλ, w, ρ_rim, logρₐ))

Return the `SVector` of the logits of the fractions of the ventilation moment `M_v` below the diameters of the number logits `logits`.

The fractions are those of [`cumulative_logits`](@ref) of the integrand of
[`capacity_point`](@ref), for the state of [`table_state`](@ref).
"""
function capacity_fraction_point(params, vel, generator, logits, (logλ, w, ρ_rim, logρₐ))
    (; particles, dist) = table_state(params, generator, logλ, w, ρ_rim)
    v_i = ice_particle_terminal_velocity(vel, particles, exp(logρₐ))
    return cumulative_logits(D -> D * sqrt(D * v_i(D)), particles, dist, v_i, logits, generator)
end

"""
    cumulative_logits(k, particles, dist, v_i, logits, generator)

Return the `SVector` of the logits of the fractions of `∫ k(D) N′(D) dD` below the diameters at the number logits `logits`, see [`number_logit`](@ref).

The first and last of `logits` are those of the bounds of the liquid-ice collisions, and the
integral extends between their diameters, split at the diameters of the other logits, the
thresholds of `particles` and the breakpoints of `v_i`. Both terms of each logit
`log(F) - log(1 - F)` are sums of intervals, at least `max(1e-30 total, floatmin)`, so the first
and last logits are those of the fractions 0 and 1 at this floor.

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
    p = liquid_ice_collision_quantile(eltype(particles))
    bounds = velocity_integral_bounds(particles, dist, v_i; p)
    D_lo, D_hi = first(bounds), last(bounds)
    breaks = SA.sort(SA.SVector(particles.D_th, particles.D_gr, particles.D_cr, velocity_breakpoints(v_i)...))
    cuts = ntuple(Val(K)) do j
        (F, Q) = logistic_fractions(logits[j])
        return clamp(UT.gamma_inc_inv(dist.μ + 1, F, Q) * exp(-dist.logλ), D_lo, D_hi)
    end
    pieces = SA.SVector(ntuple(i -> interval_integral(k, n_i, cuts[i], cuts[i + 1], breaks, generator), Val(K - 1)))
    FT = eltype(pieces)
    smallest = max(FT(1e-30) * sum(pieces), floatmin(FT))
    return SA.SVector(ntuple(Val(K)) do j
        below = sum(ntuple(i -> ifelse(i < j, pieces[i], zero(eltype(pieces))), Val(K - 1)))
        above = sum(ntuple(i -> ifelse(i ≥ j, pieces[i], zero(eltype(pieces))), Val(K - 1)))
        return log(max(below, smallest)) - log(max(above, smallest))
    end)
end

"""
    interval_integral(k, n_i, a, b, breaks, generator)

Return `∫ k(D) n_i(D) dD` from `a` to `b` with the rule `generator`, split at the diameters of `breaks` between `a` and `b`, see [`cumulative_logits`](@ref).
"""
@inline interval_integral(k::FK, n_i::FN, a, b, breaks, generator) where {FK, FN} =
    integrate(D -> n_i(D) * k(D), (a, map(D -> clamp(D, a, b), Tuple(breaks))..., b), generator)

"""
    number_logit(μ, ξ)

Return the logit `log(P) - log(1 - P)` of the number fraction `P = P(μ + 1, ξ)` of the ice below the scaled diameter `ξ = D λ`, see [`UT.gamma_inc`](@ref).

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

Return the fraction `F = 1 / (1 + exp(-z))` of the logit `z` limited to `±50`, and `1 - F`.
"""
@inline function logistic_fractions(z)
    z = clamp(z, -50, 50)
    return (inv(1 + exp(-z)), inv(1 + exp(z)))
end

@inline function fraction_between(z_a, z_b)
    (F_a, Q_a) = logistic_fractions(z_a)
    (F_b, Q_b) = logistic_fractions(z_b)
    return max(ifelse(z_b < 0, F_b - F_a, Q_a - Q_b), zero(F_a))
end

@inline fixed_cells(t, first_axis, x) =
    map((d, x_d) -> axis_cell(t, first_axis + d - 1, x_d), ntuple(identity, Val(length(x))), x)

"""
    wet_set_scan_intervals

The number of intervals of the scan of [`tabulated_wet_set`](@ref).
"""
const wet_set_scan_intervals = 24

"""
    tabulated_wet_set(tables, state, key, cloud, rain, coeffs, logT)

Return the freeze and shed partition of the liquid-ice collisions per unit ice number from the wet-set tables.

The result is a `NamedTuple` with the frozen mass `cloud_frozen`, `rain_frozen`, the shed mass
`cloud_shed`, `rain_shed`, and the rime volume of the frozen mass `cloud_rime`, `rain_rime`, of
cloud and rain, all per unit ice number, see [Wet-set tables](@ref).

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
    (; groups) = tables
    μ = FT(primal_value(state.dist.μ))
    logλ = FT(primal_value(state.dist.logλ))
    (_, w, ρ_rim, logρₐ) = key
    (_, w_p, ρ_rim_p, logρₐ_p) = map(@inline(x -> oftype(μ, primal_value(x))), key)
    t_c, t_r, t_cap = groups.cloud_collected, groups.rain_collected, groups.capacity
    ρ′_local = state.particles.params.ρ_rim_local
    L_c_p, L_r_p, c₁_p, c₂_p = map(@inline(x -> oftype(μ, primal_value(x))), (cloud.mass, rain.mass, c₁, c₂))

    # The log ratio of the collected mass to the maximum freezing rate of a particle of diameter `exp(logD)`
    cells_c_p = fixed_cells(t_c, 2, (w_p, ρ_rim_p, logρₐ_p, FT(primal_value(cloud.logx̄))))
    cells_r_p = fixed_cells(t_r, 2, (w_p, ρ_rim_p, logρₐ_p, FT(primal_value(rain.logx̄))))
    cells_cap_p = fixed_cells(t_cap, 2, (w_p, ρ_rim_p, logρₐ_p))
    smallest = floatmin(FT)
    @inline function log_ratio(logD)
        cell_D = axis_cell(t_c, 1, logD)
        u =
            L_c_p * exp(interpolate_or(cloud.present, 0, t_c, (cell_D, cells_c_p...))[1]) +
            L_r_p * exp(interpolate_or(rain.present, 0, t_r, (cell_D, cells_r_p...))[1])
        m = c₁_p * exp(logD) + c₂_p * exp(interpolate(t_cap, (cell_D, cells_cap_p...))[1])
        return log(max(u, smallest)) - log(m)
    end
    # The integrands of the partition per unit log D at a particle of diameter `exp(logD)`, up to a
    # constant factor, with the local rime density at the impact scale of the mass it collects
    cells_c = fixed_cells(t_c, 2, (w, ρ_rim, logρₐ, cloud.logx̄))
    cells_r = fixed_cells(t_r, 2, (w, ρ_rim, logρₐ, rain.logx̄))
    cells_cap = fixed_cells(t_cap, 2, (w, ρ_rim, logρₐ))
    @inline function scan_point(logD)
        cell_D = axis_cell(t_c, 1, logD)
        (log_u_c, log_κ_c) = interpolate_or(cloud.present, 0, t_c, (cell_D, cells_c...))
        (log_u_r, log_κ_r) = interpolate_or(rain.present, 0, t_r, (cell_D, cells_r...))
        u_c = cloud.mass * exp(log_u_c)
        u_r = rain.mass * exp(log_u_r)
        m = c₁ * exp(logD) + c₂ * exp(interpolate(t_cap, (cell_D, cells_cap...))[1])
        f = min(UT.guarded_quotient(m, u_c + u_r), one(m))
        ρ′_c = ρ′_local(exp(log_κ_c - logT))
        ρ′_r = ρ′_local(exp(log_κ_r - logT))
        # D N′(D), up to a constant factor, since dD = D d(log D)
        g = exp((μ + 1) * (logD + logλ) - exp(logD + logλ))
        log_ratio_at = primal_value(log(max(u_c + u_r, smallest)) - log(m))
        return (
            log_ratio_at,
            SA.SVector(
                g * u_c * f, g * u_c * (1 - f), g * u_r * f, g * u_r * (1 - f),
                g * u_c * f / ρ′_c, g * u_c / ρ′_c, g * u_c, g * u_r * f / ρ′_r, g * u_r / ρ′_r, g * u_r,
            ),
        )
    end
    (logξ_lo, logξ_hi) = lookup(groups.collision_bounds, μ)
    logD_lo, logD_hi = logξ_lo - logλ, logξ_hi - logλ
    # The scan brackets the zeros of the log ratio and sums the integrands by the trapezoid rule
    n_scan = wet_set_scan_intervals
    @inline trapezoid(i, Δx) = ifelse((i == 0) | (i == n_scan), Δx / 2, Δx)
    scan = UT.scan_zeros(
        scan_point, log_ratio, RS.RegulaFalsiMethod, logD_lo, logD_hi, Val(n_scan), Val(2), FixedIterations{FT}(),
        3,
        trapezoid,
    )
    (; roots, n_found) = scan
    log_ratio_lo = scan.f_lo
    sums = scan.total
    # The three intervals between the bounds and the zeros are wet and dry in turn, and an
    # interval past the last zero is empty
    wet_lo = log_ratio_lo > 0
    wet = (wet_lo, !wet_lo, wet_lo)
    found = (n_found ≥ 1, n_found ≥ 2)

    # The integrals over W and the dry set from the fraction tables, with the logits -Inf and Inf at the bounds
    ξ = map(@inline(r -> exp(r + logλ)), roots)
    z = map(@inline((f, ξ) -> ifelse(f, number_logit(μ, ξ), oftype(μ, Inf))), found, ξ)
    @inline function fractions(t, x, present)
        F = map(@inline((f, z) -> (v = lookup_or(present, 0, t, z, x...)[1]; ifelse(f, v, oftype(v, Inf)))), found, z)
        zs = (-oftype(F[1], Inf), F..., oftype(F[1], Inf))
        @inline on(set) =
            sum(ntuple(@inline(i -> ifelse(set[i], fraction_between(zs[i], zs[i + 1]), zero(F[1]))), Val(3)))
        return (; wet = on(wet), dry = on(map(!, wet)))
    end
    cloud_fraction = fractions(groups.cloud_fraction, (key..., cloud.logx̄), cloud.present)
    rain_fraction = fractions(groups.rain_fraction, (key..., rain.logx̄), rain.present)
    M_v = fractions(groups.capacity_fraction, key, true).wet * exp(lookup(groups.ice, key...)[1])
    # The first moment in closed form, (μ + 1) / λ times the regularized incomplete gamma function P(μ + 2, ξ)
    PQ = map(@inline(ξ -> UT.gamma_inc(μ + 2, ξ)), ξ)
    P = (zero(FT), ifelse(found[1], first(PQ[1]), one(FT)), ifelse(found[2], first(PQ[2]), one(FT)), one(FT))
    Q = (one(FT), ifelse(found[1], last(PQ[1]), zero(FT)), ifelse(found[2], last(PQ[2]), zero(FT)), zero(FT))
    @inline between(i) = max(ifelse(P[i + 1] < one(μ) / 2, P[i + 1] - P[i], Q[i] - Q[i + 1]), zero(μ))
    M_1 = sum(ntuple(@inline(i -> ifelse(wet[i], between(i), zero(μ))), Val(3))) * (μ + 1) * exp(-logλ)
    frozen_on_W = c₁ * M_1 + c₂ * M_v
    frozen = cloud_fraction.dry * cloud.M + rain_fraction.dry * rain.M + frozen_on_W
    shed = max(cloud_fraction.wet * cloud.M + rain_fraction.wet * rain.M - frozen_on_W, zero(frozen_on_W))

    (cloud_f, cloud_s, rain_f, rain_s, cloud_fB, cloud_B, cloud_u, rain_fB, rain_B, rain_u) = sums
    # The smaller of the frozen and the shed mass, divided between cloud and rain
    freezes_less = frozen < shed
    smaller = min(frozen, shed, cloud.M + rain.M)
    # A wet set between two scan points sheds in proportion to the mass collected on it
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
    @inline relative_inverse_density(fB, f, B, u) =
        UT.guarded_quotient(UT.guarded_quotient(fB, f), UT.guarded_quotient(B, u))
    return (;
        cloud_frozen,
        rain_frozen,
        cloud_shed = cloud.M - cloud_frozen,
        rain_shed = rain.M - rain_frozen,
        cloud_rime = cloud_frozen / cloud.ρ′ * relative_inverse_density(cloud_fB, cloud_f, cloud_B, cloud_u),
        rain_rime = rain_frozen / rain.ρ′ * relative_inverse_density(rain_fB, rain_f, rain_B, rain_u),
    )
end
