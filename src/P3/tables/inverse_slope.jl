"""
    slope_point(params, (logλ, w, ρ_rim))

Return the log mean particle mass [`logLdivN`](@ref) at the slope `logλ`, the rime coordinate `w` and the rime density `ρ_rim`, the value of the inverse slope table.

[`tabulated_distribution_logλ`](@ref) inverts the table along `logλ`.
"""
slope_point(params, (logλ, w, ρ_rim)) = logLdivN(P3Particles(params, rime_mass_fraction_at(w), ρ_rim), logλ)

"""
    tabulated_distribution_logλ(table, particles, ρq_ice, ρn_ice)

Return `logλ`, the log of the slope parameter of the P3 size distribution, from the inverse slope table.

The result solves the equation of [`get_distribution_logλ`](@ref), with the same values without
ice mass or number, and is a primal value.

# Arguments
- `table`: the `slope` table of a [`P3Tables`](@ref), the values of [`slope_point`](@ref) on `(logλ, w, ρ_rim)`
- `particles`: [`P3Particles`](@ref)
- `ρq_ice`: ice mass concentration [kg/m³]
- `ρn_ice`: ice number concentration [1/m³], bounded by [`bounded_ice_number`](@ref)
"""
@inline function tabulated_distribution_logλ(table, particles, ρq_ice, ρn_ice)
    particles = _primal_particles(particles)
    FT = eltype(particles)
    ρq_ice = FT(primal_value(ρq_ice))
    ρn_ice = FT(primal_value(ρn_ice))
    lo, hi = map(FT, logλ_bounds(particles.params))
    target = log(max(ρq_ice, floatmin(FT))) - log(max(ρn_ice, floatmin(FT)))

    (; data) = table
    nλ, nw, nρ = Base.tail(size(data))
    cell_w = uniform_index(FT(rime_coordinate(particles.F_rim)), FT(table.lo[2]), FT(table.hi[2]), nw)
    cell_ρ = uniform_index(particles.ρ_rim, FT(table.lo[3]), FT(table.hi[3]), nρ)

    # The column decreases with `logλ`, and the bisection keeps the target between points k_lo and k_hi
    k_lo, k_hi = 1, nλ
    for _ in 1:bisection_steps(nλ)
        k = (k_lo + k_hi) >>> 1
        above = bilinear_column(data, k, cell_w, cell_ρ) ≥ target
        k_lo = ifelse(above, k, k_lo)
        k_hi = ifelse(above, k_hi, k)
    end
    c_lo = bilinear_column(data, k_lo, cell_w, cell_ρ)
    c_hi = bilinear_column(data, k_lo + 1, cell_w, cell_ρ)
    Δlogλ = (FT(table.hi[1]) - FT(table.lo[1])) / (nλ - 1)
    s = clamp(UT.guarded_quotient(c_lo - target, c_lo - c_hi), zero(FT), one(FT))
    logλ = FT(table.lo[1]) + (k_lo - 1 + s) * Δlogλ

    slope = (c_hi - c_lo) / Δlogλ
    residual = logLdivN(particles, logλ) - target
    logλ = clamp(ifelse(slope < 0, logλ - residual / slope, logλ), lo, hi)

    logλ = ifelse(ρn_ice ≥ floatmin(FT), logλ, lo)
    return ifelse(ρq_ice ≥ floatmin(FT), logλ, hi)
end

@inline function bilinear_column(data, k, (j, a_j), (l, a_l))
    @inbounds lower = (1 - a_j) * data[1, k, j, l] + a_j * data[1, k, j + 1, l]
    @inbounds upper = (1 - a_j) * data[1, k, j, l + 1] + a_j * data[1, k, j + 1, l + 1]
    return (1 - a_l) * lower + a_l * upper
end

@inline bisection_steps(n) = 8 * sizeof(Int) - leading_zeros(Int(n) - 2)

"""
    P3Distribution(particles, tables::P3Tables, ρq_ice, ρn_ice)

Return the size distribution with the slope `logλ` from the inverse slope table of `tables`.

The number is bounded as for a `QuadratureRule`, and `logλ` follows from
[`tabulated_distribution_logλ`](@ref) instead of the solve.
"""
@inline function P3Distribution(particles::P3Particles, tables::P3Tables, ρq_ice, ρn_ice)
    ρn = bounded_ice_number(particles.params, ρq_ice, ρn_ice)
    logλ = tabulated_distribution_logλ(tables.groups.slope, particles, ρq_ice, ρn)
    return P3Distribution(particles.params.slope, ρn, logλ)
end
