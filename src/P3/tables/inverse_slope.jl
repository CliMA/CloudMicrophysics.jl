"""
    SlopePoint{P}

The point function of the inverse slope table, the log mean particle mass [`logLdivN`](@ref) at a point `(logλ, w, ρ_rim)`.

The table is the forward map of the slope solve, `log(L/N)` as a function of `logλ`, the rime
coordinate `w` (see [`rime_coordinate`](@ref)) and the rime density `ρ_rim`.
[`tabulated_distribution_logλ`](@ref) inverts it.

# Construction

    SlopePoint(params)

with `params` a [`CMP.ParametersP3`](@ref).

# Fields
$(FIELDS)
"""
struct SlopePoint{P}
    "The P3 parameters"
    params::P
end
(f::SlopePoint)((logλ, w, ρ_rim)) = logLdivN(P3Particles(f.params, rime_mass_fraction_at(w), ρ_rim), logλ)

group_table(::Val{:slope}, values, grid, ::Type{FT}) where {FT} = UniformTable(values, table_axes(grid).slope, FT)

"""
    tabulated_distribution_logλ(table, particles, ρq_ice, ρn_ice)

Return `logλ`, the log of the slope parameter of the P3 size distribution, from the inverse slope table.

The result solves the same equation as [`get_distribution_logλ`](@ref), `logLdivN(particles, logλ)
= log(ρq_ice / ρn_ice)`, in three steps:

1. The table `logLdivN` is interpolated bilinearly in the rime coordinate and the rime density of
   `particles`, which gives a decreasing function of `logλ` on the table points.
2. A bisection with a fixed number of steps finds the segment between two `logλ` points that
   holds the target, and the target is inverted linearly on the segment.
3. One Newton step on the closed form [`logLdivN`](@ref) with the slope of the segment refines
   the result, which is then kept within the bounds of the table.

As for [`get_distribution_logλ`](@ref), the result is the upper bound without ice mass, the lower
bound with mass but without number, and a primal value, constant under differentiation with
respect to the state.

# Arguments
- `table`: the `slope` table of a [`P3Tables`](@ref), the values of [`SlopePoint`](@ref) on `(logλ, w, ρ_rim)`
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

"""
    bilinear_column(data, k, (j, a_j), (l, a_l))

Return `data[1, k, :, :]` interpolated bilinearly in the cell `(j, l)` with the fractions `(a_j, a_l)`.
"""
@inline function bilinear_column(data, k, (j, a_j), (l, a_l))
    @inbounds lower = (1 - a_j) * data[1, k, j, l] + a_j * data[1, k, j + 1, l]
    @inbounds upper = (1 - a_j) * data[1, k, j, l + 1] + a_j * data[1, k, j + 1, l + 1]
    return (1 - a_l) * lower + a_l * upper
end

"""
    bisection_steps(n)

Return the number of bisection steps that reduce `n` points to one segment, `⌈log₂(n - 1)⌉`.
"""
@inline bisection_steps(n) = 8 * sizeof(Int) - leading_zeros(Int(n) - 2)

"""
    P3Distribution(particles, tables::P3Tables, ρq_ice, ρn_ice)

Return the size distribution with the slope `logλ` from the inverse slope table of `tables`.

The number is bounded as for a `QuadratureRule`, and `logλ` follows from
[`tabulated_distribution_logλ`](@ref) instead of the solve.
"""
function P3Distribution(particles::P3Particles, tables::P3Tables, ρq_ice, ρn_ice)
    ρn = bounded_ice_number(particles.params, ρq_ice, ρn_ice)
    logλ = tabulated_distribution_logλ(tables.groups.slope, particles, ρq_ice, ρn)
    return P3Distribution(particles.params.slope, ρn, logλ)
end
