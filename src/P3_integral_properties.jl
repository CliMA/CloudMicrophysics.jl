
# Quadrature rules (`QuadratureRule`, `ChebyshevGauss`, `GaussLegendre`,
# `integrate`, `subintervals`, and the `node`/`weight`/`inv_weight_fun`
# accessors) live in the `Quadrature` module. It is included before
# `Parameters` so a constructed rule can be stored on a parameter struct and
# shipped to GPU kernels rather than rebuilt in-kernel. The names are imported
# and re-exported in `P3.jl`.

"""
    integral_bounds(particles, dist; p, moment_order = 0)

Compute the integration bounds for the moments of the P3 size distribution.

The moments are

    Mⁿ = ∫_a^b Dⁿ * N′(D) dD = N₀ * ∫_a^b Dⁿ D^μ * exp(-λ * D) dD

 where `Mⁿ` is the `n`-th moment of the size distribution.
 Here `n ≡ moment_order` and `a` and `b` are the integration bounds.
 For a proper moment, `a=0` and `b=∞`. For the numerical integration, `a` and `b`
 are determined by this function.

# Arguments
- `particles`: [`P3Particles`](@ref)
- `dist`: [`P3Distribution`](@ref)
- `p`: The integration bounds are set to the `p`-th and `1-p`-th quantiles of the size distribution.
- `moment_order`: For integrands proportional to moments of the size distribution,
    `moment_order` can be used to indicate the order of the moment.
    May provide more accurate bounds; thus more accurate integration.
    Default: `moment_order = 0`.

# Returns
- `bnds`: The integration bounds (a `Tuple`), for use in numerical integration (c.f. [`integrate`](@ref)).
"""
@inline function integral_bounds(particles, (; μ, logλ); p, moment_order = 0)
    FT = eltype(particles)
    # Get reduced lower and upper bounds from quantiles
    k = μ + moment_order
    λ = exp(logλ)
    # μ == 1 here, so use the unit-μ quantile (avoids a `(z/λ)^1` runtime pow per bound)
    D_min = DT.generalized_gamma_quantile_unit_μ(k, λ, FT(p))
    D_max = DT.generalized_gamma_quantile_unit_μ(k, λ, FT(1 - p))

    # Only integrate up to the maximum diameter, `D_max`, including intermediate thresholds
    # If `F_rim` is very close to 1, `D_cr` may be greater than `D_max`, in which case it is disregarded.
    bnds = segment_boundaries(particles, D_min, D_max)
    # A breakpoint at the decay scale of the size distribution keeps each subinterval resolvable
    # at low order
    D_e = clamp(3 / λ, D_min, D_max)
    return insert_breakpoints(bnds, (D_e,))
end

"""
    insert_breakpoints(bnds, Ds)

Return the integration bounds `bnds` with the diameters `Ds` inserted, in ascending order.

Each diameter in `Ds` must lie in `[first(bnds), last(bnds)]`.
"""
@inline insert_breakpoints(bnds, Ds) = Tuple(SA.sort(SA.SVector(bnds..., Ds...)))

"""
    velocity_integral_bounds(particles, dist, v_term; p, moment_order = 0)

Compute the integration bounds for a velocity-weighted P3 integral.

This function returns the same integration bounds as [`integral_bounds`](@ref),
  but with additional breakpoints that a terminal-velocity closure `v_term` may introduce,
  see [`velocity_breakpoints`](@ref).

# Arguments
- `particles`: [`P3Particles`](@ref)
- `dist`: [`P3Distribution`](@ref)
- `v_term`: terminal-velocity closure, see [`ice_particle_terminal_velocity`](@ref)

See [`integral_bounds`](@ref) for the keyword arguments and the return value.
"""
function velocity_integral_bounds(
    particles, dist, v_term::V; p, moment_order = 0,
) where {V}
    bnds = integral_bounds(particles, dist; p, moment_order)
    D_min, D_max = first(bnds), last(bnds)
    breaks = map(D -> clamp(oftype(D_min, D), D_min, D_max), velocity_breakpoints(v_term))
    return insert_breakpoints(bnds, breaks)
end

"""
    D_m(particles, dist)

Compute the mass-weighted mean particle size [m] of a P3 size distribution.

The size is `∫ D m(D) N′(D) dD / ∫ m(D) N′(D) dD` for the [`P3Distribution`](@ref) `dist` and the
[`P3Particles`](@ref) `particles`, or zero without ice number.

# Arguments
- `particles`: [`P3Particles`](@ref)
- `dist`: [`P3Distribution`](@ref)
"""
function D_m(particles, (; ρn, μ, logλ))
    log_∫DmN′ = logmass_gamma_moment(particles, μ, logλ; n = 1)  # log(∫ D m(D) N′(D) dD / N₀)
    log_∫mN′ = logmass_gamma_moment(particles, μ, logλ; n = 0)   # log(∫ m(D) N′(D) dD / N₀)
    D_m = exp(log_∫DmN′ - log_∫mN′)
    # Compare the primal value, since a dual number with zero value and nonzero partials exceeds zero
    return ifelse(FD.value(ρn) > 0, D_m, zero(D_m))
end
