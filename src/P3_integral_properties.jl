
# Quadrature rules (`QuadratureRule`, `ChebyshevGauss`, `GaussLegendre`,
# `integrate`, `subintervals`, and the `node`/`weight`/`inv_weight_fun`
# accessors) live in the `Quadrature` module. It is included before
# `Parameters` so a constructed rule can be stored on a parameter struct and
# shipped to GPU kernels rather than rebuilt in-kernel. The names are imported
# and re-exported in `P3.jl`, so the `integrate(...; quad = ...)` API used
# throughout the P3 scheme is unchanged.

"""
    integral_bounds(state::P3State, logλ; p, moment_order = 0)

Compute the integration bounds for the P3 size distribution,

    Mⁿ = ∫_a^b Dⁿ * N′(D) dD = N₀ * ∫_a^b Dⁿ D^μ * exp(-λ * D) dD

 where `Mⁿ` is the `n`-th moment of the size distribution.
 Here `n ≡ moment_order` and `a` and `b` are the integration bounds.
 For a proper moment, `a=0` and `b=∞`. For the numerical integration, `a` and `b`
 are determined by this function.

# Arguments
- `state`: [`P3State`](@ref) object
- `logλ`: The log of the slope parameter [log(1/m)]
- `p`: The integration bounds are set to the `p`-th and `1-p`-th quantiles of the size distribution.
- `moment_order`: For integrands proportional to moments of the size distribution,
    `moment_order` can be used to indicate the order of the moment.
    May provide more accurate bounds; thus more accurate integration.
    Default: `moment_order = 0`.

# Returns
- `bnds`: The integration bounds (a `Tuple`), for use in numerical integration (c.f. [`integrate`](@ref)).
"""
@inline function integral_bounds(state::P3State{FT}, logλ; p, moment_order = 0) where {FT}
    # Get reduced lower and upper bounds from quantiles
    k = get_μ(state, logλ) + moment_order
    λ = exp(logλ)
    # μ == 1 here, so use the unit-μ quantile (avoids a `(z/λ)^1` runtime pow per bound)
    D_min = DT.generalized_gamma_quantile_unit_μ(k, λ, FT(p))
    D_max = DT.generalized_gamma_quantile_unit_μ(k, λ, FT(1 - p))

    # Only integrate up to the maximum diameter, `D_max`, including intermediate thresholds
    # If `F_rim` is very close to 1, `D_cr` may be greater than `D_max`, in which case it is disregarded.
    bnds = segment_boundaries(state, D_min, D_max)
    # `D_max` sits `log(1/p)` decay lengths into the size-distribution tail; a
    # breakpoint at the decay scale keeps each subinterval resolvable at low order
    D_e = clamp(3 / λ, D_min, D_max)
    return Tuple(SA.sort(SA.SVector(bnds..., D_e)))
end

"""
    velocity_integral_bounds(state::P3State, logλ, v_term; p, moment_order = 0)

Compute the integration bounds for a velocity-weighted P3 integral: the
mass-regime [`integral_bounds`](@ref) with the [`velocity_breakpoints`](@ref)
of the terminal-velocity closure `v_term` clamped into `[D_min, D_max]` and
re-sorted, so each breakpoint coincides with a subinterval boundary. Returns a
fixed-length tuple.
"""
function velocity_integral_bounds(state::P3State{FT}, logλ, v_term::V; p, moment_order = 0) where {FT, V}
    bnds = integral_bounds(state, logλ; p, moment_order)
    breaks = map(D -> clamp(FT(D), first(bnds), last(bnds)), velocity_breakpoints(v_term))
    return Tuple(SA.sort(SA.SVector(bnds..., breaks...)))
end

"""
    D_m(state, logλ)

Compute the mass weighted mean particle size [m]

# Arguments
 - `state`: [`P3State`](@ref) object
 - `logλ`: The log of the slope parameter [log(1/m)]
"""
function D_m(state, logλ)
    μ = get_μ(state, logλ)
    mass_weighted_moment = logmass_gamma_moment(state, μ, logλ; n = 1)
    # Same presence gate as `logN′ice` (`P3_size_distribution.jl`): `get_logN₀`
    # requires a present `ρn_ice`, and a floored-but-finite surrogate would
    # manufacture a nonzero amplitude from an absent population rather than
    # leaving it absent.
    (; ρn_ice, ρq_ice) = state
    present = FD.value(ρn_ice) > zero(FD.value(ρn_ice))
    ρn_ice_safe = ifelse(present, ρn_ice, one(ρn_ice))
    log_N₀_safe = get_logN₀(ρn_ice_safe, μ, logλ)
    log_N₀ = ifelse(present, log_N₀_safe, oftype(log_N₀_safe, -Inf))
    return exp(log_N₀ + mass_weighted_moment) / ρq_ice
end
