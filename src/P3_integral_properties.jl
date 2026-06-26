
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
function integral_bounds(state::P3State{FT}, logλ; p, moment_order = 0) where {FT}
    # Get reduced lower and upper bounds from quantiles
    k = get_μ(state, logλ) + moment_order
    D_min = DT.generalized_gamma_quantile(k, FT(1), exp(logλ), FT(p))
    D_max = DT.generalized_gamma_quantile(k, FT(1), exp(logλ), FT(1 - p))

    # Only integrate up to the maximum diameter, `D_max`, including intermediate thresholds
    # If `F_rim` is very close to 1, `D_cr` may be greater than `D_max`, in which case it is disregarded.
    return segment_boundaries(state, D_min, D_max)
end

"""
    D_m(state::P3State, logλ)

Compute the mass weighted mean particle size [m]

# Parameters
 - `state`: [`P3State`](@ref) object
 - `logλ`: The log of the slope parameter [log(1/m)]
"""
function D_m(state, logλ)
    μ = get_μ(state, logλ)
    mass_weighted_moment = logmass_gamma_moment(state, μ, logλ; n = 1)
    log_N₀ = get_logN₀(state.ρn_ice, μ, logλ)
    return exp(log_N₀ + mass_weighted_moment) / state.ρq_ice
end
