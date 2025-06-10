
"""
    ∫fdD(f, dist; [p = 1e-6], kwargs...)

Integrate the function `f` over the size distribution `dist`

!!! note "Usage"
    This function is useful for integrating functions over the size distribution `dist`.
    It is a light wrapper around `QGK.quadgk` that automatically inserts appropriate
    size distribution thresholds as integration limits.

This method calls [`∫fdD_error`](@ref), which returns both the value of the integral 
    and the estimated error. Since the error is typically not of interest, this method
    only returns the value of the integral.

# Arguments
- `f`: The function to integrate
- `dist`: The distribution object, passed to [`integral_bounds`](@ref) to set the integration limits.
- `p`: The integration bounds are set to the `p`-th and `1-p`-th quantiles of the size distribution.
    Default: `p = 1e-6` (i.e. 99.9998% of the size distribution is integrated).
- `kwargs`: Additional optional keyword arguments to pass to [`QGK.quadgk`](https://juliamath.github.io/QuadGK.jl/stable/api/#QuadGK.quadgk)
    - `rtol`: The relative tolerance for the integration, default: `rtol = sqrt(eps(FT))`
    - `atol`: The absolute tolerance for the integration, default: `atol = 0`
    - `maxevals`: The maximum number of function evaluations, default: `maxevals = 10^7`
    - `order`: The order of the quadrature rule, default: `order = 7`

# Returns
- `value`: The value of the integral

!!! note "Integral accuracy"
    To achieve highest accuracy, which can be challenging when integrating the
    [`N′ice`](@ref) function, it is recommended to increase the `order` of the 
    quadrature rule and set `rtol = 0`. Experimentally, `order = 44` has been found
    to be sufficient for most cases. 

    For convenience, passing `accurate = true` will set `rtol = 0` and `order = 44`.

# Examples

```jldoctest
julia> import CloudMicrophysics.Parameters as CMP,
              CloudMicrophysics.P3Scheme   as P3

julia> params = CMP.ParametersP3(Float64);

julia> state = P3.get_state(params; F_rim = 0.0, ρ_rim = 400.0, L_ice = 0.002, N_ice = 1000.0);

julia> logλ = P3.get_distribution_logλ(state);

julia> f(D) = D^3 * P3.N′ice(state, logλ)(D);  # Define a function to integrate

julia> P3.∫fdD(f, state, logλ; p = 0.01)  # Integrate the function
0.0008519464332296608

julia> P3.∫fdD(state, logλ; p = 0.01) do D  # Integrate with a `do`-block
           P3.ice_mass(state, D) * P3.N′ice(state, logλ)(D)
       end
0.0017027833723511712
```
"""
function ∫fdD(f, state::P3State, logλ; p = 1e-6, moment_order = 0, kwargs...)
    return ∫fdD_error(f, state, logλ; p, moment_order, kwargs...)[1]
end

"""
    ∫fdD_error(f, dist; p, [moment_order = 0], [accurate = false], kwargs...)

Integrate the function `f` over the size distribution `dist`

# Returns
- `value`: The value of the integral
- `error`: The estimated error of the integral

# Notes
See [`∫fdD`](@ref), which only returns the value of the integral and not the error, for details.
"""
function ∫fdD_error(f, state::P3State, logλ; p, moment_order = 0, accurate = false, kwargs...)
    # Get integration bounds
    bnds = integral_bounds(state, logλ; p, moment_order)
    # Use a more accurate quadrature rule if requested
    accurate && (kwargs = (; rtol = 0, order = 44, kwargs...))
    return QGK.quadgk(f, bnds...; kwargs...)
end

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
- `bnds`: The integration bounds (a `Tuple`), for use in [`QGK.quadgk`].
"""
function integral_bounds(state::P3State{FT}, logλ; p, moment_order = 0) where {FT}
    # Get mass thresholds
    (; D_th, D_gr, D_cr) = get_thresholds_ρ_g(state)
    # Get bounds from quantiles
    k = get_μ(state, logλ) + moment_order
    D_min = DT.generalized_gamma_quantile(k, FT(1), exp(logλ), FT(p))
    D_max = DT.generalized_gamma_quantile(k, FT(1), exp(logλ), FT(1 - p))

    # Only integrate up to the maximum diameter, `D_max`, including intermediate thresholds
    # If `F_rim` is very close to 1, `D_cr` may be greater than `D_max`, in which case it is disregarded.
    bnds = (D_min, filter(<(D_max), (D_th, D_gr, D_cr))..., D_max)
    return bnds
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
    log_N₀ = get_logN₀(state.N_ice, μ, logλ)
    return exp(log_N₀ + mass_weighted_moment) / state.L_ice
end
