"""
    ∫fdD(f, state::P3State; [D_max = 1], kwargs...)

Integrate the function `f` over the size distribution of the ice particles.

!!! note "Usage"
    This function is useful for integrating functions over the size distribution of the ice particles.
    It is a light wrapper around `QGK.quadgk` that automatically inserts the
    calculated size distribution thresholds (i.e. `D_th`, `D_gr`, `D_cr`) as integration limits.
    The upper integration limit is by default set to `D_max = 1 m`, with the 
    presumption that all integrals will decay to zero much quicker than this.
    Formally, integration is to infinity.

This method calls [`∫fdD_error`](@ref), which returns both the value of the integral 
    and the estimated error. Since the error is typically not of interest, this method
    only returns the value of the integral.

# Arguments
- `f`: The function to integrate
- `state`: The [`P3State`](@ref) object
- `D_max`: The maximum diameter to integrate to [m]
- `kwargs`: Additional keyword arguments to pass to [`QGK.quadgk`](https://juliamath.github.io/QuadGK.jl/stable/api/#QuadGK.quadgk)
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

```julia-repl
julia> import CloudMicrophysics.Parameters as CMP
julia> import CloudMicrophysics.P3Scheme as P3

# Get a state object
julia> params = CMP.ParametersP3(Float64)
julia> state = P3.get_state(params; F_rim = 0.0, ρ_r = 400.0)

# Define a function to integrate
julia> f(D) = D^3

# Integrate the function
julia> P3.∫fdD(f, state)
0.25

# Integrate with a `do`-block
julia> P3.∫fdD(state) do D
           P3.ice_mass(state, D)
       end
0.006392317884295288
```
"""
function ∫fdD(f, state::P3State; D_max = 1, kwargs...)
    return ∫fdD_error(f, state; D_max, kwargs...)[1]
end

"""
    ∫fdD_error(f, state::P3State; D_max = 1, kwargs...)

Integrate the function `f` over the size distribution of the ice particles.

# Returns
- `value`: The value of the integral
- `error`: The estimated error of the integral

# Notes
See [`∫fdD`](@ref), which only returns the value of the integral and not the error, for details.
"""
function ∫fdD_error(f, state::P3State; D_max = 1, accurate = false, kwargs...)
    thresholds = threshold_tuple(state)
    @assert thresholds[1] > 0
    # Only integrate up to the maximum diameter, `D_max`, including intermediate thresholds
    # If `F_rim` is very close to 1, `D_cr` may be greater than `D_max`, in which case it is disregarded.
    bnds = (0, filter(<(D_max), thresholds)..., D_max)
    # Use a more accurate quadrature rule if requested
    accurate && (kwargs = (; rtol = 0, order = 44, kwargs...))
    return QGK.quadgk(f, bnds...; kwargs...)
end


"""
    D_m(dist::P3Distribution)

Compute the mass weighted mean particle size [m]

# Parameters
 - `dist`: [`P3Distribution`](@ref) object
"""
function D_m(dist::P3Distribution{FT}) where {FT}
    (; log_λ, log_N₀, L) = dist
    L < eps(FT) && return FT(0)

    state = get_state(dist)
    (; F_rim, ρ_g, D_th, D_gr, D_cr) = state
    (; ρ_i, ρ_l, mass) = get_parameters(state)
    (; α_va, β_va) = mass

    # Calculate the mass weighted mean particle size
    μ = get_μ(state, log_λ)
    ∞ = FT(Inf)
    G = log_integrate_moment_psd
    # G_liqfrac = G(0, ∞, ρ_l * π / 6, 4, μ, log_λ)
    G_small_spherical = G(0, D_th, ρ_i * π / 6, 4, μ, log_λ)
    G_summed = if isunrimed(state)
        G_large_unrimed = G(D_th, ∞, α_va, β_va + 1, μ, log_λ)
        logsumexp((G_small_spherical, G_large_unrimed))
    else
        G_dense_nonspherical = G(D_th, D_gr, α_va, β_va + 1, μ, log_λ)
        G_graupel = G(D_gr, D_cr, ρ_g * π / 6, 4, μ, log_λ)
        G_partially_rimed =
            G(D_cr, ∞, α_va / (1 - F_rim), β_va + 1, μ, log_λ)
        logsumexp((
            G_small_spherical,
            G_dense_nonspherical,
            G_graupel,
            G_partially_rimed,
        ))
    end

    # compute `F_liq`-weighted average and normalize by `L`
    return exp(G_summed) * exp(log_N₀) / L
    # TODO: Implement `F_liq`-weighted average
    # return weighted_average(state.F_liq, exp(G_liqfrac), exp(G_summed)) * exp(log_N₀) / L
end
