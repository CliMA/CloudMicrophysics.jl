
# Callable returned by `ice_particle_terminal_velocity`: piecewise small/large-ice
# Chen 2022 velocity with optional aspect-ratio correction.
struct P3IceParticleVelocityFunctor{FT, VS, VL, S} <: Function
    v_term_small::VS
    v_term_large::VL
    D_cutoff::FT
    state::S
    use_aspect_ratio::Bool
end
@inline function (f::P3IceParticleVelocityFunctor)(D)
    vₜ = ifelse(D <= f.D_cutoff, f.v_term_small(D), f.v_term_large(D))
    return ifelse(f.use_aspect_ratio, cbrt(ϕᵢ(f.state, D)) * vₜ, vₜ)
end

"""
    ice_particle_terminal_velocity(velocity_params, ρₐ, state::P3State; [use_aspect_ratio])

Returns a single-argument function `v_term(D)` that gives the Chen 2022
terminal velocity of an ice particle of maximum dimension `D`.

# Arguments
 - `velocity_params`: A [`CMP.Chen2022VelType`](@ref)
 - `ρₐ`: Air density [kg/m³]
 - `state`: A [`P3State`](@ref)

# Keyword arguments
 - `use_aspect_ratio`: include the aspect-ratio correction (default `true`)
"""
@inline function ice_particle_terminal_velocity(
    velocity_params::CMP.Chen2022VelType, ρₐ, state::P3State; use_aspect_ratio = true,
)
    FT = typeof(ρₐ)
    (; small_ice, large_ice) = velocity_params
    D_cutoff = small_ice.cutoff
    ρᵢ = FT(916.7)  # TODO: Use parameter
    v_term_small = CO.particle_terminal_velocity(small_ice, ρₐ, ρᵢ)
    v_term_large = CO.particle_terminal_velocity(large_ice, ρₐ, ρᵢ)
    return P3IceParticleVelocityFunctor(v_term_small, v_term_large, D_cutoff, state, use_aspect_ratio)
end

struct P3NumberWeightedIntegrand{N, V} <: Function
    n::N
    v_term::V
end
@inline (f::P3NumberWeightedIntegrand)(D) = f.n(D) * f.v_term(D)

"""
    ice_terminal_velocity_number_weighted(
        velocity_params::CMP.Chen2022VelType, ρₐ, state::P3State, logλ;
        [use_aspect_ratio], [p], [quad],
    )

Return the terminal velocity of the number-weighted mean ice particle size.

# Arguments
- `velocity_params`: A [`CMP.Chen2022VelType`](@ref) with terminal velocity parameters
- `ρₐ`: Air density [kg/m³]
- `state`: A [`P3State`](@ref)
- `logλ`: The log of the slope parameter [log(1/m)]

# Keyword arguments
 - `use_aspect_ratio`: Bool flag set to `true` if we want to consider the effects
    of particle aspect ratio on its terminal velocity (default: `true`)
 - `p`: Tolerance parameter for the integral bounds. Default is 1e-6.
 - `quad`: Quadrature rule, default is `ChebyshevGauss(100)`

See also [`ice_terminal_velocity_mass_weighted`](@ref)
"""
function ice_terminal_velocity_number_weighted(
    velocity_params::CMP.Chen2022VelType, ρₐ, state::P3State, logλ;
    use_aspect_ratio = true, p = 1e-6, quad = ChebyshevGauss(100),
)
    (; ρn_ice, ρq_ice) = state
    # TODO - do we want to swicth to ϵ_numerics(FT)
    if ρn_ice < eps(one(ρn_ice)) || ρq_ice < eps(one(ρq_ice))
        return zero(ρn_ice)
    end

    v_term = ice_particle_terminal_velocity(velocity_params, ρₐ, state; use_aspect_ratio)
    n = DT.size_distribution(state, logλ)

    # ∫n(D) v(D) dD
    number_weighted_integrand = P3NumberWeightedIntegrand(n, v_term)

    bnds = integral_bounds(state, logλ; p)
    return integrate(number_weighted_integrand, bnds, quad) / ρn_ice
end

struct P3MassWeightedIntegrand{N, V, S} <: Function
    n::N
    v_term::V
    state::S
end
@inline (f::P3MassWeightedIntegrand)(D) = f.n(D) * f.v_term(D) * ice_mass(f.state, D)

"""
    ice_terminal_velocity_mass_weighted(velocity_params::CMP.Chen2022VelType, ρₐ, state::P3State, logλ; [use_aspect_ratio], [∫kwargs...])

Return the terminal velocity of the mass-weighted mean ice particle size.

# Arguments
- `velocity_params`: A [`CMP.Chen2022VelType`](@ref) with terminal velocity parameters
- `ρₐ`: Air density [kg/m³]
- `state`: A [`P3State`](@ref)
- `logλ`: The log of the slope parameter [log(1/m)]

# Keyword arguments
 - `use_aspect_ratio`: Bool flag set to `true` if we want to consider the effects
    of particle aspect ratio on its terminal velocity (default: `true`)
 - `p`: Tolerance parameter for the integral bounds. Default is 1e-6.
 - `quad`: Quadrature rule, default is `ChebyshevGauss(100)`

See also [`ice_terminal_velocity_number_weighted`](@ref)
"""
function ice_terminal_velocity_mass_weighted(
    velocity_params::CMP.Chen2022VelType, ρₐ, state::P3State, logλ;
    use_aspect_ratio = true, p = 1e-6, quad = ChebyshevGauss(100),
)
    (; ρn_ice, ρq_ice) = state
    # TODO - do we want to swicth to ϵ_numerics(FT)
    if ρn_ice < eps(one(ρn_ice)) || ρq_ice < eps(one(ρq_ice))
        return zero(ρq_ice)
    end

    v_term = ice_particle_terminal_velocity(velocity_params, ρₐ, state; use_aspect_ratio)
    n = DT.size_distribution(state, logλ)  # Number concentration at diameter D

    # ∫n(D) m(D) v(D) dD
    mass_weighted_integrand = P3MassWeightedIntegrand(n, v_term, state)

    bnds = integral_bounds(state, logλ; p)
    return integrate(mass_weighted_integrand, bnds, quad) / ρq_ice
end

"""
    ice_terminal_velocity_number_weighted_from_prognostic(
        velocity_params, ρₐ, params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, logλ; kw...
    )

Pointwise wrapper that takes the *raw prognostic* P3 ice state
(`ρq_ice`, `ρn_ice`, `ρq_rim`, `ρb_rim`) and returns the number-weighted
mean ice terminal velocity. Builds the per-cell `P3State` via
[`state_from_prognostic`](@ref), so `F_rim` is regularised to
`[0, 1 - eps(FT)]` and `ρ_rim` is clamped to `[0, 0.8 ρ_l]`.

Designed for `@.`-broadcast use from a host (CA, KiD, etc.) where the
state must be reconstructed from prognostic variables every cell.
"""
@inline function ice_terminal_velocity_number_weighted_from_prognostic(
    velocity_params, ρₐ, params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, logλ; kw...,
)
    state = state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    return ice_terminal_velocity_number_weighted(velocity_params, ρₐ, state, logλ; kw...)
end

"""
    ice_terminal_velocity_mass_weighted_from_prognostic(
        velocity_params, ρₐ, params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, logλ; kw...
    )

Mass-weighted counterpart to
[`ice_terminal_velocity_number_weighted_from_prognostic`](@ref). Builds
the per-cell `P3State` via the regularised
[`state_from_prognostic`](@ref).
"""
@inline function ice_terminal_velocity_mass_weighted_from_prognostic(
    velocity_params, ρₐ, params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, logλ; kw...,
)
    state = state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    return ice_terminal_velocity_mass_weighted(velocity_params, ρₐ, state, logλ; kw...)
end
