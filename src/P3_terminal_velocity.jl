
"""
    ice_particle_terminal_velocity(velocity_params, ρₐ, state::P3State; [use_aspect_ratio])

Returns a single-argument function `v_term(D)` that gives the Chen 2022
terminal velocity of an ice particle of maximum dimension `D`.

The size-independent coefficient work (`Chen2022_vel_coeffs` and
monodisperse-PDF construction for both small- and large-ice regimes) is
done once at call time; the returned closure only does the per-`D`
evaluation and — if `use_aspect_ratio = true` — the aspect-ratio
correction `cbrt(ϕᵢ(state, D))`. `ϕᵢ(state, D)` is O(1) because `state`
caches the regime thresholds.

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
    ρᵢ = FT(916.7)
    v_term_small = CO.particle_terminal_velocity(small_ice, ρₐ, ρᵢ)
    v_term_large = CO.particle_terminal_velocity(large_ice, ρₐ, ρᵢ)

    # Bare regime-split sedimentation velocity (no aspect-ratio correction).
    sedimentation_velocity(D) =
        D <= D_cutoff ? v_term_small(D) : v_term_large(D)
    # With aspect-ratio correction, composed on top.
    sedimentation_velocity_aspect_ratio(D) =
        cbrt(ϕᵢ(state, D)) * sedimentation_velocity(D)

    return use_aspect_ratio ? sedimentation_velocity_aspect_ratio : sedimentation_velocity
end
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
    use_aspect_ratio = true, p = 1e-6, quad = ChebyshevGauss(100)
)
    (; N_ice, L_ice) = state
    # TODO - do we want to swicth to ϵ_numerics(FT)
    if N_ice < eps(one(N_ice)) || L_ice < eps(one(L_ice))
        return zero(N_ice)
    end

    v_term = ice_particle_terminal_velocity(velocity_params, ρₐ, state; use_aspect_ratio)
    n = DT.size_distribution(state, logλ)

    # ∫n(D) v(D) dD
    number_weighted_integrand(D) = n(D) * v_term(D)

    bnds = integral_bounds(state, logλ; p)
    return integrate(number_weighted_integrand, bnds, quad) / N_ice
end

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
    (; N_ice, L_ice) = state
    # TODO - do we want to swicth to ϵ_numerics(FT)
    if N_ice < eps(one(N_ice)) || L_ice < eps(one(L_ice))
        return zero(L_ice)
    end

    v_term = ice_particle_terminal_velocity(velocity_params, ρₐ, state; use_aspect_ratio)
    n = DT.size_distribution(state, logλ)  # Number concentration at diameter D

    # ∫n(D) m(D) v(D) dD
    mass_weighted_integrand(D) = n(D) * v_term(D) * ice_mass(state, D)

    bnds = integral_bounds(state, logλ; p)
    return integrate(mass_weighted_integrand, bnds, quad) / L_ice
end

"""
    ice_terminal_velocity_number_weighted_from_prognostic(
        velocity_params, ρₐ, params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, logλ; kw...
    )

Pointwise wrapper that takes the *raw prognostic* P3 ice state
(`ρq_ice`, `ρn_ice`, `ρq_rim`, `ρb_rim`) and returns the number-weighted
mean ice terminal velocity. Builds the per-cell `P3State` via
[`get_state_from_prognostic`](@ref), so `F_rim` is regularised to
`[0, 1 - eps(FT)]` and `ρ_rim` is clamped to `[0, 0.8 ρ_l]`.

Designed for `@.`-broadcast use from a host (CA, KiD, etc.) where the
state must be reconstructed from prognostic variables every cell.
"""
@inline function ice_terminal_velocity_number_weighted_from_prognostic(
    velocity_params, ρₐ, params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, logλ; kw...,
)
    state = get_state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    return ice_terminal_velocity_number_weighted(velocity_params, ρₐ, state, logλ; kw...)
end

"""
    ice_terminal_velocity_mass_weighted_from_prognostic(
        velocity_params, ρₐ, params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, logλ; kw...
    )

Mass-weighted counterpart to
[`ice_terminal_velocity_number_weighted_from_prognostic`](@ref). Builds
the per-cell `P3State` via the regularised
[`get_state_from_prognostic`](@ref).
"""
@inline function ice_terminal_velocity_mass_weighted_from_prognostic(
    velocity_params, ρₐ, params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, logλ; kw...,
)
    state = get_state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    return ice_terminal_velocity_mass_weighted(velocity_params, ρₐ, state, logλ; kw...)
end
