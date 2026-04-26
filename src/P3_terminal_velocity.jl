
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
        [use_aspect_ratio], [∫kwargs...],
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
 - `∫kwargs...`: Additional optional keyword arguments to pass to the quadrature rule

See also [`ice_terminal_velocity_mass_weighted`](@ref)
"""
function ice_terminal_velocity_number_weighted(
    velocity_params::CMP.Chen2022VelType, ρₐ, state::P3State, logλ;
    use_aspect_ratio = true, p = 1e-6, ∫kwargs...,
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
    return integrate(number_weighted_integrand, bnds...; ∫kwargs...) / N_ice
end
function ice_terminal_velocity_number_weighted(
    velocity_params::CMP.Chen2022VelType, ρₐ, params::CMP.ParametersP3, L_ice, N_ice, F_rim, ρ_rim, logλ;
    use_aspect_ratio = true, ∫kwargs...,
)
    state = get_state(params; L_ice, N_ice, F_rim, ρ_rim)
    return ice_terminal_velocity_number_weighted(velocity_params, ρₐ, state, logλ; use_aspect_ratio, ∫kwargs...)
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
 - `∫kwargs...`: Keyword arguments passed to the quadrature rule

See also [`ice_terminal_velocity_number_weighted`](@ref)
"""
function ice_terminal_velocity_mass_weighted(
    velocity_params::CMP.Chen2022VelType, ρₐ, state::P3State, logλ;
    use_aspect_ratio = true, p = 1e-6, ∫kwargs...,
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
    return integrate(mass_weighted_integrand, bnds...; ∫kwargs...) / L_ice
end
function ice_terminal_velocity_mass_weighted(
    velocity_params::CMP.Chen2022VelType, ρₐ, params::CMP.ParametersP3, L_ice, N_ice, F_rim, ρ_rim, logλ;
    use_aspect_ratio = true, ∫kwargs...,
)
    state = get_state(params; L_ice, N_ice, F_rim, ρ_rim)
    return ice_terminal_velocity_mass_weighted(velocity_params, ρₐ, state, logλ; use_aspect_ratio, ∫kwargs...)
end
