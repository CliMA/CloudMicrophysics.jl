
"""
    ice_particle_terminal_velocity(velocity_params, ρₐ, state::P3State; [use_aspect_ratio])
    ice_particle_terminal_velocity(velocity_params, ρₐ, params::CMP.ParametersP3, F_rim, ρ_rim; [use_aspect_ratio])

Returns the terminal velocity of a single ice particle as a function of its size 
    (maximum dimension, `D`) using the Chen 2022 parametrization.

# Arguments
 - `state`: A [`P3State`](@ref)
 - `velocity_params`: A [`CMP.Chen2022VelType`](@ref) with terminal velocity parameters
 - `ρₐ`: Air density [kg/m³]

# Keyword arguments
 - `use_aspect_ratio`: Bool flag set to `true` if we want to consider the effects
    of particle aspect ratio on its terminal velocity (default: `true`)
"""
function ice_particle_terminal_velocity(
    velocity_params::CMP.Chen2022VelType, ρₐ, state_components...; use_aspect_ratio = true,
)
    function v_term(D::FT) where {FT}
        (; small_ice, large_ice) = velocity_params
        D_cutoff = small_ice.cutoff # Diameter cutoff between small and large ice regimes
        ρᵢ = FT(916.7) # ρᵢ = p3_density(p3, D, F_rim, th) # TODO: tmp
        v_term_small = CO.particle_terminal_velocity(small_ice, ρₐ, ρᵢ)
        v_term_large = CO.particle_terminal_velocity(large_ice, ρₐ, ρᵢ)
        # `state_components...` is either `state` or `(params, F_rim, ρ_rim)`
        ϕ_factor = use_aspect_ratio ? cbrt(ϕᵢ(state_components..., D)) : FT(1)
        vₜ = D <= D_cutoff ? v_term_small(D) : v_term_large(D)
        return ϕ_factor * vₜ
    end
    return v_term
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
