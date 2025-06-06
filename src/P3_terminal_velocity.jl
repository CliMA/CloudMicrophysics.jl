
"""
    ice_particle_terminal_velocity(state, velocity_params, ρₐ; [use_aspect_ratio])

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
    state::P3State, velocity_params::CMP.Chen2022VelType, ρₐ; use_aspect_ratio = true,
)
    function v_term(D::FT) where {FT}
        (; small_ice, large_ice) = velocity_params
        D_cutoff = small_ice.cutoff # Diameter cutoff between small and large ice regimes
        ρᵢ = FT(916.7) # ρᵢ = p3_density(p3, D, F_rim, th) # TODO: tmp
        v_term_small = CO.particle_terminal_velocity(small_ice, ρₐ, ρᵢ)
        v_term_large = CO.particle_terminal_velocity(large_ice, ρₐ, ρᵢ)
        ϕ_factor = use_aspect_ratio ? cbrt(ϕᵢ(state, D)) : FT(1)
        vₜ = D <= D_cutoff ? v_term_small(D) : v_term_large(D)
        return ϕ_factor * vₜ
    end
    return v_term
end

function ice_terminal_velocity_number_weighted(
    state::P3State, velocity_params::CMP.Chen2022VelType, ρₐ;
    use_aspect_ratio = true, ∫kwargs...,
)
    if state.N_ice < eps(state.N_ice) || state.L_ice < eps(state.L_ice)
        return zero(state.N_ice)
    end

    v_term = ice_particle_terminal_velocity(state, velocity_params, ρₐ; use_aspect_ratio)

    # ∫N(D) v(D) dD
    number_weighted_integrand(D) = N′ice(state, logλ) * v_term(D)

    return ∫fdD(number_weighted_integrand, state, logλ; ∫kwargs...)
end

function ice_terminal_velocity_mass_weighted(
    state::P3State, velocity_params::CMP.Chen2022VelType, ρₐ, logλ;
    use_aspect_ratio = true, ∫kwargs...,
)
    (; N_ice, L_ice) = state
    if N_ice < eps(N_ice) || L_ice < eps(L_ice)
        return zero(L_ice)
    end

    v_term = ice_particle_terminal_velocity(state, velocity_params, ρₐ; use_aspect_ratio)
    N′ = N′ice(state, logλ)  # Number concentration at diameter D
    mass = get_mass_law(state)

    # ∫N(D) m(D) v(D) dD
    mass_weighted_integrand(D) = N′(D) * v_term(D) * mass(D)

    return ∫fdD(mass_weighted_integrand, state, logλ; ∫kwargs...)
end
