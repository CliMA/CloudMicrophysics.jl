
"""
    ice_particle_terminal_velocity(state, velocity_params, ρₐ; [use_aspect_ratio])

Returns the terminal velocity of a single ice particle as a function of its size 
    (maximum dimension, `D`) using the Chen 2022 parametrization.

The first method returns a function of `D`, while the second evaluates at a specific `D`.

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

"""
    ice_terminal_velocity(dist, velocity_params, ρₐ; [use_aspect_ratio], [∫kwargs...])

Compute the mass and number weighted fall speeds for ice

See Eq. C10 of [MorrisonMilbrandt2015](@cite) and use the Chen 2022 terminal velocity scheme.

# Arguments
- `dist`: A [`P3Distribution`](@ref)
- `velocity_params`: A [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: The density of air [kg/m³]

# Keyword arguments
- `use_aspect_ratio`: Bool flag set to true if we want to consider the effects
   of particle aspect ratio on its terminal velocity (default: `true`)
- `∫kwargs`: Additional optional keyword arguments to pass to [`∫fdD`](@ref)

# Returns
- `v_n`: The number weighted fall speed
- `v_m`: The mass weighted fall speed
"""
function ice_terminal_velocity(
    dist::P3Distribution, velocity_params::CMP.Chen2022VelType, ρₐ;
    use_aspect_ratio = true, ∫kwargs...,
)
    (; state, L, N) = dist
    if N < eps(N) || L < eps(L)
        return zero(N), zero(L)
    end

    v_term = ice_particle_terminal_velocity(state, velocity_params, ρₐ; use_aspect_ratio)

    # ∫N(D) m(D) v(D) dD
    mass_weighted_integrand(D) = N′ice(dist, D) * v_term(D) * ice_mass(state, D)
    # ∫N(D) v(D) dD
    number_weighted_integrand(D) = N′ice(dist, D) * v_term(D)

    v_m = ∫fdD(mass_weighted_integrand, dist; ∫kwargs...)
    v_n = ∫fdD(number_weighted_integrand, dist; ∫kwargs...)
    return (v_n / N, v_m / L)
end
