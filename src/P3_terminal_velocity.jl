
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
        ρᵢ = FT(916.7) # ρᵢ = p3_density(p3, D, F_rim, th) # TODO: tmp
        ϕ_factor = use_aspect_ratio ? cbrt(ϕᵢ(state, D)) : FT(1)
        (ai, bi, ci) = if D <= small_ice.cutoff
            CO.Chen2022_vel_coeffs_B2(small_ice, ρₐ, ρᵢ)
        else
            CO.Chen2022_vel_coeffs_B4(large_ice, ρₐ, ρᵢ)
        end
        ϕ_factor * sum(@. sum(ai * D^bi * exp(-ci * D)))
    end
    return v_term
end

"""
    ice_terminal_velocity(dist, velocity_params, ρₐ; [use_aspect_ratio], [accurate])

Compute the mass and number weighted fall speeds for ice

See Eq. C10 of [MorrisonMilbrandt2015](@cite) and use the Chen 2022 terminal velocity scheme.

# Arguments
- `dist`: A [`P3Distribution`](@ref)
- `velocity_params`: A [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: The density of air [kg/m³]

# Keyword arguments
- `use_aspect_ratio`: Bool flag set to true if we want to consider the effects
   of particle aspect ratio on its terminal velocity (default: `true`)
- `accurate`: Set to `true` to perform a more accurate numerical integration
    see [`∫fdD`](@ref) for details. Default is `false`.

# Returns
- `v_n`: The number weighted fall speed
- `v_m`: The mass weighted fall speed
"""
function ice_terminal_velocity(
    dist::P3Distribution, velocity_params::CMP.Chen2022VelType, ρₐ;
    use_aspect_ratio = true, accurate = false,
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

    v_m = ∫fdD(mass_weighted_integrand, state; accurate)
    v_n = ∫fdD(number_weighted_integrand, state; accurate)
    return (v_n / N, v_m / L)
end
