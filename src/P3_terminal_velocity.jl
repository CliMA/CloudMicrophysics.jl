"""
    ice_particle_terminal_velocity(state, D, Chen2022, ρₐ, use_aspect_ratio)

 - state - P3State
 - D - maximum particle dimension
 - Chen2022 - a struct with terminal velocity parameters from Chen 2022
 - ρₐ - air density
 - use_aspect_ratio - Bool flag set to true if we want to consider the effects
   of particle aspect ratio on its terminal velocity (default: true)

Returns the terminal velocity of a single ice particle as a function
of its size (maximum dimension, D) using the Chen 2022 parametrization.
"""
function ice_particle_terminal_velocity(
    state::P3State,
    D::FT,
    Chen2022::CMP.Chen2022VelType,
    ρₐ::FT,
    use_aspect_ratio = true,
) where {FT}
    # TODO - tmp
    #ρᵢ = p3_density(p3, D, F_rim, th)
    ρᵢ = FT(916.7)
    if D <= Chen2022.small_ice.cutoff
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_B2(Chen2022.small_ice, ρₐ, ρᵢ)
    else
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_B4(Chen2022.large_ice, ρₐ, ρᵢ)
    end
    v = sum(@. sum(ai * D^bi * exp(-ci * D)))

    return ifelse(use_aspect_ratio, ϕᵢ(state, D)^FT(1 / 3) * v, v)
end

"""
   p3_particle_terminal_velocity(p3, D, Chen2022, ρₐ, F_rim, F_liq, th, use_aspect_ratio)

 - p3 - p3 parameters
 - D - maximum particle dimension
 - Chen2022 - struct with terminal velocity parameters from Chen 2022
 - ρₐ - air density
 - F_rim - rime mass fraction (L_rim/L_ice) [-]
 - F_liq - liquid fraction (L_liq / L_p3_tot) [-]
 - th - thresholds as calculated by thresholds()
 - use_aspect_ratio - Bool flag set to true if we want to consider the effects
   of particle aspect ratio on its terminal velocity (default: true)

Returns the terminal velocity of a single mixed-phase particle using the Chen 2022
parametrizations by computing an F_liq-weighted average of solid and liquid
phase terminal velocities, using the maximum dimension of the whole particle for both.
"""
function p3_particle_terminal_velocity(
    state::P3State,
    D::FT,
    Chen2022::CMP.Chen2022VelType,
    ρₐ::FT,
    use_aspect_ratio = true,
) where {FT}
    v_i = ice_particle_terminal_velocity(state, D, Chen2022, ρₐ, use_aspect_ratio)
    v_r = CM2.rain_particle_terminal_velocity(D, Chen2022.rain, ρₐ)
    return weighted_average(state.F_liq, v_r, v_i)
end

"""
    velocity_difference(type, Dₗ, Dᵢ, p3, Chen2022, ρₐ, F_rim, F_liq, th, aspect_ratio)

 - type - a struct containing the size distribution parameters of the particle colliding with ice
 - Dₗ - maximum dimension of the particle colliding with ice
 - Dᵢ - maximum dimension of ice particle
 - p3 - a struct containing P3 parameters
 - Chen2022 - a struct containing Chen 2022 velocity parameters
 - ρ_a - density of air
 - F_rim - rime mass fraction (L_rim/L_ice) [-]
 - F_liq - liquid fraction (L_liq / L_p3_tot) [-]
 - th - P3 particle properties thresholds
 - use_aspect_ratio - Bool flag set to true if we want to consider the effects of
   particle aspect ratio on its terminal velocity (default: true)

Returns the absolute value of the velocity difference between a mixed-phase particle and
cloud or rain drop as a function of their sizes. It uses Chen 2022 velocity
parameterization for ice and rain and assumes no sedimentation of cloud droplets.
"""
function velocity_difference(
    pdf_r::Union{
        CMP.RainParticlePDF_SB2006{FT},
        CMP.RainParticlePDF_SB2006_limited{FT},
    },
    Dₗ::FT,
    Dᵢ::FT,
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    ρₐ::FT,
    F_rim::FT,
    th,
    use_aspect_ratio = true,
) where {FT}
    # velocity difference for rain-ice collisions
    return abs(
        p3_particle_terminal_velocity(
            p3,
            Dᵢ,
            Chen2022,
            ρₐ,
            F_rim,
            F_liq,
            th,
            use_aspect_ratio,
        ) - CM2.rain_particle_terminal_velocity(Dₗ, Chen2022.rain, ρₐ),
    )
end
function velocity_difference(
    pdf_c::CMP.CloudParticlePDF_SB2006{FT},
    Dₗ::FT,
    Dᵢ::FT,
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    ρₐ::FT,
    F_rim::FT,
    F_liq::FT,
    th,
    use_aspect_ratio = true,
) where {FT}
    # velocity difference for cloud-ice collisions
    return abs(
        p3_particle_terminal_velocity(
            p3,
            Dᵢ,
            Chen2022,
            ρₐ,
            F_rim,
            F_liq,
            th,
            use_aspect_ratio,
        ),
    )
end

"""
    ice_terminal_velocity(dist, Chen2022, ρₐ, use_aspect_ratio)

 - dist - a struct with P3 distribution parameters
 - Chen2022 - a struct with terminal velocity parameters as in Chen(2022)
 - ρₐ - density of air
 - use_aspect_ratio - Bool flag set to true if we want to consider the effects
   of particle aspect ratio on its terminal velocity (default: true)

Returns the mass and number weighted fall speeds for ice following
eq C10 of Morrison and Milbrandt (2015) and using Chen 2022 terminal velocity scheme.
"""
function ice_terminal_velocity(dist::P3Distribution{FT}, Chen2022::CMP.Chen2022VelType, ρₐ::FT, use_aspect_ratio = true) where {FT}
    L = exp(log_LdN₀(dist.state, dist.log_λ) + dist.log_N₀)
    N = exp(log_NdN₀(dist.state, dist.log_λ) + dist.log_N₀)
    if N < eps(FT) || L < eps(FT)
        return FT(0), FT(0)
    end

    # Fetch the particle properties thresholds
    thresholds = threshold_tuple(dist.state)

    # ∫N(D) m(D) v(D) dD
    v_m = QGK.quadgk(
        D ->
            N′ice(dist, D) *
            ice_mass(dist.state, D) *
            p3_particle_terminal_velocity(dist.state, D, Chen2022, ρₐ, use_aspect_ratio),
        FT(0), thresholds..., FT(1), # non-Inf upper limit
        rtol = FT(1e-6),
    )[1]

    # ∫N(D) v(D) dD
    v_n = QGK.quadgk(
        D ->
            N′ice(dist, D) *
            p3_particle_terminal_velocity(dist.state, D, Chen2022, ρₐ, use_aspect_ratio),
        FT(0), thresholds..., FT(1), # non-Inf upper limit
        rtol = FT(1e-6),
    )[1]
    return (v_n / N, v_m / L)
end
