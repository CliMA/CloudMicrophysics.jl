"""
    ice_particle_terminal_velocity(D, Chen2022, ρₐ)

 - D - maximum particle dimension
 - Chen2022 - a struct with terminal velocity parameters from Chen 2022
 - ρₐ - air density

Returns the terminal velocity of a single ice particle as a function
of its size (maximum dimension) using Chen 2022 parametrization.
Needed for numerical integrals in the P3 scheme.
"""
function ice_particle_terminal_velocity(
    D::FT,
    Chen2022::CMP.Chen2022VelTypeSnowIce,
    ρₐ::FT,
) where {FT}
    if D <= Chen2022.cutoff
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_small(Chen2022, ρₐ)
    else
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_large(Chen2022, ρₐ)
    end
    return sum(@. sum(ai * D^bi * exp(-ci * D)))
end

"""
    velocity_difference(type, Dₗ, Dᵢ, p3, Chen2022, ρ_a, F_r, th)

 - type - a struct containing the size distribution parameters of the particle colliding with ice
 - Dₗ - maximum dimension of the particle colliding with ice
 - Dᵢ - maximum dimension of ice particle
 - p3 - a struct containing P3 parameters
 - Chen2022 - a struct containing Chen 2022 velocity parameters
 - ρ_a - density of air
 - F_r - rime mass fraction (L_rim/ L_ice)
 - th - P3 particle properties thresholds

Returns the absolute value of the velocity difference between an ice particle and
cloud or rain drop as a function of their sizes. It uses Chen 2022 velocity
parameterization for ice and rain and assumes no sedimentation of cloud droplets.
Needed for P3 numerical integrals
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
    ρ_a::FT,
) where {FT}
    # velocity difference for rain-ice collisions
    return abs(
        ice_particle_terminal_velocity(Dᵢ, Chen2022.snow_ice, ρ_a) -
        CM2.rain_particle_terminal_velocity(Dₗ, Chen2022.rain, ρ_a),
    )
end
function velocity_difference(
    pdf_c::CMP.CloudParticlePDF_SB2006{FT},
    Dₗ::FT,
    Dᵢ::FT,
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    ρ_a::FT,
) where {FT}
    # velocity difference for cloud-ice collisions
    return abs(ice_particle_terminal_velocity(Dᵢ, Chen2022.snow_ice, ρ_a))
end

"""
    ice_terminal_velocity_2(p3, Chen2022, q, N, ρ_r, F_r, ρ_a)

 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struch with terminal velocity parameters as in Chen(2022)
 - L - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (L_rim/q_ice)
 - ρ_a - density of air

Returns the mass and number weighted fall speeds for ice following
eq C10 of Morrison and Milbrandt (2015).
"""
function ice_terminal_velocity(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelTypeSnowIce,
    L::FT,
    N::FT,
    ρ_r::FT,
    F_r::FT,
    ρₐ::FT,
) where {FT}

    # get the particle properties thresholds
    th = thresholds(p3, ρ_r, F_r)
    # get the size distribution parameters
    (λ, N₀) = distribution_parameter_solver(p3, L, N, ρ_r, F_r)
    # get the integral limit
    D_max = get_ice_bound(p3, λ, eps(FT))

    # ∫N(D) m(D) v(D) dD
    v_m = QGK.quadgk(
        D ->
            N′ice(p3, D, λ, N₀) *
            p3_mass(p3, D, F_r, th) *
            ice_particle_terminal_velocity(D, Chen2022, ρₐ),
        FT(0),
        D_max,
    )[1]

    # ∫N(D) v(D) dD
    v_n = QGK.quadgk(
        D ->
            N′ice(p3, D, λ, N₀) *
            ice_particle_terminal_velocity(D, Chen2022, ρₐ),
        FT(0),
        D_max,
    )[1]
    return (v_n / N, v_m / L)
end
