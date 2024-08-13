"""
    ϕᵢ(p3, D, F_rim, th)

 - p3 - a struct containing P3 parameters
 - D - maximum dimension of ice particle [m]
 - F_rim - rime mass fraction (L_rim/ L_ice) [-]
 - th - P3 particle properties thresholds

Returns the aspect ratio (ϕ) for an ice particle with mass, cross-section area,
and ice density determined using the size-dependent particle property regimes
following Morrison and Milbrandt (2015). The density of nonspherical
particles is assumed to be equal to the particle mass divided by the volume of a
spherical particle with the same D_max.
"""
function ϕᵢ(p3::PSP3, D::FT, F_rim::FT, th) where {FT}
    # Let F_liq = 0 to calculate ϕᵢ for an ice particle with D
    # of the total particle but treating a mixed-phase particle
    # as an ice particle with the same D. We need this for the F_liq
    # weighted average of velocities discussed in the docs.
    # TODO? should the above comment be here or should it be in the docs?
    F_liq = FT(0)
    mᵢ = p3_mass(p3, D, F_rim, F_liq, th)
    aᵢ = p3_area(p3, D, F_rim, F_liq, th)
    ρᵢ = p3_density(p3, D, F_rim, th)

    return ifelse(D == 0, FT(0), 16 * ρᵢ^2 * aᵢ^3 / (9 * FT(π) * mᵢ^2))
end

"""
    ice_particle_terminal_velocity(p3, D, Chen2022, ρₐ, F_rim, th, use_aspect_ratio)

 - p3 - p3 parameters
 - D - maximum particle dimension
 - Chen2022 - a struct with terminal velocity parameters from Chen 2022
 - ρₐ - air density
 - F_rim - rime mass fraction (L_rim/(L_ice - L_liq)) [-]
 - th - P3 particle properties thresholds
 - use_aspect_ratio - Bool flag set to true if we want to consider the effects
   of particle aspect ratio on its terminal velocity (default: true)

Returns the terminal velocity of a single ice particle as a function
of its size (maximum dimension, D) using the Chen 2022 parametrization.
"""
function ice_particle_terminal_velocity(
    p3::PSP3,
    D::FT,
    Chen2022::CMP.Chen2022VelTypeSnowIce,
    ρₐ::FT,
    F_rim::FT,
    th,
    use_aspect_ratio = true,
) where {FT}
    if D <= Chen2022.cutoff
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_small(Chen2022, ρₐ)
    else
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_large(Chen2022, ρₐ)
    end
    v = sum(@. sum(ai * D^bi * exp(-ci * D)))

    return ifelse(use_aspect_ratio, ϕᵢ(p3, D, F_rim, th)^FT(1 / 3) * v, v)
end

"""
   p3_particle_terminal_velocity(p3, D, Chen2022, ρₐ, F_rim, F_liq, th, use_aspect_ratio)

 - p3 - p3 parameters
 - D - maximum particle dimension
 - Chen2022 - struct with terminal velocity parameters from Chen 2022
 - ρₐ - air density
 - F_rim - rime mass fraction (L_rim/(L_ice - L_liq)) [-]
 - F_liq - liquid fraction (L_liq/L_ice) [-]
 - th - thresholds as calculated by thresholds()
 - use_aspect_ratio - Bool flag set to true if we want to consider the effects
   of particle aspect ratio on its terminal velocity (default: true)

Returns the terminal velocity of a single mixed-phase particle using the Chen 2022
parametrizations by computing an F_liq-weighted average of solid and liquid
phase terminal velocities, using the maximum dimension of the whole particle for both.
"""
function p3_particle_terminal_velocity(
    p3::PSP3,
    D::FT,
    Chen2022::CMP.Chen2022VelType,
    ρₐ::FT,
    F_rim::FT,
    F_liq::FT,
    th,
    use_aspect_ratio = true,
) where {FT}
    v_i = ice_particle_terminal_velocity(
        p3,
        D,
        Chen2022.snow_ice,
        ρₐ,
        F_rim,
        th,
        use_aspect_ratio,
    )
    v_r = CM2.rain_particle_terminal_velocity(D, Chen2022.rain, ρₐ)
    v = (1 - F_liq) * v_i + F_liq * v_r
    return v
end

"""
    velocity_difference(type, Dₗ, Dᵢ, p3, Chen2022, ρₐ, F_rim, F_liq, th, aspect_ratio)

 - type - a struct containing the size distribution parameters of the particle colliding with ice
 - Dₗ - maximum dimension of the particle colliding with ice
 - Dᵢ - maximum dimension of ice particle
 - p3 - a struct containing P3 parameters
 - Chen2022 - a struct containing Chen 2022 velocity parameters
 - ρ_a - density of air
 - F_rim - rime mass fraction (L_rim/(L_ice - L_liq)) [-]
 - F_liq - liquid fraction (L_liq/L_ice) [-]
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
    ice_terminal_velocity(p3, Chen2022, L, N, ρ_r, F_rim, ρₐ, use_aspect_ratio)

 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struch with terminal velocity parameters as in Chen(2022)
 - L - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_rim - rime mass fraction (L_rim/q_ice)
 - ρₐ - density of air
 - use_aspect_ratio - Bool flag set to true if we want to consider the effects
   of particle aspect ratio on its terminal velocity (default: true)

Returns the mass and number weighted fall speeds for ice following
eq C10 of Morrison and Milbrandt (2015) and using Chen 2022 terminal velocity scheme.
"""
function ice_terminal_velocity(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    L::FT,
    N::FT,
    ρ_r::FT,
    F_rim::FT,
    F_liq::FT,
    ρₐ::FT,
    use_aspect_ratio = true,
) where {FT}
    if N < eps(FT) || L < eps(FT)
        return FT(0), FT(0)
    else
        # get the particle properties thresholds
        th = thresholds(p3, ρ_r, F_rim)
        # get the size distribution parameters
        (λ, N₀) = distribution_parameter_solver(p3, L, N, ρ_r, F_rim, F_liq)
        # get the integral limit
        D_max = get_ice_bound(p3, λ, eps(FT))

        # ∫N(D) m(D) v(D) dD
        v_m = QGK.quadgk(
            D ->
                N′ice(p3, D, λ, N₀) *
                p3_mass(p3, D, F_rim, F_liq, th) *
                p3_particle_terminal_velocity(
                    p3,
                    D,
                    Chen2022,
                    ρₐ,
                    F_rim,
                    F_liq,
                    th,
                    use_aspect_ratio,
                ),
            FT(0),
            D_max,
            rtol = FT(1e-6),
        )[1]

        # ∫N(D) v(D) dD
        v_n = QGK.quadgk(
            D ->
                N′ice(p3, D, λ, N₀) * p3_particle_terminal_velocity(
                    p3,
                    D,
                    Chen2022,
                    ρₐ,
                    F_rim,
                    F_liq,
                    th,
                    use_aspect_ratio,
                ),
            FT(0),
            D_max,
            rtol = FT(1e-6),
        )[1]
        return (v_n / N, v_m / L)
    end
end
