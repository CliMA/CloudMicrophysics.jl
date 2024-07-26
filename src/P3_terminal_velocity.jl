"""
    ϕᵢ(mᵢ, aᵢ, ρᵢ)

 - mᵢ - particle mass
 - aᵢ - particle area
 - ρᵢ - ice density

Returns the aspect ratio (ϕ) for an ice particle with given mass, area, and ice density.
"""
function ϕᵢ(mᵢ::FT, aᵢ::FT, ρᵢ::FT) where {FT}
    # TODO - add some notes on how we derived it
    # TODO - should we make use of other P3 properties like rimed fraction and volume?
    return 16 * ρᵢ^2 * aᵢ^3 / (9 * π * mᵢ^2)
end

"""
    ice_particle_terminal_velocity(D, Chen2022, ρₐ, mᵢ, aᵢ, ρᵢ)

 - D - maximum particle dimension
 - Chen2022 - a struct with terminal velocity parameters from Chen 2022
 - ρₐ - air density
 - mᵢ - particle mass
 - aᵢ - particle area
 - ρᵢ - ice density

Returns the terminal velocity of a single ice particle as a function
of its size (maximum dimension) using Chen 2022 parametrization.
Needed for numerical integrals in the P3 scheme.
"""
function ice_particle_terminal_velocity(
    D::FT,
    Chen2022::CMP.Chen2022VelTypeSnowIce,
    ρₐ::FT,
    mᵢ::FT,
    aᵢ::FT,
    ρᵢ::FT,
) where {FT}
    if D <= Chen2022.cutoff
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_small(Chen2022, ρₐ)
    else
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_large(Chen2022, ρₐ)
    end
    v = sum(@. sum(ai * D^bi * exp(-ci * D)))

    κ = FT(-1 / 6)
    return ϕᵢ(mᵢ, aᵢ, ρᵢ)^κ * v
end
# """
#     ice_particle_terminal_velocity(D, Chen2022, ρₐ)

#  - D - maximum particle dimension
#  - Chen2022 - a struct with terminal velocity parameters from Chen 2022
#  - ρₐ - air density

# Returns the terminal velocity of a single ice particle as a function
# of its size (maximum dimension) using Chen 2022 parametrization.
# Needed for numerical integrals in the P3 scheme.
# """
# function ice_particle_terminal_velocity(
#     D::FT,
#     Chen2022::CMP.Chen2022VelTypeSnowIce,
#     ρₐ::FT,
# ) where {FT}
#     if D <= Chen2022.cutoff
#         (ai, bi, ci) = CO.Chen2022_vel_coeffs_small(Chen2022, ρₐ)
#     else
#         (ai, bi, ci) = CO.Chen2022_vel_coeffs_large(Chen2022, ρₐ)
#     end
#     return sum(@. sum(ai * D^bi * exp(-ci * D)))
# end

"""
   p3_particle_terminal_velocity(p3, D, Chen2022, ρ_a, F_liq)

 - D - maximum particle dimension
 - Chen2022 - struct with terminal velocity parameters as in Chen(2022)
 - ρ_a - density of air
 - F_liq - liquid fraction (q_liq/q_i,tot)

Returns the terminal velocity of a single mixed-phase particle using the Chen 2022
parametrizations by computing an F_liq-weighted average of solid and liquid
phase terminal velocities.
"""
function p3_particle_terminal_velocity(
    p3::PSP3,
    D::FT,
    Chen2022::CMP.Chen2022VelType,
    ρ_a::FT,
    F_r::FT,
    F_liq::FT,
    th,
) where {FT}
    v_i = ice_particle_terminal_velocity(
        D,
        Chen2022.snow_ice,
        ρ_a,
        p3_mass(p3, D, F_r, FT(0), th),
        p3_area(p3, D, F_r, FT(0), th),
        p3.ρ_i,
    )
    v_r = CM2.rain_particle_terminal_velocity(D, Chen2022.rain, ρ_a)
    v = (1 - F_liq) * v_i + F_liq * v_r
    return v
end

"""
    velocity_difference(type, Dₗ, Dᵢ, Chen2022, ρ_a)

 - type - a struct containing the size distribution parameters of the particle colliding with ice
 - Dₗ - maximum dimension of the particle colliding with ice
 - Dᵢ - maximum dimension of ice particle
 - Chen2022 - a struct containing Chen 2022 velocity parameters
 - ρ_a - density of air
 - F_liq - liquid fraction (q_liq/q_i,tot)

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
    Chen2022::CMP.Chen2022VelType,
    ρ_a::FT,
    F_liq::FT,
) where {FT}
    # velocity difference for rain-ice collisions
    return abs(
        p3_particle_terminal_velocity(p3, D, Chen2022, ρₐ, F_rim, F_liq, th) -
        CM2.rain_particle_terminal_velocity(Dₗ, Chen2022.rain, ρ_a),
    )
end
function velocity_difference(
    pdf_c::CMP.CloudParticlePDF_SB2006{FT},
    Dₗ::FT,
    Dᵢ::FT,
    Chen2022::CMP.Chen2022VelType,
    ρ_a::FT,
) where {FT}
    # velocity difference for cloud-ice collisions
    return abs(p3_particle_terminal_velocity(p3, D, Chen2022, ρₐ, F_rim, F_liq, th))
end

"""
    ice_terminal_velocity(p3, Chen2022, q, N, ρ_r, F_rim, ρ_a)

 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struch with terminal velocity parameters as in Chen(2022)
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_rim - rime mass fraction (q_rim/q_i)
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - ρ_a - density of air

Returns the mass and number weighted fall speeds for ice following
eq C10 of Morrison and Milbrandt (2015).
"""
function ice_terminal_velocity(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    q::FT,
    N::FT,
    ρ_r::FT,
    F_rim::FT,
    F_liq::FT,
    ρₐ::FT,
) where {FT}

    # get the particle properties thresholds
    th = thresholds(p3, ρ_r, F_rim)
    # get the size distribution parameters
    (λ, N₀) = distribution_parameter_solver(p3, q, N, ρ_r, F_rim, F_liq)
    # get the integral limit
    D_max = get_ice_bound(p3, λ, eps(FT))

    # ∫N(D) m(D) v(D) dD
    v_m = QGK.quadgk(
        D ->
            N′ice(p3, D, λ, N₀) *
            p3_mass(p3, D, F_rim, F_liq, th) *
            p3_particle_terminal_velocity(p3, D, Chen2022, ρₐ, F_rim, F_liq, th),
        FT(0),
        D_max,
    )[1]

    # ∫N(D) v(D) dD
    v_n = QGK.quadgk(
        D ->
            N′ice(p3, D, λ, N₀) *
            p3_particle_terminal_velocity(p3, D, Chen2022, ρₐ, F_rim, F_liq, th),
        FT(0),
        D_max,
    )[1]
    return (v_n / N, v_m / q)
end
