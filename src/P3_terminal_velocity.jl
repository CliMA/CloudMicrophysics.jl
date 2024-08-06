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
    ice_particle_terminal_velocity(D, Chen2022, ρₐ, mᵢ, aᵢ, ρᵢ, aspect_ratio)
 - D - maximum particle dimension
 - Chen2022 - a struct with terminal velocity parameters from Chen 2022
 - ρₐ - air density
 - mᵢ - particle mass
 - aᵢ - particle area
 - ρᵢ - ice density
 - aspect_ratio - Boolean which determines whether or not aspect ratio is used
    - default: true

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
    mᵢ::FT,
    aᵢ::FT,
    ρᵢ::FT,
    aspect_ratio = true,
) where {FT}
    if D <= Chen2022.cutoff
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_small(Chen2022, ρₐ)
    else
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_large(Chen2022, ρₐ)
    end
    v = sum(@. sum(ai * D^bi * exp(-ci * D)))

    κ = FT(-1 / 6)
    if aspect_ratio
        return ifelse(D == 0, FT(0), ϕᵢ(mᵢ, aᵢ, ρᵢ)^κ * v)
    else
        return v
    end
end

"""
    velocity_difference(type, Dₗ, Dᵢ, p3, Chen2022, ρ_a, F_r, th, aspect_ratio)

 - type - a struct containing the size distribution parameters of the particle colliding with ice
 - Dₗ - maximum dimension of the particle colliding with ice
 - Dᵢ - maximum dimension of ice particle
 - p3 - a struct containing P3 parameters
 - Chen2022 - a struct containing Chen 2022 velocity parameters
 - ρ_a - density of air
 - F_r - rime mass fraction (L_rim/ L_ice)
 - th - P3 particle properties thresholds
 - aspect_ratio - Boolean which determines whether or not aspect ratio is used
    - default: true

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
    ρₐ::FT,
    F_r::FT,
    th,
    aspect_ratio = true,
) where {FT}
    # velocity difference for rain-ice collisions
    return abs(
        ice_particle_terminal_velocity(Dᵢ, Chen2022.snow_ice, ρₐ, p3_mass(p3, Dᵢ, F_r, th), p3_area(p3, Dᵢ, F_r, th), p3.ρ_i, aspect_ratio) -
        CM2.rain_particle_terminal_velocity(Dₗ, Chen2022.rain, ρₐ),
    )
end
function velocity_difference(
    pdf_c::CMP.CloudParticlePDF_SB2006{FT},
    Dₗ::FT,
    Dᵢ::FT,
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    ρₐ::FT,
    F_r::FT,
    th,
    aspect_ratio = true,
) where {FT}
    # velocity difference for cloud-ice collisions
    return abs(ice_particle_terminal_velocity(Dᵢ, Chen2022.snow_ice, ρₐ, p3_mass(p3, Dᵢ, F_r, th), p3_area(p3, Dᵢ, F_r, th), p3.ρ_i, aspect_ratio))
end

"""
    ice_terminal_velocity(p3, Chen2022, q, N, ρ_r, F_r, ρ_a, aspect_ratio)

 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struch with terminal velocity parameters as in Chen(2022)
 - L - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (L_rim/q_ice)
 - ρ_a - density of air
 - aspect_ratio - Boolean which determines whether or not aspect ratio is used
    - default: true

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
    aspect_ratio = true,
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
            ice_particle_terminal_velocity(D, Chen2022, ρₐ, p3_mass(p3, D, F_r, th), p3_area(p3, D, F_r, th), p3.ρ_i, aspect_ratio),
        FT(0),
        D_max,
    )[1]

    # ∫N(D) v(D) dD
    v_n = QGK.quadgk(
        D ->
            N′ice(p3, D, λ, N₀) *
            ice_particle_terminal_velocity(D, Chen2022, ρₐ, p3_mass(p3, D, F_r, th), p3_area(p3, D, F_r, th), p3.ρ_i, aspect_ratio),
        FT(0),
        D_max,
    )[1]
    return (v_n / N, v_m / L)
end
