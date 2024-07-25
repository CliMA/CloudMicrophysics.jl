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

"""
    mixed_phase_particle_velocity(p3, D, Chen2022, ρ_a, F_r, F_liq, th)

 - p3 - p3 parameters
 - Chen2022 - struct with terminal velocity parameters as in Chen(2022)
 - ρ_a - density of air
 - D - maximum particle dimension
 - F_r - rime mass fraction (q_rim/ q_i)
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - th - thresholds as calculated by thresholds()

Returns the terminal velocity of a single mixed-phase particle using the Chen 2022
parametrizations by computing an F_liq-weighted average of solid and liquid
phase terminal velocities.
"""
function mixed_phase_particle_velocity(
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
    ice_terminal_velocity(p3, Chen2022, q, N, ρ_r, F_liq, F_r, ρ_a)
 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struch with terminal velocity parameters as in Chen(2022)
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - F_r - rime mass fraction (q_rim/q_i)
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
    F_liq::FT,
    F_r::FT,
    ρ_a::FT,
) where {FT}
    # Get thresholds
    th = thresholds(p3, ρ_r, F_r)
    D_th = D_th_helper(p3)

    # Get size distribution parameters
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, F_liq, F_r)
    μ = DSD_μ(p3, λ)

    # Get integration bounds
    bound = get_ice_bound(p3, λ, eps(FT))

    v(D) = mixed_phase_particle_velocity(p3, D, Chen2022, ρ_a, F_r, F_liq, th)

    f(D) = N′ice(p3, D, λ, N_0) * v(D)

    f_m(D) = f(D) * p3_mass(p3, D, F_r, F_liq, th)

    v_n, err = QGK.quadgk(D -> f(D), FT(0), 2 * bound)

    v_m, err = QGK.quadgk(D -> f_m(D), FT(0), 2 * bound)

    return (v_n / N), (v_m / q)
end

"""
    velocity_difference(type, Dₗ, Dᵢ, p3, Chen2022, ρ_a, F_r, th)

 - type - a struct containing the size distribution parameters of the particle colliding with ice
 - Dₗ - maximum dimension of the particle colliding with ice
 - Dᵢ - maximum dimension of ice particle
 - p3 - a struct containing P3 parameters
 - Chen2022 - a struct containing Chen 2022 velocity parameters
 - ρ_a - density of air
 - F_r - rime mass fraction (q_rim/ q_i,ice)
 - F_liq - liquid fraction (q_liq/q_i,tot)
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
    F_r::FT,
    F_liq::FT,
    th,
) where {FT}
    return abs(
        mixed_phase_particle_velocity(p3, Dᵢ, Chen2022, ρ_a, F_r, F_liq, th) -
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
    F_r::FT,
    F_liq::FT,
    th,
) where {FT}
    # velocity difference for cloud collisions
    return abs(mixed_phase_particle_velocity(p3, Dᵢ, Chen2022, ρ_a, F_r, F_liq, th))
end
