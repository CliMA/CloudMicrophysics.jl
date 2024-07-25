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
    ice_terminal_velocity(p3, Chen2022, q, N, ρ_r, F_r, ρ_a)

 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struch with terminal velocity parameters as in Chen(2022)
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i)
 - ρ_a - density of air

Returns the mass and number weighted fall speeds for ice following
eq C10 of Morrison and Milbrandt (2015).
"""
function ice_terminal_velocity(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelTypeSnowIce,
    q::FT,
    N::FT,
    ρ_r::FT,
    F_r::FT,
    ρ_a::FT,
) where {FT}
    # Get the pree parameters for terminal velocities of small
    # and large particles
    small = CO.Chen2022_vel_coeffs_small(Chen2022, ρ_a)
    large = CO.Chen2022_vel_coeffs_large(Chen2022, ρ_a)
    get_p(prs, it) = (prs[1][it], prs[2][it], prs[3][it])

    # Get the thresholds for different particles regimes
    (; D_cr, D_gr, ρ_g, ρ_d) = thresholds(p3, ρ_r, F_r)
    D_th = D_th_helper(p3)
    D_ct = Chen2022.cutoff

    # Get the shape parameters of the particle size distribution
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, F_r)
    μ = DSD_μ(p3, λ)

    # TODO: Change when each value used depending on type of particle
    # TODO: or keep fixed and add to ClimaParams...?
    κ = FT(-1 / 6) #FT(1/3)
    # Redefine α_va to be in si units
    α_va = α_va_si(p3)

    aₛ(a) = a * N_0
    bₛ(b) = b + μ
    cₛ(c) = c + λ

    aₛ_m(a) = aₛ(a) * FT(π) / 6 * p3.ρ_i
    bₛ_m(b) = bₛ(b) + 3

    spheres_n(a, b, c) = (aₛ(a), bₛ(b), cₛ(c))
    spheres_m(a, b, c) = (aₛ_m(a), bₛ_m(b), cₛ(c))

    # TODO - combine different aspect ratio computations
    aₙₛ(a) = aₛ(a) * (16 * p3.ρ_i^2 * p3.γ^3 / (9 * FT(π) * α_va^2))^κ
    bₙₛ(b) = bₛ(b) + κ * (3 * p3.σ - 2 * p3.β_va)

    aₙₛ_m(a) = aₙₛ(a) * α_va
    bₙₛ_m(b) = bₙₛ(b) + p3.β_va

    non_spheres_n(a, b, c) = (aₙₛ(a), bₙₛ(b), cₛ(c))
    non_spheres_m(a, b, c) = (aₙₛ_m(a), bₙₛ_m(b), cₛ(c))

    aᵣₛ(a) = aₛ(a) * (p3.ρ_i / ρ_g)^(2 * κ)
    aᵣₛ_m(a) = aᵣₛ(a) * FT(π) / 6 * ρ_g

    rimed_n(a, b, c) = (aᵣₛ(a), bₛ(b), cₛ(c))
    rimed_m(a, b, c) = (aᵣₛ_m(a), bₛ_m(b), cₛ(c))

    v_n_D_cr(D, a, b, c) =
        a *
        N_0 *
        D^(b + μ) *
        exp((-c - λ) * D) *
        (
            16 * p3.ρ_i^2 * (F_r * π / 4 * D^2 + (1 - F_r) * p3.γ * D^p3.σ)^3 /
            (9 * π * (α_va / (1 - F_r) * D^p3.β_va)^2)
        )^κ
    v_m_D_cr(D, a, b, c) = v_n_D_cr(D, a, b, c) * (α_va / (1 - F_r) * D^p3.β_va)

    v_m = 0
    v_n = 0
    for i in 1:2
        if F_r == 0
            v_m += ∫_Γ(FT(0), D_th, spheres_m(get_p(small, i)...)...)
            v_n += ∫_Γ(FT(0), D_th, spheres_n(get_p(small, i)...)...)

            v_m += ∫_Γ(
                D_th,
                D_ct,
                Inf,
                non_spheres_m(get_p(small, i)...)...,
                non_spheres_m(get_p(large, i)...)...,
            )
            v_n += ∫_Γ(
                D_th,
                D_ct,
                Inf,
                non_spheres_n(get_p(small, i)...)...,
                non_spheres_n(get_p(large, i)...)...,
            )
        else
            # Velocity coefficients for small particles
            v_m += ∫_Γ(FT(0), D_th, spheres_m(get_p(small, i)...)...)
            v_n += ∫_Γ(FT(0), D_th, spheres_n(get_p(small, i)...)...)
            is_large = false

            # D_th to D_gr
            if !is_large && D_gr > D_ct
                v_m += ∫_Γ(
                    D_th,
                    D_ct,
                    D_gr,
                    non_spheres_m(get_p(small, i)...)...,
                    non_spheres_m(get_p(large, i)...)...,
                )
                v_n += ∫_Γ(
                    D_th,
                    D_ct,
                    D_gr,
                    non_spheres_n(get_p(small, i)...)...,
                    non_spheres_n(get_p(large, i)...)...,
                )
                # Switch to large particles
                is_large = true
            else
                v_m += ∫_Γ(D_th, D_gr, non_spheres_m(get_p(small, i)...)...)
                v_n += ∫_Γ(D_th, D_gr, non_spheres_n(get_p(small, i)...)...)
            end

            # D_gr to D_cr
            if !is_large && D_cr > D_ct
                v_m += ∫_Γ(
                    D_gr,
                    D_ct,
                    D_cr,
                    rimed_m(get_p(small, i)...)...,
                    rimed_m(get_p(large, i)...)...,
                )
                v_n += ∫_Γ(
                    D_gr,
                    D_ct,
                    D_cr,
                    rimed_n(get_p(small, i)...)...,
                    rimed_n(get_p(large, i)...)...,
                )
                # Switch to large particles
                is_large = true
            elseif is_large
                v_m += ∫_Γ(D_gr, D_cr, rimed_m(get_p(large, i)...)...)
                v_n += ∫_Γ(D_gr, D_cr, rimed_n(get_p(large, i)...)...)
            else
                v_m += ∫_Γ(D_gr, D_cr, rimed_m(get_p(small, i)...)...)
                v_n += ∫_Γ(D_gr, D_cr, rimed_n(get_p(small, i)...)...)
            end

            # D_cr to Infinity
            if !is_large
                (Im, em) =
                    QGK.quadgk(D -> v_m_D_cr(D, get_p(small, i)...), D_cr, D_ct)
                (In, en) =
                    QGK.quadgk(D -> v_n_D_cr(D, get_p(small, i)...), D_cr, D_ct)
                v_m += Im
                v_n += In

                # Switch to large particles
                (Im, em) =
                    QGK.quadgk(D -> v_m_D_cr(D, get_p(large, i)...), D_ct, Inf)
                (In, en) =
                    QGK.quadgk(D -> v_n_D_cr(D, get_p(large, i)...), D_ct, Inf)
                v_m += Im
                v_n += In
            else
                # TODO - check if it should be large or small
                (Im, em) =
                    QGK.quadgk(D -> v_m_D_cr(D, get_p(large, i)...), D_cr, Inf)
                (In, en) =
                    QGK.quadgk(D -> v_n_D_cr(D, get_p(large, i)...), D_cr, Inf)
                v_m += Im
                v_n += In
            end
        end
    end
    return (v_n / N, v_m / q)
end

"""
    velocity_difference(type, Dₗ, Dᵢ, p3, Chen2022, ρ_a, F_r, th)

 - type - a struct containing the size distribution parameters of the particle colliding with ice
 - Dₗ - maximum dimension of the particle colliding with ice
 - Dᵢ - maximum dimension of ice particle
 - p3 - a struct containing P3 parameters
 - Chen2022 - a struct containing Chen 2022 velocity parameters
 - ρ_a - density of air
 - F_r - rime mass fraction (q_rim/ q_i)
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
    th,
) where {FT}
    return abs(
        ice_particle_terminal_velocity(
            Dᵢ,
            Chen2022.snow_ice,
            ρ_a,
            p3_mass(p3, Dᵢ, F_r, th),
            p3_area(p3, Dᵢ, F_r, th),
            p3.ρ_i,
        ) - CM2.rain_particle_terminal_velocity(Dₗ, Chen2022.rain, ρ_a),
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
    th,
) where {FT}
    # velocity difference for cloud collisions
    return abs(
        ice_particle_terminal_velocity(
            Dᵢ,
            Chen2022.snow_ice,
            ρ_a,
            p3_mass(p3, Dᵢ, F_r, th),
            p3_area(p3, Dᵢ, F_r, th),
            p3.ρ_i,
        ),
    )
end
