"""
Predicted particle properties scheme (P3) for ice, which includes:
 - threshold solver
 - shape parameters solver
 - m(D) regime
 - a(D) regime

Implementation of Morrison and Milbrandt 2015 doi: 10.1175/JAS-D-14-0065.1

Note: Particle size is defined as its maximum length (i.e. max dimesion).

Changes to accomodate for liquid fraction: shape params:
    - added mass_liq
    - for p3_mass and p3_area, added a F_liq-weighted linear
    average of the original regime and the liquid part (assumed
    spherical for area)
    - for mass mixing ratios in different regimes:
        - modified each q_x such that q_xnew = (1 - F_liq) * q_x
        - added q_liq (mass mixing ratio of liquid on mixed-phase particles)
    - with nonzero F_liq, and with the new q_x formulations, the 
    q_over_N_gamma helper function in the shape parameter solver now includes
    for the mixing ratio q in the numerator q_liq as well as q_[insert regime]
        - for F_liq = 0, q_new = q_old, which is what we want
        - for nonzero F_liq, q_new = q_old + q_liq
    - so, in theory, if we want the shape parameters for the ice cores we pass
    F_liq = 0, even when it is nonzero, and if we want the shape parameters for the whole
    particle, we pass F_liq
    - ***with this approach, I'm unsure how to change the bounds we use for q_liq to
    correspond to the different cases in q_over_N_gamma***
    - ***I don't believe the μ given λ and vice versa need to be changed, but perhaps the
    get_bounds function will need some fine tuning, assuming the bounds may change given
    PSDs will be different under the liquid fraction regime***
For terminal velocity:
    - kept original terminal_velocity function with F_liq = 0 to
    calculate terminal velocity of ice cores
    - added/adding terminal_velocity_liq which uses Chen2022VelTypeRain
    - added/adding terminal_velocity_tot which computes an F_liq-weighted
    linear average of the terminal velocities of the ice core and liquid part
"""
module P3Scheme

import SpecialFunctions as SF
import QuadGK as QGK
import RootSolvers as RS
import HCubature as HC

import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics.TerminalVelocity as TV

const PSP3 = CMP.ParametersP3

export thresholds, distribution_parameter_solver

"""
    α_va_si(p3)

 - p3 - a struct with P3 scheme parameters

Returns `α_va` coefficient for the assumed particle mass(size) relation for
large unrimed ice and dense nonspherical ice, in base SI units: kg m^(-β_va).
`β_va` is another coefficient of the mass(size) relation.
From measurements of mass grown by vapor diffusion and aggregation
in midlatitude cirrus by Brown and Francis (1995)
doi: 10.1175/1520-0426(1995)012<0410:IMOTIW>2.0.CO;2
"""
α_va_si(p3::PSP3{FT}) where {FT} = p3.α_va * 10^(6 * p3.β_va - 3)

"""
    D_th_helper(p3)

 - p3 - a struct with P3 scheme parameters

Returns the critical size separating spherical and nonspherical ice, in meters.
Eq. 8 in Morrison and Milbrandt (2015).
"""
D_th_helper(p3::PSP3{FT}) where {FT} =
    (FT(π) * p3.ρ_i / 6 / α_va_si(p3))^(1 / (p3.β_va - 3))

"""
    D_cr_helper(p3, F_r, ρ_g)

 - p3 - a struct with P3 scheme parameters
 - F_r - rime mass fraction (q_rim/q_i) [-]
 - ρ_g - is the effective density of a spherical graupel particle [kg/m^3]

Returns the size of equal mass for graupel and partially rimed ice, in meters.
Eq. 14 in Morrison and Milbrandt (2015).
"""
function D_cr_helper(p3::PSP3{FT}, F_r::FT, ρ_g::FT) where {FT}
    α_va = α_va_si(p3)
    return (1 / (1 - F_r) * 6 * α_va / FT(π) / ρ_g)^(1 / (3 - p3.β_va))
end

"""
    D_gr_helper(p3, ρ_g)

 - p3 - a struct with P3 scheme parameters
 - ρ_g - is the effective density of a spherical graupel particle [kg/m^3]

Returns the size of equal mass for graupel and unrimed ice, in meters.
Eq. 15 in Morrison and Milbrandt (2015).
"""
function D_gr_helper(p3::PSP3{FT}, ρ_g::FT) where {FT}
    α_va = α_va_si(p3)
    return (6 * α_va / FT(π) / ρ_g)^(1 / (3 - p3.β_va))
end

"""
    ρ_g_helper(ρ_r, F_r, ρ_d)

 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i) [-]
 - ρ_g - is the effective density of a spherical graupel particle [kg/m^3]

Returns the density of total (deposition + rime) ice mass for graupel, in kg/m3
Eq. 16 in Morrison and Milbrandt (2015).
"""
ρ_g_helper(ρ_r::FT, F_r::FT, ρ_d::FT) where {FT} = F_r * ρ_r + (1 - F_r) * ρ_d

"""
    ρ_d_helper(p3, D_cr, D_gr)

 - p3 - a struct with P3 scheme parameters
 - D_cr - is the size of equal mass for graupel and partially rimed ice, in meters
 - D_gr - the size of equal mass for graupel and unrimed ice, in meters

Returns the density of unrimed ice mass, in kg/m3
Eq. 17 in Morrison and Milbrandt (2015).
"""
function ρ_d_helper(p3::PSP3{FT}, D_cr::FT, D_gr::FT) where {FT}
    α_va = α_va_si(p3)
    β_m2 = p3.β_va - 2
    return 6 * α_va * (D_cr^β_m2 - D_gr^β_m2) / FT(π) / β_m2 /
           max(D_cr - D_gr, eps(FT))
end

"""
    thresholds(p3, ρ_r, F_r)

 - p3 - a struct with P3 scheme parameters
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i) [-]

Solves the nonlinear system consisting of D_cr, D_gr, ρ_g, ρ_d
for a given rime density and rime mass fraction.
Returns a named tuple containing:
 - D_cr - is the threshold size separating partially rimed ice and graupel [m],
 - D_gr - is the threshold size separating graupel and dense nonspherical ice [m],
 - ρ_g - is the effective density of a spherical graupel particle [kg/m3],
 - ρ_d - is the density of the unrimed portion of the particle [kg/m3],
"""
function thresholds(p3::PSP3{FT}, ρ_r::FT, F_r::FT) where {FT}

    @assert F_r >= FT(0)   # rime mass fraction must be positive ...
    @assert F_r < FT(1)    # ... and there must always be some unrimed part

    if F_r == FT(0)
        return (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0))
    else
        @assert ρ_r > FT(0)   # rime density must be positive ...
        @assert ρ_r <= p3.ρ_l # ... and as a bulk ice density can't exceed the density of water

        P3_problem(ρ_d) =
            ρ_d - ρ_d_helper(
                p3,
                D_cr_helper(p3, F_r, ρ_g_helper(ρ_r, F_r, ρ_d)),
                D_gr_helper(p3, ρ_g_helper(ρ_r, F_r, ρ_d)),
            )

        ρ_d =
            RS.find_zero(
                P3_problem,
                RS.SecantMethod(FT(0), FT(1000)),
                RS.CompactSolution(),
            ).root
        ρ_g = ρ_g_helper(ρ_r, F_r, ρ_d)

        return (;
            D_cr = D_cr_helper(p3, F_r, ρ_g),
            D_gr = D_gr_helper(p3, ρ_g),
            ρ_g,
            ρ_d,
        )
    end
end

"""
    mass_(p3, D, ρ, F_r)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension [m]
 - ρ - bulk ice density (ρ_i for small ice, ρ_g for graupel) [kg/m3]
 - F_r - rime mass fraction [q_rim/q_i]

Returns mass as a function of size for differen particle regimes [kg]
"""
# for spherical ice (small ice or completely rimed ice)
mass_s(D::FT, ρ::FT) where {FT} = FT(π) / 6 * ρ * D^3
# for large nonspherical ice (used for unrimed and dense types)
mass_nl(p3::PSP3, D::FT) where {FT} = α_va_si(p3) * D^p3.β_va
# for partially rimed ice
mass_r(p3::PSP3, D::FT, F_r::FT) where {FT} =
    α_va_si(p3) / (1 - F_r) * D^p3.β_va
# for liquid
mass_liq(p3::PSP3, D::FT) where {FT} = (FT(π) / 6) * p3.ρ_l * D^3

"""
    p3_mass(p3, D, F_r, F_liq, th)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension
 - F_r - rime mass fraction (q_rim/q_i)
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - th - P3Scheme nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns mass(D) regime, used to create figures for the docs page.
"""
function p3_mass(
    p3::PSP3,
    D::FT,
    F_r::FT,
    F_liq::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT}
    D_th = D_th_helper(p3)
    if D_th > D
        return (1 - F_liq) * mass_s(D, p3.ρ_i) + F_liq * mass_liq(p3, D)         # small spherical ice
    elseif F_r == 0
        return (1 - F_liq) * mass_nl(p3, D) + F_liq * mass_liq(p3, D)            # large nonspherical unrimed ice
    elseif th.D_gr > D >= D_th
        return (1 - F_liq) * mass_nl(p3, D) + F_liq * mass_liq(p3, D)            # dense nonspherical ice
    elseif th.D_cr > D >= th.D_gr
        return (1 - F_liq) * mass_s(D, th.ρ_g) + F_liq * mass_liq(p3, D)         # graupel
    elseif D >= th.D_cr
        return (1 - F_liq) * mass_r(p3, D, F_r) + F_liq * mass_liq(p3, D)        # partially rimed ice
    end
end

"""
    A_(p3, D)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension

Returns particle projected area as a function of size for different particle regimes
"""
# for spherical particles
A_s(D::FT) where {FT} = FT(π) / 4 * D^2
# for nonspherical particles
A_ns(p3::PSP3, D::FT) where {FT} = p3.γ * D^p3.σ
# partially rimed ice
A_r(p3::PSP3, F_r::FT, D::FT) where {FT} =
    F_r * A_s(D) + (1 - F_r) * A_ns(p3, D)

"""
    p3_area(p3, D, F_r, F_liq, th)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension
 - F_r - rime mass fraction (q_rim/q_i)
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - th - P3Scheme nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns area(D), used to create figures for the documentation.
"""
function p3_area(
    p3::PSP3,
    D::FT,
    F_r::FT,
    F_liq::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT}
    # Area regime:
    if D_th_helper(p3) > D
        return A_s(D)                                         # small spherical ice
    elseif F_r == 0
        return (1 - F_liq) * A_ns(p3, D) + F_liq * A_s(D)        # large nonspherical unrimed ice
    elseif th.D_gr > D >= D_th_helper(p3)
        return (1 - F_liq) * A_ns(p3, D) + F_liq * A_s(D)        # dense nonspherical ice
    elseif th.D_cr > D >= th.D_gr
        return A_s(D)                                         # graupel
    elseif D >= th.D_cr
        return (1 - F_liq) * A_r(p3, F_r, D) + F_liq * A_s(D)    # partially rimed ice
    else
        throw("D not in range")
    end
end

# Some wrappers to cast types from SF.gamma
# (which returns Float64 even when the input is Float32)
Γ(a::FT, z::FT) where {FT} = FT(SF.gamma(a, z))
Γ(a::FT) where {FT} = FT(SF.gamma(a))

"""
    μ_to_λ(μ)

 - μ - parameter for gamma distribution of N′

Returns corresponding λ to given μ value
"""
function μ_to_λ(p3::PSP3, μ::FT) where {FT}
    return ((μ + p3.c) / p3.a)^(1 / p3.b)
end

"""
    DSD_μ_approx(p3, q, N, ρ_r, F_r, F_liq)

 - p3 - a struct with P3 scheme parameters
 - q - mass mixing ratio
 - N - total ice number concentration [1/m3]
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i)
 - F_liq - liquid fraction (q_liq/q_i,tot)

Returns the approximated shape parameter μ for a given q and N value
"""
function DSD_μ_approx(
    p3::PSP3,
    q::FT,
    N::FT,
    ρ_r::FT,
    F_r::FT,
    F_liq::FT,
) where {FT}
    # Get thresholds for given F_r, ρ_r
    th = thresholds(p3, ρ_r, F_r)

    # Get min and max lambda values
    λ_0 = μ_to_λ(p3, FT(0))
    λ_6 = μ_to_λ(p3, p3.μ_max)

    # Get corresponding q/N values at given F_r
    q_over_N_min = log(q_over_N_gamma(p3, F_liq, F_r, log(λ_0), FT(0), th))
    q_over_N_max = log(q_over_N_gamma(p3, F_liq, F_r, log(λ_6), p3.μ_max, th))

    # Return approximation between them
    μ = (p3.μ_max / (q_over_N_max - q_over_N_min)) * (log(q / N) - q_over_N_min)

    # Clip approximation between 0 and 6
    return min(p3.μ_max, max(FT(0), μ))
end

"""
    DSD_μ(p3, λ)

- p3 - a struct with P3 scheme parameters
- λ - slope parameter for gamma distribution of N′ [1/m]

Returns the shape parameter μ for a given λ value
Eq. 3 in Morrison and Milbrandt (2015).
"""
function DSD_μ(p3::PSP3, λ::FT) where {FT}
    #@assert λ > FT(0)
    return min(p3.μ_max, max(FT(0), p3.a * λ^p3.b - p3.c))
end

"""
    DSD_N₀(p3, μ, N, λ)

 - p3 - a struct with P3 scheme parameters
 - N - total ice number concentration [1/m3]
 - μ - shape parameter of N′ gamma distribution
 - λ - slope parameter for gamma distribution of N′ [1/m]

Returns the shape parameter N₀ from Eq. 2 in Morrison and Milbrandt (2015).
"""
function DSD_N₀(p3::PSP3, N::FT, λ::FT) where {FT}
    μ = DSD_μ(p3, λ)
    return N / Γ(1 + μ) * λ^(1 + μ)
end

"""
    ∫_Γ(x₀, x_end, c1, c2, c3)

 - x₀ - lower bound
 - x_end - upper bound
 - c1, c2, c3 - respective constants

f(D, c1, c2, c3) = c1 * D ^ (c2) * exp(-c3 * D)

Integrates f(D, c1, c2, c3) dD from x₀ to x_end
"""
function ∫_Γ(x₀::FT, x_end::FT, c1::FT, c2::FT, c3::FT) where {FT}
    if x_end == Inf
        return c1 * c3^(-c2 - 1) * (Γ(1 + c2, x₀ * c3))
    elseif x₀ == 0
        return c1 * c3^(-c2 - 1) * (Γ(1 + c2) - Γ(1 + c2, x_end * c3))
    else
        return c1 * c3^(-c2 - 1) * (Γ(1 + c2, x₀ * c3) - Γ(1 + c2, x_end * c3))
    end
end

"""
    ∫_Γ(x₀, xₘ, x_end, c1, c2, c3, c4, c5, c6)

 - x₀ - lower bound
 - xₘ - switch point
 - x_end - upper bound
 - c1, c2, c3 - respective constants for the first part of the integral
 - c4, c5, c6 - respective constants for the second part of the integral

f(D, c1, c2, c3) = c1 * D ^ (c2) * exp(-c3 * D)

Integrates f(D, c1, c2, c3) dD from x₀ to xₘ and f(D, c4, c5, c6) dD from xₘ to x_end
"""
function ∫_Γ(
    x₀::FT,
    xₘ::FT,
    x_end::FT,
    c1::FT,
    c2::FT,
    c3::FT,
    c4::FT,
    c5::FT,
    c6::FT,
) where {FT}
    return ∫_Γ(x₀, xₘ, c1, c2, c3) + ∫_Γ(xₘ, x_end, c4, c5, c6)
end

"""
    q_(p3, ρ, F_r, F_liq, λ, μ, D_min, D_max)

 - p3 - a struct with P3 scheme parameters
 - ρ - bulk ice density (ρ_i for small ice, ρ_g for graupel) [kg/m^3]
 - F_r - rime mass fraction [q_rim/q_i]
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - μ - shape parameter of N′ gamma distribution
 - λ - slope parameter of N′ gamma distribution
 - D_min - minimum bound for regime
 - D_max - maximum bound for regime (if not specified, then infinity)

 Returns ice mass density for a given m(D) regime
"""
# small, spherical ice or graupel (completely rimed, spherical)
# D_min = 0, D_max = D_th, ρ = ρᵢ
# or
# q_rim > 0 and D_min = D_gr, D_max = D_cr, ρ = ρ_g
function q_s(
    p3::PSP3,
    F_liq::FT,
    ρ::FT,
    μ::FT,
    λ::FT,
    D_min::FT,
    D_max::FT,
) where {FT}
    return ∫_Γ(D_min, D_max, (1 - F_liq) * FT(π) / 6 * ρ, μ + 3, λ)
end
# q_rim = 0 and D_min = D_th, D_max = inf
function q_rz(p3::PSP3, F_liq::FT, μ::FT, λ::FT, D_min::FT) where {FT}
    return ∫_Γ(D_min, FT(Inf), (1 - F_liq) * α_va_si(p3), μ + p3.β_va, λ)
end
# q_rim > 0 and D_min = D_th and D_max = D_gr
function q_n(p3::PSP3, F_liq::FT, μ::FT, λ::FT, D_min::FT, D_max::FT) where {FT}
    return ∫_Γ(D_min, D_max, (1 - F_liq) * α_va_si(p3), μ + p3.β_va, λ)
end
# partially rimed ice or large unrimed ice (upper bound on D is infinity)
# q_rim > 0 and D_min = D_cr, D_max = inf
function q_r(p3::PSP3, F_liq::FT, F_r::FT, μ::FT, λ::FT, D_min::FT) where {FT}
    return ∫_Γ(
        D_min,
        FT(Inf),
        (1 - F_liq) * α_va_si(p3) / (1 - F_r),
        μ + p3.β_va,
        λ,
    )
end
# F_liq != 0 (liquid mass on mixed-phase particles for D in [D_min, D_max])
function q_liq(
    p3::PSP3,
    F_liq::FT,
    μ::FT,
    λ::FT,
    D_min::FT,
    D_max::FT,
) where {FT}
    return ∫_Γ(D_min, D_max, F_liq * (FT(π) / 6) * p3.ρ_l, μ + 3, λ)
end

"""
    q_over_N_gamma(p3, F_liq, F_r, λ, th)

 - p3 - a struct with P3 scheme parameters
 - F_r - rime mass fraction [q_rim/q_i]
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - log_λ - logarithm of the slope parameter of N′ gamma distribution
 - μ - shape parameter of N′ gamma distribution
 - th - thresholds() nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns q/N for all values of D (sum over all regimes).
Eq. 5 in Morrison and Milbrandt (2015).
"""
function q_over_N_gamma(
    p3::PSP3,
    F_liq::FT,
    F_r::FT,
    log_λ::FT,
    μ::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT}

    D_th = D_th_helper(p3)
    λ = exp(log_λ)
    N = Γ(1 + μ) / (λ^(1 + μ))

    return ifelse(
        F_r == FT(0),
        (
            q_s(p3, F_liq, p3.ρ_i, μ, λ, FT(0), D_th) +
            q_rz(p3, F_liq, μ, λ, D_th) +
            q_liq(p3, F_liq, μ, λ, FT(0), FT(Inf))
        ) / N,
        (
            q_s(p3, F_liq, p3.ρ_i, μ, λ, FT(0), D_th) +
            q_n(p3, F_liq, μ, λ, D_th, th.D_gr) +
            q_s(p3, F_liq, th.ρ_g, μ, λ, th.D_gr, th.D_cr) +
            q_r(p3, F_liq, F_r, μ, λ, th.D_cr) +
            q_liq(p3, F_liq, μ, λ, FT(0), FT(Inf))
        ) / N,
    )
end

"""
    get_bounds(N, q, F_r, p3, th)

 - N - ice number concentration [1/m3]
 - q - mass mixing ratio
 - μ - shape parameter of N′ gamma distribution
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - F_r -rime mass fraction [q_rim/q_i]
 - p3 - a struct with P3 scheme parameters
 - th -  thresholds() nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

 Returns estimated guess for λ from q to be used in distribution_parameter_solver()
"""
function get_bounds(
    N::FT,
    q::FT,
    μ::FT,
    F_liq::FT,
    F_r::FT,
    p3::PSP3,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT}
    goal = q / N

    if goal >= 1e-8
        left = FT(1)
        right = FT(6 * 1e3)
        radius = FT(0.2)
    elseif goal >= 2 * 1e-9
        left = FT(6 * 1e3)
        right = FT(3 * 1e4)
        radius = FT(-0.1)
    else
        left = FT(4 * 1e4)
        right = FT(1e6)
        radius = FT(0.2)
    end

    ql = q_over_N_gamma(p3, F_liq, F_r, log(left), μ, th)
    qr = q_over_N_gamma(p3, F_liq, F_r, log(right), μ, th)

    guess =
        left * (goal / (ql))^((log(right) - log(left)) / (log(qr) - log(ql)))

    max = log(guess * exp(radius))
    min = log(guess)

    return (; min, max)
end

"""
    distrbution_parameter_solver()

 - p3 - a struct with P3 scheme parameters
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - F_r - rime mass fraction (q_rim/q_i)

Solves the nonlinear system consisting of N_0 and λ for P3 prognostic variables
Returns a named tuple containing:
 - N_0 - intercept size distribution parameter [1/m4]
 - λ - slope size distribution parameter [1/m]
"""
function distribution_parameter_solver(
    p3::PSP3{FT},
    q::FT,
    N::FT,
    ρ_r::FT,
    F_liq::FT,
    F_r::FT,
) where {FT}
    # Get the thresholds for different particles regimes
    th = thresholds(p3, ρ_r, F_r)

    # Get μ given q and N
    μ = DSD_μ_approx(p3, q, N, ρ_r, F_r, F_liq)

    # To ensure that λ is positive solve for x such that λ = exp(x)
    shape_problem(x) = q / N - q_over_N_gamma(p3, F_liq, F_r, x, μ, th)

    # Get intial guess for solver
    (; min, max) = get_bounds(N, q, μ, F_liq, F_r, p3, th)

    # Find slope parameter
    x =
        RS.find_zero(
            shape_problem,
            RS.SecantMethod(min, max),
            RS.CompactSolution(),
            RS.RelativeSolutionTolerance(eps(FT)),
            5,
        ).root

    return (; λ = exp(x), N_0 = DSD_N₀(p3, N, exp(x)))
end

"""
    terminal_velocity_(p3, Chen2022, q, N, ρ_r, F_r, ρ_a)

 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struch with terminal velocity parameters as in Chen(2022)
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i)
 - ρ_a - density of air

 Returns the mass and number weighted fall speeds
 Eq C10 of Morrison and Milbrandt (2015)
"""
function terminal_velocity(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelTypeSnowIce,
    q::FT,
    N::FT,
    ρ_r::FT,
    F_r::FT,
    ρ_a::FT,
) where {FT}
    # Get the free parameters for terminal velocities of small
    # and large particles
    small = CO.Chen2022_vel_coeffs_small(Chen2022, ρ_a)
    large = CO.Chen2022_vel_coeffs_large(Chen2022, ρ_a)
    get_p(prs, it) = (prs[1][it], prs[2][it], prs[3][it])

    # Get the thresholds for different particles regimes
    (; D_cr, D_gr, ρ_g, ρ_d) = thresholds(p3, ρ_r, F_r)
    D_th = D_th_helper(p3)
    D_ct = Chen2022.cutoff

    # Get the shape parameters of the particle size distribution
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, FT(0), F_r) # F_liq=0 to maintain terminal vel of ice part
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
    terminal_velocity_liq(p3, Chen2022, q, N, ρ_r, F_liq, F_r, ρ_a)

 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struch with terminal velocity parameters as in Chen(2022)
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - F_r - rime mass fraction (q_rim/q_i)
 - ρ_a - density of air

 Returns the mass and number weighted fall speeds
 Eq C10 of Morrison and Milbrandt (2015)
 for the liquid part of mixed-phase particles
 (For now using aspect ratio=1, κ=0 as in Microphysics2M,
 assuming spherical rain drop size...
 later we can change this, especially if we don't like how the
 F_liq-weighted terminal velocity comes out)
"""
function terminal_velocity_liq(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelTypeRain,
    q::FT,
    N::FT,
    ρ_r::FT,
    F_liq::FT,
    F_r::FT,
    ρ_a::FT,
) where {FT}
    # Get the pree parameters for terminal velocities of small particles
    # (No large particles for rain?)
    small = CO.Chen2022_vel_coeffs_small(Chen2022, ρ_a)
    get_p(prs, it) = (prs[1][it], prs[2][it], prs[3][it])

    # Get the shape parameters of the particle size distribution
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, F_liq, F_r)
    μ = DSD_μ(p3, λ)
    κ = 0

    aₛ(a) = a * N_0
    bₛ(b) = b + μ
    cₛ(c) = c + λ

    aₛ_m(a) = aₛ(a) * FT(π) / 6 * p3.ρ_l
    bₛ_m(b) = bₛ(b) + 3

    spheres_n(a, b, c) = (aₛ(a), bₛ(b), cₛ(c))
    spheres_m(a, b, c) = (aₛ_m(a), bₛ_m(b), cₛ(c))

    v_m = 0
    v_n = 0

    v(D, a, b, c) = a * D^b * exp(-c * D)

    for i in 1:3
        # TODO: fix bounds for quadgk integral
        (Im, em) =
            QGK.quadgk(D -> v(D, spheres_m(get_p(small, i)...)...), FT(0), Inf)
        (In, en) =
            QGK.quadgk(D -> v(D, spheres_m(get_p(small, i)...)...), FT(0), Inf)
        v_m += Im
        v_n += In
    end
    return (v_n / N, v_m / q)
end

"""
    terminal_velocity_tot(p3, Chen2022, q, N, ρ_r, F_liq, F_r, ρ_a)

 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struct with terminal velocity parameters as in Chen(2022)
            - pass (Chen2022VelTypeSnowIce, Chen2022VelTypeRain)
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - F_r - rime mass fraction (q_rim/q_i)
 - ρ_a - density of air

 Returns the mass and number weighted fall speeds
 Eq C10 of Morrison and Milbrandt (2015)
 for a mixed-phase particle
"""
function terminal_velocity_tot(
    p3::PSP3,
    Chen2022_ice::CMP.Chen2022VelTypeSnowIce,
    Chen2022_rain::CMP.Chen2022VelTypeRain,
    q::FT,
    N::FT,
    ρ_r::FT,
    F_liq::FT,
    F_r::FT,
    ρ_a::FT,
) where {FT}
    return (1 - F_liq) *
           terminal_velocity(p3, Chen2022_ice, q, N, ρ_r, F_r, ρ_a)[1] +
           F_liq * terminal_velocity_liq(
        p3,
        Chen2022_rain,
        q,
        N,
        ρ_r,
        F_liq,
        F_r,
        ρ_a,
    )[1],
    (1 - F_liq) * terminal_velocity(p3, Chen2022_ice, q, N, ρ_r, F_r, ρ_a)[2] +
    F_liq *
    terminal_velocity_liq(p3, Chen2022_rain, q, N, ρ_r, F_liq, F_r, ρ_a)[2]
end

"""
    D_m (p3, q, N, ρ_r, F_r, F_liq)

 - p3 - a struct with P3 scheme parameters
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - F_r - rime mass fraction (q_rim/q_i)

 Return the mass weighted mean particle size [m]
"""
function D_m(p3::PSP3, q::FT, N::FT, ρ_r::FT, F_r::FT, F_liq::FT) where {FT}
    # Get the thresholds for different particles regimes
    th = thresholds(p3, ρ_r, F_r)
    D_th = D_th_helper(p3)

    # Get the shape parameters
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, F_liq, F_r)
    μ = DSD_μ(p3, λ)

    # Redefine α_va to be in si units
    α_va = α_va_si(p3)

    # Calculate numerator
    n = 0
    if F_r == 0
        n += ∫_Γ(FT(0), D_th, π / 6 * p3.ρ_i * N_0, μ + 4, λ)
        n += ∫_Γ(D_th, Inf, α_va * N_0, μ + p3.β_va + 1, λ)
    else
        n += ∫_Γ(FT(0), D_th, π / 6 * p3.ρ_i * N_0, μ + 4, λ)
        n += ∫_Γ(D_th, th.D_gr, α_va * N_0, μ + p3.β_va + 1, λ)
        n += ∫_Γ(th.D_gr, th.D_cr, π / 6 * th.ρ_g * N_0, μ + 4, λ)
        n += ∫_Γ(th.D_cr, Inf, α_va / (1 - F_r) * N_0, μ + p3.β_va + 1, λ)
    end
    # Normalize by q
    return n / q
end

"""
    D_m_liq (p3, q, N, ρ_r, F_r, F_liq)

 - p3 - a struct with P3 scheme parameters
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - F_r - rime mass fraction (q_rim/q_i)

 Return the mass weighted mean particle size [m]
"""
function D_m_liq(p3::PSP3, q::FT, N::FT, ρ_r::FT, F_r::FT, F_liq::FT) where {FT}
    # Get the shape parameters
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, F_liq, F_r)
    μ = DSD_μ(p3, λ)
    n(D) = F_liq * N_0 * p3.ρ_l * (FT(π) / 6) * D^(μ + 4) * exp(-λ * D)
    (n, em) = QGK.quadgk(D -> n(D), FT(0), Inf)
    # Normalize by q
    return n / q
end

"""
    D_m_tot(p3, q, N, ρ_r, F_r, F_liq)

 - p3 - a struct with P3 scheme parameters
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - F_r - rime mass fraction (q_rim/q_i)

 Return the mass weighted mean particle size [m]
"""
function D_m_tot(p3::PSP3, q::FT, N::FT, ρ_r::FT, F_r::FT, F_liq::FT) where {FT}
    return (1 - F_liq) * D_m(p3, q, N, ρ_r, F_r, F_liq) +
           F_liq * D_m_liq(p3, q, N, ρ_r, F_r, F_liq)
end


"""
    N′ice(p3, D, λ, N0)

 - p3 - a struct containing P3 scheme parameters
 - D - diameter of particle 
 - λ - shape parameter of distribution 
 - N0 - shape parameter of distribution

 Returns the distribution of ice particles (assumed to be of the form 
 N'(D) = N0 * D ^ μ * exp(-λD)) at given D 
"""
function N′ice(p3::PSP3, D::FT, λ::FT, N_0::FT) where {FT}
    return N_0 * D^DSD_μ(p3, λ) * exp(-λ * D)
end

"""
    get_ice_bounds(p3, λ, tolerance)

 - p3 - a struct containing p3 parameters 
 - λ - shape parameters of ice distribution 
 - tolerance - tolerance for how much of the distribution we want to integrate over

 Returns the bound on the distribution that would guarantee that 1-tolerance 
 of the ice distribution is integrated over. This is calculated by setting 
 N_0(1 - tolerance) = ∫ N'(D) dD from 0 to bound and solving for bound. 
 This was further simplified to cancel out the N_0 from both sides. 
 The guess was calculated through a linear approximation extrapolated from 
 numerical solutions. 
"""
function get_ice_bound(p3, λ::FT, tolerance::FT) where {FT}
    ice_problem(x) =
        tolerance - Γ(1 + DSD_μ(p3, λ), FT(exp(x)) * λ) / Γ(1 + DSD_μ(p3, λ))
    guess = log(19 / 6 * (DSD_μ(p3, λ) - 1) + 39) - log(λ)
    log_ice_x =
        RS.find_zero(
            ice_problem,
            RS.SecantMethod(guess - 1, guess),
            RS.CompactSolution(),
            RS.RelativeSolutionTolerance(eps(FT)),
            5,
        ).root

    return exp(log_ice_x)
end

end
