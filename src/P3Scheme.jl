"""
Predicted particle properties scheme (P3) for ice, which includes:
 - threshold solver
 - shape parameters solver
 - m(D) regime
 - a(D) regime

Implementation of Morrison and Milbrandt 2015 doi: 10.1175/JAS-D-14-0065.1

Note: Particle size is defined as its maximum length (i.e. max dimesion).
"""
module P3Scheme

import SpecialFunctions as SF
import QuadGK as QGK
import RootSolvers as RS
import HCubature as HC

import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics as CM
import CloudMicrophysics.TerminalVelocity as TV
import CloudMicrophysics.Microphysics2M as CM2

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
mass_s(D::FT, ρ::FT) where {FT <: Real} = FT(π) / 6 * ρ * D^3
# for large nonspherical ice (used for unrimed and dense types)
mass_nl(p3::PSP3, D::FT) where {FT <: Real} = α_va_si(p3) * D^p3.β_va
# for partially rimed ice
mass_r(p3::PSP3, D::FT, F_r::FT) where {FT <: Real} =
    α_va_si(p3) / (1 - F_r) * D^p3.β_va

"""
    p3_mass(p3, D, F_r, th)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension
 - F_r - rime mass fraction (q_rim/q_i)
 - th - P3Scheme nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns mass(D) regime, used to create figures for the docs page.
"""
function p3_mass(
    p3::PSP3,
    D::FT,
    F_r::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT <: Real}
    D_th = D_th_helper(p3)
    if D_th > D
        return mass_s(D, p3.ρ_i)          # small spherical ice
    elseif F_r == 0
        return mass_nl(p3, D)             # large nonspherical unrimed ice
    elseif th.D_gr > D >= D_th
        return mass_nl(p3, D)             # dense nonspherical ice
    elseif th.D_cr > D >= th.D_gr
        return mass_s(D, th.ρ_g)          # graupel
    elseif D >= th.D_cr
        return mass_r(p3, D, F_r)         # partially rimed ice
    end
end

"""
    A_(p3, D)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension

Returns particle projected area as a function of size for different particle regimes
"""
# for spherical particles
A_s(D::FT) where {FT <: Real} = FT(π) / 4 * D^2
# for nonspherical particles
A_ns(p3::PSP3, D::FT) where {FT <: Real} = p3.γ * D^p3.σ
# partially rimed ice
A_r(p3::PSP3, F_r::FT, D::FT) where {FT <: Real} =
    F_r * A_s(D) + (1 - F_r) * A_ns(p3, D)

"""
    p3_area(p3, D, F_r, th)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension
 - F_r - rime mass fraction (q_rim/q_i)
 - th - P3Scheme nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns area(D), used to create figures for the documentation.
"""
function p3_area(
    p3::PSP3,
    D::FT,
    F_r::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT <: Real}
    # Area regime:
    if D_th_helper(p3) > D
        return A_s(D)                      # small spherical ice
    elseif F_r == 0
        return A_ns(p3, D)                 # large nonspherical unrimed ice
    elseif th.D_gr > D >= D_th_helper(p3)
        return A_ns(p3, D)                 # dense nonspherical ice
    elseif th.D_cr > D >= th.D_gr
        return A_s(D)                      # graupel
    elseif D >= th.D_cr
        return A_r(p3, F_r, D)             # partially rimed ice
    else
        throw("D not in range")
    end
end

# Some wrappers to cast types from SF.gamma
# (which returns Float64 even when the input is Float32)
Γ(a::FT, z::FT) where {FT <: Real} = FT(SF.gamma(a, z))
Γ(a::FT) where {FT <: Real} = FT(SF.gamma(a))

"""
    μ_to_λ(μ)

 - μ - parameter for gamma distribution of N′

Returns corresponding λ to given μ value
"""
function μ_to_λ(p3::PSP3, μ::FT) where {FT}
    return ((μ + p3.c) / p3.a)^(1 / p3.b)
end

"""
    DSD_μ_approx(p3, q, N, ρ_r, F_r)

 - p3 - a struct with P3 scheme parameters
 - q - mass mixing ratio
 - N - total ice number concentration [1/m3]
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i)

Returns the approximated shape parameter μ for a given q and N value
"""
function DSD_μ_approx(p3::PSP3, q::FT, N::FT, ρ_r::FT, F_r::FT) where {FT}
    # Get thresholds for given F_r, ρ_r
    th = thresholds(p3, ρ_r, F_r)

    # Get min and max lambda values
    λ_0 = μ_to_λ(p3, FT(0))
    λ_6 = μ_to_λ(p3, p3.μ_max)

    # Get corresponding q/N values at given F_r
    q_over_N_min = log(q_over_N_gamma(p3, F_r, log(λ_0), FT(0), th))
    q_over_N_max = log(q_over_N_gamma(p3, F_r, log(λ_6), p3.μ_max, th))

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
    q_(p3, ρ, F_r, λ, μ, D_min, D_max)

 - p3 - a struct with P3 scheme parameters
 - ρ - bulk ice density (ρ_i for small ice, ρ_g for graupel) [kg/m^3]
 - F_r - rime mass fraction [q_rim/q_i]
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
function q_s(p3::PSP3, ρ::FT, μ::FT, λ::FT, D_min::FT, D_max::FT) where {FT}
    return ∫_Γ(D_min, D_max, FT(π) / 6 * ρ, μ + 3, λ)
end
# q_rim = 0 and D_min = D_th, D_max = inf
function q_rz(p3::PSP3, μ::FT, λ::FT, D_min::FT) where {FT}
    return ∫_Γ(D_min, FT(Inf), α_va_si(p3), μ + p3.β_va, λ)
end
# q_rim > 0 and D_min = D_th and D_max = D_gr
function q_n(p3::PSP3, μ::FT, λ::FT, D_min::FT, D_max::FT) where {FT}
    return ∫_Γ(D_min, D_max, α_va_si(p3), μ + p3.β_va, λ)
end
# partially rimed ice or large unrimed ice (upper bound on D is infinity)
# q_rim > 0 and D_min = D_cr, D_max = inf
function q_r(p3::PSP3, F_r::FT, μ::FT, λ::FT, D_min::FT) where {FT}
    return ∫_Γ(D_min, FT(Inf), α_va_si(p3) / (1 - F_r), μ + p3.β_va, λ)
end

"""
    q_over_N_gamma(p3, F_r, λ, th)

 - p3 - a struct with P3 scheme parameters
 - F_r - rime mass fraction [q_rim/q_i]
 - log_λ - logarithm of the slope parameter of N′ gamma distribution
 - μ - shape parameter of N′ gamma distribution
 - th - thresholds() nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns q/N for all values of D (sum over all regimes).
Eq. 5 in Morrison and Milbrandt (2015).
"""
function q_over_N_gamma(
    p3::PSP3,
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
        (q_s(p3, p3.ρ_i, μ, λ, FT(0), D_th) + q_rz(p3, μ, λ, D_th)) / N,
        (
            q_s(p3, p3.ρ_i, μ, λ, FT(0), D_th) +
            q_n(p3, μ, λ, D_th, th.D_gr) +
            q_s(p3, th.ρ_g, μ, λ, th.D_gr, th.D_cr) +
            q_r(p3, F_r, μ, λ, th.D_cr)
        ) / N,
    )
end

"""
    get_bounds(N, q, F_r, p3, th)

 - N - ice number concentration [1/m3]
 - q - mass mixing ratio
 - μ - shape parameter of N′ gamma distribution
 - F_r -rime mass fraction [q_rim/q_i]
 - p3 - a struct with P3 scheme parameters
 - th -  thresholds() nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

 Returns estimated guess for λ from q to be used in distribution_parameter_solver()
"""
function get_bounds(
    N::FT,
    q::FT,
    μ::FT,
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

    ql = q_over_N_gamma(p3, F_r, log(left), μ, th)
    qr = q_over_N_gamma(p3, F_r, log(right), μ, th)

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
    F_r::FT,
) where {FT}
    # Get the thresholds for different particles regimes
    th = thresholds(p3, ρ_r, F_r)

    # Get μ given q and N
    μ = DSD_μ_approx(p3, q, N, ρ_r, F_r)

    # To ensure that λ is positive solve for x such that λ = exp(x)
    shape_problem(x) = q / N - q_over_N_gamma(p3, F_r, x, μ, th)

    # Get intial guess for solver
    (; min, max) = get_bounds(N, q, μ, F_r, p3, th)

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
    terminal_velocity(p3, Chen2022, q, N, ρ_r, F_r, ρ_a)

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
    D_m (p3, q, N, ρ_r, F_r)

 - p3 - a struct with P3 scheme parameters
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i)

 Return the mass weighted mean particle size [m]
"""
function D_m(p3::PSP3, q::FT, N::FT, ρ_r::FT, F_r::FT) where {FT}
    # Get the thresholds for different particles regimes
    th = thresholds(p3, ρ_r, F_r)
    D_th = D_th_helper(p3)

    # Get the shape parameters
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, F_r)
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
    N′_ice(D, p3, λ, N0)

 - D - diameter of particle 
 - p3 - a struct containing P3 scheme parameters
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

"""
    integration_bounds(p3, tolerance, λ, N_0, Nᵢ, qᵣ, Nᵣ)
    
    - p3 - a struct containing P3 Scheme parameters 
    - tolerance - tolerance to which distributions need to be evaluated to
    - λ - shape parameter of ice distribution
    - N_0 - intercept size distribution of ice 
    - Nᵢ - number mixing ratio of ice
    - q - mass mizing ratio of colliding species
    - N - number mixing ratio of colliding species
   
    Returns the bounds over which to integrate rain, cloud, and ice distributions to ensure 
    coverage or more than (1 - tolerance) of each distribution
   """
function integration_bounds(
    pdf::Union{
        CMP.RainParticlePDF_SB2006{FT},
        CMP.RainParticlePDF_SB2006_limited{FT},
        CMP.CloudParticlePDF_SB2006{FT},
    },
    p3::PSP3,
    tolerance::FT,
    λ::FT,
    q::FT,
    N::FT,
    ρ_a::FT,
) where {FT}
    colliding_x = CM2.get_distribution_bound(pdf, q, N, ρ_a, tolerance)
    ice_bound = get_ice_bound(p3, λ, tolerance)
    return (2 * colliding_x, 2 * ice_bound)
end

"""
    a(D1, D2) 

 - D1 - maximum dimension of first particle 
 - D2 - maximum dimension of second particle 

 Returns the collision kernel (assumed to be of the form π(r1 + r2)^2) for the 
 two colliding particles
"""
function a(D1::FT, D2::FT) where {FT}
    # TODO make this more accurate for non-spherical particles 
    return π * (D1 / 2 + D2 / 2)^2
end

"""
    vel_diff(type, D, Dᵢ, p3, Chen2022, ρ_a, F_r, th)

 - type - defines what is colliding with ice ("rain" or "cloud")
 - D - maximum dimension of colliding particle 
 - Dᵢ - maximum dimension of ice particle 
 - p3 - a struct containing P3 parameters
 - Chen2022 - a struct containing Chen 2022 velocity parameters 
 - ρ_a - density of air 
 - F_r - rime mass fraction (q_rim/ q_i)
 - th - thresholds as calculated by thresholds()

 Returns the corresponding velocity difference of colliding particles depending on type
"""
function vel_diff(
    pdf_r::Union{
        CMP.RainParticlePDF_SB2006{FT},
        CMP.RainParticlePDF_SB2006_limited{FT},
    },
    D::FT,
    Dᵢ::FT,
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    ρ_a::FT,
    F_r::FT,
    th,
) where {FT}
    return abs(
        TV.velocity_chen(
            Dᵢ,
            Chen2022.snow_ice,
            ρ_a,
            p3_mass(p3, D, F_r, th),
            p3_area(p3, D, F_r, th),
            p3.ρ_i,
        ) - TV.velocity_chen(D, Chen2022.rain, ρ_a),
    )
end

# velocity difference for cloud collisions
function vel_diff(
    pdf_c::CMP.CloudParticlePDF_SB2006{FT},
    D::FT,
    Dᵢ::FT,
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    ρ_a::FT,
    F_r::FT,
    th,
) where {FT}
    return abs(
        TV.velocity_chen(
            Dᵢ,
            Chen2022.snow_ice,
            ρ_a,
            p3_mass(p3, D, F_r, th),
            p3_area(p3, D, F_r, th),
            p3.ρ_i,
        ),
    )
end

"""
    ice_collisions(type, p3, Chen2022, ρ_a, F_r, qᵣ, qᵢ, Nᵣ, Nᵢ, ρ, ρ_r, E_ri)

 - type - defines what is colliding with the ice ("cloud" or "rain")
 - p3 - a struct with P3 scheme parameters 
 - Chen2022 - a struct with terminal velocity parameters as in Chen (2022) 
 - qᵢ - mass mixing ratio of ice 
 - Nᵢ - number mixing ratio of ice 
 - q_c - mass mixing ratio of colliding species 
 - N_c - number mixing ratio of colliding species
 - ρ_a - density of air 
 - F_r - rime mass fraction (q_rim/ q_i) 
 - ρ_r - rime density (q_rim/B_rim) 
 - T - temperature (in K)
 - E_ci - collision efficiency between ice and colliding species 

 Returns the rate of collisions between cloud ice and rain 
 Equivalent to the measure of QRCOL in Morrison and Mildbrandt (2015)
"""
function ice_collisions(
    pdf::Union{
        CMP.RainParticlePDF_SB2006{FT},
        CMP.RainParticlePDF_SB2006_limited{FT},
        CMP.CloudParticlePDF_SB2006{FT},
    },
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    qᵢ::FT,
    Nᵢ::FT,
    q_c::FT,
    N_c::FT,
    ρ_a::FT,
    F_r::FT,
    ρ_r::FT,
    T::FT,
    E_ci = FT(1),
) where {FT}
    ρ_l = p3.ρ_l
    th = thresholds(p3, ρ_r, F_r)
    (λ, N_0) = distribution_parameter_solver(p3, qᵢ, Nᵢ, ρ_r, F_r)
    (colliding_bound, ice_bound) =
        integration_bounds(pdf, p3, eps(FT), λ, q_c * ρ_a, N_c, ρ_a)

    if T > p3.T_freeze
        f_warm(D_c, Dᵢ) =
            E_ci / ρ_a *
            CM2.particle_size_distribution(pdf, D_c, q_c, ρ_a, N_c) *
            N′ice(p3, Dᵢ, λ, N_0) *
            a(D_c, Dᵢ) *
            p3_mass(p3, Dᵢ, F_r, th) *
            vel_diff(pdf, D_c, Dᵢ, p3, Chen2022, ρ_a, F_r, th)
        (dqdt, error) = HC.hcubature(
            d -> f_warm(d[1], d[2]),
            (FT(0), FT(0)),
            (colliding_bound, ice_bound),
        )
    else
        f_cold(D_c, Dᵢ) =
            E_ci / ρ_a *
            CM2.particle_size_distribution(pdf, D_c, q_c, ρ_a, N_c) *
            N′ice(p3, Dᵢ, λ, N_0) *
            a(D_c, Dᵢ) *
            mass_s(D_c, ρ_l) *
            vel_diff(pdf, D_c, Dᵢ, p3, Chen2022, ρ_a, F_r, th)
        (dqdt, error) = HC.hcubature(
            d -> f_cold(d[1], d[2]),
            (FT(0), FT(0)),
            (colliding_bound, ice_bound),
        )
    end

    return dqdt
end

"""
    dmdt_mass(p3, D, F_r, th) 

 - p3 - a struct containing p3 parameters 
 - D - maximum dimension of the particle 
 - F_r - rime mass fraction (q_rim/ q_i) 
 - th - thresholds as calculated by thresholds()

Returns the value equivalent to dm(D)/dt * 4 / D for each P3 regime
    4 / D comes from dD/dt
"""
function dmdt_mass(p3, D, F_r, th)
    D_th = D_th_helper(p3)
    if D_th > D
        return 2 * π * p3.ρ_i * D
    elseif F_r == 0
        return 4 * α_va_si(p3) * p3.β_va * D^(p3.β_va - 2)
    elseif th.D_gr > D >= D_th
        return 4 * α_va_si(p3) * p3.β_va * D^(p3.β_va - 2)
    elseif th.D_cr > D >= th.D_gr
        return 2 * π * th.ρ_g * D
    elseif D >= th.D_cr
        return 4 * α_va_si(p3) / (1 - F_r) * p3.β_va * D^(p3.β_va - 2)
    end
end

"""
    p3_melt(p3, Chen2022, aps, tps, q, N, ρ, T, ρ_a, F_r, ρ_r)

 - p3 - a struct containing p3 parameters 
 - Chen2022 - struct containing Chen 2022 velocity parameters 
 - aps - air properties
 - tps - thermodynamics parameters
 - q - mass mixing ratio of ice 
 - N - number mixing ratio of ice 
 - T - temperature (K) 
 - ρ_a - air density
 - F_r - rime mass fraction (q_rim/ q_i)
 - ρ_r - rime density (q_rim/B_rim) 

 Returns the calculated melting rate of ice
 Equivalent to the measure of QIMLT in Morrison and Mildbrandt (2015) 
"""
function p3_melt(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q::FT,
    N::FT,
    T::FT,
    ρ_a::FT,
    F_r::FT,
    ρ_r::FT,
) where {FT}
    # Get constants
    (; ν_air, D_vapor, K_therm) = aps
    a = p3.vent_a
    b = p3.vent_b
    L_f = TD.latent_heat_fusion(tps, T)
    T_freeze = p3.T_freeze
    N_sc = ν_air / D_vapor
    ρ_l = p3.ρ_l

    # Get distribution values  
    th = thresholds(p3, ρ_r, F_r)
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, F_r)

    # Get bound 
    ice_bound = get_ice_bound(p3, λ, eps(FT))

    # Define function pieces
    N_re(D) =           # TODO: What is this for non-spherical particles?
        D * TV.velocity_chen(
            D,
            Chen2022.snow_ice,
            ρ_a,
            p3_mass(p3, D, F_r, th),
            p3_area(p3, D, F_r, th),
            p3.ρ_i,
        ) / ν_air
    F(D) = a + b * N_sc^(1 / 3) * N_re(D)^(1 / 2)
    dmdt(D) =
        dmdt_mass(p3, D, F_r, th) / ρ_l * K_therm / L_f * (T - T_freeze) * F(D)
    f(D) = 1 / (2 * ρ_a) * dmdt(D) * N′ice(p3, D, λ, N_0)

    # Integrate 
    (dqdt, error) = QGK.quadgk(d -> f(d), FT(0), 2 * ice_bound)

    return dqdt
end

"""
    p3_het_freezing(mass, tps, q, N, T, ρ_a, qᵥ, aero_type)

 - mass - true if calculating change in mass, false for change in number
 - tps - thermodynamics parameters
 - q - mass mixing ratio of rain
 - N - number mixing ratio of rain
 - T - temperature in K 
 - ρ_a - density of air 
 - qᵥ - mixing ratio of water vapor 
 - aero_type - type of aerosols present 

 Returns the rate of hetergoeneous freezing within rain or cloud water 
    If mass false corresponds to NRHET in Morrison and Mildbrandt (2015)
    If mass true corresponds to QRHET in Morrison and Mildbrandt (2015)
"""
function p3_rain_het_freezing(
    mass::Bool,
    pdf_r::CMP.RainParticlePDF_SB2006{FT},
    p3::PSP3,
    tps::TDP.ThermodynamicsParameters{FT},
    q::FT,
    N::FT,
    T::FT,
    ρ_a::FT,
    qᵥ::FT,
    aero_type,
) where {FT}
    ρ_w = p3.ρ_l

    Rₐ = TD.gas_constant_air(tps, TD.PhasePartition(qᵥ))
    R_v = TD.Parameters.R_v(tps)
    e = qᵥ * ρ_a * R_v / Rₐ

    a_w = CO.a_w_eT(tps, e, T)
    a_w_ice = CO.a_w_ice(tps, T)
    Δa_w = a_w - a_w_ice
    J_immersion = CM.HetIceNucleation.ABIFM_J(aero_type, Δa_w)

    bound = CM2.get_distribution_bound(pdf_r, q, N, ρ_a, eps(FT))
    if mass
        f_rain_mass(D) =
            J_immersion *
            CM2.particle_size_distribution(pdf_r, D, q * ρ_a, ρ_a, N) *
            mass_s(D, ρ_w) *
            π *
            D^2
        dqdt, = QGK.quadgk(d -> f_rain_mass(d), FT(0), 2 * bound)
        return dqdt
    else
        f_rain(D) =
            J_immersion *
            CM2.particle_size_distribution(pdf_r, D, q * ρ_a, ρ_a, N) *
            π *
            D^2
        dNdt, = QGK.quadgk(d -> f_rain(d), FT(0), 2 * bound)
        return dNdt
    end
end

end
