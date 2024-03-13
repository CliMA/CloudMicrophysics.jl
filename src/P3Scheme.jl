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

import RootSolvers as RS
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO

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
    integrate(a, b, c1, c2, c3)

 - a - lower bound 
 - b - upper bound 
 - c1, c2, c3 - respective constants 

 Integrates the function c1 * D ^ (c2) * exp(-c3 * D) dD from a to b 
    Returns the result
"""
function integrate(a::FT, b::FT, c1::FT, c2::FT, c3::FT) where{FT}
    if b == Inf
        return c1 * c3 ^ (-c2 - 1) * (Γ(1 + c2, a * c3))
    elseif a == 0 
        return c1 * c3 ^ (-c2 - 1) * (Γ(1 + c2) - Γ(1 + c2, b * c3))
    else 
        return c1 * c3 ^ (-c2 - 1) * (Γ(1 + c2, a * c3) - Γ(1 + c2, b * c3))
    end
end

"""
    get_coeffs(p3, th) 

 - p3 - a struct with P3 scheme parameters
 - th - thresholds tuple as returned by thresholds()
 - F_r - rime mass fraction [q_rim/q_i]

Returns the coefficients for m(D), a(D), and the respective powers of D 
    Where the indices are as follows: 
        1 - small, spherical ice 
        2 - large, unrimed ice 
        3 - dense, nonspherical ice 
        4 - graupel 
        5 - partially rimed ice 
        6 - second half of partially rimed ice (only for a)
"""
function get_coeffs(p3::PSP3, th, F_r::FT) where{FT}
    α_va = α_va_si(p3)
    m = [π / 6 * p3.ρ_i, α_va, α_va, π / 6 * th.ρ_g, α_va / (1- F_r)]
    m_power = [FT(3), p3.β_va, p3.β_va, 3, p3.β_va]
    a = [π/4, p3.γ, p3.γ, π/4, F_r * π/4, (1 - F_r) * p3.γ] 
    a_power = [FT(2), p3.σ, p3.σ, FT(2), FT(2), p3.σ]
    return (; m, m_power, a, a_power)
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
    return integrate(D_min, D_max, FT(π) / 6 * ρ, μ + 3, λ)
end
# q_rim = 0 and D_min = D_th, D_max = inf
function q_rz(p3::PSP3, μ::FT, λ::FT, D_min::FT) where {FT}
    return integrate(D_min, FT(Inf), α_va_si(p3), μ + p3.β_va, λ)
end
# q_rim > 0 and D_min = D_th and D_max = D_gr
function q_n(p3::PSP3, μ::FT, λ::FT, D_min::FT, D_max::FT) where {FT}
    return integrate(D_min, D_max, α_va_si(p3), μ + p3.β_va, λ)
end
# partially rimed ice or large unrimed ice (upper bound on D is infinity)
# q_rim > 0 and D_min = D_cr, D_max = inf
function q_r(p3::PSP3, F_r::FT, μ::FT, λ::FT, D_min::FT) where {FT}
    return integrate(D_min, FT(Inf), α_va_si(p3) / (1 - F_r), μ + p3.β_va, λ)
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
    ϕ_coeff(m, a) 

 - m - coefficient of mass
 - a - coefficient of mass 

 Returns the coefficient of aspect ratio (ignoring powers of D)
"""
function ϕ_coeff(p3::PSP3, m::FT, a::FT) where{FT}
    return 16 * p3.ρ_i ^ 2 * a ^ 3 / (9 * π^2 * m^2)
end

"""
    terminal_velocity_mass(p3, Chen2022, q, N, ρ_r, F_r)

 - p3 - a struct with P3 scheme parameters
 - Chen 2022 - a struch with terminal velocity parameters as in Chen(2022)
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i)
 - ρ_a - density of air 

 Returns the mass (total)-weighted fall speed
 Eq C10 of Morrison and Milbrandt (2015)
"""
function terminal_velocity_mass(p3::PSP3, Chen2022::CMP.Chen2022VelTypeSnowIce, q::FT, N::FT, ρ_r::FT, F_r::FT, ρ_a::FT) where{FT}

    # Get the thresholds for different particles regimes
    th = thresholds(p3, ρ_r, F_r)
    D_th = D_th_helper(p3)
    cutoff = FT(0.000625) # TO be added to the struct

    # Get the shape parameters
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, F_r)
    μ = DSD_μ(p3, λ)

    # Get the ai, bi, ci constants (in si units) for velocity calculations
    (ai, bi, ci) = CO.Chen2022_vel_coeffs_small(Chen2022, ρ_a)

    κ = FT(-1/6) #FT(1/3)

    # Redefine α_va to be in si units
    α_va = α_va_si(p3)

    # TODO Update the velocity to use a different formulation for D > 0.625 mm 
    v = 0
    for i in 1:2
        if F_r == 0 
            v += integrate(FT(0), D_th, π / 6 * p3.ρ_i * ai[i] * N_0, bi[i] + μ + 3, ci[i] + λ)
            v += integrate(D_th, cutoff, α_va * ai[i] * N_0 * (16 * p3.ρ_i ^ 2 * p3.γ^3/(9 * π * α_va ^ 2)) ^ κ, bi[i] + μ + p3.β_va + κ * (3 * p3.σ - 2 * p3.β_va), ci[i] + λ)
            
            # Get velocity coefficients for large particles
            (ai, bi, ci) = CO.Chen2022_vel_coeffs_large(Chen2022, ρ_a)
            v += integrate(cutoff, Inf, α_va * ai[i] * N_0 * (16 * p3.ρ_i ^ 2 * p3.γ^3/(9 * π * α_va ^ 2)) ^ κ, bi[i] + μ + p3.β_va + κ * (3 * p3.σ - 2 * p3.β_va), ci[i] + λ)
        else 
            large = false
            v += integrate(FT(0), D_th, π / 6 * p3.ρ_i * ai[i] * N_0, bi[i] + μ + 3, ci[i] + λ)
            
            # D_th to D_gr
            if !large && th.D_gr > cutoff 
                v += integrate(D_th, cutoff, α_va * ai[i] * N_0 * (16 * p3.ρ_i ^ 2 * p3.γ^3/(9 * π * α_va ^ 2)) ^ κ, bi[i] + μ + p3.β_va + κ * (3 * p3.σ - 2 * p3.β_va), ci[i] + λ)
                
                # large particles 
                (ai, bi, ci) = CO.Chen2022_vel_coeffs_large(Chen2022, ρ_a)
                large = true 

                v += integrate(cutoff, th.D_gr, α_va * ai[i] * N_0 * (16 * p3.ρ_i ^ 2 * p3.γ^3/(9 * π * α_va ^ 2)) ^ κ, bi[i] + μ + p3.β_va + κ * (3 * p3.σ - 2 * p3.β_va), ci[i] + λ)
            else
                v += integrate(D_th, th.D_gr, α_va * ai[i] * N_0 * (16 * p3.ρ_i ^ 2 * p3.γ^3/(9 * π * α_va ^ 2)) ^ κ, bi[i] + μ + p3.β_va + κ * (3 * p3.σ - 2 * p3.β_va), ci[i] + λ)
            end

            # D_gr to D_cr
            if !large && th.D_cr > cutoff 
                v += integrate(th.D_gr, cutoff, π / 6 * th.ρ_g * ai[i] * N_0 * (p3.ρ_i / th.ρ_g)^(2 * κ), bi[i] + μ + 3, ci[i] + λ)
                    
                # large particles 
                (ai, bi, ci) = CO.Chen2022_vel_coeffs_large(Chen2022, ρ_a)
                large = true 

                v += integrate(cutoff, th.D_cr, π / 6 * th.ρ_g * ai[i] * N_0 * (p3.ρ_i / th.ρ_g)^(2 * κ), bi[i] + μ + 3, ci[i] + λ)
            else 
                v += integrate(th.D_gr, th.D_cr, π / 6 * th.ρ_g * ai[i] * N_0 * (p3.ρ_i / th.ρ_g)^(2 * κ), bi[i] + μ + 3, ci[i] + λ)
            end
            
            # D_cr to Infinity
            if !large
                # approximating sigma as 2 (for closed form integration) (this overestimates A LOT)
                v += integrate(th.D_cr, cutoff, 1/(1 - F_r) * α_va * ai[i] * N_0 * (16 * p3.ρ_i ^ 2 * (F_r * π / 4 + (1 - F_r) * p3.γ)^3 / (9 * π * (1/(1 - F_r) * α_va) ^ 2)) ^ (κ), bi[i] + μ + p3.β_va + κ*(6 - 2 * p3.β_va), ci[i] + λ)
                    
                # large particles 
                (ai, bi, ci) = CO.Chen2022_vel_coeffs_large(Chen2022, ρ_a)
                large = true

                v += integrate(cutoff, Inf, 1/(1 - F_r) * α_va * ai[i] * N_0 * (16 * p3.ρ_i ^ 2 * (F_r * π / 4 + (1 - F_r) * p3.γ)^3 / (9 * π * (1/(1 - F_r) * α_va) ^ 2)) ^ (κ), bi[i] + μ + p3.β_va + κ*(6 - 2 * p3.β_va), ci[i] + λ)
            else 
                # approximating sigma as 2 (for closed form integration) (this overestimates A LOT)
                v += integrate(th.D_cr, Inf, 1/(1 - F_r) * α_va * ai[i] * N_0 * (16 * p3.ρ_i ^ 2 * (F_r * π / 4 + (1 - F_r) * p3.γ)^3 / (9 * π * (1/(1 - F_r) * α_va) ^ 2)) ^ (κ), bi[i] + μ + p3.β_va + κ*(6 - 2 * p3.β_va), ci[i] + λ)
            end
        end
    end

    return v / q
end

"""
    terminal_velocity_number(p3, Chen2022, q, N, ρ_r, F_r)

 - p3 - a struct with P3 scheme parameters
 - Chen 2022 - a struch with terminal velocity parameters as in Chen(2022)
 - q - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i)
 - ρ_a - density of air 

 Returns the number (total)-weighted fall speed
 Eq C11 of Morrison and Milbrandt (2015)
"""
function terminal_velocity_number(p3::PSP3, Chen2022::CMP.Chen2022VelTypeSnowIce, q::FT, N::FT, ρ_r::FT, F_r::FT, ρ_a::FT) where{FT}
    # Get the thresholds for different particles regimes
    th = thresholds(p3, ρ_r, F_r)
    D_th = D_th_helper(p3)

    # Get the shape parameters
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, F_r)
    μ = DSD_μ(p3, λ)

    # Get the ai, bi, ci constants (in si units) for velocity calculations
    (; As, Bs, Cs, Es, Fs, Gs) = Chen2022

    bi = [Bs + ρ_a * Cs, Bs + ρ_a * Cs]
    ai = [Es * ρ_a^As * 10^(3 * bi[1]), Fs * ρ_a^As * 10^(3 * bi[2])]
    ci = [0, Gs * 10^3]

    κ = FT(1/3)

    # Redefine α_va to be in si units
    α_va = α_va_si(p3)

    # TODO Update the velocity to use a different formulation for D > 0.625 mm 
    v = 0
    for i in 1:2
        if F_r == 0 
            v += integrate(FT(0), D_th, ai[i] * N_0, bi[i] + μ, ci[i] + λ)
            v += integrate(D_th, Inf, ai[i] * N_0 * (16 * p3.ρ_i ^ 2 * p3.γ^3/(9 * π * α_va ^ 2)) ^ κ, bi[i] + μ + κ * (3 * p3.σ - 2 * p3.β_va), ci[i] + λ)
        else 
            v += integrate(FT(0), D_th, ai[i] * N_0, bi[i] + μ, ci[i] + λ)
            v += integrate(D_th, th.D_gr, ai[i] * N_0 * (16 * p3.ρ_i ^ 2 * p3.γ^3/(9 * π * α_va ^ 2)) ^ κ, bi[i] + μ + κ * (3 * p3.σ - 2 * p3.β_va), ci[i] + λ)
            v += integrate(th.D_gr, th.D_cr, ai[i] * N_0 * (p3.ρ_i / th.ρ_g)^(2 * κ), bi[i] + μ, ci[i] + λ)
            v += (
                integrate(th.D_cr, Inf, ai[i] * N_0 * (p3.γ * (1-F_r)) ^ (3 * κ), bi[i] + μ + κ*(3 * p3.σ), ci[i] + λ) + 
                integrate(th.D_cr, Inf, ai[i] * N_0 * (p3.γ ^ 2 * π * F_r * 3/4 * (1-F_r)^2) ^ κ, bi[i] + μ + κ*(2 + 2 * p3.σ), ci[i] + λ) + 
                integrate(th.D_cr, Inf, ai[i] * N_0 * (p3.γ * π ^ 2 * F_r ^ 2 * 3/16 * (1-F_r)) ^ κ, bi[i] + μ + κ*(4 + p3.σ), ci[i] + λ) + 
                integrate(th.D_cr, Inf, ai[i] * N_0 * (F_r^3 * π^3 / 64) ^ κ, bi[i] + μ + 6 * κ, ci[i] + λ)
            )
        end
    end

    return v / N
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
        n += integrate(FT(0), D_th, π / 6 * p3.ρ_i * N_0, μ + 4, λ)
        n += integrate(D_th, Inf, α_va * N_0, μ + p3.β_va + 1, λ)
    else 
        n += integrate(FT(0), D_th, π / 6 * p3.ρ_i * N_0, μ + 4, λ) 
        n += integrate(D_th, th.D_gr, α_va * N_0, μ + p3.β_va + 1, λ)
        n += integrate(th.D_gr, th.D_cr, π / 6 * th.ρ_g * N_0, μ + 4, λ) 
        n += integrate(th.D_cr, Inf, α_va / (1 - F_r) * N_0, μ + p3.β_va + 1, λ)
    end
    
    # Normalize by q
    return n/q

end

end
