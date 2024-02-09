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
import CLIMAParameters as CP
import CloudMicrophysics.Parameters as CMP

const PSP3 = CMP.ParametersP3

export thresholds, p3_mass, distribution_parameter_solver

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
    if D_th_helper(p3) > D
        return mass_s(D, p3.ρ_i)          # small spherical ice
    end
    if F_r == 0
        return mass_nl(p3, D)             # large nonspherical unrimed ice
    end
    if th.D_gr > D >= D_th_helper(p3)
        return mass_nl(p3, D)             # dense nonspherical ice
    end
    if th.D_cr > D >= th.D_gr
        return mass_s(D, th.ρ_g)          # graupel
    end
    if D >= th.D_cr
        return mass_r(p3, D, F_r)         # partially rimed ice
    end

    # TODO - would something like this be better?
    #return ifelse(D_th_helper(p3) > D, mass_s(D, p3.ρ_i),
    #       ifelse(F_r == 0, mass_nl(p3, D),
    #       ifelse(th.D_gr > D >= D_th_helper(p3), mass_nl(p3, D),
    #       ifelse(th.D_cr > D >= th.D_gr, mass_s(D, th.ρ_g),
    #           mass_r(p3, D, F_r)))))

end

# Some wrappers to cast types from SF.gamma
# (which returns Float64 even when the input is Float32)
Γ(a::FT, z::FT) where {FT <: Real} = FT(SF.gamma(a, z))
Γ(a::FT) where {FT <: Real} = FT(SF.gamma(a))

"""
    DSD_μ(p3, λ)

 - p3 - a struct with P3 scheme parameters
 - λ - slope parameter for gamma distribution of N′ [1/m]

Returns the shape parameter μ for a given λ value
Eq. 3 in Morrison and Milbrandt (2015).
"""
function DSD_μ(p3::PSP3, λ::FT) where {FT}
    #println("λ = ", λ)
    #@assert λ > FT(0)
    return min(p3.μ_max, max(FT(0), p3.a * λ^p3.b - p3.c))
end

"""
    DSD_N₀(p3, N, λ)
 - p3 - a struct with P3 scheme parameters
 - N - total ice number concentration [1/m3]
 - λ - slope parameter [1/m]

Returns the shape parameter N₀ from Eq. 2 in Morrison and Milbrandt (2015).
"""
function DSD_N₀(p3::PSP3, N::FT, λ::FT) where {FT}
    μ = DSD_μ(p3, λ)
    return N / Γ(1 + μ) * λ^(1 + μ)
end

"""
    q_(p3, ρ, F_r, N_0, λ, μ, D_min, D_max)

 - p3 - a struct with P3 scheme parameters
 - ρ - bulk ice density (ρ_i for small ice, ρ_g for graupel) [kg/m^3]
 - F_r - rime mass fraction [q_rim/q_i]
 - N_0 - intercept parameter of N′ gamma distribution
 - λ - slope parameter of N′ gamma distribution
 - D_min - minimum bound for regime
 - D_max - maximum bound for regime (if not specified, then infinity)

 Returns ice mass density for a given m(D) regime
"""
# small, spherical ice or graupel (completely rimed, spherical)
# D_min = 0, D_max = D_th, ρ = ρᵢ
# or
# q_rim > 0 and D_min = D_gr, D_max = D_cr, ρ = ρ_g
function q_s(p3::PSP3, ρ::FT, N_0::FT, λ::FT, D_min::FT, D_max::FT) where {FT}
    x = DSD_μ(p3, λ) + 4
    return FT(π) / 6 * ρ * N_0 / λ^x * (Γ(x, λ * D_min) - Γ(x, λ * D_max))
end
# q_rim = 0 and D_min = D_th, D_max = inf
function q_rz(p3::PSP3, N_0::FT, λ::FT, D_min::FT) where {FT}
    x = DSD_μ(p3, λ) + p3.β_va + 1
    return α_va_si(p3) * N_0 / λ^x * (Γ(x, λ * D_min))
end
# q_rim > 0 and D_min = D_th and D_max = D_gr
function q_n(p3::PSP3, N_0::FT, λ::FT, D_min::FT, D_max::FT) where {FT}
    x = DSD_μ(p3, λ) + p3.β_va + 1
    return α_va_si(p3) * N_0 / λ^x * (Γ(x, λ * D_min) - Γ(x, λ * D_max))
end
# partially rimed ice or large unrimed ice (upper bound on D is infinity)
# q_rim > 0 and D_min = D_cr, D_max = inf
function q_r(p3::PSP3, F_r::FT, N_0::FT, λ::FT, D_min::FT) where {FT}
    x = DSD_μ(p3, λ) + p3.β_va + 1
    return α_va_si(p3) * N_0 / (1 - F_r) / λ^x * (Γ(x, λ * D_min))
end

"""
    q_gamma(p3, F_r, N, λ, th)

 - p3 - a struct with P3 scheme parameters
 - F_r - rime mass fraction [q_rim/q_i]
 - N - ice number concentration [1/m3]
 - log_λ - logarithm of the slope parameter of N′ gamma distribution
 - th - thresholds() nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns ice mass density for all values of D (sum over all regimes).
Eq. 5 in Morrison and Milbrandt (2015).
"""
function q_gamma(
    p3::PSP3,
    F_r::FT,
    N::FT,
    log_λ::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT}

    D_th = D_th_helper(p3)
    λ = exp(log_λ)
    N_0 = DSD_N₀(p3, N, λ)
    
    return ifelse(
        F_r == FT(0),
        q_s(p3, p3.ρ_i, N_0, λ, FT(0), D_th) + q_rz(p3, N_0, λ, D_th),
        q_s(p3, p3.ρ_i, N_0, λ, FT(0), D_th) +
        q_n(p3, N_0, λ, D_th, th.D_gr) +
        q_s(p3, th.ρ_g, N_0, λ, th.D_gr, th.D_cr) +
        q_r(p3, F_r, N_0, λ, th.D_cr),
    )
end

"""
    get_bounds()

"""
function get_bounds(N::FT, q::FT, F_r::FT, p3::PSP3, th) where{FT}
    leftpt = FT(1)
    rightpt = FT(1e6)
    radius = FT(0.8)

    q_left = q_gamma(p3, F_r, N, log(leftpt), th)
    q_right = q_gamma(p3, F_r, N, log(rightpt), th)

    guess = (q/q_left)^((log(rightpt) - log(leftpt))/(log(q_right) - log(q_left)))

    max = log(guess * exp(radius))
    min = log(guess)

    if guess < FT(2.5 * 1e4)
        min = log(FT(20000))
        max = log(FT(50000))
    end
    return (min = min, max = max)
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

    # To ensure that λ is positive solve for x such that λ = exp(x)
    # We divide by N̂ to deal with large N₀ values for Float32
    N̂ = FT(1e20)
    shape_problem(x) = q / N̂ - q_gamma(p3, F_r, N / N̂, x, th)

    # Get intial guess for solver 
    (min, max) = get_bounds(N, q, F_r, p3, th)
    #min = log(20000)
    #max = log(50000)
    # Find slope parameter
    sol =
        RS.find_zero(
            shape_problem,
            RS.SecantMethod(FT(min), FT(max)),
            RS.CompactSolution(),
            RS.RelativeSolutionTolerance(eps(FT)),
            10,
        )
    x=sol.root
    converge = sol.converged
    #= x =
        RS.find_zero(
            shape_problem,
            RS.SecantMethod(FT(log(4 * 1e4)), FT(log(5 * 1e4))),
            RS.CompactSolution(),
            RS.RelativeSolutionTolerance(eps(FT)),
            10,
        ).root
 =#
    return (; λ = exp(x), N_0 = DSD_N₀(p3, N, exp(x)), converged = converge)
end

end
