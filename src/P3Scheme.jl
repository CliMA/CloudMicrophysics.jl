"""
Predicted particle properties scheme (P3) for ice, which includes:
 - threshold solver
 - m(D) regime
 - TODO

Implementation of Morrison and Milbrandt 2015 doi: 10.1175/JAS-D-14-0065.1

Note: Particle size is defined as its maximum length (i.e. max dimesion).
"""
module P3Scheme

import Integrals as IN
import RootSolvers as RS
import CLIMAParameters as CP
import CloudMicrophysics.Parameters as CMP
import SpecialFunctions as SF


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
α_va_si(p3::PSP3{FT}) where {FT} = FT(p3.α_va * 10^(6 * p3.β_va - 3))

"""
    D_th_helper(p3)

 - p3 - a struct with P3 scheme parameters

Returns the critical size separating spherical and nonspherical ice, in meters.
Eq. 8 in Morrison and Milbrandt (2015).
"""
D_th_helper(p3::PSP3) = (π * p3.ρ_i / 6 / α_va_si(p3))^(1 / (p3.β_va - 3))

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
    return (1 / (1 - F_r) * 6 * α_va / π / ρ_g)^(1 / (3 - p3.β_va))
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
    return (6 * α_va / π / ρ_g)^(1 / (3 - p3.β_va))
end

"""
    ρ_g_helper(ρ_r, F_r, ρ_d)

 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i) [-]
 - ρ_g - is the effective density of a spherical graupel particle [kg/m^3]

Returns the density of total (deposition + rime) ice mass for graupel, in kg/m3
Eq. 16 in Morrison and Milbrandt (2015).
"""
ρ_g_helper(ρ_r::FT, F_r::FT, ρ_d::FT) where {FT} =
    FT(F_r * ρ_r + (1 - F_r) * ρ_d)

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
    β_m2 = p3.β_va - FT(2)
    return FT(
        6 * α_va * (D_cr^β_m2 - D_gr^β_m2) / π / β_m2 /
        max(D_cr - D_gr, eps(FT)),
    )
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

    @assert ρ_r > FT(0)   # rime density must be positive ...
    @assert ρ_r <= p3.ρ_l # ... and as a bulk ice density can't exceed the density of water
    @assert F_r > FT(0)   # rime mass fraction must be positive ...
    @assert F_r < FT(1)   # ... and there must always be some unrimed part

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
mass_r(p3::PSP3, D::FT, F_r::FT) where {FT <: Real} = α_va_si(p3) / (1 - F_r) * D^p3.β_va

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
end

"""
    μ_calc(λ)

 - λ - slope parameter for gamma distribution of N′

 Returns the shape parameter (μ) corresponding to the given λ value
 Eq. 3 in Morrison and Milbrandt (2015).
"""
function μ_calc(λ::FT) where {FT}
    μ = 0.00191 * λ^(0.8) - 2
    if μ < 0
        μ = FT(0)
    elseif μ > 6
        μ = FT(6)
    end
    return μ
end

"""
    q_(p3, ρ, F_r, N_0, λ, μ, D_min, D_max)

 - p3 - a struct with P3 scheme parameters
 - ρ - bulk ice density (ρ_i for small ice, ρ_g for graupel) [kg/m^3]
 - F_r - rime mass fraction [q_rim/q_i]
 - N_0 - intercept parameter of N′ gamma distribution
 - λ - slope parameter of N′ gamma distribution 
 - μ - shape parameter of N′ gamma distribution
 - D_min - minimum bound for regime 
 - D_max - maximum bound for regime (if not specified, then infinity)

 Returns the prognostic mass mixing ratio for a given regime
"""
# small, spherical ice or graupel (completely rimed, spherical)
function q_s(ρ::FT, N_0::FT, λ::FT, μ::FT, D_min::FT, D_max::FT) where {FT} 
    return (π/6 * ρ * N_0) * λ^(-1 * (μ + 4)) * (SF.gamma(μ + 4, λ*D_min) - SF.gamma(μ + 4, λ*D_max)) 
end
# dense, non-spherical ice
function q_n(p3::PSP3, N_0::FT, λ::FT, μ::FT, D_min::FT, D_max::FT) where {FT} 
    return (p3.α_va * N_0) * (λ)^(-1 * (μ + p3.β_va + 1)) * (SF.gamma(μ + p3.β_va + 1, λ*D_min) - SF.gamma(μ + p3.β_va + 1, λ*D_max))
end
# partially rimed ice or large unrimed ice
function q_r(p3::PSP3, F_r::FT, N_0::FT, λ::FT, μ::FT, D_min::FT) where{FT} 
    return (p3.α_va * N_0 /(1 - F_r)) * (λ)^(-1 * (μ + p3.β_va + 1)) * (SF.gamma(μ + p3.β_va + 1) + (SF.gamma(μ + p3.β_va + 1, λ*D_min) - (μ + p3.β_va)*SF.gamma(μ + p3.β_va)))
end

"""
    q_gamma(p3, F_r, N_0, λ, μ, th)

 - p3 - a struct with P3 scheme parameters 
 - F_r - rime mass fraction [q_rim/q_i]
 - N_0 - intercept parameter of N′ gamma distribution 
 - λ - slope parameter of N′ gamma distribution 
 - μ - shape parameter of N′ gamma distribution 
 - th - thresholds() nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

 Returns prognostic mass mixing ratio for all values of D 
    Eq. 5 in Morrison and Milbrandt (2015).
"""
function q_gamma(p3::PSP3, F_r::FT, N_0::FT, λ::FT, μ::FT, th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0))) where{FT}
    D_th = D_th_helper(p3) 

    println("q_s = ", q_s(p3.ρ_i, N_0, λ, μ, FT(0), D_th))
    println("q_n = ", q_n(p3, N_0, λ, μ, D_th, th.D_gr))
    println("q_s 2 = ", q_s(th.ρ_g, N_0, λ, μ, th.D_gr, th.D_cr))
    println("q_r = ", q_r(p3, F_r, N_0, λ, μ, th.D_cr))

    if F_r == 0
        return q_s(p3.ρ_i, N_0, λ, μ, FT(0), D_th) + q_r(p3, F_r, N_0, λ, μ, D_th)
    else 
        return q_s(p3.ρ_i, N_0, λ, μ, FT(0), D_th) + q_n(p3, N_0, λ, μ, D_th, th.D_gr) + q_s(th.ρ_g, N_0, λ, μ, th.D_gr, th.D_cr) + q_r(p3, F_r, N_0, λ, μ, th.D_cr)
    end
end

"""
    N′(D, p)
    
 - D - maximum particle dimension
 - p - a tuple containing N_0, λ, μ (intrcept, slope, and shape parameters for N′ respectively)
 
 Returns the value of N′ 
 Eq. 2 in Morrison and Milbrandt (2015).   
"""
N′(D, p) = p.N_0 * D ^ (p.μ) * exp(-p.λ * D)

"""
"""
function N_0_helper(N::FT, λ::FT, μ::FT) where{FT}
    return N/(λ^(-1 - μ) * SF.gamma(1 + μ))
end

"""
    distrbution_parameter_solver()

 - p3 - a struct with P3 scheme parameters
 - q - mass mixing ratio 
 - N - number mixing ratio 
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i)

Solves the nonlinear system consisting of N_0 and λ
for a given number mixing ratio (N) and mass mixing ratio (q).
    Returns a named tuple containing:
     - N_0 - size distribution parameter related to N and q
     - λ - size distribution parameter related to N and q [m^-1]
"""
function distribution_parameter_solver(p3::PSP3{FT}, q::FT, N::FT, ρ_r::FT, F_r::FT) where {FT}

    th = thresholds(p3, ρ_r, F_r)
    
    shape_problem(λ) = q - q_gamma(p3, F_r, N_0_helper(N, λ, μ_calc(λ)), λ, μ_calc(λ), th)  

    λ = RS.find_zero(shape_problem, RS.SecantMethod(FT(100), FT(5000)), RS.CompactSolution(), RS.RelativeSolutionTolerance(1e-3), 10,).root 

    return(;
        λ = λ,
        N_0 = N_0_helper(N, λ, μ_calc(λ))
    )
end


"""
    N_helper(N_0, λ)

 - N_0 - intercept parameter of N′ gamma distribution
 - λ - slope parameter of N′ gamma distribution 

Returns the prognostic number mixing ratio 
Eq. 4 in Morrison and Milbrandt (2015).
"""
function N_helper(N_0::FT, λ::FT) where {FT}
    problem = IN.IntegralProblem(N′, 0, Inf, (N_0 = N_0, λ = λ, μ = μ_calc(λ)))
    sol = IN.solve(problem, IN.HCubatureJL(), reltol = 1e-10, abstol = 1e-10)
    return FT(sol.u)
end

"""
    q_helper(N_0, λ)

 - p3 - a struct with P3 scheme parameters
 - N_0 - intercept parameter of N′ gamma distribution
 - λ - slope parameter of N′ gamma distribution 
 - F_r - rime mass fraction (q_rim/q_i)
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]

Returns the prognostic mass mixing ratio
Eq. 5 in Morrison and Milbrandt (2015).
"""
function q_helper(p3::PSP3{FT}, N_0::FT, λ::FT, F_r::FT, ρ_r::FT) where {FT}
    th = thresholds(p3, ρ_r, F_r)
    q′(D, p) = p3_mass(p3, D, F_r, th) * N′(D, p)
    problem = IN.IntegralProblem(q′, 0, Inf, (N_0 = N_0, λ = λ, μ = μ_calc(λ)))
    sol = IN.solve(problem, IN.HCubatureJL(), reltol = 1e-10, abstol = 1e-10)
    return FT(sol.u)
end

end
