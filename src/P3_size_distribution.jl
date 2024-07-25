"""
    μ_to_λ(μ)

 - μ - parameter for gamma distribution of N′

Returns corresponding λ to given μ value
"""
function μ_to_λ(p3::PSP3, μ::FT) where {FT}
    return ((μ + p3.c) / p3.a)^(1 / p3.b)
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
    get_ice_bound(p3, λ, tolerance)

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
function get_ice_bound(p3::PSP3, λ::FT, tolerance::FT) where {FT}
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
