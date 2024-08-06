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
    guess = log(FT(19) / 6 * (DSD_μ(p3, λ) - 1) + 39) - log(λ)
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
    L_(p3, ρ, F_rim, F_liq, λ, μ, D_min, D_max)

 - p3 - a struct with P3 scheme parameters
 - ρ - bulk ice density (ρ_i for small ice, ρ_g for graupel) [kg/m^3]
 - F_rim - rime mass fraction (L_rim/(L_i - L_liq)) [-]
 - F_liq - liquid fraction (L_liq/L_i,tot):
    - zero if solving for ice core shape parameters
 - μ - shape parameter of N′ gamma distribution
 - λ - slope parameter of N′ gamma distribution
 - D_min - minimum bound for regime
 - D_max - maximum bound for regime (if not specified, then infinity)

 Returns ice mass density for a given m(D) regime
"""
# small, spherical ice or graupel (completely rimed, spherical)
# D_min = 0, D_max = D_th, ρ = ρᵢ
# or
# L_rim > 0 and D_min = D_gr, D_max = D_cr, ρ = ρ_g
function L_s(
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
# L_rim = 0 and D_min = D_th, D_max = inf
function L_rz(p3::PSP3, F_liq::FT, μ::FT, λ::FT, D_min::FT) where {FT}
    return ∫_Γ(D_min, FT(Inf), (1 - F_liq) * α_va_si(p3), μ + p3.β_va, λ)
end
# L_rim > 0 and D_min = D_th and D_max = D_gr
function L_n(p3::PSP3, F_liq::FT, μ::FT, λ::FT, D_min::FT, D_max::FT) where {FT}
    return ∫_Γ(D_min, D_max, (1 - F_liq) * α_va_si(p3), μ + p3.β_va, λ)
end
# partially rimed ice or large unrimed ice (upper bound on D is infinity)
# L_rim > 0 and D_min = D_cr, D_max = inf
function L_r(p3::PSP3, F_rim::FT, F_liq::FT, μ::FT, λ::FT, D_min::FT) where {FT}
    return ∫_Γ(
        D_min,
        FT(Inf),
        (1 - F_liq) * α_va_si(p3) / (1 - F_rim),
        μ + p3.β_va,
        λ,
    )
end
# F_liq != 0 (liquid mass on mixed-phase particles for D in [D_min, D_max])
function L_liq(p3::PSP3, F_liq::FT, μ::FT, λ::FT) where {FT}
    return ∫_Γ(FT(0), FT(Inf), F_liq * (FT(π) / 6) * p3.ρ_l, μ + 3, λ)
end

"""
    L_over_N_gamma(p3, F_rim, F_liq, log_λ, th)

 - p3 - a struct with P3 scheme parameters
 - F_rim - rime mass fraction (L_rim/(L_i - L_liq)) [-]
 - F_liq - liquid fraction (L_liq/L_i,tot):
    - zero if solving for ice core shape parameters
 - log_λ - logarithm of the slope parameter of N′ gamma distribution
 - μ - shape parameter of N′ gamma distribution
 - th - thresholds() nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns L/N for all values of D (sum over all regimes).
Eq. 5 in Morrison and Milbrandt (2015).
"""
function L_over_N_gamma(
    p3::PSP3,
    F_rim::FT,
    F_liq::FT,
    log_λ::FT,
    μ::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT}

    D_th = D_th_helper(p3)
    λ = exp(log_λ)
    N = Γ(1 + μ) / (λ^(1 + μ))

    return ifelse(
        F_rim == FT(0),
        (
            L_s(p3, F_liq, p3.ρ_i, μ, λ, FT(0), D_th) +
            L_rz(p3, F_liq, μ, λ, D_th) +
            L_liq(p3, F_liq, μ, λ)
        ) / N,
        (
            L_s(p3, F_liq, p3.ρ_i, μ, λ, FT(0), D_th) +
            L_n(p3, F_liq, μ, λ, D_th, th.D_gr) +
            L_s(p3, F_liq, th.ρ_g, μ, λ, th.D_gr, th.D_cr) +
            L_r(p3, F_rim, F_liq, μ, λ, th.D_cr) +
            L_liq(p3, F_liq, μ, λ)
        ) / N,
    )
end

"""
    DSD_μ_approx(p3, L, N, ρ_r, F_rim, F_liq)

 - p3 - a struct with P3 scheme parameters
 - L - mass mixing ratio
 - N - total ice number concentration [1/m3]
 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_rim - rime mass fraction (L_rim/(L_i - L_liq)) [-]
 - F_liq - liquid fraction (L_liq/L_i,tot):
    - zero if solving for ice core shape parameters

Returns the approximated shape parameter μ for a given L and N value
"""
function DSD_μ_approx(
    p3::PSP3,
    L::FT,
    N::FT,
    ρ_r::FT,
    F_rim::FT,
    F_liq::FT,
) where {FT}
    # Get thresholds for given F_rim, ρ_r
    th = thresholds(p3, ρ_r, F_rim)

    # Get min and max lambda values
    λ_0 = μ_to_λ(p3, FT(0))
    λ_6 = μ_to_λ(p3, p3.μ_max)

    # Get corresponding L/N values at given F_rim
    L_over_N_min = log(L_over_N_gamma(p3, F_rim, F_liq, log(λ_0), FT(0), th))
    L_over_N_max = log(L_over_N_gamma(p3, F_rim, F_liq, log(λ_6), p3.μ_max, th))

    # Return approximation between them
    μ = (p3.μ_max / (L_over_N_max - L_over_N_min)) * (log(L / N) - L_over_N_min)

    # Clip approximation between 0 and 6
    return min(p3.μ_max, max(FT(0), μ))
end

"""
    get_bounds(N, L, F_rim, F_liq, p3, th)

 - N - ice number concentration [1/m3]
 - L - mass mixing ratio
 - μ - shape parameter of N′ gamma distribution
 - F_rim - rime mass fraction (L_rim/(L_i - L_liq)) [-]
 - F_liq - liquid fraction (L_liq/L_i,tot):
    - zero if solving for ice core shape parameters
 - p3 - a struct with P3 scheme parameters
 - th -  thresholds() nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

 Returns estimated guess for λ from L to be used in distribution_parameter_solver()
"""
function get_bounds(
    N::FT,
    L::FT,
    μ::FT,
    F_rim::FT,
    F_liq::FT,
    p3::PSP3,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT}
    goal = L / N

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

    Ll = L_over_N_gamma(p3, F_rim, F_liq, log(left), μ, th)
    Lr = L_over_N_gamma(p3, F_rim, F_liq, log(right), μ, th)

    guess =
        left * (goal / (Ll))^((log(right) - log(left)) / (log(Lr) - log(Ll)))

    max = log(guess * exp(radius))
    min = log(guess)

    return (; min, max)
end

"""
    distrbution_parameter_solver()

 - p3 - a struct with P3 scheme parameters
 - L - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_rim - rime mass fraction (L_rim/(L_i - L_liq)) [-]
 - F_liq - liquid fraction (L_liq/L_i,tot):
    - zero if solving for ice core shape parameters

Solves the nonlinear system consisting of N_0 and λ for P3 prognostic variables
Returns a named tuple containing:
 - N_0 - intercept size distribution parameter [1/m4]
 - λ - slope size distribution parameter [1/m]
"""
function distribution_parameter_solver(
    p3::PSP3{FT},
    L::FT,
    N::FT,
    ρ_r::FT,
    F_rim::FT,
    F_liq::FT,
) where {FT}
    # Get the thresholds for different particles regimes
    th = thresholds(p3, ρ_r, F_rim)

    # Get μ given L and N
    μ = DSD_μ_approx(p3, L, N, ρ_r, F_rim, F_liq)

    # To ensure that λ is positive solve for x such that λ = exp(x)
    shape_problem(x) = L / N - L_over_N_gamma(p3, F_rim, F_liq, x, μ, th)

    # Get intial guess for solver
    (; min, max) = get_bounds(N, L, μ, F_rim, F_liq, p3, th)

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
    D_m (p3, L, N, ρ_r, F_rim, F_liq)

 - p3 - a struct with P3 scheme parameters
 - L - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_rim - rime mass fraction (L_rim/(L_i - L_liq)) [-]
 - F_liq - liquid fraction (L_liq/L_i,tot):
    - zero if solving for ice core
    mass weighted mean particle size

 Return the mass weighted mean particle size [m]
"""
function D_m(p3::PSP3, L::FT, N::FT, ρ_r::FT, F_rim::FT, F_liq::FT) where {FT}
    # Get the thresholds for different particles regimes
    th = thresholds(p3, ρ_r, F_rim)
    D_th = D_th_helper(p3)

    # Get the shape parameters
    (λ, N_0) = distribution_parameter_solver(p3, L, N, ρ_r, F_rim, F_liq)
    μ = DSD_μ(p3, λ)

    # Redefine α_va to be in si units
    α_va = α_va_si(p3)

    # Calculate numerator
    n_nl = 0
    if F_rim == 0
        n_nl += ∫_Γ(FT(0), D_th, π / 6 * p3.ρ_i * N_0, μ + 4, λ)
        n_nl += ∫_Γ(D_th, FT(Inf), α_va * N_0, μ + p3.β_va + 1, λ)
        n_l = ∫_Γ(FT(0), FT(Inf), N_0 * p3.ρ_l * (FT(π) / 6), μ + 4, λ)
    else
        n_nl += ∫_Γ(FT(0), D_th, π / 6 * p3.ρ_i * N_0, μ + 4, λ)
        n_nl += ∫_Γ(D_th, th.D_gr, α_va * N_0, μ + p3.β_va + 1, λ)
        n_nl += ∫_Γ(th.D_gr, th.D_cr, π / 6 * th.ρ_g * N_0, μ + 4, λ)
        n_nl +=
            ∫_Γ(th.D_cr, FT(Inf), α_va / (1 - F_rim) * N_0, μ + p3.β_va + 1, λ)
        n_l = ∫_Γ(FT(0), FT(Inf), N_0 * p3.ρ_l * (FT(π) / 6), μ + 4, λ)
    end


    # F_liq-weighted average:
    n = (1 - F_liq) * n_nl + F_liq * n_l

    # Normalize by q
    return n / L
end
