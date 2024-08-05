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
    return ifelse(D < eps(FT), FT(0), N_0 * D^DSD_μ(p3, λ) * exp(-λ * D))
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
    if λ < FT(2)
        return (5 \ λ)
    else
        ice_problem(x) =
            tolerance - Γ(1 + DSD_μ(p3, λ), exp(x) * λ) / Γ(1 + DSD_μ(p3, λ))
        low_guess = log(eps(FT))
        high_guess = log(FT(4))
        guess = log(FT(19) / 6 * (DSD_μ(p3, λ) - 1) + 39) - log(λ)
        guess = Base.min(high_guess, Base.max(low_guess, guess))
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

"""
    L_(p3, ρ, F_rim, λ, μ, D_min, D_max)

 - p3 - a struct with P3 scheme parameters
 - ρ - bulk ice density (ρ_i for small ice, ρ_g for graupel) [kg/m^3]
 - F_rim - rime mass fraction [L_rim/L_ice]
 - μ - shape parameter of N′ gamma distribution
 - λ - slope parameter of N′ gamma distribution
 - D_min - minimum bound for regime
 - D_max - maximum bound for regime (if not specified, then infinity)

 Returns ice content for a given m(D) regime
"""
# small, spherical ice or graupel (completely rimed, spherical)
# D_min = 0, D_max = D_th, ρ = ρᵢ
# or
# L_rim > 0 and D_min = D_gr, D_max = D_cr, ρ = ρ_g
function L_s(p3::PSP3, ρ::FT, μ::FT, λ::FT, D_min::FT, D_max::FT) where {FT}
    return ∫_Γ(D_min, D_max, FT(π) / 6 * ρ, μ + 3, λ)
end
# L_rim = 0 and D_min = D_th, D_max = inf
function L_rz(p3::PSP3, μ::FT, λ::FT, D_min::FT) where {FT}
    return ∫_Γ(D_min, FT(Inf), α_va_si(p3), μ + p3.β_va, λ)
end
# L_rim > 0 and D_min = D_th and D_max = D_gr
function L_n(p3::PSP3, μ::FT, λ::FT, D_min::FT, D_max::FT) where {FT}
    return ∫_Γ(D_min, D_max, α_va_si(p3), μ + p3.β_va, λ)
end
# partially rimed ice or large unrimed ice (upper bound on D is infinity)
# L_rim > 0 and D_min = D_cr, D_max = inf
function L_r(p3::PSP3, F_rim::FT, μ::FT, λ::FT, D_min::FT) where {FT}
    return ∫_Γ(D_min, FT(Inf), α_va_si(p3) / (1 - F_rim), μ + p3.β_va, λ)
end

"""
    L_over_N_gamma(p3, F_rim, λ, th)

 - p3 - a struct with P3 scheme parameters
 - F_rim - rime mass fraction [L_rim/L_ice]
 - log_λ - logarithm of the slope parameter of N′ gamma distribution
 - μ - shape parameter of N′ gamma distribution
 - th - thresholds() nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns L/N for all values of D (sum over all regimes).
Eq. 5 in Morrison and Milbrandt (2015).
"""
function L_over_N_gamma(
    p3::PSP3,
    F_rim::FT,
    log_λ::FT,
    μ::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT}
    if exp(log_λ) < eps(FT)
        throw(AssertionError("(1 / λ) too large; un-physical PSD"))
    else
        D_th = D_th_helper(p3)
        λ = exp(log_λ)
        N = Γ(1 + μ) / (λ^(1 + μ))

        L = 0
        if F_rim == FT(0)
            L += L_s(p3, p3.ρ_i, μ, λ, FT(0), D_th)
            L += L_rz(p3, μ, λ, D_th)
        else
            L += L_s(p3, p3.ρ_i, μ, λ, FT(0), D_th)
            L += L_n(p3, μ, λ, D_th, th.D_gr)
            L += L_s(p3, th.ρ_g, μ, λ, th.D_gr, th.D_cr)
            L += L_r(p3, F_rim, μ, λ, th.D_cr)
        end
        # if N = 0, return some form of "mean mass" based on
        # mean D = 1 / λ
        # TODO - come up with a better way to treat this case in the
        # context of the general solver architecture
        return ifelse(N < eps(FT), p3_mass(p3, (1 / λ), F_rim, th), L / N)
    end
end

"""
    DSD_μ_approx(p3, L_ice, N, ρ_r, F_rim)

 - p3 - a struct with P3 scheme parameters
 - L_ice - ice content [kg/m3]
 - N - total ice number concentration [1/m3]
 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_rim - rime mass fraction (L_rim/L_ice)

Returns the approximated shape parameter μ for a given q and N value
"""
function DSD_μ_approx(p3::PSP3, L::FT, N::FT, ρ_r::FT, F_rim::FT) where {FT}
    if L < FT(0)
        throw("negative mass")
    elseif N < eps(FT)
        return FT(0)
    else
        # Get thresholds for given F_rim, ρ_r
        th = thresholds(p3, ρ_r, F_rim)

        # Get min and max lambda values
        λ_0 = μ_to_λ(p3, FT(0))
        λ_6 = μ_to_λ(p3, p3.μ_max)

        # Get corresponding L/N values at given F_rim
        L_over_N_min = L_over_N_gamma(p3, F_rim, log(λ_0), FT(0), th)
        L_over_N_max = L_over_N_gamma(p3, F_rim, log(λ_6), p3.μ_max, th)

        # Return approximation between them
        μ =
            (p3.μ_max / (log(L_over_N_max) - log(L_over_N_min))) *
            (log(L / N) - log(L_over_N_min))

        # Clip approximation between 0 and 6
        return min(p3.μ_max, max(FT(0), μ))
    end
end

"""
    get_bounds(N, L, F_rim, p3, th)

 - N - ice number concentration [1/m3]
 - L - ice content [kg/m3]
 - μ - shape parameter of N′ gamma distribution
 - F_rim -rime mass fraction [L_rim/L_ice]
 - p3 - a struct with P3 scheme parameters
 - th -  thresholds() nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

 Returns estimated guess for λ from L to be used in distribution_parameter_solver()
"""
function get_bounds(
    N::FT,
    L::FT,
    μ::FT,
    F_rim::FT,
    p3::PSP3,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT}
    if N < eps(FT) || L < eps(FT)
        return (; min = FT(0), max = eps(FT))
        # TODO - if N = 0, then set guess to something meaningful
        # (i.e. set min, max = 0, eps(FT) to search for λ in this range)
        # ALSO - if L = 0 -> return meaningful min, max
    else
        goal = ifelse(L < eps(FT), FT(0), L / N)
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

        Ll = L_over_N_gamma(p3, F_rim, log(left), μ, th)
        Lr = L_over_N_gamma(p3, F_rim, log(right), μ, th)

        low_guess = eps(FT)
        high_guess = FT(1e5)
        guess =
            left *
            (goal / (Ll))^((log(right) - log(left)) / (log(Lr) - log(Ll)))
        guess = Base.min(high_guess, Base.max(low_guess, guess))
        
        max = log(guess * exp(radius))
        min = log(guess)

        return (; min, max)
    end
end

"""
    distrbution_parameter_solver()

 - p3 - a struct with P3 scheme parameters
 - L - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_rim - rime mass fraction (L_rim/L_ice)

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
) where {FT}
    # Get the thresholds for different particles regimes
    th = thresholds(p3, ρ_r, F_rim)

    # Get μ given L and N
    μ = DSD_μ_approx(p3, L, N, ρ_r, F_rim)

    # To ensure that λ is positive solve for x such that λ = exp(x)
    shape_problem(x) = L / N - L_over_N_gamma(p3, F_rim, x, μ, th)

    # Get intial guess for solver
    (; min, max) = get_bounds(N, L, μ, F_rim, p3, th)
    # Find slope parameter
    x =
        RS.find_zero(
            shape_problem,
            RS.SecantMethod(min, max),
            RS.CompactSolution(),
            RS.RelativeSolutionTolerance(eps(FT)),
            5,
        ).root

    params = ifelse(
        N < eps(FT),
        (; λ = FT(1e4), N_0 = FT(0)),
        (; λ = exp(x), N_0 = DSD_N₀(p3, N, exp(x))),
    )

    return params
end

"""
    D_m (p3, L, N, ρ_r, F_rim)

 - p3 - a struct with P3 scheme parameters
 - L - mass mixing ratio
 - N - number mixing ratio
 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_rim - rime mass fraction (L_rim/L_ice)

 Return the mass weighted mean particle size [m]
"""
function D_m(p3::PSP3, L::FT, N::FT, ρ_r::FT, F_rim::FT) where {FT}
    # Get the thresholds for different particles regimes
    th = thresholds(p3, ρ_r, F_rim)
    D_th = D_th_helper(p3)

    # Get the shape parameters
    (λ, N_0) = distribution_parameter_solver(p3, L, N, ρ_r, F_rim)
    μ = DSD_μ(p3, λ)

    # Redefine α_va to be in si units
    α_va = α_va_si(p3)

    # Calculate numerator
    n = 0
    if F_rim == 0
        n += ∫_Γ(FT(0), D_th, π / 6 * p3.ρ_i * N_0, μ + 4, λ)
        n += ∫_Γ(D_th, Inf, α_va * N_0, μ + p3.β_va + 1, λ)
    else
        n += ∫_Γ(FT(0), D_th, π / 6 * p3.ρ_i * N_0, μ + 4, λ)
        n += ∫_Γ(D_th, th.D_gr, α_va * N_0, μ + p3.β_va + 1, λ)
        n += ∫_Γ(th.D_gr, th.D_cr, π / 6 * th.ρ_g * N_0, μ + 4, λ)
        n += ∫_Γ(th.D_cr, Inf, α_va / (1 - F_rim) * N_0, μ + p3.β_va + 1, λ)
    end
    # Normalize by L
    return ifelse(L < eps(FT), FT(0), n / L)
end
