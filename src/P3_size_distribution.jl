
"""
    log_N′ice(dist, D)

Compute the log of the ice particle number concentration at diameter `D` given the distribution `dist`
"""
function log_N′ice(dist::P3Distribution, D)
    (; log_N₀, log_λ) = dist
    μ = get_μ(dist)
    return log_N₀ + μ * log(D) - exp(log_λ) * D
end

"""
    N′ice(dist, D)

Compute the ice particle number concentration at diameter `D` given the distribution `dist`
"""
N′ice(dist::P3Distribution, D) = exp(log_N′ice(dist, D))



"""
    log_integrate_moment_psd(x, y, a, b, μ, log_λ)

Computes the log of the integral of the moment of the PSD from x to y

i.e. integral of the form
    ``∫_x^y (aD^b) D^μ e^{-λD} dD``
"""
function log_integrate_moment_psd(x, y, a, b, μ, log_λ)
    b_μ_1 = b + μ + 1
    (_, q_x) = SF.gamma_inc(b_μ_1, exp(log_λ) * x)
    (_, q_y) = SF.gamma_inc(b_μ_1, exp(log_λ) * y)
    return log(a) - b_μ_1 * log_λ + SF.loggamma(b_μ_1) + log(abs(q_x - q_y))
end

"""
    get_μ(dist)
    get_μ(state, log_λ)
    get_μ(params, log_λ)
    
Compute the slope parameter μ
"""
get_μ(dist::P3Distribution)  = get_μ(get_state(dist), dist.log_λ)
get_μ(state::P3State, log_λ) = get_μ(get_parameters(state), log_λ)
get_μ(params::PSP3, log_λ)   = get_μ(params.slope, log_λ)

get_μ((; a, b, c, μ_max)::SlopePowerLaw, log_λ) = clamp(a * exp(log_λ)^b - c, 0, μ_max)
get_μ((; μ)::SlopeConstant, _) = μ

"""
    log_LdN₀(state, log_λ)

Compute `log(L/N₀)` given the `state` and `log_λ`
"""
function log_LdN₀(state::P3State{FT}, log_λ) where {FT}
    (; F_rim, ρ_g, D_th, D_gr, D_cr) = state
    (; ρ_i, mass) = get_parameters(state)
    (; α_va, β_va) = mass
    
    μ = get_μ(state, log_λ)
    G = log_integrate_moment_psd
    G_small_spherical = G(FT(0), D_th, ρ_i * FT(π) / 6, 3, μ, log_λ)
    # G_liqfrac = G(FT(0), FT(Inf), ρ_l * π / 6, 3, μ, log_λ)  # TODO: Implement liquid fraction (need to do weighted average)
    if isunrimed(state)
        # L_rim = 0
        G_large_unrimed = G(D_th, FT(Inf), α_va, β_va, μ, log_λ)
        return logsumexp((G_small_spherical, G_large_unrimed))
    else
        # L_rim > 0
        G_dense_nonspherical =  G(D_th, D_gr,     α_va,                 β_va,   μ, log_λ)
        G_graupel =             G(D_gr, D_cr,     ρ_g * FT(π) / 6,      3,      μ, log_λ)
        G_partially_rimed =     G(D_cr, FT(Inf),  α_va / (1 - F_rim),   β_va,   μ, log_λ)
        return logsumexp((G_small_spherical, G_dense_nonspherical, G_graupel, G_partially_rimed))
    end
end

"""
    log_NdN₀(state, log_λ)
    log_NdN₀(slope, log_λ)

Compute `log(N/N₀)` given either the `state` or the `slope` parameterization, and `log_λ`

Note: This function is equivalent to `log_integrate_moment_psd(0, Inf, 1, 0, μ, log_λ)`
"""
log_NdN₀(state::P3State, log_λ) = log_NdN₀(get_parameters(state).slope, log_λ)
function log_NdN₀(slope::SlopeLaw, log_λ)
    μ = get_μ(slope, log_λ)
    return SF.loggamma(μ + 1) - (μ + 1) * log_λ
end

"""
    log_LdN(state, log_λ)

Compute `log(L/N)` given the `state` and `log_λ`
"""
log_LdN(state::P3State, log_λ) = log_LdN₀(state, log_λ) - log_NdN₀(state, log_λ)

"""
    get_log_N₀(state, log_N, log_λ)

Compute `log(N₀)` given the `state`, `log_N`, and `log_λ`
"""
function get_log_N₀(state::P3State, log_N, log_λ)
    μ = get_μ(state, log_λ)
    return log_N - SF.loggamma(μ + 1) + (μ + 1) * log_λ
end

function get_log_λ_initial_bounds(
    state::P3State{FT},
    target_log_mass_per_particle::FT,
) where {FT}

    # Choose λ range based on mass per particle
    if target_log_mass_per_particle >= log(FT(1e-8))
        λ_min = FT(1)        # Small particles
        λ_max = FT(6e3)
        search_width = FT(0.2)
    elseif target_log_mass_per_particle >= log(FT(2e-9))
        λ_min = FT(6e3)      # Medium particles
        λ_max = FT(3e4)
        search_width = FT(-0.1)
    else
        λ_min = FT(4e4)      # Large particles
        λ_max = FT(1e6)
        search_width = FT(0.2)
    end

    # Convert λ bounds to log space
    log_λ_min = log(λ_min)
    log_λ_max = log(λ_max)

    # Get mass per particle at bounds in log space
    mass_per_particle_at_min = log_LdN(state, log_λ_min)
    mass_per_particle_at_max = log_LdN(state, log_λ_max)

    # Linear interpolation in log space to estimate λ
    log_λ_guess = simple_linear_interpolation(
        mass_per_particle_at_min, log_λ_min,
        mass_per_particle_at_max, log_λ_max,
        target_log_mass_per_particle
    )

    # Define search interval around guess
    search_min = log_λ_guess - search_width
    search_max = log_λ_guess + search_width

    return (; min = search_min, max = search_max)
end

"""
    simple_linear_interpolation(x₁, y₁, x₂, y₂, x_target)

Simple linear interpolation between two points `(x₁, y₁)` and `(x₂, y₂)` to estimate the value at `x_target`
"""
function simple_linear_interpolation(x₁, y₁, x₂, y₂, x_target)
    return y₁ + (y₂ - y₁) * (x_target - x₁) / (x₂ - x₁)
end

function distribution_parameter_solver(
    state::P3State{FT},
    log_L::FT,
    log_N::FT,
    search_bounds = (; min = log(FT(10)), max = log(FT(1e7)))
) where {FT}
    target_log_LdN = log_L - log_N

    shape_problem(log_λ) = log_LdN(state, log_λ) - target_log_LdN


    # TODO: I think bound is wrong -- compare with `get_bounds`
    # search_bounds = get_log_λ_initial_bounds(state, target_log_LdN)

    # Find slope parameter
    sol = RS.find_zero(
        shape_problem,
        RS.SecantMethod(search_bounds.min, search_bounds.max),
        RS.VerboseSolution(),
        # RS.CompactSolution(),
        RS.RelativeSolutionTolerance(eps(FT)),
        50,
    )
    log_λ = sol.root

    log_N₀ = get_log_N₀(state, log_N, log_λ)

    return P3Distribution(state, log_λ, log_N₀)
end

"""
    distribution_parameter_solver_all_solutions(p3, L, N, ρ_r, F_rim, F_liq)

Find all solutions for λ
"""
function distribution_parameter_solver_all_solutions(state::P3State{FT}, log_L, log_N) where {FT}
    # Find bounds by evaluating function incrementally, then apply root finding with bounds above and below zero-point
    target_log_LdN = log_L - log_N
    
    shape_problem(log_λ) = log_LdN(state, log_λ) - target_log_LdN

    Δλ = 0.01
    λs = 10.0 .^ (2.0:Δλ:6.0)
    λ_bnds = Tuple[]
    # Loop over λs and find where shape_problem changes sign
    for i in 1:length(λs)-1
        if shape_problem(log(λs[i])) * shape_problem(log(λs[i+1])) < 0
            push!(λ_bnds, (λs[i], λs[i+1]))
        end
    end
    
    # Apply root finding with bounds above and below zero-point
    dists = P3Distribution[]
    for (λ_min, λ_max) in λ_bnds
        log_λ_min, log_λ_max = log(λ_min), log(λ_max)
        dist = distribution_parameter_solver(state, log_L, log_N, (; min = log_λ_min, max = log_λ_max))
        push!(dists, dist)
    end
    return dists
end

"""
    D_m (p3, L, N, ρ_r, F_rim)

 - p3 - a struct with P3 scheme parameters
 - L - ice mass content [kg/m3]
 - N - number concentration [1/m3]
 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_rim - rime mass fraction (L_rim / L_ice) [-]
 - F_liq - liquid fraction (L_liq / L_p3_tot) [-]:
    - zero if solving for ice core D_m
 - p3 - a struct with P3 scheme parameters

 Return the mass weighted mean particle size [m]
"""
# function D_m(p3::PSP3, L::FT, N::FT, ρ_r::FT, F_rim::FT, F_liq::FT) where {FT}
function D_m(state::P3State, log_L::FT, log_N::FT) where {FT}
    exp(L) < eps(FT) && return FT(0)

    (; F_rim, ρ_g, D_th, D_gr, D_cr) = state
    (; ρ_i, ρ_l, mass) = get_parameters(state)
    (; α_va, β_va) = mass
    
    # Get the shape parameters
    (; log_λ, log_N₀) = dist = distribution_parameter_solver(state, log_L, log_N)
    
    μ = get_μ(state, log_λ)

    # Calculate numerator
    G = log_integrate_moment_psd
    n_nl = 0
    ∞ = FT(Inf)
    G_small_spherical = G(0, D_th, ρ_i * π / 6, 4, μ, log_λ)
    G_liqfrac = G(0, ∞, ρ_l * π / 6, 4, μ, log_λ)
    G_summed = if isunrimed(state)
        G_large_unrimed = G(D_th, ∞, α_va, β_va + 1, μ, log_λ)
        logsumexp((G_small_spherical, G_large_unrimed))
    else
        G_dense_nonspherical = G(D_th, D_gr, α_va, β_va + 1, μ, log_λ)
        G_graupel = G(D_gr, D_cr, ρ_g * π / 6, 4, μ, log_λ)
        G_partially_rimed = G(D_cr, ∞, α_va / (1 - F_rim), β_va + 1, μ, log_λ)
        logsumexp((G_small_spherical, G_dense_nonspherical, G_graupel, G_partially_rimed))
    end

    # compute F_liq-weighted average and normalize by L
    return weighted_average(F_liq, exp(G_summed), exp(G_liqfrac)) * exp(log_N₀) / exp(log_L)
end
