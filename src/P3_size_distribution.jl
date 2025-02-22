
### --------------------------------- ###
### ----- P3Distribution STRUCT ----- ###
### --------------------------------- ###

"""
    P3Distribution{FT}

Distribution of ice particles in size.

# Fields
$(FIELDS)
"""
@kwdef struct P3Distribution{FT}
    "Particle state, see [`P3State`](@ref)"
    state::P3State{FT}

    # Integral solution parameters
    "The mass concentration [kg/m³]"
    L::FT
    "The number concentration [1/m³]"
    N::FT

    # Distribution parameters
    "Logarithm of the slope parameter"
    log_λ::FT
    "Logarithm of the intercept parameter"
    log_N₀::FT
end

"""
    get_state(dist::P3Distribution)

Return the particle state from a [`P3Distribution`](@ref) object.
"""
get_state(dist::P3Distribution) = dist.state

Base.eltype(::P3Distribution{FT}) where {FT} = FT
Base.broadcastable(state::P3Distribution) = tuple(state)

"""
    get_parameters(dist::P3Distribution)

Return the parameters from a [`P3Distribution`](@ref) object.
"""
get_parameters(dist::P3Distribution) = get_parameters(get_state(dist))

"""
    isunrimed(dist::P3Distribution)

Return `true` if the particle state associated with a [`P3Distribution`](@ref) object is unrimed, `false` otherwise.
"""
isunrimed(dist::P3Distribution) = isunrimed(get_state(dist))


"""
    log_N′ice(dist, D)

Compute the log of the ice particle number concentration at diameter `D` 
    given the distribution `dist`
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


### ------------------------------------------------ ###
### ----- Obtaining P3 distribution parameters ----- ###
### ------------------------------------------------ ###


"""
    log_integrate_moment_psd(D₁, D₂, a, b, μ, log_λ)

Computes the log of the integral of the moment of the PSD from `D₁` to `D₂`

i.e. integral of the form
    ``∫_{D₁}^{D₂} (aD^b) D^μ e^{-λD} dD``
"""
function log_integrate_moment_psd(D₁, D₂, a, b, μ, log_λ)
    b_μ_1 = b + μ + 1
    (_, q_D₁) = SF.gamma_inc(b_μ_1, exp(log_λ) * D₁)
    (_, q_D₂) = SF.gamma_inc(b_μ_1, exp(log_λ) * D₂)
    return log(a) - b_μ_1 * log_λ + SF.loggamma(b_μ_1) + log(abs(q_D₁ - q_D₂))
end

"""
    get_μ(dist)
    get_μ(state, log_λ)
    get_μ(params, log_λ)
    
Compute the slope parameter μ

# Arguments
- `dist`: [`P3Distribution`](@ref) object
- `state`: [`P3State`](@ref) object
- `params`: [`CMP.ParametersP3`](@ref) object
- `log_λ`: The log of the slope parameter [log(1/m)]
"""
get_μ(dist::P3Distribution) = get_μ(get_state(dist), dist.log_λ)
get_μ(state::P3State, log_λ) = get_μ(get_parameters(state), log_λ)
get_μ(params::PSP3, log_λ) = get_μ(params.slope, log_λ)

get_μ((; a, b, c, μ_max)::CMP.SlopePowerLaw, log_λ) =
    clamp(a * exp(log_λ)^b - c, 0, μ_max)
get_μ((; μ)::CMP.SlopeConstant, _) = μ

"""
    log∫DⁿmN′dD(state, log_λ; n = 0)

Compute `log(∫_0^∞ Dⁿ m(D) N′(D) dD)` given the `state` and `log_λ`.
    This is the log of the `n`-th moment of the mass-weighted PSD.

# Arguments
- `state::P3State`: [`P3State`](@ref) object
- `log_λ`: The log of the slope parameter [log(1/m)]
- `n`: The order of the moment [dimensionless]

# Note:
- For `n = 0`, this evaluates to `log(L/N₀)`, see [`log_L_div_N₀`](@ref)
- For `n = 1`, this evaluates to the (unnormalized) mass-weighted mean particle size, see [`D_m`](@ref)
"""
function log∫DⁿmN′dD(state::P3State{FT}, log_λ; n = 0) where {FT}
    @assert n ≥ 0 "The moment order must be non-negative"
    (; F_rim, ρ_g, D_th, D_gr, D_cr) = state
    (; ρ_i, mass) = get_parameters(state)
    (; α_va, β_va) = mass

    μ = get_μ(state, log_λ)
    ∞ = FT(Inf)
    G = log_integrate_moment_psd
    # G_liqfrac = G(FT(0), FT(Inf), ρ_l * FT(π) / 6, 3 + n, μ, log_λ)  # TODO: Implement liquid fraction (need to do weighted average)
    G_small_spherical = G(0, D_th, ρ_i * π / 6, 3 + n, μ, log_λ)
    if isunrimed(state)
        G_large_unrimed = G(D_th, ∞, α_va, β_va + n, μ, log_λ)
        return LogExpFunctions.logsumexp((G_small_spherical, G_large_unrimed))
    else
        G_dense_nonspherical = G(D_th, D_gr, α_va, β_va + n, μ, log_λ)
        G_graupel = G(D_gr, D_cr, ρ_g * π / 6, 3 + n, μ, log_λ)
        G_partially_rimed = G(D_cr, ∞, α_va / (1 - F_rim), β_va + n, μ, log_λ)
        return LogExpFunctions.logsumexp((
            G_small_spherical,
            G_dense_nonspherical,
            G_graupel,
            G_partially_rimed,
        ))
    end
end


"""
    log_L_div_N₀(state, log_λ)

Compute `log(L/N₀)` given the `state` and `log_λ`
"""
log_L_div_N₀(state::P3State, log_λ) = log∫DⁿmN′dD(state, log_λ; n = 0)

"""
    log_N_div_N₀(state, log_λ)
    log_N_div_N₀(slope, log_λ)

Compute `log(N/N₀)` given either the `state` or the `slope` parameterization, and `log_λ`

Note: This function is equivalent to `log_integrate_moment_psd(0, Inf, 1, 0, μ, log_λ)`
"""
log_N_div_N₀(state::P3State, log_λ) = log_N_div_N₀(get_parameters(state), log_λ)
log_N_div_N₀(params::PSP3, log_λ) = log_N_div_N₀(params.slope, log_λ)
function log_N_div_N₀(slope::CMP.SlopeLaw, log_λ)
    μ = get_μ(slope, log_λ)
    return SF.loggamma(μ + 1) - (μ + 1) * log_λ
end

"""
    log_L_div_N(state, log_λ)

Compute `log(L/N)` given the `state` (or `params`) and `log_λ`

# Arguments
- `state::P3State`: [`P3State`](@ref) object, or
- `params::PSP3`: [`CMP.ParametersP3`](@ref) object, and
- `log_λ`: The log of the slope parameter [log(1/m)]
"""
log_L_div_N(state::P3State, log_λ) = log_L_div_N₀(state, log_λ) - log_N_div_N₀(state, log_λ)

"""
    get_log_N₀(state; N, log_λ)
    get_log_N₀(params; N, log_λ)

Compute `log(N₀)` given the `state`, `N`, and `log_λ`

# Arguments
- `state::P3State`: [`P3State`](@ref) object, or
- `params::PSP3`: [`CMP.ParametersP3`](@ref) object

# Keyword arguments
- `N`: The number concentration [1/m³]
- `log_λ`: The log of the slope parameter [log(1/m)]
"""
function get_log_N₀(state::P3State; N, log_λ)
    get_log_N₀(get_parameters(state); N, log_λ)
end
function get_log_N₀(params::PSP3; N, log_λ)
    μ = get_μ(params, log_λ)
    return log(N) - SF.loggamma(μ + 1) + (μ + 1) * log_λ
end

"""
    get_distribution_parameters(state; L, N, [log_λ_min, log_λ_max])

Solve for the distribution parameters given the state, and the mass (`L`) and number (`N`) concentrations.

The assumed distribution is of the form

```math
N(D) = N₀ D^μ e^{-λD}
```
where `N(D)` is the number concentration at diameter `D` and `μ` is the slope parameter.
    The slope parameter is parameterized, e.g. [`CMP.SlopePowerLaw`](@ref) or [`CMP.SlopeConstant`](@ref).

This algorithm solves for `log_λ = log(λ)` and `log_N₀ = log(N₀)` 
    given `L` and `N` by solving the equations:

```math
\\begin{align*}
\\log(L) &= \\log ∫_0^∞ m(D) N(D)\\ \\mathrm{d}D, \\\\
\\log(N) &= \\log ∫_0^∞ N(D)\\ \\mathrm{d}D, \\\\
\\end{align*}
```
where `m(D)` is the mass of a particle at diameter `D` (see [`ice_mass`](@ref)).
    The procedure is decribed in detail in [the P3 docs](@ref "Parameterizations for the slope parameter \$μ\$").


# Arguments
- `state`: [`P3State`](@ref) object

# Keyword arguments
- `L`: The mass concentration [kg/m³]
- `N`: The number concentration [1/m³]
- `search_bound_min`: The minimum value of the search bounds [log(1/m)], default is `log(1e1)`
- `search_bound_max`: The maximum value of the search bounds [log(1/m)], default is `log(1e7)`

# Returns
- [`P3Distribution`](@ref) object with the distribution parameters `log_λ` and `log_N₀`

# Examples

```jldoctest
julia> import CloudMicrophysics.Parameters as CMP,
              CloudMicrophysics.P3Scheme   as P3

julia> params = CMP.ParametersP3(Float64);

julia> state = P3.get_state(params; F_rim = 0.0, ρ_r = 400.0);

julia> dist = P3.get_distribution_parameters(state; L = 1e-3, N = 1e3)
P3Distribution{Float64}
├── state: is unrimed
├── log_λ = 5.4897008376530385 [log(1/m)]
└── log_N₀ = 12.397456116635176 [log(1/m^3)]
```
"""
function get_distribution_parameters(
    state::P3State{FT};
    L,
    N,
    log_λ_min = log(1e1),
    log_λ_max = log(1e7),
) where {FT}
    target_log_LdN = log(L) - log(N)

    shape_problem(log_λ) = log_L_div_N(state, log_λ) - target_log_LdN

    # Find slope parameter
    sol = RS.find_zero(
        shape_problem,
        RS.SecantMethod(FT(log_λ_min), FT(log_λ_max)),
        RS.CompactSolution(),
        RS.RelativeSolutionTolerance(eps(FT)),
        50,
    )
    log_λ = sol.root

    log_N₀ = get_log_N₀(state; N, log_λ)

    return P3Distribution(; state, L, N, log_λ, log_N₀)
end

"""
    get_distribution_parameters_all_solutions(state; L, N)

Find all solutions for `log_λ` given the `state` ([`P3State`](@ref)), `L`, and `N`.

!!! note "Usage"
    This function is experimental, and usually only relevant for the
    [`SlopePowerLaw`](@ref) parameterization, which can have multiple solutions
    for `log_λ` for a given `log_L` and `log_N`.
"""
function get_distribution_parameters_all_solutions(
    state::P3State{FT};
    L,
    N,
) where {FT}
    # Find bounds by evaluating function incrementally, then apply root finding with bounds above and below zero-point
    target_log_LdN = log(L) - log(N)

    shape_problem(log_λ) = log_L_div_N(state, log_λ) - target_log_LdN

    Δλ = 0.01
    λs = 10.0 .^ (2.0:Δλ:6.0)
    λ_bnds = Tuple[]
    # Loop over λs and find where shape_problem changes sign
    for i in 1:(length(λs) - 1)
        if shape_problem(log(λs[i])) * shape_problem(log(λs[i + 1])) < 0
            push!(λ_bnds, (λs[i], λs[i + 1]))
        end
    end

    # Apply root finding with bounds above and below zero-point
    dists = P3Distribution[]
    for (λ_min, λ_max) in λ_bnds
        log_λ_min, log_λ_max = log(λ_min), log(λ_max)
        dist = get_distribution_parameters(state; L, N, log_λ_min, log_λ_max)
        push!(dists, dist)
    end
    return dists
end

### ----------------- ###
### ----- UTILS ----- ###
### ----------------- ###

function Base.show(io::IO, p3s::P3Distribution{FT}) where {FT}
    rimed = isunrimed(p3s) ? "unrimed" : "rimed"
    println(io, "P3Distribution{$FT}")
    println(io, "├── state: is $rimed")
    println(io, "├── log_λ = $(p3s.log_λ) [log(1/m)]")
    println(io, "└── log_N₀ = $(p3s.log_N₀) [log(1/m^3)]")
end
