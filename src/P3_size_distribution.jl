
"""
    logN′ice(state, logλ)

Compute the log of the ice particle number concentration at diameter `D` 
    given the distribution `dist`
"""
function logN′ice(state::P3State, logλ)
    μ = get_μ(state, logλ)
    log_N₀ = get_logN₀(N_ice, μ, logλ)
    return function logN′(D)
        logD = log(D)
        log_N₀ + μ * logD - exp(logλ + logD)
    end
end

"""
    N′ice(dist, D)

Compute the ice particle number concentration at diameter `D` given the distribution `dist`
"""
N′ice(state::P3State, logλ) = exp ∘ logN′ice(state, logλ)


### ------------------------------------------------ ###
### ----- Obtaining P3 distribution parameters ----- ###
### ------------------------------------------------ ###

"""
    loggamma_inc_moment(D₁, D₂, μ, logλ, [k = 0], [scale = 1])

Compute `log(Iᵏ)` where `Iᵏ` is the following integral:

    ``I^k = ∫_{D₁}^{D₂} G(D) D^k dD``

 ``G(D) ≡ D^μ e^{-λD}`` is the (unnormalized) gamma kernel, and `k` is an arbitrary exponent.

 If `scale` is provided, `log(scale ⋅ Iᵏ)` is returned.

 With appropriate scaling, we can compute useful quantities like:
 - the `k`-th moment of the ice PSD,
    ``M^k = N₀ I^k``
 - combined power law and moment weighted integrals,
    ``∫_{D₁}^{D₂} (aD^b) D^n K(D) dD ≡ a I^(b + n)``

# Arguments
 - `D₁`: The minimum diameter [`m`]
 - `D₂`: The maximum diameter [`m`]
 - `μ`: The PSD shape parameter [`-`]
 - `logλ`: The log of the PSD slope parameter [`log(1/m)`]
 - `k`: An arbitrary exponent [`-`], default is `0`
 - `scale`: The scale factor [`-`], default is `1`

# Extended help
 ## Implementation details
 We can write `∫_D₁^D₂ G(D) D^k dD`, where `G(D) = D^μ e^{-λD}` as:
    `∫_D₁^∞ G(D) D^k dD - ∫_D₂^∞ G(D) D^k dD`
 with the transformation `x = λD`, and `z = μ+k+1`, each term can be written as:
    `∫_{Dᵢ}^∞ G(D) D^k dD = ∫_{λDᵢ}^∞ x^z e^{-x} dx / λ^z = Γ(z, λDᵢ) / λ^z`
 where `Γ(z, λDᵢ) = q ⋅ Γ(z)` and `q` is the incomplete gamma function ratio given by
    `(_, q) = SF.gamma_inc(z, x)`.
 This means that the integral `∫_{Dᵢ}^∞ G(D) D^k dD` is computed as:
    `Γ(z) ⋅ q / λ^z`
 The full integral from `D₁` to `D₂` is then:
    `Γ(z) ⋅ (q_D₁ - q_D₂) / λ^z`
 In log-space, this is:
    `- z log(λ) + logΓ(z) + log(q_D₁ - q_D₂)`
 
"""
function loggamma_inc_moment(D₁, D₂, μ, logλ, k = 0, scale = 1)
    z = k + μ + 1
    (_, q_D₁) = SF.gamma_inc(z, exp(logλ) * D₁)
    (_, q_D₂) = SF.gamma_inc(z, exp(logλ) * D₂)
    return -z * logλ + SF.loggamma(z) + log(q_D₁ - q_D₂) + log(scale)
end

"""
    loggamma_moment(μ, logλ; [k = 0], [scale = 1])

Compute `log(scale ⋅ ∫_0^∞ G(D) D^k dD)`, 
 where `G(D) ≡ D^μ e^{-λD}` is the (unnormalized) gamma kernel, 
 `k` is an arbitrary exponent, and `scale` is a scale factor.

# Arguments
 - `μ`: The PSD shape parameter [`-`]
 - `logλ`: The log of the PSD slope parameter [`log(1/m)`]

# Keyword arguments
- `k`: An arbitrary exponent [`-`], default is `0`
- `scale`: The scale factor [`-`], default is `1`.

# Extended help
 ## Implementation details
 If we write `D̃ = λD`, then the integral can be rewritten as:
 `(1/λ)^(μ+k+1) ∫_0^∞ D̃^(μ+k+1) e^{-D̃} D̃^k dD̃ = Γ(z) / λ^z`
 where `z = μ + k + 1`.

 Thus, the log of the integral is:
 `logΓ(z) - z * logλ + log(scale)`
"""
function loggamma_moment(μ, logλ; k = 0, scale = 1)
    z = k + μ + 1
    return SF.loggamma(z) - z * logλ + log(scale)
end

"""
    get_μ(slope::CMP.SlopeLaw, logλ)
    get_μ(state::P3State, logλ)
    
Compute the slope parameter μ

# Arguments
- `slope`: [`CMP.SlopeLaw`](@ref) object, or
- `state`: [`P3State`](@ref) object, or
- `params`: [`CMP.ParametersP3`](@ref) object
- `logλ`: The log of the slope parameter [log(1/m)]
"""
get_μ((; a, b, c, μ_max)::CMP.SlopePowerLaw, logλ) = clamp(a * exp(logλ)^b - c, 0, μ_max)
get_μ((; μ)::CMP.SlopeConstant, logλ...) = μ
get_μ((; params)::P3State, logλ) = get_μ(params.slope, logλ)

"""
    log∫mass_gamma_moment(state, logλ; [n=0])

Compute `log(∫_0^∞ Dⁿ m(D) N′(D) dD)` given the `state` and `logλ`.
    This is the log of the `n`-th moment of the mass-weighted PSD.

# Arguments
- `state`: [`P3State`](@ref) object
- `logλ`: The log of the slope parameter [log(1/m)]
- `n`: The order of the moment [dimensionless]

# Note:
- For `n = 0`, this evaluates to `log(L/N₀)`, see [`logLdN₀`](@ref)
- For `n = 1`, this evaluates to the (unnormalized) mass-weighted mean particle size, see [`D_m`](@ref)
"""
@inline function log∫mass_gamma_moment(state::P3State, logλ; n = 0)
    mass_law = get_mass_law(state)
    μ = get_μ(state, logλ)
    
    return LogExpFunctions.logsumexp(
        let (D_min, D_max) = CO.get_bounds(piece), (a, b) = CO.get_coefficients(piece)
            loggamma_inc_moment(D_min, D_max, μ, logλ, b + n, a)
        end for piece in mass_law.pieces
    )
end

"""
    logLdivN₀(state, logλ)

Compute `log(L/N₀)` given the parameters that define the state, and `logλ`
"""
logLdivN₀(state::P3State, logλ) = log∫mass_gamma_moment(state, logλ; n = 0)

"""
    logNdivN₀(state, logλ)

Compute `log(N/N₀)` given the `state`, and `logλ`
"""
logNdivN₀(state::P3State, logλ) = loggamma_moment(0, get_μ(state, logλ), logλ)

"""
    logLdivN(state, logλ)

Compute `log(L/N)` given the `state` and `logλ`

# Arguments
- `state`: [`P3State`](@ref) object
- `logλ`: The log of the slope parameter [log(1/m)]
"""
logLdivN(state::P3State, logλ) = logLdivN₀(state, logλ) - logNdivN₀(state, logλ)

"""
    get_logN₀(N_ice, μ, logλ)

Compute `log(N₀)` given the `state`, `N`, and `logλ`,

        N  = N₀ ∫ G(D) dD
    log N₀ = log N - log(∫G(D) dD) 
           = log(N) - log( ∫D^μ e^{-λD} dD )
           = log(N) - M⁰

# Arguments
- `N_ice`: The number concentration [1/m³]
- `μ`: The slope parameter
- `logλ`: The log of the slope parameter [log(1/m)]

# Keyword arguments
- `N`: The number concentration [1/m³]
- `logλ`: The log of the slope parameter [log(1/m)]
"""
get_logN₀(N_ice, μ, logλ) = log(N_ice) - loggamma_moment(0, μ, logλ)

"""
    get_distribution_logλ(params, L_ice, N_ice, F_rim, ρ_rim; [logλ_min, logλ_max])
    get_distribution_logλ(state; [logλ_min, logλ_max])

Solve for the distribution parameters given the state, and the mass (`L`) and number (`N`) concentrations.

The assumed distribution is of the form

```math
N′(D) = N₀ D^μ e^{-λD}
```
where `N′(D)` is the number concentration at diameter `D` and `μ` is the slope parameter.
    The slope parameter is parameterized, e.g. [`CMP.SlopePowerLaw`](@ref) or [`CMP.SlopeConstant`](@ref).

This algorithm solves for `logλ = log(λ)` and `log_N₀ = log(N₀)` 
    given `L_ice` and `N_ice` by solving the equations:

```math
\\begin{align*}
\\log(L) &= \\log ∫_0^∞ m(D) N′(D)\\ \\mathrm{d}D, \\\\
\\log(N) &= \\log ∫_0^∞ N′(D)\\ \\mathrm{d}D, \\\\
\\end{align*}
```
where `m(D)` is the mass of a particle at diameter `D` (see [`ice_mass`](@ref)).
    The procedure is decribed in detail in [the P3 docs](@ref "Parameterizations for the slope parameter \$μ\$").


# Arguments
- `state`: [`P3State`](@ref) object
- `params`: [`CMP.ParametersP3`](@ref) object
- `L_ice`: The mass concentration [kg/m³]
- `N_ice`: The number concentration [1/m³]
- `F_rim`: The rime mass fraction
- `ρ_rim`: The rime density

# Keyword arguments
- `logλ_min`: The minimum value of the search bounds [log(1/m)], default is `log(1e1)`
- `logλ_max`: The maximum value of the search bounds [log(1/m)], default is `log(1e7)`
"""
function get_distribution_logλ(
    state::P3State{FT}; logλ_min = log(1e1), logλ_max = log(1e7)
) where {FT}
    target_log_LdN = log(state.L_ice) - log(state.N_ice)

    shape_problem(logλ) = logLdN(state, logλ) - target_log_LdN

    # Find slope parameter
    sol = RS.find_zero(
        shape_problem,
        RS.SecantMethod(FT(logλ_min), FT(logλ_max)),
        RS.CompactSolution(),
        RS.RelativeSolutionTolerance(eps(FT)),
        50,
    )
    return sol.root  # logλ
end
function get_distribution_logλ(params::CMP.ParametersP3, L_ice, N_ice, F_rim, ρ_rim; kwargs...)
    state = get_state(params; L_ice, N_ice, F_rim, ρ_rim)
    return get_distribution_logλ(state; kwargs...)
end

"""
    get_distribution_parameters_all_solutions(state; L, N)

Find all solutions for `logλ` given the `state` ([`P3State`](@ref)), `L`, and `N`.

!!! note "Usage"
    This function is experimental, and usually only relevant for the
    [`SlopePowerLaw`](@ref) parameterization, which can have multiple solutions
    for `logλ` for a given `log_L` and `log_N`.
"""
function get_distribution_parameters_all_solutions(state::P3State)
    # Find bounds by evaluating function incrementally, then apply root finding with bounds above and below zero-point
    target_log_LdN = log(state.L_ice) - log(state.N_ice)

    shape_problem(logλ) = logLdN(state, logλ) - target_log_LdN

    Δλ = 0.01
    λs = 10.0 .^ (2.0:Δλ:6.0)
    logλ_bnds = Tuple[]
    # Loop over λs and find where shape_problem changes sign
    for i in 1:(length(λs) - 1)
        if shape_problem(log(λs[i])) * shape_problem(log(λs[i + 1])) < 0
            push!(logλ_bnds, (log(λs[i]), log(λs[i + 1])))
        end
    end

    # Apply root finding with bounds above and below zero-point
    logλs = [get_distribution_logλ(state; logλ_min, logλ_max) for (logλ_min, logλ_max) in logλ_bnds]
    return logλs
end
