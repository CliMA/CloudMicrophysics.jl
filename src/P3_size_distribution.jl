import CloudMicrophysics.DistributionTools: size_distribution

# Callable returned by `logN′ice`: evaluates `log(N′(D))` for a fixed state and slope.
# We store `λ = exp(logλ)` (computed once when the functor is built) rather than `logλ`
# so the slope term is a single multiply `λ * D` in the quadrature hot loop, instead of
# `exp(logλ + logD)` (a transcendental per evaluation). `logD = log(D)` is still needed
# for the `μ * logD` term.
struct P3LogNumberFunctor{FT} <: Function
    log_N₀::FT
    μ::FT
    λ::FT
end
@inline function (f::P3LogNumberFunctor)(D)
    logD = log(D)
    return f.log_N₀ + f.μ * logD - f.λ * D
end

"""
    logN′ice(state, logλ)

Return a callable that computes `log(N′(D))`, the log of the ice particle number
concentration at diameter `D`, for the [`P3State`](@ref) `state` and log-slope `logλ`.
"""
function logN′ice(state::P3State, logλ)
    μ = get_μ(state, logλ)
    (; ρn_ice) = state
    # Without ice number, the distribution is zero (log(N₀) = -Inf). The substitute `ρn_ice = 1`
    # keeps `get_logN₀` and its ForwardDiff partials finite in the discarded branch.
    present = FD.value(ρn_ice) > zero(FD.value(ρn_ice))
    ρn_ice_safe = ifelse(present, ρn_ice, one(ρn_ice))
    log_N₀_safe = get_logN₀(ρn_ice_safe, μ, logλ)
    log_N₀ = ifelse(present, log_N₀_safe, oftype(log_N₀_safe, -Inf))
    # Promote to a common type: differentiating w.r.t. the ice number makes
    # `log_N₀` a `Dual` while `μ` (a function of the fixed `logλ`) stays a plain
    # float, and `P3LogNumberFunctor` stores both in a single field type.
    λ = exp(logλ)
    return P3LogNumberFunctor(promote(log_N₀, μ, λ)...)
end

# Callable returned by `size_distribution`: `n(D) = exp(logN′(D))`.
struct P3SizeDistributionFunctor{F} <: Function
    logN′::F
end
@inline (f::P3SizeDistributionFunctor)(D) = exp(f.logN′(D))

"""
    size_distribution(state::P3State, logλ)

Return `n(D)`, a function that computes the size distribution for ice particles at diameter `D`

# Arguments
- `state`: The [`P3State`](@ref)
- `logλ`: The log of the slope parameter [log(1/m)]
"""
DT.size_distribution(state::P3State, logλ) = P3SizeDistributionFunctor(logN′ice(state, logλ))

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
    `(_, q) = UT.gamma_inc(z, x)`.
 This means that the integral `∫_{Dᵢ}^∞ G(D) D^k dD` is computed as:
    `Γ(z) ⋅ q / λ^z`
 The full integral from `D₁` to `D₂` is then:
    `Γ(z) ⋅ (q_D₁ - q_D₂) / λ^z`
 In log-space, this is:
    `- z log(λ) + logΓ(z) + log(q_D₁ - q_D₂)`

See also [`gamma_inc_moment`](@ref)
"""
function loggamma_inc_moment(D₁, D₂, μ, logλ, k = 0, scale = 1)
    FT = UT.promote_typeof(D₁, D₂, μ, logλ)
    D₁ < D₂ || return log(FT(0))  # return log(0) if D₁ ≥ D₂
    z = k + μ + 1
    # `λ⋅D ≡ xexpy(D, logλ) ≡ D * exp(logλ)` (numerically stable)
    x1 = LogExpFunctions.xexpy(D₁, logλ)
    x2 = LogExpFunctions.xexpy(D₂, logλ)
    (p1, q1) = UT.gamma_inc(z, x1)
    (p2, q2) = UT.gamma_inc(z, x2)
    Δq = x2 < z + 1 ? p2 - p1 : q1 - q2
    # `Δq ≥ 0` analytically but can round below zero. The floor keeps `log(Δq)` and its
    # derivative finite.
    Δq = max(Δq, floatmin(FT))
    return -z * logλ + SF.loggamma(z) + log(Δq) + log(FT(scale))
end

"""
    gamma_inc_moment(D₁, D₂, p, α)

`∫_{D₁}^{D₂} D^p e^{-α D} dD = α^{-(p+1)} Γ(p+1) [Q(p+1,αD₁) - Q(p+1,αD₂)]`
with `Q` the regularized upper incomplete gamma.

Returns `0` if `D₂ ≤ D₁`, and `NaN` if `α ≤ 0`

See also [`loggamma_inc_moment`](@ref)
"""
@inline function gamma_inc_moment(D₁, D₂, p, α)
    FT = float(promote_type(typeof(D₁), typeof(D₂), typeof(α)))
    D₂ > D₁ || return zero(FT)
    α > 0 || return FT(NaN)
    z = p + 1
    x1 = α * D₁
    x2 = α * D₂
    (p1, q1) = UT.gamma_inc(z, x1)
    (p2, q2) = UT.gamma_inc(z, x2)
    Δq = x2 < z + 1 ? p2 - p1 : q1 - q2
    Δq = max(Δq, zero(FT))
    return SF.gamma(z) * Δq / α^z
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

The implementation follows the same logic as [`loggamma_inc_moment`](@ref),
    but with `D₁ = 0` and `D₂ = ∞`, which implies `q_D₁ = 1` and `q_D₂ = 0`.
"""
function loggamma_moment(μ, logλ; k = 0, scale = 1)
    FT = eltype(μ)
    z = k + μ + 1
    return -z * logλ + SF.loggamma(z) + log(FT(scale))
end

"""
    get_μ(slope::CMP.SlopeLaw, logλ)
    get_μ(state::P3State, logλ)
    
Compute the shape parameter μ

# Arguments
- `slope`: [`CMP.SlopeLaw`](@ref) object, or
- `state`: [`P3State`](@ref) object
- `logλ`: The log of the slope parameter [log(1/m)]
"""
get_μ((; a, b, c, μ_max)::CMP.SlopePowerLaw, logλ) = clamp(a * exp(logλ)^b - c, 0, μ_max)
function get_μ((; a, b, c, μ_max, κ)::CMP.SmoothSlopePowerLaw, logλ)
    p = a * exp(logλ)^b - c
    M₀ = LogExpFunctions.log1pexp(κ * p) / κ
    return μ_max - LogExpFunctions.log1pexp(κ * (μ_max - M₀)) / κ
end
get_μ((; μ)::CMP.SlopeConstant, logλ...) = μ
get_μ((; params)::P3State, logλ) = get_μ(params.slope, logλ)

"""
    logmass_gamma_moment(state::P3State, μ, logλ; n = 0)

Compute `log(∫_0^∞ Dⁿ m(D) N′(D) dD)` given the `state`, `μ`, and `logλ`.
    This is the log of the `n`-th moment of the mass-weighted PSD.

# Arguments
- `state`: [`P3State`](@ref) object
- `μ`: The shape parameter [`-`]
- `logλ`: The log of the slope parameter [log(1/m)]

# Keyword arguments
- `n`: The order of the moment, default is `0`

# Note:
- For `n = 0`, this evaluates to `log(L/N₀)`
- For `n = 1`, this evaluates to the (unnormalized) mass-weighted mean particle size, see [`D_m`](@ref)
"""
function logmass_gamma_moment(state::P3State, μ, logλ; n = 0)
    bnds = segment_boundaries(state)
    moments = UU.unrolled_map(subintervals(bnds)) do (D_lo, D_hi)
        (a, b) = ice_mass_coeffs(state, (D_lo + D_hi) / 2)
        loggamma_inc_moment(D_lo, D_hi, μ, logλ, b + n, a)
    end
    return UT.unrolled_logsumexp(moments)
end

"""
    logLdivN(state, logλ)

Compute `log(L/N)` given the `state` and `logλ`

# Arguments
- `state`: [`P3State`](@ref) object
- `logλ`: The log of the slope parameter [log(1/m)]
"""
function logLdivN(state::P3State, logλ)
    μ = get_μ(state, logλ)
    logLdivN₀ = logmass_gamma_moment(state, μ, logλ; n = 0)
    logNdivN₀ = loggamma_moment(μ, logλ; k = 0)
    return logLdivN₀ - logNdivN₀
end

"""
    get_logN₀(N_ice, μ, logλ)

Compute `log(N₀)` given `N_ice`, `μ`, and `logλ`,

        N  = N₀ ∫ G(D) dD
    log N₀ = log N - log(∫G(D) dD)
           = log(N) - log( ∫D^μ e^{-λD} dD )
           = log(N) - log(M⁰)

# Arguments
- `N_ice`: The number concentration [1/m³]
- `μ`: The shape parameter [`-`]
- `logλ`: The log of the slope parameter [log(1/m)]

Requires `N_ice > 0`.
"""
function get_logN₀(N_ice, μ, logλ)
    logNdivN₀ = loggamma_moment(μ, logλ; k = 0)
    return log(N_ice) - logNdivN₀
end

"""
    FixedIterations{FT}()

A `RootSolvers.AbstractTolerance` whose convergence predicate is always `false`,
so the bracketing solver never exits early and always runs the full iteration
budget. This makes the iteration count independent of the input, eliminating
warp divergence from data-dependent early-exit on the GPU (at the cost of the
warm-start speedup - a tighter initial bracket improves accuracy but not the
iteration count). Each caller sets its own iteration budget.
"""
struct FixedIterations{FT} <: RS.AbstractTolerance{FT} end
@inline (::FixedIterations)(x1, x2, y) = false

"""
    _mean_mass_target_logλ(state, a, b, target_mass)

Solve `log(target_mass) = log(a) + loggamma(b+μ+1) - loggamma(μ+1) - b·logλ` for `logλ`, where
`μ = get_μ(state, logλ)`, by a fixed-point iteration with a fixed step count.

`(a, b)` are the coefficients of the ice mass power law `a D^b` in the regime that contains
`target_mass`. The solution is exact when the whole distribution lies in that regime.
"""
@inline function _mean_mass_target_logλ(state, a, b, target_mass)
    FT = eltype(state)
    log_target = log(target_mass)
    logλ = (log(a) + SF.loggamma(b + 1) - log_target) / b  # μ = 0 seed
    for _ in 1:3
        μ = get_μ(state, logλ)
        logλ = (log(a) + SF.loggamma(b + μ + 1) - SF.loggamma(μ + 1) - log_target) / b
    end
    return logλ
end

"""
    _derived_logλ_bracket(state)

Default search bounds `(logλ_min, logλ_max)` of [`get_distribution_logλ`](@ref).

`logλ_max` is the slope at which the mean particle mass equals
[`ice_mean_particle_mass_min`](@ref). Such small particles are small spherical ice, so
`logλ_max` does not depend on `F_rim` or `ρ_rim`. `logλ_min = 2`.
"""
@inline function _derived_logλ_bracket(state::P3State)
    FT = eltype(state)
    p3 = state.params
    mass_min = ice_mean_particle_mass_min(p3)
    (a_small, b_small) = ice_mass_coeffs(state, zero(FT))
    logλ_max = _mean_mass_target_logλ(state, a_small, b_small, mass_min)
    margin = FT(100) * eps(FT)  # rounding margin
    return (FT(2), logλ_max + margin)
end

"""
    get_distribution_logλ(state, [logλ_guess, logλ_min, logλ_max])

Solve for the distribution parameters given the state, and the mass (`L`) and number (`N`) concentrations.

The assumed distribution is of the form

```math
N′(D) = N₀ D^μ e^{-λD}
```
where `N′(D)` is the number concentration at diameter `D` and `μ` is the shape parameter.
    The shape parameter is parameterized, e.g. [`CMP.SlopePowerLaw`](@ref) or [`CMP.SlopeConstant`](@ref).

This algorithm solves for `logλ = log(λ)`
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
- `state`: The [`P3State`](@ref)
- `logλ_guess`: Optional initial guess, which narrows the search bounds but does not change
  the iteration count
- `logλ_min`, `logλ_max`: The search bounds [log(1/m)]. By default, `logλ_min = 2`, and
  `logλ_max` is the slope at which the mean particle mass equals the mass of a newly nucleated
  crystal.
"""
function get_distribution_logλ(state, logλ_guess = nothing, logλ_min = nothing, logλ_max = nothing)
    FT = eltype(state)
    (; ρn_ice, ρq_ice) = state
    lo, hi = if logλ_min === nothing || logλ_max === nothing
        (dlo, dhi) = _derived_logλ_bracket(state)
        (logλ_min === nothing ? dlo : FT(logλ_min), logλ_max === nothing ? dhi : FT(logλ_max))
    else
        (FT(logλ_min), FT(logλ_max))
    end
    # Without ice mass, the limit of vanishing mean particle mass is the smallest particles, `hi`
    if ρq_ice <= 0 && ρn_ice > 0
        return hi
    end

    # Two logs instead of the log of the ratio, which can underflow to zero in Float32
    target_log_LdN = log(ρq_ice) - log(ρn_ice)

    shape_problem(logλ) = logLdivN(state, logλ) - target_log_LdN
    f_lo, f_hi = shape_problem(lo), shape_problem(hi)
    if !isfinite(f_lo) || !isfinite(f_hi) || f_lo * f_hi > 0
        # `logLdivN` decreases in `logλ`, so an infinite target selects a bound by its sign:
        # vanishing mean mass (-Inf) gives `hi`, no number (+Inf) gives `lo`, and 0/0 (NaN) gives `hi`
        if !isfinite(target_log_LdN)
            isnan(target_log_LdN) && return hi
            return target_log_LdN > 0 ? lo : hi
        end
        return abs(f_lo) ≤ abs(f_hi) ? lo : hi
    end
    (lo, f_lo, hi, f_hi) =
        _narrow_bracket(shape_problem, lo, f_lo, hi, f_hi, logλ_guess)

    # A fixed iteration count (no early exit) keeps GPU warps convergent. The counts reach the
    # rounding floor of `logLdivN` at each precision; `test/p3_tests.jl` asserts the residual.
    maxiters = FT === Float32 ? 10 : 12
    sol = RS.find_zero(
        shape_problem,
        RS.BrentsMethod(lo, hi),
        RS.CompactSolution(),
        FixedIterations{FT}(),
        maxiters,
    )
    return clamp(sol.root, lo, hi)  # logλ, within the search bounds
end

"""
    get_distribution_logλ_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, args...)

Compute `log(λ)` for P3, using prognostic ice variables directly.
Trailing `args...` are forwarded to [`get_distribution_logλ`](@ref).

The P3 variables `F_rim` and `ρ_rim` are computed in a regularised way
"""
function get_distribution_logλ_from_prognostic(
    params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, args...,
)
    state = state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    return get_distribution_logλ(state, args...)
end

@inline _narrow_bracket(_sp, lo, f_lo, hi, f_hi, ::Nothing) = (lo, f_lo, hi, f_hi)
@inline function _narrow_bracket(shape_problem, lo, f_lo, hi, f_hi, p::Real)
    p_ = oftype(lo, p)
    valid = isfinite(p_) & (lo < p_ < hi)
    p_clean = ifelse(valid, p_, lo)
    f_p = shape_problem(p_clean)
    valid &= isfinite(f_p)

    left = valid & (f_lo * f_p < 0)
    right = valid & !left

    new_hi = ifelse(left, p_clean, hi)
    new_f_hi = ifelse(left, f_p, f_hi)
    new_lo = ifelse(right, p_clean, lo)
    new_f_lo = ifelse(right, f_p, f_lo)

    return (new_lo, new_f_lo, new_hi, new_f_hi)
end

"""
    get_distribution_logλ_all_solutions(state)

Find all solutions for `logλ` given the `state` ([`P3State`](@ref)), `L`, and `N`.

!!! note "Usage"
    This function is experimental, and usually only relevant for the
    [`CMP.SlopePowerLaw`](@ref) parameterization, which can have multiple solutions
    for `logλ` for a given `log_L` and `log_N`.
"""
function get_distribution_logλ_all_solutions(state::P3State)
    # Find bounds by evaluating function incrementally, then apply root finding with bounds above and below zero-point
    target_log_LdN = log(state.ρq_ice) - log(state.ρn_ice)

    shape_problem(logλ) = logLdivN(state, logλ) - target_log_LdN

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
    logλs = [get_distribution_logλ(state, nothing, logλ_min, logλ_max) for (logλ_min, logλ_max) in logλ_bnds]
    return logλs
end
