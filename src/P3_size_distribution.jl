import CloudMicrophysics.DistributionTools: size_distribution

"""
    logN′ice(state, logλ)

Compute the log of the ice particle number concentration at diameter `D` given the distribution `dist`
"""
function logN′ice(state::P3State, logλ)
    μ = get_μ(state, logλ)
    log_N₀ = get_logN₀(state.ρn_ice, μ, logλ)
    return function logN′(D)
        logD = log(D)
        log_N₀ + μ * logD - exp(logλ + logD)
    end
end

"""
    size_distribution(state::P3State, logλ)

Return `n(D)`, a function that computes the size distribution for ice particles at diameter `D`

# Arguments
- `state`: The [`P3State`](@ref)
- `logλ`: The log of the slope parameter [log(1/m)]
"""
DT.size_distribution(state::P3State, logλ) = exp ∘ logN′ice(state, logλ)

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

 !!! note "Log/linear twin"
     [`gamma_inc_moment`](@ref) (just below) is the **linear-space**
     sibling of this function: it returns the same incomplete-gamma
     partial moment `Γ(z)(q_D₁−q_D₂)/λ^z` *without* the `log`. The
     closed-form rain-inner collision split (`closed_rain_inner_NM`,
     `P3_processes.jl`) needs the linear form because it accumulates
     **sign-alternating** addends (the `±sgn` split pieces and the
     `vᵢ_part − v_part` subtraction) that cannot be done under a `log`
     / `logsumexp`. The two are kept adjacent so the pair stays in
     sync (P1-2). `Common.Chen2022_exponential_pdf` is the *full-domain*
     `[0,∞)` 1M analogue (complete, not incomplete, gamma) — not a
     drop-in for the closed form, which needs *partial* `[D₁,D₂]`
     moments because of the velocity-crossover `D*` split.

"""
function loggamma_inc_moment(D₁, D₂, μ, logλ, k = 0, scale = 1)
    FT = float(UT.promote_typeof(D₁, D₂, μ, logλ))
    D₁ < D₂ || return log(FT(0))  # return log(0) if D₁ ≥ D₂
    z = k + μ + 1
    # NOTE: We use `LogExpFunctions.xexpy(D, logλ)` to compute `λD = D * exp(logλ)`.
    # When `logλ` is large, `exp(logλ) = Inf`, so the naive product `D * exp(logλ)`
    # yields `0 * Inf = NaN` when `D = 0`. `xexpy` correctly returns `0` in that case.
    (_, q_D₁) = SF.gamma_inc(z, LogExpFunctions.xexpy(D₁, logλ))
    (_, q_D₂) = SF.gamma_inc(z, LogExpFunctions.xexpy(D₂, logλ))
    return -z * logλ + SF.loggamma(z) + log(q_D₁ - q_D₂) + log(FT(scale))
end

"""
    gamma_inc_moment(D₁, D₂, p, α)

`∫_{D₁}^{D₂} D^p e^{-α D} dD = α^{-(p+1)} Γ(p+1) [Q(p+1,αD₁) − Q(p+1,αD₂)]`
with `Q` the regularized upper incomplete gamma.

The **linear-space twin** of [`loggamma_inc_moment`](@ref) (same
`Γ(z)(q_{D₁}−q_{D₂})/λ^z`, `z = p+1` reduction, without the `log` and
without `loggamma_inc_moment`'s `logλ`/`xexpy` overflow guard, which is
unneeded here — the closed-form slopes `α0=1/D̄r`, `αⱼ=α0+cⱼ` are plain
positive scalars well inside range). Used by the closed-form rain-inner
collision split (`closed_rain_inner_NM`, `P3_processes.jl`), whose
**sign-alternating** addends cannot be accumulated under a `log` /
`logsumexp` — see the "Log/linear twin" note on
[`loggamma_inc_moment`](@ref). Keep the two adjacent and in sync (P1-2).

`Common.Chen2022_exponential_pdf` is the *full-domain* `[0,∞)` 1M
analogue (complete gamma); it is **not** a drop-in here because the
velocity-crossover `D*` split needs *partial* `[D₁,D₂]` (incomplete)
moments.

Returns `0` if `D₂ ≤ D₁`. Returns `NaN` if `α ≤ 0` (a mis-set slope);
for the dispatched `(SB2006, Chen)` bundle `α > 0` is an invariant
(`D̄r > 0`, Chen rain `cⱼ > 0`) so that branch is unreachable, and
`get_liquid_integrals_rain_closed` folds any non-finite result to `0`
as a defensive belt-and-suspenders (P1-1). `z = p+1 > 0` is likewise
guaranteed (Chen rain `Bⱼ > 0`, `bounds_r` lower > 0) so `SF.gamma(z)`
is always finite (P2-2). `p+1` is a Float64 constant — the `Dual`
enters only via `gamma_inc`'s **2nd** argument (AD-safe, A2).
"""
@inline function gamma_inc_moment(D₁, D₂, p, α)
    # AD-generic working type: promotes Float64 (e.g. a Float64 D₁/D₂
    # quadrature endpoint) and a `Dual` slope `α` (Dual via state/PSD
    # params) to the Dual so the partial moment differentiates. `p` is
    # a plain Float64/Int exponent (z = p+1 stays a non-Dual constant —
    # the A2-safe `gamma_inc` 2nd-arg-only pattern; never a Dual 1st arg).
    T = float(promote_type(typeof(D₁), typeof(D₂), typeof(α)))
    D₂ > D₁ || return zero(T)
    α > 0 || return T(NaN)
    z = p + 1  # > 0: Chen rain Bⱼ>0 & bounds_r lower>0 ⇒ SF.gamma(z) finite
    (_, q1) = SF.gamma_inc(z, α * D₁)
    (_, q2) = SF.gamma_inc(z, α * D₂)
    # Float64 * Dual / Dual auto-promotes to the Dual `T` (the `zero(T)`
    # / `T(NaN)` early returns are why `T` is computed explicitly).
    return SF.gamma(z) * (q1 - q2) / α^z
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
    logmass_gamma_moment(state, logλ; [n=0])

Compute `log(∫_0^∞ Dⁿ m(D) N′(D) dD)` given the `state` and `logλ`.
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

Compute `log(N₀)` given the `state`, `N`, and `logλ`,

        N  = N₀ ∫ G(D) dD
    log N₀ = log N - log(∫G(D) dD) 
           = log(N) - log( ∫D^μ e^{-λD} dD )
           = log(N) - M⁰

# Arguments
- `N_ice`: The number concentration [1/m³]
- `μ`: The shape parameter [`-`]
- `logλ`: The log of the slope parameter [log(1/m)]
"""
function get_logN₀(N_ice, μ, logλ)
    logNdivN₀ = loggamma_moment(μ, logλ; k = 0)
    logN₀ = log(N_ice) - logNdivN₀
    return logN₀
end

"""
    get_distribution_logλ(state, [logλ_guess, logλ_min, logλ_max])

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
- `state`: The [`P3State`](@ref)
- `logλ_guess`: Optional warm-start seed. When `logλ_guess` is finite and
    strictly inside `[logλ_min, logλ_max]`, one extra evaluation of the
    residual at the guess is used to **halve the Brent bracket** by sign:
    whichever endpoint sits on the same side of the root as `f(logλ_guess)`
    is replaced by the guess, and Brent then converges from the tighter
    bracket. If the guess is `nothing` / non-finite / out-of-bracket, or
    evaluates to a non-finite residual, the solver falls back to the full
    bracket with no behaviour change. Default: `nothing`.

    **Caveat — monotonicity depends on the `μ(λ)` parameterization.**

    - Under [`CMP.SlopeConstant`](@ref) (μ fixed, independent of λ),
      `log(L/N)(logλ)` is a strictly decreasing function of `logλ` —
      larger slope ⇒ smaller mean mass. Brent finds the unique root on
      any valid bracket and the sign-based warm-start narrowing is
      exact.
    - Under [`CMP.SlopePowerLaw`](@ref) (μ = clamp(a·λ^b - c, 0, μ_max)
      — piecewise flat, rising, flat), the residual is *not* globally
      monotonic and the same target `L/N` can have several roots (see
      `docs/src/plots/P3SlopeParameterizations.jl` for a three-root
      example). In that regime the halved bracket still contains at
      least one root, but it may not be the root closest to the guess;
      the warm-start is heuristic rather than exact. For smooth cell
      evolution (step-to-step `logλ` continuity) this is usually fine
      because the guess and the current root stay on the same branch,
      but callers that need all roots should use
      `get_distribution_logλ_all_solutions`.
- `logλ_min`: The minimum value of the search bounds [log(1/m)], default is `log(1e1)`
- `logλ_max`: The maximum value of the search bounds [log(1/m)], default is `log(1e7)`
"""
function get_distribution_logλ(state, logλ_guess = nothing, logλ_min = 2, logλ_max = 17)
    FT = eltype(state)
    ϵₘ = UT.ϵ_numerics_2M_M(FT)
    ϵₙ = UT.ϵ_numerics_2M_N(FT)
    (; ρn_ice, ρq_ice) = state
    (ρn_ice < ϵₙ || ρq_ice < ϵₘ) && return log(zero(ρq_ice))
    target_log_LdN = log(ρq_ice) - log(ρn_ice)

    shape_problem(logλ) = logLdivN(state, logλ) - target_log_LdN

    # Cold bracket. Brent converges to *some* root inside a valid bracket.
    # Whether that root is unique over the full bracket depends on the
    # `μ(λ)` parameterization: monotone for `SlopeConstant`, possibly
    # multi-valued for `SlopePowerLaw` (see docstring caveat). If either
    # endpoint residual is non-finite or the endpoints share a sign, the
    # target `L/N` lies outside the representable PSD range — saturate at
    # the nearer endpoint rather than failing.
    lo, hi = FT(logλ_min), FT(logλ_max)
    f_lo, f_hi = shape_problem(lo), shape_problem(hi)
    if !isfinite(f_lo) || !isfinite(f_hi) || f_lo * f_hi > 0
        return abs(f_lo) ≤ abs(f_hi) ? lo : hi
    end

    # Optional hot-start: one probe at the prior step's `logλ` halves the
    # bracket by sign, and Brent takes over from there. Any further
    # narrowing is Brent's job — its interpolating iteration is more
    # efficient per f-eval than any fixed-offset probe we could add here.
    # Non-monotonic states may have the halved bracket land on a different
    # root than the guess; see the docstring caveat.
    (lo, f_lo, hi, f_hi) =
        _narrow_bracket(shape_problem, lo, f_lo, hi, f_hi, logλ_guess)

    sol = RS.find_zero(
        shape_problem,
        RS.BrentsMethod(lo, hi),
        RS.CompactSolution(),
        RS.RelativeSolutionTolerance(eps(FT)),
        100,
    )
    return sol.root  # logλ
end

"""
    get_distribution_logλ_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)

Compute `log(λ)` for P3, using prognostic ice variables directly

The P3 variables `F_rim` and `ρ_rim` are computed in a regularised way
"""
function get_distribution_logλ_from_prognostic(
    params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, args...,
)
    state = state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    return get_distribution_logλ(state, args...)
end

# One-probe narrowing of a valid bracket `[lo, hi]` using a probe point `p`.
# No-op if `p` is not usable (nothing / non-finite / outside the bracket) or
# if `f(p)` is non-finite. Otherwise, replace whichever endpoint sits on the
# same side of the root as `p` — a single sign check on `f_lo * f_p`.
@inline _narrow_bracket(_sp, lo, f_lo, hi, f_hi, ::Nothing) = (lo, f_lo, hi, f_hi)
@inline function _narrow_bracket(shape_problem, lo, f_lo, hi, f_hi, p::Real)
    p_ = oftype(lo, p)  # match FT of the bracket for type stability
    (isfinite(p_) && lo < p_ < hi) || return (lo, f_lo, hi, f_hi)
    f_p = shape_problem(p_)
    isfinite(f_p) || return (lo, f_lo, hi, f_hi)
    return f_lo * f_p < 0 ? (lo, f_lo, p_, f_p) : (p_, f_p, hi, f_hi)
end

"""
    get_distribution_logλ_all_solutions(state)

Find all solutions for `logλ` given the `state` ([`P3State`](@ref)), `L`, and `N`.

!!! note "Usage"
    This function is experimental, and usually only relevant for the
    [`SlopePowerLaw`](@ref) parameterization, which can have multiple solutions
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
