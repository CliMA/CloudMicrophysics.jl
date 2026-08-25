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
    # `get_logN₀` requires a present `ρn_ice`; a benign substitute is used
    # before calling it, and the result is replaced with the doctrine-correct
    # absent answer (`log(N₀) = -Inf`, an identically-zero distribution)
    # afterward, matching `guarded_quotient`'s own construction so nothing
    # non-finite is built from a live `ForwardDiff` partial at the common
    # `ρn_ice = 0` state. A floored-but-finite surrogate (`floatmin(FT)`, as
    # used elsewhere on this branch) is wrong here specifically: unlike a
    # target that a downstream root-find bracket or `logsumexp` only reads by
    # sign or magnitude-order, `log_N₀` is consumed directly as the
    # distribution's amplitude, so a finite surrogate manufactures a real,
    # if minuscule, nonzero population instead of leaving it absent.
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
    # Cancellation, not smallness: `Δq` is a difference of two regularized
    # incomplete gamma evaluations and is nonnegative analytically, but can
    # round slightly below zero. Clamped to `floatmin(FT)` rather than to
    # zero: `log` of an exact-zero `Dual` gives a zero-valued partial times an
    # infinite derivative, `0 * Inf = NaN`, where `floatmin(FT)` keeps
    # `log`'s derivative finite (`1/floatmin(FT)` does not overflow at either
    # precision) while remaining negligible next to any populated segment's
    # moment in the `unrolled_logsumexp` this feeds.
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

# Analytic ForwardDiff partials of `M(D₁, D₂, p, α) = ∫_{D₁}^{D₂} Dᵖ e^{-αD} dD`:
# `∂M/∂D₁ = -D₁ᵖ e^{-αD₁}`, `∂M/∂D₂ = D₂ᵖ e^{-αD₂}`, `∂M/∂α = -M(D₁, D₂, p+1, α)`.
# The value is taken from the plain-`Real` method on the argument values.
@inline function gamma_inc_moment(D₁, D₂, p, α::FD.Dual{T}) where {T}
    v₁ = FD.value(T, D₁)
    v₂ = FD.value(T, D₂)
    vα = FD.value(T, α)
    M = gamma_inc_moment(v₁, v₂, p, vα)
    ∂D₁ = -v₁^p * exp(-vα * v₁)
    ∂D₂ = v₂^p * exp(-vα * v₂)
    ∂α = -gamma_inc_moment(v₁, v₂, p + 1, vα)
    Z = zero(FD.partials(α))
    part = ∂D₁ * _moment_partials(T, Z, D₁) + ∂D₂ * _moment_partials(T, Z, D₂) + ∂α * FD.partials(α)
    return FD.Dual{T}(M, part)
end

# Partials of an argument with respect to tag `T`; a non-`Dual` argument contributes `Z`.
@inline _moment_partials(::Type{T}, Z, x::FD.Dual{T}) where {T} = FD.partials(x)
@inline _moment_partials(::Type{T}, Z, x) where {T} = Z

"""
    gamma_inc_Q_chain(zf0, x, invΓ)

Return `(Q(zf0,x), Q(zf0+1,x), ..., Q(zf0+5,x))`, the regularized upper
incomplete gamma at six consecutive orders starting at `zf0`, via at most one
`gamma_inc` evaluation at `zf0` and the upward recurrence
`Q(p+1,x) = Q(p,x) + x^p e^{-x} / Γ(p+1)` for the remaining five orders.
`invΓ` holds `1/Γ(zf0+1), ..., 1/Γ(zf0+5)` (from
[`gamma_inc_moment_channel_setup`](@ref)); the boundary term's `x^p` shares
one `log(x)` across all five steps.

At `zf0 == 1` no `gamma_inc` evaluation is needed at all: `Q(1,x) = e^{-x}`
exactly, which the recurrence's own boundary factor already computes. The ice
channel runs at `p0 = 0` and so takes that path; the rain channels carry a
non-integer `zf0` and take the general one.
"""
@inline function gamma_inc_Q_chain(zf0, x, invΓ)
    logx = log(x)
    ex = exp(-x)
    # Q(1, x) = e^{-x} exactly, so at zf0 == 1 the chain's own boundary term already holds the
    # answer and the iterative evaluation is redundant. The ice channel runs at p0 = 0, hence
    # zf0 = 1, so this is the ice path; the rain channels carry zf0 = b + 1 with b a Chen exponent
    # and take the general branch. The predicate is a CHANNEL parameter, uniform across the cells
    # in a warp, so it does not divide a warp the way a state-dependent gate would.
    q0 = zf0 == oneunit(zf0) ? ex : UT.gamma_inc(zf0, x)[2]
    q1 = q0 + exp(zf0 * logx) * ex * invΓ[1]
    q2 = q1 + exp((zf0 + 1) * logx) * ex * invΓ[2]
    q3 = q2 + exp((zf0 + 2) * logx) * ex * invΓ[3]
    q4 = q3 + exp((zf0 + 3) * logx) * ex * invΓ[4]
    q5 = q4 + exp((zf0 + 4) * logx) * ex * invΓ[5]
    return (q0, q1, q2, q3, q4, q5)
end

"""
    gamma_inc_moment_channel_setup(p0, α, D_min, D_max)

Precompute, for one collision-rate channel (six consecutive moment orders
`p0, ..., p0 + 5` at rate `α`), the pieces of the crossing-split moment that
do not depend on the crossing diameter: the `Q`-chain reciprocal-Γ
denominators, the six per-order scale factors `Γ(p0+k+1) / α^(p0+k+1)`, and
the `Q`-chains at the two fixed endpoints `D_min`, `D_max`. Pass the result to
[`gamma_inc_moment_channel_finish`](@ref) once the crossing diameter is known.

An `Integer` order `p0` routes `α^(p0+1)` through `Base.power_by_squaring`
instead of the general real-exponent path, and `Γ(p0+1)` through `UT.fac`
(`= p0!`) instead of `SF.gamma`, since `SF.gamma` on small integers reads a
host-memory factorial table.
"""
@inline function gamma_inc_moment_channel_setup(p0, α, D_min, D_max)
    z0 = p0 + 1
    FT = float(promote_type(typeof(z0), typeof(α), typeof(D_min), typeof(D_max)))
    return _gamma_inc_moment_channel_setup(FT(z0), z0, FT(SF.gamma(z0)), α, D_min, D_max)
end
@inline function gamma_inc_moment_channel_setup(p0::Integer, α, D_min, D_max)
    z0 = p0 + 1
    FT = float(promote_type(typeof(α), typeof(D_min), typeof(D_max)))
    return _gamma_inc_moment_channel_setup(FT(z0), z0, FT(UT.fac(p0)), α, D_min, D_max)
end
@inline function _gamma_inc_moment_channel_setup(zf0, zpow0, Γ0, α, D_min, D_max)
    FT = float(promote_type(typeof(zf0), typeof(Γ0), typeof(α), typeof(D_min), typeof(D_max)))
    if !(α > 0)
        nan = FT(NaN)
        nan6 = (nan, nan, nan, nan, nan, nan)
        return (; zf0, α, invΓ = (nan, nan, nan, nan, nan), scale = nan6, Q_min = nan6, Q_max = nan6)
    end
    # Γ(z0+k) = (z0+k-1) Γ(z0+k-1), shared by the Q-chain's boundary-term
    # denominators (invΓ) and the per-order scale factors below.
    Γ1 = zf0 * Γ0
    Γ2 = (zf0 + 1) * Γ1
    Γ3 = (zf0 + 2) * Γ2
    Γ4 = (zf0 + 3) * Γ3
    Γ5 = (zf0 + 4) * Γ4
    invΓ = (inv(Γ1), inv(Γ2), inv(Γ3), inv(Γ4), inv(Γ5))

    invα = inv(α)
    scale0 = Γ0 / α^zpow0
    scale1 = scale0 * zf0 * invα
    scale2 = scale1 * (zf0 + 1) * invα
    scale3 = scale2 * (zf0 + 2) * invα
    scale4 = scale3 * (zf0 + 3) * invα
    scale5 = scale4 * (zf0 + 4) * invα
    scale = (scale0, scale1, scale2, scale3, scale4, scale5)

    Q_min = gamma_inc_Q_chain(zf0, α * D_min, invΓ)
    Q_max = gamma_inc_Q_chain(zf0, α * D_max, invΓ)
    return (; zf0, α, invΓ, scale, Q_min, Q_max)
end

"""
    gamma_inc_moment_channel_finish(setup, D_min, Dstar, D_max)

Combine a channel's outer-diameter-independent
[`gamma_inc_moment_channel_setup`](@ref) with a fresh `Q`-chain evaluation at
the crossing diameter `Dstar`, returning the six consecutive-order
crossing-split moment pairs
`(∫_{D_min}^{Dstar} D^p e^{-α D} dD, ∫_{Dstar}^{D_max} D^p e^{-α D} dD)`.
"""
@inline function gamma_inc_moment_channel_finish(setup, D_min, Dstar, D_max)
    (; zf0, α, invΓ, scale, Q_min, Q_max) = setup
    FT = typeof(scale[1])
    Q_star = gamma_inc_Q_chain(zf0, α * Dstar, invΓ)
    z0 = zero(FT)
    return ntuple(Val(6)) do k
        m_lo = D_min < Dstar ? scale[k] * max(Q_min[k] - Q_star[k], z0) : z0
        m_hi = Dstar < D_max ? scale[k] * max(Q_star[k] - Q_max[k], z0) : z0
        (m_lo, m_hi)
    end
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

Requires `N_ice` present (strictly positive in its value lane): `log(N_ice)`
is `-Inf` at zero and undefined below it. A caller that may see an absent
`N_ice` (including a zero-valued `ForwardDiff.Dual`) gates on presence before
calling, rather than this function inventing a value at absence.
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
iteration count). The budget is a property of the problem being solved, not of
this type, so each caller sets its own; the one in
[`get_distribution_logλ`](@ref) is measured, the others are not.
"""
struct FixedIterations{FT} <: RS.AbstractTolerance{FT} end
@inline (::FixedIterations)(x1, x2, y) = false

"""
    _mean_mass_target_logλ(state, a, b, target_mass)

Solve `log(target_mass) = log(a) + loggamma(b+μ+1) - loggamma(μ+1) - b·logλ` for `logλ`, where
`μ = get_μ(state, logλ)`, by a fixed-point iteration with a fixed step count.

`(a, b)` are the ice mass power-law coefficients (`a·D^b`) of whichever regime `target_mass`
falls in; the caller selects them (see [`_derived_logλ_bracket`](@ref)). The untruncated
gamma-moment ratio this inverts is the exact mean mass of a gamma-distributed population
confined to one power-law regime, and an approximation where the true generating distribution
spans more than one.

Three iterations, fixed rather than convergence-gated so the loop stays warp-convergent on GPU
(the same reasoning [`FixedIterations`](@ref) states for the outer bracketed solve): measured to
reach both precisions' rounding floor by the second iteration and bit-stable by the third, at
the derivation battery's own worst case.
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

Derive [`get_distribution_logλ`](@ref)'s default upper search bound `logλ_max` from the existing
mean-particle-mass floor ([`ice_mean_particle_mass_min`](@ref), already read by the ice number
adjustment); the lower bound `logλ_min` stays the literal `2`.

**The window's two edges are different kinds of quantity, not two instances of the same one, and
the asymmetry below is measured, not an oversight - do not re-symmetrize it on tidiness
grounds.** `ice_mean_particle_mass_min` is a PHYSICAL FLOOR: no ice population's mean mass sits
below it, so a bound derived from it is exact where it applies - [`ice_mass`](@ref) is a
piecewise power law, `a·D^b`, and at this smallest mean mass the relevant regime is `D → 0`, the
small-spherical-ice coefficients, which do not depend on the state's own `F_rim`/`ρ_rim`
(measured constant across every combination tried) and reproduce the production solve to machine
precision, every state tried. `ice_mean_particle_mass_max` is instead the ice number
adjustment's own REGULARIZATION TARGET - the value that process relaxes states TOWARD, not a
ceiling states respect. A production pin census (the gen-2 record, `219,086` populated cells)
found real graupel/hail cells with mean mass up to `0.0011` kg, over a hundred times that target,
and a bracket bound derived from it regressed pinning from `0.00` to `0.17` percent on the same
record. The literal `2` is retained rather than merely left alone: two independent censuses (the
ratio-form landing's own, and this one) show `0.00` percent pinning at this bound on real
records, and the largest observed production mean mass stays an order below the bound's own
`~13` cm implied mean size.

`logλ_max` is exact at this end (machine precision, every state in the derivation battery), so
the small margin added is a floating-point/iteration-precision cushion, not a correction for a
measured approximation error.
"""
@inline function _derived_logλ_bracket(state::P3State)
    FT = eltype(state)
    p3 = state.params
    mass_min = ice_mean_particle_mass_min(p3)
    (a_small, b_small) = ice_mass_coeffs(state, zero(FT))
    logλ_max = _mean_mass_target_logλ(state, a_small, b_small, mass_min)
    margin = FT(100) * eps(FT)  # a floating-point/iteration-precision cushion; this end's
    # derivation is exact, not an approximation - see the docstring above.
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
- `logλ_guess`: Optional initial guess used to narrow the search bracket (see
  [`_narrow_bracket`](@ref)) before the fixed-iteration solve runs; it does not change the
  iteration count itself. Within the 2M+P3 substep march, the caller passes the first
  substep's shape as this guess and refreshes the shape internally on every later substep of
  the same step. This argument is expected to narrow to a pure warm-start seed for the next
  timestep once the march owns its shape state internally; its contract is documented here so
  that transition does not require a signature change.
- `logλ_min`: The minimum value of the search bounds [log(1/m)]. Defaults to a value derived
  from `state` by [`_derived_logλ_bracket`](@ref) rather than a literal.
- `logλ_max`: The maximum value of the search bounds [log(1/m)]. Defaults likewise.
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
    # The target is the log of the MEAN PARTICLE MASS, so it is formed as ONE
    # log of the ratio and needs no floor on either operand. The mean mass is
    # well scaled across the whole physical range (1e-14 to 1e-3 kg) even when
    # both operands are individually tiny, which is exactly the regime where
    # splitting it into two logs forced a floor.
    #
    # The floors this replaces were `max(ρq_ice, eps(FT))` and
    # `max(ρn_ice, eps(FT))`, and the mass one was a defect rather than a
    # safeguard: `ρq_ice` runs 1e-12 to 1e-6 kg/m³ and so STRADDLES
    # `eps(Float32)` = 1.19e-7, while `ρn_ice` runs 1e3 to 1e6 m⁻³ and sits ten
    # orders above it. At the measured epicentre state the mass floor was
    # 39,480x larger than the mass and won, shifting the Float32 target by
    # 10.58 in the log; no logλ in [2, 17] could match it, the endpoints failed
    # to straddle, and the no-bracket fallback returned `logλ_min`. Worse than
    # the visible pin: below `eps(Float32)` the floored target became
    # INDEPENDENT OF THE MASS, identical for ρq_ice = 1e-12 and 1e-45 alike, so
    # the solve stopped seeing the ice mass at all and returned plausible-looking
    # wrong values. Measured across the physical grid, the ratio form makes
    # Float32 reproduce Float64 to five decimals where the floored form was off
    # by 4.8 to 11.7 in the target.
    #
    # A mass-free, number-carrying population (`ρq_ice ≤ 0`, `ρn_ice > 0`) returns
    # `logλ_max` directly, ahead of the no-bracket fallback below. `logLdivN` is strictly
    # decreasing in `logλ` (measured over the whole bracket), so the physical limit as mean
    # particle mass approaches zero is `logλ_max` (`hi`), the smallest particles the bracket
    # admits. Previously, `target_log_LdN = log(0 / ρn_ice) = -Inf` reached the no-bracket
    # branch below, where a magnitude tie-break returned `lo` instead - the largest
    # representable particle, about 13 cm in mean size, and a 77,734x error in
    # `ice_terminal_velocity_number_weighted` at a measured production census state. The
    # host's number-keyed presence mask does not exclude `ρq_ice = 0`, so this state is
    # reached in production (31.9% of ice-bearing cell-frames over the gen-3 record). The
    # literal `[2, 17]` bracket this landed against is now the DEFAULT VALUE's own derivation
    # target, not the value itself: see [`_derived_logλ_bracket`](@ref) for the per-state
    # bounds that replaced it; a prior measurement established that `[2, 17]` had margin to
    # spare before that replacement.
    if ρq_ice <= 0 && ρn_ice > 0
        return hi
    end

    # A pathological but populated pair (mean mass far below the nucleation mass, e.g. a
    # ratio of ~2e-45 kg) produces a target no `logλ` in `[lo, hi]` can match, at either
    # precision - a state with no valid shape rather than a precision failure. Making such
    # states impossible belongs to the ice minimum-quantities doctrine, not to this solver.
    #
    # Split into two logs here, not `log(ρq_ice / ρn_ice)`: forming the quotient first can
    # underflow to exactly zero in Float32 when `ρq_ice` is a positive subnormal divided by
    # a normal `ρn_ice`, giving `log(0) = -Inf` even though `ρq_ice` itself is nonzero. This
    # is NOT the floored split form the comment above retired - there is no `eps(FT)` floor
    # on either operand here, so the defect that motivated the ratio form (an oversized,
    # physically-unscaled floor overriding the real value) does not apply; `log` of any
    # positive float, subnormal included, is already finite. Matches the pre-fix target in
    # [`get_distribution_logλ_all_solutions`](@ref). The affected band sits 23 orders of
    # magnitude below the mass of a single 1 μm ice crystal, so no physical state reaches
    # it either way; this is a robustness fix, expected to move ordinary states by at most
    # a rounding-level amount from computing two logs instead of one, not a trajectory
    # change.
    target_log_LdN = log(ρq_ice) - log(ρn_ice)

    shape_problem(logλ) = logLdivN(state, logλ) - target_log_LdN
    f_lo, f_hi = shape_problem(lo), shape_problem(hi)
    if !isfinite(f_lo) || !isfinite(f_hi) || f_lo * f_hi > 0
        # A non-finite `target_log_LdN` has a known direction and is read off directly
        # rather than compared against another non-finite residual. Previously,
        # `abs(f_lo) ≤ abs(f_hi) ? lo : hi` compared `Inf ≤ Inf` (true regardless of which
        # end is correct) or, for `NaN`, `NaN ≤ NaN` (false, the right answer by accident).
        # Both are replaced by reading the target's sign:
        #   `target = -Inf` (mean mass → 0) → `hi`, by the same monotonicity as the early
        #     return above. The split-log form above no longer reaches this from a merely
        #     subnormal `ρq_ice` (that case is now a large finite target, handled by the
        #     ordinary bracketing path instead) - kept as a defensive branch rather than
        #     removed, since no exhaustive check has shown it unreachable.
        #   `target = +Inf` (mean mass → ∞, e.g. `ρn_ice = 0` with `ρq_ice > 0`) → `lo`, the
        #     opposite end by the same monotonicity.
        #   `target = NaN` (`ρq_ice = ρn_ice = 0`) → `hi`, matching the prior return but now
        #     explicit rather than an `NaN ≤ NaN` accident.
        #
        # A finite but unbracketable target (the populated-but-pathological case above) is
        # unaffected: `f_lo` and `f_hi` are both finite there, the `abs` comparison is a
        # distance measure between two representable residuals, and the returned bound is
        # the nearest representable one for a state with no shape.
        if !isfinite(target_log_LdN)
            isnan(target_log_LdN) && return hi
            return target_log_LdN > 0 ? lo : hi
        end
        return abs(f_lo) ≤ abs(f_hi) ? lo : hi
    end
    (lo, f_lo, hi, f_hi) =
        _narrow_bracket(shape_problem, lo, f_lo, hi, f_hi, logλ_guess)

    # Fixed iteration count (no early-exit) keeps GPU warps convergent. The
    # budget is the count at which each precision reaches its OWN rounding
    # floor, so it is a derived quantity rather than two tuned numbers. Measured
    # over 366 physical states spanning `F_rim` in [0, 0.99], `ρ_rim` in
    # [50, 900] kg/m³ and mean particle mass in [1e-14, 1e-3] kg, with the error
    # taken as the relative error in the diagnosed mean particle mass at the
    # returned root against a Float64-bisected reference:
    #
    #     iterations |   6      |   8      |  10      |  12      |  24
    #     Float32    | 1.1e0    | 7.4e-2   | 1.28e-3  | 1.28e-3  | 1.28e-3
    #     Float64    | 1.1e0    | 7.5e-2   | 3.9e-6   | 2.1e-14  | 1.4e-14
    #
    # The plateau is the rounding floor of `logLdivN` itself, so both slope laws
    # reach the same one (Float32 worst case 1.27994e-3 hard, 1.27896e-3
    # smoothed) and no budget can go below it. Before this measurement the Float32
    # budget was 8, where the worst error was 7.4e-2 with the hard
    # `SlopePowerLaw` and 2.3e-2 with `SmoothSlopePowerLaw`, worst at heavily
    # rimed small ice (`F_rim = 0.9`, mean mass 1e-9 kg), and mean particle mass
    # sets terminal velocity and every size-dependent rate. That budget
    # conflated representable precision with convergence rate: at 8 iterations
    # the error is ITERATION-limited at BOTH precisions, Float64 being no better
    # than Float32 there, so Brent needs the same count either way until the
    # floor is reached, and Float32, whose floor is the higher one, had been
    # given the FEWER iterations.
    #
    # The residual is asserted directly in `test/p3_tests.jl`. It is NOT
    # constrained by that file's `N ≈ ∫N′ dD` integral checks: `logN₀` is
    # derived from the returned `logλ`, so those integrals balance for any
    # root, converged or not.
    #
    # The per-cell GPU cost of the two extra iterations is UNMEASURED. Measure it
    # rather than assuming it negligible, since this solve runs on every cell of
    # every step.
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
