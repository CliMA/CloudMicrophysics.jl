import CloudMicrophysics.DistributionTools: size_distribution

# Callable returned by `logN′ice`: evaluates `log(N′(D))` for a fixed distribution.
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
    logN′ice(dist)

Return a callable that computes `log(N′(D))` for the [`P3Distribution`](@ref) `dist`.

`N′(D)` is the ice particle number concentration at diameter `D`.

# Arguments
- `dist`: [`P3Distribution`](@ref)
"""
function logN′ice((; μ, logλ, logN₀))
    # Promote to a common type: differentiating w.r.t. the ice number makes
    # `logN₀` a `Dual` while `μ` and `logλ` stay plain floats, and
    # `P3LogNumberFunctor` stores all three in a single field type.
    return P3LogNumberFunctor(promote(logN₀, μ, exp(logλ))...)
end

# Callable returned by `size_distribution`: `n(D) = exp(logN′(D))`.
struct P3SizeDistributionFunctor{F} <: Function
    logN′::F
end
@inline (f::P3SizeDistributionFunctor)(D) = exp(f.logN′(D))

"""
    size_distribution(dist::P3Distribution)

Return `n(D)`, the ice size distribution of the [`P3Distribution`](@ref) `dist` at diameter `D`.
"""
DT.size_distribution(dist::P3Distribution) = P3SizeDistributionFunctor(logN′ice(dist))

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
"""
function loggamma_inc_moment(D₁, D₂, μ, logλ, k = 0, scale = 1)
    FT = UT.promote_typeof(D₁, D₂, μ, logλ)
    D₁ < D₂ || return FT(-Inf)  # an empty segment
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
    loggamma_moment(μ, logλ; [k = 0], [scale = 1])

Compute `log(scale ⋅ ∫_0^∞ G(D) D^k dD)`.

`G(D) ≡ D^μ e^{-λD}` is the (unnormalized) gamma kernel, `k` is an arbitrary exponent, and `scale`
is a scale factor.

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
    get_μ(slope, logλ)

Compute the shape parameter μ.

# Arguments
- `slope`: [`CMP.SlopeLaw`](@ref)
- `logλ`: The log of the slope parameter [log(1/m)]
"""
get_μ((; a, b, c, μ_max)::CMP.SlopePowerLaw, logλ) = clamp(a * exp(logλ)^b - c, 0, μ_max)
function get_μ((; a, b, c, μ_max, κ)::CMP.SmoothSlopePowerLaw, logλ)
    p = a * exp(logλ)^b - c
    M₀ = LogExpFunctions.log1pexp(κ * p) / κ
    return μ_max - LogExpFunctions.log1pexp(κ * (μ_max - M₀)) / κ
end
get_μ((; μ)::CMP.SlopeConstant, logλ...) = μ

"""
    logmass_gamma_moment(particles, μ, logλ; n = 0)

Compute `log(∫_0^∞ Dⁿ m(D) N′(D) dD)` given the `particles`, `μ`, and `logλ`.

This is the log of the `n`-th moment of the mass-weighted PSD.

# Arguments
- `particles`: [`P3Particles`](@ref)
- `μ`: The shape parameter [`-`]
- `logλ`: The log of the slope parameter [log(1/m)]

# Keyword arguments
- `n`: The order of the moment, default is `0`

# Note:
- For `n = 0`, this evaluates to `log(L/N₀)`
- For `n = 1`, this evaluates to the (unnormalized) mass-weighted mean particle size, see [`D_m`](@ref)
"""
function logmass_gamma_moment(particles, μ, logλ; n = 0)
    bnds = segment_boundaries(particles)
    moments = UU.unrolled_map(subintervals(bnds)) do (D_lo, D_hi)
        (a, b) = ice_mass_coeffs(particles, (D_lo + D_hi) / 2)
        loggamma_inc_moment(D_lo, D_hi, μ, logλ, b + n, a)
    end
    return UT.unrolled_logsumexp(moments)
end

"""
    logLdivN(particles, logλ)

Compute `log(L/N)` given the `particles` and `logλ`.

# Arguments
- `particles`: [`P3Particles`](@ref)
- `logλ`: The log of the slope parameter [log(1/m)]
"""
function logLdivN(particles, logλ)
    μ = get_μ(particles.params.slope, logλ)
    logLdivN₀ = logmass_gamma_moment(particles, μ, logλ; n = 0)
    logNdivN₀ = loggamma_moment(μ, logλ; k = 0)
    return logLdivN₀ - logNdivN₀
end

"""
    get_logN₀(N_ice, μ, logλ)

Compute `log(N₀)` given `N_ice`, `μ`, and `logλ`.

The intercept follows from

        N  = N₀ ∫ G(D) dD
    log N₀ = log N - log(∫G(D) dD)
           = log(N) - log( ∫D^μ e^{-λD} dD )
           = log(N) - log(M⁰)

# Arguments
- `N_ice`: The number concentration [1/m³]
- `μ`: The shape parameter [`-`]
- `logλ`: The log of the slope parameter [log(1/m)]

Without ice number, `log(N₀) = -Inf`, with zero derivatives.
"""
function get_logN₀(N_ice, μ, logλ)
    logN₀ = log(N_ice) - loggamma_moment(μ, logλ; k = 0)
    return ifelse(N_ice ≥ floatmin(typeof(N_ice)), logN₀, oftype(logN₀, -Inf))
end

"""
    FixedIterations{FT}()

A `RootSolvers.AbstractTolerance` whose convergence predicate is always `false`.

The solver stops at the iteration budget instead of at a tolerance. This bounds
the cost of a GPU warp by the budget for every input. `BrentsMethod` still returns
early when the residual is exactly zero, which leaves the result unchanged. A tighter
initial bracket improves accuracy but not the budget. Each caller sets its own budget.
"""
struct FixedIterations{FT} <: RS.AbstractTolerance{FT} end
@inline (::FixedIterations)(x1, x2, y) = false

"""
    _logλ_max(params)

Return the largest `logλ` of the shape solve.

It is the slope at which the mean particle mass equals the mass of a newly nucleated crystal, a
solid-ice sphere of diameter `D_nuc`.

Small ice is spherical, `m(D) = (π/6) ρ_i D³`, so the mean particle mass is
`(π/6) ρ_i Γ(μ + 4) / Γ(μ + 1) / λ³`. Equating it with `(π/6) ρ_i D_nuc³` gives

```math
λ_{max} D_{nuc} = \\sqrt[3]{(μ + 1)(μ + 2)(μ + 3)}, \\qquad μ = μ(λ_{max}).
```

One update from `μ = 0` solves this exactly while `μ` is constant at such small sizes, as it is
for the default slope law, [`CMP.SmoothSlopePowerLaw`](@ref).

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
"""
@inline function _logλ_max((; slope, D_nuc))
    logλ(μ) = log(cbrt((μ + 1) * (μ + 2) * (μ + 3)) / D_nuc)
    return logλ(get_μ(slope, logλ(zero(D_nuc))))
end

"""
    logλ_bounds(params)

Return the smallest and largest `logλ` [log(1/m)] of the slope solve.

The smallest is 2, and the largest is the slope at which the mean particle mass equals the mass of a
newly nucleated crystal.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
"""
@inline logλ_bounds(params) = (2 * one(params.D_nuc), _logλ_max(params))

"""
    get_distribution_logλ(particles, ρq_ice, ρn_ice)

Solve for `logλ`, the log of the slope parameter of the P3 size distribution.

The solve uses the ice mass `ρq_ice`, the bounded ice number `ρn_ice`, and the mass relation of the
[`P3Particles`](@ref) `particles`.

The assumed distribution is of the form

```math
N′(D) = N₀ D^μ e^{-λD}
```
where `N′(D)` is the number concentration at diameter `D` and `μ` is the shape parameter.
    The shape parameter is parameterized, e.g. [`CMP.SlopePowerLaw`](@ref) or [`CMP.SlopeConstant`](@ref).

This algorithm solves for `logλ = log(λ)`
    given `ρq_ice` and `ρn_ice` by solving the equations:

```math
\\begin{align*}
\\log(L) &= \\log ∫_0^∞ m(D) N′(D)\\ \\mathrm{d}D, \\\\
\\log(N) &= \\log ∫_0^∞ N′(D)\\ \\mathrm{d}D, \\\\
\\end{align*}
```
where `m(D)` is the mass of a particle at diameter `D` (see [`ice_mass`](@ref)).
    The procedure is described in detail in [the P3 docs](@ref "Parameterizations for the slope parameter \$μ\$").

The solve uses the primal (as opposed to `ForwardDiff.Dual`) values of its arguments and returns
a primal value, so `logλ` is constant under differentiation with respect to the state.

# Arguments
- `particles`: The [`P3Particles`](@ref)
- `ρq_ice`: The ice mass concentration [kg/m³]
- `ρn_ice`: The ice number concentration [1/m³], bounded by [`bounded_ice_number`](@ref)

The search bounds are [`logλ_bounds`](@ref). Without ice mass, the result is the upper bound, the
smallest particles. With mass but without number, it is the lower bound, the largest particles.
"""
function get_distribution_logλ(particles, ρq_ice, ρn_ice)
    particles = _primal_particles(particles)
    FT = eltype(particles)
    ρq_ice = FT(FD.value(ρq_ice))
    ρn_ice = FT(FD.value(ρn_ice))
    lo, hi = map(FT, logλ_bounds(particles.params))
    ρq_ice ≥ floatmin(FT) || return hi
    ρn_ice ≥ floatmin(FT) || return lo

    # Two logs instead of the log of the ratio, which can underflow to zero in Float32
    target_log_LdN = log(ρq_ice) - log(ρn_ice)

    shape_problem(logλ) = logLdivN(particles, logλ) - target_log_LdN
    f_lo, f_hi = shape_problem(lo), shape_problem(hi)
    if !isfinite(f_lo) || !isfinite(f_hi) || f_lo * f_hi > 0
        return abs(f_lo) ≤ abs(f_hi) ? lo : hi
    end

    # A fixed budget instead of a tolerance bounds the GPU cost. 18 iterations reach the rounding
    # floor of `logLdivN` in both precisions, including the flat stretch and the knee of the
    # `μ` ramp; `test/p3_tests.jl` asserts the residual.
    maxiters = 18
    sol = RS.find_zero(
        shape_problem,
        RS.BrentsMethod(lo, hi),
        RS.CompactSolution(),
        FixedIterations{FT}(),
        maxiters,
    )
    return sol.root  # logλ
end
