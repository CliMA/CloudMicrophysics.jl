"""
    Utilities

Lightweight numerical utility functions shared across CloudMicrophysics modules.
Contains pure numerical operations with no physics dependencies.
"""
module Utilities

import UnrolledUtilities as UU
import SpecialFunctions as SF
import ForwardDiff as FD

export clamp_to_nonneg, ϵ_numerics, ϵ_numerics_2M_M, ϵ_numerics_2M_N, ϵ_numerics_P3_B
export promote_typeof
export fac

"""
    promote_typeof(args...)

The common promoted type of the arguments' types.

Use it to type early returns and fallback values from all the arguments the
main-path result derives from. Typing them from a single argument
(`FT = eltype(q_tot)`-style) makes the function's return a union when a
caller mixes plain floats with `ForwardDiff.Dual`s (or float widths) across
arguments — non-concrete, heap-boxed, and silent.
"""
@inline promote_typeof(args...) = Base.promote_typeof(args...)
export unrolled_logsumexp
export sgs_weight_function, rime_mass_fraction, rime_density
export gamma_inc, gamma_inc_inv

"""
    gamma_inc(a, x)

Fast, fixed-iteration approximation of the incomplete gamma function
(replacing `SpecialFunctions.gamma_inc`).

Returns `(P(a, x), Q(a, x))` where `P` is the lower regularized incomplete gamma
and `Q` is the upper regularized incomplete gamma function.

Optimized for GPU execution and low compile times. The primary domain split
(`x < a + 1`) is an `if ... else` branch, which is cheap when neighbouring
grid points fall on the same side (spatial coherence in atmospheric fields),
and the per-iteration guards in the continued fraction use branchless `ifelse`.
The series / continued-fraction loops execute a fixed number of iterations
(20 for `Float32`, 30 for `Float64`), completely eliminating warp divergence
while guaranteeing excellent physical accuracy.

# Automatic differentiation
Analytic rules are provided for the derivative with respect to `x`
(`∂P/∂x = x^{a-1} e^{-x} / Γ(a)`, `∂Q/∂x = -∂P/∂x`), so AD does not trace the
internal iteration. The derivative with respect to the shape parameter `a` is
**not** implemented: passing an `a` that depends on the differentiation variable
(nonzero partials) raises an error rather than silently returning a wrong (zero)
gradient. See `_assert_const_shape`.
"""
@inline function gamma_inc(a::Real, x::Real)
    FT = float(promote_type(typeof(a), typeof(x)))
    return _gamma_inc(FT(a), FT(x))
end

# `P, Q` as `Dual`s carrying the analytic x-derivative. `val_a` is the (plain)
# shape parameter. Shared by the `gamma_inc(::Real, ::Dual)` and
# `gamma_inc(::Dual, ::Dual)` overloads (the latter after the shape guard).
@inline function _gamma_inc_dx(val_a::Real, x::FD.Dual)
    val_x = FD.value(x)
    P, Q = gamma_inc(val_a, val_x)
    T = FD.tagtype(typeof(x))
    deriv = val_x > 0 ? exp((val_a - 1) * log(val_x) - val_x - SF.loggamma(val_a)) : zero(val_x)
    return (FD.Dual{T}(P, deriv * FD.partials(x)), FD.Dual{T}(Q, -deriv * FD.partials(x)))
end

@inline gamma_inc(a::Real, x::FD.Dual) = _gamma_inc_dx(a, x)

@inline function gamma_inc(a::FD.Dual, x::FD.Dual)
    _assert_const_shape(a)
    return _gamma_inc_dx(FD.value(a), x)
end

@inline function gamma_inc(a::FD.Dual, x::Real)
    _assert_const_shape(a)
    P, Q = gamma_inc(FD.value(a), x)
    T = FD.tagtype(typeof(a))
    z = zero(FD.partials(a))
    return (FD.Dual{T}(P, z), FD.Dual{T}(Q, z))
end

# `lΓa = SF.loggamma(a)` is passed in so callers that evaluate this repeatedly for a
# fixed `a` (e.g. the Halley loop in `_gamma_inc_inv`) compute the expensive loggamma
# once instead of every iteration. The 2-arg form computes it for one-off callers.
@inline _gamma_inc(a::FT, x::FT) where {FT <: Real} = _gamma_inc(a, x, SF.loggamma(a))
@inline function _gamma_inc(a::FT, x::FT, lΓa::FT) where {FT <: Real}
    if x <= 0
        return (zero(FT), one(FT))
    elseif isinf(x)
        return (one(FT), zero(FT))
    end

    # factor = x^a * e^-x / Gamma(a)
    # Using loggamma for numerical stability
    factor = exp(a * log(x) - x - lΓa)
    maxiters = FT === Float32 ? 20 : 30

    if x < a + 1
        # Series expansion for P(a, x)
        term = one(FT) / a
        sum_P = term
        for k in 1:maxiters
            term *= x / (a + k)
            sum_P += term
        end
        P = factor * sum_P
        P = clamp(P, zero(FT), one(FT))
        return (P, one(FT) - P)
    else
        # Continued fraction (Lentz's method) for Q(a, x)
        tiny = FT(1e-30)

        # k = 1
        b_1 = x + 1 - a
        c = b_1 + 1 / tiny
        d = 1 / b_1
        h = d

        for k in 1:maxiters
            a_k = -FT(k) * (FT(k) - a)
            b_k = x + 2 * k + 1 - a

            d_tmp = b_k + a_k * d
            d = ifelse(abs(d_tmp) < tiny, tiny, d_tmp)

            c_tmp = b_k + a_k / c
            c = ifelse(abs(c_tmp) < tiny, tiny, c_tmp)

            d = 1 / d
            delta = c * d
            h *= delta
        end
        Q = factor * h
        Q = clamp(Q, zero(FT), one(FT))
        return (one(FT) - Q, Q)
    end
end

# Guard used by the AD rules for `gamma_inc` / `gamma_inc_inv`: the derivative
# with respect to the shape parameter `a` is not implemented (only the `x`- /
# `p`-derivative is). To avoid silently returning a wrong (zero) gradient, error
# if `a` actually depends on the differentiation variable. An `a` that merely
# happens to be a `Dual` with zero partials (e.g. promoted from a constant) is
# allowed and treated as a constant.
@inline function _assert_const_shape(a::FD.Dual)
    iszero(FD.partials(a)) || error(
        "gamma_inc/gamma_inc_inv: differentiation with respect to the shape " *
        "parameter `a` is not supported (only the x-/p-derivative is implemented).",
    )
    return nothing
end

"""
    gamma_inc_inv(a, p, q)

Fast, GPU-friendly inverse of `gamma_inc` using Halley's method.
Finds `x` such that `P(a, x) = p` and `Q(a, x) = q`.

# Automatic differentiation
The derivative with respect to `p` is provided analytically via the inverse
function theorem (`dx/dp = 1 / (∂P/∂x) = Γ(a) / (x^{a-1} e^{-x})`), so AD does
not trace Halley's iteration. As with [`gamma_inc`](@ref), differentiating with
respect to the shape parameter `a` is not supported and errors (see
`_assert_const_shape`).
"""
@inline function gamma_inc_inv(a::Real, p::Real, q::Real)
    FT = float(promote_type(typeof(a), typeof(p), typeof(q)))
    return _gamma_inc_inv(FT(a), FT(p), FT(q))
end

# `x` as a `Dual` carrying the analytic p-derivative. `val_a` is the (plain)
# shape parameter. Shared by the `gamma_inc_inv(::Real, ::Dual, ::Dual)` and
# `gamma_inc_inv(::Dual, ::Dual, ::Dual)` overloads (the latter after the guard).
@inline function _gamma_inc_inv_dp(val_a::Real, p::FD.Dual, q::FD.Dual)
    val_p = FD.value(p)
    val_q = FD.value(q)
    x_val = gamma_inc_inv(val_a, val_p, val_q)
    T = FD.tagtype(typeof(p))
    dP_dx = exp((val_a - 1) * log(x_val) - x_val - SF.loggamma(val_a))
    dx_dp = dP_dx > 0 ? inv(dP_dx) : zero(x_val)
    return FD.Dual{T}(x_val, dx_dp * FD.partials(p))
end

@inline gamma_inc_inv(a::Real, p::FD.Dual, q::FD.Dual) = _gamma_inc_inv_dp(a, p, q)

@inline function gamma_inc_inv(a::FD.Dual, p::FD.Dual, q::FD.Dual)
    _assert_const_shape(a)
    return _gamma_inc_inv_dp(FD.value(a), p, q)
end

@inline function gamma_inc_inv(a::FD.Dual, p::Real, q::Real)
    _assert_const_shape(a)
    x_val = gamma_inc_inv(FD.value(a), p, q)
    T = FD.tagtype(typeof(a))
    return FD.Dual{T}(x_val, zero(FD.partials(a)))
end

@inline function _gamma_inc_inv(a::FT, p::FT, q::FT) where {FT <: Real}
    if p <= 0
        return zero(FT)
    elseif q <= 0
        return FT(Inf)
    end

    # Initial guess
    if p < 0.5
        x = FT((p * SF.gamma(a + 1))^(1 / a))
    else
        x = FT(a - log(q))
    end

    # Halley's method
    # Use Q-q residual when p > 0.5 to avoid catastrophic cancellation
    use_q = p > FT(0.5)
    # loggamma(a) is loop-invariant, so compute once and thread it into `_gamma_inc`
    # and the derivative below (was recomputed ~2×/iteration, up to ~30×/call).
    lΓa = SF.loggamma(a)
    for i in 1:15
        P, Q = _gamma_inc(a, x, lΓa)
        f = use_q ? Q - q : P - p
        # Derivative: dP/dx = x^(a-1) e^-x / Gamma(a), dQ/dx = -dP/dx
        fprime = exp((a - 1) * log(x) - x - lΓa)
        # When using Q-q, the derivative is -fprime
        fprime = use_q ? -fprime : fprime
        if fprime == 0
            break
        end
        # f'' / f' = (a - 1 - x) / x  (same sign regardless of residual choice)
        fprime2_over_fprime = (a - 1 - x) / x

        # Halley step
        step = f / (fprime * (one(FT) - FT(0.5) * f * fprime2_over_fprime / fprime))

        # Protect against taking a step that makes x <= 0
        if x - step <= 0
            step = FT(0.5) * x
        end

        x = FT(x - step)
        if abs(step) < eps(FT) * x
            break
        end
    end
    return x
end

"""
    clamp_to_nonneg(x)

Clamp values to be non-negative.
Compatible with dual numbers (AD) and GPUs.

# Arguments
- `x`: value to clamp

# Returns
- `max(zero(x), x)`
"""
@inline clamp_to_nonneg(x) = max(zero(x), x)

"""
    fac(n)

Integer factorial `n!`, valid for `0 ≤ n ≤ 20`.
"""
@inline fac(n) = prod(1:n; init = one(n))

"""
    ϵ_numerics(FT)

Smallest number that is different than zero for the purpose of microphysics
computations. Returns `cbrt(floatmin(FT))` to avoid underflow issues.
"""
@inline ϵ_numerics(FT) = cbrt(floatmin(FT))

"""
    ϵ_numerics_2M_M(FT)

Numerical epsilon for 2-moment mass calculations.
"""
@inline ϵ_numerics_2M_M(FT) = eps(FT)

"""
    ϵ_numerics_2M_N(FT)

Numerical epsilon for 2-moment number calculations.
"""
@inline ϵ_numerics_2M_N(FT) = eps(FT)

"""
    ϵ_numerics_P3_B(FT)

Numerical epsilon for P3 bulk microphysics mass calculations
    relating to rim volume, B_rim.
"""
@inline ϵ_numerics_P3_B(FT) = eps(FT)


"""
    unrolled_logsumexp(x)

Compute `log(sum(exp, x))` in a statically unrolled fashion.

This method uses [`UnrolledUtilities`](https://github.com/CliMA/UnrolledUtilities.jl)
to produce fully unrolled code with no dynamic dispatch or reductions,
making it transparent to GPU compilers.

The standard shift-by-max trick is used for numerical stability.

Note: This code is two-pass (find max, then sum shifted exponentials). 
LogExpFunctions.jl implements a one-pass version, but is not unrolled,
so may result in more complicated GPU code.

## Extended help

Other implementation options were considered, detailed below, and may be revisited in the future.
For now, this implementation is sufficient.

### Naive implementation

This is the most straightforward implementation, but it is not numerically stable.

```julia
log(UU.unrolled_sum(exp, x))
```

### One-pass unrolled implementation

This is reaches into LogExpFunctions.jl internals,

```julia
FT = eltype(x)
return LogExpFunctions._logsumexp_onepass_result(
    UU.unrolled_reduce(LogExpFunctions._logsumexp_onepass_op, x, (FT(-Inf), zero(FT)))
)
```

### Dispatch-wrapper for reduce

Pass a wrapper to compile to unrolled reduce
```julia
# Note: This is a sketch, not tested
struct UnrolledWrapper{T}
    x::T
end
Base.iterate(w::UnrolledWrapper) = iterate(w.x)
Base.iterate(w::UnrolledWrapper, state) = iterate(w.x, state)
Base.length(w::UnrolledWrapper) = length(w.x)
Base.eltype(w::UnrolledWrapper) = eltype(w.x)
Base.reduce(op, w::UnrolledWrapper) = UU.unrolled_reduce(op, w.x)  # use unrolled reduce
# ... then call:
LogExpFunctions.logsumexp(UnrolledWrapper(x))
```
"""
function unrolled_logsumexp(x)
    # Find the maximum (ps: if any element is NaN, then xmax = NaN)
    xmax = UU.unrolled_maximum(x)

    # Handle non-finite values: if xmax is +Inf or -Inf or NaN, return it directly
    # (avoids Inf - Inf = NaN and x - NaN = NaN in the shifted exponentials below)
    isfinite(xmax) || return xmax

    # Sum shifted exponentials
    shifted_exp(xi) = exp(xi - xmax)
    s = UU.unrolled_sum(shifted_exp, x)

    return xmax + log(s)
end


"""
    sgs_weight_function(a, a_half)

Smooth, monotonic weight function `w(a)` ranging from 0 to 1.

Used as the interpolation weight in regularised divisions so the result stays
well-defined and smoothly blends toward zero when the denominator is small.

Key properties:
- `w(a) = 0` for `a ≤ 0`.
- `w(a) = 1` for `a ≥ 1`.
- `w(a_half) = 0.5`.
- Continuously differentiable; derivatives vanish at `a = 0` and `a = 1` so
  the blend is smooth.
- Grows very rapidly near `a_half`, very slowly elsewhere.

Construction: for `a ∈ (0, 1)`, a bounded sigmoid is built by composing `tanh`
with the inverse of a slower `tanh`; midpoint control `(a_half, 0.5)` is
enforced by pre-transforming `a` via `1 - (1 - a)^k`.

Mirrors the `sgs_weight_function` in `ClimaAtmos.jl/src/utils/variable_manipulations.jl`.

# Arguments
- `a`: the input variable (often approximated as an area fraction `ρa / ρ`).
- `a_half`: value of `a` at which the weight equals 0.5 (controls the
  transition point of the sigmoid).

# Returns
- `w(a)` in `[0, 1]`.
"""
@inline function sgs_weight_function(a, a_half)
    if a < 0
        zero(a)
    elseif a > min(1, 42 * a_half)   # autodiff generates NaNs when a is large
        one(a)
    elseif 4 * a < eps(typeof(a))
        # 1 - a rounds to 1, making atanh(-1) = -Inf: the value is 0 either
        # way, but autodiff generates NaNs (mirrors the upper guard)
        zero(a)
    else
        (1 + tanh(2 * atanh(1 - 2 * (1 - a)^(-1 / log2(1 - a_half))))) / 2
    end
end

"""
    _regularised_ratio(numerator, denominator, half, ϵ)

Compute `numerator / denominator` with `sgs_weight_function`-based
regularisation so the result stays finite when `denominator` is zero
or very small.

Returns `weight(denominator) * numerator / denominator`, falling to
zero when `denominator` is below machine precision.
"""
@inline function _regularised_ratio(
    numerator, denominator,
    half = eps(typeof(denominator)),
    ϵ = eps(typeof(denominator))^2,
)
    weight = sgs_weight_function(denominator, half)
    # zero of the promoted type: a single-argument zero makes the return a
    # union when numerator and denominator mix plain floats with Duals
    z = zero(promote_typeof(numerator, denominator))
    return ifelse(denominator < ϵ, z, weight * numerator / denominator)
end

"""
    rime_mass_fraction(q_rim, q_ice, q_ice_half)

Regularised rime mass fraction `F_rim = q_rim / q_ice` that stays finite when
`q_ice` is zero or very small, with the result clamped to `[0, 1]` via
`min(q_rim, q_ice)`.

# Arguments
- `q_rim`: rime specific mass `[kg rim / kg air]`.
- `q_ice`: total ice specific mass `[kg ice / kg air]`.
- `q_ice_half`: value of `q_ice` at which the blending weight equals 0.5
  (default `eps(typeof(q_ice))`).
"""
@inline rime_mass_fraction(q_rim, q_ice, kw...) =
    _regularised_ratio(min(q_rim, q_ice), q_ice, kw...)

"""
    rime_density(q_rim, b_rim, b_rim_half)

Regularised rime density `ρ_rim = q_rim / b_rim` that stays finite when
`b_rim` is zero or very small.

# Arguments
- `q_rim`: rime specific mass `[kg rim / kg air]`.
- `b_rim`: rime specific volume `[m³ rim / kg air]`.
- `b_rim_half`: value of `b_rim` at which the blending weight equals 0.5
  (default `eps(typeof(b_rim))`).
"""
@inline rime_density(q_rim, b_rim, kw...) = _regularised_ratio(q_rim, b_rim, kw...)


end # module
