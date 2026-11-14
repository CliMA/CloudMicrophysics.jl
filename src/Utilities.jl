"""
    Utilities

Lightweight numerical utility functions shared across CloudMicrophysics
modules, plus regularised helpers for derived quantities such as the rime
mass fraction and rime density.
"""
module Utilities

import UnrolledUtilities as UU
import SpecialFunctions as SF
import ForwardDiff as FD
import RootSolvers as RS
import StaticArrays as SA

export clamp_to_nonneg, ϵ_numerics
export promote_typeof, return_type_of, find_zeros, scan_zeros
export fac

"""
    promote_typeof(args...)

The common promoted type of the arguments' types.

Similar to `promote_type`, but accepts concrete variables rather than
their types. Useful to guarantee a single, concrete return type for functions
that accept mixed input types (e.g. `Float32` and `ForwardDiff.Dual{Float32}`)
"""
@inline promote_typeof(args...) = Base.promote_typeof(args...)

"""
    fast_div(a, b)

`a / b`, which for `Float32` carries the LLVM fast-math flag `afn` ("approximate functions").

On a CPU the flag changes nothing: the division stays the correctly rounded IEEE division. On
a GPU it selects the approximate division of the hardware (PTX `div.approx.f32`, error up to
2 ulp for a divisor of magnitude between `2^-126` and `2^126`; a finite divisor above `2^126`,
about `8.5e37`, gives `0`, or `NaN` for an infinite numerator) in place of the correctly rounded one, which CUDA evaluates in a
subroutine that costs several times a multiplication and up to ten times one when an operand
or the result is zero, subnormal or infinite. None of the value-changing fast-math assumptions
(no `NaN`, no `Inf`, reassociation) is made, so infinite time scales and `ifelse` guards keep
their meaning. Every other number type (`Float64`, dual numbers, arrays) uses `/`, except
that a static array is divided element by element.
"""
@inline fast_div(a::Float32, b::Float32) =
    Base.llvmcall("%3 = fdiv afn float %0, %1\nret float %3", Float32, Tuple{Float32, Float32}, a, b)
@inline fast_div(a::Real, b::Real) = _fast_div(promote(a, b)...)
@inline fast_div(a::SA.StaticArray, b::Real) = map(x -> fast_div(x, b), a)
@inline fast_div(a, b) = a / b
@inline _fast_div(a::Float32, b::Float32) = fast_div(a, b)
@inline _fast_div(a, b) = a / b

"""
    fast_pow(x, y)

`x^y` for a base `x ≥ 0`, evaluated for two `Float32` arguments with `Base.FastMath.pow_fast`.

On a GPU that is `exp2(y * log2(x))` in single precision, in place of the `Float32` `^` of
CUDA, which is evaluated in double precision. The relative error is a few
`eps(Float32) * |y log2(x)|`. Every other combination of types uses `^`.
"""
@inline fast_pow(x::Float32, y::Float32) = Base.FastMath.pow_fast(x, y)
@inline fast_pow(x, y) = x^y

"""
    fast_sqrt(x)

`sqrt(x)` for `x ≥ 0`. For `Float32` it is the LLVM square root with the fast-math flag `afn`,
without the domain check of `Base.sqrt`: a negative argument gives `NaN` instead of an error.

On a CPU the result is the correctly rounded square root. On a GPU it is the approximate one
of the hardware (PTX `sqrt.approx.f32`), and the kernel no longer carries the error path of
the domain check. Every other number type uses `sqrt`.
"""
@inline fast_sqrt(x::Float32) = Base.llvmcall(
    (
        """
        declare float @llvm.sqrt.f32(float)
        define float @entry(float %0) alwaysinline {
            %2 = call afn float @llvm.sqrt.f32(float %0)
            ret float %2
        }
        """,
        "entry",
    ), Float32, Tuple{Float32}, x)
@inline fast_sqrt(x) = sqrt(x)

"""
    @approx_ops expr

Replace, in `expr`, every `a / b` by [`fast_div`](@ref)`(a, b)`, every `sqrt(x)` by
[`fast_sqrt`](@ref)`(x)` and every `x^y` whose exponent is not an integer literal by
[`fast_pow`](@ref)`(x, y)`.

The macro is applied to the definitions of the functions evaluated at every grid point on a
GPU. It changes `Float32` evaluations only: on a GPU to operations accurate to a few units in
the last place, on a CPU only through `fast_pow` and the missing domain check of `fast_sqrt`.
Unlike `@fastmath` it assumes nothing about `NaN` or `Inf` and does not reorder operations.
"""
macro approx_ops(expr)
    return esc(_approx_ops(expr))
end
_approx_ops(x) = x
function _approx_ops(ex::Expr)
    args = map(_approx_ops, ex.args)
    if ex.head === :call && length(args) == 3 && args[1] === :/
        return Expr(:call, fast_div, args[2], args[3])
    elseif ex.head === :call && length(args) == 3 && args[1] === :^ && !(args[3] isa Integer)
        return Expr(:call, fast_pow, args[2], args[3])
    elseif ex.head === :call && length(args) == 2 && args[1] === :sqrt
        return Expr(:call, fast_sqrt, args[2])
    end
    return Expr(ex.head, args...)
end

"""
    return_type_of(f, args...)

Return the type that `f(args...)` returns, as inferred by the compiler, without calling `f`.

The type is concrete for a type-stable `f`, which every integrand must be.
"""
@inline return_type_of(f::F, args...) where {F} = Core.Compiler.return_type(f, typeof(args))

"""
    find_zeros(f, method, x_lo, x_hi, ::Val{n_scan}, ::Val{max_roots}, tol, maxiters)

Find up to `max_roots` zeros of `f` in `[x_lo, x_hi]`, in ascending order.

`f` is evaluated at `n_scan + 1` evenly spaced points, and the first `max_roots` sign changes between
adjacent points are recorded. After the scan, each recorded sign change is refined by
`RootSolvers.find_zero` with the two-point `method`, which starts from the residuals of the scan. A
scan point at which `f` is zero is a zero. Two zeros between adjacent points are missed.

# Arguments
- `f`: function of one argument
- `method`: two-point method type of RootSolvers, e.g. `RegulaFalsiMethod`
- `x_lo`, `x_hi`: bounds of the search
- `n_scan`: number of scan intervals
- `max_roots`: largest number of zeros returned
- `tol`, `maxiters`: tolerance and iteration limit of `RootSolvers.find_zero`

# Returns
- `(; roots, n_found, f_lo)`: an `NTuple{max_roots}` of the zeros, with `x_lo` in place of each zero not
  found, the number of zeros found, and `f(x_lo)`
"""
function find_zeros(
    f::F, ::Type{M}, x_lo, x_hi, n_scan::Val, max_roots::Val, tol, maxiters,
) where {F, M}
    no_summand = x -> (f(x), nothing)
    (; roots, n_found, f_lo) =
        scan_zeros(no_summand, f, M, x_lo, x_hi, n_scan, max_roots, tol, maxiters, (i, Δx) -> nothing)
    return (; roots, n_found, f_lo)
end

"""
    scan_zeros(g, f, method, x_lo, x_hi, ::Val{n_scan}, ::Val{max_roots}, tol, maxiters, weight)

Find up to `max_roots` zeros of `f` in `[x_lo, x_hi]` as [`find_zeros`](@ref) does, and compute a weighted sum of a second function `h` over the scan points.

At each of the `n_scan + 1` scan points `x`, `g(x)` returns `(f(x), h(x))`, so that one evaluation
gives the residual of the scan and the summand. The first element of `g(x)` must equal `f(x)`, and
`f` alone refines the sign changes. The sum is `Σᵢ weight(i, Δx) h(x_lo + i Δx)` over the scan
points `i = 0, …, n_scan` of spacing `Δx`, and it is `nothing` if `h` returns `nothing`.

# Returns
- `(; roots, n_found, f_lo, total)`: as [`find_zeros`](@ref), and the weighted sum `total` of `h`
"""
@inline function scan_zeros(
    g::G, f::F, ::Type{M}, x_lo, x_hi, ::Val{n_scan}, ::Val{max_roots}, tol, maxiters, weight::W,
) where {G, F, M, n_scan, max_roots, W}
    Δx = (x_hi - x_lo) / n_scan
    (f_lo, summand) = g(x_lo)
    total = _weighted(weight(0, Δx), summand)
    # Each bracket is `(x_a, x_b, f(x_a), f(x_b))`, with `x_a = x_b = x_lo` for a zero not found
    brackets = ntuple(_ -> (x_lo, x_lo, f_lo, f_lo), Val(max_roots))
    exact = ntuple(_ -> false, Val(max_roots))
    n_found = 0
    x_prev, y_prev = x_lo, f_lo
    for i in 1:n_scan
        x = x_lo + i * Δx
        (y, summand) = g(x)
        total = _add(total, _weighted(weight(i, Δx), summand))
        record = (iszero(y) | (y * y_prev < 0)) & (n_found < max_roots)
        # No entry has the index 0, so an unrecorded point changes no bracket
        k = ifelse(record, n_found + 1, 0)
        brackets = _replace_entry(brackets, (x_prev, x, y_prev, y), k)
        exact = _replace_entry(exact, iszero(y), k)
        n_found += record
        x_prev, y_prev = x, y
    end
    found = ntuple(identity, Val(max_roots)) .≤ n_found
    roots = map(
        @inline((bracket, at_point, is_found) -> _refined_zero(f, M, bracket, at_point, is_found, x_lo, tol, maxiters)),
        brackets,
        exact,
        found,
    )
    return (; roots, n_found, f_lo, total)
end

@inline function _refined_zero(
    f::F, ::Type{M}, (x_a, x_b, y_a, y_b), at_point, is_found, x_lo, tol, maxiters,
) where {F, M}
    is_found || return x_lo
    at_point && return x_b
    return @inline(RS.find_zero(f, M, x_a, x_b, y_a, y_b, RS.CompactSolution(), tol, maxiters)).root
end

@inline _weighted(w, summand) = w * summand
@inline _weighted(w, ::Nothing) = nothing
@inline _add(a, b) = a + b
@inline _add(::Nothing, ::Nothing) = nothing

# The tuple `t` with entry `i` replaced by `v`, without the bounds check of `Base.setindex`
@inline _replace_entry(t::NTuple{N}, v, i) where {N} = ntuple(k -> ifelse(k == i, v, t[k]), Val(N))
export unrolled_logsumexp
export rime_mass_fraction, rime_density
export guarded_quotient, nearest_admissible_b
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
    # The derivative vanishes as `x → ∞`: set the correct limit.
    deriv =
        (val_x > 0) && isfinite(val_x) ?  # ternary to avoid log(x) when x < 0 (which throws)
        exp((val_a - 1) * log(val_x) - val_x - SF.loggamma(val_a)) : zero(val_x)
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

# `loggamma_a = SF.loggamma(a)` is passed in so callers that evaluate this repeatedly for a
# fixed `a` (e.g. the Halley loop in `_gamma_inc_inv`) compute the expensive loggamma
# once instead of every iteration. The 2-arg form computes it for one-off callers.
@inline _gamma_inc(a::FT, x::FT) where {FT <: Real} = _gamma_inc(a, x, SF.loggamma(a))
@inline function _gamma_inc(a::FT, x::FT, loggamma_a::FT) where {FT <: Real}
    if x <= 0
        return (zero(FT), one(FT))
    elseif isinf(x)
        return (one(FT), zero(FT))
    end

    # factor = x^a * e^-x / Gamma(a)
    # Using loggamma for numerical stability
    factor = exp(a * log(x) - x - loggamma_a)
    # For Float32, 14 iters gives 0.029% error at the 99th percentile, 0.47% max error
    maxiters = FT === Float32 ? 14 : 30

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
        "gamma_inc/gamma_inc_inv: differentiation with respect to the shape \
        parameter `a` is not supported (only the x-/p-derivative is implemented).",
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
    if p < FT(0.5)
        x = FT((p * SF.gamma(a + 1))^(1 / a))
    else
        x = FT(a - log(q))
    end

    # Halley's method
    # Use Q-q residual when p > 0.5 to avoid catastrophic cancellation
    use_q = p > FT(0.5)
    # loggamma(a) is loop-invariant, so compute once and thread it into
    # `_gamma_inc` and the derivative below rather than recomputing per call.
    loggamma_a = SF.loggamma(a)
    for i in 1:15
        P, Q = _gamma_inc(a, x, loggamma_a)
        f = use_q ? Q - q : P - p
        # Derivative: dP/dx = x^(a-1) e^-x / Gamma(a), dQ/dx = -dP/dx
        fprime = exp((a - 1) * log(x) - x - loggamma_a)
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

# ---------------------------------------------------------------------------
# Small-number handling: two distinct roles, kept deliberately separate.
#
# 1. Domain sanitization (NOT a physical threshold): `clamp_to_nonneg` and the
#    `max(x, ϵ)` floors on arguments of `cbrt`/`log`/`^`. These only keep an
#    expression finite (e.g. so the *discarded* branch of a branchless
#    `ifelse(cond, rate, 0)` cannot produce Inf/NaN, which matters for
#    ForwardDiff gradients). They do not decide whether a tracer is "present".
#
# 2. Physical thresholds: `ϵ_numerics` below, used as `q > ϵ` gates that treat
#    a tracer as absent. There is exactly one threshold, always called from here
#    so the value has a single source of truth (no bare literals/zeros scattered
#    across the schemes).
#
# The mass and number thresholds of the 1-moment and 2-moment schemes are the same
# value, `cbrt(floatmin(FT))` (≈2.3e-13 @ f32, ≈2.8e-103 @ f64). Its cube is
# `floatmin`, so products of up to three small contents in the rate expressions
# stay normal floats. It is many orders of magnitude below any physically relevant
# content in both precisions, so the rates do not depend on the precision.
# ---------------------------------------------------------------------------

"""
    clamp_to_nonneg(x)

Clamp values to be non-negative. Compatible with dual numbers (AD) and GPUs.

This is **domain sanitization**, not a physical threshold: use it to keep an
argument of `cbrt`/`log`/`^` finite, not to decide whether a tracer is present
(use the [`ϵ_numerics`](@ref) gate for that).

# Arguments
- `x`: value to clamp

# Returns
- `max(zero(x), x)`
"""
@inline clamp_to_nonneg(x) = max(zero(x), x)

"""
    fac(n)

Integer factorial `n!`. Throws a `DomainError` unless `0 ≤ n ≤ 20`, the range
representable in `Int64`.
"""
@inline function fac(n)
    # 21! overflows Int64, so restrict to 0 ≤ n ≤ 20
    0 ≤ n ≤ 20 || throw(DomainError(n, "fac(n) is defined for 0 ≤ n ≤ 20"))
    return prod(1:n; init = one(n))
end

"""
    ϵ_numerics(FT)
    ϵ_numerics(x::Real)

Physical smallness threshold for the **1-moment** and **2-moment** schemes.

Below this threshold a tracer is treated as absent.

Returns `cbrt(floatmin(float(FT)))`, a floatmin-derived floor
that also keeps `cbrt`/`^` arguments above `floatmin`.

Can be called in three ways:
- With a floating-point type (e.g., `Float32` or `Float64`), returning the threshold for that type.
- With a concrete floating-point value, returning the threshold for its type.
- With an integer type or value, returning the threshold for its float counterpart (e.g., `Int32` → `Float32`).
"""
@inline ϵ_numerics(FT) = cbrt(floatmin(float(FT)))
# The two common float types return a value computed once, on the CPU. A GPU compiler does not
# fold `cbrt` of a constant, so the generic method would evaluate it in every kernel and keep
# the result, and every `max(x, ϵ)` of a constant `x`, in a register.
const _ϵ_numerics_f32 = cbrt(floatmin(Float32))
const _ϵ_numerics_f64 = cbrt(floatmin(Float64))
@inline ϵ_numerics(::Type{Float32}) = _ϵ_numerics_f32
@inline ϵ_numerics(::Type{Float64}) = _ϵ_numerics_f64
@inline ϵ_numerics(x::Real) = ϵ_numerics(typeof(x))

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
    _regularised_ratio(numerator, denominator, ϵ)

Compute `numerator / denominator`, returning zero when `denominator` is below `ϵ`.
"""
@inline function _regularised_ratio(numerator, denominator, ϵ)
    denominator_safe = ifelse(denominator > ϵ, denominator, one(denominator))
    val = numerator / denominator_safe
    return ifelse(denominator > ϵ, val, zero(val))
end

"""
    rime_mass_fraction(q_rim, q_ice, ϵ = 1e-24)

Regularised rime mass fraction `F_rim = q_rim / q_ice` (or zero if `q_ice < ϵ`), clamped to `[0, 1]`.

# Arguments
- `q_rim`: rime mass concentration `[kg rim / m³ air]`.
- `q_ice`: total ice mass concentration `[kg ice / m³ air]`.
- `ϵ`: cutoff on `q_ice` below which the fraction is zero. By default, `1e-24 kg/m³`.

# Extended help

The cutoff is independent of `eps(typeof(q_ice))`, which spans nine orders of magnitude
between `Float32` and `Float64` and would place the result at a different physical scale in each.
It sits below every admissible ice content, so a small `q_ice` paired with a
proportionally small `q_rim` returns its true fraction.
"""
@inline rime_mass_fraction(q_rim, q_ice, ϵ = oftype(q_ice, 1e-24)) =
    clamp(_regularised_ratio(q_rim, q_ice, ϵ), zero(q_ice), one(q_ice))

"""
    rime_density(q_rim, b_rim, ϵ = 1e-30)

Regularised rime density `ρ_rim = q_rim / b_rim` (or zero if `b_rim < ϵ`).

# Arguments
- `q_rim`: rime mass concentration `[kg rim / m³ air]`.
- `b_rim`: rime volume concentration `[m³ rim / m³ air]`.
- `ϵ`: cutoff on `b_rim` below which the density is zero. By default, `1e-30 m³/m³`.

# Extended help

The cutoff sits below every admissible rime volume and well above the smallest `Float32`
normal, `1.2e-38`, so a small `b_rim` paired with a proportionally small `q_rim` returns
its true density. Independent of `eps(typeof(b_rim))` for the reason given in
[`rime_mass_fraction`](@ref).
"""
@inline rime_density(q_rim, b_rim, ϵ = oftype(b_rim, 1e-30)) =
    _regularised_ratio(q_rim, b_rim, ϵ)

"""
    guarded_quotient(a, b, absent = zero(a / one(b)))

`a / b` when `b` is present, `absent` otherwise, with correct derivatives in both branches.

# Extended help

Presence is tested on `FD.value(b)`, not on `b`: `ForwardDiff` orders `Dual` numbers
lexicographically, so a `Dual` with a zero value and nonzero partials compares greater than
zero, and a `b > zero(b)` predicate would take the present branch. The test is at exactly
zero rather than a floor on `b`, whose derivative underflows at `b = 0` instead of
reporting an absent quantity.

`absent` defaults to zero; pass it where the absent-state answer differs, 
such as `Inf` for a bound with nothing to constrain it.
"""
@inline function guarded_quotient(a, b, absent = zero(a / one(b)))
    present = FD.value(b) > zero(FD.value(b))
    bs = ifelse(present, b, one(b))
    return ifelse(present, a / bs, absent)
end

"""
    nearest_admissible_b(q, b_trial, ρ_min, ρ_max)
    nearest_admissible_b(q, b_trial, (ρ_min, ρ_max))

The value nearest `b_trial` whose ratio `q / b` lies in `[ρ_min, ρ_max]`, for a nonnegative `q`. 

Where the ratio is already in range that is `b_trial` itself; 
where it is not, the nearer endpoint of `[q / ρ_max, q / ρ_min]`; 
and where `q` is not positive, zero.

# Extended help

`(q, b)` is projected onto the cone `{q, b ≥ 0 : ρ_min b ≤ q ≤ ρ_max b}` as a pair.
Clamping the two independently can leave a positive `q` whose `b` has clipped to zero, or the reverse.
"""
@inline function nearest_admissible_b(q, b_trial, ρ_min, ρ_max)
    return ifelse(q > 0, clamp(b_trial, q / ρ_max, q / ρ_min), zero(q))
end
@inline nearest_admissible_b(q, b_trial, (ρ_min, ρ_max)) =
    nearest_admissible_b(q, b_trial, ρ_min, ρ_max)

end # module
