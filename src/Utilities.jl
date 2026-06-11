"""
    Utilities

Lightweight numerical utility functions shared across CloudMicrophysics modules.
Contains pure numerical operations with no physics dependencies.
"""
module Utilities

import UnrolledUtilities as UU

export clamp_to_nonneg, ϵ_numerics, ϵ_numerics_2M_M, ϵ_numerics_2M_N, ϵ_numerics_P3_B
export unrolled_logsumexp
export sgs_weight_function, rime_mass_fraction, rime_density

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
    return ifelse(
        denominator < ϵ, zero(numerator),
        weight * numerator / denominator,
    )
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
