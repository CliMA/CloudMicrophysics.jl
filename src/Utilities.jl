"""
    Utilities

Lightweight numerical utility functions shared across CloudMicrophysics modules.
Contains pure numerical operations with no physics dependencies.
"""
module Utilities

import UnrolledUtilities as UU

export clamp_to_nonneg, ϵ_numerics, ϵ_numerics_2M_M, ϵ_numerics_2M_N, ϵ_numerics_P3_B
export unrolled_logsumexp

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


end # module
