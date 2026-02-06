"""
    Utilities

Lightweight numerical utility functions shared across CloudMicrophysics modules.
Contains pure numerical operations with no physics dependencies.
"""
module Utilities

export clamp_to_nonneg, ϵ_numerics, ϵ_numerics_2M_M, ϵ_numerics_2M_N

"""
    clamp_to_nonneg(x)

Clamp values to be non-negative in an automatic differentiation (AD) compatible way.
Uses `ifelse` to preserve dual number types during AD.

# Arguments
- `x`: value to clamp

# Returns
- `x` if `x ≥ 0`, otherwise `0*x` (preserving the type of x)
"""
@inline clamp_to_nonneg(x) = ifelse(x < 0, 0 * x, x)

"""
    ϵ_numerics(FT)

Smallest number that is different than zero for the purpose of microphysics
computations. Returns `sqrt(floatmin(FT))` to avoid underflow issues.
"""
@inline ϵ_numerics(FT) = sqrt(floatmin(FT))

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

end # module
