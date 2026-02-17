"""
    Utilities

Lightweight numerical utility functions shared across CloudMicrophysics modules.
Contains pure numerical operations with no physics dependencies.
"""
module Utilities

export clamp_to_nonneg, ϵ_numerics, ϵ_numerics_2M_M, ϵ_numerics_2M_N

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

end # module
