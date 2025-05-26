"""
    DistributionTools

A module containing tools for working with size distributions.

Currently, it contains tools for working with the generalized gamma distribution
and the exponential distribution.

Mostly used for the 2-moment microphysics.
"""
module DistributionTools

import SpecialFunctions as SF

"""
    generalized_gamma_quantile(ν, μ, B, Y)

Calculate the quantile (inverse cumulative distribution function) for a 
    generalized gamma distribution parameterized in the form:

    g(x) = A ⋅ x^ν ⋅ exp(-B ⋅ x^μ)

# Arguments
- `ν, μ, B`: The PDF parameters
- `Y`: The probability level (0 ≤ Y ≤ 1) for which to compute the quantile

# Returns
- `x`: The value x such that P(X ≤ x) = Y
"""
function generalized_gamma_quantile(ν, μ, B, Y)
    # Compute the inverse of the regularized incomplete gamma function
    z = SF.gamma_inc_inv((ν + 1) / μ, Y, 1 - Y)
    return (z / B)^(1 / μ)
end

"""
    generalized_gamma_cdf(ν, μ, B, x)

Calculate the cumulative distribution function (CDF) for a generalized gamma distribution
parameterized in the form:

    g(x) = A * x^ν * exp(-B * x^μ)

The CDF gives the probability P(X ≤ x).

# Arguments
 - `ν, μ, B`: The PDF parameters
 - `x`: The value at which to evaluate the CDF

# Returns
 - `p`: The probability P(X ≤ x)
"""
function generalized_gamma_cdf(ν, μ, B, x)
    # Check input validity
    μ > 0 || throw(DomainError(μ, "Parameter μ must be positive"))
    B > 0 || throw(DomainError(B, "Parameter B must be positive"))
    # Handle edge cases
    x ≤ 0 && return zero(x)

    # Compute the regularized incomplete gamma function
    p, _ = SF.gamma_inc((ν + 1) / μ, (B * x)^μ)
    return p
end

"""
   exponential_cdf(D_mean, D)

Calculate the cumulative distribution function (CDF) for an exponential distribution
parameterized in the form:

   n(D) = N₀ * exp(-D / D_mean)

where N₀ is a normalizing constant such that the total probability is 1.

# Arguments
- `D_mean`: The mean value of the distribution
- `D`: The point at which to evaluate the CDF (must be ≥ 0)

# Returns
- `p`: The probability P(X ≤ D)
"""
function exponential_cdf(D_mean, D)
    # Check input validity
    D_mean > 0 || throw(DomainError(D_mean, "Mean parameter must be positive"))
    # Handle edge cases
    D < 0 && return zero(D)
    # Calculate CDF: P(X ≤ D) = 1 - exp(-D/D_mean)
    return 1 - exp(-D / D_mean)
end

"""
    exponential_quantile(D_mean, Y)

Calculate the quantile (inverse cumulative distribution function) for an exponential distribution
parameterized in the form:

    n(D) = N₀ * exp(-D / D_mean)

where N₀ is a normalizing constant such that the total probability is 1.

# Arguments
- `Y`: The probability level (0 ≤ Y ≤ 1) for which to compute the quantile
- `D_mean`: The mean value of the distribution

# Returns
- `D`: The value D such that P(X ≤ D) = Y
"""
function exponential_quantile(D_mean, Y)
    # Check input validity
    (0 ≤ Y ≤ 1) || throw(DomainError(Y, "Probability Y must be in [0,1]"))
    D_mean > 0 || throw(DomainError(D_mean, "Mean parameter must be positive"))
    # Calculate quantile: x = -D_mean * ln(1-Y)
    return -D_mean * log(1 - Y)
end

end # module DistributionTools
