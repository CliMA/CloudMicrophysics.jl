
"""
    integrate(f, a, b, quad = ChebyshevGauss(100))

 Approximate the definite integral ∫ₐᵇ f(x) dx using Chebyshev-Gauss quadrature of the first kind.

# Mathematical Background

 This method transforms the integral to the standard interval [-1, 1] and applies 
 Chebyshev-Gauss quadrature. The transformation is:
 
     y = (2x - (a+b)) / (b-a)    →    x = (b-a)y/2 + (a+b)/2
 
 with Jacobian dx/dy = (b-a)/2.
 
 The integral becomes:

     ∫ₐᵇ f(x) dx = (b-a)/2 ∫₋₁¹ f(x(y)) dy
 
 Using the Chebyshev-Gauss quadrature identity:

     ∫₋₁¹ g(y) dy = ∫₋₁¹ [g(y)√(1-y²)] / √(1-y²) dy ≈ (π/n) ∑ᵢ₌₁ⁿ g(yᵢ)√(1-yᵢ²)
 
 where yᵢ = cos((2i-1)π/(2n)) are the Chebyshev nodes of the first kind.
 
 # Final Formula

     ∫ₐᵇ f(x) dx ≈ (b-a)π/(2n) ∑ᵢ₌₁ⁿ f(xᵢ)√(1-yᵢ²)
 
 where:
 - yᵢ = cos((2i-1)π/(2n)) for i = 1, ..., n
 - xᵢ = (b-a)yᵢ/2 + (a+b)/2 
 - All quadrature weights are equal: π/n

# Arguments
 - `f`: Function to integrate
 - `a`, `b`: Integration bounds. Note: if `a ≥ b`, or `a` or `b` is `NaN`, `zero(f(a))` is returned.

# Keyword arguments
 - `quad`: Quadrature scheme, default: `ChebyshevGauss(100)`

# Returns
 Approximation to the definite integral ∫ₐᵇ f(x) dx

# Notes
 This method achieves spectral convergence for smooth functions and is particularly 
 effective for analytic functions. The √(1-y²) weighting factor helps handle 
 functions with mild singularities at the interval endpoints.
 Ref: https://en.wikipedia.org/wiki/Chebyshev–Gauss_quadrature
"""
function integrate(f, a, b; quad = ChebyshevGauss(100))
    FT = eltype(float(a))
    # Pre-compute transformation parameters
    scale_factor = (b - a) / 2
    shift = (a + b) / 2

    # Compute integral using Chebyshev-Gauss quadrature
    result = zero(f(a))
    a < b || return result  # return early unless a < b
    (; n) = quad
    for i in 1:n
        # Node on [-1, 1] interval
        y = node(quad, FT(i), n)

        # Node on [a, b] interval
        x = scale_factor * y + shift

        # Total weight √(1 - y²) * wᵢ
        w = inv_weight_fun(quad, y) * weight(quad, FT(i), n)

        # Accumulate: f(x) * √(1 - y²) * wᵢ
        result += f(x) * w
    end

    return scale_factor * result
end

"""
    integrate(f, bnds...; quad = ChebyshevGauss(100))

Integrate the function `f` over each subinterval of the integration bounds, `bnds`.

# Arguments
 - `f`: Function to integrate
 - `bnds`: A tuple of bounds, `(a, b, c, d, ...)`
 - `quad`: Quadrature scheme, default: `ChebyshevGauss(100)`

 The integral is computed as the sum of the integrals over each subinterval,
 `(a, b), (b, c), (c, d), ...`.
"""
function integrate(f, bnds...; quad = ChebyshevGauss(100))
    # compute integral over each subinterval (a, b), (b, c), (c, d), ...
    return sum(integrate(f, a, b; quad) for (a, b) in zip(Base.front(bnds), Base.tail(bnds)))
end

"""
    ChebyshevGauss(n)

Quadrature scheme for Chebyshev-Gauss quadrature of the first kind.

# Arguments
 - `n`: Number of quadrature points

# Available methods

    node(::ChebyshevGauss, i::FT, n) where {FT}
    weight(::ChebyshevGauss, i::FT, n) where {FT}
    inv_weight_fun(::ChebyshevGauss, y)

- `node(quad, i, n)`: Return the `i`-th node of the `n`-point Chebyshev-Gauss quadrature scheme.
- `weight(quad, i, n)`: Return the `i`-th weight of the `n`-point Chebyshev-Gauss quadrature scheme.
- `inv_weight_fun(quad, y)`: Return the inverse of the weight function `w(x)`.



# Mathematical Background

A method to approximate the value of integrals of the kind

    ∫_{-1}^{1} f(x) w(x) dx ≈ ∑_1^n f(x_i) wᵢ(x_i)

where `w(x) = 1 / √(1 - x^2)` is the weight function, `x_i` are the nodes, and `wᵢ` are the weights.

If we are interested in only the integral of `f(x)`, we can instead integrate `g(x) = f(x) / w(x)`,

    ∫_{-1}^{1} g(x) w(x) dx ≈ ∑_1^n g(x_i) wᵢ(x_i) = ∑_1^n f(x_i) wᵢ(x_i) / w(x_i)

# References

- https://en.wikipedia.org/wiki/Chebyshev–Gauss_quadrature
- https://en.wikipedia.org/wiki/Gaussian_quadrature
"""
struct ChebyshevGauss
    n::Int
end
Base.broadcastable(quad::ChebyshevGauss) = (quad,)

@inline node(::ChebyshevGauss, i::FT, n) where {FT} = cospi((2float(FT(i)) - 1) / (2n))
@inline weight(::ChebyshevGauss, i::FT, n) where {FT} = float(FT(π)) / n
@inline inv_weight_fun(::ChebyshevGauss, y) = √(1 - y^2)

"""
    integral_bounds(state::P3State, logλ; p, moment_order = 0)

Compute the integration bounds for the P3 size distribution,

    Mⁿ = ∫_a^b Dⁿ * N′(D) dD = N₀ * ∫_a^b Dⁿ D^μ * exp(-λ * D) dD

 where `Mⁿ` is the `n`-th moment of the size distribution.
 Here `n ≡ moment_order` and `a` and `b` are the integration bounds.
 For a proper moment, `a=0` and `b=∞`. For the numerical integration, `a` and `b`
 are determined by this function.

# Arguments
- `state`: [`P3State`](@ref) object
- `logλ`: The log of the slope parameter [log(1/m)]
- `p`: The integration bounds are set to the `p`-th and `1-p`-th quantiles of the size distribution.
- `moment_order`: For integrands proportional to moments of the size distribution, 
    `moment_order` can be used to indicate the order of the moment. 
    May provide more accurate bounds; thus more accurate integration.
    Default: `moment_order = 0`.

# Returns
- `bnds`: The integration bounds (a `Tuple`), for use in numerical integration (c.f. [`integrate`](@ref)).
"""
function integral_bounds(state::P3State{FT}, logλ; p, moment_order = 0) where {FT}
    # Get reduced lower and upper bounds from quantiles
    k = get_μ(state, logλ) + moment_order
    D_min = DT.generalized_gamma_quantile(k, FT(1), exp(logλ), FT(p))
    D_max = DT.generalized_gamma_quantile(k, FT(1), exp(logλ), FT(1 - p))

    # Only integrate up to the maximum diameter, `D_max`, including intermediate thresholds
    # If `F_rim` is very close to 1, `D_cr` may be greater than `D_max`, in which case it is disregarded.
    thresholds = get_bounded_thresholds(state, D_min, D_max)
    return thresholds
end

"""
    D_m(state::P3State, logλ)

Compute the mass weighted mean particle size [m]

# Parameters
 - `state`: [`P3State`](@ref) object
 - `logλ`: The log of the slope parameter [log(1/m)]
"""
function D_m(state, logλ)
    μ = get_μ(state, logλ)
    mass_weighted_moment = logmass_gamma_moment(state, μ, logλ; n = 1)
    log_N₀ = get_logN₀(state.N_ice, μ, logλ)
    return exp(log_N₀ + mass_weighted_moment) / state.L_ice
end
