"""
    Quadrature

Numerical quadrature rules for the P3 size-distribution integrals.

Two rules are provided, both plugging into the same [`integrate`](@ref)
machinery:

- [`ChebyshevGauss`](@ref): Chebyshev-Gauss of the first kind. Nodes/weights
  have a closed form, so the rule is cheap to build for *any* order and carries
  only `n::Int`.
- [`GaussLegendre`](@ref): Gauss-Legendre. Nodes/weights have **no** closed
  form; they are computed once at construction (host side, `O(n)` via
  FastGaussQuadrature) and stored in the rule as `SVector`s. The rule is
  `isbits` and ships to GPU kernels unchanged — only *construction* is
  host-side.

The generalized Gauss-Laguerre rule [`GaussLaguerre`](@ref) integrates over `[0, ∞)` against a
normalized gamma weight, for the drop freezing integrals over a size distribution.

This module lives early in the include order (before `Parameters`) so that the
constructed quadrature object can be stored on a parameter struct and shipped to
the device, rather than being rebuilt inside a GPU kernel.
"""
module Quadrature

import StaticArrays as SA
import UnrolledUtilities as UU
import FastGaussQuadrature as FGQ
import ..Utilities as UT

export QuadratureRule, ChebyshevGauss, GaussLegendre, GaussLaguerre, integrate

abstract type QuadratureRule end

"""
    integrate(f, a, b, quad)

 Approximate the definite integral ∫ₐᵇ f(x) dx using the quadrature rule `quad`.

# Mathematical Background

 The integral is transformed to the standard interval [-1, 1]:

     y = (2x - (a+b)) / (b-a)    →    x = (b-a)y/2 + (a+b)/2

 with Jacobian dx/dy = (b-a)/2, so that

     ∫ₐᵇ f(x) dx = (b-a)/2 ∫₋₁¹ f(x(y)) dy
                 ≈ (b-a)/2 ∑ᵢ f(xᵢ) · inv_weight_fun(quad, yᵢ) · weight(quad, i)

 where `yᵢ = node(quad, i)`, `xᵢ = (b-a)yᵢ/2 + (a+b)/2`, and the
 `inv_weight_fun`/`weight` factors are rule-specific (see [`ChebyshevGauss`](@ref)
 and [`GaussLegendre`](@ref)).

# Arguments
 - `f`: Function to integrate
 - `a`, `b`: Integration bounds. Note: if `a ≥ b`, or `a` or `b` is `NaN`, zero is returned.
 - `quad`: Quadrature scheme (a `QuadratureRule`).

# Returns
 Approximation to the definite integral ∫ₐᵇ f(x) dx
"""
@inline function integrate(f::F, a::T, b::T, quad::QuadratureRule) where {F, T}
    # Pre-compute transformation parameters
    scale_factor = (b - a) / 2
    shift = (a + b) / 2

    # Compute integral using the quadrature rule
    result = zero(UT.return_type_of(f, a))
    a < b || return result  # return early unless a < b
    (; n) = quad
    for i in 1:n
        # Node on [-1, 1] interval
        y = node(quad, i)

        # Node on [a, b] interval
        x = scale_factor * y + shift

        # Total weight: rule-specific weight-function factor × wᵢ
        w = inv_weight_fun(quad, y) * weight(quad, i)

        # Accumulate: f(x) * inv_weight_fun(y) * wᵢ
        result += f(x) * w
    end

    return scale_factor * result
end

"""
    subintervals(bnds::NTuple)

Pair adjacent elements of a flat boundary tuple into the consecutive
subintervals they delimit:

    subintervals((a, b, c, d)) ≡ ((a, b), (b, c), (c, d))

`Base.front` and `Base.tail` are zero-cost on fixed-length tuples, so
the helper compiles to a direct tuple construction with no runtime
dispatch.
"""
@inline subintervals(bnds::NTuple) =
    UU.unrolled_map(tuple, Base.front(bnds), Base.tail(bnds))

"""
    integrate(f, bnds, quad)

Integrate the function `f` over each subinterval of the integration bounds, `bnds`.

# Arguments
 - `f`: Function to integrate
 - `bnds`: A tuple of bounds, `(a, b, c, d, ...)`
 - `quad`: Quadrature scheme (a `QuadratureRule`).

 The integral is computed as the sum of the integrals over each subinterval,
 `(a, b), (b, c), (c, d), ...` — see [`subintervals`](@ref).
"""
@inline function integrate(f::F, bnds::NTuple{N, T}, quad::QuadratureRule) where {F, N, T}
    result = zero(UT.return_type_of(f, bnds[1]))
    @inbounds for i in 1:(N - 1)
        result += integrate(f, bnds[i], bnds[i + 1], quad)
    end
    return result
end

# ---------------------------------------------------------------------------
# Chebyshev-Gauss quadrature (closed-form nodes/weights)
# ---------------------------------------------------------------------------

"""
    ChebyshevGauss(n)
    ChebyshevGauss(FT, n)

Quadrature scheme for Chebyshev-Gauss quadrature of the first kind.

# Arguments
 - `FT`: float type of the nodes and weights; default `Float64`
 - `n`: Number of quadrature points

# Available methods

    node(::ChebyshevGauss, i)
    weight(::ChebyshevGauss, i)
    inv_weight_fun(::ChebyshevGauss, y)

- `node(quad, i)`: Return the `i`-th node of the Chebyshev-Gauss quadrature scheme.
- `weight(quad, i)`: Return the `i`-th weight of the Chebyshev-Gauss quadrature scheme.
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
struct ChebyshevGauss{FT} <: QuadratureRule
    n::Int
end
ChebyshevGauss(n::Int) = ChebyshevGauss(Float64, n)
ChebyshevGauss(::Type{FT}, n::Int) where {FT} = ChebyshevGauss{FT}(n)
Base.broadcastable(quad::ChebyshevGauss) = (quad,)

@inline node((; n)::ChebyshevGauss{FT}, i) where {FT} = cospi((2 * FT(i) - 1) / (2n))
@inline weight((; n)::ChebyshevGauss{FT}, i) where {FT} = FT(π) / n
@inline inv_weight_fun(::ChebyshevGauss, y) = √(1 - y^2)

# ---------------------------------------------------------------------------
# Gauss-Legendre quadrature (nodes/weights computed once at construction)
# ---------------------------------------------------------------------------

"""
    GaussLegendre(n)
    GaussLegendre(FT, n)

Quadrature scheme for Gauss-Legendre quadrature, plugging into the same
[`integrate`](@ref) machinery as [`ChebyshevGauss`](@ref).

Unlike `ChebyshevGauss` (whose nodes/weights have a closed form), the
Gauss-Legendre nodes/weights have **no closed form**. They are computed **once,
at construction**, in `Float64` via `FastGaussQuadrature.gausslegendre` (an
`O(n)` host-side routine) and converted to `FT`. The constructed rule stores the
nodes/weights as `SVector{N, FT}`, so per-`integrate` access is a static lookup
— identical hot-path cost to `ChebyshevGauss`. Construction is intended to run
**once at initialization** (e.g. when building a parameter struct), never inside
a GPU kernel.

Arbitrary orders `n ≥ 1` are supported (the order is the type parameter `N`).

# GPU / type-stability

`GaussLegendre{FT, N}` stores `SVector{N, FT}` nodes/weights, so it is `isbits`
(GPU-safe) and fully concrete. The order `N` is a type parameter, fixed at
construction; the object is built host-side once and shipped to the device, so
`N` never ripples a runtime value into a kernel.

# Accuracy vs `ChebyshevGauss`

At matched `n`, Gauss-Legendre is more accurate on the smooth P3
size-distribution integrals (e.g. ~20× lower error than `ChebyshevGauss(40)` on
the dominant ice-rain collision integral, verified against an adaptive QuadGK
reference) at equal cost, and it is the default rule on the P3 parameters.
`ice_self_collection` places the diagonal `D_1 = D_2` on a subinterval boundary,
so its integrand is smooth on each subinterval and Gauss-Legendre converges
geometrically there too.

# Available methods

    node(::GaussLegendre, i)
    weight(::GaussLegendre, i)
    inv_weight_fun(::GaussLegendre, y)   # ≡ 1 (Legendre weight is 1 on [-1, 1])

# References

- https://en.wikipedia.org/wiki/Gauss–Legendre_quadrature
- A. Townsend, *FastGaussQuadrature.jl* (MIT); Bogaert (2014) O(n) nodes/weights.
"""
struct GaussLegendre{FT, N} <: QuadratureRule
    n::Int
    nodes::SA.SVector{N, FT}
    weights::SA.SVector{N, FT}
end
Base.broadcastable(quad::GaussLegendre) = (quad,)

@inline node(q::GaussLegendre, i) = @inbounds q.nodes[i]
@inline weight(q::GaussLegendre, i) = @inbounds q.weights[i]
@inline inv_weight_fun(::GaussLegendre, y) = one(y)

# Host-side, one-shot construction path (see the `GaussLegendre` docstring);
# not meant to be called inside a GPU kernel.
function GaussLegendre(::Type{FT}, n::Int) where {FT}
    n ≥ 1 || error("GaussLegendre: order n=$n must be ≥ 1")
    # Compute in Float64 for accuracy, then convert to FT.
    nodes64, weights64 = FGQ.gausslegendre(n)
    nodes = SA.SVector{n, FT}(ntuple(i -> FT(nodes64[i]), n))
    weights = SA.SVector{n, FT}(ntuple(i -> FT(weights64[i]), n))
    return GaussLegendre{FT, n}(n, nodes, weights)
end
GaussLegendre(n::Int) = GaussLegendre(Float64, n)

"""
    GaussLaguerre{FT, N}
    GaussLaguerre(FT, n; α = 0)

The `n`-point generalized Gauss-Laguerre rule on `[0, ∞)`,

    ∫₀^∞ tᵅ exp(-t) / Γ(α + 1) g(t) dt ≈ Σᵢ wᵢ g(tᵢ),

whose weight function is normalized, so the weights sum to one. The nodes and weights are
computed once at construction with FastGaussQuadrature and stored as `SVector`s.

# Fields
- `α`: exponent of the weight function [-]
- `nodes`, `weights`: the nodes `tᵢ` and the normalized weights `wᵢ`
"""
struct GaussLaguerre{FT, N}
    α::FT
    nodes::SA.SVector{N, FT}
    weights::SA.SVector{N, FT}
end
Base.broadcastable(quad::GaussLaguerre) = (quad,)

function GaussLaguerre(::Type{FT}, n::Int; α = 0) where {FT}
    n ≥ 1 || error("GaussLaguerre: order n=$n must be ≥ 1")
    # Compute in Float64 for accuracy, then convert to FT
    nodes64, weights64 = FGQ.gausslaguerre(n, Float64(α))
    w_sum = sum(weights64)
    nodes = SA.SVector{n, FT}(ntuple(i -> FT(nodes64[i]), n))
    weights = SA.SVector{n, FT}(ntuple(i -> FT(weights64[i] / w_sum), n))
    return GaussLaguerre{FT, n}(FT(α), nodes, weights)
end

end # module Quadrature
