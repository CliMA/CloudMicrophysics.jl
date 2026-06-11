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

This module lives early in the include order (before `Parameters`) so that the
constructed quadrature object can be stored on a parameter struct and shipped to
the device, rather than being rebuilt inside a GPU kernel.
"""
module Quadrature

import StaticArrays as SA
import UnrolledUtilities as UU
import FastGaussQuadrature as FGQ

export QuadratureRule, ChebyshevGauss, GaussLegendre, integrate

abstract type QuadratureRule end

"""
    integrate(f, a, b, quad = ChebyshevGauss(100))

 Approximate the definite integral ∫ₐᵇ f(x) dx using the quadrature rule `quad`.

# Mathematical Background

 The integral is transformed to the standard interval [-1, 1]:

     y = (2x - (a+b)) / (b-a)    →    x = (b-a)y/2 + (a+b)/2

 with Jacobian dx/dy = (b-a)/2, so that

     ∫ₐᵇ f(x) dx = (b-a)/2 ∫₋₁¹ f(x(y)) dy
                 ≈ (b-a)/2 ∑ᵢ f(xᵢ) · inv_weight_fun(quad, yᵢ) · weight(quad, i, n)

 where `yᵢ = node(quad, i, n)`, `xᵢ = (b-a)yᵢ/2 + (a+b)/2`, and the
 `inv_weight_fun`/`weight` factors are rule-specific (see [`ChebyshevGauss`](@ref)
 and [`GaussLegendre`](@ref)).

# Arguments
 - `f`: Function to integrate
 - `a`, `b`: Integration bounds. Note: if `a ≥ b`, or `a` or `b` is `NaN`, `zero(f(a))` is returned.

# Keyword arguments
 - `quad`: Quadrature scheme, default: `ChebyshevGauss(100)`

# Returns
 Approximation to the definite integral ∫ₐᵇ f(x) dx
"""
@inline function integrate(f::F, a::T, b::T, quad::QuadratureRule = ChebyshevGauss(100)) where {F, T}
    FT = eltype(float(a))
    # Pre-compute transformation parameters
    scale_factor = (b - a) / 2
    shift = (a + b) / 2

    # Compute integral using the quadrature rule
    result = zero(f(a))
    a < b || return result  # return early unless a < b
    (; n) = quad
    for i in 1:n
        # Node on [-1, 1] interval
        y = node(quad, FT(i), n)

        # Node on [a, b] interval
        x = scale_factor * y + shift

        # Total weight: rule-specific weight-function factor × wᵢ
        w = inv_weight_fun(quad, y) * weight(quad, FT(i), n)

        # Accumulate: f(x) * inv_weight_fun(y) * wᵢ
        result += f(x) * w
    end

    return scale_factor * result
end

"""
    subintervals(bnds::NTuple{N, T}) -> NTuple{N-1, NTuple{2, T}}

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
    integrate(f, bnds, quad = ChebyshevGauss(100))

Integrate the function `f` over each subinterval of the integration bounds, `bnds`.

# Arguments
 - `f`: Function to integrate
 - `bnds`: A tuple of bounds, `(a, b, c, d, ...)`
 - `quad`: Quadrature scheme, default: `ChebyshevGauss(100)`

 The integral is computed as the sum of the integrals over each subinterval,
 `(a, b), (b, c), (c, d), ...` — see [`subintervals`](@ref).
"""
@inline function integrate(
    f::F, bnds::NTuple{N, T}, quad::QuadratureRule = ChebyshevGauss(100),
) where {F, N, T}
    pairs = UU.unrolled_map(subintervals(bnds)) do (a, b)
        integrate(f, a, b, quad)
    end
    return UU.unrolled_sum(pairs)
end

# ---------------------------------------------------------------------------
# Chebyshev-Gauss quadrature (closed-form nodes/weights)
# ---------------------------------------------------------------------------

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
struct ChebyshevGauss <: QuadratureRule
    n::Int
end
Base.broadcastable(quad::ChebyshevGauss) = (quad,)

@inline node(::ChebyshevGauss, i::FT, n) where {FT} = cospi((2float(FT(i)) - 1) / (2n))
@inline weight(::ChebyshevGauss, i::FT, n) where {FT} = float(FT(π)) / n
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
`40` is the ClimaAtmos production `quadrature_order`.

# GPU / type-stability

`GaussLegendre{FT, N}` stores `SVector{N, FT}` nodes/weights, so it is `isbits`
(GPU-safe) and fully concrete. The order `N` is a type parameter, fixed at
construction; the object is built host-side once and shipped to the device, so
`N` never ripples a runtime value into a kernel.

# Accuracy vs `ChebyshevGauss`

At matched `n`, Gauss-Legendre is substantially more accurate on the smooth P3
size-distribution integrals (e.g. ~20× lower error than `ChebyshevGauss(40)` on
the dominant ice-rain collision integral, verified against an adaptive QuadGK
reference) at equal cost. It is **not** uniformly better on the cusp-limited
`ice_self_collection` diagonal — there both schemes are quadrature-limited at low
`n` and the structural remedy is a diagonal split, not the scheme. The
production default therefore remains `ChebyshevGauss`; switch a call site to
`GaussLegendre` deliberately where the integrand is smooth.

# Available methods

    node(::GaussLegendre, i, n)
    weight(::GaussLegendre, i, n)
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

@inline node(q::GaussLegendre, i::FT, n) where {FT} = @inbounds q.nodes[Int(i)]
@inline weight(q::GaussLegendre, i::FT, n) where {FT} = @inbounds q.weights[Int(i)]
@inline inv_weight_fun(::GaussLegendre, y) = one(y)

"""
    GaussLegendre(FT, n)

Construct the `n`-point Gauss-Legendre rule with element type `FT`.

Nodes/weights are computed in `Float64` (`FastGaussQuadrature.gausslegendre`)
and converted to `FT`. This is the host-side, one-shot construction path; it is
**not** meant to be called inside a GPU kernel.
"""
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
    build_quadrature(FT, quadrature_order)

Select and **construct** the quadrature rule for the P3 size-distribution
integrals from a single `quadrature_order` knob, in element type `FT`. This is
the host-side, one-shot builder; the returned object is `isbits` and stored on a
parameter struct for reuse in the (GPU) hot loop.

Gauss-Legendre is preferred for the orders where it is meaningfully more
accurate than Chebyshev-Gauss on the smooth P3 integrands (≈20× lower error on
the dominant ice-rain collision integral at matched `n`; see [`GaussLegendre`](@ref)),
namely `quadrature_order ∈ {16, 32, 40, 64}` (incl. the ClimaAtmos production
order 40). Any other order falls back to [`ChebyshevGauss`](@ref), preserving the
default behaviour for non-preferred orders.
"""
function build_quadrature(::Type{FT}, quadrature_order::Int) where {FT}
    return if quadrature_order in (16, 32, 40, 64)
        GaussLegendre(FT, quadrature_order)
    else
        ChebyshevGauss(quadrature_order)
    end
end

end # module Quadrature
