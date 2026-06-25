"""
    ADSpecialFunctions

AD-compatible, GPU-lean reimplementations of the special functions that
CloudMicrophysics needs: `gamma`, `loggamma`, `gamma_inc`, `gamma_inc_inv`,
`beta_inc`, `erf`, `erfc`.

Each is pure generic arithmetic (no `SpecialFunctions` internals and no `mpfr`
fallback, with only a couple of value branches), so `ForwardDiff.Dual` numbers
flow through directly (no derivative rules needed) and GPU compilation stays
small. The algorithm and term count are chosen at compile time by `precision(T)`,
which returns `SinglePrecision()` or `DoublePrecision()` and forwards through
`ForwardDiff.Dual` to its value type. A `Dual{…,Float32}` therefore runs the
single-precision algorithm in `Dual` arithmetic.

The functions target the argument ranges CloudMicrophysics actually uses
(gamma-PSD shape/moment orders and `λD` magnitudes, aerosol-activation
`erf`/`erfc`), not the full domain of a general-purpose library.

Accuracy (vs `SpecialFunctions`, measured over those ranges): `Float64`
results agree to ~1e-13 relative and analytic derivatives to ~1e-15;
`Float32` results agree to ~1e-5 relative in the normal range. A result
whose true value falls below `floatmin(Float32)` is returned as `0`, which
is the only representable `Float32` answer and matches the flush-to-zero
behaviour of typical GPU `Float32` arithmetic.
"""
module ADSpecialFunctions

import ForwardDiff

export precision, SinglePrecision, DoublePrecision

# ---------------------------------------------------------------------------
# Precision markers: a compile-time choice of algorithm/term-count.
# ---------------------------------------------------------------------------
abstract type Precision end
struct SinglePrecision <: Precision end
struct DoublePrecision <: Precision end

@inline precision(::Type{Float32}) = SinglePrecision()
@inline precision(::Type{Float64}) = DoublePrecision()
# Dual forwards to its value type V (the underlying float decides the algorithm).
@inline precision(::Type{ForwardDiff.Dual{T, V, N}}) where {T, V, N} = precision(V)
@inline precision(::Type{T}) where {T <: Real} = DoublePrecision()
@inline precision(x) = precision(typeof(x))

@inline _underlying(::SinglePrecision) = Float32
@inline _underlying(::DoublePrecision) = Float64

# Lentz continued-fraction floor (FPMIN). Chosen so `tiny²` stays well above
# `floatmin`: a derivative through `1/d` is `-d'/d²`, so an unclamped `d` near
# the floor would otherwise underflow `d²` and make the derivative non-finite.
@inline _fpmin(::SinglePrecision) = 1.0f-18
@inline _fpmin(::DoublePrecision) = 1.0e-30

# Compile-time iteration limits that bound the series/continued-fraction loops.
@inline _gamma_maxit(::SinglePrecision) = 40
@inline _gamma_maxit(::DoublePrecision) = 300
@inline _beta_maxit(::SinglePrecision) = 60
@inline _beta_maxit(::DoublePrecision) = 350
@inline _erf_maxit(::SinglePrecision) = 35
@inline _erf_maxit(::DoublePrecision) = 80
# With the convergence break, these are safety bounds; Halley needs ~4 steps.
@inline _inv_maxit(::SinglePrecision) = 12
@inline _inv_maxit(::DoublePrecision) = 20

# ---------------------------------------------------------------------------
# loggamma / gamma  (Lanczos g = 7), valid for real arguments via reflection.
# ---------------------------------------------------------------------------
const _LANCZOS = (
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
)

@inline function _loggamma_pos(z::T) where {T}
    # z > 0
    g = T(7)
    x = z - one(T)
    a = T(_LANCZOS[1])
    @inbounds for i in 2:9
        a += T(_LANCZOS[i]) / (x + T(i - 1))
    end
    t = x + g + T(1) / 2
    (x + one(T) / 2) * log(t) - t + T(0.5) * log(T(2) * T(π)) + log(a)
end

"""
    loggamma(z)

Compute `log(Γ(z))` for real `z` (via reflection for `z < 0.5`).
"""
@inline function loggamma(z::T) where {T}
    if z < T(1) / 2
        # reflection: Γ(z)Γ(1-z) = π / sin(πz)  ⇒  logΓ(z) = log(π/|sin(πz)|) - logΓ(1-z)
        return log(T(π) / abs(sin(T(π) * z))) - _loggamma_pos(one(T) - z)
    else
        return _loggamma_pos(z)
    end
end
@inline loggamma(z::Integer) = loggamma(float(z))

"""
    gamma(z)

Compute `Γ(z)` for real `z`.
"""
@inline function gamma(z::T) where {T}
    if z < T(1) / 2
        return T(π) / (sin(T(π) * z) * exp(_loggamma_pos(one(T) - z)))
    else
        return exp(_loggamma_pos(z))
    end
end
@inline gamma(z::Integer) = gamma(float(z))

# ---------------------------------------------------------------------------
# Regularized incomplete gamma  P(a,x), Q(a,x)
#   series for P when x < a+1, continued fraction for Q otherwise.
# ---------------------------------------------------------------------------
@inline function _P_series(a::T, x::T, n::Int) where {T}
    tol = eps(_underlying(precision(T)))
    ap = a
    del = one(T) / a
    s = del
    @inbounds for _ in 1:n
        ap += one(T)
        del *= x / ap
        s += del
        abs(del) <= abs(s) * tol && break
    end
    s * exp(-x + a * log(x) - _loggamma_pos(a))
end

@inline function _Q_cf(a::T, x::T, n::Int) where {T}
    tol = eps(_underlying(precision(T)))
    tiny = T(_fpmin(precision(T)))
    b = x + one(T) - a
    c = one(T) / tiny
    d = one(T) / b
    h = d
    @inbounds for i in 1:n
        an = -T(i) * (T(i) - a)
        b += T(2)
        d = an * d + b
        abs(d) < tiny && (d = tiny)
        c = b + an / c
        abs(c) < tiny && (c = tiny)
        d = one(T) / d
        del = d * c
        h *= del
        abs(del - one(T)) <= tol && break
    end
    exp(-x + a * log(x) - _loggamma_pos(a)) * h
end

"""
    gamma_inc(a, x)

Compute the regularized lower/upper incomplete gamma functions
`(P(a,x), Q(a,x))` (matches `SpecialFunctions.gamma_inc`). `a > 0`, `x ≥ 0`.
"""
@inline function gamma_inc(a, x)
    T = promote_type(typeof(a), typeof(x))
    _gamma_inc(T(a), T(x), precision(T))
end
@inline function _gamma_inc(a::T, x::T, p::Precision) where {T}
    n = _gamma_maxit(p)
    x <= zero(T) && return (zero(T), one(T))
    if x < a + one(T)
        P = _P_series(a, x, n)
        return (P, one(T) - P)
    else
        Q = _Q_cf(a, x, n)
        return (one(T) - Q, Q)
    end
end

# ---------------------------------------------------------------------------
# Regularized incomplete beta  I_x(a,b)  (continued fraction, Lentz).
# ---------------------------------------------------------------------------
@inline function _betacf(a::T, b::T, x::T, n::Int) where {T}
    tol = eps(_underlying(precision(T)))
    tiny = T(_fpmin(precision(T)))
    qab = a + b
    qap = a + one(T)
    qam = a - one(T)
    c = one(T)
    d = one(T) - qab * x / qap
    abs(d) < tiny && (d = tiny)
    d = one(T) / d
    h = d
    @inbounds for m in 1:n
        m2 = T(2m)
        aa = T(m) * (b - T(m)) * x / ((qam + m2) * (a + m2))
        d = one(T) + aa * d
        abs(d) < tiny && (d = tiny)
        c = one(T) + aa / c
        abs(c) < tiny && (c = tiny)
        d = one(T) / d
        h *= d * c
        aa = -(a + T(m)) * (qab + T(m)) * x / ((a + m2) * (qap + m2))
        d = one(T) + aa * d
        abs(d) < tiny && (d = tiny)
        c = one(T) + aa / c
        abs(c) < tiny && (c = tiny)
        d = one(T) / d
        del = d * c
        h *= del
        abs(del - one(T)) <= tol && break
    end
    h
end

"""
    beta_inc(a, b, x)

Compute the regularized incomplete beta `(I_x(a,b), 1 - I_x(a,b))`
(matches `SpecialFunctions.beta_inc`). `a,b > 0`, `0 ≤ x ≤ 1`.
"""
@inline function beta_inc(a, b, x)
    T = promote_type(typeof(a), typeof(b), typeof(x))
    _beta_inc(T(a), T(b), T(x), precision(T))
end
@inline function _beta_inc(a::T, b::T, x::T, p::Precision) where {T}
    x <= zero(T) && return (zero(T), one(T))
    x >= one(T) && return (one(T), zero(T))
    n = _beta_maxit(p)
    bt = exp(
        loggamma(a + b) - loggamma(a) - loggamma(b) +
        a * log(x) + b * log(one(T) - x),
    )
    if x < (a + one(T)) / (a + b + T(2))
        I = bt * _betacf(a, b, x, n) / a
        return (I, one(T) - I)
    else
        Ic = bt * _betacf(b, a, one(T) - x, n) / b
        return (one(T) - Ic, Ic)
    end
end

# ---------------------------------------------------------------------------
# Inverse regularized lower incomplete gamma: x with P(a,x) = p.
#   value via Halley on the AD-able P; AD flows through the iteration, so
#   ∂x/∂p AND ∂x/∂a are obtained for free (no IFT / no FD-of-a needed).
# ---------------------------------------------------------------------------
# Rational normal-quantile approx (Acklam) — feeds only the initial guess.
@inline function _norm_quantile(p::T) where {T}
    a1 = T(-3.969683028665376e+01);
    a2 = T(2.209460984245205e+02)
    a3 = T(-2.759285104469687e+02);
    a4 = T(1.383577518672690e+02)
    a5 = T(-3.066479806614716e+01);
    a6 = T(2.506628277459239e+00)
    b1 = T(-5.447609879822406e+01);
    b2 = T(1.615858368580409e+02)
    b3 = T(-1.556989798598866e+02);
    b4 = T(6.680131188771972e+01)
    b5 = T(-1.328068155288572e+01)
    c1 = T(-7.784894002430293e-03);
    c2 = T(-3.223964580411365e-01)
    c3 = T(-2.400758277161838e+00);
    c4 = T(-2.549732539343734e+00)
    c5 = T(4.374664141464968e+00);
    c6 = T(2.938163982698783e+00)
    d1 = T(7.784695709041462e-03);
    d2 = T(3.224671290700398e-01)
    d3 = T(2.445134137142996e+00);
    d4 = T(3.754408661907416e+00)
    plow = T(0.02425)
    if p < plow
        q = sqrt(-T(2) * log(p))
        return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
               ((((d1 * q + d2) * q + d3) * q + d4) * q + one(T))
    elseif p <= one(T) - plow
        q = p - T(0.5);
        r = q * q
        return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
               (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + one(T))
    else
        q = sqrt(-T(2) * log(one(T) - p))
        return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
               ((((d1 * q + d2) * q + d3) * q + d4) * q + one(T))
    end
end

@inline function _gii_guess(a::T, p::T, lga::T) where {T}
    if a > one(T)
        t = _norm_quantile(p)
        g = a * (one(T) - one(T) / (T(9) * a) + t / (T(3) * sqrt(a)))^3
        return g > zero(T) ? g : (p * a * exp(lga))^(one(T) / a)
    else
        return (p * a * exp(lga))^(one(T) / a)   # P(a,x) ≈ x^a/Γ(a+1) for small x
    end
end

"""
    gamma_inc_inv(a, p)

Invert the regularized lower incomplete gamma: return `x` with `P(a,x) = p`
(value matches `SpecialFunctions.gamma_inc_inv(a, p, 1 - p)`). `a > 0`,
`0 < p < 1`. AD-able in both `a` and `p`.
"""
@inline function gamma_inc_inv(a, p)
    T = promote_type(typeof(a), typeof(p))
    _gamma_inc_inv(T(a), T(p), precision(T))
end
@inline function _gamma_inc_inv(a::T, p::T, pr::Precision) where {T}
    tol = eps(_underlying(precision(T)))
    lga = _loggamma_pos(a)
    x = _gii_guess(a, p, lga)
    @inbounds for _ in 1:_inv_maxit(pr)
        f = _gamma_inc(a, x, pr)[1] - p
        logfp = (a - one(T)) * log(x) - x - lga
        fp = exp(logfp)                 # P'(a,x) > 0
        dx = f / fp
        fr = ((a - one(T)) / x - one(T))   # f''/f'
        step = dx / (one(T) - T(0.5) * dx * fr)   # Halley step
        x -= step
        # Halley converges cubically, so a relative step of √eps leaves a
        # remaining error of ~eps^(3/2) (machine precision); √eps also sits
        # above the inner `gamma_inc` noise floor, so the break fires reliably.
        abs(step) <= sqrt(tol) * abs(x) && break
    end
    x
end

# ---------------------------------------------------------------------------
# erf / erfc  (series near 0, continued fraction for the tail; odd symmetry).
# ---------------------------------------------------------------------------
@inline function _erf_series(x::T, n::Int) where {T}
    tol = eps(_underlying(precision(T)))
    x2 = x * x
    term = x
    s = term
    @inbounds for k in 1:n
        term *= -x2 / T(k)
        add = term / T(2k + 1)
        s += add
        abs(add) <= abs(s) * tol && break
    end
    T(2) / sqrt(T(π)) * s
end
@inline function _erfc_cf(x::T, n::Int) where {T}   # x > 0
    tol = eps(_underlying(precision(T)))
    tiny = T(_fpmin(precision(T)))
    h = x;
    c = x;
    d = zero(T)
    @inbounds for i in 1:n
        a = T(i) / T(2)
        d = x + a * d
        abs(d) < tiny && (d = tiny)
        c = x + a / c
        abs(c) < tiny && (c = tiny)
        d = one(T) / d
        del = c * d
        h *= del
        abs(del - one(T)) <= tol && break
    end
    exp(-x * x) / (sqrt(T(π)) * h)
end

"""
    erf(x)

Compute the error function `erf(x)` for real `x`.
"""
@inline erf(x::T) where {T <: Real} = _erf(x, precision(T))
@inline erf(x::Integer) = erf(float(x))
@inline function _erf(x::T, p::Precision) where {T}
    if abs(x) < T(2)
        return _erf_series(x, _erf_maxit(p))
    else
        e = _erfc_cf(abs(x), _erf_maxit(p))
        return x < zero(T) ? (e - one(T)) : (one(T) - e)
    end
end

"""
    erfc(x)

Compute the complementary error function `erfc(x) = 1 - erf(x)` for real `x`.
"""
# erfc via the continued fraction for the tail (avoids 1 - erf cancellation).
@inline function _erfc(x::T, p::Precision) where {T}
    n = _erf_maxit(p)
    if x >= T(2)
        return _erfc_cf(x, n)
    elseif x <= -T(2)
        return T(2) - _erfc_cf(-x, n)
    else
        return one(T) - _erf_series(x, n)
    end
end
@inline erfc(x::T) where {T <: Real} = _erfc(x, precision(T))
@inline erfc(x::Integer) = erfc(float(x))

end # module
