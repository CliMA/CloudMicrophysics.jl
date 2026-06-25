import Test as TT

import CloudMicrophysics.ADSpecialFunctions as ADSF
import SpecialFunctions as SF
import ForwardDiff as FD

# Combined absolute/relative error. A reference below `floor` (e.g. a `Float32`
# underflow) only has `≈ 0` as a representable answer, so compare absolutely.
relerr(v, r; floor = 0.0) =
    abs(r) <= floor ? abs(v) : abs(v - r) / max(abs(r), eps(typeof(float(r))))

TT.@testset "ADSpecialFunctions: precision dispatch" begin
    TT.@test ADSF.precision(Float32) === ADSF.SinglePrecision()
    TT.@test ADSF.precision(Float64) === ADSF.DoublePrecision()
    TT.@test ADSF.precision(1.0f0) === ADSF.SinglePrecision()
    TT.@test ADSF.precision(1.0) === ADSF.DoublePrecision()
    # A Dual forwards to its value type: Dual{…,Float32} runs the F32 algorithm.
    TT.@test ADSF.precision(FD.Dual{:t}(1.0f0, 1.0f0)) === ADSF.SinglePrecision()
    TT.@test ADSF.precision(FD.Dual{:t}(1.0, 1.0)) === ADSF.DoublePrecision()
end

TT.@testset "ADSpecialFunctions: Float64 correctness vs SpecialFunctions" begin
    # loggamma / gamma (positive and, via reflection, negative arguments)
    for z in (0.6, 1.0, 2.5, 7.3, 30.0, 123.0)
        TT.@test ADSF.loggamma(z) ≈ SF.loggamma(z) atol = 1e-11 rtol = 1e-12
        TT.@test ADSF.gamma(z) ≈ SF.gamma(z) rtol = 1e-12
    end
    for z in (-0.6, -1.3, -2.7, -4.4)   # avoid the integer poles
        TT.@test ADSF.gamma(z) ≈ SF.gamma(z) rtol = 1e-11
    end

    # incomplete gamma P, Q over a wide (a, x) grid
    for a in 10.0 .^ range(-1, 2, length = 12), x in 10.0 .^ range(-3, 2, length = 12)
        P, Q = ADSF.gamma_inc(a, x)
        Pr, Qr = SF.gamma_inc(a, x)
        TT.@test relerr(P, Pr; floor = 1e-290) < 1e-11
        TT.@test relerr(Q, Qr; floor = 1e-290) < 1e-11
        TT.@test P + Q ≈ 1
    end

    # incomplete beta
    for a in (0.5, 2.0, 8.0), b in (0.7, 3.0, 12.0), x in (0.05, 0.3, 0.6, 0.95)
        I, Ic = ADSF.beta_inc(a, b, x)
        Ir, Icr = SF.beta_inc(a, b, x)
        TT.@test relerr(I, Ir; floor = 1e-290) < 1e-11
        TT.@test I + Ic ≈ 1
    end

    # inverse incomplete gamma
    for a in (0.7, 2.0, 5.0, 20.0), p in (1e-4, 0.05, 0.3, 0.7, 0.99)
        TT.@test ADSF.gamma_inc_inv(a, p) ≈ SF.gamma_inc_inv(a, p, 1 - p) rtol = 1e-9
    end

    # erf / erfc (incl. the large-|x| tail where 1 - erf would cancel)
    for x in (-6.0, -2.5, -1.0, -0.1, 0.3, 1.5, 3.0, 7.0)
        TT.@test ADSF.erf(x) ≈ SF.erf(x) rtol = 1e-12
        TT.@test ADSF.erfc(x) ≈ SF.erfc(x) rtol = 1e-11
    end
    TT.@test ADSF.erf(0.0) == 0.0
end

TT.@testset "ADSpecialFunctions: AD vs analytic derivatives" begin
    # d/dz loggamma = digamma ; d/dz gamma = gamma * digamma
    for z in (0.7, 2.5, 12.0)
        TT.@test FD.derivative(ADSF.loggamma, z) ≈ SF.digamma(z) rtol = 1e-11
        TT.@test FD.derivative(ADSF.gamma, z) ≈ SF.gamma(z) * SF.digamma(z) rtol = 1e-11
    end

    # d/dx P(a,x) = x^(a-1) e^{-x} / Γ(a)   (the gamma pdf)
    for (a, x) in ((2.0, 1.0), (5.0, 3.0), (0.7, 0.2), (10.0, 8.0))
        d = FD.derivative(xx -> ADSF.gamma_inc(a, xx)[1], x)
        TT.@test d ≈ x^(a - 1) * exp(-x) / SF.gamma(a) rtol = 1e-11
    end

    # d/dx I_x(a,b) = x^(a-1)(1-x)^(b-1) / B(a,b)
    for (a, b, x) in ((2.0, 3.0, 0.4), (5.0, 2.0, 0.6))
        d = FD.derivative(xx -> ADSF.beta_inc(a, b, xx)[1], x)
        TT.@test d ≈ x^(a - 1) * (1 - x)^(b - 1) / SF.beta(a, b) rtol = 1e-11
    end

    # d/dx erf = 2/√π e^{-x²}; derivative is AD-correct at the series origin x=0
    c = 2 / sqrt(pi)
    for x in (-3.0, -0.5, 0.5, 1.5, 4.0)
        TT.@test FD.derivative(ADSF.erf, x) ≈ c * exp(-x^2) rtol = 1e-11
        TT.@test FD.derivative(ADSF.erfc, x) ≈ -c * exp(-x^2) rtol = 1e-11
    end
    TT.@test FD.derivative(ADSF.erf, 0.0) ≈ c rtol = 1e-12

    # gamma_inc_inv: d/dp x = 1/P'(a,x) = Γ(a) x^{1-a} e^{x} (inverse-fn theorem)
    for (a, p) in ((2.0, 0.3), (5.0, 0.7), (12.0, 0.5))
        x = ADSF.gamma_inc_inv(a, p)
        d = FD.derivative(pp -> ADSF.gamma_inc_inv(a, pp), p)
        TT.@test d ≈ SF.gamma(a) * x^(1 - a) * exp(x) rtol = 1e-9
    end

    # nested 2nd-order AD: d²/dx² P(a,x) = ((a-1)/x - 1) · pdf
    let a = 5.0, x = 3.0
        d2 = FD.derivative(xx -> FD.derivative(yy -> ADSF.gamma_inc(a, yy)[1], xx), x)
        pdf = x^(a - 1) * exp(-x) / SF.gamma(a)
        TT.@test d2 ≈ ((a - 1) / x - 1) * pdf rtol = 1e-10
    end
end

TT.@testset "ADSpecialFunctions: shape-parameter AD (∂/∂a, which SF cannot do)" begin
    # SpecialFunctions has no ∂/∂a rule and `gamma_inc_inv` has no `Dual` method
    # at all; here AD produces both. Reference: central finite difference.
    fd(f, a; h = 1e-6) = (f(a + h) - f(a - h)) / (2h)
    for (a, x) in ((2.0, 1.0), (5.0, 3.0), (8.0, 5.0))
        ad = FD.derivative(aa -> ADSF.gamma_inc(aa, x)[1], a)
        TT.@test ad ≈ fd(aa -> ADSF.gamma_inc(aa, x)[1], a) rtol = 1e-5
    end
    for (a, p) in ((2.0, 0.3), (5.0, 0.7))
        ad = FD.derivative(aa -> ADSF.gamma_inc_inv(aa, p), a)
        TT.@test ad ≈ fd(aa -> ADSF.gamma_inc_inv(aa, p), a) rtol = 1e-5
    end
    # `gamma_inc_inv` is genuinely differentiable here (SF errors on a Dual)
    TT.@test FD.derivative(aa -> ADSF.gamma_inc_inv(aa, 0.4), 3.0) isa Float64
end

TT.@testset "ADSpecialFunctions: Float32 correctness and eltype" begin
    # eltype is preserved (no silent promotion to Float64)
    TT.@test ADSF.gamma_inc(2.0f0, 1.0f0) isa Tuple{Float32, Float32}
    TT.@test ADSF.gamma_inc_inv(2.0f0, 0.3f0) isa Float32
    TT.@test ADSF.beta_inc(2.0f0, 3.0f0, 0.4f0) isa Tuple{Float32, Float32}
    TT.@test ADSF.erf(1.0f0) isa Float32
    TT.@test ADSF.gamma(2.5f0) isa Float32
    # a Dual{…,Float32} stays Float32-valued
    TT.@test ADSF.gamma_inc(FD.Dual{:t}(2.0f0, 1.0f0), 1.0f0)[1] isa FD.Dual{:t, Float32, 1}

    # Float32 accuracy in the normal range (true value ≥ floatmin(Float32))
    fmin = floatmin(Float32)
    for a in Float32.(10 .^ range(-1, 1.4, length = 8)),
        x in Float32.(10 .^ range(-1.5, 1.6, length = 10))

        r = SF.gamma_inc(Float64(a), Float64(x))[1]
        r < fmin && continue
        TT.@test relerr(Float64(ADSF.gamma_inc(a, x)[1]), r) < 1e-4
    end
end

TT.@testset "ADSpecialFunctions: type stability and allocations" begin
    TT.@inferred ADSF.gamma_inc(2.0, 1.0)
    TT.@inferred ADSF.gamma_inc(2.0f0, 1.0f0)
    TT.@inferred ADSF.gamma_inc_inv(2.0, 0.3)
    TT.@inferred ADSF.beta_inc(2.0, 3.0, 0.4)
    TT.@inferred ADSF.erf(1.0)
    TT.@inferred ADSF.erfc(1.0)
    TT.@inferred ADSF.loggamma(2.5)
    TT.@inferred ADSF.gamma(2.5)

    # scalar, branch-light kernels allocate nothing
    g() = ADSF.gamma_inc(5.0, 3.0)
    h() = ADSF.gamma_inc_inv(5.0, 0.3)
    g();
    h()
    TT.@test (@allocated g()) == 0
    TT.@test (@allocated h()) == 0
end
