import Test as TT
import ForwardDiff as FD
import SpecialFunctions as SF
import CloudMicrophysics.Utilities as UT

# Direct CPU validation of the fast incomplete-gamma approximations
# (`UT.gamma_inc` / `UT.gamma_inc_inv`) against `SpecialFunctions`, and of their
# analytic AD rules. The GPU counterpart lives in `test/gpu_tests.jl`.

# Finite-difference derivative check (see testing_and_validation.md §"AD compatibility tests").
# The AD rule must always return a finite value; the finite-difference reference
# is only accurate enough to compare against at Float64 (a Float32 central
# difference with step `√eps(Float32) ≈ 3e-4` is too noisy), matching the
# Float64-only comparison used in `test/p3_tests.jl`.
function check_derivative(f, x; rtol, atol)
    ad = FD.derivative(f, x)
    TT.@test isfinite(ad)
    if typeof(x) == Float64
        ε = sqrt(eps(typeof(x)))
        fd = (f(x + ε) - f(x - ε)) / (2ε)
        TT.@test isapprox(ad, fd; rtol, atol)
    end
end

function test_gamma_inc(FT)
    TT.@testset "gamma_inc / gamma_inc_inv [FT=$FT]" begin
        # `p` values are dyadic so `q = 1 - p` is exact and `SF.gamma_inc_inv`
        # (which requires `p + q == 1`) accepts the reference unchanged.
        avals = FT[1, 1.5, 2, 2.5, 3.5, 5, 7.5]
        xvals = FT[0.1, 0.5, 1, 2.5, 5, 8, 12]
        pvals = FT[0.03125, 0.125, 0.25, 0.5, 0.75, 0.875, 0.96875]

        # native-FT approximation is less precise than the Float64-backed SF
        atol_PQ = FT == Float32 ? FT(2e-5) : FT(1e-6)
        rtol_inv = FT == Float32 ? FT(2e-4) : FT(1e-5)

        TT.@testset "accuracy vs SpecialFunctions" begin
            for a in avals, x in xvals
                P_sf, Q_sf = SF.gamma_inc(a, x)
                P_ut, Q_ut = UT.gamma_inc(a, x)
                TT.@test isapprox(P_ut, P_sf; atol = atol_PQ)
                TT.@test isapprox(Q_ut, Q_sf; atol = atol_PQ)
            end
            for a in avals, p in pvals
                q = FT(1) - p
                TT.@test isapprox(UT.gamma_inc_inv(a, p, q), SF.gamma_inc_inv(a, p, q);
                    rtol = rtol_inv, atol = rtol_inv)
            end
        end

        TT.@testset "analytic AD rules (x- and p-derivative)" begin
            # ∂P/∂x = x^(a-1) e^-x / Γ(a); the inverse uses dx/dp = 1 / (∂P/∂x)
            for a in avals, x in xvals
                check_derivative(x -> UT.gamma_inc(a, x)[1], x; rtol = FT(1e-3), atol = FT(1e-5))
            end
            for a in avals, p in pvals
                check_derivative(p -> UT.gamma_inc_inv(a, p, FT(1) - p), p;
                    rtol = FT(1e-3), atol = FT(1e-4))
            end
        end

        TT.@testset "shape-parameter (`a`) derivative is rejected" begin
            # Differentiating w.r.t. the shape parameter is unsupported and must
            # error (not silently return a zero gradient).
            TT.@test_throws ErrorException FD.derivative(a -> UT.gamma_inc(a, FT(3))[1], FT(2.5))
            TT.@test_throws ErrorException FD.derivative(a -> UT.gamma_inc_inv(a, FT(0.4), FT(0.6)), FT(3))
            # A Dual `a` with zero partials (a promoted constant) is allowed.
            Tag = typeof(FD.Tag(identity, FT))
            a0 = FD.Dual{Tag}(FT(2.5), zero(FT))
            TT.@test UT.gamma_inc(a0, FT(3))[1] isa FD.Dual
            x1 = FD.Dual{Tag}(FT(3), one(FT))
            TT.@test FD.partials(UT.gamma_inc(a0, x1)[1])[1] != 0  # x-derivative still flows
        end
    end
end

TT.@testset "Incomplete gamma utilities" begin
    for FT in (Float32, Float64)
        test_gamma_inc(FT)
    end
end
nothing
