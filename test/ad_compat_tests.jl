using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.DistributionTools as DT
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ForwardDiff as FD
import SpecialFunctions as SF

# ForwardDiff compatibility of the pointwise 2M+P3 path: `bulk_microphysics_tendencies`
# must be differentiable w.r.t. the 8 prognostic species (with `logλ` held fixed,
# matching the substepping semantics). Differentiating *through* the `logλ` shape
# solve is out of scope here — it additionally requires a `∂/∂a` rule for the
# forward `SF.gamma_inc`.

function test_ad_compatibility(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    D(v, ∂) = FD.Dual{:ad_test}(FT(v), FT(∂))

    @testset "P3State construction with Dual prognostics ($FT)" begin
        st = P3.state_from_prognostic(p3, FT(1e-4), FT(1e4), FT(2e-5), FT(4e-8))
        std = P3.state_from_prognostic(p3, D(1e-4, 1), D(1e4, 0), D(2e-5, 0), D(4e-8, 0))
        @test eltype(std) <: FD.Dual
        # primal values are unchanged by differentiation
        @test FD.value(std.ρ_g) == st.ρ_g
        @test FD.value(std.D_gr) == st.D_gr
        @test FD.value(std.D_cr) == st.D_cr
        @test FD.value(std.D_th) == st.D_th
        # params-only field is a true constant; rime-derived fields carry sensitivity
        @test iszero(FD.partials(std.D_th))
        @test !iszero(FD.partials(std.ρ_g))
        @test !iszero(FD.partials(std.D_cr))
        # partial seeding (derivative w.r.t. a single prognostic) also constructs
        st1 = P3.state_from_prognostic(p3, D(1e-4, 1), FT(1e4), FT(2e-5), FT(4e-8))
        @test eltype(st1) <: FD.Dual
        # the `F_rim = 0` sentinel branch stays intact under Duals
        st0 = P3.state_from_prognostic(p3, D(1e-4, 1), D(1e4, 0), D(0, 0), D(0, 0))
        @test FD.value(st0.D_gr) == FT(Inf) && FD.value(st0.D_cr) == FT(Inf)
    end

    @testset "regularised ratios stay NaN-free across tiny denominators ($FT)" begin
        # below ~eps(FT)/4 the sgs_weight_function sigmoid hits atanh(-1):
        # the weight is 0 either way, but the partials were NaN. Sweep
        # denominators across that band for both regularised ratios.
        for denom in (eps(FT)^2, eps(FT) / 8, eps(FT), sqrt(eps(FT)), FT(1e-9))
            std = P3.state_from_prognostic(p3, D(denom, 1), D(10, 0), D(denom / 10, 1), D(denom / 10, 1))
            # the regularised ratios must always be differentiable
            @test all(isfinite, FD.partials(std.F_rim))
            @test all(isfinite, FD.partials(std.ρ_rim))
            # the cached thresholds are NaN/Inf sentinels in degenerate rime
            # regimes (gated downstream); they need finite partials only where
            # their value is finite
            for field in (std.ρ_g, std.D_gr, std.D_cr)
                isfinite(FD.value(field)) && @test all(isfinite, FD.partials(field))
            end
        end
    end

    @testset "Γ_incl accepts mixed argument types ($FT)" begin
        # the rain-evaporation path calls Γ_incl(params_float, dual): Microphysics2M ~:808
        g = CM2.Γ_incl(FT(-0.25), D(0.5, 1))
        @test g isa FD.Dual
        @test FD.value(g) ≈ CM2.Γ_incl(FT(-0.25), FT(0.5))
        @test CM2.Γ_incl(FT(-0.25), FT(0.5)) isa FT
    end

    @testset "BMT 2M+P3 Jacobian w.r.t. the 8 species ($FT)" begin
        function rhs(x, ρ, T, q_tot, logλ)
            t = BMT.bulk_microphysics_tendencies(
                BMT.Microphysics2Moment(), mp, tps, ρ, T, q_tot,
                x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], logλ)
            return [t.dq_lcl_dt, t.dn_lcl_dt, t.dq_rai_dt, t.dn_rai_dt,
                t.dq_ice_dt, t.dn_ice_dt, t.dq_rim_dt, t.db_rim_dt]
        end
        function consistent_logλ(ρ, x)
            st = P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8])
            return P3.get_distribution_logλ(st)
        end
        # x = [q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim]; interior
        # states (all species nonzero where the regime is active)
        regimes = (
            (; name = "warm rain", ρ = FT(1.05), T = FT(288), q_tot = FT(0.015),
                x = FT[4e-4, 8e7, 2.1e-3, 5e4, 0, 0, 0, 0], logλ = FT(-Inf)),
            (; name = "mixed phase", ρ = FT(0.78), T = FT(273.5), q_tot = FT(0.009),
                x = FT[2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8], logλ = nothing),
            (; name = "ice heavy", ρ = FT(0.45), T = FT(233), q_tot = FT(0.003),
                x = FT[1e-6, 1e6, 1e-12, 1e-2, 8e-4, 5e5, 5e-4, 9e-7], logλ = nothing),
            # sub-threshold ice with b_rim in the regularised-ratio band that
            # previously produced NaN partials via sgs_weight_function
            (; name = "cloud edge", ρ = FT(0.7), T = FT(263), q_tot = FT(0.005),
                x = FT[1e-5, 1e7, 1e-6, 1e3, 3e-8, 30, 1e-8, 2.5e-11], logλ = nothing),
        )
        for r in regimes
            logλ = isnothing(r.logλ) ? consistent_logλ(r.ρ, r.x) : r.logλ
            f = x -> rhs(x, r.ρ, r.T, r.q_tot, logλ)
            v₀ = f(r.x)
            J = FD.jacobian(f, r.x)
            @test all(isfinite, J)
            @test f(r.x) == v₀  # differentiation does not perturb the primal
        end

        # Jacobian against central finite differences (Float64 only — FD
        # truncation in Float32 is not meaningful at these magnitudes)
        if FT == Float64
            r = regimes[2]
            logλ = consistent_logλ(r.ρ, r.x)
            f = x -> rhs(x, r.ρ, r.T, r.q_tot, logλ)
            J = FD.jacobian(f, r.x)
            J_fd = similar(J)
            for j in 1:8
                h = 1e-6 * r.x[j]
                xp = copy(r.x);
                xp[j] += h
                xm = copy(r.x);
                xm[j] -= h
                J_fd[:, j] = (f(xp) - f(xm)) / 2h
            end
            # per-row scales: number rows dwarf mass rows by ~10 orders of
            # magnitude, so a single global scale would leave the mass rows
            # unconstrained
            for i in 1:8
                scale = max(maximum(abs, J[i, :]), maximum(abs, J_fd[i, :]))
                iszero(scale) && continue
                @test maximum(abs, J[i, :] - J_fd[i, :]) / scale < 1e-5
            end
        end

        # Boundary: SB2006 autoconversion Φ_au(τ) ∝ τ^0.7 has a vertical
        # tangent at exactly zero rain with cloud present; the ϵ-gate in
        # `autoconversion` keeps the Jacobian finite there
        x_boundary = FT[1e-6, 1e6, 0, 0, 8e-4, 5e5, 5e-4, 9e-7]
        logλ_b = consistent_logλ(FT(0.45), x_boundary)
        f_b = x -> rhs(x, FT(0.45), FT(233), FT(0.003), logλ_b)
        @test all(isfinite, f_b(x_boundary))
        @test all(isfinite, FD.jacobian(f_b, x_boundary))
    end
end

@testset "gamma_inc_inv Dual rule vs finite differences" begin
    # SF kernels are Float64-only internally; validate the rule in Float64
    D64(v, ∂) = FD.Dual{:ad_test}(Float64(v), Float64(∂))
    for (a, p) in ((0.5, 0.1), (1.3, 0.1), (2.5, 1e-5), (2.5, 0.5), (7.0, 1 - 1e-5), (25.0, 0.9))
        x = SF.gamma_inc_inv(a, p, 1 - p)
        # p-slot (callers pass (Y, 1-Y), so q carries the opposite partial)
        g = DT.gamma_inc_inv(D64(a, 0), D64(p, 1), D64(1 - p, -1))
        h = p * (1 - p) * 1e-4
        ∂p_fd = (SF.gamma_inc_inv(a, p + h, 1 - p - h) - SF.gamma_inc_inv(a, p - h, 1 - p + h)) / 2h
        @test FD.value(g) == x  # primal exactly the SF result
        @test FD.partials(g)[1] ≈ ∂p_fd rtol = 1e-5
        # a-slot
        ga = DT.gamma_inc_inv(D64(a, 1), D64(p, 0), D64(1 - p, 0))
        ha = 1e-5 * a
        ∂a_fd = (SF.gamma_inc_inv(a + ha, p, 1 - p) - SF.gamma_inc_inv(a - ha, p, 1 - p)) / 2ha
        @test FD.partials(ga)[1] ≈ ∂a_fd rtol = 1e-4
        # q-slot symmetry: seeding only q composes to the negative p-derivative
        gq = DT.gamma_inc_inv(D64(a, 0), D64(p, 0), D64(1 - p, 1))
        @test FD.partials(gq)[1] ≈ -FD.partials(g)[1] / 2 atol = abs(FD.partials(g)[1]) * 1e-10
        # quantile through the public wrapper stays consistent
        @test DT.gamma_inc_inv(a, p, 1 - p) == x
    end

    # a-slot deep in the upper tail: P saturates at the resolution of one(x),
    # so the rule must difference Q there (the repo's upper integration bounds
    # use q as small as eps())
    for q in (1e-8, 1e-13, eps())
        a = 3.0
        ga = DT.gamma_inc_inv(D64(a, 1), D64(1 - q, 0), D64(q, 0))
        ha = 1e-5 * a
        ∂a_fd = (SF.gamma_inc_inv(a + ha, 1 - q, q) - SF.gamma_inc_inv(a - ha, 1 - q, q)) / 2ha
        @test FD.partials(ga)[1] ≈ ∂a_fd rtol = 1e-4
    end

    # mixed Dual/plain arguments promote inside the public wrapper
    gm = DT.gamma_inc_inv(2.5, D64(0.5, 1), D64(0.5, -1))
    @test gm isa FD.Dual
    @test FD.value(gm) == SF.gamma_inc_inv(2.5, 0.5, 0.5)
end

@testset "get_distribution_logλ_from_prognostic" begin
    for FT in (Float32, Float64)
        p3 = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true).ice.scheme
        logλ = P3.get_distribution_logλ_from_prognostic(p3, FT(1e-4), FT(1e4), FT(2e-5), FT(4e-8))
        st = P3.state_from_prognostic(p3, FT(1e-4), FT(1e4), FT(2e-5), FT(4e-8))
        @test isfinite(logλ)
        @test logλ == P3.get_distribution_logλ(st)
    end
end

test_ad_compatibility(Float64)
test_ad_compatibility(Float32)
