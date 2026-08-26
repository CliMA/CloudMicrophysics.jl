using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Common as CO
import CloudMicrophysics.DistributionTools as DT
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Microphysics1M as CM1
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
        # the `F_rim = 0` branch stays intact under Duals
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
            # the cached thresholds are NaN/Inf markers in degenerate rime
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

    @testset "rain_evaporation is concretely typed for mixed arguments ($FT)" begin
        # the early return (no rain / supersaturated) must have the same type
        # as the main path when species are Duals over a plain-float q_tot —
        # a union here heap-boxes every call in the Jacobian hot loop
        sb = mp.warm_rain.seifert_beheng
        aps = mp.warm_rain.air_properties
        z = D(0, 0)
        # subsaturated (main path) and supersaturated (early return) states
        for (q_tot, T) in ((FT(0.005), FT(288)), (FT(0.02), FT(288)))
            t = @inferred CM2.rain_evaporation(
                sb, aps, tps, q_tot, D(2e-4, 1), z, D(1e-4, 1), z, FT(1.05), D(4e4, 1), T,
            )
            @test all(v -> v isa FD.Dual, values(t))
        end
    end

    @testset "terminal velocity promotes mixed coefficient tuples ($FT)" begin
        # A Dual air or ice density makes only some Chen 2022 coefficients Dual
        vel = mp.ice.terminal_velocity
        ρₐ = D(1.2, 1)
        for v_term in (
            CO.particle_terminal_velocity(vel.rain, ρₐ),
            CO.particle_terminal_velocity(vel.small_ice, ρₐ, FT(916.7)),
            CO.particle_terminal_velocity(vel.large_ice, ρₐ, FT(916.7)),
            CO.particle_terminal_velocity(vel.small_ice, ρₐ, D(916.7, 0)),
        )
            vt = v_term(FT(1e-3))
            @test vt isa FD.Dual
            @test isfinite(FD.value(vt))
        end
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

# ForwardDiff compatibility of the 1-moment process rates. The 1M kernels take
# parameters (plain floats) and state alongside each other, so a `::FT ...
# where {FT}` signature would bind both to one type and reject a `Dual`
# state; the rates are duck-typed instead, with `zero(rate)` fallbacks and
# value-keyed `ϵ_numerics` gates.
function test_ad_compatibility_1M(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics1MParams(FT)
    mpS = CMP.Microphysics1MParams(FT; snow_autoconversion = CMP.WithSupersaturation())
    chen = CMP.Chen2022VelType(FT)
    rain, snow = mp.precip.rain, mp.precip.snow
    ρ0, T_cold, T_warm = FT(1.2), FT(270), FT(280)
    state(q) = (; q_tot = q + FT(5e-3), q_lcl = q, q_icl = q, q_rai = q, q_sno = q)

    # (name, f) pairs differentiated w.r.t. a single scalar
    rates_q = (
        (
            "accretion lcl-rai",
            q -> CM1.accretion(CMP.CloudLiquidRainAccretion(), mp, tps, state(q), (; ρ = ρ0, T = T_cold)),
        ),
        (
            "accretion icl-sno",
            q -> CM1.accretion(CMP.CloudIceSnowAccretion(), mp, tps, state(q), (; ρ = ρ0, T = T_cold)),
        ),
        (
            "accretion rain sink",
            q -> CM1.accretion_rain_sink(CMP.CloudIceRainAccretion(), mp, tps, state(q), (; ρ = ρ0, T = T_cold)),
        ),
        (
            "accretion snow-rain",
            q ->
                CM1.accretion_snow_rain(CMP.RainSnowAccretion(), mp, tps, state(q), (; ρ = ρ0, T = T_cold)).S_rai_sno,
        ),
        ("snow melt", q -> CM1.conv_q_sno_to_q_rai(CMP.SnowMelt(), mp, tps, state(q), (; ρ = ρ0, T = T_warm))),
        ("cloud ice melt", q -> CM1.conv_q_icl_to_q_lcl(CMP.CloudIceMelt(), mp, tps, state(q), (; ρ = ρ0, T = T_warm))),
        (
            "rain autoconversion",
            q -> CM1.conv_q_lcl_to_q_rai(CMP.Kessler1M(), mp, tps, state(q), (; ρ = ρ0, T = T_cold)),
        ),
        (
            "snow autoconversion",
            q -> CM1.conv_q_icl_to_q_sno(CMP.NoSupersaturation(), mp, tps, state(q), (; ρ = ρ0, T = T_cold)),
        ),
        (
            "snow autoconversion supersat",
            q -> CM1.conv_q_icl_to_q_sno(CMP.WithSupersaturation(), mpS, tps, state(q), (; ρ = ρ0, T = T_cold)),
        ),
        ("v_t rain blk1m", q -> CM1.terminal_velocity(rain, mp.terminal_velocity.rain, ρ0, q)),
        ("v_t snow blk1m", q -> CM1.terminal_velocity(snow, mp.terminal_velocity.snow, ρ0, q)),
        ("v_t rain Chen", q -> CM1.terminal_velocity(rain, chen.rain, ρ0, q)),
        ("v_t snow Chen", q -> CM1.terminal_velocity(snow, chen.large_ice, ρ0, q)),
        ("v_t snow Oblate", q -> CM1.terminal_velocity(snow, chen.large_ice, ρ0, q, CM1.Oblate())),
        ("v_t snow Prolate", q -> CM1.terminal_velocity(snow, chen.large_ice, ρ0, q, CM1.Prolate())),
        ("lambda_inverse", q -> CM1.lambda_inverse(rain.pdf, rain.mass, q, ρ0)),
    )

    @testset "1M ForwardDiff w.r.t. q, $name ($FT)" for (name, f) in rates_q
        q0 = FT(1e-4)
        d = FD.derivative(f, q0)
        @test isfinite(d)
        # value must be unchanged by seeding a Dual
        @test FD.value(f(FD.Dual{:ad_1m}(q0, one(FT)))) ≈ f(q0) rtol = 10 * eps(FT)
        # Agreement with a central difference. `cbrt(eps)` is the step that
        # balances truncation against round-off. The tolerance is set by the
        # difference quotient's accuracy, not the derivative's: where |f'|
        # is large relative to |f| (the Chen velocities, f' ~ 8e3 with
        # f ~ O(1)) cancellation leaves ~0.2% in Float32. An AD defect (a
        # dropped term or a wrong chain rule) is an O(1) relative error, so
        # 1% still catches every failure mode this test exists to catch.
        # The atol scales with the observed magnitudes so the comparison
        # stays meaningful for small derivatives.
        h = q0 * cbrt(eps(FT))
        fd = (f(q0 + h) - f(q0 - h)) / (2h)
        rtol_fd = FT === Float32 ? FT(1e-2) : FT(1e-4)
        @test d ≈ fd rtol = rtol_fd atol = rtol_fd * max(abs(d), abs(fd))
    end

    # The Dual now arrives through a *different* argument than above, which
    # a single-argument `eltype(q)` would type incorrectly.
    rates_ρ = (
        ("v_t rain blk1m", r -> CM1.terminal_velocity(rain, mp.terminal_velocity.rain, r, FT(1e-4))),
        ("v_t rain Chen", r -> CM1.terminal_velocity(rain, chen.rain, r, FT(1e-4))),
        ("v_t snow Chen", r -> CM1.terminal_velocity(snow, chen.large_ice, r, FT(1e-4))),
        ("v_t snow Oblate", r -> CM1.terminal_velocity(snow, chen.large_ice, r, FT(1e-4), CM1.Oblate())),
        ("lambda_inverse", r -> CM1.lambda_inverse(rain.pdf, rain.mass, FT(1e-4), r)),
        ("get_v0", r -> CM1.get_v0(mp.terminal_velocity.rain, r)),
    )
    @testset "1M ForwardDiff w.r.t. ρ, $name ($FT)" for (name, f) in rates_ρ
        @test isfinite(FD.derivative(f, ρ0))
    end
    @testset "1M ForwardDiff w.r.t. ρ, sign and guards ($FT)" begin
        # denser air slows the fall speed
        @test FD.derivative(r -> CM1.terminal_velocity(rain, mp.terminal_velocity.rain, r, FT(1e-4)), ρ0) < 0
        # the ρ ≥ ρw density guard must not poison the derivative
        ρw = mp.terminal_velocity.rain.ρw
        @test isfinite(FD.derivative(r -> CM1.get_v0(mp.terminal_velocity.rain, r), ρw))
        @test isfinite(FD.derivative(r -> CM1.get_v0(mp.terminal_velocity.rain, r), 2 * ρw))
    end

    # Process rates must stay concretely typed with Dual state and plain
    # thermo (kernel-level argument mixes are swept in return_type_tests.jl).
    @testset "1M rates concretely typed under Dual state ($FT)" begin
        DT = FD.Dual{:ad_1m, FT, 1}
        for (opt, f) in (
            (CMP.CloudLiquidRainAccretion(), CM1.accretion),
            (CMP.RainSnowAccretion(), CM1.accretion_snow_rain),
            (CMP.CloudIceRainAccretion(), CM1.accretion_rain_sink),
            (CMP.RainEvaporation(), CM1.conv_q_rai_to_q_vap),
            (CMP.DepositionAndSublimation(), CM1.conv_q_sno_to_q_vap),
            (CMP.SnowMelt(), CM1.conv_q_sno_to_q_rai),
            (CMP.CloudIceMelt(), CM1.conv_q_icl_to_q_lcl),
            (CMP.Kessler1M(), CM1.conv_q_lcl_to_q_rai),
        )
            rts = Base.return_types(
                f,
                Tuple{
                    typeof(opt),
                    typeof(mp),
                    typeof(tps),
                    NamedTuple{(:q_tot, :q_lcl, :q_icl, :q_rai, :q_sno), NTuple{5, DT}},
                    NamedTuple{(:ρ, :T), NTuple{2, FT}},
                },
            )
            @test length(rts) == 1 && isconcretetype(rts[1])
        end
    end
end

test_ad_compatibility_1M(Float64)
test_ad_compatibility_1M(Float32)
