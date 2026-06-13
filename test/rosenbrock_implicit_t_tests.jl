using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI
import StaticArrays as SA
import StaticArrays: SVector

# `RosenbrockAverageImplicitT` substeps the raw instantaneous pointwise tendency
# AUGMENTED with the local temperature `T` as a final state component, so the
# species+temperature Jacobian carries the latent-heating row and the
# Clausius-Clapeyron column `∂f/∂T`. The species-only `RosenbrockAverage`
# advances `T` EXPLICITLY between substeps, leaving that feedback operator-split;
# at coarse substeps the split loop rings about the saturation limit. These tests
# check: (1) the augmented tendency integrated explicitly reproduces the
# explicit-T `RosenbrockAverage` trajectory (so the latent-heat bookkeeping is
# the same physics, expressed as a rate); (2) implicit-T converges to a fine
# augmented-explicit reference under substep refinement; (3) implicit-T removes
# the operator-split saturation ringing where the explicit-T mode rings.

# Count sign flips in a sequence; a single monotone crossing is 1 flip, a true
# crossing-and-return (ringing) shows as >1.
function _count_flips(seq)
    s = sign.(seq)
    nz = filter(!iszero, s)
    isempty(nz) && return 0
    return count(i -> nz[i] != nz[i - 1], 2:length(nz))
end

#####
##### 2M+P3 implicit-T
#####

# Fine augmented-explicit reference: the SAME augmented RHS the implicit-T mode
# differentiates, integrated forward-Euler to convergence. This is the
# self-consistent accuracy reference for the implicit-T solve.
function _aug_explicit_2m(g, x0, T, Δt, nsub)
    FT = eltype(x0)
    h = Δt / FT(nsub)
    x = BMT.MicroState2MP3T{FT}(x0..., T)
    for _ in 1:nsub
        f = g(x)
        x = max.(x .+ h .* f, zero(FT))
    end
    return SVector{8, FT}(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]), x[9]
end

# Explicit-T `RosenbrockAverage`-style reference trajectory (forward-Euler
# species with the explicit between-substep T update of the species-only entry).
function _explicit_T_traj_2m(mp, tps, ρ, T, q_tot, x0, logλ, Δt, nsub)
    FT = eltype(x0)
    h = Δt / FT(nsub)
    cp_d = TDI.TD.Parameters.cp_d(tps)
    x = SVector{8, FT}(x0...)
    Tsub = T
    for _ in 1:nsub
        g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ)
        f = g(x)
        xp = x
        x = max.(x .+ h .* f, zero(FT))
        Δ = x - xp
        Ts = max(FT(150), Tsub)
        Tsub += (TDI.Lᵥ(tps, Ts) * (Δ[1] + Δ[3]) + TDI.Lₛ(tps, Ts) * Δ[5]) / cp_d
    end
    return x, Tsub
end

function test_implicit_t_2m(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    consistent_logλ(ρ, x) = P3.get_distribution_logλ(
        P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]),
    )

    function impl_step(x0, ρ, T, q_tot, logλ, Δt, nsub)
        t = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverageImplicitT(), BMT.Microphysics2Moment(), mp, tps,
            ρ, T, q_tot, x0..., logλ, Δt, nsub,
        )
        return SVector{8, FT}(x0...) .+ Δt .* SVector(values(t)...)
    end

    # x = [q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim]
    regimes = (
        (; name = "warm rain", ρ = FT(1.05), T = FT(288), q_tot = FT(0.015),
            x = FT[4e-4, 8e7, 2.1e-3, 5e4, 0, 0, 0, 0], logλ = FT(-Inf), tol = 0.01),
        (; name = "mixed phase", ρ = FT(0.78), T = FT(273.5), q_tot = FT(0.009),
            x = FT[2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8], logλ = nothing, tol = 0.25),
        (; name = "ice sublimation", ρ = FT(0.45), T = FT(253), q_tot = FT(4e-4),
            x = FT[0, 0, 0, 0, 8e-4, 5e5, 5e-4, 9e-7], logλ = nothing, tol = 0.02),
    )
    floors = SVector{8, FT}(1e-9, 1e-1, 1e-9, 1e-1, 1e-9, 1e-1, 1e-8, 1e-6)
    err_metric(x, x_ref, x0) = maximum(abs.(x .- x_ref) ./ (abs.(x0) .+ abs.(x_ref) .+ floors))

    @testset "2M augmented-explicit reproduces the explicit-T trajectory ($FT)" begin
        # The augmented tendency integrated forward-Euler must match the
        # species-only explicit-T `RosenbrockAverage`-style trajectory wherever
        # the positivity floor does not clamp the increment: the dT/dt is the
        # rate form of the same latent-heat update. At a small step (Δt = 1,
        # refined nsub) the floor does not bite and they agree to roundoff
        # (with floor events the unfloored augmented `T` legitimately diverges
        # from the floored-increment explicit update, the only difference).
        Δt = FT(1)
        for r in regimes
            logλ = isnothing(r.logλ) ? consistent_logλ(r.ρ, r.x) : r.logλ
            g = BMT.Instantaneous2MP3TTendency(mp, tps, r.ρ, r.q_tot, logλ)
            xa, Ta = _explicit_T_traj_2m(mp, tps, r.ρ, r.T, r.q_tot, r.x, logλ, Δt, 64)
            xb, Tb = _aug_explicit_2m(g, r.x, r.T, Δt, 64)
            x0 = SVector{8, FT}(r.x...)
            @test err_metric(xb, xa, x0) < sqrt(eps(FT))
            @test abs(Ta - Tb) < sqrt(eps(FT)) * abs(Ta)
        end
    end

    @testset "2M implicit-T converges to the fine augmented-explicit reference ($FT)" begin
        Δt = FT(10)
        for r in regimes
            logλ = isnothing(r.logλ) ? consistent_logλ(r.ρ, r.x) : r.logλ
            g = BMT.Instantaneous2MP3TTendency(mp, tps, r.ρ, r.q_tot, logλ)
            x_ref, _ = _aug_explicit_2m(g, r.x, r.T, Δt, 2048)
            x0 = SVector{8, FT}(r.x...)
            err(x) = err_metric(x, x_ref, x0)
            errs = [err(impl_step(r.x, r.ρ, r.T, r.q_tot, logλ, Δt, n)) for n in (1, 4, 16)]
            @test all(isfinite, errs)
            # accuracy improves under substep refinement...
            @test errs[3] ≤ max(errs[1], FT(1e-3)) * (1 + sqrt(eps(FT)))
            # ...to within a regime-calibrated distance of the reference
            @test errs[3] < (FT == Float64 ? r.tol : 2 * r.tol)
        end
    end

    @testset "2M implicit-T removes the operator-split saturation ringing ($FT)" begin
        # band_spinup warm spin-up: the documented in-spec coarse-substep case
        # where the explicit-T mode rings the liquid saturation (issue-008
        # ringing probe: 3 decaying s_liq flips at nsub 4). q_vap convention
        # q_tot - q_lcl - q_rai - q_ice; s_liq = q_vap - q_sat_liq.
        x_band = FT[1e-13, 1e2, 0, 0, 0, 0, 0, 0]
        ρ = FT(1.0);
        T = FT(288);
        q_tot = FT(0.02);
        Δt = FT(60);
        nsub = 4
        qvap(x) = max(zero(FT), q_tot - x[1] - x[3] - x[5])
        sliq(Tk, x) = qvap(x) - TDI.saturation_vapor_specific_content_over_liquid(tps, Tk, ρ)

        h = Δt / FT(nsub)
        # explicit-T trajectory (species-only update + between-substep T)
        function flips_explicit()
            x = BMT.MicroState2MP3{FT}(x_band...);
            Tsub = T;
            cp_d = TDI.TD.Parameters.cp_d(tps)
            seq = FT[]
            for _ in 1:nsub
                push!(seq, sliq(Tsub, x))
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, Tsub, q_tot, FT(-Inf))
                f = g(x);
                xp = x
                J = BMT.FD.jacobian(g, x);
                z = BMT._rosenbrock_channel_mask(x)
                x = all(isfinite, J) ? BMT._rosenbrock_update(x, f, J, z, h) : max.(x .+ h .* f, 0)
                Δ = x - xp;
                Ts = max(FT(150), Tsub)
                Tsub += (TDI.Lᵥ(tps, Ts) * (Δ[1] + Δ[3]) + TDI.Lₛ(tps, Ts) * Δ[5]) / cp_d
            end
            push!(seq, sliq(Tsub, x));
            return _count_flips(seq)
        end
        # implicit-T trajectory (T inside the solve)
        function flips_implicit()
            g = BMT.Instantaneous2MP3TTendency(mp, tps, ρ, q_tot, FT(-Inf))
            x = BMT.MicroState2MP3T{FT}(x_band..., T);
            seq = FT[]
            for _ in 1:nsub
                push!(seq, sliq(x[9], x))
                f = g(x);
                J = BMT.FD.jacobian(g, x);
                z = BMT._rosenbrock_channel_mask(x)
                x = all(isfinite, J) ? BMT._rosenbrock_update(x, f, J, z, h) : max.(x .+ h .* f, 0)
            end
            push!(seq, sliq(x[9], x));
            return _count_flips(seq)
        end
        ne = flips_explicit()
        ni = flips_implicit()
        # explicit-T rings (multiple flips); implicit-T is monotone (≤ 1 crossing)
        @test ne ≥ 2
        @test ni ≤ 1
    end

    @testset "2M implicit-T degenerate and finite/non-negative ($FT)" begin
        # all-zero state: degenerate gate -> explicit substeps -> exactly zero
        t0 = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverageImplicitT(), BMT.Microphysics2Moment(), mp, tps,
            FT(1), FT(273), FT(0),
            FT(0), FT(0), FT(0), FT(0), FT(0), FT(0), FT(0), FT(0),
            FT(-Inf), FT(60), 4,
        )
        @test all(iszero, values(t0))
        # nsub defaults to 1 and accepts the trailing-argument form
        r = regimes[2]
        logλ = consistent_logλ(r.ρ, r.x)
        t1 = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverageImplicitT(), BMT.Microphysics2Moment(), mp, tps,
            r.ρ, r.T, r.q_tot, r.x..., logλ, FT(60),
        )
        @test all(isfinite, values(t1))
        # hostile stress state stays finite and non-negative
        x_stress = FT[1e-6, 1e6, 1e-12, 1e-2, 8e-4, 5e5, 5e-4, 9e-7]
        logλ_s = consistent_logλ(FT(0.45), x_stress)
        for nsub in (1, 2, 8), Δt in (FT(60), FT(300))
            t = BMT.bulk_microphysics_tendencies(
                BMT.RosenbrockAverageImplicitT(), BMT.Microphysics2Moment(), mp, tps,
                FT(0.45), FT(233), FT(0.003), x_stress..., logλ_s, Δt, nsub,
            )
            x1 = SVector{8, FT}(x_stress...) .+ Δt .* SVector(values(t)...)
            @test all(isfinite, x1)
            tol = eps(FT) .* (abs.(SVector{8, FT}(x_stress...)) .+ Δt .* abs.(SVector(values(t)...)))
            @test all(x1 .>= -tol)
        end
    end
end

test_implicit_t_2m(Float64)
test_implicit_t_2m(Float32)

#####
##### 1M implicit-T
#####

function _aug_explicit_1m(g, x0, T, Δt, nsub)
    FT = eltype(x0)
    h = Δt / FT(nsub)
    x = BMT.MicroState1MT{FT}(x0..., T)
    for _ in 1:nsub
        f = g(x)
        x = max.(x .+ h .* f, zero(FT))
    end
    return SVector{4, FT}(x[1], x[2], x[3], x[4]), x[5]
end

function _explicit_T_traj_1m(mp, tps, ρ, T, q_tot, x0, Δt, nsub)
    FT = eltype(x0)
    h = Δt / FT(nsub)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)
    x = SVector{4, FT}(x0...)
    Tsub = T
    for _ in 1:nsub
        g = BMT.Raw1MTendency(mp, tps, ρ, Tsub, q_tot)
        f = g(x)
        xp = x
        x = max.(x .+ h .* f, zero(FT))
        Δ = x - xp
        Tsub += Lv_over_cp * (Δ[1] + Δ[3]) + Ls_over_cp * (Δ[2] + Δ[4])
    end
    return x, Tsub
end

function test_implicit_t_1m(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics1MParams(FT)
    T_frz = TDI.T_freeze(tps)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)

    function impl_step(x0, ρ, T, q_tot, Δt, nsub)
        t = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverageImplicitT(), BMT.Microphysics1Moment(), mp, tps,
            ρ, T, q_tot, x0..., Δt, nsub,
        )
        return SVector{4, FT}(x0...) .+ Δt .* SVector(values(t)...)
    end

    regimes = (
        (; name = "warm rain", ρ = FT(1.0), T = T_frz + FT(17), q_tot = FT(0.018),
            x = FT[2e-3, 0, 5e-4, 0], tol = 0.01),
        (; name = "mixed warm", ρ = FT(1.2), T = T_frz + FT(5), q_tot = FT(0.012),
            x = FT[5e-4, 2e-4, 3e-4, 3e-4], tol = 0.05),
        (; name = "mixed cold", ρ = FT(1.2), T = T_frz - FT(10), q_tot = FT(0.012),
            x = FT[3e-4, 5e-4, 2e-4, 4e-4], tol = 0.1),
        (; name = "snow cold depo", ρ = FT(1.2), T = T_frz - FT(15), q_tot = FT(0.008),
            x = FT[0, 0, 0, 1e-4], tol = 0.1),
    )
    floor = FT(1e-9)
    err_metric(x, x_ref, x0) = maximum(abs.(x .- x_ref) ./ (abs.(x0) .+ abs.(x_ref) .+ floor))

    @testset "1M augmented-explicit reproduces the explicit-T trajectory ($FT)" begin
        # The 1M warm-rain dynamics are mild enough that even at Δt = 20 the
        # floor does not bite, so the augmented-explicit trajectory matches the
        # explicit-T `RosenbrockAverage`-style trajectory to roundoff.
        Δt = FT(20)
        for r in regimes
            g = BMT.Raw1MTTendency(mp, tps, r.ρ, r.q_tot, Lv_over_cp, Ls_over_cp)
            xa, Ta = _explicit_T_traj_1m(mp, tps, r.ρ, r.T, r.q_tot, r.x, Δt, 64)
            xb, Tb = _aug_explicit_1m(g, r.x, r.T, Δt, 64)
            x0 = SVector{4, FT}(r.x...)
            @test err_metric(xb, xa, x0) < sqrt(eps(FT))
            @test abs(Ta - Tb) < sqrt(eps(FT)) * abs(Ta)
        end
    end

    @testset "1M implicit-T converges to the fine augmented-explicit reference ($FT)" begin
        Δt = FT(20)
        for r in regimes
            g = BMT.Raw1MTTendency(mp, tps, r.ρ, r.q_tot, Lv_over_cp, Ls_over_cp)
            x_ref, _ = _aug_explicit_1m(g, r.x, r.T, Δt, 4096)
            x0 = SVector{4, FT}(r.x...)
            err(x) = err_metric(x, x_ref, x0)
            errs = [err(impl_step(r.x, r.ρ, r.T, r.q_tot, Δt, n)) for n in (1, 4, 16)]
            @test all(isfinite, errs)
            @test errs[3] ≤ max(errs[1], FT(1e-3)) * (1 + sqrt(eps(FT)))
            @test errs[3] < (FT == Float64 ? r.tol : 2 * r.tol)
        end
    end

    @testset "1M implicit-T removes the operator-split saturation ringing ($FT)" begin
        # Warm, strongly supersaturated spin-up at a coarse substep: the
        # explicit-T mode rings the liquid saturation (the latent-heat T
        # overshoot reverses the supersaturation sign); implicit-T is monotone.
        x0 = FT[1e-6, 0, 1e-5, 0]
        ρ = FT(1.0);
        T = T_frz + FT(8);
        q_tot = FT(0.025);
        Δt = FT(120);
        nsub = 4
        qvap(x) = max(zero(FT), q_tot - x[1] - x[2] - x[3] - x[4])
        sliq(Tk, x) = qvap(x) - TDI.saturation_vapor_specific_content_over_liquid(tps, Tk, ρ)
        h = Δt / FT(nsub)
        function flips_explicit()
            x = BMT.MicroState1M{FT}(x0...);
            Tsub = T;
            seq = FT[]
            for _ in 1:nsub
                push!(seq, sliq(Tsub, x))
                g = BMT.Raw1MTendency(mp, tps, ρ, Tsub, q_tot)
                f = g(x);
                xp = x
                J = BMT.FD.jacobian(g, x);
                z = BMT._rosenbrock_channel_mask(x)
                x = all(isfinite, J) ? BMT._rosenbrock_update(x, f, J, z, h) : max.(x .+ h .* f, 0)
                Δ = x - xp
                Tsub += Lv_over_cp * (Δ[1] + Δ[3]) + Ls_over_cp * (Δ[2] + Δ[4])
            end
            push!(seq, sliq(Tsub, x));
            return _count_flips(seq)
        end
        function flips_implicit()
            g = BMT.Raw1MTTendency(mp, tps, ρ, q_tot, Lv_over_cp, Ls_over_cp)
            x = BMT.MicroState1MT{FT}(x0..., T);
            seq = FT[]
            for _ in 1:nsub
                push!(seq, sliq(x[5], x))
                f = g(x);
                J = BMT.FD.jacobian(g, x);
                z = BMT._rosenbrock_channel_mask(x)
                x = all(isfinite, J) ? BMT._rosenbrock_update(x, f, J, z, h) : max.(x .+ h .* f, 0)
            end
            push!(seq, sliq(x[5], x));
            return _count_flips(seq)
        end
        @test flips_explicit() ≥ 2
        @test flips_implicit() ≤ 1
    end

    @testset "1M implicit-T degenerate and finite/non-negative ($FT)" begin
        t0 = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverageImplicitT(), BMT.Microphysics1Moment(), mp, tps,
            FT(1), FT(273), FT(0), FT(0), FT(0), FT(0), FT(0), FT(60), 4,
        )
        @test all(iszero, values(t0))
        t1 = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverageImplicitT(), BMT.Microphysics1Moment(), mp, tps,
            FT(1.2), FT(278), FT(0.012), FT(5e-4), FT(2e-4), FT(3e-4), FT(3e-4), FT(60),
        )
        @test all(isfinite, values(t1))
        x_stress = FT[2e-3, 1e-3, 2e-3, 3e-3]
        for nsub in (1, 2, 8), Δt in (FT(60), FT(300))
            for (ρ, T, q_tot) in (
                (FT(1.2), T_frz + FT(8), FT(0.02)),
                (FT(0.6), T_frz - FT(12), FT(0.005)),
            )
                t = BMT.bulk_microphysics_tendencies(
                    BMT.RosenbrockAverageImplicitT(), BMT.Microphysics1Moment(), mp, tps,
                    ρ, T, q_tot, x_stress..., Δt, nsub,
                )
                x1 = SVector{4, FT}(x_stress...) .+ Δt .* SVector(values(t)...)
                @test all(isfinite, x1)
                tol = eps(FT) .* (abs.(SVector{4, FT}(x_stress...)) .+ Δt .* abs.(SVector(values(t)...)))
                @test all(x1 .>= -tol)
            end
        end
    end
end

test_implicit_t_1m(Float64)
test_implicit_t_1m(Float32)

# Allocation + JET checks on the implicit-T hot call (compiler-version sensitive,
# like the other perf assertions; see performance_tests.jl, rosenbrock_mode_tests.jl).
if VERSION >= v"1.12"
    import JET
    @testset "2M implicit-T allocations and inference" begin
        FT = Float64
        tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
        mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
        p3 = mp.ice.scheme
        st = P3.state_from_prognostic(
            p3, FT(0.78) * FT(1e-4), FT(0.78) * FT(2e5),
            FT(0.78) * FT(4e-5), FT(0.78) * FT(6e-8),
        )
        logλ = P3.get_distribution_logλ(st)
        call() = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverageImplicitT(), BMT.Microphysics2Moment(), mp, tps,
            FT(0.78), FT(273.5), FT(0.009),
            FT(2e-4), FT(5e7), FT(1e-4), FT(4e4), FT(1e-4), FT(2e5), FT(4e-5), FT(6e-8),
            logλ, FT(60), 4,
        )
        call()
        @test (@allocated call()) == 0
    end

    @testset "1M implicit-T allocations and inference" begin
        for FT in (Float64, Float32)
            tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
            mp = CMP.Microphysics1MParams(FT)
            call() = BMT.bulk_microphysics_tendencies(
                BMT.RosenbrockAverageImplicitT(), BMT.Microphysics1Moment(), mp, tps,
                FT(1.2), FT(278), FT(0.012),
                FT(5e-4), FT(2e-4), FT(3e-4), FT(3e-4), FT(20), 4,
            )
            call()
            @test (@allocated call()) == 0
            rep = JET.report_call(
                BMT.bulk_microphysics_tendencies,
                typeof.((
                    BMT.RosenbrockAverageImplicitT(), BMT.Microphysics1Moment(), mp, tps,
                    FT(1.2), FT(278), FT(0.012),
                    FT(5e-4), FT(2e-4), FT(3e-4), FT(3e-4), FT(20), 4,
                )),
            )
            @test isempty(JET.get_reports(rep))
        end
    end
end
