using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.ThermodynamicsInterface as TDI
import StaticArrays: SVector

# The 1M `RosenbrockAverage` substeps the raw instantaneous pointwise 1M
# tendency with a linearized-implicit (Rosenbrock-Euler) update, linearizing
# with the exact `ForwardDiff` Jacobian — a separate option from the hand-built
# `LinearizedAverage` (donor-cell-modified system matrix). The accuracy
# reference is a finely-resolved forward-Euler integration of the same raw
# tendency, with the identical constant-latent-heat T update and frozen q_tot.

# Fine explicit reference: integrate the raw 1M tendency with `nsub`
# forward-Euler substeps, updating T from latent heating like both averaged
# schemes (constant L_v/L_s, liquid+rain on L_v, ice+snow on L_s).
function explicit_reference_1m(mp, tps, ρ, T, q_tot, x0, Δt, nsub)
    FT = typeof(q_tot)
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
        Tsub +=
            Lv_over_cp * ((x[1] - xp[1]) + (x[3] - xp[3])) +
            Ls_over_cp * ((x[2] - xp[2]) + (x[4] - xp[4]))
    end
    return x
end

function test_rosenbrock_1m_mode(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics1MParams(FT)

    function ros_step(x0, ρ, T, q_tot, Δt, nsub)
        t = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverage(), BMT.Microphysics1Moment(), mp, tps,
            ρ, T, q_tot, x0..., Δt, nsub,
        )
        return SVector{4, FT}(x0...) .+ Δt .* SVector(values(t)...)
    end
    function lin_step(x0, ρ, T, q_tot, Δt, nsub)
        t = BMT.bulk_microphysics_tendencies(
            BMT.LinearizedAverage(), BMT.Microphysics1Moment(), mp, tps,
            ρ, T, q_tot, x0..., Δt, nsub,
        )
        return SVector{4, FT}(x0...) .+ Δt .* SVector(values(t)...)
    end

    # x = [q_lcl, q_icl, q_rai, q_sno]. `tol` is the nsub=16 max-relative error
    # vs a 4096-substep reference with ~2x headroom over the measured value
    # (warm 1.3e-3, mixed-warm 5.9e-3, mixed-cold 4.0e-2, snow 3.9e-2; identical
    # F64/F32). The cold regimes plateau higher: the deposition/freezing arms
    # are near-non-smooth at the regime boundary, so the single-linearization
    # error does not vanish under refinement as cleanly as in the warm regimes.
    T_frz = TDI.T_freeze(tps)
    regimes = (
        (; name = "warm rain", ρ = FT(1.0), T = T_frz + FT(17), q_tot = FT(0.018),
            x = FT[2e-3, 0, 5e-4, 0], tol = 0.005),
        (; name = "mixed warm", ρ = FT(1.2), T = T_frz + FT(5), q_tot = FT(0.012),
            x = FT[5e-4, 2e-4, 3e-4, 3e-4], tol = 0.015),
        (; name = "mixed cold", ρ = FT(1.2), T = T_frz - FT(10), q_tot = FT(0.012),
            x = FT[3e-4, 5e-4, 2e-4, 4e-4], tol = 0.08),
        (; name = "snow cold depo", ρ = FT(1.2), T = T_frz - FT(15), q_tot = FT(0.008),
            x = FT[0, 0, 0, 1e-4], tol = 0.08),
    )

    # Per-species scale with a physical floor so a collapsing species (reference
    # depleted to ~0) cannot dominate a plain relative metric.
    floor = FT(1e-9)
    err_metric(x, x_ref, x0) =
        maximum(abs.(x .- x_ref) ./ (abs.(x0) .+ abs.(x_ref) .+ floor))

    @testset "1M RosenbrockAverage vs fine explicit reference ($FT)" begin
        Δt = FT(20)
        for r in regimes
            x_ref = explicit_reference_1m(mp, tps, r.ρ, r.T, r.q_tot, r.x, Δt, 4096)
            x0 = SVector{4, FT}(r.x...)
            err(x) = err_metric(x, x_ref, x0)
            errs = [err(ros_step(r.x, r.ρ, r.T, r.q_tot, Δt, n)) for n in (1, 4, 16)]
            @test all(isfinite, errs)
            # accuracy improves under substep refinement (slack for the one
            # stiff-channel coarse-nsub non-monotonicity documented in the
            # AD-vs-manual study)...
            @test errs[3] ≤ max(errs[1], FT(1e-3)) * (1 + sqrt(eps(FT)))
            # ...to within a regime-calibrated distance of the reference
            @test errs[3] < (FT == Float64 ? r.tol : 2 * r.tol)
        end
    end

    @testset "1M RosenbrockAverage vs LinearizedAverage ($FT)" begin
        # The two options solve different linearizations of the same raw
        # tendency; at refined substepping both track the same fine reference,
        # so their realized states must agree to the reference tolerance.
        Δt = FT(20)
        for r in regimes
            x_ref = explicit_reference_1m(mp, tps, r.ρ, r.T, r.q_tot, r.x, Δt, 4096)
            x0 = SVector{4, FT}(r.x...)
            x_ros = ros_step(r.x, r.ρ, r.T, r.q_tot, Δt, 16)
            x_lin = lin_step(r.x, r.ρ, r.T, r.q_tot, Δt, 16)
            @test all(isfinite, x_lin)
            # both within the regime tolerance of the shared reference
            @test err_metric(x_ros, x_ref, x0) < (FT == Float64 ? r.tol : 2 * r.tol)
            @test err_metric(x_lin, x_ref, x0) < (FT == Float64 ? 2 * r.tol : 4 * r.tol)
        end
    end

    @testset "1M degenerate and trivial states ($FT)" begin
        # all-zero state: degenerate channel gate -> explicit substeps -> zero
        t0 = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverage(), BMT.Microphysics1Moment(), mp, tps,
            FT(1), FT(273), FT(0), FT(0), FT(0), FT(0), FT(0), FT(60), 4,
        )
        @test all(iszero, values(t0))
        # nsub defaults to 1 and accepts the trailing-argument form
        t1 = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverage(), BMT.Microphysics1Moment(), mp, tps,
            FT(1.2), FT(278), FT(0.012), FT(5e-4), FT(2e-4), FT(3e-4), FT(3e-4), FT(60),
        )
        @test all(isfinite, values(t1))
    end

    @testset "1M substeps stay finite and non-negative ($FT)" begin
        # hostile stress: strongly out-of-equilibrium all-species state, large Δt
        x_stress = FT[2e-3, 1e-3, 2e-3, 3e-3]
        for nsub in (1, 2, 8), Δt in (FT(60), FT(300))
            for (ρ, T, q_tot) in (
                (FT(1.2), T_frz + FT(8), FT(0.02)),
                (FT(0.6), T_frz - FT(12), FT(0.005)),
            )
                t = BMT.bulk_microphysics_tendencies(
                    BMT.RosenbrockAverage(), BMT.Microphysics1Moment(), mp, tps,
                    ρ, T, q_tot, x_stress..., Δt, nsub,
                )
                x1 = SVector{4, FT}(x_stress...) .+ Δt .* SVector(values(t)...)
                @test all(isfinite, x1)
                # non-negative up to the roundoff of the host-side
                # x + Δt * t reconstruction (internally floored at zero)
                tol =
                    eps(FT) .*
                    (abs.(SVector{4, FT}(x_stress...)) .+ Δt .* abs.(SVector(values(t)...)))
                @test all(x1 .>= -tol)
            end
        end
    end

    @testset "1M near-empty channels take the explicit path ($FT)" begin
        # condensed masses in (eps, 1e-10) produce finite but enormous Jacobian
        # rows; the channel mask must route them to forward Euler so the result
        # tracks the explicit reference rather than fabricating phantom mass.
        x_band = FT[1e-13, 0, 1e-3, 0]
        Δt = FT(60)
        ρ = FT(1.0)
        T = T_frz + FT(10)
        q_tot = FT(0.02)
        x_ref = explicit_reference_1m(mp, tps, ρ, T, q_tot, x_band, Δt, 4096)
        x16 = ros_step(x_band, ρ, T, q_tot, Δt, 16)
        @test all(isfinite, x16)
        # band channel stays bounded near the reference (no phantom growth)
        @test x16[1] < max(FT(10) * x_ref[1], FT(1e-9))
    end
end

test_rosenbrock_1m_mode(Float64)
test_rosenbrock_1m_mode(Float32)

# Allocation + JET checks on the hot call (compiler-version sensitive, like the
# other perf assertions; see performance_tests.jl and rosenbrock_mode_tests.jl).
if VERSION >= v"1.12"
    import JET
    @testset "1M RosenbrockAverage allocations and inference" begin
        for FT in (Float64, Float32)
            tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
            mp = CMP.Microphysics1MParams(FT)
            call() = BMT.bulk_microphysics_tendencies(
                BMT.RosenbrockAverage(), BMT.Microphysics1Moment(), mp, tps,
                FT(1.2), FT(278), FT(0.012),
                FT(5e-4), FT(2e-4), FT(3e-4), FT(3e-4), FT(20), 4,
            )
            call()
            @test (@allocated call()) == 0
            rep = JET.report_call(
                BMT.bulk_microphysics_tendencies,
                typeof.((
                    BMT.RosenbrockAverage(), BMT.Microphysics1Moment(), mp, tps,
                    FT(1.2), FT(278), FT(0.012),
                    FT(5e-4), FT(2e-4), FT(3e-4), FT(3e-4), FT(20), 4,
                )),
            )
            @test isempty(JET.get_reports(rep))
        end
    end
end
