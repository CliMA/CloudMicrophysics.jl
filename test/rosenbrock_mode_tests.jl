using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.HetIceNucleation as CM_HetIce
import CloudMicrophysics.HomIceNucleation as CM_HomIce
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Common as CO
import CloudMicrophysics.DistributionTools as DT
import CloudMicrophysics.MicrophysicsNonEq as CMNonEq
import CloudMicrophysics.Utilities as UT
import StaticArrays: SVector
import ForwardDiff as FD

# The process-function convention bundles the microphysical state and the thermodynamic state
# into two contexts, which is what `_substep_context` builds in production. These assemble the
# same two from the flat argument order these tests are written around, so a test exercises the
# production entry point rather than a signature kept alive for its benefit.
#
# `w` and `p` are zero: they exist on the production `thermo` for the activation slot, and no
# function these helpers are used with reads either.
# The eight tendencies the host applies, selected by name from the march's fixed-shape carrier.
# The carrier also holds `dn_lcl_activation_dt`, whose number source is already inside
# `dn_lcl_dt`, and `extras`, which is a `NamedTuple` and not a number at all.
const _APPLIED_2MP3 = (
    :dq_lcl_dt, :dn_lcl_dt, :dq_rai_dt, :dn_rai_dt,
    :dq_ice_dt, :dn_ice_dt, :dq_rim_dt, :db_rim_dt,
)
_applied_2mp3(t) = map(k -> getproperty(t, k), _APPLIED_2MP3)

_micro_2mp3(q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim) =
    (; q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
_thermo_2mp3(ρ, T, logλ) = (; ρ, T, w = zero(ρ), p = zero(ρ), logλ)

# PORTED FROM the campaign's `test/rosenbrock_mode_tests.jl` (`campaign/he/p3-density-floors`
# in this project's CM clone; ~5400 lines, ~949 assertions, the campaign's largest single
# test asset). See the manifest/report for the file-by-file account of what changed and why;
# the mechanical repairs applied throughout this file are:
#
#   - `BMT._per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot, q_lcl, n_lcl, ..., logλ)` and
#     `BMT._per_process_2mp3(...)` (a long positional argument list) were renamed to
#     `p3_2m_process_rates(mp, tps, micro, thermo)`, which takes two context `NamedTuple`s
#     instead. Rather than hand-rewriting each of the ~50 call sites below into the new
#     argument shape, `_per_process_2mp3_and_riming`/`_per_process_2mp3` are defined locally
#     just below as TEST SCAFFOLDING - thin wrappers, not part of any public API, that build
#     those two `NamedTuple`s and call `BMT.p3_2m_process_rates`, reproducing the campaign's
#     own calling convention exactly. Every call site below therefore only needed `BMT.`
#     stripped from the name, not its argument list rewritten; a reader should not mistake
#     either wrapper for a live `BulkMicrophysicsTendencies` function (see their own
#     docstrings, just below the imports).
#   - Two `pp` slot names changed with the nucleation baseline: `bigg_immersion` ->
#     `immersion_freezing`, `f23_deposition` -> `ice_deposition`. The one test built
#     specifically around the F23 slot name (`test_f23_shift_has_one_definition`) is
#     reframed below: `mp.ice.ice_nucleation` now defaults to `ExponentialSupercoolingINP`
#     (Cooper), not `Frostenberg2023`, so what that test actually exercises at its chosen
#     state is the DEFAULT closure's deposition-nucleation primal consistency, not F23
#     specifically, despite its inherited name - the primal-level invariant it checks
#     (`sum(values(pp)) == g(x)`, and the slot fires) holds regardless of which target
#     spectrum is active, so the test is kept and renamed rather than dropped.
#   - `bulk_microphysics_tendencies`'s Rosenbrock-mode return is a FIXED-SHAPE carrier that has
#     grown twice: a trailing `extras` `NamedTuple` the campaign's carrier did not have, and
#     then `dn_lcl_activation_dt` beside it. Selecting positionally does not survive that.
#     `Base.front` drops ONE trailing field, so once there were two it handed nine values to an
#     eight-element `SVector` and every such site raised a `DimensionMismatch`. The eight
#     APPLIED tendencies are therefore selected BY NAME through `_applied_2mp3` below, which is
#     what the host does in its own `_select_2mp3_tendency`, so these tests no longer care how
#     many further fields the carrier grows.
#   - `rosenbrock_donor()` does not appear in this file (verified by search) so the D4
#     preset-deletion ruling needed no repair here, unlike in rosenbrock_framework_tests.jl.
#   - `MicroState2MP3T`, `Temperature2MP3Tendency`, `TemperatureCoupledJacobian` and
#     `Verbose2MP3Tendency` are ALL current names already, used correctly throughout by the
#     campaign version (this campaign file postdates the T-as-state / 9x9 work), so the
#     temperature-coupled and verbose-primal call sites below needed no repair.
#
# THE τ_act FIELD-ACCESS BUG IS FIXED, as of the tip this file was last checked against.
# `_jacobian_2mp3_manual` used to read `mp.ice.inp_depletion_model.τ_act` off the default
# `NIceProxyDepletion()`, a fieldless marker struct with no such field, which threw
# unconditionally whenever `ManualJacobian`/`TemperatureCoupledJacobian` computed a
# Jacobian - most of this file's testsets, since it is fundamentally a manual-Jacobian
# entry-by-entry validation suite. It now evaluates `CM_HetIce.delivery_rate` (the same
# state-dependent rate the primal uses) instead of reading a stored constant. The branch is
# being restacked continuously as defects surface, so treat this as a snapshot, not a
# standing guarantee - re-check `src/BMT_2mp3_jacobian.jl`'s deposition-nucleation entry
# before relying on it.
#
# A SEPARATE, DEEPER issue surfaced while fixing two of this file's OWN direct reads of the
# same now-removed field (`test_timescale_number_couplings`, one `J[6,6]` expected-value
# formula; `test_freezing_is_differentiable_with_no_pathway_open`, both fixed in place using
# `CM_HetIce.delivery_rate` the same way the Jacobian now does): `cloud_freezing_rate`
# itself no longer carries ANY ice-nucleating-particle budget or cap on the heterogeneous
# coefficient - "the immersion coefficient is Bigg alone" (its own docstring) - which is a
# further design step past the nucleation-baseline ruling this file's other repairs were
# scoped to. That broke the PREMISE, not just the call signature, of two whole test
# functions built around the old capped-Bigg design (`J_cap`, `J_het_uncapped`,
# `immersion_limit_rate`/`τ_act`/`n_active` keywords `cloud_freezing_rate` no longer
# accepts): `test_cloud_freezing_is_uniform_with_rain` and
# `test_freezing_is_differentiable_with_no_pathway_open`, both REPAIRED with several
# testsets DROPPED outright (the ones that tested the cap's own existence/scaling/behavior)
# and the rest adapted or left alone - see each function's own header comment for the exact
# account. This is worth naming as its own category of risk for the rest of this large
# file: a source design can move further than the one rename/signature-change this port was
# scoped to catch, and the only thing that reliably catches it is reading each testset's
# actual claims against current source, which was not done exhaustively across all ~40
# functions here (see the report to the team lead).

"""
    _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ, w = zero(ρ), p = zero(ρ))

TEST SCAFFOLDING, NOT PART OF THE PUBLIC API: this file's own local reproduction of the
campaign's long-positional calling convention over
[`p3_2m_process_rates`](@ref CloudMicrophysics.BulkMicrophysicsTendencies.p3_2m_process_rates),
which takes two context `NamedTuple`s (`micro`, `thermo`) instead. Defined here purely so
the ported call sites below, all written against the old convention, needed no
argument-shape rewriting - only `BMT.` stripped from the name. A reader landing on a call
site below should not read the name as a live `BulkMicrophysicsTendencies` function; the
one function this file actually calls is `BMT.p3_2m_process_rates`, right here.
"""
function _per_process_2mp3_and_riming(
    mp, tps, ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    w = zero(ρ), p = zero(ρ),
)
    micro = (; q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    thermo = (; ρ, T, w, p, logλ)
    return BMT.p3_2m_process_rates(mp, tps, micro, thermo)
end

"""
    _per_process_2mp3(mp, tps, ρ, T, q_tot, ...)

TEST SCAFFOLDING, NOT PART OF THE PUBLIC API (see [`_per_process_2mp3_and_riming`](@ref)):
the campaign's plain (no riming-split) form, the per-process breakdown alone.
"""
_per_process_2mp3(args...) = first(_per_process_2mp3_and_riming(args...))


# `RosenbrockAverage` substeps the raw instantaneous pointwise 2M+P3 tendency
# with a linearized-implicit (Rosenbrock-Euler) update. The reference for
# accuracy tests is a finely-resolved forward-Euler integration of the same
# raw tendency (identical T update and frozen logλ/q_tot semantics), so these
# tests use the unlimited exact configuration (no increment limiter) that
# converges to that reference under substep refinement.

function explicit_reference(mp, tps, ρ, T, q_tot, x0, logλ, Δt, nsub)
    FT = typeof(q_tot)
    h = Δt / FT(nsub)
    x = SVector{8, FT}(x0...)
    Tsub = T
    cp_d = TDI.TD.Parameters.cp_d(tps)
    for _ in 1:nsub
        g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ)
        f = g(x)
        xp = x
        x = max.(x .+ h .* f, zero(FT))
        Ts = max(FT(150), Tsub)
        Tsub +=
            (
                TDI.Lᵥ(tps, Ts) * ((x[1] - xp[1]) + (x[3] - xp[3])) +
                TDI.Lₛ(tps, Ts) * (x[5] - xp[5])
            ) / cp_d
    end
    return x
end

function test_rosenbrock_mode(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme

    # exact Jacobian, no increment limiter: the configuration that converges to
    # the unlimited forward-Euler reference under substep refinement
    mode = BMT.RosenbrockAverage(
        jacobian = BMT.ExactJacobian(),
        growth = BMT.ImplicitGrowth(),
        limiter = BMT.NoLimiter(),
    )

    function consistent_logλ(ρ, x)
        st = P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8])
        return P3.get_distribution_logλ(st)
    end
    function step(x0, ρ, T, q_tot, logλ, Δt, nsub)
        t = BMT.bulk_microphysics_tendencies(
            mode, BMT.Microphysics2Moment(), mp, tps,
            ρ, T, q_tot, x0..., logλ, Δt, nsub,
        )
        return SVector{8, FT}(x0...) .+ Δt .* SVector{8, FT}(_applied_2mp3(t)...)
    end

    # x = [q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim]
    regimes = (
        (; ρ = FT(1.05), T = FT(288), q_tot = FT(0.015),  # warm rain
            x = FT[4e-4, 8e7, 2.1e-3, 5e4, 0, 0, 0, 0], logλ = FT(-Inf), tol = 0.002),
        (; ρ = FT(0.78), T = FT(273.5), q_tot = FT(0.009),  # mixed phase
            x = FT[2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8], logλ = nothing, tol = 0.07),
        (; ρ = FT(0.45), T = FT(253), q_tot = FT(4e-4),  # ice sublimation
            x = FT[0, 0, 0, 0, 8e-4, 5e5, 5e-4, 9e-7], logλ = nothing, tol = 0.012),
    )

    # Per-species scales with physical floors so a collapsing species cannot
    # dominate a plain relative metric. The rime species (q_rim, b_rim) carry
    # larger floors: once ice has melted away, single precision leaves an
    # O(1e-7) coupled volume/mass residual.
    floors = SVector{8, FT}(1e-9, 1e-1, 1e-9, 1e-1, 1e-9, 1e-1, 1e-8, 1e-6)
    err_metric(x, x_ref, x0) = maximum(abs.(x .- x_ref) ./ (abs.(x0) .+ abs.(x_ref) .+ floors))

    @testset "RosenbrockAverage vs fine explicit reference ($FT)" begin
        Δt = FT(10)
        for r in regimes
            logλ = isnothing(r.logλ) ? consistent_logλ(r.ρ, r.x) : r.logλ
            x_ref = explicit_reference(mp, tps, r.ρ, r.T, r.q_tot, r.x, logλ, Δt, 2048)
            x0 = SVector{8, FT}(r.x...)
            err(x) = err_metric(x, x_ref, x0)
            errs = [err(step(r.x, r.ρ, r.T, r.q_tot, logλ, Δt, n)) for n in (1, 4, 16)]
            @test all(isfinite, errs)
            # accuracy improves under substep refinement...
            @test errs[3] ≤ max(errs[1], FT(1e-3)) * (1 + sqrt(eps(FT)))
            # ...to within a regime-calibrated distance of the reference
            @test errs[3] < (FT == Float64 ? r.tol : 2 * r.tol)
        end
    end

    @testset "rosenbrock_manual() vs fine explicit reference ($FT)" begin
        # Riming-active mixed-phase state: the frozen quadrature rates of the
        # manual Jacobian bound its accuracy relative to the exact Jacobian,
        # and the donor-linearized collision sinks keep refinement monotone
        manual = BMT.RosenbrockAverage(
            BMT.ManualJacobian(), BMT.ExplicitGrowthDiagonal(), BMT.NoLimiter(),
        )
        Δt = FT(10)
        r = (; ρ = FT(0.78), T = TDI.T_freeze(tps) - FT(8), q_tot = FT(0.006),
            x = FT[2e-4, 5e7, 2e-4, 8e4, 3e-4, 4e5, 1e-4, 1.5e-7])
        logλ = consistent_logλ(r.ρ, r.x)
        x_ref = explicit_reference(mp, tps, r.ρ, r.T, r.q_tot, r.x, logλ, Δt, 2048)
        x0 = SVector{8, FT}(r.x...)
        errs = map((1, 4, 16)) do n
            t = BMT.bulk_microphysics_tendencies(
                manual, BMT.Microphysics2Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., logλ, Δt, n,
            )
            err_metric(x0 .+ Δt .* SVector{8, FT}(_applied_2mp3(t)...), x_ref, x0)
        end
        @test all(isfinite, errs)
        @test errs[3] ≤ errs[1] * (1 + sqrt(eps(FT)))
        @test errs[3] < (FT == Float64 ? FT(0.5) : FT(1))
    end

    @testset "degenerate and trivial states ($FT)" begin
        # all-zero state: near-empty species mask -> explicit substeps -> exactly zero
        t0 = BMT.bulk_microphysics_tendencies(
            mode, BMT.Microphysics2Moment(), mp, tps,
            FT(1), FT(273), FT(0),
            FT(0), FT(0), FT(0), FT(0), FT(0), FT(0), FT(0), FT(0),
            FT(-Inf), FT(60), 4,
        )
        @test all(iszero, _applied_2mp3(t0))
        # nsub defaults to 1 and accepts the trailing-argument form
        r = regimes[2]
        logλ = consistent_logλ(r.ρ, r.x)
        t1 = BMT.bulk_microphysics_tendencies(
            mode, BMT.Microphysics2Moment(), mp, tps,
            r.ρ, r.T, r.q_tot, r.x..., logλ, FT(60),
        )
        @test all(isfinite, _applied_2mp3(t1))
    end

    @testset "substeps stay finite and non-negative ($FT)" begin
        # strongly supersaturated over liquid at 233 K (fast
        # condensation-freezing cascade), near-empty rain alongside ice
        x_stress = FT[1e-6, 1e6, 1e-12, 1e-2, 8e-4, 5e5, 5e-4, 9e-7]
        logλ = consistent_logλ(FT(0.45), x_stress)
        for nsub in (1, 2, 8), Δt in (FT(60), FT(300))
            t = BMT.bulk_microphysics_tendencies(
                mode, BMT.Microphysics2Moment(), mp, tps,
                FT(0.45), FT(233), FT(0.003), x_stress..., logλ, Δt, nsub,
            )
            x1 = SVector{8, FT}(x_stress...) .+ Δt .* SVector{8, FT}(_applied_2mp3(t)...)
            @test all(isfinite, x1)
            # non-negative up to the roundoff of the host-side x + Δt * t
            # reconstruction (internally the state is floored at zero)
            tol = eps(FT) .* (abs.(SVector{8, FT}(x_stress...)) .+ Δt .* abs.(SVector{8, FT}(_applied_2mp3(t)...)))
            @test all(x1 .>= -tol)
        end
    end

    @testset "near-empty species take the explicit path ($FT)" begin
        # condensed masses in (eps, 1e-10) produce finite but enormous Jacobian
        # rows; the species mask routes RAIN and ICE to forward Euler so the
        # result tracks the explicit reference and droplet number stays bounded.
        #
        # LIQUID is no longer masked, and this state is why. It is 92 percent
        # supersaturated over liquid, so droplet activation relaxes into the
        # empty liquid species with `1/τ_act` of order 100 per second; an
        # explicit step of `h = 3.75` s overshoots the relaxation target by
        # `h/τ_act`, and the masked path returned 4.2e9 droplets per kg against
        # a target of 1e7. The implicit path with the `-1/τ_act` diagonal
        # returns the target, which is what the bound below now reads.
        x_band = FT[1e-13, 1e2, 0, 0, 0, 0, 0, 0]
        Δt = FT(60)
        x_ref = explicit_reference(mp, tps, FT(1), FT(288), FT(0.02), x_band, FT(-Inf), Δt, 2048)
        x16 = step(x_band, FT(1), FT(288), FT(0.02), FT(-Inf), Δt, 16)
        @test all(isfinite, x16)
        @test x16[2] < 10 * x_ref[2]
        # a near-empty species alongside an active ice species leaves the
        # active species' implicit update bounded
        x_mixb = FT[1e-13, 1e2, 0, 0, 8e-4, 5e5, 5e-4, 9e-7]
        logλ_m = consistent_logλ(FT(0.45), x_mixb)
        xm = step(x_mixb, FT(0.45), FT(253), FT(4e-4), logλ_m, FT(10), 16)
        @test all(isfinite, xm)
        @test xm[2] < FT(1e6)
    end
end

test_rosenbrock_mode(Float64)
test_rosenbrock_mode(Float32)

# Allocation check on the hot call (compiler-version sensitive, like the other
# perf assertions; see performance_tests.jl)
if VERSION >= v"1.12"
    @testset "RosenbrockAverage allocations" begin
        FT = Float64
        tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
        mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
        p3 = mp.ice.scheme
        st = P3.state_from_prognostic(
            p3,
            FT(0.78) * FT(1e-4),
            FT(0.78) * FT(2e5),
            FT(0.78) * FT(4e-5),
            FT(0.78) * FT(6e-8),
        )
        logλ = P3.get_distribution_logλ(st)
        call() = BMT.bulk_microphysics_tendencies(
            BMT.rosenbrock_exact(), BMT.Microphysics2Moment(), mp, tps,
            FT(0.78), FT(273.5), FT(0.009),
            FT(2e-4), FT(5e7), FT(1e-4), FT(4e4), FT(1e-4), FT(2e5), FT(4e-5), FT(6e-8),
            logλ, FT(60), 4,
        )
        call()
        # On <= 1.11 the differentiated path allocates (the known
        # inference-depth limit behind the other >= 1.12 perf assertions), hence
        # the version restriction on this testset.
        @test (@allocated call()) == 0
    end
end

# The manual Jacobian used the raw `g.ρ` while `_per_process_2mp3` clamps its own
# inputs, so a host-delivered non-positive air density reached a `log(ρ q / N)` in the
# condensation timescale and threw. Inside a GPU kernel a throw aborts the whole kernel
# and masks the state that caused it, so the Jacobian must stay finite instead.
# A non-negative clamp is NOT enough: at ρ = 0 the timescale is NaN, which merely moves
# the failure from a throw to a non-finite Jacobian. The floor has to be positive.
function test_jacobian_density_floor(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    mode = BMT.rosenbrock_manual()
    # the manual Jacobian dispatches on MicroState2MP3, not a bare SVector
    x = BMT.MicroState2MP3{FT}(4e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8)
    logλ(ρ) = P3.get_distribution_logλ(
        P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]))

    @testset "manual Jacobian is finite at non-positive air density" begin
        for ρ_bad in FT[-1, FT(-0.245), -floatmin(FT), FT(0)]
            # logλ is diagnosed on the clamped density the entry would use
            g = BMT.Instantaneous2MP3Tendency(
                mp, tps, ρ_bad, FT(273.5), FT(0.009), logλ(max(ρ_bad, FT(1e-4))))
            f, J = BMT._tendency_and_jacobian(mode.jacobian, g, x)
            # BOTH must be finite. The primal used to go non-finite here too, because
            # `_per_process_2mp3` clamped the density to zero rather than flooring it,
            # and the Jacobian then inherited the NaN through the per-process rates -
            # so flooring inside the Jacobian alone was not enough.
            @test all(isfinite, f)
            @test all(isfinite, J)
            # And the PRIMAL as the substep actually calls it. `Instantaneous2MP3Tendency`
            # routes through the full `bulk_microphysics_tendencies` entry, NOT through
            # `_per_process_2mp3`, so flooring the per-process path alone left this one
            # NaN. `_tendency_and_jacobian` returns the per-process sum and so cannot
            # catch it.
            @test all(isfinite, g(x))
        end
    end

    @testset "the density floor is inert at physical densities" begin
        # The realised box density minimum is ~1e-2, well above the floor, so the
        # reference trajectory must not move.
        for ρ_ok in FT[FT(0.011), FT(0.05), FT(0.5), FT(1.0)]
            g = BMT.Instantaneous2MP3Tendency(
                mp, tps, ρ_ok, FT(273.5), FT(0.009), logλ(ρ_ok))
            f₁, J₁ = BMT._tendency_and_jacobian(mode.jacobian, g, x)
            gc = BMT.Instantaneous2MP3Tendency(
                mp, tps, max(ρ_ok, FT(1e-4)), FT(273.5), FT(0.009), logλ(ρ_ok))
            f₂, J₂ = BMT._tendency_and_jacobian(mode.jacobian, gc, x)
            @test J₁ == J₂
            @test f₁ == f₂
        end
    end
end

test_jacobian_density_floor(Float64)
test_jacobian_density_floor(Float32)

# Both explicit fallback branches span the whole substep width, so a stiff phase-change
# rate can convert more mass than the cell holds; the latent release from that then carries
# the substep temperature outside the thermodynamic domain, i.e. the fallback produces
# states worse than the ones it was reached to rescue. Bounding by the cell's total water is
# exact (the condensate total is linear in the increment) and inert on an admissible step.
function test_fallback_water_bound(FT)
    @testset "the water bound is inert on an admissible increment" begin
        x = BMT.MicroState2MP3{FT}(4e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8)
        d = BMT.MicroState2MP3{FT}(1e-6, 1e4, 1e-6, 1e2, 1e-6, 1e2, 1e-7, 1e-10)
        @test BMT._water_bounded_increment(x, d, FT(0.02)) == d
    end

    @testset "the water bound clips an increment that exceeds total water" begin
        x = BMT.MicroState2MP3{FT}(1e-4, 1e7, 1e-4, 1e4, 1e-4, 1e5, 0, 0)
        # an increment that would mint ice far beyond the cell's total water
        d = BMT.MicroState2MP3{FT}(0, 0, 0, 0, FT(0.9), 0, 0, 0)
        q_tot = FT(0.016)
        db = BMT._water_bounded_increment(x, d, q_tot)
        @test BMT._condensate_total(x .+ db) <= q_tot * (1 + 10 * eps(FT))
        @test all(0 .<= Tuple(db ./ ifelse.(iszero.(Tuple(d)), one(FT), Tuple(d))) .<= 1)
    end

    @testset "q_rim is not double-counted in the condensate total" begin
        # rime is a FRACTION of the ice content, so adding q_rim would overstate condensate
        x = BMT.MicroState2MP3{FT}(1e-4, 1e7, 2e-4, 1e4, 3e-4, 1e5, 2e-4, 5e-7)
        @test BMT._condensate_total(x) == x.q_lcl + x.q_rai + x.q_ice
    end

    @testset "no budget supplied leaves the bare step unchanged" begin
        x = BMT.MicroState2MP3{FT}(1e-4, 1e7, 1e-4, 1e4, 1e-4, 1e5, 0, 0)
        d = BMT.MicroState2MP3{FT}(0, 0, 0, 0, FT(0.9), 0, 0, 0)
        @test BMT._bounded_explicit_step(x, d, nothing) === d
    end
end

test_fallback_water_bound(Float64)
test_fallback_water_bound(Float32)

# The branch guards in the collision and threshold paths test `iszero` on quantities that can be
# a differentiated zero: a `ForwardDiff.Dual` whose VALUE is zero but whose partials are not.
# `iszero(::Dual)` requires zero partials, so such a quantity fails the test and selects the
# opposite branch from the plain-float primal - which is how a `0/0` ratio guard gets skipped and
# the exact-AD evaluation's own value lane goes non-finite while the primal is fine.
#
# The guards read `FD.value(...)`, which is exactly inert on plain floats, so a passing primal suite
# proves nothing about them. What distinguishes them is the DUAL evaluation at these corners.
#
# Two things are deliberately NOT asserted, because neither is guaranteed and asserting them would
# reject a correct change (both were tried and both failed against working code):
#
#   * NOT bit-identity of the dual value lane. The dual arithmetic reassociates, so the value lane
#     differs from the primal in the last digits - measured relative differences of 1e-9 to 1e-4 on
#     the `n_lcl` tendency, which is a sum of large cancelling terms. The guards' claim is that the
#     dual takes the SAME BRANCH, and a wrong branch shows up as an order-of-magnitude or sign
#     difference, so an `≈` with a loose tolerance is the assertion that actually tests the claim.
#   * NOT `all(isfinite, J)`. A non-finite exact-AD Jacobian at a degenerate corner is EXPECTED and
#     handled: the substep routes to a forward-Euler update on the primal tendency. Measured, `J` is
#     entirely non-finite at the zero-number-with-mass corner. Asserting otherwise contradicts the
#     design of the fallback.
#
# So the guarantee under test is: the primal is finite, the dual value lane agrees with it to within
# reassociation, and **the substep output is finite** whatever the Jacobian does.
function test_value_lane_branch_guards(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    ρ, T, q_tot = FT(0.78), FT(268), FT(0.009)

    # Each of these puts a different guarded quantity at a differentiated zero.
    states = (
        # rime mass with zero rime volume -> ρ_rim = 0, the B_rim guard
        ("zero rime volume", BMT.MicroState2MP3{FT}(2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 0)),
        # ice mass with zero ice number -> the collision-ratio guards
        ("ice mass, no number", BMT.MicroState2MP3{FT}(2e-4, 5e7, 1e-4, 4e4, 1e-4, 0, 0, 0)),
        # rain mass with zero rain number -> N₀r = 0, the rain-activity guards
        ("rain mass, no number", BMT.MicroState2MP3{FT}(2e-4, 5e7, 1e-4, 0, 1e-4, 2e5, 4e-5, 6e-8)),
        # unrimed ice -> F_rim = 0, the isunrimed threshold selector
        ("unrimed ice", BMT.MicroState2MP3{FT}(2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 0, 0)),
    )

    # One OUTER testset with the per-state ones nested, so a failure names the corner without
    # aborting the remaining corners - a top-level testset throws at its end and ends the file.
    @testset "value-lane branch guards ($FT)" begin
        for (name, x) in states
            @testset "$name" begin
                logλ = P3.get_distribution_logλ(
                    P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]))
                @test isfinite(logλ)
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
                f_primal = g(x)
                f_dual, _J = BMT._tendency_and_jacobian(BMT.ExactJacobian(), g, x)
                @test all(isfinite, f_primal)
                # same BRANCH, not the same rounding: a missed guard is orders of magnitude out
                @test all(
                    isapprox.(Tuple(f_dual), Tuple(f_primal); rtol = FT(1e-3), atol = FT(0)),
                )
                # the guarantee the fallback provides, whatever the Jacobian does
                out = BMT.bulk_microphysics_tendencies(
                    BMT.rosenbrock_exact(), BMT.Microphysics2Moment(), mp, tps, ρ, T, q_tot,
                    x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], logλ, FT(2), 1)
                @test all(isfinite, _applied_2mp3(out))
            end
        end
    end
end

test_value_lane_branch_guards(Float64)
test_value_lane_branch_guards(Float32)

# The accepted implicit increment used to be the ONLY unbounded path out of `_rosenbrock_update`,
# and it is the path the damage takes. Measured on the box: the 2M+P3 substep emitted
# `dq_ice_dt = 2.9e13` in a single cell, minting ice 16 orders beyond that cell's total water,
# through an ACCEPTED increment - because `_rosenbrock_system` equilibrates, so a near-singular
# system's diagonal spike scales out of `_solve_increment_acceptable` and the huge increment passes.
# The rejected branch, which was already bounded, was taken on 0 of 1206 sampled states.
#
# Two things are asserted, and they are the two the change actually guarantees:
#   1. the POSTCONDITION - after the update, condensate never exceeds the cell's total water;
#   2. INERTNESS - on ordinary states the update is bit-identical to the unbounded form, so no
#      trajectory and no sedimentation baseline moves where the bound does not bind.
# Callers that pass no budget (the 1M and temperature-coupled substeps) are inert by construction,
# since `_bounded_explicit_step(x, d, nothing) === d`; that is asserted too.
function test_accepted_increment_is_water_bounded(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    mode = BMT.rosenbrock_manual()
    h = FT(2)
    # The substep takes the rime density bounds explicitly, as production does: the
    # rime-pair projection divides by them, so there is no bounds-free entry to call.
    ρ_min, ρ_max = P3.rime_density_bounds(p3)

    states = (
        ("warm rain", FT(1.05), FT(288), FT(0.015),
            FT[4e-4, 8e7, 2.1e-3, 5e4, 0, 0, 0, 0]),
        ("mixed phase", FT(0.78), FT(273.5), FT(0.009),
            FT[2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8]),
        ("the 56.16 h runaway cell", FT(0.675), FT(278), FT(0.0063),
            FT[4.7e-4, 6.2e5, 1.1e-3, 2.8e3, 5.2e-5, 56.1, 4.2e-5, 6.4e-8]),
        ("giant ice, one crystal", FT(0.5), FT(250), FT(0.006),
            FT[1e-6, 1e-3, 1e-6, 1e-3, 5e-3, 1.0, 4e-3, 1e-5]),
    )

    @testset "accepted increment is water-bounded ($FT)" begin
        for (name, ρ, T, q_tot, xv) in states
            @testset "$name" begin
                x = BMT.MicroState2MP3{FT}(xv...)
                logλ = P3.get_distribution_logλ(P3.state_from_prognostic(
                    p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]))
                isfinite(logλ) || return
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
                f, J_raw = BMT._tendency_and_jacobian(mode.jacobian, g, x)
                J = BMT._apply_growth(mode.growth, J_raw)
                (all(isfinite, f) && all(isfinite, J)) || return
                z = BMT._species_mask(mode.jacobian, mode.growth)(x)

                bounded = first(BMT._rosenbrock_update_diag(x, f, J, z, h, q_tot, ρ_min, ρ_max))
                unbounded = first(BMT._rosenbrock_update_diag(x, f, J, z, h, nothing, ρ_min, ρ_max))

                # 1. the postcondition: condensate never exceeds the available total water
                @test BMT._condensate_total(bounded) <= q_tot * (1 + 64 * eps(FT))
                # 2. inertness where it should not bind: these states are all admissible
                @test all(Tuple(bounded) .=== Tuple(unbounded))
            end
        end

        @testset "no budget leaves the update untouched" begin
            # the 1M and temperature-coupled substeps call without `q_tot`
            x = BMT.MicroState2MP3{FT}(2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8)
            ρ, T, q_tot = FT(0.78), FT(273.5), FT(0.009)
            logλ = P3.get_distribution_logλ(P3.state_from_prognostic(
                p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]))
            g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
            f, J_raw = BMT._tendency_and_jacobian(mode.jacobian, g, x)
            J = BMT._apply_growth(mode.growth, J_raw)
            z = BMT._species_mask(mode.jacobian, mode.growth)(x)
            no_budget = first(BMT._rosenbrock_update_diag(x, f, J, z, h, nothing, ρ_min, ρ_max))
            # a budget far above any condensate this state can reach cannot bind
            slack = first(BMT._rosenbrock_update_diag(x, f, J, z, h, FT(1), ρ_min, ρ_max))
            @test all(Tuple(no_budget) .=== Tuple(slack))
        end

        # A synthetic increment that mints far beyond the budget MUST be clipped, whichever branch
        # produced it. This is the case the box hit and no sampled physical state reproduces.
        @testset "a minting increment is clipped" begin
            x = BMT.MicroState2MP3{FT}(1e-4, 1e7, 1e-4, 1e4, 1e-4, 1e5, 0, 0)
            q_tot = FT(0.016)
            d = BMT.MicroState2MP3{FT}(0, 0, 0, 0, FT(2.9e13), FT(2.9e18), 0, 0)
            db = BMT._water_bounded_increment(x, d, q_tot)
            @test BMT._condensate_total(x .+ db) <= q_tot * (1 + 64 * eps(FT))
            @test db.q_ice < d.q_ice
        end
    end
end

test_accepted_increment_is_water_bounded(Float64)
test_accepted_increment_is_water_bounded(Float32)

# The ice deposition rate must be exactly zero when there is no ice population to deposit onto.
#
# `ice_deposition_timescale` caps its return at `ICE_DEP_TIMESCALE_MAX` to stay finite as the
# capacitance integral vanishes. The unbounded quotient diverges there, so the cap makes `deficit / τ`
# finite on a state with no particles: measured `(qᵥ - qᵥ_sat_ice) / (ICE_DEP_TIMESCALE_MAX ⋅ Γᵢ)`
# exactly, 3.9e-13 kg/kg/s at 220 K. It arrives with no number source, so it manufactures
# mass-without-number states, and it applies to every ice-free supersaturated cell continuously.
#
# The rate is only half of it: the manual Jacobian and the temperature-coupled Jacobian each build
# their own timescale from their own clamped state, so each needs the same gate or the substep
# linearizes a rate that is identically zero. The magnitudes involved cannot bite numerically
# (−1/(τ_max·Γᵢ) ≈ −5e-11 against an implicit diagonal of I/h = 0.5 at h = 2 s); the claim is f/J
# consistency on the `rosenbrock_manual` path, not stability.
#
# What is asserted here is exactly what the gate guarantees, and no more:
#   - at the empty ice state the whole `ice_depsub` slot is identically zero, at every temperature;
#   - the ice-deposition entries of both Jacobians are exactly zero there, and the same derivative
#     evaluated at the capped timescale is NOT zero, so the zero is the gate's doing;
#   - the `Microphysics2Moment` entry and the per-process decomposition return the same ice mass,
#     rime mass and rime volume there, which is the `psum == full` invariant of
#     `rosenbrock_verbose_tests` at the state where the gate acts;
#   - the degeneracy predicate is TRUE there and FALSE at a populated state, which is a
#     branch-identity claim - it is why populated states cannot have moved. Bit-identity against a
#     pre-fix baseline is NOT asserted, because there is no baseline in this suite to assert it against.
function test_degenerate_ice_deposition_is_gated(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    ρ = FT(0.6)
    q_tot = FT(4e-3)
    Ts = FT[220, 233, 245, 253, 258, 261, 263, 265, 268, 273]

    @testset "degenerate ice deposition is gated [FT=$FT]" begin
        empty_state = P3.state_from_prognostic(p3, FT(0), FT(0), FT(0), FT(0))
        logλ₀ = P3.get_distribution_logλ(empty_state)

        @testset "the empty ice state deposits and sublimates nothing" begin
            for T in Ts
                τ = P3.ice_deposition_timescale(
                    mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps,
                    T, ρ, empty_state, logλ₀; quad = mp.ice.quad)
                @test P3.ice_deposition_is_degenerate(τ)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
                    FT(0), FT(0), FT(0), FT(0),   # liquid, rain
                    FT(0), FT(0), FT(0), FT(0),   # ICE-FREE
                    logλ₀)
                # every slot, not only q_ice: q_rim/b_rim are formed from the sublimation branch
                @test all(iszero, Tuple(pp.ice_depsub))
            end
        end

        # A populated state must NOT be degenerate - otherwise the gate would silence real deposition.
        # This is the inertness half of the claim, and it is the half a finiteness assertion would miss.
        @testset "a populated ice state is not degenerate and still deposits" begin
            q_ice, n_ice = FT(1e-5), FT(1e4)
            st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, FT(0), FT(0))
            logλ = P3.get_distribution_logλ(st)
            for T in FT[233, 253, 263]
                τ = P3.ice_deposition_timescale(
                    mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps,
                    T, ρ, st, logλ; quad = mp.ice.quad)
                @test isfinite(τ) && τ > 0
                @test !P3.ice_deposition_is_degenerate(τ)
                @test τ < P3.ICE_DEP_TIMESCALE_MAX(FT)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
                    FT(0), FT(0), FT(0), FT(0),
                    q_ice, n_ice, FT(0), FT(0),
                    logλ)
                # supersaturated with respect to ice at all three temperatures, so growth is positive
                @test pp.ice_depsub.q_ice > 0
                @test isfinite(pp.ice_depsub.q_ice)
                # The Jacobian keeps the deposition self-damping it linearizes: the predicate is
                # false here, so the `ifelse` selects the same branch as before the gate and the
                # whole matrix is bit-identical to the ungated arithmetic by construction.
                x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, q_ice, n_ice, 0, 0)
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
                J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
                @test J[5, 5] < 0 && isfinite(J[5, 5])
                y = BMT.MicroState2MP3T{FT}(0, 0, 0, 0, q_ice, n_ice, 0, 0, T)
                gT = BMT.Temperature2MP3Tendency(mp, tps, ρ, q_tot, logλ)
                ctx = BMT._phase_relaxation_context(mp, tps,
                    _micro_2mp3(q_tot, FT(0), FT(0), FT(0), FT(0), q_ice, n_ice, FT(0), FT(0)),
                    _thermo_2mp3(ρ, T, logλ))
                @test !P3.ice_deposition_is_degenerate(ctx.τ_i)
                J9 = BMT._jacobian_2mp3t_manual(gT, y, pp, rs, ctx)
                @test J9[5, 5] < 0 && isfinite(J9[5, 5])
                # temperature column: −∂q_sat_ice/∂T / τ_i, the saturation shift of a real
                # deposition rate
                @test J9[5, 9] < 0 && isfinite(J9[5, 9])
            end
        end

        # f/J consistency. `_jacobian_2mp3_manual` and `_phase_relaxation_context` each build their
        # own `τ_i` from their own clamped state, so neither reads the gated rate; without the same
        # gate the Jacobian carries a deposition self-derivative against a rate of exactly zero.
        #
        # At the all-zero state every other contribution to row 5 vanishes with its donor
        # (`immersion_freezing` needs cloud liquid, `rain_freezing` needs rain, `ice_melting` needs ice),
        # so those entries are the ice-deposition closure alone. Row 6 is not asserted: its
        # deposition pathway keys on the already-gated rate, and `nice_nice` carries the F23 and
        # number-adjustment relaxations, which are not this closure's.
        @testset "neither Jacobian carries an ice-deposition derivative at the empty state" begin
            x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, 0, 0, 0, 0)
            cp_v = TDI.TD.Parameters.cp_v(tps)
            dcp_dliq = TDI.TD.Parameters.cp_l(tps) - cp_v
            dcp_dice = TDI.TD.Parameters.cp_i(tps) - cp_v
            for T in Ts
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ₀)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot, Tuple(x)..., logλ₀)
                J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
                # rows are (q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim); row 5 is the
                # ice-mass row, columns 1/3/5 its liquid, rain and ice donors
                @test J[5, 1] == 0 && J[5, 3] == 0 && J[5, 5] == 0

                # The zero is the gate's, not the arithmetic's: the same closed-form derivative
                # evaluated at the capped timescale is nonzero on both branches of `_condevap_derivs`
                # (the vapor branch carries −1/(τ_max·Γᵢ), the limited branch the same self term at
                # q_limit = 0).
                τ = P3.ice_deposition_timescale(
                    mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps,
                    T, ρ, empty_state, logλ₀; quad = mp.ice.quad)
                Lₛ_T = TDI.Lₛ(tps, T)
                qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
                dqsi_dT = CMNonEq.dqcld_dT(qᵥ_sat_ice, Lₛ_T, TDI.Rᵥ(tps), T)
                cp_air = TDI.cpₘ(tps, q_tot, FT(0), FT(0))
                Γᵢ = CMNonEq.gamma_helper(Lₛ_T, cp_air, dqsi_dT)
                ci = BMT._condevap_derivs(τ, q_tot - qᵥ_sat_ice, Γᵢ, cp_air, Lₛ_T, dqsi_dT,
                    FT(0), true, dcp_dliq, dcp_dice)
                @test ci.∂s_ice != 0

                # the temperature-coupled 9x9: the species block, the bare-relaxation correction
                # and the temperature column all drop with it
                y = BMT.MicroState2MP3T{FT}(Tuple(x)..., T)
                gT = BMT.Temperature2MP3Tendency(mp, tps, ρ, q_tot, logλ₀)
                ctx = BMT._phase_relaxation_context(mp, tps,
                    _micro_2mp3(q_tot, Tuple(x)...),
                    _thermo_2mp3(ρ, T, logλ₀))
                @test P3.ice_deposition_is_degenerate(ctx.τ_i)
                J9 = BMT._jacobian_2mp3t_manual(gT, y, pp, rs, ctx)
                @test J9[5, 1] == 0 && J9[5, 3] == 0 && J9[5, 5] == 0
                @test J9[5, 9] == 0
                # the temperature tendency's own deposition correction is `(Γᵢ − 1)·pp.ice_depsub`,
                # so it is already zero wherever the slot is gated
                @test all(iszero, Tuple((ctx.Γᵢ - 1) * pp.ice_depsub))
            end
        end

        # The `Microphysics2Moment` entry computes the same rate a third time, and it is what
        # `Instantaneous2MP3Tendency` evaluates - so the exact-AD substep and the per-process
        # decomposition read different rates unless it is gated too. `rosenbrock_verbose_tests`
        # asserts `psum == full` on three populated regimes; this is the same invariant at the state
        # where the gate acts.
        @testset "the entry and the per-process decomposition agree at the empty state" begin
            x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, 0, 0, 0, 0)
            for T in Ts
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ₀)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot, Tuple(x)..., logλ₀)
                @test g(x).q_ice == sum(values(pp)).q_ice
                @test g(x).q_rim == sum(values(pp)).q_rim
                @test g(x).b_rim == sum(values(pp)).b_rim
            end
        end
    end
end

test_degenerate_ice_deposition_is_gated(Float64)
test_degenerate_ice_deposition_is_gated(Float32)

# The rain-freezing branch must draw on every nucleation pathway and respect how long the drop
# then takes to actually become ice.
#
# `liquid_freezing_rate` gives the Bigg (1953) NUCLEATION rate and nothing else: one pathway, and a
# drop that has nucleated is treated as ice in the same instant. `rain_freezing_rate` adds the Koop
# (2000) homogeneous pathway in PARALLEL (coefficients add, through one shared PSD treatment) and
# composes the result in SERIES with `drop_freezing_heat_timescale` and
# `drop_freezing_dendrite_timescale`.
#
# What is asserted here is exactly what the composition guarantees, and the guarantees are
# narrower than one might expect from a limiter:
#   - `τ_heat` is finite, non-negative, decreasing in supercooling and increasing in drop mass;
#   - `τ_heat` carries the FULL Musil wet-growth balance, conduction plus evaporation at the
#     melting-point surface, so it shortens with ambient dryness. The enhancement over the
#     conduction-only denominator is bracketed at three supercoolings, the humidity dependence is
#     asserted monotone, and the clamp that stops a wetter-than-surface ambient turning evaporation
#     into a heat SOURCE is asserted to return the conduction-only answer exactly. Nothing past
#     `ΔT*` moves: the numerator is exactly zero there and `===` says so;
#   - it is EXACTLY ZERO at and beyond the recalescence threshold ΔT*, where the drop's own cold
#     absorbs all of `Lf`. ΔT* is bisected on the actual `TDI.Lf` in the setup rather than quoted,
#     because `Lf` is temperature dependent and that moves the threshold well below `Lf/c_w`;
#   - beyond ΔT* only `drop_freezing_dendrite_timescale` is left, and it is what reaches the deep
#     supercooling where the 2.9e13 kg/kg/s mint was measured. BIT-IDENTITY IS NOT ASSERTED
#     ANYWHERE, and the two-stage form's claim to it is withdrawn: `τ_dend` is positive on every
#     populated state, so the rate is never bit-identical to the unlimited one except where that
#     rate is already exactly zero. At weak supercooling what holds is a MEASURED relative
#     perturbation, and the tests take their tolerance from it;
#   - the composition is PER SIZE, inside the PSD integrals, by a fixed-node Gauss-Laguerre rule
#     whose `exp(-u)` weight is the exponential PSD itself. So there is no single reduction
#     factor any more and `1/(1 + (τ_heat + τ_dend)/τ_nuc)` is NOT asserted. What is asserted:
#     the rule is exact in the zero-limiting limit (integrands `u³` and `u⁶`, recovering the
#     analytic moments), its shipped table satisfies the defining degree-`2n-1` exactness, the
#     production order is converged against an independently built 32-node rule over the band
#     where the composition binds, each moment's reduction lies in `(0, 1]` and falls with
#     supercooling ACROSS THE BAND, and the MASS reduction is below the NUMBER reduction
#     everywhere. Monotonicity is claimed on the band only, and the restriction is measured, not
#     convenient: past `ΔT*` the numerator saturates (`τ_dend` carries no temperature) while the
#     denominator stalls at the Koop cold-edge clamp, so the ratio rises across `ΔT*` and falls
#     again once the extrapolated Bigg coefficient overtakes the clamp. The mean-drop form has
#     the identical shape; its own test hid it by sorting by the ratio rather than by `ΔT`;
#   - the mean frozen drop mass `∂ₜq/∂ₜn` is NOT preserved, and its motion is the second thing
#     the per-size form fixes: it starts near the volume-selective `20 x̄`, descends monotonically
#     through `x̄` inside the band, and stays below it past `ΔT*`, with no branch producing the
#     transition. The mean-drop form held `20 x̄` at every temperature, including where every drop
#     is freezing and the answer must be `x̄`;
#   - `∂ₜn_frz|nuc == n/τ_nuc` exactly for the SB2006 exponential rain PSD, which is what makes
#     evaluating `τ_nuc` at the mean-mass drop exact for the number moment rather than an
#     approximation;
#   - the per-process sum still equals the full entry at a state where the limiter binds;
#   - the Koop pathway is EXACTLY off warmer than its fitted window - a hard zero, not the edge
#     value - so `J_bigg + J_koop === J_bigg` there and adding it cannot perturb ordinary
#     mixed-phase behaviour at all. That exactness claim IS made, unlike the τ_dend one, because
#     the warm-side treatment is a hard zero rather than a small number;
#   - the shipped Koop coefficients give J in SI at both window edges, so a stray factor of 1e6
#     in either direction fails a test rather than hiding;
#   - the wrapper never throws, at any temperature. `homogeneous_J_cubic` raises a DomainError
#     outside its window, and a DomainError inside a GPU kernel is the crash family this whole
#     campaign has been chasing.
#
# DELIBERATELY NOT ASSERTED:
#   - that the limited RATE is bounded by the donor over a step. It is not, and requiring it would
#     be the wrong claim: a relaxation with τ_eff ≪ dt necessarily has rate ⋅ dt ≫ q. Boundedness
#     is asserted where it actually lives, by driving the production substep entry and bounding the
#     realised increment - see "the vertex state is donor-bounded THROUGH the substep".
#   - that `τ_heat > 0` at the vertex/mint state. It is exactly zero there, by construction, and
#     `τ_dend` is what carries the reduction on those states.
# An `n`-point Gauss-Laguerre rule, built here rather than imported, so the convergence testset
# compares the shipped table against an INDEPENDENT construction and not against a longer copy of
# itself. Nodes are the roots of the Laguerre polynomial `Lₙ`, found by Newton iteration from the
# standard `gaulag` initial guesses, with `wᵢ = xᵢ / ((n+1)² Lₙ₊₁(xᵢ)²)`. Float64 throughout; the
# rate converts at use, exactly as it does for the production table.
function _gauss_laguerre_reference(n::Int)
    # (Lₖ(x), Lₖ₋₁(x)) by the three-term recurrence k Lₖ = (2k−1−x) Lₖ₋₁ − (k−1) Lₖ₋₂,
    # iteratively: the same recurrence written recursively is exponential in k.
    function Lpair(k, x)
        p1, p0 = 1.0 - x, 1.0
        k == 0 && return (1.0, 0.0)
        for j in 2:k
            p1, p0 = ((2j - 1 - x) * p1 - (j - 1) * p0) / j, p1
        end
        return (p1, p0)
    end
    xs = zeros(n)
    ws = zeros(n)
    for i in 1:n
        z =
            i == 1 ? 3 / (1 + 2.4n) :
            i == 2 ? xs[1] + 15 / (1 + 2.5n) :
            xs[i - 1] + ((1 + 2.55 * (i - 2)) / (1.9 * (i - 2))) * (xs[i - 1] - xs[i - 2])
        for _ in 1:100
            p, pm1 = Lpair(n, z)
            dp = n * (p - pm1) / z          # Lₙ′(x) = n (Lₙ(x) − Lₙ₋₁(x)) / x
            dz = p / dp
            z -= dz
            abs(dz) < 1e-14 * abs(z) && break
        end
        xs[i] = z
        ws[i] = z / ((n + 1)^2 * first(Lpair(n + 1, z))^2)
    end
    return (; nodes = Tuple(xs), weights = Tuple(ws))
end

function test_rain_freezing_is_heat_limited(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    pdf_r = mp.ice.rain_pdf
    evap = mp.warm_rain.seifert_beheng.evap
    aps = mp.warm_rain.air_properties
    T_frz = TDI.T_freeze(tps)
    c_w = TDI.cp_l(tps)
    (; ρw) = pdf_r

    # ΔT*, where the recalescence numerator `Lf(T_frz - ΔT) - c_w ΔT` reaches zero. Bisected on
    # the code's own `TDI.Lf`, so the threshold the tests use is the one the source computes.
    resid(ΔT) = TDI.Lf(tps, T_frz - ΔT) - c_w * ΔT
    lo, hi = FT(0), FT(150)
    for _ in 1:80
        mid = (lo + hi) / 2
        resid(mid) > 0 ? (lo = mid) : (hi = mid)
    end
    # `hi` is the bracket side on which `resid ≤ 0`, so the numerator's `max(⋅, 0)` is already
    # clamped there. Taking the midpoint instead would land a hair on either side and make the
    # "exactly zero" assertions below a coin flip at Float32.
    ΔT_star = hi

    ρ = FT(0.9)
    q_rai = FT(2e-4)
    # `N_rai` chosen so that x̄ = ρ q / N lands on a millimetre drop and on a much smaller one. The
    # SB2006 limited PDF clamps N₀ and λ as well as x̄, so the second target is NOT reached - its
    # N₀ clamp binds and returns a ~0.4 mm drop rather than the ~0.1 mm one requested. Every
    # assertion therefore uses the ACHIEVED x̄, never the target, and the two states are only
    # claimed to be a larger and a smaller drop. `diag/rainfrz_tauheat_magnitudes.jl` reaches the
    # true drizzle end by calling the per-drop function directly.
    Ns = (FT(ρ * q_rai / 8e-7), FT(ρ * q_rai / 5e-10))
    x̄(N) = CM2.pdf_rain_parameters(pdf_r, q_rai, ρ, N).xr_mean
    v_drop(N) = evap.α * x̄(N)^evap.β * sqrt(evap.ρ0 / ρ)
    # The ambient humidity the wet-growth balance's evaporative term reads. A LIQUID-SATURATED
    # environment is the reference the band factors are quoted for, so the sweeps use it; the
    # humidity dependence itself, and the clamp that guards it, get their own testset.
    qᵥ_sat_liq(T) = TDI.p2q(tps, T, ρ, TDI.saturation_vapor_pressure_over_liquid(tps, T))
    τ_h(N, T) = CM_HetIce.drop_freezing_heat_timescale(
        p3.vent, aps, tps, T, ρ, qᵥ_sat_liq(T), x̄(N), ρw, v_drop(N))
    τ_d(N) = CM_HetIce.drop_freezing_dendrite_timescale(x̄(N), ρw)
    hom = mp.ice.homogeneous
    # Bigg-only, the pre-Koop nucleation rate
    raw(N, T) = CM_HetIce.liquid_freezing_rate(mp.ice.rain_freezing, pdf_r, tps, q_rai, ρ, N, T)
    lim(N, T) = CM_HetIce.rain_freezing_rate(
        mp.ice.rain_freezing, hom, p3.vent, aps, tps, evap, pdf_r, q_rai, ρ, N, T, qᵥ_sat_liq(T))
    # the UNLIMITED rate at the summed coefficient - the correct denominator for the reduction,
    # since the limiter divides that and not the Bigg-only rate
    function nuc(N, T)
        l = lim(N, T)
        return CM_HetIce._liquid_freezing_rate_from_J(
            pdf_r, l.J_bigg + l.J_koop, q_rai, ρ, N, T, T_frz)
    end
    J_koop_at(T) = CM_HetIce.homogeneous_freezing_rate_coefficient(hom, tps, T)
    # the composed post-nucleation delay, relative to nucleation - the only knob in the reduction
    ratio_of(N, T) = (τ_h(N, T) + τ_d(N)) / lim(N, T).τ_nuc

    @testset "rain freezing is heat-dissipation limited [FT=$FT]" begin
        @testset "recalescence threshold" begin
            # Strictly inside the mixed-phase range and far from the naive constant-Lf estimate.
            @test ΔT_star > 40
            @test ΔT_star < 70
            @test ΔT_star < FT(0.8) * TDI.Lf(tps, T_frz) / c_w
        end

        @testset "τ_heat, drop mass $(round(Float64(x̄(N)), sigdigits = 3)) kg" for N in Ns
            # Below ΔT*: positive, finite, and strictly decreasing in supercooling, because both
            # the heat left to shed and the conduction driving temperature move the same way.
            below = FT[5, 10, 15, 20, 25, 30, FT(0.9) * ΔT_star]
            τs = [τ_h(N, T_frz - ΔT) for ΔT in below]
            @test all(isfinite, τs)
            @test all(>(0), τs)
            @test issorted(τs; rev = true)

            # At and beyond ΔT*: τ_heat is exactly zero, and τ_dend is then the ONLY stage left
            # bounding the conversion. NOT bit-identity to the raw rate - that was the two-stage
            # form's property and it is deliberately gone, because passing the raw Bigg rate
            # through at 82 K of supercooling is the defect this stage exists to remove.
            for ΔT in FT[ΔT_star, ΔT_star + 1, 60, 82, 100]
                T = T_frz - ΔT
                @test τ_h(N, T) === zero(FT)
                @test τ_d(N) > 0
                r, l = raw(N, T), lim(N, T)
                @test l.τ_heat === zero(FT)
                @test l.τ_dend === τ_d(N)
                @test l.∂ₜn_frz < r.∂ₜn_frz
                @test l.∂ₜq_frz < r.∂ₜq_frz
                @test all(isfinite, (l.∂ₜn_frz, l.∂ₜq_frz))
            end

            # Above freezing there is nothing to limit and nothing to divide by.
            @test τ_h(N, T_frz) === zero(FT)
            @test τ_h(N, T_frz + 10) === zero(FT)
        end

        @testset "τ_heat grows with drop mass" begin
            # A bigger drop holds more latent heat behind less surface per unit mass, so it takes
            # longer to freeze through at the same supercooling.
            for ΔT in FT[10, 20, 30]
                T = T_frz - ΔT
                @test τ_h(Ns[1], T) > τ_h(Ns[2], T)  # Ns[1] is the millimetre drop
            end
        end

        @testset "the wet-growth balance carries its evaporative term" begin
            # The melting-point drop surface is warmer AND wetter than the environment, so a
            # freezing drop evaporates and each evaporated kilogram exports Lᵥ, about seven times
            # Lf. That channel is a second, parallel heat sink and it SHORTENS the freezing time.
            #
            # `Δρᵥ = 0` is the exact conduction-only limit of the same expression, reached by
            # handing the function the surface's own vapour density as the ambient one. So the
            # enhancement can be measured against the landed form without keeping a second copy
            # of it: the ratio below is the full Musil denominator over the conduction-only one,
            # up to the Prandtl-for-Schmidt swap, which is asserted separately to be inert.
            τ_h_at(N, T, qᵥ) = CM_HetIce.drop_freezing_heat_timescale(
                p3.vent, aps, tps, T, ρ, qᵥ, x̄(N), ρw, v_drop(N))
            ρᵥ_sfc = ρ * TDI.p2q(tps, T_frz, ρ,
                TDI.saturation_vapor_pressure_over_liquid(tps, T_frz))
            qᵥ_no_evap = ρᵥ_sfc / ρ   # Δρᵥ = 0 exactly: the conduction-only limit

            # The magnitudes are pre-registered rather than recomputed from the formula: a
            # liquid-saturated environment gives roughly 1.6 / 1.4 / 1.3 at ΔT = 10 / 25 / 40.
            # Brackets, not point values, because they carry the shipped constants.
            for (ΔT, lo, hi) in ((FT(10), FT(1.50), FT(1.70)),
                (FT(25), FT(1.30), FT(1.50)),
                (FT(40), FT(1.20), FT(1.40)))
                T = T_frz - ΔT
                enh = τ_h_at(Ns[1], T, qᵥ_no_evap) / τ_h(Ns[1], T)
                @test lo < enh < hi
            end

            # The physical humidity dependence the conduction-only form could not represent:
            # drier air freezes drops faster, monotonically, at every supercooling below ΔT*.
            for ΔT in FT[10, 25, 40]
                T = T_frz - ΔT
                τs = [τ_h_at(Ns[1], T, FT(f) * qᵥ_no_evap) for f in (0, 0.25, 0.5, 0.75, 1)]
                @test issorted(τs)
                @test all(>(0), τs)
                @test all(isfinite, τs)
            end

            # THE CLAMP. `Δρᵥ` is clamped at zero so that an ambient wetter than the melting-point
            # surface cannot turn evaporation into a heat SOURCE and shorten the freezing time
            # without limit. It takes a liquid supersaturation ratio above 2 at ΔT = 10 K to
            # reach, so it is unreachable in cloud and present for corrupt `q_tot` only - the
            # assertion is that the clamped branch returns the conduction-only answer exactly,
            # not something merely finite.
            for ΔT in FT[10, 25, 40], f in FT[1.001, 2, 100, 1e6]
                T = T_frz - ΔT
                @test τ_h_at(Ns[1], T, f * qᵥ_no_evap) === τ_h_at(Ns[1], T, qᵥ_no_evap)
            end

            # THE VERTEX MUST NOT MOVE, and it cannot: past ΔT* the numerator `max(Lf - c_w ΔT, 0)`
            # is exactly zero, so every denominator gives exactly zero. Asserted with `===` across
            # a humidity range spanning the clamp, because "the deep-supercooled falsifiers are
            # untouched by this commit" is the claim the whole change is landed against.
            for ΔT in FT[ΔT_star, ΔT_star + 1, 60, 82, 100], N in Ns
                T = T_frz - ΔT
                for qᵥ in FT[0, 1e-6, 1e-3, 1] .* one(FT)
                    @test τ_h_at(N, T, qᵥ) === zero(FT)
                end
                # and the whole composed rate is bit-identical between a bone-dry and a
                # saturated ambient, since τ_heat is the only place the humidity enters
                dry = CM_HetIce.rain_freezing_rate(mp.ice.rain_freezing, hom, p3.vent, aps,
                    tps, evap, pdf_r, q_rai, ρ, N, T, zero(FT))
                wet = CM_HetIce.rain_freezing_rate(mp.ice.rain_freezing, hom, p3.vent, aps,
                    tps, evap, pdf_r, q_rai, ρ, N, T, qᵥ_no_evap)
                @test dry.∂ₜn_frz === wet.∂ₜn_frz
                @test dry.∂ₜq_frz === wet.∂ₜq_frz
            end

            # Retiring the Schmidt-for-Prandtl borrow is a correctness-of-form change with no
            # measurable magnitude, and that is worth pinning: the shipped ν_air and D_vapor give
            # N_sc = 0.708 against a material Prandtl number of 0.71, so the ventilation factors
            # differ by less than 0.1 % in the cube root. If a constant moves so that this stops
            # being true, the two-ventilation-factor split starts to matter and should be reviewed
            # rather than silently carried.
            N_sc = aps.ν_air / aps.D_vapor
            N_pr = CM_HetIce.PRANDTL_NUMBER_AIR(FT)
            @test FT(0.6) < N_pr < FT(0.8)
            @test abs(cbrt(N_pr) / cbrt(N_sc) - 1) < FT(1e-3)
        end

        @testset "τ_dend is D/v_dend and grows with drop mass" begin
            v_dend = CM_HetIce.DENDRITE_GROWTH_VELOCITY(FT)
            @test v_dend > 0
            rtol = FT == Float64 ? FT(1e-12) : FT(1e-5)
            for N in Ns
                D = cbrt(6 * x̄(N) / (FT(π) * ρw))
                @test isapprox(τ_d(N), D / v_dend; rtol)
                @test τ_d(N) > 0
                # milliseconds for a raindrop, so h/τ_dend ≫ 1 at any model timestep
                @test τ_d(N) < FT(0.1)
            end
            @test τ_d(Ns[1]) > τ_d(Ns[2])  # Ns[1] is the larger drop
            # temperature-independent, unlike τ_heat: the kinetic plateau carries no ΔT
            @test τ_d(Ns[1]) === τ_d(Ns[1])
        end

        @testset "both post-nucleation stages vanish with the population" begin
            # An empty rain state gives x̄ = 0, where `x_drop/D` is 0/0 if written naively.
            @test CM_HetIce.drop_freezing_heat_timescale(
                p3.vent, aps, tps, FT(250), ρ, qᵥ_sat_liq(FT(250)),
                zero(FT), ρw, zero(FT)) === zero(FT)
            # ...and with no vapour in the air either, where the evaporative term is largest
            @test CM_HetIce.drop_freezing_heat_timescale(
                p3.vent, aps, tps, FT(250), ρ, zero(FT), zero(FT), ρw, zero(FT)) === zero(FT)
            @test CM_HetIce.drop_freezing_dendrite_timescale(zero(FT), ρw) === zero(FT)
            for (q_s, n_s) in ((FT(0), FT(0)), (FT(2e-4), FT(0)), (FT(0), FT(5000)))
                N_s = n_s * ρ
                l = CM_HetIce.rain_freezing_rate(
                    mp.ice.rain_freezing, hom, p3.vent, aps, tps, evap, pdf_r, q_s, ρ, N_s,
                    FT(250), qᵥ_sat_liq(FT(250)))
                @test isfinite(l.τ_heat)
                @test l.∂ₜn_frz === zero(FT)
                @test l.∂ₜq_frz === zero(FT)
            end
        end

        @testset "the mean-mass drop reproduces the PSD number rate exactly" begin
            # ∂ₜn_frz|nuc = J (π/6) M_D³ = J π n Dr_mean³ = n/τ_nuc, because x̄ = π ρw Dr_mean³ for
            # the SB2006 exponential PSD. This identity is what makes evaluating τ_nuc at the
            # mean-mass drop exact for the number moment; if the PSD convention moves, it breaks.
            rtol = FT == Float64 ? FT(1e-12) : FT(1e-5)
            for N in Ns, ΔT in FT[10, 25, 40]
                T = T_frz - ΔT
                n = N / ρ
                # `nuc`, not `raw`: τ_nuc is built from J_bigg + J_koop, and by ΔT = 40 K the
                # homogeneous pathway dominates, so comparing against Bigg alone would be wrong.
                @test isapprox(nuc(N, T).∂ₜn_frz, n / lim(N, T).τ_nuc; rtol)
            end
        end

        @testset "weak supercooling leaves the NUMBER rate alone, but not the mass rate" begin
            # The mean-drop form's claim was that at weak supercooling the whole rate is
            # essentially unchanged, with the perturbation measured from τ_heat/τ_nuc at the mean
            # drop. Per size that survives for the NUMBER moment and fails for the MASS moment,
            # and the failure is physics rather than an error: ∂ₜq is weighted by D⁶, so it is
            # carried by drops several times the PSD scale, and those cross into the heat-limited
            # regime far warmer than the mean drop does. Measured, the mass rate is already about
            # 9 % below the unlimited one at ΔT = 5 K, where the mean-drop factor is 1 to four
            # decimal places. That effect was previously invisible, not absent.
            for N in Ns, ΔT in FT[5, 10]
                T = T_frz - ΔT
                ratio = ratio_of(N, T)
                @test ratio < FT(0.05)  # the mean-drop diagnostic: still inert this warm
                # τ_dend must be subdominant here, or the third stage is reaching too warm
                @test τ_d(N) < FT(0.01) * τ_h(N, T)
                fn = lim(N, T).∂ₜn_frz / nuc(N, T).∂ₜn_frz
                fq = lim(N, T).∂ₜq_frz / nuc(N, T).∂ₜq_frz
                @test FT(0.75) < fn ≤ 1
                @test FT(0.35) < fq < fn   # strictly more limited, at every state
            end
        end

        @testset "the per-size reduction is bounded and ordered between the moments" begin
            # The scalar identity `reduction = 1/(1 + (τ_heat+τ_dend)/τ_nuc)` is GONE: with the
            # composition inside the integrals there is no single factor, and the two moments no
            # longer share one. What replaces it:
            #   - each moment's reduction is in (0, 1]: `r(D) = 1/(τ_nuc + τ_freeze)` is strictly
            #     below the nucleation rate `1/τ_nuc` at every size, at every temperature;
            #   - the MASS reduction is below the NUMBER reduction at every supercooling, because
            #     the mass integrand's D⁶ weight puts it on the heat-limited large drops;
            #   - both fall monotonically with supercooling ACROSS THE BAND, and the restriction
            #     is measured rather than convenient - see the next testset.
            grid = FT[5, 10, 15, 20, 25, 30, 40, FT(0.95) * ΔT_star, 60, 82]
            for N in Ns
                fs = map(grid) do ΔT
                    T = T_frz - ΔT
                    l, u = lim(N, T), nuc(N, T)
                    (l.∂ₜn_frz / u.∂ₜn_frz, l.∂ₜq_frz / u.∂ₜq_frz)
                end
                for (fn, fq) in fs
                    @test 0 < fn ≤ 1
                    @test 0 < fq ≤ 1
                    @test fq < fn
                end
                # Monotone on the band, where the denominator still responds to temperature.
                band = fs[1:findfirst(==(FT(40)), grid)]
                @test issorted([fn for (fn, _) in band]; rev = true)
                @test issorted([fq for (_, fq) in band]; rev = true)
            end
        end

        @testset "past ΔT* the conversion saturates and the reduction stops tracking ΔT" begin
            # MEASURED, and the reason the monotonicity above is stated on the band only. The
            # reduction is a ratio to the UNLIMITED rate at the summed coefficient, and past ΔT*
            # neither part of that ratio behaves like a function of ΔT:
            #
            #   - the numerator SATURATES. With τ_heat exactly zero only τ_dend is left, and it
            #     carries no temperature, so the per-size rate is the same at ΔT = 60 and 82
            #     despite J_bigg growing by six orders between them. That is the whole point of
            #     the third stage. NOT bit-identical, and the residual is DERIVED rather than
            #     tolerated - see the tolerance below.
            #   - the denominator STALLS and then resumes. J_koop is clamped at its cold window
            #     edge from ΔT ≈ 43 K, and until the extrapolated J_bigg overtakes it near 70 K the
            #     summed coefficient barely moves. So across ΔT*, where the numerator jumps as
            #     τ_heat vanishes, the ratio RISES: measured 2.2e-14 at 0.95 ΔT* against 9.4e-13 at
            #     60 K, then falling again to 8.7e-14 at 82 K once Bigg dominates the sum.
            #
            # None of this is a property of the per-size composition - the mean-drop form has the
            # identical shape, and its own test never showed it because that test sorted BY THE
            # RATIO rather than by ΔT. Pinned here so the carve-out above is a measurement.
            # THE TOLERANCE IS DERIVED FROM THE QUADRATURE, NOT CHOSEN. Past ΔT* the heat stage
            # is gone, so at every node `r(D) = 1/(τ_nuc(D) + τ_dend(D))` and the ONLY remaining
            # ΔT dependence is `τ_nuc`. The 60-vs-82 gap is therefore bounded by the largest
            # `τ_nuc/τ_dend` over the nodes, taken at the WARMER state where `τ_nuc` is larger;
            # and since that ratio goes as `D⁻⁴` the maximum sits at the SMALLEST node,
            # `D₁ = u₁ Dr_mean`. So the test PREDICTS its own residual instead of absorbing it.
            #
            # Measured (`diag/saturation_residual_check.jl`, job 6941483, both precisions): the
            # bound over-predicts the gap by a STATE-INDEPENDENT factor 1.38, across two states
            # whose residuals differ by four orders - 2.5e-9 at the millimetre drop and 4.6e-5 at
            # fresh drizzle. The factor 2 below is therefore ~45 % headroom everywhere, not a
            # guess.
            #
            # A LITERAL CANNOT DO THIS JOB, and that is what retiring the clamp cascade exposed.
            # The mean-mass window lets a fresh-drizzle state follow `L/N` down to a mean size the
            # cascade used to pin at `λ_max`, and `τ_nuc ∝ D⁻³` at the smallest node turned a
            # 1.8e-9 residual into 3.3e-5. No single literal covers both: the value that admits
            # drizzle is five orders looser than the millimetre state deserves, and the value that
            # pins the millimetre state fails on drizzle. The derived form is also tighter than
            # the literal it replaces where it matters - 4.9e-9 against 1e-6 at the millimetre
            # drop, a 200× stronger assertion - and it follows the PSD if the window's ends move.
            #
            # Calibrated on the NUMBER moment, which is where the residual lives. The mass moment
            # sits two orders inside the same bound (2.3e-7 against 4.6e-5 at drizzle), because
            # its `D³` weight moves the integral off the smallest node - the same asymmetry the
            # quadrature-convergence testset measures.
            u₁ = FT(CM_HetIce.rain_freezing_quadrature().nodes[1])
            v_dend = CM_HetIce.DENDRITE_GROWTH_VELOCITY(FT)
            for N in Ns
                deep = map(T -> lim(N, T), (T_frz - FT(60), T_frz - FT(82)))
                @test deep[1].τ_heat === zero(FT)
                @test deep[2].τ_heat === zero(FT)
                D₁ = CM2.pdf_rain_parameters(pdf_r, q_rai, ρ, N).Dr_mean * u₁
                J₆₀ = deep[1].J_bigg + deep[1].J_koop
                τ_nuc₁ = 1 / (J₆₀ * FT(π) / 6 * D₁^3)
                rtol_sat = 2 * τ_nuc₁ / (D₁ / v_dend)
                # J-independent past ΔT*, to the smallest node's own nucleation residual
                @test isapprox(deep[1].∂ₜn_frz, deep[2].∂ₜn_frz; rtol = rtol_sat)
                @test isapprox(deep[1].∂ₜq_frz, deep[2].∂ₜq_frz; rtol = rtol_sat)
                # ...while the unlimited rate is still growing steeply between the two
                @test nuc(N, T_frz - FT(82)).∂ₜn_frz > FT(10) * nuc(N, T_frz - FT(60)).∂ₜn_frz
                # and the rise across ΔT* is in the denominator's stall, not in the rate
                near = lim(N, T_frz - FT(0.95) * ΔT_star)
                @test near.τ_heat > 0
                @test deep[1].∂ₜn_frz > near.∂ₜn_frz
            end
        end

        @testset "the quadrature is exact in the zero-limiting limit" begin
            # THE load-bearing assertion of the per-size form, and the sibling of the mean-drop
            # form's exactness identity. Turn the two later stages off and the integrands become
            # exactly `u³` and `u⁶` against the exponential PSD's own `exp(-u)` weight, so an
            # n-node Gauss-Laguerre rule reproduces the analytic moments `J (π/6) M_D³` and
            # `J ρw (π/6)² M_D⁶` to round-off for every n ≥ 4. Driven through the rule itself
            # rather than through `rain_freezing_rate`, because `τ_dend` is positive on every
            # populated state and cannot be switched off from outside.
            (; nodes, weights) = CM_HetIce.rain_freezing_quadrature()
            @test length(nodes) == CM_HetIce.RAIN_FREEZING_QUADRATURE_ORDER
            @test length(weights) == length(nodes)

            # The defining property of the n-point rule: ∫₀^∞ exp(-u) u^k du = k! exactly for
            # every k ≤ 2n-1. This is what certifies the hard-coded table; a mistyped digit or a
            # node from the wrong order fails it long before any physics test would.
            nq = length(nodes)
            for k in 0:(2 * nq - 1)
                approx = sum(w * u^k for (u, w) in zip(nodes, weights))
                @test isapprox(approx, factorial(big(k)); rtol = 1e-10)
            end

            # ...and the same statement in the moments the rate actually forms.
            rtol = FT == Float64 ? FT(1e-10) : FT(1e-3)
            for N in Ns
                Dr = CM2.pdf_rain_parameters(pdf_r, q_rai, ρ, N).Dr_mean
                n = N / ρ
                for k in (3, 6)
                    quadmom = sum(FT(w) * (Dr * FT(u))^k for (u, w) in zip(nodes, weights)) * n
                    @test isapprox(quadmom, n * FT(factorial(k)) * Dr^k; rtol)
                end
            end
        end

        @testset "the quadrature order is converged where the composition binds" begin
            # `RAIN_FREEZING_QUADRATURE_ORDER` is 8, and the count is a measurement rather than a
            # preference. Compared against a much finer rule built here, the production rule must
            # be converged over the band where the per-size composition changes the answer. It is
            # NOT converged past ΔT*, where the number integrand's peak drops below the smallest
            # node; that is documented on the constant, the error is a one-signed under-estimate,
            # and it is inert through the implicit update, so the reference sweep stops at 30 K
            # rather than pretending otherwise.
            fine = _gauss_laguerre_reference(32)
            for N in Ns, ΔT in FT[5, 10, 15, 20, 25, 30]
                T = T_frz - ΔT
                l = lim(N, T)
                f = CM_HetIce.rain_freezing_rate(mp.ice.rain_freezing, hom, p3.vent, aps, tps,
                    evap, pdf_r, q_rai, ρ, N, T, qᵥ_sat_liq(T); quad = fine)
                @test isapprox(l.∂ₜq_frz, f.∂ₜq_frz; rtol = FT(0.05))
                @test isapprox(l.∂ₜn_frz, f.∂ₜn_frz; rtol = FT(0.15))
            end
        end

        @testset "the event mass slides from 20 x̄ toward and below x̄" begin
            # THE second caveat the per-size composition exists to fix, and the reason it is worth
            # the quadrature. `∂ₜq_frz/∂ₜn_frz` is the mean mass converted per freezing event.
            #
            # Unlimited, it is exactly 20 x̄: freezing is volume-selective, `∫ x r n` weights D⁶
            # against `∫ r n`'s D³, and `M_D⁶/M_D³ = 120 Dr³` for the exponential PSD. That is the
            # right answer when freezing is RARE - only the biggest drops go - and the wrong one
            # when every drop is freezing, where the converted mass has to relax toward x̄. The
            # mean-drop form divided both moments by one factor and so preserved 20 x̄ at every
            # temperature, including the deep limit where it is most wrong.
            #
            # Normalised by `π ρw Dr_mean³` rather than by `xr_mean`: the two agree only when the
            # SB2006 limited PDF's clamps are slack, and the analytic `20` is exact in the first.
            for N in Ns
                Dr = CM2.pdf_rain_parameters(pdf_r, q_rai, ρ, N).Dr_mean
                x_scale = FT(π) * ρw * Dr^3

                # the unlimited rate really does carry 20 x̄, which is what the descent starts from
                u = nuc(N, T_frz - FT(10))
                @test isapprox(u.∂ₜq_frz / u.∂ₜn_frz / x_scale, FT(20);
                    rtol = FT == Float64 ? FT(1e-10) : FT(1e-3))

                band = FT[5, 10, 15, 20, 25, 30, 35, 40]
                ev = map(band) do ΔT
                    l = lim(N, T_frz - ΔT)
                    l.∂ₜq_frz / l.∂ₜn_frz / x_scale
                end
                # strictly descending across the band, with no branch or switch producing it
                @test issorted(ev; rev = true)
                @test all(>(0), ev)
                # the two endpoints of the descent, bracketed: near the volume-selective 20 at the
                # warm end, and already below x̄ by 40 K
                @test FT(15) < ev[1] < FT(20)
                @test ev[end] < 1
                # it crosses x̄ inside the band rather than at an endpoint, which is the
                # every-drop-freezing transition the composition is supposed to produce
                @test any(>(1), ev) && any(<(1), ev)

                # ...and past ΔT* it stays well below x̄ instead of returning to 20
                deep = lim(N, T_frz - FT(82))
                @test deep.∂ₜq_frz / deep.∂ₜn_frz / x_scale < FT(0.5)
            end
        end

        @testset "the Koop pathway carries SI units and its own window" begin
            (; Δa_w_min, Δa_w_max) = hom
            @test Δa_w_min < Δa_w_max

            # UNITS. `homogeneous_J_cubic` returns 10^(logJ + 6); the +6 is the cm⁻³ -> m⁻³
            # conversion of Koop's original fit, so the value is ALREADY SI and must not be
            # converted again. These bracket the two window edges evaluated from the shipped
            # coefficients (4.2e2 and 2.9e24 m⁻³ s⁻¹). A stray extra 1e6 lands at 4.2e8 / 2.9e30
            # and fails both; a missing one lands at 4.2e-4 / 2.9e18 and fails both. That is the
            # whole point of hardcoding a magnitude here rather than restating the formula.
            J_lo = CM_HomIce.homogeneous_J_cubic(hom, Δa_w_min)
            J_hi = CM_HomIce.homogeneous_J_cubic(hom, Δa_w_max)
            @test FT(1e2) < J_lo < FT(1e4)
            @test FT(1e23) < J_hi < FT(1e26)

            # monotone increasing in Δa_w across the fitted window
            grid = range(Δa_w_min, Δa_w_max; length = 17)
            @test issorted([CM_HomIce.homogeneous_J_cubic(hom, FT(a)) for a in grid])

            # The wrapper must NEVER throw. `homogeneous_J_cubic` raises a DomainError outside the
            # window, which inside a GPU kernel is an abort of exactly the family this campaign has
            # been chasing, so the wrapper clamps before calling rather than guarding after.
            for T in FT(150):FT(5):FT(300)
                @test isfinite(J_koop_at(T))
                @test J_koop_at(T) ≥ 0
            end
        end

        @testset "the Koop window is entered and clamped where Δa_w says" begin
            (; Δa_w_max) = hom
            # Where the pathway switches on, measured on the code's own saturation vapour
            # pressures rather than assumed from a hand-computed Δa_w(T).
            lo_k, hi_k = FT(1), FT(80)
            for _ in 1:60
                mid = (lo_k + hi_k) / 2
                J_koop_at(T_frz - mid) > 0 ? (hi_k = mid) : (lo_k = mid)
            end
            ΔT_koop_on = hi_k
            # Plausibility, and a real check: the window bounds are on Δa_w, so a misread of
            # Δa_w_min/Δa_w_max moves this, while a units error does not (that is what the
            # magnitude assertions above are for).
            @test ΔT_koop_on > 25
            @test ΔT_koop_on < 36

            # Zero strictly warmer than the window, positive strictly colder.
            @test J_koop_at(T_frz - FT(0.95) * ΔT_koop_on) === zero(FT)
            @test J_koop_at(T_frz - FT(1.05) * ΔT_koop_on) > 0

            # Cold side: CLAMPED at the upper edge, so J stops responding to T entirely. This is
            # inert by construction - see the τ_dend comparison in the next testset.
            J_clamped = CM_HomIce.homogeneous_J_cubic(hom, Δa_w_max)
            for T in (FT(220), FT(210), FT(200), FT(191.15), FT(170))
                @test J_koop_at(T) === J_clamped
            end
        end

        @testset "warm inertness: the Koop pathway is exactly off above the window" begin
            # This exactness claim CAN be made, unlike the τ_dend one, because the warm-side
            # treatment is a hard zero rather than a small number: J_bigg + 0 === J_bigg, so the
            # summed-coefficient nucleation rate is bit-identical to the Bigg-only rate and adding
            # the pathway cannot perturb ordinary mixed-phase behaviour at all.
            for N in Ns, ΔT in FT[5, 10, 15, 20, 25]
                T = T_frz - ΔT
                l = lim(N, T)
                @test l.J_koop === zero(FT)
                @test l.J_bigg + l.J_koop === l.J_bigg
                @test nuc(N, T).∂ₜn_frz === raw(N, T).∂ₜn_frz
                @test nuc(N, T).∂ₜq_frz === raw(N, T).∂ₜq_frz
            end
        end

        @testset "the homogeneous pathway owns the cold end, not Bigg's extrapolation" begin
            # The point of the change: past the handoff the scheme stops depending on Bigg being
            # accidentally right outside its roughly 0 to 40 K fitted range.
            @test J_koop_at(FT(243)) < lim(Ns[1], FT(243)).J_bigg   # Bigg still leads at -30 C
            @test J_koop_at(FT(233)) > lim(Ns[1], FT(233)).J_bigg   # Koop leads by -40 C

            # Both pathways multiply the SAME V_drop through the same PSD moments, so the
            # crossover between them is independent of drop size. Asserted, because a drop-size
            # dependence here would mean the two are not sharing one PSD treatment.
            for N in Ns
                @test (J_koop_at(FT(233)) > lim(N, FT(233)).J_bigg) === true
            end

            # And the homogeneous term ALONE is already fast enough that the conversion stages own
            # τ_eff: Bigg could be deleted at this temperature without changing the answer.
            for N in Ns
                T = FT(236)
                τ_nuc_koop_only = 1 / (J_koop_at(T) * x̄(N) / ρw)
                @test τ_nuc_koop_only < FT(0.01) * (τ_h(N, T) + τ_d(N))
            end
        end

        @testset "the states of record are reduced by τ_dend alone" begin
            # The mint (trigger cell of the deterministic reproducer) and the vertex-like state
            # from `diag/rainfrz_cap_magnitudes.jl` both sit colder than T* = T_frz - ΔT*, so
            # τ_heat is exactly zero on them and the two-stage form passed the raw Bigg rate
            # through untouched. τ_dend is the stage that reaches these states, and this is the
            # falsifier for it: many orders of reduction, from the dendrite time alone.
            for (ρ_s, T_s, q_s, n_s) in (
                (FT(0.83713067), FT(191.148), FT(1.03934624e-4), FT(208.75424)),
                (FT(0.66), FT(191.15), FT(2e-4), FT(5000)),
            )
                @test T_s < T_frz - ΔT_star
                N_s = n_s * ρ_s
                r = CM_HetIce.liquid_freezing_rate(
                    mp.ice.rain_freezing, pdf_r, tps, q_s, ρ_s, N_s, T_s)
                qᵥ_s = TDI.p2q(tps, T_s, ρ_s,
                    TDI.saturation_vapor_pressure_over_liquid(tps, T_s))
                l = CM_HetIce.rain_freezing_rate(
                    mp.ice.rain_freezing, hom, p3.vent, aps, tps, evap, pdf_r, q_s, ρ_s, N_s,
                    T_s, qᵥ_s)

                @test l.τ_heat === zero(FT)              # no heat left to export this cold
                @test l.τ_dend > 0                        # ...so this is the only stage acting
                @test l.τ_nuc + l.τ_heat + l.τ_dend ≥ l.τ_dend

                # Adding the Koop pathway cannot disturb this result, and not by luck. Here
                # (τ_heat + τ_dend)/τ_nuc ≫ 1, so the limited rate tends to f_psd/τ_freeze - the
                # nucleation coefficient CANCELS out of the answer entirely. Whatever J the two
                # pathways sum to, the dendrite stage sets the conversion. Asserted rather than
                # argued, since it is the reason the vertex falsifier survived the change.
                @test (l.τ_heat + l.τ_dend) / l.τ_nuc > FT(1e6)
                @test l.J_koop > 0                        # Koop IS active this cold...
                @test l.J_koop < l.J_bigg                 # ...and still below the extrapolation
                @test all(isfinite, (l.∂ₜn_frz, l.∂ₜq_frz, l.τ_heat, l.τ_dend))
                @test l.∂ₜn_frz > 0
                @test l.∂ₜq_frz > 0
                # at least ten orders of magnitude off the raw Bigg rate, both moments
                @test l.∂ₜn_frz < FT(1e-10) * r.∂ₜn_frz
                @test l.∂ₜq_frz < FT(1e-10) * r.∂ₜq_frz
                # The mean frozen drop mass is NOT preserved here, and that is the point of the
                # per-size form: the mean-drop factor divided both moments by one number and kept
                # the volume-selective 20 x̄ even where every drop is freezing. Past ΔT* the
                # surviving stage is τ_dend = D/v_dend, so the smallest drops convert fastest and
                # the event mass collapses. Asserted as an inequality, not a value: the deep-limit
                # number moment is the one the fixed-node rule under-resolves.
                @test l.∂ₜq_frz / l.∂ₜn_frz < FT(0.5) * r.∂ₜq_frz / r.∂ₜn_frz
            end
        end

        @testset "the vertex state is donor-bounded THROUGH the substep" begin
            # The boundedness does not live in the rate and must not be asserted there: past ΔT*
            # the limited rate is still orders above what the donor can supply over a step, which
            # is correct for a relaxation with τ_eff ≪ dt. It lives in the linearized-implicit
            # update, whose increment for a decay q/τ_eff is -h(q/τ_eff)/(1 + h/τ_eff) - bounded
            # by q at every h. So this drives the PRODUCTION substep entry, rosenbrock_manual()
            # (ManualJacobian + ExplicitGrowthDiagonal + EndStateSaturationAdjustment), at
            # dt = 2 s and nsub = 1, exactly as the box runs it.
            #
            # The mode choice is load-bearing, not incidental. `_species_mask(::ManualJacobian, _)`
            # is `_full_species_mask`, so the EMPTY ice species stays implicit and the damping
            # reaches the ice row. Under `ExactJacobian + ImplicitGrowth` the near-empty mask would
            # route ice (q_ice < 1e-10) to forward Euler instead, where the increment is h ⋅ f and
            # no timescale can bound it. A bound demonstrated on that mode would not transfer.
            #
            # Cloud liquid is set to zero so the ice gain is attributable to rain freezing rather
            # than to the cloud immersion branch. The remaining ice sources at an ice-free
            # supersaturated state - F23 deposition nucleation, and any vapour the step condenses -
            # are MEASURED from the per-process decomposition and added to the budget rather than
            # absorbed into a tolerance.
            manual = BMT.rosenbrock_manual()
            dt = FT(2)
            ρ_v = FT(0.66)
            T_v = FT(191.15)
            q_rai_v = FT(2e-4)
            n_rai_v = FT(5000)
            q_tot_v = FT(0.008295026)   # the trigger cell's, per diag/rainfrz_cap_magnitudes.jl
            x_v = FT[0, 0, q_rai_v, n_rai_v, 0, 0, 0, 0]   # ice-free, cloud-free
            st_v = P3.state_from_prognostic(p3, FT(0), FT(0), FT(0), FT(0))
            logλ_v = P3.get_distribution_logλ(st_v)

            # The state must genuinely be the pathological one, or the bound proves nothing: the
            # UNLIMITED rate has to overshoot the whole cell's water by many orders here.
            r_v = CM_HetIce.liquid_freezing_rate(
                mp.ice.rain_freezing, pdf_r, tps, q_rai_v, ρ_v, n_rai_v * ρ_v, T_v)
            @test r_v.∂ₜq_frz * dt > FT(1e6) * q_tot_v

            t = BMT.bulk_microphysics_tendencies(
                manual, BMT.Microphysics2Moment(), mp, tps,
                ρ_v, T_v, q_tot_v, x_v..., logλ_v, dt, 1,
            )
            rates = SVector{8, FT}(_applied_2mp3(t)...)
            @test all(isfinite, _applied_2mp3(t))
            @test all(isfinite, rates)

            Δ = dt .* rates
            Δq_ice = Δ[5]

            # (1) The conservation bound the mint violated by sixteen orders: ice cannot gain more
            # condensate than the cell's total water.
            @test Δq_ice ≤ q_tot_v * (1 + FT(1e-3))

            # (2) The tighter, rain-attributable bound, against a budget measured not assumed.
            pp_v = _per_process_2mp3(mp, tps, ρ_v, T_v, q_tot_v, x_v..., logλ_v)
            other =
                max(pp_v.ice_deposition.q_ice, zero(FT)) +
                max(pp_v.ice_depsub.q_ice, zero(FT)) +
                max(pp_v.immersion_freezing.q_ice, zero(FT)) +
                max(pp_v.cloud_condevap.q_lcl, zero(FT))
            @test Δq_ice ≤ (q_rai_v + other * dt) * (1 + FT(0.05))

            # (3) Rain is drained, not multiplied, and the state stays non-negative.
            @test Δ[3] ≤ 0
            @test -Δ[3] ≤ q_rai_v * (1 + FT(1e-3))
            @test -Δ[4] ≤ n_rai_v * (1 + FT(1e-3))
            x1 = SVector{8, FT}(x_v...) .+ Δ
            tol = eps(FT) .* (abs.(SVector{8, FT}(x_v...)) .+ abs.(Δ))
            @test all(x1 .≥ -tol)

            # (4) τ_eff ≪ dt here, so the step should convert nearly all of the rain. Stated as a
            # lower bound too: a limiter that quietly switched the process off would also satisfy
            # every upper bound above, and that would be a different defect.
            @test -Δ[3] ≥ FT(0.5) * q_rai_v
        end

        @testset "the per-process sum equals the full entry where the limiter binds" begin
            # Two copies of the rate exist - `_per_process_2mp3` (`psum` below, via
            # `_per_process_2mp3_and_riming`) and the `Microphysics2Moment` entry (`full`, via
            # `Instantaneous2MP3Tendency` -> `_instantaneous_2mp3_tendency` ->
            # `bulk_microphysics_tendencies(Microphysics2Moment(), ...)`) - and they must move
            # together or `f` and the linear post-solve attribution split. Card #23's stage-1
            # commit made `_per_process_2mp3_and_riming`'s condensation/deposition BARE while the
            # `Microphysics2Moment()` dispatch stayed folded (0-FOLDSCOPE's original scoping,
            # shared with the 1M twin) - a real, measured divergence at this exact state for one
            # night's worth of this landing's history, re-derived rather than patched at the time.
            # A later commit (this landing's dispatch-scoped consistency fix) made the
            # `Microphysics2Moment()` dispatch ALSO bare - its OWN 2M-only call sites, no code
            # shared with the 1M dispatch - which makes `full` and `psum` bit-identical again, by
            # the SAME exact-equality invariant this test originally asserted. The intervening
            # reconciliation formula is gone because the thing it reconciled is gone.
            #
            # Evaluated at 25 K of supercooling with millimetre rain, where the limiter is
            # asserted to bind, so the comparison actually exercises the new code on both sides.
            T = T_frz - FT(25)
            x = BMT.MicroState2MP3{FT}(3e-4, 5e7, q_rai, Ns[1] / ρ, 1e-5, 1e4, 0, 0)
            st = P3.state_from_prognostic(p3, ρ * x.q_ice, ρ * x.n_ice, FT(0), FT(0))
            logλ = P3.get_distribution_logλ(st)
            q_tot = FT(0.008)

            pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot, Tuple(x)..., logλ)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
            full = Tuple(g(x))
            psum = Tuple(sum(values(pp)))

            # the limiter must actually be engaged here, or this testset proves nothing
            l = lim(Ns[1], T)
            @test l.τ_heat > 0
            @test l.∂ₜn_frz < FT(0.9) * raw(Ns[1], T).∂ₜn_frz

            @test all(isfinite, full)
            @test all(isfinite, psum)
            # Compared on a scale rather than bit-for-bit: the two assemblies sum in different
            # orders, so only agreement to round-off is claimed.
            scale = maximum(abs, psum) + eps(FT)
            rtol = FT == Float64 ? FT(1e-10) : FT(1e-4)
            @test maximum(abs.(full .- psum)) ≤ rtol * scale
        end
    end
end

test_rain_freezing_is_heat_limited(Float64)
test_rain_freezing_is_heat_limited(Float32)

# Cloud droplets get the SAME freezing composition rain gets, evaluated at cloud sizes.
#
# The principle (discussion 04, Haakon): cloud drops and rain drops are not different substances,
# they are different sizes. Process physics is written as size-dependent functions of the drop
# population, and tendencies differ because sizes differ, never because of the chosen
# categorization. So `cloud_freezing_rate` is `rain_freezing_rate`'s composition with the cloud
# PSD and the Stokes fall speed substituted, and the two share
# `_composed_liquid_freezing_moments` literally.
#
# REPAIRED FROM the campaign version, which validated a TRANSITIONAL design (an INP-budget cap on
# the heterogeneous coefficient, via `immersion_limit_rate`/`τ_act`/`n_active`, competing with an
# uncapped `J_het_uncapped` the return `NamedTuple` carried) that has been further superseded on
# this branch: the current `cloud_freezing_rate` (`src/IceNucleation.jl`) is Bigg alone, with "no
# ice-nucleating-particle budget above it" (its own docstring), and its return carries `J_het`
# (now always the uncapped Bigg coefficient - there is no capped/uncapped distinction left to
# carry) and no `J_cap`/`J_het_uncapped` fields at all. Four of the campaign's testsets tested
# that cap mechanism directly and are DROPPED, not repaired, because there is nothing left in the
# default path for them to test: "warm of the Koop window the whole block is unchanged" (compared
# against the removed `min(Bigg, cap)` block), "the budget applies to the heterogeneous
# coefficient, never the summed rate" and "the sub-235 K exemption emerges structurally, with no
# threshold in the code" (both about the cap's own existence and scaling). One testset is ADAPTED
# rather than dropped: "the kinetic series BINDS once the homogeneous pathway takes over" compared
# `J_koop` against `J_het_uncapped`; since `J_het` IS that same uncapped coefficient now, the
# comparison is kept with `J_het` substituted directly. Every other testset below (the PSD-shape
# and quadrature-exactness checks, the kinetic-series-is-inert check, the Koop-window continuity
# check, the degenerate-state check, the cloud-vs-rain sizing check) never referenced the cap
# mechanism and needed only `cld()`'s call signature repaired (`τ_act`/`n_active` dropped, per the
# current `cloud_freezing_rate(opt, hom, vent, aps, tps, pdf_c, q, ρ, N, T, qᵥ; quad)` signature -
# see the report to the team lead for the full account, including why this went beyond the
# `τ_act` field-access bug alone).
function test_cloud_freezing_is_uniform_with_rain(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    pdf_c = mp.ice.cloud_pdf
    aps = mp.warm_rain.air_properties
    hom = mp.ice.homogeneous
    T_frz = TDI.T_freeze(tps)
    (; ρw) = pdf_c

    ρ = FT(0.9)
    q_lcl = FT(1e-3)
    # Three cloud sizes spanning the range the brief names, 4e-12 to 5e-10 kg.
    Ns = (FT(ρ * q_lcl / 4e-12), FT(ρ * q_lcl / 5e-11), FT(ρ * q_lcl / 5e-10))
    qᵥ_sat_liq(T) = TDI.p2q(tps, T, ρ, TDI.saturation_vapor_pressure_over_liquid(tps, T))

    cld(N, T) = CM_HetIce.cloud_freezing_rate(
        mp.ice.rain_freezing, hom, p3.vent, aps, tps,
        pdf_c, q_lcl, ρ, N, T, qᵥ_sat_liq(T),
    )
    # PSD scale, taken from λc rather than from ρq/N so the SB2006 clamps cannot invalidate it.
    psd(N) = CM2.pdf_cloud_parameters(pdf_c, q_lcl, ρ, N)
    x_scale(N) = FT(π) / 3 * ρw / psd(N).λc     # = (π/6) ρw ⟨D³⟩, the PSD's mean droplet mass

    @testset "cloud freezing is the same composition as rain [FT=$FT]" begin
        @testset "the cloud rule's α is DERIVED from the shipped shape parameters" begin
            # The generalized Gauss-Laguerre weight exponent is α = (νcD + 1)/μcD - 1, which is 1
            # only because νc = μc = 1. If either shape parameter moves in ClimaParams the rule
            # silently stops integrating this PSD, so the coupling is asserted rather than
            # commented. This is the test that fails first if the table is ever mismatched.
            (; νcD, μcD) = psd(Ns[1])
            @test νcD == 5
            @test μcD == 3
            @test isapprox((νcD + 1) / μcD - 1, FT(1); rtol = FT(1e-6))

            # ...and the table itself, against its defining exactness: for weight t^α exp(-t),
            # ∫ t^α exp(-t) t^k dt = Γ(α + 1 + k), which at α = 1 is (k+1)!.
            (; nodes, weights) = CM_HetIce.cloud_freezing_quadrature()
            @test length(nodes) == CM_HetIce.CLOUD_FREEZING_QUADRATURE_ORDER
            @test length(weights) == length(nodes)
            nq = length(nodes)
            @test isapprox(sum(weights), 1.0; rtol = 1e-12)   # Γ(2) = 1: the weight normalizes
            for k in 0:(2 * nq - 1)
                approx = sum(w * u^k for (u, w) in zip(nodes, weights))
                @test isapprox(approx, factorial(big(k + 1)); rtol = 1e-10)
            end
        end

        @testset "the quadrature is exact in the zero-limiting limit" begin
            # The sibling of the rain rule's exactness identity, on the cloud PSD, and the one
            # assertion nothing else in this file can substitute for.
            rtol = FT == Float64 ? FT(1e-10) : FT(1e-4)
            (; nodes, weights) = CM_HetIce.cloud_freezing_quadrature()
            for N in Ns
                (; λc, νcD, μcD) = psd(N)
                n = N / ρ
                D_at = t -> (t / λc)^(1 / μcD)     # the composition's own node-to-diameter map
                for k in (3, 6)
                    quadmom = n * sum(FT(w) * D_at(FT(t))^k for (t, w) in zip(nodes, weights))
                    @test isapprox(quadmom, DT.generalized_gamma_Mⁿ(νcD, μcD, λc, n, k); rtol)
                end
            end
        end

        @testset "the kinetic series is inert where the heterogeneous pathway drives freezing" begin
            # Discussion 04's expectation, measured rather than assumed: at cloud sizes τ_heat is
            # milliseconds and τ_nuc is seconds or longer, so uniformity costs nothing warm of the
            # Koop window. Stated as the retained fraction the composition actually applies.
            for N in Ns, ΔT in FT[5, 10, 15, 20, 25, 30]
                T = T_frz - ΔT
                l = cld(N, T)
                @test l.J_koop === zero(FT)             # the window has not opened yet
                retained = l.τ_nuc / (l.τ_nuc + l.τ_heat + l.τ_dend)
                @test retained > FT(0.995)
                @test l.τ_heat < FT(1.5)
                @test l.τ_dend < FT(1e-3)
            end
            # tighter at the small end, where the whole series is four orders below τ_nuc
            for ΔT in FT[5, 10, 15, 20, 25, 30]
                l = cld(Ns[1], T_frz - ΔT)
                @test l.τ_nuc / (l.τ_nuc + l.τ_heat + l.τ_dend) > FT(0.9999)
            end
        end

        @testset "the kinetic series BINDS once the homogeneous pathway takes over" begin
            # The half of discussion 04's expectation that measurement contradicts, kept as an
            # assertion because it is a real and reportable consequence of uniformity rather than
            # a defect. `J_het` is the uncapped Bigg coefficient in the current source (see the
            # file header), which is exactly what this comparison needs.
            T = FT(233)
            for N in Ns
                l = cld(N, T)
                @test l.J_koop > 0
                @test l.J_koop > FT(1e6) * l.J_het   # homogeneous dominates Bigg by orders
                retained = l.τ_nuc / (l.τ_nuc + l.τ_heat + l.τ_dend)
                @test retained < FT(0.01)
                @test l.τ_nuc < l.τ_heat                       # nucleation is no longer the bottleneck
            end
            # ...and it is still inert THROUGH a 2 s step, which is why this is a rate statement
            # and not a stability one.
            for N in Ns
                l = cld(N, T)
                hoverτ = FT(2) / (l.τ_heat + l.τ_dend)
                @test hoverτ / (1 + hoverτ) > FT(0.9)
            end
        end

        @testset "the composition is continuous across the Koop window opening" begin
            # The window is a hard zero on the warm side, so the SUM has a kink but no jump.
            N = Ns[1]
            ΔTs = FT.(range(28, 34; length = 25))
            rates = [cld(N, T_frz - ΔT).∂ₜn_frz for ΔT in ΔTs]
            @test all(isfinite, rates)
            @test issorted(rates)                          # monotone increasing with supercooling
            # no step larger than a factor of a few between adjacent samples of a 0.25 K grid
            for i in 2:length(rates)
                @test rates[i] ≤ FT(20) * max(rates[i - 1], floatmin(FT))
            end
        end

        @testset "degenerate states stay exactly zero and finite" begin
            for (q_s, n_s) in ((FT(0), FT(0)), (FT(1e-3), FT(0)), (FT(0), FT(1e8)))
                l = CM_HetIce.cloud_freezing_rate(
                    mp.ice.rain_freezing, hom, p3.vent, aps, tps,
                    pdf_c, q_s, ρ, n_s, FT(250), qᵥ_sat_liq(FT(250)),
                )
                @test l.∂ₜn_frz === zero(FT)
                @test l.∂ₜq_frz === zero(FT)
                @test isfinite(l.J_het)
                @test isfinite(l.J_koop)
            end
            # at and above freezing the rate is exactly zero
            for T in (FT(280), T_frz)
                l = cld(Ns[1], T)
                @test l.∂ₜn_frz === zero(FT)
                @test l.∂ₜq_frz === zero(FT)
            end
            # 2 K of supercooling: freezing is available, paired, and small
            let l = cld(Ns[1], T_frz - FT(2))
                @test l.∂ₜn_frz > 0
                @test l.∂ₜq_frz > 0          # paired: no number without the mass to carry it
                @test l.∂ₜn_frz < FT(1e-6)   # and small, as the rolloff toward ΔT = 0 requires
            end
        end

        @testset "rain and cloud go through one composition, differing only by size" begin
            # The principle, asserted structurally: the mean converted mass per event scales with
            # the population's own size, and the cloud answer sits far below the rain answer at
            # the same temperature purely because the droplets are smaller. Nothing in the code
            # branches on the category.
            T = T_frz - FT(20)
            for N in Ns
                l = cld(N, T)
                @test l.∂ₜn_frz > 0
                @test l.∂ₜq_frz > 0
                # event mass in units of the PSD's own mean droplet mass: volume-selective, so
                # above 1 and below the exponential PSD's 20, since the cloud gamma is narrower
                ev = l.∂ₜq_frz / l.∂ₜn_frz / x_scale(N)
                @test 1 < ev < 20
            end
            # larger droplets freeze more mass per event, at the same temperature
            evs = map(N -> cld(N, T).∂ₜq_frz / cld(N, T).∂ₜn_frz, Ns)
            @test issorted(evs)   # Ns is ordered smallest droplet first
        end
    end
end

test_cloud_freezing_is_uniform_with_rain(Float64)
test_cloud_freezing_is_uniform_with_rain(Float32)
# The liquid twin of the testset above, and the same claim: cloud condensation must be exactly zero
# when there is no droplet population to condense onto.
#
# `cloud_condensation_timescale` caps its return at `CLOUD_COND_TIMESCALE_MAX` as the capacitance
# integral vanishes, so `sat_excess / τ` stays finite on a state with no droplets. The degeneracy is
# reached through a switch rather than a limit: `log_pdf_cloud_parameters_mass` returns
# `logA = -Inf` for `q_lcl < eps(FT)`, the diameter moment underflows to zero, and the quotient goes
# to the cap exactly. The rate that leaks is a pure mass source - the number slot is zero by
# construction ("number neglected") - so it manufactures liquid mass-without-number states.
#
# What is asserted is exactly what the gate guarantees:
#   - at the droplet-free state the whole `cloud_condevap` slot is zero while the cell is
#     supersaturated over liquid, i.e. on the growth branch, where the leak lives;
#   - the same rate evaluated at the capped timescale is NOT zero, so the zero is the gate's doing;
#   - the cloud-mass row of both Jacobians is zero there, including the temperature column;
#   - the entry and the per-process decomposition return the same cloud mass there (`psum == full`);
#   - a populated droplet state is NOT degenerate and still condenses, which is the inertness half.
# Bit-identity against a pre-fix baseline is not asserted; this suite has no such baseline.
function test_degenerate_cloud_condensation_is_gated(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    sb = mp.warm_rain.seifert_beheng
    aps = mp.warm_rain.air_properties
    ρ = FT(0.6)
    q_tot = FT(4e-3)
    # cold enough that q_tot exceeds saturation over liquid at every one of them by at least a
    # factor of two, which the setup below asserts rather than assumes; 263 K is deliberately
    # excluded, where q_sat_liq ≈ 3.9e-3 leaves no margin against q_tot
    Ts = FT[220, 233, 245, 253]

    τ_cond(T, q_lcl, n_lcl) = CM2.cloud_condensation_timescale(
        sb.pdf_c, aps, tps, T, ρ, q_lcl,
        ρ * CM2.number_bounded_by_mass_limits(
            (; x_min = sb.pdf_c.xc_min, x_max = sb.pdf_c.xc_max), q_lcl, n_lcl))

    @testset "degenerate cloud condensation is gated [FT=$FT]" begin
        empty_state = P3.state_from_prognostic(p3, FT(0), FT(0), FT(0), FT(0))
        logλ₀ = P3.get_distribution_logλ(empty_state)

        @testset "the droplet-free state condenses and evaporates nothing" begin
            for T in Ts
                qᵥ_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
                @test q_tot > qᵥ_sat_liq   # on the growth branch, which is where the leak is
                τ = τ_cond(T, FT(0), FT(0))
                @test CM2.cloud_condensation_is_degenerate(τ)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
                    FT(0), FT(0), FT(0), FT(0),   # DROPLET-FREE, rain-free
                    FT(0), FT(0), FT(0), FT(0),   # ice-free
                    logλ₀)
                @test all(iszero, Tuple(pp.cloud_condevap))

                # The zero is the gate's, not the arithmetic's: the same conversion evaluated at
                # the capped timescale mints cloud liquid at `sat_excess / (τ_max·Γₗ)`.
                # Through the shared inner function with the timescale this test computed.
                # `CloudLiquidFormation` is a zero-field marker under the process-function
                # convention and the timescale now reaches the conversion through `mp`, so a
                # test that wants a SPECIFIC timescale calls what the public method forwards
                # to rather than inventing a parameter set to carry it.
                leak = CMNonEq._conv_q_vap_to_q_lcl_const(
                    τ, tps,
                    (; q_tot, q_lcl = FT(0), q_icl = FT(0), q_rai = FT(0), q_sno = FT(0)),
                    (; ρ, T))
                @test leak > 0
            end
        end

        # Both Jacobians build their own τ_l from their own clamped state, so neither reads the
        # gated rate. At the all-zero state every other contribution to the cloud-mass row vanishes
        # with its donor, but droplet activation does NOT: it is a nucleation source, it fires on
        # a droplet-free supersaturated cell by construction, and its paired mass writes those
        # three columns through the vapour brake. So what identifies "no condensation derivative
        # here" is no longer a zero - it is that the row is EXACTLY the activation pairing, which
        # is a stronger assertion because it also pins the pairing in the linearization.
        @testset "neither Jacobian carries a condensation derivative at the droplet-free state" begin
            x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, 0, 0, 0, 0)
            x_seed = FT(CM2.activation_droplet_mass(mp.warm_rain.seifert_beheng.pdf_c))
            for T in Ts
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ₀)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot, Tuple(x)..., logλ₀)
                J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
                # row 1 is the cloud-mass row, row 2 the droplet-number row; columns 1/3/5 are
                # their liquid, rain and ice donors
                @test J[1, 1] === x_seed * J[2, 1]
                @test J[1, 3] === x_seed * J[2, 3]
                @test J[1, 5] === x_seed * J[2, 5]

                y = BMT.MicroState2MP3T{FT}(Tuple(x)..., T)
                gT = BMT.Temperature2MP3Tendency(mp, tps, ρ, q_tot, logλ₀)
                ctx = BMT._phase_relaxation_context(mp, tps,
                    _micro_2mp3(q_tot, Tuple(x)...),
                    _thermo_2mp3(ρ, T, logλ₀))
                @test CM2.cloud_condensation_is_degenerate(ctx.τ_l)
                J9 = BMT._jacobian_2mp3t_manual(gT, y, pp, rs, ctx)
                @test J9[1, 1] === x_seed * J9[2, 1]
                @test J9[1, 3] === x_seed * J9[2, 3]
                @test J9[1, 5] === x_seed * J9[2, 5]
                @test J9[1, 9] == 0
                # the temperature tendency's own condensation correction is `(Γₗ − 1)·cloud_condevap`,
                # so it is already zero wherever the slot is gated
                @test all(iszero, Tuple((ctx.Γₗ - 1) * pp.cloud_condevap))
            end
        end

        @testset "the entry and the per-process decomposition agree at the droplet-free state" begin
            x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, 0, 0, 0, 0)
            for T in Ts
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ₀)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot, Tuple(x)..., logλ₀)
                @test g(x).q_lcl == sum(values(pp)).q_lcl
            end
        end

        # A populated droplet state must NOT be degenerate, or the gate would silence real
        # condensation. The predicate is false there, so the `ifelse` selects the same branch as
        # before the gate and the arithmetic is untouched by construction.
        @testset "a populated droplet state is not degenerate and still condenses" begin
            q_lcl, n_lcl = FT(1e-4), FT(1e8)
            for T in FT[233, 245, 253]
                τ = τ_cond(T, q_lcl, n_lcl)
                @test isfinite(τ) && τ > 0
                @test !CM2.cloud_condensation_is_degenerate(τ)
                @test τ < CM2.CLOUD_COND_TIMESCALE_MAX(FT)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
                    q_lcl, n_lcl, FT(0), FT(0),
                    FT(0), FT(0), FT(0), FT(0),
                    logλ₀)
                @test pp.cloud_condevap.q_lcl > 0
                @test isfinite(pp.cloud_condevap.q_lcl)

                x = BMT.MicroState2MP3{FT}(q_lcl, n_lcl, 0, 0, 0, 0, 0, 0)
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ₀)
                J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
                @test isfinite(J[1, 1])
                # card #26's own closed form, at Γ=1 (card #23) with dcp_dliq=dcp_dice=0: the
                # Γ-coupling terms in `_condevap_derivs` vanish identically, leaving
                # `∂s_liq = -1/τ + rate·(1/(3 q_lcl)) = (1/τ)·(sat_excess/(3 q_lcl) - 1)` exactly
                # (τ > 0 always) - so `sat_excess/(3 q_lcl) > 1` is not the design note's own
                # leading-order APPROXIMATION anymore (that estimate's own error came from a
                # Γ-coupling term this call site no longer carries at all), it is the EXACT sign
                # criterion for the isolated condensation entry. `immersion_freezing` (freezing) is a
                # SEPARATE process contributing to the same `lcl_lcl` slot; where it is negligible
                # (`pp.immersion_freezing.q_lcl` below `qmin`-scale) `J[1,1]` reduces to the isolated
                # entry and the criterion governs the FULL entry too - true at T=245/253 here
                # (ratio 10.0/7.1, both `>1`, both precisions, `|J[1,1]|` 0.5-0.8, comfortably off
                # zero). This is a regression guard for card #26 itself: if the `dlog_τ_dq_liq`
                # term is ever dropped again, `J[1,1]` reverts to `-1/τ < 0` at these states and
                # this assertion catches it, where an `isfinite`-only check would not.
                bigg_negligible = T != FT(233)   # see the T=233 branch below
                if bigg_negligible
                    sat_excess_l =
                        TDI.q_vap(q_tot, q_lcl, FT(0)) -
                        TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
                    crit = sat_excess_l / (3 * q_lcl)
                    @test sign(J[1, 1]) == (crit > 1 ? 1 : -1)
                else
                    # T=233 K: `pp.immersion_freezing.q_lcl` is real here (immersion freezing is cold-
                    # temperature-active, not negligible like at 245/253) and DOMINATES `lcl_lcl`
                    # (measured: J[1,1] = -789.76 against ExactJacobian's -288.03, a ~2.7x gap) -
                    # card #16's own already-documented bigg-immersion donor-diagonal Jacobian
                    # understatement, pre-existing and untouched by cards #23/#26. `J[1,1]` stays
                    # strongly negative regardless of what the condensation entry alone would do,
                    # so the sign assertion here tests bigg immersion's sign, not condensation's.
                    @test pp.immersion_freezing.q_lcl < 0   # freezing genuinely active at this T
                    @test J[1, 1] < 0
                end
                y = BMT.MicroState2MP3T{FT}(Tuple(x)..., T)
                gT = BMT.Temperature2MP3Tendency(mp, tps, ρ, q_tot, logλ₀)
                ctx = BMT._phase_relaxation_context(mp, tps,
                    _micro_2mp3(q_tot, Tuple(x)...),
                    _thermo_2mp3(ρ, T, logλ₀))
                @test !CM2.cloud_condensation_is_degenerate(ctx.τ_l)
                J9 = BMT._jacobian_2mp3t_manual(gT, y, pp, rs, ctx)
                # J9's species block IS J8 now (card #23's dispatch-scoped consistency fix), so
                # J9[1,1] carries the identical state-dependent sign for the identical reason.
                @test J9[1, 1] == J[1, 1] && isfinite(J9[1, 1])
                # temperature column: −∂q_sat_liq/∂T / τ_l, the saturation shift term alone - this
                # has no Γ or τ(q) dependence (it is not part of `_condevap_derivs` at all, see
                # `_jacobian_2mp3t_manual`'s `t1`), so its sign is fixed by `∂q_sat_liq/∂T > 0`
                # regardless of the entry above, unaffected by either card.
                @test J9[1, 9] < 0 && isfinite(J9[1, 9])
            end
        end

        # The trigger is the droplet NUMBER, and it is precision independent. It used to be
        # `q_lcl < eps(FT)`, which gated a real droplet population out of condensation whenever
        # its mass content was small - 1.19e-7 kg/kg at Float32 is a loading the box visits with
        # 10^6 droplets per cubic metre present, 2.2e-16 at Float64 is not, so one state had two
        # different physics in the two precisions. Recorded here as the replacement assertion.
        @testset "the trigger is the droplet number, not the mass, at both precisions" begin
            q_trace, n_trace = FT(1e-9), FT(1e3)
            for T in FT[233, 253]
                τ = τ_cond(T, q_trace, n_trace)
                @test !CM2.cloud_condensation_is_degenerate(τ)
                @test isfinite(τ) && τ > 0
                pp = _per_process_2mp3(mp, tps, ρ, T, q_tot,
                    q_trace, n_trace, FT(0), FT(0),
                    FT(0), FT(0), FT(0), FT(0),
                    logλ₀)
                @test pp.cloud_condevap.q_lcl > 0
            end
        end
    end
end

test_degenerate_cloud_condensation_is_gated(Float64)
test_degenerate_cloud_condensation_is_gated(Float32)

# The other half of the same gate: a droplet population that exists must be able to grow.
#
# The states below are the ones the box actually occupied for more than a day - droplets present
# in numbers up to 7e6 per m³, air 20 to 26 percent supersaturated over liquid, cloud liquid mass
# content 0.0 - and the closure was blind to all of them because its degeneracy switch keyed on
# the mass. `τ_cond` sat at `CLOUD_COND_TIMESCALE_MAX` exactly at every droplet number from 1e4 to
# 7e6 per m³, so the rate was either the 1e10 leak or, once gated, exactly zero. Number was
# supplied and mass never followed; the number adjustment then drained the unsupported number, and
# the cold start could not bootstrap a liquid phase at all.
#
# What is asserted:
#   - `τ_cond` is finite and far below the cap at every probe state, and the slot condenses;
#   - the timescale falls with droplet number, i.e. the population is being seen rather than
#     merely admitted;
#   - the Float32 cliff case bootstraps: 7e6 droplets per m³ carrying a 1 μm seed each is
#     q = 2.6e-8 kg/kg, below `eps(Float32)`, and used to be gated at single precision only;
#   - the genuinely empty state is untouched (its own testset above), and a populated ordinary
#     state is unchanged, which `test_populated_cloud_psd_is_unchanged` pins bit-for-bit.
function test_massfree_droplet_population_condenses(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    sb = mp.warm_rain.seifert_beheng
    aps = mp.warm_rain.air_properties
    ρ = FT(0.6)
    T = FT(290)

    # `sat_excess` is what the entry supplies from the state, and it is what decides the zero-mass
    # arm: a droplet number with no mass yet is a real population while the air is supersaturated.
    # These states all are, so it is passed positive here for the same reason the entry passes it.
    τ_cond(q_lcl, n_lcl, sat_excess = one(FT)) = CM2.cloud_condensation_timescale(
        sb.pdf_c, aps, tps, T, ρ, q_lcl,
        ρ * CM2.number_bounded_by_mass_limits(
            (; x_min = sb.pdf_c.xc_min, x_max = sb.pdf_c.xc_max), q_lcl, n_lcl, sat_excess))

    @testset "a mass-free droplet population condenses [FT=$FT]" begin
        empty_state = P3.state_from_prognostic(p3, FT(0), FT(0), FT(0), FT(0))
        logλ₀ = P3.get_distribution_logλ(empty_state)
        qᵥ_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)

        # the probe's own droplet numbers, in per m³, carried as specific numbers
        Ns = FT[1e4, 1e6, 7e6]
        τs = FT[]
        for S in FT[1.20, 1.26]
            q_tot = S * qᵥ_sat_liq
            for N in Ns
                n_lcl = N / ρ
                τ = τ_cond(FT(0), n_lcl)
                @test isfinite(τ) && τ > 0
                @test !CM2.cloud_condensation_is_degenerate(τ)
                @test τ < CM2.CLOUD_COND_TIMESCALE_MAX(FT) / 1000
                pp = _per_process_2mp3(mp, tps, ρ, T, q_tot,
                    FT(0), n_lcl, FT(0), FT(0),
                    FT(0), FT(0), FT(0), FT(0),
                    logλ₀)
                @test pp.cloud_condevap.q_lcl > 0
                @test isfinite(pp.cloud_condevap.q_lcl)
                S == FT(1.20) && push!(τs, τ)
            end
        end
        # more droplets, more surface, faster relaxation - the population is seen, not just admitted
        @test issorted(τs, rev = true)

        # The Float32 cliff: a 1 μm seed on the measured 7e6 droplets per m³ is below eps(Float32),
        # so this state condensed at Float64 and was gated at Float32 before the closure change.
        q_seed = FT(sb.pdf_c.xc_min) * FT(7e6) / ρ
        @test q_seed < eps(Float32)
        τ_seed = τ_cond(q_seed, FT(7e6) / ρ)
        @test !CM2.cloud_condensation_is_degenerate(τ_seed)
        pp_seed = _per_process_2mp3(mp, tps, ρ, T, FT(1.26) * qᵥ_sat_liq,
            q_seed, FT(7e6) / ρ, FT(0), FT(0),
            FT(0), FT(0), FT(0), FT(0),
            logλ₀)
        @test pp_seed.cloud_condevap.q_lcl > 0

        # The Jacobian follows the rate: a live relaxation carries a negative self-derivative.
        x = BMT.MicroState2MP3{FT}(0, FT(7e6) / ρ, 0, 0, 0, 0, 0, 0)
        g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, FT(1.26) * qᵥ_sat_liq, logλ₀)
        (pp, rs) = _per_process_2mp3_and_riming(mp, tps, ρ, T, FT(1.26) * qᵥ_sat_liq,
            Tuple(x)..., logλ₀)
        J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
        @test J[1, 1] < 0 && isfinite(J[1, 1])
        @test all(isfinite, J)
    end
end

test_massfree_droplet_population_condenses(Float64)
test_massfree_droplet_population_condenses(Float32)

# The zero-mass droplet arm: a droplet number carrying no mass yet is retained while the air is
# supersaturated and drained only under subsaturation, which is deactivation. The corner is not
# hypothetical - the positivity clamp `max(x + Δx, 0)` writes exact zeros into masses while leaving
# numbers untouched, so the state reaches it on its own, and a census of the convecting record
# found 33.5 percent of populated cells sitting there. Draining on the mass alone emptied every one
# of them before its condensation had a step to run.
function test_zero_mass_droplet_arm(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    sb = mp.warm_rain.seifert_beheng
    aps = mp.warm_rain.air_properties
    pdf_c = sb.pdf_c
    win = (; x_min = pdf_c.xc_min, x_max = pdf_c.xc_max)
    par = (; sb.numadj.τ, win...)
    τ = FT(sb.numadj.τ)
    ρ, T = FT(0.6), FT(290)
    n = FT(7e6) / ρ

    @testset "the zero-mass droplet arm is decided by saturation [FT=$FT]" begin
        wet, dry = FT(1e-4), FT(-1e-4)

        # (a) supersaturated and mass-free: the population is kept, exactly, and nothing relaxes
        @test CM2.number_bounded_by_mass_limits(win, FT(0), n, wet) === n
        @test CM2.number_tendency_from_mass_limits(par, FT(0), n, wet) === zero(FT)

        # (b) subsaturated and mass-free: deactivation, at the full relaxation rate
        @test CM2.number_bounded_by_mass_limits(win, FT(0), n, dry) === zero(FT)
        @test CM2.number_tendency_from_mass_limits(par, FT(0), n, dry) === -n / τ
        # exactly at saturation the excess is zero and the arm drains: the retention test is
        # strict, so a cell in equilibrium with no mass is not held open indefinitely
        @test CM2.number_tendency_from_mass_limits(par, FT(0), n, zero(FT)) === -n / τ

        # (c) the two functions remain one object. This identity is why the bounded form carries
        # the excess at all: rates must see a retained population and must NOT see a draining one.
        for s in FT[-1e-4, 0, 1e-4], q in FT[0, 1e-9, 1e-4], m in FT[0, 1e2, n, 1e12]
            @test CM2.number_tendency_from_mass_limits(par, q, m, s) ===
                  (CM2.number_bounded_by_mass_limits(win, q, m, s) - m) / τ
        end

        # (d) BIT-IDENTICAL wherever there is mass: the argument may not perturb a populated cell,
        # and the default may not perturb a caller that does not pass it.
        for q in FT[1e-9, 1e-6, 1e-4, 1e-3], m in FT[1e2, n, 1e12]
            base_b = CM2.number_bounded_by_mass_limits(win, q, m)
            base_t = CM2.number_tendency_from_mass_limits(par, q, m)
            for s in FT[-1e-4, 0, 1e-4, 1e3]
                @test CM2.number_bounded_by_mass_limits(win, q, m, s) === base_b
                @test CM2.number_tendency_from_mass_limits(par, q, m, s) === base_t
            end
        end

        # (e) the retained population is a population: it has a finite condensation timescale at
        # the activation-size floor, and the entry grows mass onto it
        qᵥ_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
        N_kept = ρ * CM2.number_bounded_by_mass_limits(win, FT(0), n, wet)
        τ_l = CM2.cloud_condensation_timescale(pdf_c, aps, tps, T, ρ, FT(0), N_kept)
        @test isfinite(τ_l) && !CM2.cloud_condensation_is_degenerate(τ_l)

        # (f) and the drained one is NOT seen, which is the reason the bounded form is gated too:
        # the closure evaluated at a positive number and zero mass in subsaturated air returns an
        # EVAPORATION rate, and evaporating mass the cell does not have is a vapor source out of
        # nothing. Asserted through the entry, where the sign is decided from the state.
        logλ₀ = P3.get_distribution_logλ(
            P3.state_from_prognostic(mp.ice.scheme, FT(0), FT(0), FT(0), FT(0)))
        for (label, q_tot, want_positive) in
            (("supersaturated", FT(1.26) * qᵥ_sat_liq, true),
            ("subsaturated", FT(0.5) * qᵥ_sat_liq, false))
            pp = _per_process_2mp3(mp, tps, ρ, T, q_tot,
                FT(0), n, FT(0), FT(0), FT(0), FT(0), FT(0), FT(0), logλ₀)
            if want_positive
                @test pp.cloud_condevap.q_lcl > 0
                @test pp.cloud_numadj.n_lcl == 0        # kept, not relaxed
            else
                @test pp.cloud_condevap.q_lcl == 0      # no evaporation of absent mass
                @test pp.cloud_numadj.n_lcl ≈ -n / τ rtol = sqrt(eps(FT))
            end
        end

        # (g) f/J consistency of the arm, which is the entry that would otherwise be wrong by
        # `-1/τ` against a rate of exactly zero
        @test BMT._numadj_derivs(FT, FT(0), n, win.x_min, win.x_max, τ, wet)[2] === zero(FT)
        @test BMT._numadj_derivs(FT, FT(0), n, win.x_min, win.x_max, τ, dry)[2] === -1 / τ

        # (h) the branch is constant under differentiation: the excess is a state-dependent
        # quantity, so a switch differentiated through it would put the derivative of a step
        # function into the Jacobian. Seeding it must change nothing at all.
        seeded = FD.Dual{Nothing}(wet, one(FT))
        @test CM2.number_tendency_from_mass_limits(par, FT(0), n, seeded) === zero(FT)
        @test FD.value(CM2.number_bounded_by_mass_limits(win, FT(0), n, seeded)) === n

        # (i) rain and ice are untouched: they take the default and drain at zero mass exactly as
        # before. Their doctrine is a separate decision and this commit does not pre-empt it.
        rain = (; sb.numadj.τ, x_min = sb.pdf_r.xr_min, x_max = sb.pdf_r.xr_max)
        @test CM2.number_tendency_from_mass_limits(rain, FT(0), n) === -n / τ
        @test CM2.number_bounded_by_mass_limits(
            (; x_min = sb.pdf_r.xr_min, x_max = sb.pdf_r.xr_max), FT(0), n) === zero(FT)
    end
end

test_zero_mass_droplet_arm(Float64)
test_zero_mass_droplet_arm(Float32)

# The inertness half, asserted bit-for-bit rather than by tolerance: the mean-droplet-mass floor
# cannot touch a state whose mean droplet mass is above it, and `number_bounded_by_mass_limits`
# holds every number-adjusted state there by construction (`N ≤ ρ q / xc_min`).
function test_populated_cloud_psd_is_unchanged(FT)
    sb = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true).warm_rain.seifert_beheng
    pdf_c = sb.pdf_c
    ρ = FT(0.6)

    # the pre-change arithmetic, kept inline so the assertion does not depend on the code it checks
    function reference_logAB(q, N)
        safe_q = max(q, eps(FT))
        safe_N = max(N, eps(FT))
        (; νc, μc, loggamma_z1, loggamma_z2) = pdf_c
        logx̄ = log(ρ * safe_q / safe_N)
        z1 = (νc + 1) / μc
        logB = -μc * (logx̄ + loggamma_z1 - loggamma_z2)
        return (log(μc) + log(safe_N) + z1 * logB - loggamma_z1, logB)
    end

    @testset "the cloud PSD is bit-identical on populated states [FT=$FT]" begin
        for q in FT[1e-5, 1e-4, 1e-3], N in FT[1e7, 1e8, 5e8]
            x̄ = ρ * q / N
            @test x̄ > pdf_c.xc_min           # the floor cannot bind here
            logA, logB = CM2.log_pdf_cloud_parameters_mass(pdf_c, q, ρ, N)
            refA, refB = reference_logAB(q, N)
            @test logA === refA
            @test logB === refB
        end
        # and the genuinely empty state keeps its degenerate answer exactly
        logA, logB = CM2.log_pdf_cloud_parameters_mass(pdf_c, FT(0), ρ, FT(0))
        @test logA == -Inf
        @test logB == Inf
    end
end

test_populated_cloud_psd_is_unchanged(Float64)
test_populated_cloud_psd_is_unchanged(Float32)

# Sublimation must take the crystals with the mass, at every ice loading it removes mass from.
#
# The ice-number sublimation pathway is `∂ₜn_ice_dep = n_ice·(∂ₜq_ice_dep/q_ice)`, and its quotient
# was guarded by `q_ice > ϵₘ` where the only hazard is 0/0. Below the guard the false arm returned
# ZERO rather than declining to sublimate, so the mass drained at the full capacitance rate and
# every crystal stayed behind - a number-without-mass manufacturer, and (because `ϵₘ = eps(FT)`) one
# that fires at Float32 on a physical ice loading and never at Float64.
#
# The invariant asserted is the one the fix is for: under sublimation the two moments lose the same
# FRACTION per second, so the last of the mass takes the last of the number with it. On the
# mass-limited branch - the branch a trace ice loading in a subsaturated cell always takes - that
# fraction is exactly −1/(τ_dep·Γᵢ) for both.
function test_trace_ice_sublimation_removes_number(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    ρ = FT(0.6)
    # 5e-8 kg/kg is below eps(Float32) = 1.19e-7 and above eps(Float64) = 2.2e-16, so before the
    # fix this one state had different phase-change physics in the two precisions
    q_ice, n_ice = FT(5e-8), FT(1e4)
    # subsaturated over ice at every temperature below (the tightest is 233 K, where
    # q_sat_ice ≈ 2.0e-4), asserted in the loop rather than assumed
    q_tot = FT(1e-4)
    Ts = FT[233, 245, 253]

    @testset "trace ice sublimates its number with its mass [FT=$FT]" begin
        @test q_ice < eps(Float32)
        st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, FT(0), FT(0))
        logλ = P3.get_distribution_logλ(st)

        for T in Ts
            qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
            @test q_tot < qᵥ_sat_ice   # subsaturated: the sublimation branch is the live one
            τ = P3.ice_deposition_timescale(
                mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps,
                T, ρ, st, logλ; quad = mp.ice.quad)
            # the slot must be live, or this measures the degeneracy gate instead
            @test !P3.ice_deposition_is_degenerate(τ)

            pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
                FT(0), FT(0), FT(0), FT(0),   # droplet-free, rain-free
                q_ice, n_ice, FT(0), FT(0),
                logλ)
            @test pp.ice_depsub.q_ice < 0
            @test pp.ice_depsub.n_ice < 0   # exactly zero at Float32 before the fix
            # equal fractional loss ⇒ the mean particle mass is untouched by sublimation
            @test pp.ice_depsub.n_ice / n_ice ≈ pp.ice_depsub.q_ice / q_ice rtol = sqrt(eps(FT))

            # the entry forms the same rate a second time, and it is what
            # `Instantaneous2MP3Tendency` evaluates - both now bare (card #23's dispatch-scoped
            # consistency fix made `Microphysics2Moment()`'s dispatch bare too, its own 2M-only
            # call sites, no code shared with the 1M dispatch), so this is bit-identical again.
            x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, q_ice, n_ice, 0, 0)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
            @test g(x).n_ice == sum(values(pp)).n_ice

            # J must run the pathway f runs. Column 1 of the ice-number row was the sublimation
            # pathway's OWN cp-coupling cross-term before card #23 - `_condevap_derivs`'s
            # `∂s_liq` on the ice call's limited branch, `-q_ice·dinv_dliq`, nonzero only through
            # `dinv_dliq`'s dependence on `dcp_dliq` (the Γ-coupling). Card #23 zeroes `dcp_dliq`
            # at this call site (bare rate), so `dinv_dliq = 0` and the entry is now EXACTLY zero -
            # matching what the bare-by-construction 9x9 always computed here (see below). It was
            # nonzero, not zero, at Float32 before card #23; that history stays in the git log,
            # not as a live assertion of a value the fix retires.
            J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
            @test J[6, 1] == 0 && isfinite(J[6, 1])
            @test all(isfinite, J)
            # the self-damping is the fractional mass loss, so it is negative and bounded
            @test J[6, 6] < 0 && isfinite(J[6, 6])

            y = BMT.MicroState2MP3T{FT}(Tuple(x)..., T)
            gT = BMT.Temperature2MP3Tendency(mp, tps, ρ, q_tot, logλ)
            ctx = BMT._phase_relaxation_context(mp, tps,
                _micro_2mp3(q_tot, Tuple(x)...),
                _thermo_2mp3(ρ, T, logλ))
            J9 = BMT._jacobian_2mp3t_manual(gT, y, pp, rs, ctx)
            # J9[6,1] is STILL zero, but for a simpler reason now: `_jacobian_2mp3_manual` is
            # bare-consistent by construction (card #23), so `_jacobian_2mp3t_manual`'s species
            # block reuses it directly with no folded-vs-bare correction left to apply (the
            # `Δ6` cancellation this comment used to describe no longer exists - removed along
            # with the rest of that now-dead machinery, see the landing commit). Both sides are
            # zero because the underlying entry is zero, not because two nonzero numbers cancel.
            @test J9[6, 1] == 0 && isfinite(J9[6, 1])
            # card #23: `_jacobian_2mp3t_manual`'s species block IS `_jacobian_2mp3_manual`'s
            # output now (no per-entry Γ correction remains), so J9[6,6] agrees with J[6,6]
            # exactly rather than through the fractional-loss correction this assertion used to
            # need - verified bit-identical at this exact state before writing this assertion.
            @test J9[6, 6] == J[6, 6]
            @test J9[6, 6] < 0 && isfinite(J9[6, 6])
            @test all(isfinite, J9)
        end

        # Inertness where the guard never bound: a populated ice state took the true quotient
        # before and after, so the proportionality is the same statement there.
        @testset "a populated ice state is unaffected" begin
            q_ice_p, n_ice_p = FT(1e-5), FT(1e4)
            st_p = P3.state_from_prognostic(p3, ρ * q_ice_p, ρ * n_ice_p, FT(0), FT(0))
            logλ_p = P3.get_distribution_logλ(st_p)
            for T in Ts
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
                    FT(0), FT(0), FT(0), FT(0),
                    q_ice_p, n_ice_p, FT(0), FT(0),
                    logλ_p)
                @test pp.ice_depsub.q_ice < 0
                @test pp.ice_depsub.n_ice < 0
                @test pp.ice_depsub.n_ice / n_ice_p ≈ pp.ice_depsub.q_ice / q_ice_p rtol =
                    sqrt(eps(FT))
            end
        end

        # Deposition still leaves the number alone: only the sublimation branch carries number.
        @testset "the deposition branch carries no number" begin
            q_tot_super = FT(4e-3)
            for T in Ts
                qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
                @test q_tot_super > qᵥ_sat_ice
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot_super,
                    FT(0), FT(0), FT(0), FT(0),
                    q_ice, n_ice, FT(0), FT(0),
                    logλ)
                @test pp.ice_depsub.q_ice > 0
                @test pp.ice_depsub.n_ice == 0
            end
        end
    end
end

test_trace_ice_sublimation_removes_number(Float64)
test_trace_ice_sublimation_removes_number(Float32)

# The ice number adjustment's lower mean-mass bound used to be an independent literal, `1e-12` kg,
# annotated "~10 μm crystal" in one of its three copies. It is a 12.772 μm solid-ice sphere, 2.083x
# the mass of the crystal the scheme actually nucleates at (`ρ_i (π/6) (10 μm)³ = 4.7998e-13` kg),
# so a freshly nucleated population read as `n > q / x_min` and the guard relaxed its number toward
# 48% of what nucleation had just supplied, conserving the mass. `x_min` is now DERIVED from the
# nucleation diameter, which is what these tests hold.
#
# What is asserted is exactly what the derivation guarantees:
#   - the three former copies return one value, and it is the nucleation mass;
#   - the adjustment is EXACTLY inert on a population at the nucleation mean mass, which is the
#     property the old literal did not have - and the old literal's non-inertness is asserted too,
#     inline, so a future edit that quietly restores it fails here rather than passing vacuously;
#   - the melt number rate carries no mean-mass floor; the nucleation size enters melting
#     through `ice_melt_fraction_limit`, from the same shared constant.
# What is NOT asserted: that anything about ice onset improves. That is a trajectory claim and it
# belongs to the ensemble, not to a unit test.
function test_ice_numadj_bound_is_the_nucleation_mass(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    numadj = BMT._ice_numadj_params(p3)
    (; m_nuc) = CMP.ice_seed(p3)

    @testset "the ice numadj lower bound is the nucleation mass [FT=$FT]" begin
        @testset "one constant, derived from the nucleation diameter" begin
            @test (2 * CMP.ice_seed(p3).r_nuc) == FT(2e-6)
            @test m_nuc ≈ p3.ρ_i * FT(π) / 6 * (2 * CMP.ice_seed(p3).r_nuc)^3
            @test m_nuc ≈ FT(3.8399e-15) rtol = FT(1e-4)
            # the three sites that used to carry `1e-12` independently
            @test numadj.x_min == m_nuc
            @test P3.ice_mean_particle_mass_min(p3) == m_nuc
            # and the bound is strictly below the old literal, by the factor the nascent
            # diameter sets: a 2 μm starter rather than the 10 μm the literal was chosen against
            @test numadj.x_min < FT(1e-12)
            @test FT(1e-12) / numadj.x_min ≈ FT(260.42) rtol = FT(1e-3)
        end

        # `q_ice` is kept well above `eps(Float32)` so this exercises the x_min CLAMP arm and not
        # the empty arm, which is a separate defect with its own testset.
        q_ice = FT(1e-5)
        n_nuc = q_ice / m_nuc      # a population made entirely of fresh crystals
        @test q_ice > eps(Float32)

        @testset "the adjustment is exactly inert at the nucleation mean mass" begin
            @test CM2.number_tendency_from_mass_limits(numadj, q_ice, n_nuc) == 0
            # one crystal per kg coarser: still interior, still inert
            @test CM2.number_tendency_from_mass_limits(numadj, q_ice, n_nuc / 2) == 0
            # finer than anything the scheme can make: the guard is supposed to act there
            @test CM2.number_tendency_from_mass_limits(numadj, q_ice, 2 * n_nuc) < 0
        end

        @testset "the old literal was NOT inert there" begin
            old = (; τ = numadj.τ, x_min = FT(1e-12), x_max = numadj.x_max)
            ∂ₜn_old = CM2.number_tendency_from_mass_limits(old, q_ice, n_nuc)
            @test ∂ₜn_old < 0
            @test ∂ₜn_old ≈ (q_ice / FT(1e-12) - n_nuc) / numadj.τ
            # the target it relaxed toward, as a fraction of the number just nucleated
            @test (q_ice / FT(1e-12)) / n_nuc ≈ FT(3.8399e-3) rtol = FT(1e-2)
        end

        @testset "the production per-process slot is inert there too" begin
            tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
            ρ = FT(0.6)
            q_tot = FT(4e-3)
            state = P3.state_from_prognostic(p3, q_ice * ρ, n_nuc * ρ, FT(0), FT(0))
            logλ = P3.get_distribution_logλ(state)
            pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, FT(245), q_tot,
                FT(0), FT(0), FT(0), FT(0),          # liquid, rain
                q_ice, n_nuc, FT(0), FT(0),          # fresh pristine ice
                logλ)
            @test pp.ice_numadj.n_ice == 0
            @test all(iszero, Tuple(pp.ice_numadj))
        end

        @testset "the melt number rate carries no mean-mass floor" begin
            # `dNdt = ρn_ice * melt_frac`: the zero-mass state melts no number (previously
            # `dLdt / m_nuc` through the floored mean mass), and the nucleation size enters
            # through `ice_melt_fraction_limit` instead.
            tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
            vel = mp.ice.terminal_velocity
            aps = mp.warm_rain.air_properties
            ρₐ = FT(1.2)
            Nᵢ = FT(2e5) * ρₐ
            state = P3.P3State(p3, FT(1e-4) * ρₐ, Nᵢ, FT(0.8), FT(800))
            logλ = P3.get_distribution_logλ(state)
            state₀ = P3.P3State(p3, FT(0), Nᵢ, FT(0.8), FT(800))
            rate = P3.ice_melt(vel, aps, tps, FT(273.16), ρₐ, state₀, logλ; quad = mp.ice.quad)
            @test rate.melt_frac == 0
            @test rate.dNdt == 0
            @test isfinite(rate.dLdt)
            # the nucleation size still decides the bound, from the same shared constant
            r_nuc = (2 * CMP.ice_seed(p3).r_nuc) / 2
            lim = P3.ice_melt_fraction_limit(aps, tps, p3, FT(273.16))
            L_f = TDI.Lf(tps, FT(273.16))
            ΔT = FT(273.16) - p3.T_freeze
            @test lim.inv_τ ≈ 3 * aps.K_therm * ΔT / (p3.ρ_i * r_nuc^2 * L_f) rtol =
                sqrt(eps(FT))
        end
    end
end

test_ice_numadj_bound_is_the_nucleation_mass(Float64)
test_ice_numadj_bound_is_the_nucleation_mass(Float32)

# The number adjustment's "no mass -> no particles" arm was coded as `q < ϵ_numerics_2M_M(FT)`,
# i.e. `eps(FT)`, which at Float32 is 1.1920929e-7 kg/kg. Over that band the target number is zero,
# so `∂ₜn = -n/τ` runs at full strength - carried implicitly by `_numadj_derivs`, hence
# unconditionally stable at h = 2 s - while the mass is untouched. That is a mass-without-number
# manufacturer on the low-mass side of cloud, rain and ice, live at Float32 only. The arm is now
# `q > 0`, which is what the in-function comment ("when q == 0") and the docstring formula always
# described, and the clamp already produced a zero target at q == 0 on its own.
#
# The change also repairs a documented invariant that was false: `number_bounded_by_mass_limits` is
# advertised as returning the `n_target` of `number_tendency_from_mass_limits`, and inside the band
# it returned `n` while the tendency targeted zero - so the process rates were evaluated at a
# population the adjustment was simultaneously destroying. The equality is now asserted exactly.
#
# What is NOT asserted: that any census fraction moves, or that any run survives. Both need the box.
function test_numadj_presence_test_is_zero_mass(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    sb = mp.warm_rain.seifert_beheng
    numadj_ice = BMT._ice_numadj_params(p3)
    species = (
        ("cloud", (; τ = sb.numadj.τ, x_min = sb.pdf_c.xc_min, x_max = sb.pdf_c.xc_max)),
        ("rain", (; τ = sb.numadj.τ, x_min = sb.pdf_r.xr_min, x_max = sb.pdf_r.xr_max)),
        ("ice", numadj_ice),
    )
    # inside eps(Float32) at both precisions, so the Float32 run exercises the old band and the
    # Float64 run checks that the same code path gives the same answer there
    q_band = FT(5e-8)
    ϵ_f32 = FT(eps(Float32))

    @testset "the numadj presence test is zero mass, not eps [FT=$FT]" begin
        @testset "$name: a trace population keeps its number" for (name, p) in species
            n_interior = q_band / sqrt(p.x_min * p.x_max)   # mean mass inside [x_min, x_max]
            @test p.x_min < q_band / n_interior < p.x_max
            @test CM2.number_tendency_from_mass_limits(p, q_band, n_interior) == 0
            # the old gate at the same state, computed inline so this cannot pass vacuously
            old_target = ifelse(q_band < ϵ_f32, zero(FT), n_interior)
            @test (old_target - n_interior) / p.τ ≈ -n_interior / p.τ
            @test -n_interior / p.τ < 0

            # the bounds still act inside the band - it is a presence test, not a bypass
            @test CM2.number_tendency_from_mass_limits(p, q_band, 4 * q_band / p.x_min) < 0
            @test CM2.number_tendency_from_mass_limits(p, q_band, q_band / p.x_max / 4) > 0
        end

        @testset "$name: zero and negative mass still empty the number" for (name, p) in species
            n = FT(1e4)
            @test CM2.number_tendency_from_mass_limits(p, zero(FT), n) == -n / p.τ
            @test CM2.number_tendency_from_mass_limits(p, -one(FT) * FT(1e-9), n) == -n / p.τ
            @test CM2.number_bounded_by_mass_limits(p, zero(FT), n) == 0
        end

        # the invariant the docstring claims and the band used to break
        @testset "$name: the bounded number IS the tendency's target" for (name, p) in species
            for q in FT[0, 1e-30, 1e-12, 5e-8, 1e-7, 1e-5, 1e-3, 1e-1]
                for n in FT[0, 1e-6, 1, 1e4, 1e10]
                    ∂ₜn = CM2.number_tendency_from_mass_limits(p, q, n)
                    nb = CM2.number_bounded_by_mass_limits(p, q, n)
                    @test ∂ₜn == (nb - n) / p.τ
                    @test isfinite(∂ₜn) && isfinite(nb)
                end
            end
        end

        @testset "$name: f and J agree on the band" for (name, p) in species
            for (q, n) in ((q_band, q_band / sqrt(p.x_min * p.x_max)),
                (q_band, 4 * q_band / p.x_min), (zero(FT), FT(1e4)))
                (∂q, ∂n) = BMT._numadj_derivs(FT, q, n, p.x_min, p.x_max, p.τ)
                @test ∂q == 0
                h = sqrt(eps(FT)) * max(n, one(FT))
                fd =
                    (
                        CM2.number_tendency_from_mass_limits(p, q, n + h) -
                        CM2.number_tendency_from_mass_limits(p, q, n - h)
                    ) / (2h)
                @test isapprox(∂n, fd; rtol = FT(1e-3), atol = FT(1e-6) / p.τ)
            end
        end

        @testset "the production ice slot no longer empties a trace population" begin
            ρ = FT(0.6)
            q_tot = FT(4e-3)
            n_ice = q_band / sqrt(numadj_ice.x_min * numadj_ice.x_max)
            state = P3.state_from_prognostic(p3, q_band * ρ, n_ice * ρ, FT(0), FT(0))
            logλ = P3.get_distribution_logλ(state)
            pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, FT(245), q_tot,
                FT(0), FT(0), FT(0), FT(0),          # liquid, rain
                q_band, n_ice, FT(0), FT(0),         # trace ice, consistent number
                logλ)
            @test pp.ice_numadj.n_ice == 0
            # the adjustment is number-only by construction, so an active arm can only ever remove
            # number and never the mass that goes with it; pin that so a "mass return path" cannot
            # be added here by accident instead of in a paired term
            @test all(iszero, Tuple(pp.ice_numadj))
            pp_empty = _per_process_2mp3(mp, tps, ρ, FT(245), q_tot,
                FT(0), FT(0), FT(0), FT(0),
                FT(0), n_ice, FT(0), FT(0),          # number without ANY mass
                logλ)
            @test pp_empty.ice_numadj.n_ice == -n_ice / numadj_ice.τ
            @test pp_empty.ice_numadj.q_ice == 0
        end
    end
end

test_numadj_presence_test_is_zero_mass(Float64)
test_numadj_presence_test_is_zero_mass(Float32)

# The whole mixed-phase ice block - liquid-ice collisions, aggregation and melting - sat behind
# `if q_ice > ϵ_numerics_2M_M(FT)`. That is `eps(FT)`, 1.1920929e-7 kg/kg at Float32, used as an
# absolute gate on a mass. Below it the melting rate was not small, it was identically zero at every
# temperature including 300 K, and it went to zero discontinuously as `q_ice` crossed the threshold
# from above, mid-melt. The gate is now the presence of both ice moments.
#
# The tests are written so that the Float64 run is the CONTROL for the Float32 one: the states are
# chosen around `eps(Float32)`, where Float64 always ran the block and Float32 never did. Asserting
# the same answers at both precisions is the precision-split repair stated as a test.
#
# What is NOT asserted: that trace ice removal changes IWP, the census, or any survival fraction.
# The audit's own arithmetic puts the IWP bias below a percent; the reason to fix this is the
# discontinuity and the precision split, not the bias.
function test_mixed_phase_gate_is_ice_presence(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    vel = mp.ice.terminal_velocity
    aps = mp.warm_rain.air_properties
    quad = mp.ice.quad
    T_freeze = TDI.T_freeze(tps)
    ρ = FT(0.9)
    q_tot = FT(6e-3)
    T_warm = T_freeze + FT(5)
    ϵ_f32 = FT(eps(Float32))
    x̄ = FT(1e-11)   # mean ice particle mass held fixed across the sweep [kg]

    ppcall(T, q_ice, n_ice, logλ) = _per_process_2mp3(mp, tps, ρ, T, q_tot,
        FT(0), FT(0), FT(0), FT(0), q_ice, n_ice, FT(0), FT(0), logλ)
    mkstate(q_ice, n_ice) = P3.state_from_prognostic(p3, q_ice * ρ, n_ice * ρ, FT(0), FT(0))

    @testset "the mixed-phase gate is ice presence [FT=$FT]" begin
        @testset "the predicate names exactly the degenerate states" begin
            @test P3.ice_population_is_present(mkstate(FT(1e-4), FT(1e-4) / x̄))
            @test P3.ice_population_is_present(mkstate(FT(1e-9), FT(1e-9) / x̄))
            @test !P3.ice_population_is_present(mkstate(FT(0), FT(0)))
            @test !P3.ice_population_is_present(mkstate(FT(5e-8), FT(0)))   # mass, no number
            @test !P3.ice_population_is_present(mkstate(FT(0), FT(1e4)))    # number, no mass
        end

        @testset "presence scales with the nucleation mass" begin
            # The mass conjunct is a quotient, `ρq_ice / m_nuc > ρn_ice`: a
            # population averaging less than one nucleated crystal's mass per particle reads
            # absent - the subnormal-mass, populated-number crash state among them - and the
            # product form is avoided because it underflows Float32 at trace number.
            (; m_nuc) = CMP.ice_seed(p3)
            mkvol(L, N) = P3.state_from_prognostic(p3, L, N, FT(0), FT(0))
            N = FT(1e4)
            @test !P3.ice_population_is_present(mkvol(nextfloat(zero(FT)), FT(0.55)))
            @test !P3.ice_population_is_present(mkvol(FT(0.5) * N * m_nuc, N))
            @test P3.ice_population_is_present(mkvol(2 * N * m_nuc, N))
            # one crystal's mass per particle: decided by the predicate's own rounded
            # quotient, deterministically
            q_bnd = N * m_nuc
            @test P3.ice_population_is_present(mkvol(q_bnd, N)) == (q_bnd / m_nuc > N)
            # the mass conjunct alone reads positive at any mass; the number conjunct keeps
            # mass-without-number absent
            @test !P3.ice_population_is_present(mkvol(FT(5e-8), FT(0)))
            # the trace-number regime where the product form underflows Float32 stays decided
            # by the mass
            @test P3.ice_population_is_present(mkvol(FT(1e-20), FT(1e-33)))
        end

        @testset "trace ice in warm air melts, at the rate ice_melt gives" begin
            q_ice = ϵ_f32 / 2
            n_ice = q_ice / x̄
            state = mkstate(q_ice, n_ice)
            logλ = P3.get_distribution_logλ(state)
            melt = P3.ice_melt(vel, aps, tps, T_warm, ρ, state, logλ; quad)
            @test melt.dLdt > 0            # the rate was always there; the gate hid it
            pp = ppcall(T_warm, q_ice, n_ice, logλ)
            @test pp.ice_melting.q_ice == -melt.dLdt / ρ
            @test pp.ice_melting.q_rai == melt.dLdt / ρ
            @test pp.ice_melting.n_ice == -melt.dNdt / ρ
            # the removal timescale is physical, not a trickle: the trace population is gone in
            # well under an hour of warm air rather than never
            @test 0 < q_ice / (melt.dLdt / ρ) < 3600
        end

        @testset "melting is continuous across the old threshold" begin
            rates = map((FT(0.5), FT(0.99), FT(1.01), FT(2))) do f
                q_ice = f * ϵ_f32
                state = mkstate(q_ice, q_ice / x̄)
                logλ = P3.get_distribution_logλ(state)
                pp = ppcall(T_warm, q_ice, q_ice / x̄, logλ)
                (f, -pp.ice_melting.q_ice)
            end
            @test all(r -> r[2] > 0, rates)
            # the loading doubles from 0.5 to 0.99 and again from 1.01 to 2, so the rate rises
            @test rates[1][2] < rates[2][2]
            @test rates[3][2] < rates[4][2]
            # and the step across the OLD threshold, 0.99 -> 1.01, is the 2% the loading changed
            # by rather than the 100% cliff the epsilon gate produced. Asserted on |jump| and not
            # on monotonicity: the two are 2% apart and the F32 shape solve's own residual is a
            # few 1e-3, so a monotonicity test on that pair would be tighter than the solver.
            jump = (rates[3][2] - rates[2][2]) / rates[2][2]
            @test abs(jump) < FT(0.1)
        end

        @testset "healthy loadings are untouched" begin
            for (q_ice, T) in ((FT(1e-4), T_warm), (FT(1e-4), T_freeze - FT(10)),
                (FT(1e-3), T_freeze - FT(20)))
                n_ice = q_ice / x̄
                state = mkstate(q_ice, n_ice)
                @test P3.ice_population_is_present(state)
                logλ = P3.get_distribution_logλ(state)
                pp = ppcall(T, q_ice, n_ice, logλ)
                @test all(isfinite, Tuple(pp.ice_melting))
                @test all(isfinite, Tuple(pp.ice_aggregation))
                @test all(isfinite, Tuple(pp.liquid_ice_collision))
                @test pp.ice_aggregation.n_ice <= 0
            end
        end

        @testset "a degenerate ice state still gets exactly zero, and stays finite" begin
            for (q_ice, n_ice) in ((FT(0), FT(0)), (FT(5e-8), FT(0)), (FT(0), FT(1e4)))
                logλ = P3.get_distribution_logλ(mkstate(q_ice, n_ice))
                for T in (T_warm, T_freeze - FT(15))
                    pp = ppcall(T, q_ice, n_ice, logλ)
                    @test all(iszero, Tuple(pp.ice_melting))
                    @test all(iszero, Tuple(pp.ice_aggregation))
                    @test all(iszero, Tuple(pp.liquid_ice_collision))
                end
            end
        end

        @testset "the crash-corpus killing and trace states get exactly zero melt" begin
            # The round-5 killing state's own recorded values (a different air density than
            # this function's shared `ρ`, so it carries its own local one) and the archived
            # trace state, with the crash-corpus's own documented `n_ice` companion (its exact
            # `n_ice` was not recorded; 1.0 #/kg specific matches
            # `notes/corpus-test-design.md`). Both temperatures are above `T_freeze`, so melt
            # would otherwise fire. `pp.ice_melting` isolates melt's own contribution from the
            # unrelated orphan-mass-drain and number-adjustment processes, which legitimately
            # act on these same low-mass states (see `test_ice_orphan_doctrine` below) and are
            # out of this claim's scope.
            for (label, T, ρ_state, ρq_ice, ρn_ice) in (
                ("killing_state", FT(291.95975), FT(1.0475167), FT(1.0e-45), FT(0.55014247)),
                ("trace_state", FT(294.5), FT(1.0475167),
                    FT(5.5e-36) * FT(1.0475167), FT(1.0) * FT(1.0475167)),
            )
                @testset "$label" begin
                    state = P3.state_from_prognostic(p3, ρq_ice, ρn_ice, FT(0), FT(0))
                    @test !P3.ice_population_is_present(state)
                    logλ = P3.get_distribution_logλ(state)
                    pp = _per_process_2mp3(mp, tps, ρ_state, T, q_tot,
                        FT(0), FT(0), FT(0), FT(0), ρq_ice / ρ_state, ρn_ice / ρ_state,
                        FT(0), FT(0), logλ)
                    @test all(iszero, Tuple(pp.ice_melting))
                end
            end
        end

        # the quadrature now runs on states it never saw: the shape solve floors both moments with
        # the same epsilon, so a small enough loading pins logλ near its bracket floor and the
        # integration bounds stretch out. Nothing may become non-finite there. The
        # `q_ice / 1e-13` rows sit below the nucleation mean mass, so the mixed-phase trio is
        # exactly zero for them and the finiteness assertion covers the remaining processes.
        @testset "the newly admitted states stay finite" begin
            for q_ice in FT[1e-7, 1e-9, 1e-12, 1e-15, 1e-20]
                for n_ice in (q_ice / x̄, q_ice / FT(1e-5), q_ice / FT(1e-13))
                    logλ = P3.get_distribution_logλ(mkstate(q_ice, n_ice))
                    @test isfinite(logλ)
                    for T in (T_warm, T_freeze - FT(15))
                        pp = ppcall(T, q_ice, n_ice, logλ)
                        @test all(v -> all(isfinite, Tuple(v)), values(pp))
                    end
                end
            end
        end

        # one copy gated and the other not is the f/J split the deposition gate hit; the entry and
        # the per-process decomposition have to agree at the states the gate newly admits
        @testset "the entry and the per-process decomposition agree on the trace states" begin
            for q_ice in FT[1e-8, 5e-8, 1e-7]
                n_ice = q_ice / x̄
                logλ = P3.get_distribution_logλ(mkstate(q_ice, n_ice))
                x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, q_ice, n_ice, 0, 0)
                for T in (T_warm, T_freeze - FT(15))
                    g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
                    psum = sum(values(ppcall(T, q_ice, n_ice, logλ)))
                    @test g(x).q_ice ≈ psum.q_ice rtol = sqrt(eps(FT))
                    @test g(x).n_ice ≈ psum.n_ice rtol = sqrt(eps(FT))
                    @test g(x).q_rai ≈ psum.q_rai rtol = sqrt(eps(FT))
                end
            end
        end
    end
end

test_mixed_phase_gate_is_ice_presence(Float64)
test_mixed_phase_gate_is_ice_presence(Float32)

# Riming is the one condensate-to-condensate transfer whose manual-Jacobian column did not
# sum to zero: the block carried the cloud and rain donor sinks and no ice receiver, so a
# linearization of an exactly conservative primal created condensate inside the solve.
# Autoconversion, accretion, Bigg immersion, rain freezing and melt all carry both sides.
#
# What is asserted, and at which strength:
#   - `_riming_jacobian_block` balances EXACTLY. The receivers are built as the negatives of
#     the sink terms, so `receivers == -sink` is a bitwise identity and is tested as one.
#   - the assembled 8x8 balances to ACCUMULATION ROUNDING only. The block's entries are added
#     into accumulators that already hold other processes, so the column sum of the assembled
#     matrix cannot be bitwise zero and is bounded against the column's own scale instead.
#   - the block touches exactly TEN entries (card #17 added the tenth: `_riming_jacobian_block`
#     itself still only ever writes nine, but wet growth's own self-term (`BIWET`, restored on
#     `brim_brim`/`J[8,8]` directly rather than through the block, since it is a relaxation
#     toward a fixed density rather than a donor transfer) is zeroed in `rs_zero` alongside the
#     donor split, so it shows up as a tenth entry that moves between `J` and `J0` exactly like
#     the other nine). That IS bitwise: every other accumulator is the same expression whether
#     the split is passed or zeroed.
function test_riming_jacobian_conserves_condensate(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    q_floor = FT(TDI.TD.Parameters.q_min(tps))
    T_frz = TDI.T_freeze(tps)
    rs_zero = (; ∂ₜq_lcl_frz = zero(FT), ∂ₜq_lcl_shd = zero(FT), ∂ₜq_rai_frz = zero(FT),
        ∂ₜb_rim_lcl = zero(FT), ∂ₜb_rim_rai = zero(FT), f_shd = zero(FT))
    # rows of the 8x8 that carry condensed water; `q_rim` is a fraction of `q_ice` and does
    # not enter the total (see `_condensate_total`)
    cond_rows = (1, 3, 5)
    # strong riming: a healthy supercooled cloud and a healthy rimed ice population, at three
    # supercoolings and with and without rain to exercise both donor columns
    # `q_tot` is set a little above the cell's own liquid saturation so the cloud is
    # maintained rather than evaporating; the supercooling and the rime fraction (0.5) are
    # what make riming the dominant condensate-to-condensate transfer.
    states = (
        (; ρ = FT(0.9), ΔT = FT(8), q_tot = FT(4.5e-3),
            x = FT[1e-3, 1e8, 2e-4, 5e3, 2e-4, 3e5, 1e-4, 2.5e-7]),
        (; ρ = FT(0.78), ΔT = FT(15), q_tot = FT(3.8e-3),
            x = FT[8e-4, 6e7, 5e-4, 2e3, 5e-4, 2e5, 2.5e-4, 5e-7]),
        (; ρ = FT(0.6), ΔT = FT(25), q_tot = FT(1.9e-3),
            x = FT[4e-4, 4e7, 0, 0, 3e-4, 1e5, 1.5e-4, 3e-7]),
    )

    @testset "the riming Jacobian block conserves condensate ($FT)" begin
        for s in states
            x = BMT.MicroState2MP3{FT}(s.x...)
            T = T_frz - s.ΔT
            st = P3.state_from_prognostic(p3, s.ρ * x.q_ice, s.ρ * x.n_ice,
                s.ρ * x.q_rim, s.ρ * x.b_rim)
            logλ = P3.get_distribution_logλ(st)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, s.ρ, T, s.q_tot, logλ)
            pp, rs = _per_process_2mp3_and_riming(mp, tps, s.ρ, T, s.q_tot,
                Tuple(x)..., logλ)

            # the transfer is active, or the test asserts nothing
            @test rs.∂ₜq_lcl_frz < 0
            @test pp.liquid_ice_collision.q_ice > 0

            # the split reproduces the aggregated rates it was cut from, so the Jacobian and
            # the tendency are reading one transfer and not two; each comparison is scaled by
            # its own largest operand, since the three aggregated slots differ by orders
            atol_c(vals...) = 8 * eps(FT) * maximum(abs, vals)
            @test isapprox(rs.∂ₜq_lcl_frz + rs.∂ₜq_lcl_shd,
                pp.liquid_ice_collision.q_lcl;
                atol = atol_c(rs.∂ₜq_lcl_frz, rs.∂ₜq_lcl_shd))
            @test isapprox(rs.∂ₜq_rai_frz - rs.∂ₜq_lcl_shd,
                pp.liquid_ice_collision.q_rai;
                atol = atol_c(rs.∂ₜq_rai_frz, rs.∂ₜq_lcl_shd))
            @test isapprox(-(rs.∂ₜq_lcl_frz + rs.∂ₜq_rai_frz),
                pp.liquid_ice_collision.q_ice;
                atol = atol_c(rs.∂ₜq_lcl_frz, rs.∂ₜq_rai_frz))

            dlcl = 1 / max(q_floor, x.q_lcl)
            drai = 1 / max(q_floor, x.q_rai)
            rb = BMT._riming_jacobian_block(rs, dlcl, drai)

            # EXACT: the cloud column's receivers are the sink negated, so the three condensate
            # contributions cancel bitwise when the receivers are summed first. This is the
            # invariant the block exists for.
            @test (rb.ice_lcl + rb.rai_lcl) + rb.lcl_lcl == 0
            @test rb.ice_rai + rb.rai_rai == 0
            # each donor loses and each receiver gains
            @test rb.lcl_lcl <= 0 && rb.rai_rai <= 0
            @test rb.ice_lcl >= 0 && rb.rai_lcl >= 0 && rb.ice_rai >= 0
            @test rb.rai_rai <= 0

            J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
            J0 = BMT._jacobian_2mp3_manual(g, x, pp, rs_zero)
            @test all(isfinite, J)

            # EXACT: the block writes nine entries via `_riming_jacobian_block`, plus the two
            # wet-growth self-terms the shed-flux driver puts on the rime rows directly, `(8, 8)`
            # for BIWET (card #17 - a relaxation toward a fixed density, not a donor transfer, so
            # it bypasses the block) and `(7, 7)` for QIWET, which is driven by the same `f_shd`
            # and so also differs between the two calls. Eleven entries total, no others.
            #
            # The melting densification's `(8, 7)` is deliberately NOT here: it is driven by
            # `melt_frac`, not by the riming state, so it is the same expression in both calls and
            # listing it would skip a comparison that must be made.
            written = ((1, 1), (5, 1), (3, 1), (3, 3), (5, 3), (7, 1), (7, 3), (8, 1), (8, 3),
                (7, 7), (8, 8))
            for i in 1:8, j in 1:8
                (i, j) in written && continue
                @test J[i, j] == J0[i, j]
            end

            # ROUNDING: the assembled condensate column sums differ between the two matrices
            # by the block's contribution, which is zero up to the accumulation
            for col in (1, 3)
                Δ = sum(J[r, col] - J0[r, col] for r in cond_rows)
                scale = maximum(abs(J[r, col]) for r in cond_rows)
                @test abs(Δ) <= 8 * eps(FT) * scale
            end

            # ...and it is a real cancellation, not an already-zero term: the donor-only form
            # the block replaced left the cloud column short by the whole riming sink
            residue_before = min(pp.liquid_ice_collision.q_lcl, zero(FT)) * dlcl
            @test abs(residue_before) > 0
            @test abs(sum(J[r, 1] - J0[r, 1] for r in cond_rows)) <
                  FT(1e-3) * abs(residue_before)
        end

        # The rain donor's own sink is unconditional, which the aggregated net rate was not:
        # `∂ₜq_r = -QRFRZ + QCSHD` is positive wherever shedding exceeds rain freezing, and
        # `min(., 0)` then dropped the rain diagonal entirely. Asserted on a synthetic split
        # rather than a sampled state, because whether shedding dominates is empirical.
        rs_shed = (; ∂ₜq_lcl_frz = FT(-1e-6), ∂ₜq_lcl_shd = FT(-3e-6),
            ∂ₜq_rai_frz = FT(-1e-6), ∂ₜb_rim_lcl = FT(1e-12), ∂ₜb_rim_rai = FT(1e-12))
        @test rs_shed.∂ₜq_rai_frz - rs_shed.∂ₜq_lcl_shd > 0  # the aggregated rain rate is a source
        rb_s = BMT._riming_jacobian_block(rs_shed, FT(1e3), FT(1e3))
        @test rb_s.rai_rai < 0
        @test (rb_s.ice_lcl + rb_s.rai_lcl) + rb_s.lcl_lcl == 0
        @test rb_s.ice_rai + rb_s.rai_rai == 0
        # and the block is identically zero on a zero split
        @test all(iszero, values(BMT._riming_jacobian_block(rs_zero, FT(1e3), FT(1e3))))
    end

    # Wherever the collision rate is zero the block contributes nothing, so the matrix is the
    # one the donor-only form produced, bit for bit.
    @testset "the riming block is bit-inert where the transfer is off ($FT)" begin
        x = BMT.MicroState2MP3{FT}(1e-3, 1e8, 2e-4, 5e3, 0, 0, 0, 0)  # ice-free: gated off
        T = T_frz - FT(8)
        g = BMT.Instantaneous2MP3Tendency(mp, tps, FT(0.9), T, FT(6e-3), FT(-Inf))
        pp, rs = _per_process_2mp3_and_riming(mp, tps, FT(0.9), T, FT(6e-3),
            Tuple(x)..., FT(-Inf))
        @test all(iszero, values(rs))
        @test all(iszero, Tuple(pp.liquid_ice_collision))
        @test BMT._jacobian_2mp3_manual(g, x, pp, rs) ==
              BMT._jacobian_2mp3_manual(g, x, pp, rs_zero)
    end

    # End to end on the production preset at the box's own step: the substep must leave a
    # physical state, and its condensate must stay inside the cell's total water. That bound
    # is what the mint violated by sixteen orders.
    @testset "the accepted increment stays inside the water budget while riming ($FT)" begin
        manual = BMT.rosenbrock_manual()
        Δt = FT(2)
        for s in states
            T = T_frz - s.ΔT
            st = P3.state_from_prognostic(p3, s.ρ * s.x[5], s.ρ * s.x[6],
                s.ρ * s.x[7], s.ρ * s.x[8])
            logλ = P3.get_distribution_logλ(st)
            t = BMT.bulk_microphysics_tendencies(
                manual, BMT.Microphysics2Moment(), mp, tps,
                s.ρ, T, s.q_tot, s.x..., logλ, Δt, 1,
            )
            x0 = SVector{8, FT}(s.x...)
            x1 = x0 .+ Δt .* SVector{8, FT}(_applied_2mp3(t)...)
            @test all(isfinite, x1)
            @test all(x1 .>= 0)
            # `_water_bounded_increment` bounds a RISING condensate total by `q_tot`, so this
            # is the guarantee as written and not a tighter one
            @test x1[1] + x1[3] + x1[5] <= max(x0[1] + x0[3] + x0[5], s.q_tot)
        end
    end
end

test_riming_jacobian_conserves_condensate(Float64)
test_riming_jacobian_conserves_condensate(Float32)

# Ice aggregation was the only process in the per-process breakdown carrying a non-zero rate
# with no representation of any kind in the manual Jacobian, so a quadratic ice-number sink
# was integrated as a bare `h·f`. The entry is closed form: at frozen `logλ` the rate is
# exactly second order in `n_ice`, so the derivative is `2·rate/n_ice`.
function test_ice_aggregation_jacobian_entry(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    n_floor = FT(TDI.TD.Parameters.q_min(tps))
    T_frz = TDI.T_freeze(tps)
    zero_state = BMT.MicroState2MP3{FT}(0, 0, 0, 0, 0, 0, 0, 0)
    states = (
        (; ρ = FT(0.9), ΔT = FT(10), q_tot = FT(3e-3),
            x = FT[0, 0, 0, 0, 1e-4, 1e5, 0, 0]),
        (; ρ = FT(0.6), ΔT = FT(25), q_tot = FT(1.5e-3),
            x = FT[0, 0, 0, 0, 1e-3, 1e6, 5e-4, 1e-6]),
    )

    @testset "ice aggregation carries its own quadratic diagonal ($FT)" begin
        for s in states
            x = BMT.MicroState2MP3{FT}(s.x...)
            T = T_frz - s.ΔT
            st = P3.state_from_prognostic(p3, s.ρ * x.q_ice, s.ρ * x.n_ice,
                s.ρ * x.q_rim, s.ρ * x.b_rim)
            logλ = P3.get_distribution_logλ(st)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, s.ρ, T, s.q_tot, logλ)
            pp, rs = _per_process_2mp3_and_riming(mp, tps, s.ρ, T, s.q_tot,
                Tuple(x)..., logλ)

            # aggregation is active and is a number sink
            @test pp.ice_aggregation.n_ice < 0

            J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
            # the same matrix with the rate removed from the breakdown: the Jacobian reads
            # `pp` slot by slot, so this isolates the new entry exactly
            pp0 = merge(pp, (; ice_aggregation = zero_state))
            J0 = BMT._jacobian_2mp3_manual(g, x, pp0, rs)

            # EXACT: nothing but the ice-number diagonal moves
            for i in 1:8, j in 1:8
                (i, j) == (6, 6) && continue
                @test J[i, j] == J0[i, j]
            end
            # The entry is the exact second-order donor derivative. The comparison DIFFERENCES two
            # sums that share their sublimation, nucleation and number-adjustment terms, so the
            # precision it can reach is set by the magnitude of that common part and not by the small
            # term being extracted: `J0[6, 6]` runs to ~1e-2 at the rimed state while the entry is
            # ~1e-5, and 8 ulps of the entry asks the subtraction for about 1 part in 1e19 of the
            # numbers it subtracts. Measured (both precisions): the disagreement is 0.16 ulps of
            # `J0[6, 6]` and 159 ulps of `expected`.
            #
            # Where nothing cancels - `J0[6, 6] == 0`, the unrimed state - the two paths agree to
            # within a rounding of the reciprocal, and that is what pins the closed form: the
            # Jacobian multiplies the rate by a reciprocal it shares across the row where this
            # divides, and the two round differently in the last place for some values of the rate
            # and identically for others, so the equality is bit-for-bit or one ulp off depending on
            # the quadrature order the rate was integrated at. Two ulps of the entry admits that and
            # nothing else; a wrong exponent is off by a factor of two. Elsewhere the tolerance is 8
            # ulps of the quantity actually differenced. Do not tighten this back to
            # `abs(expected)`: it held that way only while the aggregation rate ran ~1/E_stick
            # larger than the sticking efficiency now makes it, so the old form was scaled to the
            # wrong quantity throughout and merely had room to spare.
            expected = 2 * pp.ice_aggregation.n_ice / max(n_floor, x.n_ice)
            if J0[6, 6] == 0
                @test isapprox(J[6, 6] - J0[6, 6], expected; atol = 2 * eps(abs(expected)))
            else
                @test isapprox(J[6, 6] - J0[6, 6], expected;
                    atol = 8 * eps(FT) * abs(J0[6, 6]))
            end
            # it damps, and `ExplicitGrowthDiagonal` keeps a damping diagonal
            @test J[6, 6] < J0[6, 6]
            @test BMT._apply_growth(BMT.ExplicitGrowthDiagonal(), J)[6, 6] == J[6, 6]

            # MEASURED, not a theorem: at these ordinary populations the entry is far below
            # the `1/h = 0.5` implicit diagonal of the box's 2 s step, so it cannot be the
            # reason a single step ran away. It bites only where `n_ice` is already large.
            @test abs(expected) * FT(2) < 1
        end
    end

    @testset "the aggregation entry is absent where the rate is ($FT)" begin
        # ice-free: the whole ice block is gated off, so the entry must be exactly zero
        x = BMT.MicroState2MP3{FT}(1e-4, 1e8, 0, 0, 0, 0, 0, 0)
        T = T_frz - FT(10)
        g = BMT.Instantaneous2MP3Tendency(mp, tps, FT(0.9), T, FT(3e-3), FT(-Inf))
        pp, rs = _per_process_2mp3_and_riming(mp, tps, FT(0.9), T, FT(3e-3),
            Tuple(x)..., FT(-Inf))
        @test pp.ice_aggregation.n_ice == 0
        J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
        pp0 = merge(pp, (; ice_aggregation = zero_state))
        @test J == BMT._jacobian_2mp3_manual(g, x, pp0, rs)
    end
end

test_ice_aggregation_jacobian_entry(Float64)
test_ice_aggregation_jacobian_entry(Float32)

# Rain breakup's donor linearization carried `+rate_br/n_rai` unconditionally, but the rate is
# exponential in a mean diameter that falls with its own donor, so the true derivative changes
# sign. Summed with self-collection the pair went positive for every mean volume drop diameter
# above `Deq = 0.9` mm, and `ExplicitGrowthDiagonal` then deleted the rain-number diagonal -
# self-collection's genuine damping with it.
function test_rain_number_diagonal_sign(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    brek = mp.warm_rain.seifert_beheng.brek
    ρw = mp.warm_rain.seifert_beheng.pdf_r.ρw
    zero_state = BMT.MicroState2MP3{FT}(0, 0, 0, 0, 0, 0, 0, 0)
    ρ = FT(1)
    T = FT(290)          # warm and ice-free: rain evaporation, self-collection and breakup only
    q_rai = FT(1e-3)
    # pick `n_rai` from the target mean volume diameter, `x̄ = π ρw D³/6`, spanning the
    # equilibrium diameter `Deq = 0.9` mm and the sign change of the breakup derivative at
    # `3/κbr = 1.30` mm, up to the limited PDF's largest admissible drop
    Drs = FT[3e-4, 6e-4, 8e-4, 1.0e-3, 1.4e-3, 2.0e-3]

    @testset "the rain-number diagonal is a damping at every drop size ($FT)" begin
        @test brek.Deq ≈ FT(9e-4)
        for Dr in Drs
            xr = FT(π) * ρw * Dr^3 / 6
            n_rai = ρ * q_rai / xr
            q_tot = FT(0.014)
            x = BMT.MicroState2MP3{FT}(2e-4, 5e7, q_rai, n_rai, 0, 0, 0, 0)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, FT(-Inf))
            pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
                Tuple(x)..., FT(-Inf))
            J = BMT._jacobian_2mp3_manual(g, x, pp, rs)

            # self-collection is a sink and breakup a source, at every size in the band
            @test pp.rain_selfcol.n_rai < 0
            @test pp.rain_breakup.n_rai >= 0

            # the diagonal damps, so the growth treatment leaves it alone. Before the gate the
            # pair went positive above `Deq` and the whole entry was deleted.
            @test J[4, 4] < 0
            @test BMT._apply_growth(BMT.ExplicitGrowthDiagonal(), J)[4, 4] == J[4, 4]

            # EXACT: above `Deq` the breakup rate no longer enters the diagonal at all, so
            # removing it from the breakdown changes nothing
            pp0 = merge(pp, (; rain_breakup = zero_state))
            J0 = BMT._jacobian_2mp3_manual(g, x, pp0, rs)
            above = pp.rain_selfcol.n_rai + pp.rain_breakup.n_rai > 0
            if above
                @test Dr > brek.Deq
                @test J[4, 4] == J0[4, 4]
            else
                # below `Deq` the pair is still a damping and is carried unchanged
                @test J[4, 4] >= J0[4, 4]
            end
        end
    end

    @testset "the growth treatment no longer acts on the manual matrix ($FT)" begin
        # MEASURED at these states, not proved: with the rain-number diagonal fixed, no
        # diagonal of the manual 2M+P3 Jacobian is positive, so `ExplicitGrowthDiagonal` is
        # inert on the production path rather than silently deleting one entry.
        T_frz = TDI.T_freeze(tps)
        cases = (
            (; ρ = FT(1), T = FT(290), q_tot = FT(0.014), logλ = FT(-Inf),
                x = FT[2e-4, 5e7, 1e-3, 2.4e2, 0, 0, 0, 0]),
            (; ρ = FT(0.9), T = T_frz - FT(8), q_tot = FT(4.5e-3), logλ = nothing,
                x = FT[1e-3, 1e8, 2e-4, 5e3, 2e-4, 3e5, 1e-4, 2.5e-7]),
            (; ρ = FT(0.6), T = T_frz - FT(25), q_tot = FT(1.9e-3), logλ = nothing,
                x = FT[4e-4, 4e7, 0, 0, 3e-4, 1e5, 1.5e-4, 3e-7]),
        )
        for c in cases
            x = BMT.MicroState2MP3{FT}(c.x...)
            logλ = if isnothing(c.logλ)
                st = P3.state_from_prognostic(mp.ice.scheme, c.ρ * x.q_ice, c.ρ * x.n_ice,
                    c.ρ * x.q_rim, c.ρ * x.b_rim)
                P3.get_distribution_logλ(st)
            else
                c.logλ
            end
            g = BMT.Instantaneous2MP3Tendency(mp, tps, c.ρ, c.T, c.q_tot, logλ)
            pp, rs = _per_process_2mp3_and_riming(mp, tps, c.ρ, c.T, c.q_tot,
                Tuple(x)..., logλ)
            J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
            @test BMT._apply_growth(BMT.ExplicitGrowthDiagonal(), J) == J
        end
    end
end

test_rain_number_diagonal_sign(Float64)
test_rain_number_diagonal_sign(Float32)

# The rime-density ratchet: melting must remove rime along the ray through the origin.
#
# `ρ_rim = q_rim/b_rim` is a derived quotient of two prognostic moments, so its physicality is a
# property of every `(dq_rim, db_rim)` pair. A removal pair preserves the bulk quotient exactly if
# and only if it exits the `(q_rim, b_rim)` plane along a ray through the origin; a removal at any
# other implied density moves the quotient AWAY from that density, without bound.
#
# The replaced melt drain removed volume at `state.ρ_rim`, the derived quotient CLAMPED to `ρ_i` at
# state construction. On a bulk sitting at `ρ_i(1 + δ)` the clamp makes the drain's implied density
# `ρ_i ≠ ρ_bulk`, so removing all but a fraction `f` of the rime leaves
#
#     ρ(f)/ρ_i = f(1 + δ) / (f − δ(1 − f)),        i.e.  δ' ≈ δ/f  for small δ,
#
# and `b_rim` reaches zero FIRST at `f* = δ/(1 + δ)`, where the positivity clamp strands an orphan
# `q_rim > 0` with `b_rim = 0` and a formally infinite quotient. Near-complete melt-out is the
# normal fate of a sedimenting rimed particle, so this is a one-signed amplifier with gain `1/f` per
# melting episode - the consumer clamp is not merely failing to fix the prognostics, it IS the
# mechanism that manufactures them.
#
# Control here is the replaced arithmetic written out; treatment is the ray form, both in isolation
# and on the production substep path.
function test_rime_melt_drain_is_ray_form(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    ρ_i = p3.ρ_i
    ρ = FT(0.8)
    n_ice = FT(1e4)

    # The replaced drain: rime mass leaves at `F_rim`, rime VOLUME at the clamped derived quotient.
    # Both ratios are taken from `state_from_prognostic`, so this is the regularisation the code ran.
    function control_drain(q_ice, q_rim, b_rim, Δq_melt)
        st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, ρ * q_rim, ρ * b_rim)
        dq = Δq_melt * st.F_rim
        db = st.ρ_rim > 0 ? Δq_melt * st.F_rim / st.ρ_rim : zero(FT)
        return (max(q_rim - dq, zero(FT)), max(b_rim - db, zero(FT)))
    end
    # The ray form: one fractional loss, both moments.
    function ray_drain(q_ice, q_rim, b_rim, Δq_melt)
        frac = Δq_melt / max(q_ice, floatmin(FT))
        return (max(q_rim - q_rim * frac, zero(FT)), max(b_rim - b_rim * frac, zero(FT)))
    end

    @testset "the melt drain is ray-form [FT=$FT]" begin
        # --- closed form: the control reproduces the ratchet, the treatment does not have one
        @testset "closed form" begin
            q_ice0, q_rim0 = FT(1e-4), FT(9e-5)   # F_rim = 0.9, away from the F_rim < 1 clamp
            # δ resolved by the precision: at Float32, eps = 1.2e-7 and ρ_i(1 + δ) has to hold δ
            δs = FT == Float64 ? FT[1e-6, 1e-4, 1e-3] : FT[1e-4, 1e-3]
            for δ in δs, f in FT[1e-1, 1e-2, 1e-3]
                b_rim0 = q_rim0 / (ρ_i * (1 + δ))
                Δq = (1 - f) * q_ice0            # melt all but the fraction f of the ice mass
                @test Δq > 0
                # Removing all but a fraction `f` in ONE bite is a cancellation of that severity,
                # so a quotient that survives it exactly still carries `O(eps/f)` of rounding.
                # That is the floor these assertions measure against, and it is not the ratchet:
                # it is unbiased, it is not amplified (the next drain re-anchors on the current
                # pair), and at `f = 1e-3` it is ten orders below the excursion the control makes.
                tol = 32 * eps(FT) / f

                qc, bc = control_drain(q_ice0, q_rim0, b_rim0, Δq)
                f★ = δ / (1 + δ)                 # where the control's b_rim reaches zero
                if f > f★ * (1 + 64 * eps(FT))
                    # the ratchet, in closed form. The right-hand side cancels a second time as
                    # `f` approaches `f★`, hence the `(1 + δ/f)` factor on the tolerance.
                    @test bc > 0
                    # the denominator `f − δ(1−f)` cancels a second time as `f` approaches `f★`,
                    # by the factor `f/(f − δ(1−f))` - a thousandfold at δ = f = 1e-3, which is
                    # the state the design's own 1001·ρ_i number comes from
                    cancel = f / abs(f - δ * (1 - f))
                    @test qc / bc ≈ ρ_i * f * (1 + δ) / (f - δ * (1 - f)) rtol =
                        2 * tol * max(one(FT), cancel)
                    δ′ = qc / bc / ρ_i - 1
                    @test δ′ > δ                 # one-signed: the excursion only ever grows
                    # the small-δ gain, asserted only where the approximation it is is valid
                    f > 4 * f★ && @test δ′ ≈ δ / f rtol = FT(0.05) + 4 * δ / f
                else
                    # past f★ the volume moment is truncated first and the pair is orphaned
                    @test bc == 0 && qc > 0
                end

                qr, br = ray_drain(q_ice0, q_rim0, b_rim0, Δq)
                @test br > 0                     # the pair never orphans: both scale by `f`
                @test qr / br ≈ q_rim0 / b_rim0 rtol = tol
                @test qr ≈ f * q_rim0 rtol = tol
                @test br ≈ f * b_rim0 rtol = tol
            end
        end

        # --- the three numbers the design was written against, at the precision they were computed
        if FT == Float64
            @testset "the design's arithmetic" begin
                q_ice0, q_rim0 = FT(1e-4), FT(9e-5)
                ratchet(δ, f) = begin
                    b0 = q_rim0 / (ρ_i * (1 + δ))
                    q, b = control_drain(q_ice0, q_rim0, b0, (1 - f) * q_ice0)
                    (q, b)
                end
                q, b = ratchet(1e-6, 1e-3)       # δ = 1e-6 melted to f = 1e-3 gives δ′ = 1e-3
                @test q / b / ρ_i - 1 ≈ 1e-3 rtol = 2e-3
                q, b = ratchet(1e-3, 1e-3)       # δ = 1e-3 melted to f = 1e-3 gives ρ = 1001 ρ_i
                @test q / b / ρ_i ≈ 1001 rtol = 1e-2
                q, b = ratchet(0.0185, 0.018)    # δ = 0.0185 cannot survive a melt past f = 0.018
                @test b == 0 && q > 0
            end
        end

        # --- a melting episode integrated in many small steps
        @testset "a melting episode" begin
            δ = FT == Float64 ? FT(1e-5) : FT(1e-4)
            nstep = 20
            f_end = FT(1e-3)
            step_fac = f_end^(1 / FT(nstep))     # equal fractional bites of the ice mass
            q_ice0, q_rim0 = FT(1e-4), FT(9e-5)
            b_rim0 = q_rim0 / (ρ_i * (1 + δ))

            qc, bc, qr, br, q_ice = q_rim0, b_rim0, q_rim0, b_rim0, q_ice0
            for _ in 1:nstep
                Δq = q_ice * (1 - step_fac)
                qc, bc = control_drain(q_ice, qc, bc, Δq)
                qr, br = ray_drain(q_ice, qr, br, Δq)
                q_ice -= Δq
            end
            @test q_ice ≈ f_end * q_ice0 rtol = 64 * eps(FT)

            # control: the excursion has been amplified by ~1/f, or the pair orphaned outright
            if bc == 0
                @test qc > 0                     # orphaned: mass left with no volume
            else
                δ_ctrl = qc / bc / ρ_i - 1
                @test δ_ctrl > 100 * δ           # the gain is 1/f_end = 1000; 100 is the margin
                @test δ_ctrl ≈ δ / f_end rtol = FT(0.3)
            end
            # treatment: the quotient is the one it started with, to ulps, and the pair survives
            @test br > 0 && qr > 0
            @test qr / br ≈ ρ_i * (1 + δ) rtol = 4 * nstep * eps(FT)
            # ... and it is still inside the physical interval, which is the point
            @test qr / br ≤ ρ_i * (1 + δ) * (1 + 4 * nstep * eps(FT))
        end

        # --- the production path: `_per_process_2mp3`, the entry, and the Jacobian
        @testset "the production substep" begin
            T_freeze = TDI.TD.Parameters.T_freeze(tps)
            T = T_freeze + FT(5)                  # melting live
            q_tot = FT(5e-3)
            q_ice, q_rim = FT(1e-4), FT(9e-5)
            for δ in FT[0, 1e-3, 2e-2]            # at, just above, and well above solid ice
                b_rim = q_rim / (ρ_i * (1 + δ))
                st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, ρ * q_rim, ρ * b_rim)
                logλ = P3.get_distribution_logλ(st)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
                    FT(0), FT(0), FT(0), FT(0),   # droplet-free, rain-free: melt is the only
                    q_ice, n_ice, q_rim, b_rim,   # process writing the rime pair here
                    logλ)
                # the slot must be live, or this measures nothing
                @test pp.ice_melting.q_ice < 0
                @test pp.ice_melting.q_rim < 0 && pp.ice_melting.b_rim < 0

                # The drain is NOT along the ray: melting densifies the rime toward solid
                # ice, so the volume drains faster than the mass by exactly ρ_i / ρ_rim and the
                # drain's implied density is ρ_rim² / ρ_i rather than the bulk quotient itself.
                # At δ = 0 the two coincide, which is the degenerate case the sweep includes.
                ρ_rim_state = q_rim / b_rim
                @test pp.ice_melting.q_rim / pp.ice_melting.b_rim ≈
                      ρ_rim_state^2 / ρ_i rtol = 8 * eps(FT)
                # equivalently: one fractional mass loss, and a volume loss ρ_i / ρ_rim times it
                @test pp.ice_melting.b_rim / b_rim ≈
                      (ρ_i / ρ_rim_state) * (pp.ice_melting.q_rim / q_rim) rtol = 8 * eps(FT)
                # and it is the ice mass's own fractional loss
                @test pp.ice_melting.q_rim / q_rim ≈ pp.ice_melting.q_ice / q_ice rtol =
                    8 * eps(FT)

                # the entry forms the same rates a second time, and it is what
                # `Instantaneous2MP3Tendency` evaluates
                x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, q_ice, n_ice, q_rim, b_rim)
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
                psum = sum(values(pp))
                @test g(x).q_rim ≈ psum.q_rim rtol = 8 * eps(FT)
                @test g(x).b_rim ≈ psum.b_rim rtol = 8 * eps(FT)

                # The two rime diagonals are no longer the same fractional loss, because the
                # melting densification reshapes the pair rather than rescaling it: the rime MASS
                # drains along the ray so its self-derivative is -melt_frac, while the rime VOLUME
                # carries -2 melt_frac ρ_i / ρ_rim. Sublimation and the orphan drain are ray-form
                # and add the SAME term to both rows, so their whole difference is the melting one
                # and asserting the difference measures the densification alone. Above freezing
                # the entry suppresses ice DEPOSITION but not sublimation, which is why J[7, 7]
                # itself is the sum of two fractional losses below.
                J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
                melt_frac = -pp.ice_melting.q_ice / q_ice
                @test J[8, 8] - J[7, 7] ≈
                      melt_frac * (1 - 2 * ρ_i / ρ_rim_state) rtol = 64 * eps(FT)
                @test J[8, 8] < J[7, 7] < 0 && isfinite(J[8, 8])
                drains = pp.ice_melting.q_ice + min(pp.ice_depsub.q_ice, zero(FT))
                @test J[7, 7] ≈ drains / q_ice rtol = 8 * eps(FT)
                @test all(isfinite, J)
            end
        end

        # --- the ICE-FREE family: the drains' fractional loss must vanish in the partials too
        #
        # The fraction is `∂ₜq_ice_sink / q_ice`, and forming it as `rate / max(q_ice, floatmin)`
        # gives the right VALUE at an ice-free state and NON-FINITE DERIVATIVES: the quotient rule
        # weights `∂q_ice` by `−rate/q_ice²`, and `floatmin(FT)^2` underflows to exactly zero. The
        # value lane looks healthy throughout, which is why a band sweep over rimed states cannot
        # see it - the ice-free family is outside the band by construction. Caught in the wild by
        # `ad_compat_tests`' warm-rain regime, where both ice gates are shut and every ice moment
        # is exactly zero.
        #
        # This is the falsifier for the presence gate: at states with no ice the rime rows must be
        # exactly zero and EVERY partial of the entry must be finite. It is the same defect family
        # as a bracketing solver returning a finite root with NaN derivatives - a guard sized for
        # the value is blind to the partials.
        @testset "ice-free states leave finite partials" begin
            # The predicate itself, at the operation level. A `Dual` whose value is zero and whose
            # partials are not must read as ABSENT: `q > zero(q)` reads it as present and then runs
            # the division the guard exists to avoid, which is how this defect survived a first fix.
            # `UT.guarded_quotient` is the one helper both lanes now share; these assertions are
            # what keep the two old idioms - floored denominator, naked select - red.
            let rate = FD.Dual{Nothing}(FT(0), FT(1), FT(0)), q = FD.Dual{Nothing}(FT(0), FT(0), FT(1))

                # THE DISCRIMINATOR, asserted rather than inferred. The two predicates disagree
                # on exactly this Dual, and that disagreement is the whole defect: the value lane
                # reads ABSENT, the Dual comparison reads PRESENT. Any guard built on the second
                # runs the division it exists to avoid, however the quotient is then written - so
                # the two-sided select (`b_safe = ifelse(active, b, one(b))`) is NOT sufficient on
                # its own either, because it tests `active` the same way. Reading `FD.value` is
                # what makes it work; the safe denominator is belt to that braces.
                @test !(FD.value(q) > zero(FD.value(q)))   # value lane: absent
                @test q > zero(q)                          # Dual comparison: present - the trap

                f = UT.guarded_quotient(rate, q)
                @test iszero(FD.value(f))
                @test all(iszero, FD.partials(f))
                # the naked select, written out here as the CONTROL: with the same predicate the
                # guard leaks, which is what a bare `ifelse(q > 0, ...)` reduces to
                @test !all(isfinite, FD.partials(ifelse(q > zero(q), rate / q, zero(rate / q))))
                # and it is still the exact quotient wherever the population is there
                qp = FD.Dual{Nothing}(FT(2), FT(0), FT(1))
                @test FD.value(UT.guarded_quotient(rate, qp)) == 0
                rp = FD.Dual{Nothing}(FT(3), FT(1), FT(0))
                @test FD.value(UT.guarded_quotient(rp, qp)) == FT(1.5)
                @test all(isfinite, FD.partials(UT.guarded_quotient(rp, qp)))
                # the explicit absent value, which is what makes one helper serve both lanes: a
                # budget with nothing to apportion does not bind, so its cap is Inf, not zero
                @test UT.guarded_quotient(rate, q, oftype(FD.value(q), Inf)) == FT(Inf)
            end
            function rhs(x, ρ, T, q_tot, logλ)
                t = BMT.bulk_microphysics_tendencies(
                    BMT.Microphysics2Moment(), mp, tps, ρ, T, q_tot,
                    x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], logλ)
                return [t.dq_lcl_dt, t.dn_lcl_dt, t.dq_rai_dt, t.dn_rai_dt,
                    t.dq_ice_dt, t.dn_ice_dt, t.dq_rim_dt, t.db_rim_dt]
            end
            # `ad_compat_tests`' own warm-rain regime, plus two colder ice-free states so the
            # sublimation branch is live at one of them and the deposition branch at the other -
            # a single warm state would leave the gate untested on the branch that actually runs
            regimes = (
                (; name = "warm rain", ρ = FT(1.05), T = FT(288), q_tot = FT(0.015)),
                (; name = "cold subsaturated", ρ = FT(0.6), T = FT(245), q_tot = FT(1e-4)),
                (; name = "cold supersaturated", ρ = FT(0.6), T = FT(245), q_tot = FT(4e-3)),
            )
            for r in regimes
                x0 = FT[4e-4, 8e7, 2.1e-3, 5e4, 0, 0, 0, 0]   # every ice moment exactly zero
                @test all(iszero, x0[5:8])
                f = x -> rhs(x, r.ρ, r.T, r.q_tot, FT(-Inf))
                v = f(x0)
                J = FD.jacobian(f, x0)
                @test all(isfinite, v)
                @test all(isfinite, J)                        # all 64 partials
                # Where no rime SOURCE is active - above freezing, so every freezing pathway is
                # off - the rime rows are exactly zero in value and in every derivative. Below
                # freezing they are NOT, and that is correct physics rather than a leak: Bigg
                # immersion freezes cloud droplets into fully rimed embryo graupel, so an
                # ice-free supercooled state with liquid present legitimately creates rime.
                if r.T > TDI.TD.Parameters.T_freeze(tps)
                    @test v[7] == 0 && v[8] == 0
                    @test all(iszero, J[7, :]) && all(iszero, J[8, :])
                else
                    @test isfinite(v[7]) && isfinite(v[8])
                end
            end
        end

        # --- inert where the replaced form was already correct: a consistent interior state
        @testset "inert on a consistent state" begin
            T = TDI.TD.Parameters.T_freeze(tps) + FT(5)
            q_tot = FT(5e-3)
            q_ice, n_ice_c = FT(1e-4), FT(1e4)
            # a physically ordinary rime density, well inside [ρ′(1), ρ_i] and off both the
            # `ρ_rim ≤ ρ_i` clamp and the rime-volume presence cutoff
            for (F_rim, ρ_rim) in ((FT(0.3), FT(400)), (FT(0.9), FT(700)))
                q_rim = F_rim * q_ice
                b_rim = q_rim / ρ_rim
                st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice_c, ρ * q_rim, ρ * b_rim)
                @test st.ρ_rim ≈ ρ_rim rtol = 8 * eps(FT)   # neither clamp nor taper binds
                logλ = P3.get_distribution_logλ(st)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
                    FT(0), FT(0), FT(0), FT(0), q_ice, n_ice_c, q_rim, b_rim, logλ)
                # the rime MASS still drains along the ray, so the ray form predicts it exactly
                ray_q = -(-pp.ice_melting.q_ice) * st.F_rim
                ray_b = -(-pp.ice_melting.q_ice) * st.F_rim / st.ρ_rim
                @test pp.ice_melting.q_rim ≈ ray_q rtol = 32 * eps(FT)
                # the rime VOLUME does not: the densification drains it ρ_i / ρ_rim times faster,
                # which is where this state's rime density sits below solid ice
                @test pp.ice_melting.b_rim ≈
                      ray_b * p3.ρ_i / st.ρ_rim rtol = 32 * eps(FT)
                @test pp.ice_melting.b_rim < ray_b < 0
            end
        end
    end
end

# the per-process decomposition sums to the entry in every slot, at an F23-active state
function test_f23_shift_has_one_definition(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    ρ, T = FT(0.8), FT(245)
    q_tot = FT(4e-3)              # supersaturated over ice at 245 K, so F23 is live
    q_ice, n_ice = FT(1e-6), FT(1e3)
    x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, q_ice, n_ice, 0, 0)
    st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, FT(0), FT(0))
    logλ = P3.get_distribution_logλ(st)

    @testset "the F23 INP shift has one definition [FT=$FT]" begin
        pp = _per_process_2mp3(mp, tps, ρ, T, q_tot,
            FT(0), FT(0), FT(0), FT(0), q_ice, n_ice, FT(0), FT(0), logλ)
        @test pp.ice_deposition.n_ice > 0
        g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
        psum = sum(values(pp))
        for f in (:q_lcl, :n_lcl, :q_rai, :n_rai, :q_ice, :n_ice, :q_rim, :b_rim)
            @test getfield(g(x), f) ≈ getfield(psum, f) rtol = 16 * eps(FT)
        end
    end
end

test_f23_shift_has_one_definition(Float64)
test_f23_shift_has_one_definition(Float32)

test_rime_melt_drain_is_ray_form(Float64)
test_rime_melt_drain_is_ray_form(Float32)

# Sublimation removes rime along the same ray melting does.
#
# The sublimation rim drain had the melt drain's structure and so its defect: mass left at
# `F_rim`, volume at the regularised and CLAMPED `state.ρ_rim`, and the `ifelse(ρ_rim > 0)` guard
# drained mass while leaving volume behind whenever the taper zeroed the quotient. Reusing the
# `sub_frac_ice` the number pathway already forms makes the drain's implied density the bulk
# quotient at every state, so sublimation - like the mean particle mass under the (q, n) pairing -
# leaves the rime density alone, and the last of the rime mass takes the last of the rime volume.
function test_rime_sublimation_drain_is_ray_form(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    ρ_i = p3.ρ_i
    ρ = FT(0.6)
    q_ice, n_ice = FT(1e-5), FT(1e4)
    q_tot = FT(1e-4)              # subsaturated over ice at every T below
    Ts = FT[233, 245, 253]

    @testset "the sublimation drain is ray-form [FT=$FT]" begin
        for δ in FT[0, 1e-3, 2e-2], T in Ts
            q_rim = FT(0.9) * q_ice
            b_rim = q_rim / (ρ_i * (1 + δ))
            st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, ρ * q_rim, ρ * b_rim)
            logλ = P3.get_distribution_logλ(st)
            qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
            @test q_tot < qᵥ_sat_ice     # the sublimation branch is the live one
            τ = P3.ice_deposition_timescale(
                mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps,
                T, ρ, st, logλ; quad = mp.ice.quad)
            @test !P3.ice_deposition_is_degenerate(τ)

            pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
                FT(0), FT(0), FT(0), FT(0), q_ice, n_ice, q_rim, b_rim, logλ)
            @test pp.ice_depsub.q_ice < 0
            @test pp.ice_depsub.q_rim < 0 && pp.ice_depsub.b_rim < 0

            # one fractional loss, all three of ice mass, rime mass and rime volume
            @test pp.ice_depsub.q_rim / q_rim ≈ pp.ice_depsub.q_ice / q_ice rtol = 8 * eps(FT)
            @test pp.ice_depsub.b_rim / b_rim ≈ pp.ice_depsub.q_ice / q_ice rtol = 8 * eps(FT)
            # so the drain's implied density is the bulk quotient, at and above the ρ_i clamp
            @test pp.ice_depsub.q_rim / pp.ice_depsub.b_rim ≈ q_rim / b_rim rtol = 8 * eps(FT)

            # the pair reaches zero TOGETHER: the step that empties the rime mass empties the
            # rime volume, so the positivity clamp truncates the pair jointly instead of
            # stranding one moment. (The two crossings coincide because the factor is shared.)
            h_q = -q_rim / pp.ice_depsub.q_rim
            h_b = -b_rim / pp.ice_depsub.b_rim
            @test h_q ≈ h_b rtol = 8 * eps(FT)

            # the entry forms the same rates, and it is what the substep evaluates
            x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, q_ice, n_ice, q_rim, b_rim)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
            psum = sum(values(pp))
            @test g(x).q_rim ≈ psum.q_rim rtol = 8 * eps(FT)
            @test g(x).b_rim ≈ psum.b_rim rtol = 8 * eps(FT)

            # J carries the drain f runs, exactly and equally on both rim rows. Melting is off
            # below freezing, so the sublimation drain owns these two diagonals here.
            J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
            @test J[7, 7] == J[8, 8]
            @test J[7, 7] < 0 && isfinite(J[7, 7])
            @test J[7, 7] ≈ pp.ice_depsub.q_ice / q_ice rtol = 8 * eps(FT)
            @test all(isfinite, J)
        end

        # depositional growth is pristine: it carries no rime at all, in either moment
        @testset "the deposition branch carries no rime" begin
            q_tot_super = FT(4e-3)
            q_rim = FT(0.9) * q_ice
            b_rim = q_rim / ρ_i
            st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, ρ * q_rim, ρ * b_rim)
            logλ = P3.get_distribution_logλ(st)
            for T in Ts
                @test q_tot_super > TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot_super,
                    FT(0), FT(0), FT(0), FT(0), q_ice, n_ice, q_rim, b_rim, logλ)
                @test pp.ice_depsub.q_ice > 0
                @test pp.ice_depsub.q_rim == 0 && pp.ice_depsub.b_rim == 0
            end
        end

        # inert where the replaced form was already correct
        @testset "inert on a consistent state" begin
            T = FT(245)
            for (F_rim, ρ_rim) in ((FT(0.3), FT(400)), (FT(0.9), FT(700)))
                q_rim = F_rim * q_ice
                b_rim = q_rim / ρ_rim
                st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, ρ * q_rim, ρ * b_rim)
                @test st.ρ_rim ≈ ρ_rim rtol = 8 * eps(FT)   # neither clamp nor taper binds
                logλ = P3.get_distribution_logλ(st)
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
                    FT(0), FT(0), FT(0), FT(0), q_ice, n_ice, q_rim, b_rim, logλ)
                old_q = pp.ice_depsub.q_ice * st.F_rim
                old_b = pp.ice_depsub.q_ice * st.F_rim / st.ρ_rim
                @test pp.ice_depsub.q_rim ≈ old_q rtol = 32 * eps(FT)
                @test pp.ice_depsub.b_rim ≈ old_b rtol = 32 * eps(FT)
            end
        end
    end
end

test_rime_sublimation_drain_is_ray_form(Float64)
test_rime_sublimation_drain_is_ray_form(Float32)

"""
    test_biwet_diagonal_restored(FT)

Card #17: the manual Jacobian's `brim_brim` entry (`J[8,8]`) previously carried no
contribution from wet growth (`BIWET`, [`CMP3.bulk_liquid_ice_collision_sources`](@ref)) -
only the ray-form melt/sublimation/orphan drains, each contributing identically to both
`rim_rim` (`J[7,7]`) and `brim_brim`. Every OTHER riming-Jacobian test in this file
deliberately zeros cloud/rain input to isolate those ray-form drains alone, which also
zeros `f_shd` (the collision integral has nothing to collect) - so none of them exercise
this diagonal. Here cloud AND rain collection are both active, `f_shd > 0`, and `BIWET`'s
own self-term `-f_shd/τ_wet` is exact and closed-form with `state`/`logλ`/the quadrature
rates held fixed (only `B_rim`'s explicit factor varies - see
`notes/card17-biwet-decomposition.md`), and is now restored on the diagonal.

Two things matter beyond the closed-form match: `J[7,7] != J[8,8]` here, CORRECTLY, not
accidentally - wet growth relaxes the `(L_rim, B_rim)` pair TOWARD a fixed density
endpoint rather than draining it along a fixed ray, so mass and volume respond
differently and the ray-form symmetry every sibling test exercises does not apply. And
the fix is not asserted to close the AD gap exactly - the geometry-mediated remainder
(BCCOL/BRCOL's own dependence on `F_rim`/`ρ_rim` through `state`) is a measured, documented
Tier-3 limitation, so only that the fix moves the manual entry strictly closer to the AD
one is checked here.
"""
function test_biwet_diagonal_restored(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    τ_wet = p3.τ_wet
    ρ, T, q_tot = FT(0.8), FT(268.0), FT(6e-3)
    q_lcl, n_lcl = FT(1e-3), FT(1e8)
    q_rai, n_rai = FT(5e-4), FT(2e5)
    q_ice, n_ice = FT(2e-4), FT(5e5)
    F_rim, ρ_rim = FT(0.4), FT(500.0)
    q_rim = F_rim / (1 - F_rim) * q_ice
    b_rim = q_rim / ρ_rim

    @testset "BIWET's self-term is on the diagonal [FT=$FT]" begin
        st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, ρ * q_rim, ρ * b_rim)
        logλ = P3.get_distribution_logλ(st)
        pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot,
            q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ)

        # the state must actually be wet-growing, or this measures nothing
        @test rs.f_shd > 0 && rs.f_shd <= 1

        x = BMT.MicroState2MP3{FT}(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
        g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
        J_m = BMT._jacobian_2mp3_manual(g, x, pp, rs)
        _, J_e = BMT._tendency_and_jacobian(BMT.ExactJacobian(), g, x)
        J_e88 = Matrix(J_e)[8, 8]

        # (1) the restored entry is exactly the closed form
        @test J_m[8, 8] ≈ -rs.f_shd / τ_wet rtol = 8 * eps(FT)

        # (2) both rime rows now take the SAME wet-growth self-term, because the shed-flux
        # driver replaced the wet-growth indicator and drives QIWET and BIWET alike. At this
        # subfreezing state melting is inactive, so nothing else separates the two rows and the
        # diagonals coincide; the asymmetry this clause used to assert belonged to the retired
        # indicator. What still distinguishes BIWET is clause (1), which pins its closed form.
        @test J_m[7, 7] ≈ J_m[8, 8] rtol = 8 * eps(FT)

        # (3) the fix moves the manual entry TOWARD the AD entry rather than overshooting
        # or moving the wrong way; the residual is the documented geometry remainder, not
        # asserted to vanish here
        old_gap = abs(J_e88 - zero(FT))
        new_gap = abs(J_e88 - J_m[8, 8])
        @test new_gap < old_gap
        @test all(isfinite, Matrix(J_m))
    end
end

test_biwet_diagonal_restored(Float64)
test_biwet_diagonal_restored(Float32)

# Composition falsifier: the rime-density interval must be forward-invariant under the PRODUCTION
# substep, not merely under each process in isolation.
#
# Each source pair deposits at an implied density inside `[ρ′_rim(1), ρ_i]` and each sink now
# removes along the ray through the origin, so process by process the mediant property confines
# the bulk quotient to that interval. What that argument does NOT cover is the composition: one
# Rosenbrock increment mixes all of them through a linear solve whose row weights are positive
# only in the diagonally dominant regime, and the positivity clamp, the water bound and the
# saturation adjustment all act afterwards. This drives random reachable states through
# `rosenbrock_manual()` - the box's own configuration - and asserts what the design claims:
#
#   1. the quotient stays inside the interval, up to an accumulated per-step ulp tolerance;
#   2. no orphan pair is ever produced, in either direction.
#
# The tolerance is derived rather than tuned: two moments each carrying `O(1)` ulps of rounding per
# step give `O(2K)` ulps of relative drift in their quotient over `K` steps, and states sitting AT
# `ρ_i` - which freezing and wet growth both drive them to - can only leave by that much. A failure
# outside it is a real composition excursion and the design's own pre-registered branch: it would
# mean the drains are pair-consistent and the SOLVE is not.
function test_rime_density_interval_is_forward_invariant(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    ρ_min, ρ_max = P3.rime_density_bounds(p3)
    mode = BMT.rosenbrock_manual()
    Δt = FT(2)          # the production box step
    nsub = 1            # the production substep count
    K = 20              # steps per trajectory
    tol = 16 * K * eps(FT)

    # A deterministic LCG rather than `Random`, so the sweep is reproducible without depending on
    # the RNG stream of a particular Julia version.
    seed = UInt64(20260727)
    nextrand() = (seed = (0x5851f42d4c957f2d * seed + 0x14057b7ef767814f); FT(seed >> 11) / FT(2^53))
    logu(lo, hi) = exp(log(lo) + nextrand() * (log(hi) - log(lo)))
    u(lo, hi) = lo + nextrand() * (hi - lo)

    n_traj = 48
    worst_hi, worst_lo, n_orphan, n_live = zero(FT), FT(Inf), 0, 0
    # the largest rime mass ever stranded on an ice-free cell, for the residue bound below
    resid = zero(FT)
    ϵₘ = eps(FT)
    # the rimed-cell threshold of the event census: below this a cell carries no rime worth naming
    Q_RIM_RIMED = FT(1e-9)
    resid_by_number = false

    @testset "the rime-density interval is forward-invariant [FT=$FT]" begin
        @test ρ_min < ρ_max
        for traj in 1:n_traj
            ρ = FT(u(0.4, 1.2))
            T = FT(u(220, 300))
            q_lcl = FT(logu(1e-8, 2e-3))
            n_lcl = FT(logu(1e6, 5e8))
            q_rai = FT(logu(1e-8, 2e-3))
            n_rai = FT(logu(1e2, 1e6))
            q_ice = FT(logu(1e-8, 2e-3))
            n_ice = FT(logu(1e2, 1e6))
            # total water carries the condensate plus positive vapour, so the states are
            # reachable rather than already violating the water budget
            q_tot = q_lcl + q_rai + q_ice + FT(logu(1e-5, 2e-2))
            # a REACHABLE rime state: inside the interval by construction, which is the premise
            # forward invariance is asserted from
            F_rim = FT(u(0, 0.98))
            ρ_rim = FT(u(ρ_min, ρ_max))
            q_rim = F_rim * q_ice
            b_rim = q_rim / ρ_rim

            x = FT[q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim]
            for step in 1:K
                st = P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8])
                logλ = P3.get_distribution_logλ(st)
                t = BMT.bulk_microphysics_tendencies(
                    mode, BMT.Microphysics2Moment(), mp, tps,
                    ρ, T, q_tot, x..., logλ, Δt, nsub,
                )
                rates = collect(FT, _applied_2mp3(t))
                all(isfinite, rates) || break     # a non-finite substep is a different defect
                x = max.(x .+ Δt .* rates, zero(FT))
                # The host's own grid-mean rime constraint, applied here because this sweep claims
                # to test states the MODEL reaches: `enforce_grid_mean_microphysics_constraints!`
                # (ClimaAtmos, mass_flux_closures.jl) zeroes BOTH rime moments where there is no
                # ice mass or number. Without it the sweep carries states the model never does -
                # specifically the rounding residue a near-total drain leaves behind, whose
                # quotient is noise. See the ice-free residue assertion below.
                if x[5] <= ϵₘ || x[6] <= ϵₘ
                    if x[7] > resid
                        resid = x[7]
                        resid_by_number = x[5] > ϵₘ   # the NUMBER condition fired, not the mass
                    end
                    x[7] = zero(FT)
                    x[8] = zero(FT)
                end

                # the presence half of the invariant: never one moment without the other
                orphan = (x[7] > 0) != (x[8] > 0)
                n_orphan += orphan
                @test !orphan

                if x[7] > 0 && x[8] > 0
                    n_live += 1
                    ρ_bulk = x[7] / x[8]
                    worst_hi = max(worst_hi, ρ_bulk / ρ_max)
                    worst_lo = min(worst_lo, ρ_bulk / ρ_min)
                    @test ρ_bulk <= ρ_max * (1 + tol)
                    # Rime increments arrive along fixed-density rays, one per process, and each
                    # of them preserves the density interval. Their sum need not: at this state
                    # the net increment is a drain of effective density 710.7 acting on a pair at
                    # 539.2, and removing material denser than the pair leaves the remainder
                    # below the interval's lower bound for a single step, after which the
                    # interval holds again. The rime pair floor now projects the summed
                    # increment back onto the interval at the positivity floor, so this state
                    # no longer breaches it; the marker is retired per its own stated exit
                    # condition.
                    @test ρ_bulk >= ρ_min * (1 - tol)
                end
            end
        end
        # THE REGRESSION CASE: integration-tip trajectory 47, the same-state control that sits
        # over the ceiling at Float32 (921.6) and inside it at Float64 (900.3) on the identical
        # draw. Stepped here with the host constraint in force, which is what the model does, the
        # pair is zeroed with the ice and no out-of-interval state survives at either precision.
        @testset "the trajectory-47 state" begin
            ρ47, T47, q_tot47 = FT(0.8467667), FT(284.1068), FT(0.00075612945)
            x = FT[1.0124526e-4, 1.6376661e8, 5.0111406e-4, 41951.746,
                2.4247203e-7, 4729.9385, 8.541345e-8, 9.5069896e-11]
            @test ρ_min < x[7] / x[8] < ρ_max  # starting inside the interval
            # It ENTERS the newly admitted trace-ice band rather than starting in it: q_ice begins
            # at 2.03x eps(Float32) and drops below within one step. That entry is the precondition
            # the regression case exists to exercise, so it is asserted after the loop rather than
            # before it - asserting it on the initial state was simply wrong.
            entered_band = false
            for _ in 1:K
                st = P3.state_from_prognostic(p3, ρ47 * x[5], ρ47 * x[6], ρ47 * x[7], ρ47 * x[8])
                t = BMT.bulk_microphysics_tendencies(
                    mode, BMT.Microphysics2Moment(), mp, tps, ρ47, T47, q_tot47,
                    x..., P3.get_distribution_logλ(st), Δt, nsub)
                r = collect(FT, _applied_2mp3(t))
                all(isfinite, r) || break
                x = max.(x .+ Δt .* r, zero(FT))
                entered_band |= x[5] < eps(Float32)
                if x[5] <= ϵₘ || x[6] <= ϵₘ
                    # the residue the drain leaves is negligible in MASS, which is the claim that
                    # matters; its quotient is noise and the host zeroes the pair. The bound is the
                    # sweep's own live worst residue (below, re-measured the same way), not an
                    # independent guess - the eps/f floor scales with how coarse the step is at the
                    # crossing, which scales with the melt rate, so a rate change moves this floor.
                    # Re-measured after the melt-rate correction of 2026-08-06: trajectory 47's
                    # own worst residue is 1.11836655e-8 (Float32, first crossing step; Float64
                    # stays at 3.23e-12), against the sweep-wide worst of 1.0281e-7 (Float32)
                    # reported below - old bound 1e-9 was sized to the pre-fix rate.
                    @test x[7] < FT(2e-7)
                    x[7] = zero(FT)
                    x[8] = zero(FT)
                end
                @test (x[7] > 0) == (x[8] > 0)
                x[7] > 0 && @test ρ_min * (1 - tol) <= x[7] / x[8] <= ρ_max * (1 + tol)
            end
            @test entered_band   # or the case is not exercising the regime it was filed under
        end

        # the sweep has to have exercised rimed states, or it asserts nothing
        @test n_live > 100
        # A near-total drain leaves a rounding residue, and the RATIO of two such residues carries
        # no information: the ray form preserves the quotient to O(eps/f) when a step removes all
        # but a fraction f, so at f -> 0 there is no quotient to preserve. Measured on the
        # integration tip's trajectory 47 at Float32: the drain's implied density matched the bulk
        # to six digits (898.732 against 898.733) and still left q_rim = 4.45e-16 on a cell with
        # q_ice = 0, whose ratio to b_rim = 4.83e-19 is 921.6 - above solid ice by 0.5 percent,
        # which is exactly the eps/f floor at f = 3e-5. That state is not a rime population and the
        # host zeroes it; what has to be bounded is the residue's MASS, not its quotient.
        #
        # The residue is REPORTED rather than asserted, because it is a property of the HOST's
        # rule and not of this scheme, and the measurement says the obvious bound is false. At
        # Float32 the worst residue is 1.0281e-7 kg/kg, about 103 times `Q_RIM_RIMED`, the
        # threshold at which the event census counts a cell as rimed at all; at Float64 it is
        # 9.81e-12, still several orders above `eps`. (Re-measured 2026-08-06 after the melt-rate
        # correction moved the residue's magnitude; the mechanism and the earlier 2.5e-8 / 1.7e-12
        # pre-fix figures are otherwise unchanged from when this was first written.)
        #
        # `resid_by_number` reports WHICH host condition fired, and it is `false` at both
        # precisions: the ice MASS condition, not the number one. So the residue is rime left on a
        # cell whose ice mass has reached `eps`, and it is not bounded by that ice mass here,
        # because this sweep models only ONE of the host's two rime rules - the zeroing - and not
        # its `min(rho_q_rim, rho_q_ice)` bound. That is a gap in this harness rather than in the
        # scheme, and it sits with the host pair-breaker already flagged; a first draft of this
        # comment blamed the NUMBER condition and the instrument refuted it.
        #
        # Asserting a bound on the host's rule from here would be asserting what this branch cannot
        # guarantee, and loosening one until it passed would be tuning. Two earlier bounds were
        # both wrong: `1e-3 * eps(FT)`, off by seven orders because it came from a single
        # trajectory's numbers rather than from a scale, and `Q_RIM_RIMED`, which assumed the rule
        # only ever fires on trace mass.
        @info "rime-density composition sweep [FT=$FT]" n_live n_orphan worst_hi worst_lo tol resid resid_by_number
    end
end

test_rime_density_interval_is_forward_invariant(Float64)
test_rime_density_interval_is_forward_invariant(Float32)


# The two Tier-1 relaxation timescales are functions of the populations they relax, and the manual
# Jacobian differentiates them as such: `1/τ_i ∝ n_ice` from the P3 capacitance integral at frozen
# `logλ` (the number enters `log N₀` additively and the quadrature bounds are quantiles of the
# shape alone) and `1/τ_l ∝ N_lcl^{2/3} q_lcl^{1/3}` from the droplet diameter moment (the exponent
# is `(ν_c + 2) − μ_c (ν_cD + 2)/μ_cD = 2/3` for every SB2006 shape pair). Both were previously
# differentiated as constants, which left the ice-mass row without an ice-number entry, the
# cloud-mass row without a droplet-number entry, and the sublimation number pathway - degree-2
# homogeneous in `n_ice`, since its own rate carries a factor of `n_ice` through `1/τ_i` - with half
# its self-derivative.
#
# Every entry is checked against ForwardDiff differentiating the RATES at the consumer,
# `_per_process_2mp3`, so the degeneracy gates, the mean-mass bound and the capacitance quadrature
# are all inside the differentiated function; the manual entries are built by the recipe
# (homogeneity degree times rate over donor), not by symbolic differentiation of the closure. The
# values the code carried BEFORE are asserted to be wrong at the same states, so the testset fails
# if the entries are ever dropped back to zero or halved as well as if they are mis-derived.
function test_timescale_number_couplings(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    pdf_c = mp.warm_rain.seifert_beheng.pdf_c
    ρ = FT(0.9)
    tol = sqrt(eps(FT))
    # the Tier-1 donor floor, matching `_jacobian_2mp3_manual`'s own `n_floor = q_floor`; needed
    # here to subtract the ice-aggregation contribution to the ice-number diagonal
    n_floor = FT(TDI.TD.Parameters.q_min(tps))

    # The phase-change slots alone, differentiated in ONE state variable, with the value lane read
    # from THE SAME evaluation. Two reasons for the hand-seeded dual rather than `FD.jacobian`:
    # restricting the differentiated function to these slots is what makes the comparison valid
    # (the full tendency's ice-mass row also carries the mixed-phase collision transfer, whose
    # ice-side receivers the Tier-3 recipe deliberately leaves explicit), and carrying the value
    # lane is what separates the derivative from the primal-versus-dual evaluation difference
    # measured below. `Nothing` is a fine perturbation tag here: nothing differentiates through
    # these calls a second time.
    #
    # THE CONSEQUENCE, which is what a comparison against a `J` entry has to respect: these lanes
    # are the phase-change slots and ONLY those, so `der_*` is never comparable to a `J` entry that
    # other processes also write - the difference has to be subtracted first. It bit the ice-NUMBER
    # diagonal and not the ice-mass row, because ice aggregation is a number-only sink: it writes
    # slot 6 alone, so it lands on `J[6, 6]` and touches nothing this helper covers. Widening the
    # helper to include it would be the wrong repair, since it is a two-body collision rather than
    # a phase change and the temperature-coupled assertions below rely on exactly that distinction
    # when they scale by `Γᵢ`.
    part(u) = u isa FD.Dual ? FD.partials(u, 1) : zero(FT)
    function slot_lanes(T, q_tot, x, logλ, k)
        y = ntuple(i -> FD.Dual{Nothing}(FT(x[i]), FT(i == k)), 8)
        pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot, y..., logλ)
        s = Tuple(pp.cloud_condevap + pp.ice_depsub)
        return (map(FD.value, s), map(part, s))
    end

    logλ_of(x) = P3.get_distribution_logλ(
        P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]))
    # `s` is the vapor saturation ratio over ice asked for; the state's own condensate is added on
    # top, so the loading never moves the supersaturation the state is meant to sit at
    q_tot_at(T, x, s) =
        FT(s) * TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ) + x[1] + x[3] + x[5]
    # the same construction over liquid, which is how the cloud closure's growth branch is reached:
    # supersaturation over ice is not supersaturation over liquid at these temperatures
    q_tot_liq(T, x, s) =
        FT(s) * TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ) + x[1] + x[3] + x[5]

    # (q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim). Ordinary loadings: 0.1 g/kg of cloud
    # at 1e8 kg⁻¹, drizzle, and 1e-5 kg/kg of ice at 1e4 kg⁻¹, i.e. a 1e-9 kg mean crystal - inside
    # the ice number-adjustment bounds, so `_numadj_derivs` contributes nothing to row 6.
    ordinary = BMT.MicroState2MP3{FT}(1e-4, 1e8, 1e-5, 1e3, 1e-5, 1e4, 0, 0)
    # fresh crystals: a tenth of the mass spread over ten times the number, so the mean crystal is
    # 1e-11 kg (about 30 μm) - the small-crystal regime the coupling is largest in, and the one
    # transport and the number adjustment keep creating. Kept above eps(Float32) so the mixed-phase
    # block runs at both precisions.
    fresh = BMT.MicroState2MP3{FT}(1e-4, 1e8, 1e-5, 1e3, 1e-6, 1e5, 0, 0)
    # the same ordinary ice with half its mass rime at 400 kg/m³, so the thresholds are live
    rimed = BMT.MicroState2MP3{FT}(1e-4, 1e8, 1e-5, 1e3, 1e-5, 1e4, 5e-6, 5e-6 / 400)
    T = FT(253)

    @testset "the relaxation timescales carry their population couplings [FT=$FT]" begin
        for (label, x) in (("ordinary", ordinary), ("fresh ice", fresh), ("rimed", rimed))
            logλ = logλ_of(x)
            # the droplet number must be inside the mean-mass bound, or the condensation entry is
            # correctly zero and the agreement test below asserts nothing
            @test CM2.number_bounded_by_mass_limits(
                (; x_min = pdf_c.xc_min, x_max = pdf_c.xc_max), x[1], x[2]) == x[2]
            # the third case reaches the cloud closure's GROWTH branch; the first two leave it on
            # the evaporation-limited branch, since ice supersaturation at 253 K is liquid
            # subsaturation
            cases = (
                ("ice deposition", q_tot_at(T, x, FT(1.2)), 1),
                ("ice sublimation", q_tot_at(T, x, FT(0.8)), -1),
                ("liquid growth", q_tot_liq(T, x, FT(1.1)), 1),
            )
            for (branch, q_tot, ice_sign) in cases
                pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot, Tuple(x)..., logλ)
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
                J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
                (val_n, der_n) = slot_lanes(T, q_tot, x, logλ, 6)   # ∂/∂n_ice
                (val_c, der_c) = slot_lanes(T, q_tot, x, logλ, 2)   # ∂/∂n_lcl

                @testset "$label, $branch" begin
                    # both closures must be LIVE, or every entry below is a gate's zero
                    @test pp.ice_depsub.q_ice != 0
                    @test pp.cloud_condevap.q_lcl != 0
                    @test sign(pp.ice_depsub.q_ice) == ice_sign
                    if branch == "liquid growth"
                        @test pp.cloud_condevap.q_lcl > 0   # the vapor branch, by construction
                    end

                    # MEASURED, and it is why the ice entry is compared at two tolerances: the
                    # DIFFERENTIATED evaluation of an UNRIMED state does not reproduce the primal
                    # Seeding n_ice alone reproduces the primal rate EXACTLY, which is what makes
                    # the derivative comparisons below statements about the derivative and nothing
                    # else. Seeding the FULL state does not, at an unrimed state: `P3State` selects
                    # `D_gr`/`D_cr` on `iszero(F_rim)`, and with `q_rim` perturbed `F_rim` is a
                    # `Dual` whose value is zero and whose partials are not, which is not `iszero`,
                    # so that pass takes the RIMED threshold branch and integrates the capacitance
                    # over different segment boundaries (7.4e-6 at the ordinary state, 1.1e-6 at
                    # the fresh one, nothing at the rimed one - job 6925581, and measured again by
                    # `diag/taucoup_conditioning.jl`). That is a property of the exact-AD path,
                    # which differentiates exactly that full-state function, and it is pre-existing:
                    # `isunrimed` documents the identical trap one function over. It is left as a
                    # measurement rather than an assertion here, so that fixing it upstream does
                    # not fail this testset.
                    @test val_n[5] == pp.ice_depsub.q_ice

                    # (a) the ice-mass row's ice-number entry, identically zero before. Exact
                    # linearity first, both lanes from ONE evaluation.
                    @test der_n[5] ≈ val_n[5] / x[6] rtol = tol
                    @test der_n[5] != 0
                    # the manual entry is that same form on the PRIMAL rate, which is the rate `f`
                    # is built from, so f and J stay consistent with each other
                    @test J[5, 6] ≈ pp.ice_depsub.q_ice / x[6] rtol = tol
                    # and the two therefore agree, at the tolerance of the division alone
                    @test J[5, 6] ≈ der_n[5] rtol = tol
                    # SIGN, asserted explicitly so it is not later "fixed": on the deposition
                    # branch this is a POSITIVE OFF-DIAGONAL, because more crystals deposit more
                    # vapor. It is a growth coupling and belongs in the linearization;
                    # `ExplicitGrowthDiagonal` removes positive diagonals only.
                    @test sign(J[5, 6]) == sign(pp.ice_depsub.q_ice)

                    # (b) the cloud-mass row's droplet-number entry, identically zero before. The
                    # condensation closure is closed-form with no quadrature and no threshold
                    # branch, so its two evaluations agree and one tolerance does.
                    @test val_c[1] ≈ pp.cloud_condevap.q_lcl rtol = tol
                    @test der_c[1] ≈ 2 * val_c[1] / (3 * x[2]) rtol = tol
                    @test der_c[1] != 0
                    @test J[1, 2] ≈ der_c[1] rtol = tol
                    @test sign(J[1, 2]) == sign(pp.cloud_condevap.q_lcl)
                    # the 2/3 is the diameter moment's homogeneity degree and nothing else: a
                    # degree-1 reading would be 1.5x this entry
                    @test !isapprox(der_c[1], val_c[1] / x[2]; rtol = FT(0.1))

                    # (c) the sublimation number pathway's self-derivative, halved before
                    if ice_sign < 0
                        @test pp.ice_depsub.n_ice < 0
                        # Row 6 column 6 is a SUM, so every process reaching it has to be
                        # accounted for or the comparison is against the wrong quantity. These
                        # three are zero at this state and asserted so.
                        @test pp.ice_deposition.n_ice == 0
                        @test pp.ice_numadj.n_ice == 0
                        @test pp.ice_melting.n_ice == 0
                        # Ice aggregation is the fourth, and it CANNOT be asserted zero: it is a
                        # number sink driven by the ice population itself, so it is live wherever
                        # there is ice, which is every state this testset uses. It is measured and
                        # subtracted instead. Its Jacobian entry is the one `_jacobian_2mp3_manual`
                        # adds as `2·pp.ice_aggregation.n_ice·dnice`, verified in its own testset
                        # above by zeroing the slot.
                        agg = 2 * pp.ice_aggregation.n_ice / max(n_floor, x[6])
                        @test pp.ice_aggregation.n_ice < 0
                        @test agg != 0
                        # BOTH factors of two on this diagonal are HOMOGENEITY DEGREES, and that is
                        # the thing to carry forward when a fifth process is added here. The
                        # sublimation pathway is `∂ₜn = n·∂ₜq/q` with `∂ₜq` itself LINEAR in `n`
                        # (the timescale carries the population), so differentiating in `n` hits
                        # both factors and gives `2(n/q)·(∂ₜq/n)`. Aggregation is QUADRATIC in `n`
                        # (self-collection is a two-body rate), so `∂/∂n` of `−c·n²` is `2·rate/n`.
                        # Same exponent, different reason. Anything added to this diagonal needs
                        # its own degree worked out and subtracted here, or these assertions will
                        # fail exactly as they did when aggregation was overlooked - the manual
                        # entry was right and agreed with ForwardDiff of the composed tendency to
                        # the last bit on the rimed state, while the assertion was comparing a
                        # composed diagonal against one pathway.
                        @test der_n[6] ≈ 2 * val_n[5] / x[5] rtol = tol
                        # the degree-1 reading, which is what the code carried
                        @test !isapprox(der_n[6], val_n[5] / x[5]; rtol = FT(0.1))
                        @test J[6, 6] ≈ 2 * pp.ice_depsub.q_ice / x[5] + agg rtol = tol
                        # `der_n` differentiates `cloud_condevap + ice_depsub` alone, by
                        # construction in `slot_lanes`, so the aggregation term has to come off
                        # before the two are comparable at all
                        @test J[6, 6] - agg ≈ der_n[6] rtol = tol
                        # the pathway ties the two entries together: ∂ₜn = n·∂ₜq/q with ∂ₜq itself
                        # linear in n gives ∂/∂n = 2(n/q)·(∂ₜq/n)
                        @test J[6, 6] - agg ≈ 2 * (x[6] / x[5]) * J[5, 6] rtol = tol
                        # sublimation damps its own number, and the doubling keeps it bounded by
                        # 2/(τ_i·Γᵢ) rather than making it stiff. Aggregation damps too, so the
                        # sum is negative whether or not it is separated out.
                        @test J[6, 6] < 0 && isfinite(J[6, 6])
                        @test J[6, 6] - agg < 0
                    else
                        # deposition carries no number, so the pathway is off and its
                        # self-derivative with it; row 6 column 6 is then F23 and the number
                        # adjustment, neither of which this testset owns
                        @test pp.ice_depsub.n_ice == 0
                        @test der_n[6] == 0
                    end

                    # the rim rows keep NO number coupling: the rim drain's other three donors are
                    # explicit, so linearizing it in one of them alone would be the inconsistency,
                    # not the omission
                    @test J[7, 6] == 0 && J[8, 6] == 0
                    @test all(isfinite, J)

                    # card #23: `_jacobian_2mp3_manual` (`J` above) linearizes the BARE
                    # condensation/deposition rate directly now (Γ=1 at both its call sites), so
                    # the 9x9's species block - `_jacobian_2mp3t_manual`'s own `J8c`, built by
                    # calling `_jacobian_2mp3_manual` with these SAME `g`/`x`/`pp`/`rs` - is no
                    # longer a Γ-SCALED version of `J`, it IS `J`: both compute the identical
                    # closed form from the identical bare primal. Before this card, `J` was
                    # folded and the 9x9 corrected it by Γ per entry (the relationship these three
                    # assertions used to check); now there is nothing left to correct, so the
                    # entries agree exactly rather than up to a Γ factor.
                    y = BMT.MicroState2MP3T{FT}(Tuple(x)..., T)
                    gT = BMT.Temperature2MP3Tendency(mp, tps, ρ, q_tot, logλ)
                    ctx = BMT._phase_relaxation_context(mp, tps,
                        _micro_2mp3(q_tot, Tuple(x)...),
                        _thermo_2mp3(ρ, T, logλ))
                    J9 = BMT._jacobian_2mp3t_manual(gT, y, pp, rs, ctx)
                    @test J9[5, 6] == J[5, 6]
                    @test J9[1, 2] == J[1, 2]
                    if ice_sign < 0
                        # Same argument as [5,6]/[1,2] above, restated for the composed diagonal:
                        # `J9[6,6]` and `J[6,6]` are now the same computation on the same bare
                        # rate, so they agree exactly - `agg` (ice aggregation, live at every
                        # state with ice, no Γ-coupling of its own) drops out of the comparison
                        # rather than needing to be subtracted before scaling.
                        @test J9[6, 6] == J[6, 6]
                    end
                    @test all(isfinite, J9)
                end
            end
        end

        # An entry must be zero exactly where its rate is gated to zero (the f/J consistency
        # doctrine). ForwardDiff agrees at each of these states for the same reason the manual
        # entry does - the gate is an `ifelse` inside the differentiated function - so the two are
        # checked together rather than the manual value alone.
        @testset "the entries vanish exactly where their rates are gated" begin
            empty = BMT.MicroState2MP3{FT}(0, 0, 0, 0, 0, 0, 0, 0)
            logλ₀ = logλ_of(BMT.MicroState2MP3{FT}(0, 0, 0, 0, 1e-5, 1e4, 0, 0))
            q_tot = q_tot_at(T, empty, FT(1.3))

            # ice-free and droplet-free: both capacitance integrals underflow, both slots are
            # gated, and neither number coupling may survive
            pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot, Tuple(empty)..., logλ₀)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ₀)
            @test all(iszero, Tuple(pp.ice_depsub))
            @test all(iszero, Tuple(pp.cloud_condevap))
            J = BMT._jacobian_2mp3_manual(g, empty, pp, rs)
            (_, der_n) = slot_lanes(T, q_tot, empty, logλ₀, 6)
            (_, der_c) = slot_lanes(T, q_tot, empty, logλ₀, 2)
            # Both entries vanish here. J[5, 6] because the deposition slot is gated, and
            # J[1, 2] because every contributor to it is gated as well: the air is droplet-free
            # and supersaturated over liquid, which is where droplet activation would fire, and
            # this parameter set carries no aerosol, so the activation rate is identically zero
            # by the activation-off contract of the default construction. The state where that
            # entry is written by activation is asserted where activation arrives, against a
            # parameter set that supplies an aerosol.
            @test J[5, 6] == 0
            @test J[1, 2] == 0
            @test der_n[5] == 0 && der_c[1] == 0
            # Row 6 column 6 is NOT zero here and must not be asserted to be: this state is
            # supersaturated over ice, so the ice number adjustment relaxes an empty number at
            # −1/τ_numadj and deposition nucleation relaxes toward its target at the state-
            # dependent delivery rate (`CM_HetIce.delivery_rate`, no longer a stored constant -
            # see the file header). Both are their own entries, not this closure's; pinning row 6
            # to exactly their sum is the assertion that says the deposition pathway contributes
            # nothing.
            @test pp.ice_deposition.n_ice > 0
            S_i_here = TDI.supersaturation_over_ice(tps, q_tot, empty.q_lcl + empty.q_rai, empty.q_ice, ρ, T)
            inv_τ_dep_here = CM_HetIce.delivery_rate(mp.ice.ice_nucleation, mp, tps, T, S_i_here)
            @test J[6, 6] ≈
                  -inv_τ_dep_here - 1 / BMT._ice_numadj_params(p3).τ rtol =
                4 * eps(FT)
            @test der_n[6] == 0

            # above freezing with supersaturation over ice the entry caps deposition at zero, so
            # the ice-number coupling of a rate that is identically zero goes with it
            T_warm = TDI.T_freeze(tps) + FT(5)
            x = ordinary
            logλ = logλ_of(x)
            q_warm = q_tot_at(T_warm, x, FT(1.2))
            pp_w, rs_w = _per_process_2mp3_and_riming(mp, tps, ρ, T_warm, q_warm, Tuple(x)..., logλ)
            g_w = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T_warm, q_warm, logλ)
            @test pp_w.ice_depsub.q_ice == 0
            J_w = BMT._jacobian_2mp3_manual(g_w, x, pp_w, rs_w)
            @test J_w[5, 6] == 0
            @test slot_lanes(T_warm, q_warm, x, logλ, 6)[2][5] == 0

            # where the mean-mass bound clamps the droplet number the rate reads a number the
            # state no longer sets, so the primal has NO n_lcl sensitivity and the entry must
            # vanish with it. ForwardDiff returns exactly zero here, which is the check that this
            # is the clamp's own structure and not a threshold this test invented.
            n_hi = 10 * ordinary[1] / pdf_c.xc_min
            clamped = BMT.MicroState2MP3{FT}(
                ordinary[1], n_hi, ordinary[3], ordinary[4],
                ordinary[5], ordinary[6], ordinary[7], ordinary[8])
            logλ_c = logλ_of(clamped)
            q_c = q_tot_at(T, clamped, FT(1.2))
            pp_c, rs_c = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_c, Tuple(clamped)..., logλ_c)
            g_c = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_c, logλ_c)
            @test CM2.number_bounded_by_mass_limits(
                (; x_min = pdf_c.xc_min, x_max = pdf_c.xc_max), clamped[1], n_hi) < n_hi
            @test pp_c.cloud_condevap.q_lcl != 0     # the rate is live, only its donor is pinned
            @test BMT._jacobian_2mp3_manual(g_c, clamped, pp_c, rs_c)[1, 2] == 0
            @test slot_lanes(T, q_c, clamped, logλ_c, 2)[2][1] == 0
        end

        # Structure: the new couplings write only into the rows their own closure writes. On an
        # ice-and-vapor state every warm-rain and mixed-phase process is inactive for want of a
        # donor, so the two number columns are the phase-change closure alone.
        @testset "the couplings are confined to the phase-change rows" begin
            x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, 1e-5, 1e4, 0, 0)
            logλ = logλ_of(x)
            q_tot = q_tot_at(T, x, FT(0.8))
            pp, rs = _per_process_2mp3_and_riming(mp, tps, ρ, T, q_tot, Tuple(x)..., logλ)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
            J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
            @test pp.ice_depsub.q_ice < 0
            for i in 1:8
                # Column 2 is empty apart from the cloud number adjustment's own diagonal, which
                # relaxes an empty droplet number and is not this closure's entry: with no
                # droplets the condensation slot is gated and every warm-rain number pathway
                # vanishes with its donor.
                if i != 2
                    @test J[i, 2] == 0
                end
                # column 6 carries the deposition mass row and the sublimation number row, and
                # nothing else - the rim drain's number coupling is the deliberate omission
                if i != 5 && i != 6
                    @test J[i, 6] == 0
                end
            end
            @test J[5, 6] != 0 && J[6, 6] != 0
        end
    end
end

test_timescale_number_couplings(Float64)
test_timescale_number_couplings(Float32)
# `UT.guarded_quotient` and the trap it exists for, at the operation level.
#
# `ForwardDiff` orders `Dual`s LEXICOGRAPHICALLY. On a `Dual` whose value is zero and whose
# partials are live the value lane reads ABSENT and the `Dual` comparison reads PRESENT, so any
# guard built on the second runs the division it exists to avoid - however the quotient is then
# written. The three controls below keep the naked form, the half-fixed form and the magnitude
# floor red forever, which is the only way this stays fixed: each of them has the right VALUE and
# non-finite DERIVATIVES, so every value-checking suite in the set passes them.
function test_guarded_quotient_reads_the_value_lane(FT)
    @testset "guarded_quotient reads the value lane [FT=$FT]" begin
        a = FD.Dual{Nothing}(FT(1), FT(1), FT(0))
        b = FD.Dual{Nothing}(FT(0), FT(0), FT(1))     # value zero, LIVE seed

        # THE DISCRIMINATOR: the two predicates disagree on exactly this Dual.
        @test !(FD.value(b) > zero(FD.value(b)))      # value lane: absent
        @test b > zero(b)                             # Dual comparison: present - the trap

        g = UT.guarded_quotient(a, b)
        @test iszero(FD.value(g))
        @test all(iszero, FD.partials(g))

        # CONTROL 1, the naked select - what `ifelse(b > 0, a/b, zero(a/b))` reduces to
        @test !all(isfinite, FD.partials(ifelse(b > zero(b), a / b, zero(a / b))))
        # CONTROL 2, the HALF-FIXED form: a benign denominator selected on the same Dual predicate
        let bs = ifelse(b > zero(b), b, one(b))
            @test !all(isfinite, FD.partials(ifelse(b > zero(b), a / bs, zero(a / bs))))
        end
        # CONTROL 3, the magnitude floor: floatmin² underflows, so ∂b is weighted by -a/0
        @test !all(isfinite, FD.partials(a / max(b, floatmin(FT))))

        # a DEAD zero must also give exact zero, and the present case must be the exact quotient
        bd = FD.Dual{Nothing}(FT(0), FT(0), FT(0))
        @test iszero(FD.value(UT.guarded_quotient(a, bd)))
        @test all(iszero, FD.partials(UT.guarded_quotient(a, bd)))
        bp = FD.Dual{Nothing}(FT(2), FT(0), FT(1))
        @test UT.guarded_quotient(a, bp) === a / bp            # inert wherever present
        @test all(isfinite, FD.partials(UT.guarded_quotient(a, bp)))
        # plain (non-Dual) arguments go through unchanged
        @test UT.guarded_quotient(FT(3), FT(2)) === FT(1.5)
        @test UT.guarded_quotient(FT(3), FT(0)) === FT(0)
        # the non-zero `absent` answer: a budget with nothing to apportion does not bind
        @test UT.guarded_quotient(a, b, oftype(a / one(b), Inf)) == FT(Inf)
        @test all(iszero, FD.partials(UT.guarded_quotient(a, b, oftype(a / one(b), Inf))))
    end
end

test_guarded_quotient_reads_the_value_lane(Float64)
test_guarded_quotient_reads_the_value_lane(Float32)
# The freezing composition must be differentiable where NO NUCLEATION PATHWAY IS OPEN, and this
# testset exists because my own degenerate-state falsifiers checked VALUES and never partials.
#
# The per-size form writes the per-drop rate as `1/(1/(J (π/6) D³) + τ_freeze)`, which is the
# right shape at an overflowed `J` but has a guarded value and UNGUARDED partials at a vanishing
# one: `1/0` is `Inf`, so `r` is correctly zero, while `d(1/x)/dx = -∂x/x²` is `0/0`. THE
# POPULATION IS ABSENT is the physical state family that reaches it: `D` is zero, so
# `J (π/6) D³` is zero, and the rain fall-speed law `α x^β` additionally has an infinite
# derivative at `x = 0`.
#
# REPAIRED FROM the campaign version, which also carried a SECOND family - "THE INP BUDGET IS
# SPENT", `immersion_limit_rate` returning exactly zero once `n_active ≥ INPC(T)/ρ` on the
# capped heterogeneous coefficient - that no longer applies: `cloud_freezing_rate` is Bigg alone
# now, with no INP budget or cap above it at all (see
# `test_cloud_freezing_is_uniform_with_rain`'s header for the fuller account of the same
# removed mechanism). Three testsets that tested that budget/cap directly are DROPPED, not
# repaired: "cloud: a spent INP budget gives exactly zero, partials included" and "the cap is a
# presence gate, not a magnitude floor" test the removed cap's own behavior end to end; the
# `immersion_limit_rate`/`n_active` precondition check inside "the ad_compat cloud-edge state,
# through the production entry" is dropped with it, keeping that testset's own AD-safety check
# (`edge_rhs`/`FD.jacobian`), which does not depend on the cap. The two population-presence
# testsets ("both categories: an absent population...", "above the freezing gate, both
# categories...") are REPAIRED, not dropped: they test the presence/temperature gate, which is
# unaffected, and needed only the removed `τ_act` keyword dropped from their
# `cloud_freezing_rate` calls.
#
# Gated on PRESENCE rather than floored in magnitude: a floor would make the rate small instead
# of absent and would hand the state a constant's partials, which is a wrong derivative rather
# than a missing one. Asserted at both precisions, on values and partials together.
function test_freezing_is_differentiable_with_no_pathway_open(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    pdf_c = mp.ice.cloud_pdf
    pdf_r = mp.ice.rain_pdf
    evap = mp.warm_rain.seifert_beheng.evap
    aps = mp.warm_rain.air_properties
    hom = mp.ice.homogeneous
    T_frz = TDI.T_freeze(tps)

    allfinite(x) = isfinite(FD.value(x)) && all(isfinite, FD.partials(x))
    exactzero(x) = iszero(FD.value(x)) && all(iszero, FD.partials(x))
    seed(v, i, n) = FD.Dual{:pathway}(v, ntuple(k -> FT(k == i), n)...)

    @testset "freezing differentiates where no pathway is open [FT=$FT]" begin
        @testset "the predicate must read the VALUE LANE (operation level)" begin
            # THE DISCRIMINATOR, asserted rather than inferred, and the reason a guard that looks
            # right can be wrong. `ForwardDiff` orders `Dual`s LEXICOGRAPHICALLY, so on a Dual
            # whose value is zero and whose partials are live the two predicates DISAGREE: the
            # value lane reads absent, the Dual comparison reads present. Any guard built on the
            # second runs the division it exists to avoid, however the quotient is then written.
            let a = FD.Dual{Nothing}(FT(1), FT(1), FT(0)), b = FD.Dual{Nothing}(FT(0), FT(0), FT(1))
                @test !(FD.value(b) > zero(FD.value(b)))   # value lane: absent
                @test b > zero(b)                          # Dual comparison: present - the trap

                g = UT.guarded_quotient(a, b)
                @test iszero(FD.value(g))
                @test all(iszero, FD.partials(g))

                # CONTROL 1, the naked select: same predicate, so the guard leaks. This is what a
                # bare `ifelse(b > 0, a/b, zero(a/b))` reduces to, and it must stay red forever.
                @test !all(isfinite, FD.partials(ifelse(b > zero(b), a / b, zero(a / b))))
                # CONTROL 2, the HALF-FIXED form: substituting a benign denominator is not enough
                # on its own, because the substitution is selected on the same Dual predicate.
                let bs = ifelse(b > zero(b), b, one(b))
                    @test !all(isfinite, FD.partials(ifelse(b > zero(b), a / bs, zero(a / bs))))
                end
                # CONTROL 3, the magnitude floor: right value, non-finite derivatives, because the
                # quotient rule weights ∂b by -a/b² and floatmin² underflows to exactly zero.
                @test !all(isfinite, FD.partials(a / max(b, floatmin(FT))))

                # ...and it is still the exact quotient wherever the denominator is present
                bp = FD.Dual{Nothing}(FT(2), FT(0), FT(1))
                @test FD.value(UT.guarded_quotient(a, bp)) == FT(0.5)
                @test all(isfinite, FD.partials(UT.guarded_quotient(a, bp)))
                @test UT.guarded_quotient(a, bp) === a / bp     # inert where present
                # the non-zero `absent` answer, for a budget with nothing to apportion
                @test UT.guarded_quotient(a, b, oftype(a / one(b), Inf)) == FT(Inf)
            end
        end

        @testset "a LIVE zero coefficient through the shared composition" begin
            # The state my first fix would have failed on, and the reason the guard is not merely
            # a select: `J` with value zero and LIVE partials. It is reachable whenever the INP
            # budget lands exactly at `n_active`, where `max(zero, INPC - n_active)` returns the
            # live branch rather than the constant zero. Driven at the operation level because a
            # trajectory state that hits the tie exactly is measure-zero to construct.
            ρ, T = FT(0.9), FT(263)
            qᵥ = TDI.p2q(tps, T, ρ, TDI.saturation_vapor_pressure_over_liquid(tps, T))
            Dr = FT(3e-4)
            n = FT(250)
            for J in (FD.Dual{Nothing}(FT(0), FT(0), FT(1)),   # live zero  <- the trap
                FD.Dual{Nothing}(FT(0), FT(0), FT(0)),         # dead zero
                FD.Dual{Nothing}(FT(1e5), FT(1), FT(0)))       # present
                m = CM_HetIce._composed_liquid_freezing_moments(
                    p3.vent, aps, tps, T, ρ, qᵥ, pdf_r.ρw, J, n,
                    u -> Dr * u, (D, x) -> evap.α * x^evap.β * sqrt(evap.ρ0 / ρ),
                    CM_HetIce.rain_freezing_quadrature(),
                )
                @test allfinite(m.∂ₜn_frz)
                @test allfinite(m.∂ₜq_frz)
                iszero(FD.value(J)) && @test exactzero(m.∂ₜn_frz)
                iszero(FD.value(J)) && @test exactzero(m.∂ₜq_frz)
            end
            # and a LIVE zero DIAMETER, where the rain fall-speed law α x^β additionally has an
            # infinite derivative - the second hazard the node substitution covers
            for Dr0 in (FD.Dual{Nothing}(FT(0), FT(0), FT(1)), FD.Dual{Nothing}(FT(0), FT(0), FT(0)))
                m = CM_HetIce._composed_liquid_freezing_moments(
                    p3.vent, aps, tps, T, ρ, qᵥ, pdf_r.ρw,
                    FD.Dual{Nothing}(FT(1e5), FT(1), FT(0)), n,
                    u -> Dr0 * u, (D, x) -> evap.α * x^evap.β * sqrt(evap.ρ0 / ρ),
                    CM_HetIce.rain_freezing_quadrature(),
                )
                @test exactzero(m.∂ₜn_frz)
                @test exactzero(m.∂ₜq_frz)
            end
        end

        @testset "the ad_compat cloud-edge state, through the production entry" begin
            # The regime verbatim, including the `n_active = n_ice` the call site passes and my
            # earlier probe defaulted to zero - which is exactly why this went unnoticed.
            ρ, T, q_tot = FT(0.7), FT(263), FT(0.005)
            x0 = FT[1e-5, 1e7, 1e-6, 1e3, 3e-8, 30, 1e-8, 2.5e-11]
            st = P3.state_from_prognostic(p3, ρ * x0[5], ρ * x0[6], ρ * x0[7], ρ * x0[8])
            logλ = P3.get_distribution_logλ(st)
            function edge_rhs(x)
                t = BMT.bulk_microphysics_tendencies(
                    BMT.Microphysics2Moment(), mp, tps, ρ, T, q_tot,
                    x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], logλ)
                return [t.dq_lcl_dt, t.dn_lcl_dt, t.dq_rai_dt, t.dn_rai_dt,
                    t.dq_ice_dt, t.dn_ice_dt, t.dq_rim_dt, t.db_rim_dt]
            end
            @test all(isfinite, edge_rhs(x0))
            @test all(isfinite, FD.jacobian(edge_rhs, x0))
            @test CM_HetIce.homogeneous_freezing_rate_coefficient(hom, tps, T) === zero(FT)
        end

        @testset "both categories: an absent population gives exactly zero, partials included" begin
            ρ, T = FT(0.9), FT(250)
            qᵥ = TDI.p2q(tps, T, ρ, TDI.saturation_vapor_pressure_over_liquid(tps, T))
            for (q_s, n_s) in ((FT(0), FT(0)), (FT(2e-4), FT(0)), (FT(0), FT(5000)))
                r = CM_HetIce.rain_freezing_rate(
                    mp.ice.rain_freezing, hom, p3.vent, aps, tps, evap, pdf_r,
                    seed(q_s, 1, 2), ρ, seed(n_s * ρ, 2, 2), T, qᵥ,
                )
                @test exactzero(r.∂ₜn_frz)
                @test exactzero(r.∂ₜq_frz)
                c = CM_HetIce.cloud_freezing_rate(
                    mp.ice.rain_freezing, hom, p3.vent, aps, tps, pdf_c,
                    seed(q_s, 1, 2), ρ, seed(n_s * ρ, 2, 2), T, qᵥ,
                )
                @test exactzero(c.∂ₜn_frz)
                @test exactzero(c.∂ₜq_frz)
            end
        end

        @testset "above the freezing gate, both categories, partials included" begin
            ρ = FT(0.9)
            for T in (FT(280), T_frz)
                qᵥ = TDI.p2q(tps, T, ρ, TDI.saturation_vapor_pressure_over_liquid(tps, T))
                r = CM_HetIce.rain_freezing_rate(
                    mp.ice.rain_freezing, hom, p3.vent, aps, tps, evap, pdf_r,
                    seed(FT(2e-4), 1, 2), ρ, seed(FT(5000 * 0.9), 2, 2), T, qᵥ,
                )
                @test exactzero(r.∂ₜn_frz)
                @test exactzero(r.∂ₜq_frz)
                c = CM_HetIce.cloud_freezing_rate(
                    mp.ice.rain_freezing, hom, p3.vent, aps, tps, pdf_c,
                    seed(FT(1e-4), 1, 2), ρ, seed(FT(1e8 * 0.9), 2, 2), T, qᵥ,
                )
                @test exactzero(c.∂ₜn_frz)
                @test exactzero(c.∂ₜq_frz)
            end
        end
    end
end

test_freezing_is_differentiable_with_no_pathway_open(Float64)
test_freezing_is_differentiable_with_no_pathway_open(Float32)

# that function rely on the returned rate already being non-positive above `T_freeze`. Asserted on
# a grid in `(T, S_i)` straddling `T_freeze`, with `S_i` set through `q_tot` against the ice
# saturation specific content at each temperature. `===` rather than `==` because `min` treats
# signed zeros specially.
function test_above_freezing_deposition_limiter_is_single_copy(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    T_frz = FT(TDI.TD.Parameters.T_freeze(tps))
    ρ = FT(0.6)
    q_ice, n_ice = FT(1e-5), FT(1e4)
    st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, FT(0), FT(0))
    logλ = P3.get_distribution_logλ(st)

    ΔTs = FT[-10, -3, -1, -FT(0.1), 0, FT(0.1), 1, 3, 10]
    S_is = FT[-FT(0.5), -FT(0.05), 0, FT(0.05), FT(0.5)]

    @testset "the above-freezing deposition limiter has one copy [FT=$FT]" begin
        @testset "the limiter fires on the growth branch above freezing and nowhere else" begin
            for ΔT in ΔTs, rate in FT[-1, -eps(FT), 0, eps(FT), 1]
                T = T_frz + ΔT
                fires = CMNonEq.wet_surface_deposition_limiter(rate, tps, T)
                @test fires == (T > T_frz && rate > 0)
            end
        end

        @testset "min(rate, 0) is the identity on what the shared function returns" begin
            for ΔT in ΔTs, S_i in S_is
                T = T_frz + ΔT
                qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
                q_tot = (1 + S_i) * qᵥ_sat_ice + q_ice
                τ_dep = P3.ice_deposition_timescale(
                    mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps,
                    T, ρ, st, logλ; quad = mp.ice.quad)
                micro = (; q_tot, q_lcl = FT(0), q_icl = q_ice, q_rai = FT(0), q_sno = FT(0))
                # Through the shared inner conversion, for the same reason as the liquid
                # side above: `ConstantTimescale` is a zero-field marker and the timescale
                # reaches the conversion through `mp`, so a test wanting a SPECIFIC one calls
                # what the public method forwards to.
                rate = CMNonEq._conv_q_vap_to_q_icl_const(τ_dep, tps, micro, (; ρ, T))
                @test isfinite(rate)
                if T > T_frz
                    @test min(rate, zero(T)) === rate
                    @test rate ≤ 0
                elseif S_i > 0
                    @test rate > 0
                end
            end
        end

        @testset "the per-process decomposition and the full entry agree across T_freeze" begin
            x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, q_ice, n_ice, 0, 0)
            for ΔT in ΔTs, S_i in S_is
                T = T_frz + ΔT
                qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
                q_tot = (1 + S_i) * qᵥ_sat_ice + q_ice
                pp = _per_process_2mp3(mp, tps, ρ, T, q_tot, Tuple(x)..., logλ)
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
                full = Tuple(g(x))
                psum = Tuple(sum(values(pp)))
                @test all(isfinite, full)
                @test all(isfinite, psum)
                scale = maximum(abs, psum) + eps(FT)
                rtol = FT == Float64 ? FT(1e-10) : FT(1e-4)
                @test maximum(abs.(full .- psum)) ≤ rtol * scale
                if T > T_frz
                    @test pp.ice_depsub.q_ice ≤ 0
                end
            end
        end
    end
end

test_above_freezing_deposition_limiter_is_single_copy(Float64)
test_above_freezing_deposition_limiter_is_single_copy(Float32)

# The temperature-coupled Jacobian took the melting rate's temperature column as
# `rate / max(1e-3, T - T_freeze)`. The conduction-only integral carries no temperature
# dependence at all, so the derivative is available exactly and without the division: `ice_melt`
# returns it alongside the rate, and `_ice_melting_species` is linear in the three rates. The
# retired floor made two separate errors, and they are asserted separately below: it understated
# the column by `ΔT / 1e-3` inside the millikelvin layer, and it omitted the `L_f(T)` term
# everywhere, which is `(cp_l - cp_i)/L_f ≈ 6e-3` per K of excess.
function test_melting_temperature_column_is_exact(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    vel = mp.ice.terminal_velocity
    aps = mp.warm_rain.air_properties
    quad = mp.ice.quad
    T_frz = FT(TDI.TD.Parameters.T_freeze(tps))
    ∂L_f_∂T = FT(TDI.TD.Parameters.cp_l(tps) - TDI.TD.Parameters.cp_i(tps))
    ρ = FT(0.6)
    q_ice, n_ice = FT(1e-5), FT(1e4)
    st = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, FT(0), FT(0))
    logλ = P3.get_distribution_logλ(st)
    # supersaturated over ice above freezing, so the deposition column is inactive there and the
    # melting contribution owns the temperature column
    q_tot = FT(2e-2)
    x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, q_ice, n_ice, 0, 0)

    # Float32 cannot resolve a microkelvin excess against 273 K, so the grid is filtered rather
    # than shortened, and the smallest resolved excess is still inside the retired floor's layer.
    ΔTs = FT[1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10]
    resolved = filter(ΔT -> (T_frz + ΔT) - T_frz > 0, ΔTs)
    floor_ΔT = FT(1e-3)

    melt_at(T) = P3.ice_melt(vel, aps, tps, T, ρ, st, logλ; quad)

    @testset "the melting temperature column is exact [FT=$FT]" begin
        @test minimum(resolved) < floor_ΔT

        @testset "the derivative is the rate over the excess, corrected for L_f(T)" begin
            for ΔT in resolved
                T = T_frz + ΔT
                ΔT_eff = T - T_frz
                lf_corr = 1 - ΔT_eff * ∂L_f_∂T / TDI.Lf(tps, T)
                melt = melt_at(T)
                @test melt.dLdt > 0
                @test melt.∂dLdt_∂T ≈ melt.dLdt / ΔT_eff * lf_corr rtol = sqrt(eps(FT))
                @test melt.∂dNdt_∂T ≈ melt.dNdt / ΔT_eff * lf_corr rtol = sqrt(eps(FT))
            end
        end

        @testset "and agrees with ForwardDiff of ice_melt" begin
            for ΔT in resolved
                T = T_frz + ΔT
                melt = melt_at(T)
                dL = FD.derivative(t -> melt_at(t).dLdt, T)
                dN = FD.derivative(t -> melt_at(t).dNdt, T)
                @test isfinite(dL) && isfinite(dN)
                @test melt.∂dLdt_∂T ≈ dL rtol = sqrt(eps(FT))
                @test melt.∂dNdt_∂T ≈ dN rtol = sqrt(eps(FT))
            end
        end

        @testset "at and below freezing the rates and their derivatives are zero" begin
            for ΔT in FT[-10, -1, -1e-3, 0]
                T = T_frz + ΔT
                ctx = BMT._phase_relaxation_context(mp, tps,
                    _micro_2mp3(q_tot, Tuple(x)...),
                    _thermo_2mp3(ρ, T, logλ))
                @test all(iszero, Tuple(ctx.∂melting_∂T))
            end
        end

        @testset "the column agrees with ForwardDiff of the melting contribution" begin
            for ΔT in resolved
                T = T_frz + ΔT
                ctx = BMT._phase_relaxation_context(mp, tps,
                    _micro_2mp3(q_tot, Tuple(x)...),
                    _thermo_2mp3(ρ, T, logλ))
                ad = FD.derivative(T) do t
                    m = melt_at(t)
                    SVector{8}(Tuple(BMT._ice_melting_species(
                        ρ, p3.ρ_i, zero(FT), zero(FT),
                        m.dNdt, m.dLdt, m.melt_frac)))
                end
                col = SVector{8}(Tuple(ctx.∂melting_∂T))
                @test all(isfinite, ad)
                @test all(isfinite, col)
                @test maximum(abs.(col .- ad)) ≤ sqrt(eps(FT)) * (maximum(abs, ad) + eps(FT))
            end
        end

        @testset "the two errors the retired floor made, one at a time" begin
            for ΔT in resolved
                T = T_frz + ΔT
                ΔT_eff = T - T_frz
                lf_corr = 1 - ΔT_eff * ∂L_f_∂T / TDI.Lf(tps, T)
                pp = _per_process_2mp3(mp, tps, ρ, T, q_tot, Tuple(x)..., logλ)
                ctx = BMT._phase_relaxation_context(mp, tps,
                    _micro_2mp3(q_tot, Tuple(x)...),
                    _thermo_2mp3(ρ, T, logλ))
                linear = pp.ice_melting.q_ice / ΔT_eff
                floored = pp.ice_melting.q_ice / max(floor_ΔT, ΔT_eff)
                understated = min(ΔT_eff / floor_ΔT, one(FT))
                @test ctx.∂melting_∂T.q_ice < 0
                @test floored ≈ linear * understated rtol = sqrt(eps(FT))
                @test ctx.∂melting_∂T.q_ice ≈ linear * lf_corr rtol = sqrt(eps(FT))
                if ΔT_eff ≥ floor_ΔT
                    @test understated == 1
                end
            end
        end

        @testset "the 9x9 assembly is finite across T_freeze and carries the column" begin
            g = BMT.Temperature2MP3Tendency(mp, tps, ρ, q_tot, logλ)
            for ΔT in vcat(FT[-10, -1, 0], resolved)
                T = T_frz + ΔT
                y = BMT.MicroState2MP3T{FT}(Tuple(x)..., T)
                f9, J9 = BMT._tendency_and_jacobian(
                    BMT.TemperatureCoupledJacobian(), g, y)
                @test all(isfinite, f9)
                @test all(isfinite, J9)
                ctx = BMT._phase_relaxation_context(mp, tps,
                    _micro_2mp3(q_tot, Tuple(x)...),
                    _thermo_2mp3(ρ, T, logλ))
                # above freezing the deposition column is inactive and there is no liquid, so
                # the species rows of the temperature column are the melting contribution alone
                if T > T_frz
                    @test ctx.sat_excess_i > 0
                    for i in 1:8
                        @test J9[i, 9] == ctx.∂melting_∂T[i]
                    end
                end
            end
        end
    end
end

test_melting_temperature_column_is_exact(Float64)
test_melting_temperature_column_is_exact(Float32)
# The rain self-collection/breakup pair is the one process whose rate legitimately changes sign
# through the POPULATION'S shape rather than a thermodynamic attractor: breakup dominates above
# the equilibrium diameter and collection below it. The donor-linearization recipe applied to the
# pair written as a corrected sink therefore produces a POSITIVE diagonal entry - anti-damping -
# exactly where breakup wins. Written as the relaxation it already is, toward n_eq = q_rai/x_eq,
# the frozen-shape entry is -1/τ_eff and is non-positive by construction.
#
# This testset fails in both directions: it asserts the recipe it replaces IS positive in that
# band (otherwise the reformulation is unmotivated) and that the relaxation entry is not.
function test_rain_relaxation_jacobian_sign(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    td = CP.create_toml_dict(FT)
    win = CMP.RainParticlePDF_SB2006_windowed(td)
    mp = CMP.Microphysics2MParams(td; with_ice = true, rain_pdf = win)
    sb = mp.warm_rain.seifert_beheng
    ρ = FT(1)
    T = FT(283)
    q_tot = FT(5e-3)
    empty_state = P3.state_from_prognostic(mp.ice.scheme, FT(0), FT(0), FT(0), FT(0))
    logλ₀ = P3.get_distribution_logλ(empty_state)
    x_eq = FT(π) / 6 * win.ρw * sb.brek.Deq^3

    @testset "rain breakup does not anti-damp the Jacobian diagonal [FT=$FT]" begin
        # ratios of n_rai to n_eq: below 1 the mean drop is larger than equilibrium and breakup
        # dominates, which is the anti-damping band
        seen_positive_recipe = false
        for q_rai in FT[1e-6, 1e-5, 1e-4, 1e-3],
            r in FT[0.05, 0.2, 0.5, 0.8, 1.2, 2.0, 10.0, 100.0]

            n_rai = r * q_rai / x_eq
            x = BMT.MicroState2MP3{FT}(0, 0, q_rai, n_rai, 0, 0, 0, 0)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ₀)
            pp, rs = _per_process_2mp3_and_riming(
                mp, tps, ρ, T, q_tot, Tuple(x)..., logλ₀)
            J = BMT._jacobian_2mp3_manual(g, x, pp, rs)

            # row 4 is the rain-number row; its diagonal must never be a growth direction
            @test isfinite(J[4, 4])
            @test J[4, 4] <= 0

            # the recipe this replaces, on the same state: (sc + br)/n_rai
            recipe = (pp.rain_selfcol.n_rai + pp.rain_breakup.n_rai) / n_rai
            recipe > 0 && (seen_positive_recipe = true)
        end
        # the defect being fixed is real on this grid, not hypothetical
        @test seen_positive_recipe
    end
end

test_rain_relaxation_jacobian_sign(Float64)
test_rain_relaxation_jacobian_sign(Float32)

# THE PHANTOM CORNER: mass with no number. The mean-mass bound used to relax the number towards
# `q/x_max` there, which is positive at every positive mass, so it MANUFACTURED a population at the
# largest size the window allows - for cloud droplets 2.6e-10 kg, a 790 μm drop, four times the
# diameter at which the scheme calls something rain - and fed it to the freezing and collision
# integrals. The invention is deleted for the warm phase: the rates read the cell empty, and its
# orphan mass is removed as mass, by evaporation at the rate a population of minimum-mass droplets
# would evaporate at, and only where the air is subsaturated.
function test_phantom_corner_invents_nothing(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    sb = mp.warm_rain.seifert_beheng
    aps = mp.warm_rain.air_properties
    pdf_c = sb.pdf_c
    pdf_r = sb.pdf_r
    win_c = (; x_min = pdf_c.xc_min, x_max = pdf_c.xc_max)
    win_r = (; x_min = pdf_r.xr_min, x_max = pdf_r.xr_max)
    par_c = (; sb.numadj.τ, win_c...)
    τ = FT(sb.numadj.τ)
    ρ, T = FT(0.9), FT(285)
    q = FT(1e-6)

    @testset "the phantom corner invents no population [FT=$FT]" begin
        # (a) THE INVENTION IS GONE, both halves, both warm species, and the old value is named
        # here so the assertion is against the behaviour and not against itself
        @test CM2.number_bounded_by_mass_limits(win_c, q, FT(0); invent_from_zero = false) === zero(FT)
        @test CM2.number_tendency_from_mass_limits(par_c, q, FT(0); invent_from_zero = false) === zero(FT)
        @test CM2.number_bounded_by_mass_limits(win_r, q, FT(0); invent_from_zero = false) === zero(FT)
        # what it used to be: a population at the window's largest mass
        @test CM2.number_bounded_by_mass_limits(win_c, q, FT(0)) === q / FT(pdf_c.xc_max)
        @test CM2.number_bounded_by_mass_limits(win_c, q, FT(0)) > 0

        # (b) ICE IS UNTOUCHED - it takes the default, and its corner doctrine is its own decision
        p3 = mp.ice.scheme
        ice = BMT._ice_numadj_params(p3)
        @test CM2.number_tendency_from_mass_limits(
            (; ice.τ, ice.x_min, ice.x_max), q, FT(0)) > 0

        # (c) BIT-IDENTICAL wherever a number exists. The argument may not perturb any populated
        # cell, and the default may not perturb a caller that does not pass it.
        for qq in FT[1e-9, 1e-6, 1e-4, 1e-3], n in FT[1e2, 1e8, 1e12]
            @test CM2.number_bounded_by_mass_limits(win_c, qq, n; invent_from_zero = false) ===
                  CM2.number_bounded_by_mass_limits(win_c, qq, n)
            @test CM2.number_tendency_from_mass_limits(par_c, qq, n; invent_from_zero = false) ===
                  CM2.number_tendency_from_mass_limits(par_c, qq, n)
        end
        # and the identity that ties the two halves together survives the new arm
        for qq in FT[0, 1e-9, 1e-4], n in FT[0, 1e2, 1e8]
            @test CM2.number_tendency_from_mass_limits(par_c, qq, n; invent_from_zero = false) ===
                  (CM2.number_bounded_by_mass_limits(win_c, qq, n; invent_from_zero = false) - n) / τ
        end

        # (d) THE 790 MICRON STATE IS DEAD BY ASSERTION. Every rate that used to receive the
        # manufactured population receives zero instead, because it receives zero population.
        N_invented = ρ * CM2.number_bounded_by_mass_limits(win_c, q, FT(0))
        x̄_invented = ρ * q / N_invented
        @test x̄_invented ≈ FT(pdf_c.xc_max) rtol = sqrt(eps(FT))
        # MEASURED, and it corrects a number this campaign has been quoting: the invented droplet
        # is 79 μm across, not the 790 μm the discussion-12 brief and the ruling both carry. The
        # brief's figure is a factor of ten out - `cbrt(6·2.6e-10/(π·1000))` is 7.92e-5 m - and the
        # assertion is written against the arithmetic rather than against the quoted value, which
        # is how the discrepancy surfaced. The physical point is unchanged and is sharper stated
        # correctly: `xc_max` IS `xr_min`, so the manufactured particle sits exactly ON the
        # cloud-rain boundary, a drizzle drop conjured inside the cloud category.
        D_invented = cbrt(6 * x̄_invented / (FT(π) * FT(pdf_c.ρw)))
        @test FT(7e-5) < D_invented < FT(9e-5)
        @test FT(pdf_c.xc_max) == FT(pdf_r.xr_min)
        N_now = ρ * CM2.number_bounded_by_mass_limits(win_c, q, FT(0); invent_from_zero = false)
        @test N_now === zero(FT)
        bigg = CM_HetIce.liquid_freezing_rate(
            mp.ice.rain_freezing, pdf_c, tps, q, ρ, N_now, FT(250))
        @test bigg.∂ₜn_frz == 0 && bigg.∂ₜq_frz == 0
        acnv = CM2.autoconversion(sb.acnv, pdf_c, q, FT(0), ρ, N_now)
        @test all(iszero, (acnv.dq_lcl_dt, acnv.dq_rai_dt, acnv.dN_lcl_dt, acnv.dN_rai_dt))
        accr = CM2.accretion(sb, q, q, ρ, N_now)
        @test accr.dq_lcl_dt == 0 && accr.dq_rai_dt == 0
        # and the same rates ARE live on the invented population, so these assertions are
        # measuring the deletion and not a state where everything happens to be zero
        bigg_old = CM_HetIce.liquid_freezing_rate(
            mp.ice.rain_freezing, pdf_c, tps, q, ρ, N_invented, FT(250))
        @test bigg_old.∂ₜn_frz != 0

        # (e) THE ORPHAN DRAIN: subsaturated only, a relaxation, and derived rather than tuned
        wet, dry = FT(1e-4), FT(-1e-4)
        # `iszero` and not `=== zero(FT)`: the rate is `-q * inv_τ`, so a zero inverse timescale
        # returns NEGATIVE zero. That is numerically inert - it sums and compares as zero - but it
        # is not bit-identical to `+0.0`, and asserting identity here would be asserting the sign
        # of a zero rather than the absence of a rate.
        @test iszero(CM2.orphan_mass_drain(aps, tps, T, ρ, q, pdf_c.xc_min, pdf_c.ρw, wet))
        @test iszero(CM2.orphan_mass_drain(aps, tps, T, ρ, q, pdf_c.xc_min, pdf_c.ρw, zero(FT)))
        drain = CM2.orphan_mass_drain(aps, tps, T, ρ, q, pdf_c.xc_min, pdf_c.ρw, dry)
        @test drain < 0 && isfinite(drain)
        # linear in the mass, so the timescale is a property of the air and the minimum size and
        # not of how much orphan mass there happens to be
        inv_τ = CM2.orphan_mass_inv_timescale(aps, tps, T, ρ, pdf_c.xc_min, pdf_c.ρw, dry)
        @test drain ≈ -q * inv_τ rtol = sqrt(eps(FT))
        @test CM2.orphan_mass_drain(aps, tps, T, ρ, 3 * q, pdf_c.xc_min, pdf_c.ρw, dry) ≈
              3 * drain rtol = sqrt(eps(FT))
        # against the closure it is derived from: the same rate the condensation timescale gives
        # at the number that makes the mean mass exactly x_min
        N_xmin = ρ * q / FT(pdf_c.xc_min)
        τ_cl = CM2.cloud_condensation_timescale(pdf_c, aps, tps, T, ρ, q, N_xmin)
        @test !CM2.cloud_condensation_is_degenerate(τ_cl)
        # MEASURED, against a prediction of this testset's own that was wrong and is worth keeping
        # because it bears on the ruling's stated intent. The drain was specified as a relaxation
        # rather than an instant removal so that a one-step flicker into the corner would cost
        # `dt/τ` of the mass instead of all of it. At CLOUD minimum size that is not what the
        # derivation gives: τ_orphan is 0.42 s at 10 percent subsaturation, which is FASTER than
        # the 2 s substep, so a cloud cell flickering into the corner for one step does lose
        # essentially all of its orphan mass. That is not a defect in the derivation - a 2 μm
        # droplet really does evaporate in a fraction of a second - it is the physical rate being
        # shorter than the timestep, and the design intent survives only in the weaker sense that
        # the removal is a rate rather than a set-to-zero. At RAIN minimum size the intent holds
        # with room to spare, because the same expression scales as D_min/x_min and the two
        # species are three orders of magnitude apart.
        τ_c = 1 / inv_τ
        @test FT(0.1) < τ_c < FT(2)
        inv_τ_r = CM2.orphan_mass_inv_timescale(aps, tps, T, ρ, pdf_r.xr_min, pdf_r.ρw, dry)
        @test 1 / inv_τ_r > FT(60)
        @test inv_τ / inv_τ_r > FT(100)
        @test FT(pdf_r.xr_min) > FT(pdf_c.xc_min)

        # (f) THROUGH THE ENTRY, where the branch is decided from the state, per species
        p3 = mp.ice.scheme
        logλ₀ = P3.get_distribution_logλ(P3.state_from_prognostic(p3, FT(0), FT(0), FT(0), FT(0)))
        qᵥ_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
        for (label, q_tot, drains) in
            (("subsaturated", FT(0.5) * qᵥ_sat + q, true),
            ("supersaturated", FT(1.2) * qᵥ_sat + q, false))
            pp = _per_process_2mp3(mp, tps, ρ, T, q_tot,
                q, FT(0), q, FT(0), FT(0), FT(0), FT(0), FT(0), logλ₀)
            if drains
                @test pp.cloud_orphan.q_lcl < 0
                @test pp.rain_orphan.q_rai < 0
            else
                @test pp.cloud_orphan.q_lcl == 0
                @test pp.rain_orphan.q_rai == 0
            end
            # the orphan slot is MASS ONLY: it removes what is there, it does not create a carrier
            @test all(
                iszero,
                (pp.cloud_orphan.n_lcl, pp.cloud_orphan.q_rai,
                    pp.cloud_orphan.n_rai, pp.cloud_orphan.q_ice, pp.cloud_orphan.n_ice),
            )
            @test all(
                iszero,
                (pp.rain_orphan.q_lcl, pp.rain_orphan.n_lcl,
                    pp.rain_orphan.n_rai, pp.rain_orphan.q_ice, pp.rain_orphan.n_ice),
            )
            # and the number adjustment is not manufacturing anything alongside it
            @test pp.cloud_numadj.n_lcl == 0
            @test pp.rain_numadj.n_rai == 0
        end

        # (g) A POPULATED CELL SEES NO ORPHAN SLOT AT ALL, which is what confines the whole
        # mechanism to the corner it was written for
        pp_pop = _per_process_2mp3(mp, tps, ρ, T, FT(0.5) * qᵥ_sat + q,
            q, FT(1e8), q, FT(1e3), FT(0), FT(0), FT(0), FT(0), logλ₀)
        @test pp_pop.cloud_orphan.q_lcl == 0 && pp_pop.rain_orphan.q_rai == 0

        # (h) f/J consistency of the drain: the diagonal is the exact -1/τ_orphan and it appears
        # exactly where the rate does
        x = BMT.MicroState2MP3{FT}(q, 0, q, 0, 0, 0, 0, 0)
        q_tot_dry = FT(0.5) * qᵥ_sat + q
        g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot_dry, logλ₀)
        (pp, rs) = _per_process_2mp3_and_riming(
            mp, tps, ρ, T, q_tot_dry, Tuple(x)..., logλ₀)
        J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
        se = BMT._liquid_sat_excess(tps, ρ, T, q_tot_dry, q, q, FT(0))
        inv_c = CM2.orphan_mass_inv_timescale(aps, tps, T, ρ, pdf_c.xc_min, pdf_c.ρw, se)
        inv_r = CM2.orphan_mass_inv_timescale(aps, tps, T, ρ, pdf_r.xr_min, pdf_r.ρw, se)
        @test J[1, 1] ≈ -inv_c rtol = sqrt(eps(FT))
        @test J[3, 3] ≈ -inv_r rtol = sqrt(eps(FT))
        @test all(isfinite, J)
        # the rate it linearizes, from the same numbers
        @test pp.cloud_orphan.q_lcl ≈ -q * inv_c rtol = sqrt(eps(FT))
        @test pp.rain_orphan.q_rai ≈ -q * inv_r rtol = sqrt(eps(FT))
    end
end

test_phantom_corner_invents_nothing(Float64)
test_phantom_corner_invents_nothing(Float32)



# The bare melt fraction `dLdt / ρq_ice` is unbounded: at a state observed in a crashed
# simulation - the smallest positive Float32 subnormal of ice mass under a populated number -
# the quotient overflows to Inf, and `-0.0 * Inf` then writes NaN into both rime slots. The
# fraction is formed once in `ice_melt`, bounded by the conduction-limited melt rate of a
# nucleation-size particle (`ice_melt_fraction_limit`), and the melt number rate is
# `ρn_ice * melt_frac`. What is asserted:
#
#   - the crash state: the fraction respects the bound and every slot of the melt
#     contribution, the per-process breakdown and the assembled Jacobian is finite;
#   - a trace state far above the subnormal edge is bounded the same way;
#   - healthy states: the bound does not bind and the fraction is the quotient, bit-exact;
#   - the number-rate identity `dNdt == ρn_ice * melt_frac`, bit-exact by construction.
function test_melt_fraction_is_conduction_bounded(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    vel = mp.ice.terminal_velocity
    aps = mp.warm_rain.air_properties
    quad = mp.ice.quad
    ρ = FT(1.0475167)

    @testset "the melt fraction is conduction-bounded [FT=$FT]" begin
        @testset "a subnormal-mass state stays finite" begin
            T = FT(291.95975)
            q_ice_sub = nextfloat(zero(FT))   # subnormal mass that passes a presence test
            n_ice_sub = FT(0.525)
            lim = P3.ice_melt_fraction_limit(aps, tps, p3, T)
            @test lim.inv_τ > 0

            # the fraction at the subnormal mass: the quotient overflows, the bound decides
            mf = P3.ice_melt_fraction(aps, tps, p3, T, q_ice_sub * ρ, FT(7.3e-7))
            @test !isfinite(FT(7.3e-7) / (q_ice_sub * ρ))
            @test mf.frac == lim.inv_τ
            @test isfinite(mf.∂frac_∂T)

            # the full melt contribution and both production paths at the state
            state = P3.state_from_prognostic(
                p3, q_ice_sub * ρ, n_ice_sub * ρ, FT(0), FT(0))
            logλ = P3.get_distribution_logλ(state)
            @test isfinite(logλ)
            melt = P3.ice_melt(vel, aps, tps, T, ρ, state, logλ; quad)
            @test all(isfinite, values(melt))
            @test 0 <= melt.melt_frac <= lim.inv_τ
            species = BMT._ice_melting_species(
                ρ, p3.ρ_i, FT(0), FT(0), melt.dNdt, melt.dLdt, melt.melt_frac)
            @test all(isfinite, Tuple(species))

            q_tot = FT(6e-3)
            x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, q_ice_sub, n_ice_sub, 0, 0)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
            f, J = BMT._tendency_and_jacobian(BMT.ManualJacobian(), g, x)
            @test all(isfinite, f)
            @test all(isfinite, J)
            @test all(isfinite, g(x))
        end

        @testset "a trace state is bounded the same way" begin
            T = FT(294.5)
            q_ice_trace = FT(5.5e-36)
            n_ice_trace = FT(1e3)
            lim = P3.ice_melt_fraction_limit(aps, tps, p3, T)
            mf = P3.ice_melt_fraction(aps, tps, p3, T, q_ice_trace * ρ, FT(1e-7) * ρ)
            @test mf.frac == lim.inv_τ
            state = P3.state_from_prognostic(
                p3, q_ice_trace * ρ, n_ice_trace * ρ, FT(0), FT(0))
            logλ = P3.get_distribution_logλ(state)
            melt = P3.ice_melt(vel, aps, tps, T, ρ, state, logλ; quad)
            @test all(isfinite, values(melt))
            @test melt.melt_frac <= lim.inv_τ
            @test melt.dNdt == state.ρn_ice * melt.melt_frac
        end

        @testset "healthy states keep the quotient bit-exactly" begin
            # the healthy reference state: q_ice 1.6e-5, dq 7.3e-7 per second, an
            # ordinary 0.046 per second fraction, two orders below the bound at 0.28 K
            T = TDI.T_freeze(tps) + FT(0.28)
            ρq = FT(1.6e-5) * ρ
            dLdt = FT(7.3e-7) * ρ
            lim = P3.ice_melt_fraction_limit(aps, tps, p3, T)
            @test dLdt / ρq < lim.inv_τ
            mf = P3.ice_melt_fraction(aps, tps, p3, T, ρq, dLdt)
            @test mf.frac == dLdt / ρq
            # and the derivative follows the quotient branch there
            ∂dLdt_∂T = FT(2e-6) * ρ
            mf∂ = P3.ice_melt_fraction(aps, tps, p3, T, ρq, dLdt, ∂dLdt_∂T)
            @test mf∂.∂frac_∂T == ∂dLdt_∂T / ρq
            # the limit's analytic temperature derivative against automatic differentiation
            ∂inv_τ_ad = FD.derivative(
                t -> P3.ice_melt_fraction_limit(aps, tps, p3, t).inv_τ, T)
            @test lim.∂inv_τ_∂T ≈ ∂inv_τ_ad rtol = sqrt(eps(FT))
        end

        @testset "the quotient branch's mass partial at trace mass is recorded" begin
            # The value is bounded on both branches; the quotient branch's derivative in
            # `ρq_ice` is `-dLdt / ρq_ice²`, genuinely enormous at trace mass and free to
            # overflow Float32. Whenever it is finite it satisfies `∂frac * ρq_ice == -frac`
            # up to rounding; asserted so the behavior a differentiating consumer sees at
            # trace mass is recorded rather than incidental.
            T = TDI.T_freeze(tps) + FT(0.28)
            ρq_tr = FT(1e-30)
            dLdt_tr = FT(1e-31)
            lim = P3.ice_melt_fraction_limit(aps, tps, p3, T)
            mf = P3.ice_melt_fraction(aps, tps, p3, T, ρq_tr, dLdt_tr)
            @test mf.frac == dLdt_tr / ρq_tr    # below the bound: the quotient branch
            @test mf.frac < lim.inv_τ
            ∂frac = FD.derivative(
                ρq -> P3.ice_melt_fraction(aps, tps, p3, T, ρq, dLdt_tr).frac, ρq_tr)
            @test !isfinite(∂frac) || isapprox(∂frac * ρq_tr, -mf.frac, rtol = sqrt(eps(FT)))
        end

        @testset "the number rate is the identity at a populated state" begin
            T = TDI.T_freeze(tps) + FT(2)
            state = P3.state_from_prognostic(
                p3, FT(1e-4) * ρ, FT(2e5) * ρ, FT(0.4) * FT(1e-4) * ρ,
                FT(0.4) * FT(1e-4) * ρ / p3.ρ_i)
            logλ = P3.get_distribution_logλ(state)
            melt = P3.ice_melt(vel, aps, tps, T, ρ, state, logλ; quad)
            @test melt.dNdt == state.ρn_ice * melt.melt_frac
            @test melt.∂dNdt_∂T == state.ρn_ice * melt.∂melt_frac_∂T
            @test melt.dLdt > 0
        end
    end
end

test_melt_fraction_is_conduction_bounded(Float64)
test_melt_fraction_is_conduction_bounded(Float32)

# Orphan ice: the warm-phase corner doctrine, one phase over. The ice number adjustment now
# passes `invent_from_zero = false`, so mass with no number is no longer given a manufactured
# carrier; the mass is removed as mass by `orphan_mass_drain_ice` below ice saturation, and
# the rime pair drains along the ray through the origin with the same fraction. What is
# asserted mirrors the phantom-corner testset: the drain fires subsaturated and not at or
# above ice saturation, it writes only the ice-mass and rime slots (total water is conserved
# the way every phase change here conserves it, through the fixed `q_tot`), no number is
# invented alongside it, and the manual Jacobian carries the exact `-1/τ_orphan` diagonals.
function test_ice_orphan_doctrine(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    aps = mp.warm_rain.air_properties
    numadj_ice = BMT._ice_numadj_params(p3)
    ρ = FT(0.9)
    T = FT(260)
    q = FT(5e-8)
    q_rim = FT(0.4) * q
    b_rim = q_rim / FT(800)
    x_min_ice = P3.ice_mean_particle_mass_min(p3)
    qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
    logλ₀ = P3.get_distribution_logλ(
        P3.state_from_prognostic(p3, q * ρ, FT(0), q_rim * ρ, b_rim * ρ))

    ppcall(q_tot) = _per_process_2mp3(mp, tps, ρ, T, q_tot,
        FT(0), FT(0), FT(0), FT(0), q, FT(0), q_rim, b_rim, logλ₀)

    @testset "orphan ice is drained as mass, not carried [FT=$FT]" begin
        @testset "no number is invented for mass without number" begin
            ∂ₜn = CM2.number_tendency_from_mass_limits(
                numadj_ice, q, FT(0); invent_from_zero = false)
            @test ∂ₜn == 0
            # the invention the flag removes, computed inline so this cannot pass vacuously
            @test CM2.number_tendency_from_mass_limits(numadj_ice, q, FT(0)) > 0
        end

        @testset "the drain fires subsaturated only, on the ice-mass and rime slots" begin
            for (q_tot, drains) in
                ((FT(0.5) * qᵥ_sat_ice + q, true), (FT(1.2) * qᵥ_sat_ice + q, false))
                pp = ppcall(q_tot)
                if drains
                    @test pp.ice_orphan.q_ice < 0
                    @test pp.ice_orphan.q_rim < 0
                    @test pp.ice_orphan.b_rim < 0
                    # the ray: both rime moments and the mass lose one fraction
                    @test pp.ice_orphan.q_rim / pp.ice_orphan.q_ice ≈ q_rim / q rtol =
                        sqrt(eps(FT))
                    @test pp.ice_orphan.b_rim / pp.ice_orphan.q_ice ≈ b_rim / q rtol =
                        sqrt(eps(FT))
                else
                    @test iszero(pp.ice_orphan.q_ice)
                    @test iszero(pp.ice_orphan.q_rim)
                    @test iszero(pp.ice_orphan.b_rim)
                end
                # mass and rime only: no vapor, liquid, rain or number slot moves
                @test all(
                    iszero,
                    (pp.ice_orphan.q_lcl, pp.ice_orphan.n_lcl,
                        pp.ice_orphan.q_rai, pp.ice_orphan.n_rai, pp.ice_orphan.n_ice),
                )
                # and the adjustment is not manufacturing a population alongside it
                @test pp.ice_numadj.n_ice == 0
            end
        end

        @testset "the drain is the derived relaxation, linear in the mass" begin
            q_tot_dry = FT(0.5) * qᵥ_sat_ice + q
            se = BMT._ice_sat_excess(tps, ρ, T, q_tot_dry, FT(0), FT(0), q)
            inv_τ = CM2.orphan_mass_inv_timescale_ice(
                aps, tps, T, ρ, x_min_ice, p3.ρ_i, se)
            @test inv_τ > 0 && isfinite(inv_τ)
            pp = ppcall(q_tot_dry)
            @test pp.ice_orphan.q_ice ≈ -q * inv_τ rtol = sqrt(eps(FT))
            drain = CM2.orphan_mass_drain_ice(
                aps, tps, T, ρ, q, x_min_ice, p3.ρ_i, se)
            @test drain ≈ -q * inv_τ rtol = sqrt(eps(FT))
            @test CM2.orphan_mass_drain_ice(
                aps, tps, T, ρ, 3 * q, x_min_ice, p3.ρ_i, se) ≈ 3 * drain rtol =
                sqrt(eps(FT))
            # zero at and above ice saturation
            @test iszero(CM2.orphan_mass_drain_ice(
                aps, tps, T, ρ, q, x_min_ice, p3.ρ_i, FT(1e-4)))
            @test iszero(CM2.orphan_mass_drain_ice(
                aps, tps, T, ρ, q, x_min_ice, p3.ρ_i, zero(FT)))
        end

        @testset "a populated cell sees no orphan slot" begin
            q_tot_dry = FT(0.5) * qᵥ_sat_ice + q
            n_ice = q / (2 * x_min_ice)
            pp = _per_process_2mp3(mp, tps, ρ, T, q_tot_dry,
                FT(0), FT(0), FT(0), FT(0), q, n_ice, q_rim, b_rim, logλ₀)
            @test all(iszero, Tuple(pp.ice_orphan))
        end

        @testset "f/J consistency: the exact diagonals on q_ice and the rime pair" begin
            q_tot_dry = FT(0.5) * qᵥ_sat_ice + q
            x = BMT.MicroState2MP3{FT}(0, 0, 0, 0, q, 0, q_rim, b_rim)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot_dry, logλ₀)
            (pp, rs) = _per_process_2mp3_and_riming(
                mp, tps, ρ, T, q_tot_dry, Tuple(x)..., logλ₀)
            J = BMT._jacobian_2mp3_manual(g, x, pp, rs)
            se = BMT._ice_sat_excess(tps, ρ, T, q_tot_dry, FT(0), FT(0), q)
            inv_τ = CM2.orphan_mass_inv_timescale_ice(
                aps, tps, T, ρ, x_min_ice, p3.ρ_i, se)
            @test J[5, 5] ≈ -inv_τ rtol = sqrt(eps(FT))
            @test J[7, 7] ≈ -inv_τ rtol = sqrt(eps(FT))
            @test J[8, 8] ≈ -inv_τ rtol = sqrt(eps(FT))
            @test all(isfinite, J)
        end
    end
end

test_ice_orphan_doctrine(Float64)
test_ice_orphan_doctrine(Float32)

function test_logl_refresh(FT)
    @testset "the shape parameter is refreshed from the marched state [FT=$FT]" begin
        mp = CMP.Microphysics2MParams(FT; with_ice = true)
        ρ = FT(1)
        x = BMT.MicroState2MP3{FT}(
            FT(1e-4), FT(1e8), FT(1e-5), FT(1e4),
            FT(1e-3), FT(1e5), FT(5e-4), FT(1e-6))

        # (1) CONSISTENT AT THE ENTRY STATE. The first substep re-uses the caller's value rather
        # than re-solving, and that is only sound if the two agree: this asserts the shortcut does
        # not quietly start the march from a different distribution than a refresh would give.
        st = P3.state_from_prognostic(
            mp.ice.scheme, ρ * x.q_ice, ρ * x.n_ice, ρ * x.q_rim, ρ * x.b_rim)
        @test BMT._refreshed_logλ(mp, ρ, x) === P3.get_distribution_logλ(st)

        # (2) LIVE WHEN THE STATE MOVES. A refresh that returned the same value whatever the march
        # did would pass (1) and be worthless; this is the knob-fires half.
        moved = BMT.MicroState2MP3{FT}(
            x.q_lcl, x.n_lcl, x.q_rai, x.n_rai,
            FT(1e-3), FT(1e7), x.q_rim, x.b_rim)   # a hundredfold more, far smaller, crystals
        @test BMT._refreshed_logλ(mp, ρ, moved) != BMT._refreshed_logλ(mp, ρ, x)

        # (3) SAFE AT THE DEGENERATE STATES THE MARCH CAN REACH. This is what makes a per-substep
        # re-solve admissible at all: mid-march the ice slots are not guaranteed to describe a
        # population, and the solve must return a bracket end rather than an arbitrary value.
        for degenerate in (
            BMT.MicroState2MP3{FT}(x.q_lcl, x.n_lcl, x.q_rai, x.n_rai,
                zero(FT), FT(1e5), zero(FT), zero(FT)),   # number, no mass
            BMT.MicroState2MP3{FT}(x.q_lcl, x.n_lcl, x.q_rai, x.n_rai,
                zero(FT), zero(FT), zero(FT), zero(FT)),  # neither
            BMT.MicroState2MP3{FT}(x.q_lcl, x.n_lcl, x.q_rai, x.n_rai,
                FT(1e-3), FT(1e5), FT(-1e-9), FT(-1e-12)), # inadmissible pair
        )
            @test isfinite(BMT._refreshed_logλ(mp, ρ, degenerate))
        end
    end
end
test_logl_refresh(Float64)
test_logl_refresh(Float32)
