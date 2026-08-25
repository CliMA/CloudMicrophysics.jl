using Test

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI

include("p3_state_battery.jl")

# A3: the threshold and boundary battery. Consolidates:
#   1. the archived crash corpus (ported from the campaign's `diag/corpus_stability_test.jl`
#      and its design note `notes/corpus-test-design.md`, both in
#      `~/science/projects/p3-2mp3-stability`), evaluated on three of that harness's four
#      public surfaces (the fourth, gated melt-fraction census, is soft/informational there
#      and is not ported);
#   2. the non-positive air-density floor, the value-lane AD branch guards at differentiated
#      zero corners, the ice-deposition degeneracy gate, the ice-population presence
#      predicate, and the substep water bound - all ported from the campaign's
#      `test/rosenbrock_mode_tests.jl` (`campaign/he/p3-density-floors` in this project's CM
#      clone) and adapted to the rebuilt tip's renamed entry points
#      (`_per_process_2mp3`/`_per_process_2mp3_and_riming` -> `p3_2m_process_rates`;
#      `P3.ice_nucleation_mass` -> `CMP.ice_seed(...).m_nuc`, verified against
#      `src/P3_processes.jl`'s current `ice_population_is_present` docstring, not assumed
#      from the campaign source);
#   3. the state-battery's own corner states (`test/p3_state_battery.jl`, already at this
#      tip), run through the production substep for a broad finiteness sweep.

#####
##### 1. The archived crash corpus
#####
# Ported from `diag/corpus_stability_test.jl::build_analytic_corpus()` (campaign project,
# not this repo). Density form (ρq_*, ρn_* in kg/m3, 1/m3), matching how every reference
# state in that project's PROJECT_STATE.md and every P3State field is recorded. `logλ =
# NaN` means "compute it from the prognostic densities"; a real value means "use exactly
# this", for the one state where the trajectory's own logλ was recorded.
const CrashState = NamedTuple{
    (:name, :ρ, :T, :ρq_tot, :ρq_lcl, :ρn_lcl, :ρq_rai, :ρn_rai, :ρq_ice, :ρn_ice, :ρq_rim, :ρb_rim, :logλ),
    Tuple{
        String,
        Float64,
        Float64,
        Float64,
        Float64,
        Float64,
        Float64,
        Float64,
        Float64,
        Float64,
        Float64,
        Float64,
        Float64,
    },
}

function build_crash_corpus()
    out = CrashState[]

    # Item 1: THE KILLING STATE (round-5 stage-state dump, m4ctrl, t=61047). Exact F32
    # values from the campaign's PROJECT_STATE.md. POSITIVE CONTROL: on the pre-floor code
    # this minted NaN in (dq_rim_dt, db_rim_dt) at Float32 via a quotient overflowing
    # (ρq_ice = 1e-45 is subnormal at F32; the same state is finite-but-absurd at F64).
    push!(
        out,
        (
            name = "killing_state", ρ = 1.0475167, T = 291.95975,
            ρq_tot = 0.0140373735,
            ρq_lcl = -1.4502772e-19, ρn_lcl = -2.1273305e-5,
            ρq_rai = -5.202259e-15, ρn_rai = -0.008766984,
            ρq_ice = 1.0e-45,                 # subnormal at F32
            ρn_ice = 0.55014247,
            ρq_rim = 0.0, ρb_rim = 0.0,       # the crash's rime pair, exactly zero
            logλ = 6.297534,                  # recorded on the real trajectory, used as-is
        ),
    )

    # Item 2: THE TRACE STATE (surviving run, remnant cell). Archive gives q_ice =
    # 5.5e-36 kg/kg (specific) and T = 294.5 K; n_ice is not recorded exactly. DOCUMENTED
    # CHOICE (matching the campaign harness): ρ = 1.2 kg/m3, n_ice = 1.0 #/kg specific
    # (same order as the killing state's 0.525 #/kg).
    let ρ = 1.2
        push!(
            out,
            (
                name = "trace_state", ρ = ρ, T = 294.5,
                ρq_tot = 1.2e-3,
                ρq_lcl = 0.0, ρn_lcl = 0.0, ρq_rai = 0.0, ρn_rai = 0.0,
                ρq_ice = 5.5e-36 * ρ, ρn_ice = 1.0 * ρ,
                ρq_rim = 0.0, ρb_rim = 0.0, logλ = NaN,
            ),
        )
    end

    # Item 3: THE HEALTHY REFERENCE. Archive gives q_ice = 1.6e-5 kg/kg at T = T_freeze +
    # 0.28 K; n_ice not recorded (DOCUMENTED CHOICE: 1.6e4 #/kg specific, an ordinary
    # non-degenerate population). `T = NaN` here is this file's own sentinel (not the
    # campaign harness's), filled at evaluation time from `T_freeze(tps) + 0.28`.
    let ρ = 1.2
        push!(
            out,
            (
                name = "healthy_reference", ρ = ρ, T = NaN,
                ρq_tot = 1.2e-3,
                ρq_lcl = 0.0, ρn_lcl = 0.0, ρq_rai = 0.0, ρn_rai = 0.0,
                ρq_ice = 1.6e-5 * ρ, ρn_ice = 1.6e4 * ρ,
                ρq_rim = 0.0, ρb_rim = 0.0, logλ = NaN,
            ),
        )
    end

    # Item 4: FINGERPRINT REFERENCE STATES (the campaign's
    # `notes/stage2-fingerprint-reference.txt`), five levels, three with NEGATIVE n_ice
    # mid-stage, at t=60800 in a stage-2 trace. The archive records only (n_ice, q_ice) at
    # each level (both specific); DOCUMENTED CHOICE (matching the campaign harness): ρ =
    # 1.0 exactly, so density and specific coincide numerically; evaluated at both a
    # sub-freezing (260 K, melt inert) and a super-freezing (280 K, melt active) companion
    # temperature.
    fingerprint = [
        ("lev7", -4.9747825e-21, 1.126e-41),
        ("lev8", -6.004194e-20, 3.3e-42),
        ("lev9", 1.8072207e-18, 0.0),
        ("lev10", 1.0553206e-16, -1.8e-44),
        ("lev11", -2.594197e-17, -3.57e-43),
    ]
    for (lev, n_ice, q_ice) in fingerprint, Tc in (260.0, 280.0)
        push!(
            out,
            (
                name = "fingerprint_$(lev)_T$(Int(Tc))", ρ = 1.0, T = Tc,
                ρq_tot = 0.0,
                ρq_lcl = 0.0, ρn_lcl = 0.0, ρq_rai = 0.0, ρn_rai = 0.0,
                ρq_ice = q_ice, ρn_ice = n_ice,
                ρq_rim = 0.0, ρb_rim = 0.0, logλ = NaN,
            ),
        )
    end

    # Item 6: MASSLESS / DEGENERATE EXEMPLARS, composed from the classes seen in items
    # 1-5. Standard companion: ρ = 1.2, T = 280 K (above freezing, melt path active unless
    # the exemplar itself is what is being tested).
    ρc, Tc = 1.2, 280.0

    # 6a: n_ice > 0, q_ice exactly 0 (number without mass).
    push!(
        out,
        (
            name = "exemplar_ice_number_no_mass", ρ = ρc, T = Tc, ρq_tot = 0.0,
            ρq_lcl = 0.0, ρn_lcl = 0.0, ρq_rai = 0.0, ρn_rai = 0.0,
            ρq_ice = 0.0, ρn_ice = 1.0 * ρc,
            ρq_rim = 0.0, ρb_rim = 0.0, logλ = NaN,
        ),
    )

    # 6b: each of the eight density moments negative ALONE, magnitude 1e-6 * ρc, everything
    # else exactly zero - the round-5 dump's negative-sibling-argument class, one species
    # at a time so a failure names its own carrier.
    species_slots = (:ρq_lcl, :ρn_lcl, :ρq_rai, :ρn_rai, :ρq_ice, :ρn_ice, :ρq_rim, :ρb_rim)
    for sp in species_slots
        vals = Dict(k => 0.0 for k in species_slots)
        vals[sp] = -1.0e-6 * ρc
        push!(
            out,
            (
                name = "exemplar_negative_$(sp)", ρ = ρc, T = Tc, ρq_tot = 0.0,
                ρq_lcl = vals[:ρq_lcl], ρn_lcl = vals[:ρn_lcl],
                ρq_rai = vals[:ρq_rai], ρn_rai = vals[:ρn_rai],
                ρq_ice = vals[:ρq_ice], ρn_ice = vals[:ρn_ice],
                ρq_rim = vals[:ρq_rim], ρb_rim = vals[:ρb_rim], logλ = NaN,
            ),
        )
    end

    # 6c: THE NECESSARY-TRIO PATTERN generalized beyond ice - subnormal-at-F32 mass with
    # populated number and (for ice) an exactly-zero rime pair, one carrier species at a
    # time. The killing state is the ice instance of this; lcl/rai have no rime pair.
    fmin32 = 1.0e-45  # subnormal at Float32 (min subnormal Float32 is ~1.4e-45)
    push!(
        out,
        (
            name = "exemplar_subnormal_ice_mass_trio", ρ = ρc, T = Tc, ρq_tot = 0.0,
            ρq_lcl = 0.0, ρn_lcl = 0.0, ρq_rai = 0.0, ρn_rai = 0.0,
            ρq_ice = fmin32, ρn_ice = 1.0 * ρc, ρq_rim = 0.0, ρb_rim = 0.0, logλ = NaN,
        ),
    )
    push!(
        out,
        (
            name = "exemplar_subnormal_lcl_mass", ρ = ρc, T = Tc, ρq_tot = 0.0,
            ρq_lcl = fmin32, ρn_lcl = 1.0 * ρc, ρq_rai = 0.0, ρn_rai = 0.0,
            ρq_ice = 0.0, ρn_ice = 0.0, ρq_rim = 0.0, ρb_rim = 0.0, logλ = NaN,
        ),
    )
    push!(
        out,
        (
            name = "exemplar_subnormal_rai_mass", ρ = ρc, T = Tc, ρq_tot = 0.0,
            ρq_lcl = 0.0, ρn_lcl = 0.0, ρq_rai = fmin32, ρn_rai = 1.0 * ρc,
            ρq_ice = 0.0, ρn_ice = 0.0, ρq_rim = 0.0, ρb_rim = 0.0, logλ = NaN,
        ),
    )

    # 6d: mass without number (the ice "birth corner" analog).
    push!(
        out,
        (
            name = "exemplar_ice_mass_no_number", ρ = ρc, T = Tc, ρq_tot = 0.0,
            ρq_lcl = 0.0, ρn_lcl = 0.0, ρq_rai = 0.0, ρn_rai = 0.0,
            ρq_ice = 1.0e-5 * ρc, ρn_ice = 0.0, ρq_rim = 0.0, ρb_rim = 0.0, logλ = NaN,
        ),
    )

    # 6e: rime present with ice fully absent (an orphan-rime exemplar; the ice-process
    # block should be entirely inert here since `ice_population_is_present` requires both
    # ice moments).
    push!(
        out,
        (
            name = "exemplar_rime_no_ice", ρ = ρc, T = Tc, ρq_tot = 0.0,
            ρq_lcl = 0.0, ρn_lcl = 0.0, ρq_rai = 0.0, ρn_rai = 0.0,
            ρq_ice = 0.0, ρn_ice = 0.0,
            ρq_rim = 1.0e-5 * ρc, ρb_rim = 1.0e-8 * ρc, logλ = NaN,
        ),
    )

    return out
end

const CRASH_CORPUS = build_crash_corpus()

for FT in (Float32, Float64)
    mp = CMP.Microphysics2MParams(
        FT; with_ice = true, is_limited = true, aerosol = CMP.PrescribedAerosol(FT),
    )
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    T_freeze = TDI.T_freeze(tps)

    @testset "A3 crash corpus ($FT)" begin
        for cs in CRASH_CORPUS
            @testset "$(cs.name)" begin
                ρ = FT(cs.ρ)
                T = isnan(cs.T) ? T_freeze + FT(0.28) : FT(cs.T)
                q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim =
                    FT(cs.ρq_tot) / ρ, FT(cs.ρq_lcl) / ρ, FT(cs.ρn_lcl) / ρ,
                    FT(cs.ρq_rai) / ρ, FT(cs.ρn_rai) / ρ, FT(cs.ρq_ice) / ρ,
                    FT(cs.ρn_ice) / ρ, FT(cs.ρq_rim) / ρ, FT(cs.ρb_rim) / ρ

                state = P3.state_from_prognostic(
                    mp.ice.scheme, FT(cs.ρq_ice), FT(cs.ρn_ice), FT(cs.ρq_rim), FT(cs.ρb_rim),
                )
                logλ = isnan(cs.logλ) ? P3.get_distribution_logλ(state) : FT(cs.logλ)
                @test isfinite(logλ)

                # (A) full raw tendency
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
                x = BMT.MicroState2MP3{FT}(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
                tend = g(x)
                @test all(isfinite, tend)
                # empty-category sink check: a category with no mass and no number
                # (post-clamp) must not receive a negative rate out of nothing
                for (qv, nv, dqdt, sname) in (
                    (q_lcl, n_lcl, tend.q_lcl, "lcl"), (q_rai, n_rai, tend.q_rai, "rai"),
                    (q_ice, n_ice, tend.q_ice, "ice"),
                )
                    empty = max(qv, zero(FT)) == 0 && max(nv, zero(FT)) == 0
                    if empty && isfinite(dqdt)
                        @test dqdt >= 0
                    end
                end

                # (C) production substep, Δt = 2 s, nsub = 1 - matching the campaign box
                # timestep and the harness this corpus is ported from
                subres = BMT.bulk_microphysics_tendencies(
                    BMT.rosenbrock_manual(), BMT.Microphysics2Moment(), mp, tps,
                    ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim,
                    logλ, FT(2), 1, zero(FT), zero(FT),
                )
                # `values(subres)` ends with `extras`, a `NamedTuple` (empty in production,
                # the T1 seam per PROJECT_STATE.md), so `isfinite` cannot run over it
                # directly; `Base.front` drops that trailing field and keeps the eight
                # species tendencies plus the activation diagnostic, all real numbers
                @test all(isfinite, Base.front(values(subres)))

                # (B) isolated, ungated melt rate - bypasses whatever presence gate the
                # entry above applies; role is finiteness, sign, and exact zero when ice
                # is fully absent
                melt = P3.ice_melt(
                    mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps, T, ρ,
                    state, logλ; quad = mp.ice.quad,
                )
                @test isfinite(melt.dNdt) && isfinite(melt.dLdt)
                @test isfinite(melt.∂dNdt_∂T) && isfinite(melt.∂dLdt_∂T)
                isfinite(melt.dNdt) && @test melt.dNdt >= 0  # melting cannot un-melt
                isfinite(melt.dLdt) && @test melt.dLdt >= 0
                ice_absent = cs.ρq_ice <= 0 && cs.ρn_ice <= 0
                if ice_absent
                    @test melt.dNdt == 0 && melt.dLdt == 0
                end
            end
        end
    end
end

#####
##### 2. Boundary and gate tests ported from the campaign's rosenbrock_mode_tests.jl
#####

for FT in (Float32, Float64)
    mp = CMP.Microphysics2MParams(
        FT; with_ice = true, is_limited = true, aerosol = CMP.PrescribedAerosol(FT),
    )
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    @testset "the all-zero state returns exactly zero ($FT)" begin
        # near-empty species mask -> explicit substeps -> exactly zero; ported from the
        # campaign's "degenerate and trivial states" testset
        out = BMT.bulk_microphysics_tendencies(
            BMT.rosenbrock_manual(), BMT.Microphysics2Moment(), mp, tps,
            FT(1), FT(273), FT(0),
            FT(0), FT(0), FT(0), FT(0), FT(0), FT(0), FT(0), FT(0),
            FT(-Inf), FT(60), 4,
        )
        @test all(iszero, Base.front(values(out)))  # `extras` excluded, see the crash-corpus comment above
    end

    @testset "the air-density floor is finite and inert ($FT)" begin
        # the manual Jacobian dispatches on MicroState2MP3, not a bare SVector
        x = BMT.MicroState2MP3{FT}(FT(4e-4), FT(5e7), FT(1e-4), FT(4e4), FT(1e-4), FT(2e5), FT(4e-5), FT(6e-8))
        logλ_at(ρ) = P3.get_distribution_logλ(
            P3.state_from_prognostic(mp.ice.scheme, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]),
        )
        for ρ_bad in FT[-1, -0.245, -floatmin(FT), 0]
            g = BMT.Instantaneous2MP3Tendency(
                mp, tps, ρ_bad, FT(273.5), FT(0.009), logλ_at(max(ρ_bad, BMT.AIR_DENSITY_FLOOR)),
            )
            f, J = BMT._tendency_and_jacobian(BMT.ManualJacobian(), g, x)
            @test all(isfinite, f)
            @test all(isfinite, J)
            @test all(isfinite, g(x))
        end
        # inert at physical densities: the realized box density minimum is of order
        # 1e-2, well above the floor
        for ρ_ok in FT[0.011, 0.05, 0.5, 1.0]
            g1 = BMT.Instantaneous2MP3Tendency(mp, tps, ρ_ok, FT(273.5), FT(0.009), logλ_at(ρ_ok))
            f1, J1 = BMT._tendency_and_jacobian(BMT.ManualJacobian(), g1, x)
            g2 = BMT.Instantaneous2MP3Tendency(
                mp, tps, max(ρ_ok, BMT.AIR_DENSITY_FLOOR), FT(273.5), FT(0.009), logλ_at(ρ_ok),
            )
            f2, J2 = BMT._tendency_and_jacobian(BMT.ManualJacobian(), g2, x)
            @test J1 == J2
            @test f1 == f2
        end
    end

    @testset "value-lane branch guards at differentiated-zero corners ($FT)" begin
        # each of these puts a different guarded quantity at a differentiated zero: a
        # `ForwardDiff.Dual` whose VALUE is zero but whose partials are not, which
        # `iszero(::Dual)` (requiring zero partials) fails, selecting a different branch
        # from the plain-float primal
        ρ, T, q_tot = FT(0.78), FT(268), FT(0.009)
        corners = (
            # rime mass with zero rime volume -> ρ_rim = 0, the B_rim guard
            ("zero rime volume", BMT.MicroState2MP3{FT}(2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 0)),
            # ice mass with zero ice number -> the collision-ratio guards
            ("ice mass, no number", BMT.MicroState2MP3{FT}(2e-4, 5e7, 1e-4, 4e4, 1e-4, 0, 0, 0)),
            # rain mass with zero rain number -> the rain-activity guards
            ("rain mass, no number", BMT.MicroState2MP3{FT}(2e-4, 5e7, 1e-4, 0, 1e-4, 2e5, 4e-5, 6e-8)),
            # unrimed ice -> F_rim = 0, the isunrimed threshold selector
            ("unrimed ice", BMT.MicroState2MP3{FT}(2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 0, 0)),
        )
        for (name, x) in corners
            @testset "$name" begin
                logλ = P3.get_distribution_logλ(
                    P3.state_from_prognostic(mp.ice.scheme, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]),
                )
                @test isfinite(logλ)
                g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)
                f_primal = g(x)
                f_dual, _J = BMT._tendency_and_jacobian(BMT.ExactJacobian(), g, x)
                @test all(isfinite, f_primal)
                # same BRANCH, not the same rounding: a missed guard is orders of
                # magnitude out, not a last-digit difference (the dual arithmetic
                # reassociates, so bit-identity is not asserted)
                @test all(isapprox.(Tuple(f_dual), Tuple(f_primal); rtol = FT(1.0e-3), atol = FT(0)))
                out = BMT.bulk_microphysics_tendencies(
                    BMT.rosenbrock_exact(), BMT.Microphysics2Moment(), mp, tps, ρ, T, q_tot,
                    x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], logλ, FT(2), 1,
                )
                # `extras` again excluded, see the crash-corpus surface (C) comment above
                @test all(isfinite, Base.front(values(out)))
            end
        end
    end

    @testset "the ice-deposition degeneracy gate ($FT)" begin
        ρ, q_tot = FT(0.6), FT(4.0e-3)
        empty_state = P3.state_from_prognostic(mp.ice.scheme, FT(0), FT(0), FT(0), FT(0))
        logλ0 = P3.get_distribution_logλ(empty_state)
        for T in FT[220, 233, 245, 253, 258, 261, 263, 265, 268, 273]
            τ = P3.ice_deposition_timescale(
                mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps,
                T, ρ, empty_state, logλ0; quad = mp.ice.quad,
            )
            @test P3.ice_deposition_is_degenerate(τ)
            micro = (;
                q_tot, q_lcl = FT(0), n_lcl = FT(0), q_rai = FT(0), n_rai = FT(0),
                q_ice = FT(0), n_ice = FT(0), q_rim = FT(0), b_rim = FT(0),
            )
            thermo = (; ρ, T, w = FT(0), p = FT(0), logλ = logλ0)
            pp, _rs = BMT.p3_2m_process_rates(mp, tps, micro, thermo)
            @test all(iszero, Tuple(pp.ice_depsub))
        end

        # a populated state must NOT be degenerate, or the gate would silence real
        # deposition - the inertness half of the claim
        q_ice, n_ice = FT(1.0e-5), FT(1.0e4)
        st = P3.state_from_prognostic(mp.ice.scheme, ρ * q_ice, ρ * n_ice, FT(0), FT(0))
        logλ = P3.get_distribution_logλ(st)
        for T in FT[233, 253, 263]
            τ = P3.ice_deposition_timescale(
                mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps, T, ρ, st, logλ;
                quad = mp.ice.quad,
            )
            @test isfinite(τ) && τ > 0
            @test !P3.ice_deposition_is_degenerate(τ)
            micro = (;
                q_tot, q_lcl = FT(0), n_lcl = FT(0), q_rai = FT(0), n_rai = FT(0),
                q_ice, n_ice, q_rim = FT(0), b_rim = FT(0),
            )
            thermo = (; ρ, T, w = FT(0), p = FT(0), logλ)
            pp, _rs = BMT.p3_2m_process_rates(mp, tps, micro, thermo)
            # supersaturated with respect to ice at all three temperatures: growth is positive
            @test pp.ice_depsub.q_ice > 0 && isfinite(pp.ice_depsub.q_ice)
        end
    end

    @testset "the ice-population presence predicate ($FT)" begin
        # `ice_population_is_present`'s own current formula (`src/P3_processes.jl`):
        # `(ρq_ice / m_nuc > ρn_ice) & (ρn_ice > 0)`, verified against source rather than
        # assumed from the campaign's (now renamed) `P3.ice_nucleation_mass`
        mkstate(q_ice, n_ice) =
            P3.state_from_prognostic(mp.ice.scheme, q_ice * FT(1), n_ice * FT(1), FT(0), FT(0))
        x̄ = FT(1.0e-11)  # mean ice particle mass held fixed across the sweep [kg]
        @test P3.ice_population_is_present(mkstate(FT(1.0e-4), FT(1.0e-4) / x̄))
        @test P3.ice_population_is_present(mkstate(FT(1.0e-9), FT(1.0e-9) / x̄))
        @test !P3.ice_population_is_present(mkstate(FT(0), FT(0)))
        @test !P3.ice_population_is_present(mkstate(FT(5.0e-8), FT(0)))  # mass, no number
        @test !P3.ice_population_is_present(mkstate(FT(0), FT(1.0e4)))   # number, no mass

        m_nuc = CMP.ice_seed(mp.ice.scheme).m_nuc
        N = FT(1.0e4)
        @test !P3.ice_population_is_present(mkstate(FT(0.5) * N * m_nuc, N))  # < 1 crystal/particle
        @test P3.ice_population_is_present(mkstate(FT(2) * N * m_nuc, N))     # > 1 crystal/particle
    end
end

@testset "substep water bound" begin
    for FT in (Float32, Float64)
        @testset "$FT" begin
            @testset "inert on an admissible increment" begin
                x = BMT.MicroState2MP3{FT}(4e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8)
                d = BMT.MicroState2MP3{FT}(1e-6, 1e4, 1e-6, 1e2, 1e-6, 1e2, 1e-7, 1e-10)
                @test BMT._water_bounded_increment(x, d, FT(0.02)) == d
            end
            @testset "clips an increment that exceeds total water" begin
                x = BMT.MicroState2MP3{FT}(1e-4, 1e7, 1e-4, 1e4, 1e-4, 1e5, 0, 0)
                d = BMT.MicroState2MP3{FT}(0, 0, 0, 0, FT(0.9), 0, 0, 0)
                q_tot = FT(0.016)
                db = BMT._water_bounded_increment(x, d, q_tot)
                @test BMT._condensate_total(x .+ db) <= q_tot * (1 + 10 * eps(FT))
            end
            @testset "q_rim is not double-counted in the condensate total" begin
                x = BMT.MicroState2MP3{FT}(1e-4, 1e7, 2e-4, 1e4, 3e-4, 1e5, 2e-4, 5e-7)
                @test BMT._condensate_total(x) == x.q_lcl + x.q_rai + x.q_ice
            end
            @testset "no budget supplied leaves the bare step unchanged" begin
                x = BMT.MicroState2MP3{FT}(1e-4, 1e7, 1e-4, 1e4, 1e-4, 1e5, 0, 0)
                d = BMT.MicroState2MP3{FT}(0, 0, 0, 0, FT(0.9), 0, 0, 0)
                @test BMT._bounded_explicit_step(x, d, nothing) === d
            end
        end
    end
end

#####
##### 3. The state-battery's corner states through the production substep
#####

@testset "corner states stay finite through the production substep (Float64)" begin
    FT = Float64
    mp = CMP.Microphysics2MParams(
        FT; with_ice = true, is_limited = true, aerosol = CMP.PrescribedAerosol(FT),
    )
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    states = generate_comprehensive_states(FT; n_lhs = 0, seed = 1)  # corner_states only
    @test length(states) > 0
    for s in states
        state = P3.state_from_prognostic(mp.ice.scheme, s.ρ * s.q_ice, s.ρ * s.n_ice, s.ρ * s.q_rim, s.ρ * s.b_rim)
        logλ = P3.get_distribution_logλ(state)
        out = BMT.bulk_microphysics_tendencies(
            BMT.rosenbrock_manual(), BMT.Microphysics2Moment(), mp, tps,
            s.ρ, s.T, s.q_tot, s.q_lcl, s.n_lcl, s.q_rai, s.n_rai, s.q_ice, s.n_ice,
            s.q_rim, s.b_rim, logλ, FT(2), 1,
        )
        @test all(isfinite, Base.front(values(out)))  # `extras` excluded, see above
    end
end
