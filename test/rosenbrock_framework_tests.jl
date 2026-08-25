using Test

import JET
import BenchmarkTools as BT

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.HetIceNucleation as CM_HetIce
import CloudMicrophysics.Utilities as UT
import ForwardDiff as FD

# The unified `RosenbrockAverage{Jacobian, GrowthTreatment, TendencyLimiter}` framework
# on the 2M+P3 model: the presets (`rosenbrock_coupled`, `rosenbrock_exact`,
# `rosenbrock_manual`), the `ExactJacobian`/`ManualJacobian` contract, and the
# `_jacobian_2mp3_manual` closed-form entries validated against `ForwardDiff`.
#
# PORTED FROM the campaign's `test/rosenbrock_framework_tests.jl`
# (`campaign/he/p3-density-floors` in this project's CM clone). This file's own 1M
# section (`test_framework_1m`, the `rosenbrock_donor() ≡ LinearizedAverage()`
# equivalence and the framework-based 1M presets) is DROPPED, not repaired: the
# `RosenbrockAverage`-on-`Microphysics1Moment` framework it exercises
# (`src/BMT_1m.jl`, `_rosenbrock_average_1m`, `MicroState1M`) does not exist on this
# branch at all - it lives on `he/p3-cm10b-1m`, a sibling branch off CM-10a that has
# not merged onto this trunk (`he/p3-cmLUT-tables`). That branch carries its own
# `test/rosenbrock_1m_tests.jl` covering the same ground (checked, not redone here,
# per instruction - see the report). `rosenbrock_donor()` itself no longer exists as a
# public name anywhere on this trunk (D4 ruling: deleted, its one production use
# inlined; `RosenbrockAverage()` with bare defaults - `DonorJacobian`,
# `ImplicitGrowth`, `NoLimiter` - is the equivalent construction where a throw-check
# still needs a donor-shaped mode, used once below).
#
# STALE REFERENCES REPAIRED, beyond the 1M section and the `rosenbrock_donor()` name:
#   - The per-process evaluator is `p3_2m_process_rates` (not `_per_process_2mp3` /
#     `_per_process_2mp3_and_riming`); this file does not call it directly (it goes
#     through `_tendency_and_jacobian`/`bulk_microphysics_tendencies`), so no call site
#     needed renaming, only the header comment did.
#   - The "_jacobian_2mp3_manual Tier-1 closed-form entries match ForwardDiff" testset's
#     condensation and ice-deposition sub-blocks called
#     `CMNonEq.conv_q_vap_to_q_lcl(CMP.CloudLiquidFormation(τ), nothing, tps, ...)` /
#     `conv_q_vap_to_q_icl(CMP.ConstantTimescale(τ), nothing, ...)`. Both
#     `CloudLiquidFormation` and `ConstantTimescale` are now fieldless marker structs
#     (`src/parameters/Microphysics1MOptions.jl`) whose dispatch reads `τ` from
#     `mp.process_params.cloud_liquid_formation.τ_relax` / `.cloud_ice_formation.τ_relax`
#     on a REAL `Microphysics1MParams`, not from a constructor argument on a `nothing`
#     `mp` - `CloudLiquidFormation(τ)` itself throws (no such method on a zero-field
#     struct). DROPPED, not repaired in place: the condensation half of this comparison
#     is already superseded by the "card #23/#26" testset below, which differentiates
#     the 2M+P3's OWN bare-rate function (`_bare_rate_conv_q_vap_to_q_lcl`) directly
#     rather than borrowing the 1M scheme's folded-rate function as a stand-in
#     reference. The ice half had no such modern equivalent in this file; one is added
#     below ("ice deposition/sublimation entries match full AD"), mirroring card
#     #23/#26's pattern exactly, on the ice side, replacing both the ice condensation
#     sub-block and the ice-number-sublimation-pathway (`n_sub`) sub-block it also
#     depended on the same stale calls for.
#   - The F23 (deposition-nucleation number pathway) sub-block is DROPPED outright, not
#     repaired at the time this file was FIRST written: `_jacobian_2mp3_manual`
#     (`src/BMT_2mp3_jacobian.jl:451`) read `mp.ice.inp_depletion_model.τ_act`, and the
#     default depletion model (`NIceProxyDepletion`) is a fieldless marker struct with
#     no such field. THIS IS NOW FIXED - `_jacobian_2mp3_manual` evaluates
#     `CM_HetIce.delivery_rate(mp.ice.ice_nucleation, mp, tps, T, S_i)` (the same
#     state-dependent rate the primal uses) instead of reading a stored constant, per
#     the team lead's confirmation. A fresh "deposition-nucleation number pathway
#     matches delivery_rate" testset was written below once the fix landed, using
#     `deposition_rate`'s current five-slot interface (`inp, mp, tps, micro, thermo`).
#     The branch is being restacked continuously as defects surface; this is a snapshot
#     of the fix, not a standing guarantee - re-check before relying on it.

function test_framework_2m(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    T_frz = TDI.T_freeze(tps)

    consistent_logλ(ρ, x) =
        P3.get_distribution_logλ(P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]))

    # x = [q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim]
    ρ = FT(0.78)
    T = FT(273.5)
    q_tot = FT(0.009)
    x = FT[2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8]
    logλ = consistent_logλ(ρ, x)

    @testset "rosenbrock_exact() on Microphysics2Moment works ($FT)" begin
        for nsub in (1, 2, 8), Δt in (FT(60), FT(300))
            t = BMT.bulk_microphysics_tendencies(
                BMT.rosenbrock_exact(), BMT.Microphysics2Moment(), mp, tps,
                ρ, T, q_tot, x..., logλ, Δt, nsub,
            )
            @test all(isfinite, Base.front(values(t)))  # `extras` is a NamedTuple, not a number
        end
    end

    @testset "rosenbrock_manual() on Microphysics2Moment works ($FT)" begin
        for nsub in (1, 2, 8), Δt in (FT(60), FT(300))
            t = BMT.bulk_microphysics_tendencies(
                BMT.rosenbrock_manual(), BMT.Microphysics2Moment(), mp, tps,
                ρ, T, q_tot, x..., logλ, Δt, nsub,
            )
            @test all(isfinite, Base.front(values(t)))
        end
    end

    @testset "_numadj_derivs matches ForwardDiff of number_tendency_from_mass_limits ($FT)" begin
        # This sub-block of the campaign's "Tier-1 closed-form entries" testset is
        # unaffected by the stale-API issues above: `CM2.number_tendency_from_mass_limits`
        # is called directly, with its own current default arguments
        # (`sat_excess = 0`, `invent_from_zero = true`), across the
        # interior/low/high/empty/trace regimes for the cloud, rain, and (a literal)
        # ice-shaped bound.
        rtol = FT == Float64 ? FT(1e-9) : FT(1e-2)
        atol = FT == Float64 ? FT(1e-12) : FT(1e-6)
        qmin = UT.ϵ_numerics_2M_M(FT)
        sb = mp.warm_rain.seifert_beheng
        numadj_species = (
            (sb.pdf_c.xc_min, sb.pdf_c.xc_max, sb.numadj.τ),
            (sb.pdf_r.xr_min, sb.pdf_r.xr_max, sb.numadj.τ),
            (FT(1e-12), FT(1e-5), FT(100)),
        )
        for (x_min, x_max, τ) in numadj_species
            q = FT(2e-4)
            n_mid = q / sqrt(x_min * x_max)
            for (q_test, n_test) in (
                (q, q / x_max * FT(0.5)),  # low clamp
                (q, n_mid),                # interior
                (q, q / x_min * FT(2)),    # high clamp
                (zero(FT), n_mid),         # empty: the presence test is `q > 0`
                (qmin / 2, n_mid),         # trace but present: no longer the empty arm
            )
                manual = collect(BMT._numadj_derivs(FT, q_test, n_test, x_min, x_max, τ))
                h(v) = CM2.number_tendency_from_mass_limits((; x_min, x_max, τ), v[1], v[2])
                fd = FD.gradient(h, FT[q_test, n_test])
                # the mass coupling is treated explicitly; the implicit part is
                # the relaxation diagonal
                @test manual[1] == 0
                @test isapprox(manual[2], fd[2]; rtol, atol)
            end
        end
    end

    @testset "card #23/#26: the bare-rate condensation self-entry matches full AD ($FT)" begin
        # Card #23's falsifier ladder item 1 (unit-level exact identity), updated for the
        # LANDING's bare rate - the design-phase check in
        # notes/condensation-tau-derivative-liquid-design.md used the FOLDED CMNonEq helper as
        # ground truth; the landing's primal is bare, so the reference here differentiates
        # `_bare_rate_conv_q_vap_to_q_lcl` instead. Two checks per state:
        #   (A) the bare-rate self-entry ALONE (`_condevap_derivs` with Γ=1, dcp=0,
        #       `dlog_τ_dq_limit=0`, i.e. τ held fixed) against AD of the bare rate at fixed τ -
        #       isolates card #23's own contribution;
        #   (B) the FULL production entry (card #23 + card #26's `dlog_τ_dq_liq` together)
        #       against a FULL AD reference that differentiates the bare rate ALL THE WAY
        #       THROUGH `τ_l`'s own `q_lcl` dependence (`cloud_condensation_timescale` called
        #       inside the differentiated closure) - the true total derivative, not a
        #       leading-order approximation.
        # States span the vapor branch at three supersaturations (low/mid/the small-q_lcl,
        # high-S regime `notes/condensation-tau-derivative-liquid-design.md` section 5 names as
        # its worst measured corpus ratio) and one evaporation-LIMITED state.
        sb = mp.warm_rain.seifert_beheng
        aps = mp.warm_rain.air_properties
        rtol_id = FT == Float64 ? FT(1e-10) : FT(1e-4)
        condevap_id_states = (
            (; ρ = FT(1.1), T = FT(290.0), q_tot = FT(0.0165), q_lcl = FT(5e-4), n_lcl = FT(8e7), q_rai = FT(1e-5)),
            (; ρ = FT(1.1), T = FT(290.0), q_tot = FT(0.0180), q_lcl = FT(5e-4), n_lcl = FT(8e7), q_rai = FT(1e-5)),
            (; ρ = FT(0.95), T = FT(292.0), q_tot = FT(0.0175), q_lcl = FT(8e-5), n_lcl = FT(5e7), q_rai = FT(1e-5)),
            (; ρ = FT(1.0), T = FT(288.0), q_tot = FT(0.008), q_lcl = FT(2e-6), n_lcl = FT(1e6), q_rai = FT(1e-6)),
        )
        for st in condevap_id_states
            (; ρ, T, q_tot, q_lcl, n_lcl, q_rai) = st
            q_ice = zero(FT)
            sat_excess_arm = BMT._liquid_sat_excess(tps, ρ, T, q_tot, q_lcl, q_rai, q_ice)
            n_lcl_j = CM2.number_bounded_by_mass_limits(
                (; x_min = sb.pdf_c.xc_min, x_max = sb.pdf_c.xc_max), q_lcl, n_lcl, sat_excess_arm;
                invent_from_zero = false)
            @test n_lcl_j == n_lcl   # precondition: not on the mean-mass floor
            τ_l = CM2.cloud_condensation_timescale(sb.pdf_c, aps, tps, T, ρ, q_lcl, n_lcl_j * ρ)
            @test !CM2.cloud_condensation_is_degenerate(τ_l)   # precondition: not at the cap
            qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_ice)
            qᵥ_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
            sat_excess_l = qᵥ - qᵥ_sat_liq

            # (A) bare rate alone, τ fixed
            cl_bare = BMT._condevap_derivs(τ_l, sat_excess_l, one(FT), FT(1), FT(1), FT(1),
                q_lcl, false, zero(FT), zero(FT))
            fd_bare(q_lcl_v) = BMT._bare_rate_conv_q_vap_to_q_lcl(τ_l, tps,
                (; q_tot, q_lcl = q_lcl_v, q_icl = q_ice, q_rai, q_sno = zero(q_lcl_v)), (; ρ, T))
            @test isapprox(cl_bare.∂s_liq, FD.derivative(fd_bare, q_lcl); rtol = rtol_id)

            # (B) the full production entry, τ = τ_l(q_lcl) differentiated through
            dlog_τ_dq_liq = UT.guarded_quotient(-one(FT), 3 * q_lcl)
            cl_full = BMT._condevap_derivs(τ_l, sat_excess_l, one(FT), FT(1), FT(1), FT(1),
                q_lcl, false, zero(FT), zero(FT), dlog_τ_dq_liq)
            function fd_full(q_lcl_v)
                sat_excess_arm_v = BMT._liquid_sat_excess(tps, ρ, T, q_tot, q_lcl_v, q_rai, q_ice)
                n_lcl_j_v = CM2.number_bounded_by_mass_limits(
                    (; x_min = sb.pdf_c.xc_min, x_max = sb.pdf_c.xc_max), q_lcl_v, n_lcl, sat_excess_arm_v;
                    invent_from_zero = false)
                τ_v = CM2.cloud_condensation_timescale(sb.pdf_c, aps, tps, T, ρ, q_lcl_v, n_lcl_j_v * ρ)
                BMT._bare_rate_conv_q_vap_to_q_lcl(τ_v, tps,
                    (; q_tot, q_lcl = q_lcl_v, q_icl = q_ice, q_rai, q_sno = zero(q_lcl_v)), (; ρ, T))
            end
            @test isapprox(cl_full.∂s_liq, FD.derivative(fd_full, q_lcl); rtol = rtol_id)

            # NOT compared against the full `J8[1,1]`: `lcl_lcl` sums every process with a q_lcl
            # donor diagonal (autoconversion, accretion, activation, self-collection, ... - Tier
            # 2/3 approximations by design), not condensation alone. A first version of this test
            # made that comparison and measured only ~6e-5 ABSOLUTE agreement at a state where the
            # condensation term itself is O(0.1) - not because `cl.∂s_liq` is wrong, but because
            # it is answering a different question than "what is J[1,1]" (the same trap
            # notes/condensation-tau-derivative-liquid-design.md section 2 names and deliberately
            # does not chase). `cl.∂s_liq` is exercised directly above; a same-state cross-check
            # against `_jacobian_2mp3_manual`'s own `cl` local would need an accessor this
            # function does not expose, and is not worth adding for this alone.
        end
    end

    @testset "ice deposition/sublimation entries match full AD ($FT)" begin
        # Ice-side counterpart of "card #23/#26" above, replacing the campaign's
        # stale ice condensation and ice-number-sublimation ("n_sub") sub-blocks (see the
        # file header). Two checks:
        #   (A) the bare-rate ice self-entry (`_condevap_derivs` at fixed τ_i) against AD
        #       of `_bare_rate_conv_q_vap_to_q_icl` at fixed τ_i - the closed-form
        #       `∂s_ice` self-derivative `_jacobian_2mp3_manual`'s `ice_ice` entry uses;
        #   (B) the ice-number sublimation pathway `n_per_q * (∂s_ice - sub_frac)`
        #       (`_jacobian_2mp3_manual`'s `nice_ice`/`nice_nice` construction) against AD
        #       of the number tendency `n_ice * (∂ₜq_ice_dep / q_ice)` on the sublimation
        #       branch, at the SAME fixed τ_i (the manual entry does not differentiate
        #       through τ_i's own `n_ice` dependence in the donor cross-terms either; only
        #       its own `n_ice`-proportional self-term is carried exactly, which is a
        #       separate, already-covered claim - see `_jacobian_2mp3_manual`'s own
        #       docstring on `nice_nice`).
        # τ_i is the REAL, state-dependent `P3.ice_deposition_timescale`, not a config
        # constant: unlike the 1M scheme, the 2M+P3 primal has no fixed ice condensation
        # timescale to borrow, so τ_i is evaluated once per state and held fixed while
        # `q_ice` is differentiated, matching what `_jacobian_2mp3_manual` itself does
        # (its own `τ_i` is built from the SAME clamped state before `_condevap_derivs`
        # is called, at every state, not re-solved inside the derivative).
        rtol_id = FT == Float64 ? FT(1e-10) : FT(1e-4)
        qmin = UT.ϵ_numerics_2M_M(FT)
        ice_dep_states = (
            (; ρ = FT(0.9), T = FT(250.0), q_tot = FT(0.0008), q_lcl = FT(0.0), q_rai = FT(0.0),
                q_ice = FT(1e-4), n_ice = FT(1e5)),   # supersaturated over ice: growth branch
            (; ρ = FT(0.9), T = FT(270.0), q_tot = FT(0.0002), q_lcl = FT(0.0), q_rai = FT(0.0),
                q_ice = FT(1e-4), n_ice = FT(1e5)),   # subsaturated: sublimation branch
        )
        for st in ice_dep_states
            (; ρ, T, q_tot, q_lcl, q_rai, q_ice, n_ice) = st
            state = P3.state_from_prognostic(p3, ρ * q_ice, ρ * n_ice, zero(ρ), zero(ρ))
            logλ_st = P3.get_distribution_logλ(state)
            τ_i = P3.ice_deposition_timescale(
                mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps, T, ρ, state, logλ_st;
                quad = mp.ice.quad,
            )
            @test isfinite(τ_i) && τ_i > 0
            @test !P3.ice_deposition_is_degenerate(τ_i)  # precondition: not at the cap

            qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_ice)
            qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
            sat_excess_i = qᵥ - qᵥ_sat_ice

            # (A) bare rate alone, τ_i fixed
            ci_bare = BMT._condevap_derivs(τ_i, sat_excess_i, one(FT), FT(1), FT(1), FT(1),
                q_ice, true, zero(FT), zero(FT))
            fd_bare(q_ice_v) = BMT._bare_rate_conv_q_vap_to_q_icl(τ_i, tps,
                (; q_tot, q_lcl, q_icl = q_ice_v, q_rai, q_sno = zero(q_ice_v)), (; ρ, T))
            @test isapprox(ci_bare.∂s_ice, FD.derivative(fd_bare, q_ice); rtol = rtol_id)

            # (B) the ice-number sublimation pathway, at the same fixed τ_i
            ∂ₜq_ice_dep = fd_bare(q_ice)  # the bare rate's own value at this state
            dep_active = !(sat_excess_i > 0 && T > TDI.T_freeze(tps))
            n_sub_active = dep_active && ∂ₜq_ice_dep < 0
            if n_sub_active
                sub_frac = ∂ₜq_ice_dep / max(q_ice, floatmin(FT))
                n_per_q = n_ice / max(qmin, q_ice)
                manual_nice_ice = n_per_q * (ci_bare.∂s_ice - sub_frac)
                fd_n(q_ice_v) = begin
                    r = BMT._bare_rate_conv_q_vap_to_q_icl(τ_i, tps,
                        (; q_tot, q_lcl, q_icl = q_ice_v, q_rai, q_sno = zero(q_ice_v)), (; ρ, T))
                    ifelse(r < 0, n_ice * r / q_ice_v, zero(r))
                end
                @test isapprox(manual_nice_ice, FD.derivative(fd_n, q_ice); rtol = rtol_id)
            end
        end
        # confirm both branches (growth and sublimation) were exercised by the corpus above
        @test any(
            begin
                (; ρ, T, q_tot, q_lcl, q_rai, q_ice) = st
                qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_ice)
                qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
                (qᵥ - qᵥ_sat_ice) < 0
            end for st in ice_dep_states
        )
    end

    @testset "deposition-nucleation number pathway matches delivery_rate ($FT)" begin
        # Written fresh once the `τ_act` field-access bug (see the file header) was fixed:
        # `_jacobian_2mp3_manual`'s deposition-nucleation entry now evaluates
        # `CM_HetIce.delivery_rate(mp.ice.ice_nucleation, mp, tps, T, S_i)` (the SAME
        # state-dependent rate the primal uses) rather than reading a stored constant, so
        # this is a value-level structural check that the Jacobian entry equals what it
        # claims to, NOT an AD-derivative comparison: the entry deliberately does not carry
        # `S_i`'s own dependence on the donors (per the team lead: "the supersaturation
        # dependence is still NOT carried in the Jacobian... the block keeps the number
        # diagonal and deliberately drops the mass-donor coupling"), so differentiating
        # through `delivery_rate` would be asserting something the design does not claim.
        # An assertion about the deposition slot's behavior NEAR saturation belongs against
        # the rate itself (the boundary battery), not against this Jacobian entry.
        #
        # The all-zero state isolates the entry cleanly: with q_ice = n_ice = 0, ice
        # aggregation and the ice-number sublimation pathway are both inactive, and the
        # ice-number-adjustment's own contribution at zero mass and zero number is the
        # SEPARATE, closed-form `-1/τ_numadj` term `_ice_numadj_params` already covers
        # (`test_ice_numadj_timescale_key`, below), so `nice_nice` here is exactly the sum
        # of those two known pieces and nothing else.
        ρ_dep, q_tot_dep = FT(0.6), FT(4.0e-3)
        empty = BMT.MicroState2MP3{FT}(0, 0, 0, 0, 0, 0, 0, 0)
        for T_dep in FT[220, 233, 245, 253]  # below the Cooper gate's 258.15 K threshold
            logλ_dep = P3.get_distribution_logλ(P3.state_from_prognostic(p3, FT(0), FT(0), FT(0), FT(0)))
            g_dep = BMT.Instantaneous2MP3Tendency(mp, tps, ρ_dep, T_dep, q_tot_dep, logλ_dep)
            micro_dep = (;
                q_tot = q_tot_dep, q_lcl = FT(0), n_lcl = FT(0), q_rai = FT(0), n_rai = FT(0),
                q_ice = FT(0), n_ice = FT(0), q_rim = FT(0), b_rim = FT(0),
            )
            thermo_dep = (; ρ = ρ_dep, T = T_dep, w = FT(0), p = FT(0), logλ = logλ_dep)
            pp, rs = BMT.p3_2m_process_rates(mp, tps, micro_dep, thermo_dep)
            dep_active = pp.ice_deposition.n_ice > 0
            J = BMT._jacobian_2mp3_manual(g_dep, empty, pp, rs)
            if dep_active
                S_i = TDI.supersaturation_over_ice(tps, q_tot_dep, FT(0), FT(0), ρ_dep, T_dep)
                inv_τ_dep = CM_HetIce.delivery_rate(mp.ice.ice_nucleation, mp, tps, T_dep, S_i)
                @test inv_τ_dep > 0
                expected = -inv_τ_dep - 1 / BMT._ice_numadj_params(p3).τ
                @test isapprox(J[6, 6], expected; rtol = FT(1e-9))
            end
        end
        # confirm the gate actually opened somewhere in the sweep, or the testset is vacuous
        @test any(FT[220, 233, 245, 253]) do T_dep
            S_i = TDI.supersaturation_over_ice(tps, q_tot_dep, FT(0), FT(0), ρ_dep, T_dep)
            CM_HetIce.delivery_rate(mp.ice.ice_nucleation, mp, tps, T_dep, S_i) > 0
        end
    end

    @testset "non-Exact RosenbrockAverage on Microphysics2Moment throws ($FT)" begin
        for mode in (BMT.RosenbrockAverage(), BMT.rosenbrock_coupled())
            @test_throws ArgumentError BMT.bulk_microphysics_tendencies(
                mode, BMT.Microphysics2Moment(), mp, tps,
                ρ, T, q_tot, x..., logλ, FT(60), 4,
            )
        end
    end

    @testset "RosenbrockAverage on warm-rain-only parameters throws ($FT)" begin
        mp_warm = CMP.Microphysics2MParams(FT; with_ice = false)
        @test_throws "requires P3 ice parameters" BMT.bulk_microphysics_tendencies(
            BMT.rosenbrock_exact(), BMT.Microphysics2Moment(), mp_warm, tps,
            ρ, T, q_tot, x..., logλ, FT(60), 4,
        )
    end
end

function _framework_exact_args(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    st = P3.state_from_prognostic(
        p3, FT(0.78) * FT(1e-4), FT(0.78) * FT(2e5), FT(0.78) * FT(4e-5), FT(0.78) * FT(6e-8),
    )
    logλ = P3.get_distribution_logλ(st)
    return (
        BMT.rosenbrock_exact(), BMT.Microphysics2Moment(), mp, tps,
        FT(0.78), FT(273.5), FT(0.009),
        FT(2e-4), FT(5e7), FT(1e-4), FT(4e4), FT(1e-4), FT(2e5), FT(4e-5), FT(6e-8),
        logλ, FT(60), 4,
    )
end

"""
    test_ice_numadj_timescale_key(FT)

`ParametersP3.τ_numadj` reads the `P3_ice_number_adjustment_timescale` TOML key, and
the value reaches `_ice_numadj_params`'s number-adjustment rate.
"""
function test_ice_numadj_timescale_key(FT)
    p3 = CMP.ParametersP3(CP.create_toml_dict(FT))
    override = joinpath(pkgdir(CM), "src", "parameters", "toml", "P3_numadj_override.toml")
    p3_override = CMP.ParametersP3(CP.create_toml_dict(FT; override_file = override))

    @testset "ice number-adjustment timescale reads its own TOML key ($FT)" begin
        @test p3.τ_numadj == FT(100)
        @test p3_override.τ_numadj == FT(50)

        numadj = BMT._ice_numadj_params(p3)
        numadj_override = BMT._ice_numadj_params(p3_override)
        @test numadj.τ == FT(100)
        @test numadj_override.τ == FT(50)
        @test numadj.x_min == numadj_override.x_min
        @test numadj.x_max == numadj_override.x_max

        # Below the low bound, so the target clamps to `q / x_max` and the rate is a
        # nonzero constant divided by `τ`.
        q = FT(2e-4)
        n = q / numadj.x_max * FT(0.5)
        rate = CM2.number_tendency_from_mass_limits(numadj, q, n; invent_from_zero = false)
        rate_override = CM2.number_tendency_from_mass_limits(numadj_override, q, n; invent_from_zero = false)
        @test rate > 0
        @test isapprox(rate_override, 2 * rate; rtol = 10 * eps(FT))
    end
end

function test_framework_exact_inference(FT)
    args = _framework_exact_args(FT)
    @testset "rosenbrock_exact() inference and allocations ($FT)" begin
        @test (@inferred BMT.bulk_microphysics_tendencies(args...)) isa NamedTuple
        JET.@test_opt BMT.bulk_microphysics_tendencies(args...)
        trail = BT.@benchmark $(splat(BMT.bulk_microphysics_tendencies))($args) samples = 100 evals = 1
        @test trail.memory == 0
    end
end

test_framework_2m(Float64)
test_framework_2m(Float32)

test_ice_numadj_timescale_key(Float64)
test_ice_numadj_timescale_key(Float32)

test_framework_exact_inference(Float64)
test_framework_exact_inference(Float32)
