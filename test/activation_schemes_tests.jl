using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.AerosolModel as CMAM
import CloudMicrophysics.ThermodynamicsInterface as TDI

"""
Unit tests for the pluggable aerosol activation scheme dispatch introduced
in issue 012 (BMT shared-API refactor).

Each test exercises one scheme in isolation via the pure helper
`BMT.activation_source`. Integration with the full BMT tendency tuple
is tested separately by re-running the standard BMT P3 smoke tests with
each scheme plugged in.
"""
function test_activation_schemes(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    # A reasonable cold-cloud thermodynamic state
    ρ = FT(1.0)
    T = FT(275)
    q_lcl = FT(1e-3)           # cloudy
    q_lcl_dry = FT(1e-9)       # no cloud
    q_ice = FT(0)
    q_vap_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
    q_tot_sat = q_vap_sat + q_lcl             # S ≈ 0
    q_tot_super = FT(1.01) * q_vap_sat + q_lcl  # supersaturated
    q_tot_sub = FT(0.9) * q_vap_sat            # subsaturated, no cloud
    n_lcl_empty = FT(0)
    n_lcl_partial = FT(5e7)
    w_up = FT(0.5)
    w_dn = FT(-0.1)
    p = FT(90_000)

    @testset "NoActivation (always zero)" begin
        scheme = CMP.NoActivation()
        @test BMT.activation_source(scheme, tps, ρ, T, q_tot_sat,
            q_lcl, q_ice, n_lcl_empty, w_up, p) == FT(0)
        @test BMT.activation_source(scheme, tps, ρ, T, q_tot_super,
            q_lcl, q_ice, n_lcl_partial, w_up, p) == FT(0)
    end

    @testset "DiagnosticNc" begin
        scheme = CMP.DiagnosticNc{FT}(N_c = FT(1e8), q_thresh = FT(1e-7), τ_relax = FT(60))

        # Cloudy, no droplets yet → positive source relaxing toward N_c
        s_empty = BMT.activation_source(scheme, tps, ρ, T, q_tot_sat,
            q_lcl, q_ice, n_lcl_empty, w_up, p)
        @test s_empty > FT(0)
        @test isapprox(s_empty, FT(1e8) / FT(60); rtol = 1e-12)

        # Cloudy, partially activated → smaller positive source
        s_partial = BMT.activation_source(scheme, tps, ρ, T, q_tot_sat,
            q_lcl, q_ice, n_lcl_partial, w_up, p)
        @test s_partial > FT(0)
        @test s_partial < s_empty
        @test isapprox(s_partial, (FT(1e8) - n_lcl_partial) / FT(60); rtol = 1e-12)

        # Cloudy but already beyond target → negative source
        s_over = BMT.activation_source(scheme, tps, ρ, T, q_tot_sat,
            q_lcl, q_ice, FT(2e8), w_up, p)
        @test s_over < FT(0)

        # No cloud → target is zero, existing droplets relax out
        s_nocloud = BMT.activation_source(scheme, tps, ρ, T, q_tot_sub,
            q_lcl_dry, q_ice, n_lcl_partial, w_up, p)
        @test s_nocloud < FT(0)
        @test isapprox(s_nocloud, -n_lcl_partial / FT(60); rtol = 1e-12)

        # Updraft/pressure are ignored by this tier
        s_dn = BMT.activation_source(scheme, tps, ρ, T, q_tot_sat,
            q_lcl, q_ice, n_lcl_empty, w_dn, p)
        @test s_dn == s_empty
    end

    @testset "TwomeyActivation" begin
        scheme = CMP.TwomeyActivation{FT}(
            C = FT(1e8), k = FT(0.4), w_min = FT(0.01),
            q_thresh = FT(1e-7), τ_relax = FT(60),
        )

        # Cloudy, supersaturated, updraft → positive
        s_super = BMT.activation_source(scheme, tps, ρ, T, q_tot_super,
            q_lcl, q_ice, n_lcl_empty, w_up, p)
        @test s_super > FT(0)

        # Cloudy, saturated only (S ≈ 0) → target S^k ≈ 0 → small/negative
        # relaxation (pulls n_lcl toward 0)
        s_sat = BMT.activation_source(scheme, tps, ρ, T, q_tot_sat,
            q_lcl, q_ice, n_lcl_partial, w_up, p)
        @test s_sat <= FT(0)

        # Subsidence → gate closes, target is zero
        s_dn = BMT.activation_source(scheme, tps, ρ, T, q_tot_super,
            q_lcl, q_ice, n_lcl_partial, w_dn, p)
        @test s_dn < FT(0)      # relaxes existing droplets toward zero
        @test isapprox(s_dn, -n_lcl_partial / FT(60); rtol = 1e-12)

        # No cloud → zero target
        s_dry = BMT.activation_source(scheme, tps, ρ, T, q_tot_sub,
            q_lcl_dry, q_ice, n_lcl_empty, w_up, p)
        @test s_dry == FT(0)
    end

    @testset "FixedARGActivation" begin
        # Build a simple single-mode kappa-distribution (continental-ish)
        act_params = CMP.AerosolActivationParameters(FT)
        air_props = CMP.AirProperties(FT)
        mode = CMAM.Mode_κ(
            FT(5e-8),     # r_dry [m]
            FT(2.0),      # std
            FT(1e9),      # N [m⁻³]
            (FT(1),),     # vol mix
            (FT(1),),     # mass mix
            (FT(0),),     # molar (unused)
            (FT(0.6),),   # κ
        )
        distribution = CMAM.AerosolDistribution((mode,))
        scheme = CMP.FixedARGActivation{FT, typeof(act_params), typeof(air_props), typeof(distribution)}(
            act_params = act_params,
            air_properties = air_props,
            distribution = distribution,
            q_thresh = FT(1e-7),
            τ_relax = FT(60),
        )

        # Cloudy, supersaturated, strong updraft → kernel should produce a
        # positive tendency (or, if the input state happens to yield N_act
        # below n_lcl, zero — check for non-negativity).
        s_super = BMT.activation_source(scheme, tps, ρ, T, q_tot_super,
            q_lcl, q_ice, n_lcl_empty, w_up, p)
        @test s_super >= FT(0)

        # Subsaturated, no cloud → zero
        s_dry = BMT.activation_source(scheme, tps, ρ, T, q_tot_sub,
            q_lcl_dry, q_ice, n_lcl_empty, w_up, p)
        @test s_dry == FT(0)

        # Subsidence → zero
        s_dn = BMT.activation_source(scheme, tps, ρ, T, q_tot_super,
            q_lcl, q_ice, n_lcl_empty, w_dn, p)
        @test s_dn == FT(0)

        # Zero pressure should not throw
        s_bad = BMT.activation_source(scheme, tps, ρ, T, q_tot_super,
            q_lcl, q_ice, n_lcl_empty, w_up, FT(0))
        @test s_bad == FT(0)
    end
end

function test_activation_schemes_in_bmt(FT)
    # Build a P3+warm-rain param bundle and confirm the activation-scheme
    # threads through the 2M+P3 entry point and adds to dn_lcl_dt.
    # The scheme lives inside `mp.warm_rain.activation_scheme`, so we build
    # one parameter bundle per scheme under test.
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp_base = CMP.Microphysics2MParams(FT; with_ice = true)

    _with_activation_scheme(mp, scheme) = CMP.Microphysics2MParams(;
        warm_rain = CMP.WarmRainParams2M(;
            seifert_beheng = mp.warm_rain.seifert_beheng,
            air_properties = mp.warm_rain.air_properties,
            condevap       = mp.warm_rain.condevap,
            subdep         = mp.warm_rain.subdep,
            activation_scheme = scheme,
        ),
        ice = mp.ice,
    )

    ρ = FT(1.0)
    T = FT(275)
    q_vap_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
    q_lcl = FT(1e-3)
    q_ice = FT(0)      # no ice → skip all P3 ice processes
    q_tot = q_vap_sat + q_lcl
    n_lcl = FT(5e7)
    q_rai = FT(0)
    n_rai = FT(0)
    q_rim = FT(0)
    b_rim = FT(0)
    logλ = FT(0)
    n_ice = FT(0)

    scheme_diag = CMP.DiagnosticNc{FT}(N_c = FT(1e8), q_thresh = FT(1e-7), τ_relax = FT(60))
    mp_no   = _with_activation_scheme(mp_base, CMP.NoActivation())
    mp_diag = _with_activation_scheme(mp_base, scheme_diag)

    t_no = BMT.bulk_microphysics_tendencies(
        BMT.Microphysics2Moment(), mp_no, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai,
        q_ice, n_ice, q_rim, b_rim, logλ,
    )
    @test t_no.dn_lcl_activation_dt == FT(0)

    t_diag = BMT.bulk_microphysics_tendencies(
        BMT.Microphysics2Moment(), mp_diag, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai,
        q_ice, n_ice, q_rim, b_rim, logλ,
    )
    @test t_diag.dn_lcl_activation_dt > FT(0)
    # dn_lcl_dt with activation should equal the baseline dn_lcl_dt plus the
    # activation source (other processes unchanged).
    @test isapprox(
        t_diag.dn_lcl_dt,
        t_no.dn_lcl_dt + t_diag.dn_lcl_activation_dt;
        rtol = 1e-12,
    )
end


function test_repair_ice_state(FT)
    # Realistic mid-tropospheric air density so the per-volume convention
    # (N [1/m³], L [kg/m³]) is genuinely exercised — with ρ = 1 the two
    # unit conventions coincide numerically and unit bugs hide.
    ρ = FT(0.8)
    m_max = FT(1e-5)
    q_floor = FT(1e-14)

    @testset "repair_ice_state scalar" begin
        # Convention: N [1/m³], L [kg/m³], ρ [kg/m³]. Mean particle mass is
        # L / N [kg/particle]; healthy if ≤ m_max.

        # Healthy cell: mean-crystal mass well below m_max → no change
        N0 = FT(1e5)
        L0 = FT(1e-6)        # implied mean ≈ 1e-11 kg (fine)
        @test BMT.repair_ice_state(N0, L0, ρ) == N0

        # Pathological cell: mass present, number tiny → repair
        N1 = FT(1e-3)        # essentially no crystals; implied mean ≈ 1e-3 kg (!)
        L1 = FT(1e-6)
        repaired = BMT.repair_ice_state(N1, L1, ρ)
        @test repaired > N1
        @test isapprox(repaired, L1 / m_max; rtol = 1e-12)   # exactly the floor

        # Trace-ice cell (L below q_floor·ρ) → collapse N to zero
        # (top-of-domain phantom-number pathology). Repair is symmetric:
        # mass-without-number on one tail → impute number; number-without-mass
        # on the other tail → zero number.
        N2 = FT(0.01)
        L2 = FT(1e-16)       # < q_floor·ρ = 0.8e-14
        @test BMT.repair_ice_state(N2, L2, ρ) == FT(0)

        # Edge: L just below q_floor·ρ → N → 0; just above triggers
        # upper-mass-bound repair. Probe both branches with the SAME
        # input N so the contrast is clear.
        L_below = q_floor * ρ * FT(0.5)
        L_above = q_floor * ρ * FT(2.0)
        @test BMT.repair_ice_state(FT(1e6), L_below, ρ) == FT(0)
        @test BMT.repair_ice_state(FT(0), L_above, ρ) > FT(0)

        # Already-large N → repair doesn't reduce
        N3 = FT(1e8)
        L3 = FT(1e-5)        # implied mean = 1e-13 kg (fine)
        @test BMT.repair_ice_state(N3, L3, ρ) == N3
    end
end

@testset "Activation schemes (Float64)" begin
    test_activation_schemes(Float64)
    test_activation_schemes_in_bmt(Float64)
    test_repair_ice_state(Float64)
end

@testset "Activation schemes (Float32)" begin
    test_activation_schemes(Float32)
end
nothing
