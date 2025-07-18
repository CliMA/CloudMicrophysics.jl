import Test as TT
import BenchmarkTools as BT
import JET

import ClimaParams as CP

import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.ArtifactCalling as AFC
import CloudMicrophysics.Common as CO
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.HetIceNucleation as CMI_het
import CloudMicrophysics.HomIceNucleation as CMI_hom
import CloudMicrophysics.MicrophysicsNonEq as CMN
import CloudMicrophysics.Microphysics0M as CM0
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Nucleation as HN
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.CloudDiagnostics as CMD

function bench_press(
    type,
    foo,
    args,
    min_run_time,
    min_memory = 0.0,
    min_allocs = 0.0,
)
    @info "Testing $foo"

    # converts a string argument to tuple
    args isa String && (args = (args,))

    # Benchmark foo
    trail = BT.@benchmark $foo($args...) samples = 100 evals = 100

    show(stdout, MIME("text/plain"), trail)
    println("\n")

    TT.@test BT.minimum(trail).time < min_run_time
    TT.@test trail.memory <= min_memory
    TT.@test trail.allocs <= min_allocs

    # Test that foo is free from optimization failures
    # and unresolved method dispatches
    JET.@test_opt(foo(args...))

    # Test that the return type of foo is of correct type (Float32 vs Float64)
    TT.@test typeof(foo(args...)) == type
end

function benchmark_test(FT)

    # Artifact calling
    bench_press(
        String,
        AFC.AIDA_ice_nucleation,
        ("in05_17_aida.edf"),
        50000,
        30000,
        300,
    )

    # 0-moment microphysics
    p0m = CMP.Parameters0M(FT)
    # 1-moment microphysics
    liquid = CMP.CloudLiquid(FT)
    ice = CMP.CloudIce(FT)
    rain = CMP.Rain(FT)
    snow = CMP.Snow(FT)
    ce = CMP.CollisionEff(FT)
    # 2-moment microphysics
    override_file = joinpath(
        pkgdir(CM),
        "src",
        "parameters",
        "toml",
        "SB2006_limiters.toml",
    )
    toml_dict = CP.create_toml_dict(FT; override_file)
    sb2006 = CMP.SB2006(toml_dict)
    sb2006_no_limiters = CMP.SB2006(toml_dict, false)

    # P3 scheme
    p3 = CMP.ParametersP3(FT)
    # terminal velocity
    blk1mvel = CMP.Blk1MVelType(FT)
    sb2006vel = CMP.SB2006VelType(FT)
    ch2022 = CMP.Chen2022VelType(FT)
    # aerosol activation
    ap = CMP.AerosolActivationParameters(FT)
    # ice nucleation
    desert_dust = CMP.DesertDust(FT)
    kaolinite = CMP.Kaolinite(FT)
    ip = CMP.IceNucleationParameters(FT)
    H2SO4_prs = CMP.H2SO4SolutionParameters(FT)
    ip_frostenberg = CMP.Frostenberg2023(FT)
    # aerosol nucleation parameters
    mixed_nuc = CMP.MixedNucleationParameters(FT)
    h2so4_nuc = CMP.H2S04NucleationParameters(FT)
    organ_nuc = CMP.OrganicNucleationParameters(FT)
    # air and thermodynamics parameters
    aps = CMP.AirProperties(FT)
    wtr = CMP.WaterProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    ρ_air = FT(1.2)
    T_air = FT(280)
    q_liq = FT(5e-4)
    q_ice = FT(5e-4)
    q_tot = FT(1e-3)
    q = TDI.TD.PhasePartition(q_tot)
    q_rai = FT(1e-4)
    q_sno = FT(1e-4)
    N_liq = FT(1e8)
    N_rai = FT(1e8)
    e = FT(600)

    ρ_r = FT(400.0)
    F_rim = FT(0.95)
    F_liq = FT(0.33)
    N = FT(1e8)
    L = ρ_air * q_ice

    T_air_2 = FT(250)
    RH_2 = FT(1.5)
    T_air_cold = FT(230)
    S_ice = FT(1.2)
    dSi_dt = FT(0.05)

    w_air = FT(0.5)
    p_air = FT(1e5)
    r_aer = FT(0.05 * 1e-6)
    σ_aer = FT(2)
    N_aer = FT(100.0 * 1e6)
    M_seasalt = FT(0.058443)
    κ_seasalt = FT(1.12)
    seasalt_mode = AM.Mode_κ(
        r_aer,
        σ_aer,
        N_aer,
        (FT(1),),
        (FT(1),),
        (M_seasalt,),
        (κ_seasalt,),
    )
    aer_distr = AM.AerosolDistribution((seasalt_mode,))

    x_sulph = FT(0.1)
    Delta_a_w = FT(0.27)
    r_liq = FT(1e-6)
    V_liq = FT(4 / 3 * π * r_liq^3)
    Δt = FT(25)

    INPC = FT(1e5)

    # P3 scheme - TODO - bring them back once we optimize P3 scheme
    #bench_press(P3.thresholds, (p3, ρ_r, F_rim), 12e6, 2048, 80)
    #if FT == Float64
    #    bench_press(
    #        P3.distribution_parameter_solver,
    #        (p3, q_ice, N, ρ_r, F_r),
    #        1e5,
    #    )
    #    bench_press(
    #        P3.ice_terminal_velocity,
    #        (p3, ch2022.snow_ice, q_ice, N, ρ_r, F_r, ρ_air, false),
    #        2.5e5,
    #        800,
    #        4,
    #    )
    #    bench_press(
    #        P3.ice_terminal_velocity,
    #        (p3, ch2022.snow_ice, q_ice, N, ρ_r, F_r, ρ_air, true),
    #        2.5e5,
    #        800,
    #        4,
    #    )
    #    bench_press(P3.D_m, (p3, q_ice, N, ρ_r, F_r), 1e5)
    #end
    # P3 ice nucleation
    #bench_press(
    #    P3.het_ice_nucleation,
    #    (
    #        kaolinite,
    #        tps,
    #        TDI.TD.PhasePartition(q_tot, q_liq, q_ice),
    #        N_liq,
    #        RH_2,
    #        T_air_2,
    #        ρ_air,
    #        Δt,
    #    ),
    #    200,
    #)
    #if FT == Float64
    #    bench_press(
    #        P3.ice_melt,
    #        (p3, ch2022.snow_ice, aps, tps, L, N, T_air, ρ_air, F_rim, ρ_r, Δt),
    #        3.7e5,
    #        2e3,
    #        3,
    #    )
    #end
    #bench_press(CMI_het.P3_deposition_N_i, (ip.p3, T_air_cold), 230)
    #bench_press(CMI_het.P3_het_N_i, (ip.p3, T_air_cold, N_liq, V_liq, Δt), 230)

    # Chen 2022 terminal velocity
    bench_press(FT, CMN.terminal_velocity, (liquid, ch2022.rain, ρ_air, q_liq), 350)
    bench_press(
        FT,
        CMN.terminal_velocity,
        (ice, ch2022.small_ice, ρ_air, q_ice),
        400,
    )
    bench_press(FT, CM1.terminal_velocity, (rain, ch2022.rain, ρ_air, q_rai), 850)
    bench_press(
        FT,
        CM1.terminal_velocity,
        (snow, ch2022.large_ice, ρ_air, q_sno),
        850,
    )

    # aerosol activation
    bench_press(
        FT,
        AA.total_N_activated,
        (ap, aer_distr, aps, tps, T_air, p_air, w_air, q_tot, FT(0), FT(0)),
        1300,
    )

    # Common
    bench_press(
        FT,
        CO.H2SO4_soln_saturation_vapor_pressure,
        (H2SO4_prs, x_sulph, T_air_cold),
        50,
    )
    bench_press(FT, CO.a_w_xT, (H2SO4_prs, tps, x_sulph, T_air_cold), 230)
    bench_press(FT, CO.a_w_eT, (tps, e, T_air_cold), 230)
    bench_press(FT, CO.a_w_ice, (tps, T_air_cold), 230)

    # ice nucleation
    bench_press(
        FT,
        CMI_het.dust_activated_number_fraction,
        (desert_dust, ip.deposition, S_ice, T_air_2),
        50,
    )
    bench_press(
        FT,
        CMI_het.MohlerDepositionRate,
        (desert_dust, ip.deposition, S_ice, T_air_2, dSi_dt, N_aer),
        80,
    )
    bench_press(FT, CMI_het.deposition_J, (kaolinite, Delta_a_w), 230)
    bench_press(FT, CMI_het.ABIFM_J, (desert_dust, Delta_a_w), 230)
    bench_press(
        FT,
        CMI_het.INP_concentration_frequency,
        (ip_frostenberg, INPC, T_air_cold),
        150,
    )
    bench_press(FT, CMI_hom.homogeneous_J_cubic, (ip.homogeneous, Delta_a_w), 230)
    bench_press(FT, CMI_hom.homogeneous_J_linear, (ip.homogeneous, Delta_a_w), 230)

    # non-equilibrium
    bench_press(FT, CMN.τ_relax, (liquid,), 15)
    bench_press(FT, CMN.conv_q_vap_to_q_liq_ice, (ice, FT(2e-3), FT(1e-3)), 15)
    bench_press(
        FT,
        CMN.conv_q_vap_to_q_liq_ice_MM2015,
        (liquid, tps, FT(0.00145), FT(0), FT(0), FT(0), FT(0), FT(0.8), FT(263)),
        70,
    )


    # 0-moment
    bench_press(FT, CM0.remove_precipitation, (p0m, q), 12)
    bench_press(FT, CM0.remove_precipitation, (p0m, q_liq, q_ice), 12)

    # 1-moment
    bench_press(
        FT,
        CM1.accretion,
        (liquid, rain, blk1mvel.rain, ce, q_liq, q_rai, ρ_air),
        350,
    )
    bench_press(FT, CMD.radar_reflectivity_1M, (rain, q_rai, ρ_air), 300)

    # 2-moment
    for sb in [sb2006, sb2006_no_limiters]
        bench_press(
            @NamedTuple{au::CM2.LiqRaiRates{FT}, sc::FT},
            CM2.autoconversion_and_liquid_self_collection,
            (sb, q_liq, q_rai, ρ_air, N_liq),
            300,
        )
        bench_press(
            @NamedTuple{sc::FT, br::FT},
            CM2.rain_self_collection_and_breakup,
            (sb, q_rai, ρ_air, N_rai),
            1200,
        )
        bench_press(
            @NamedTuple{evap_rate_0::FT, evap_rate_1::FT},
            CM2.rain_evaporation,
            (sb, aps, tps, q_tot, q_liq, q_ice, q_rai, q_sno, ρ_air, N_rai, T_air),
            2000,
        )
        bench_press(
            Tuple{FT, FT},
            CM2.rain_terminal_velocity,
            (sb, sb2006vel, q_rai, ρ_air, N_rai),
            700,
        )
        bench_press(
            Tuple{FT, FT},
            CM2.rain_terminal_velocity,
            (sb, ch2022.rain, q_rai, ρ_air, N_rai),
            2200,
        )
        bench_press(
            FT,
            CMD.radar_reflectivity_2M,
            (sb, q_liq, q_rai, N_liq, N_rai, ρ_air),
            2000,
        )
        bench_press(
            FT,
            CMD.effective_radius_2M,
            (sb, q_liq, q_rai, N_liq, N_rai, ρ_air),
            2000,
        )
    end
    bench_press(
        FT,
        CMD.effective_radius_Liu_Hallet_97,
        (wtr, ρ_air, q_liq, N_liq, q_rai, N_rai),
        300,
    )
    bench_press(
        FT,
        CM2.number_increase_for_mass_limit,
        (sb2006.numadj, FT(5e-6), q_rai, ρ_air, N_rai),
        50,
    )
    bench_press(
        FT,
        CM2.number_decrease_for_mass_limit,
        (sb2006.numadj, FT(2.6e-10), q_rai, ρ_air, N_rai),
        50,
    )
    # Homogeneous Nucleation
    bench_press(
        @NamedTuple{binary_rate::FT, ternary_rate::FT},
        HN.h2so4_nucleation_rate,
        (FT(1e12), FT(1), FT(1), FT(208), h2so4_nuc),
        470,
    )
    bench_press(
        FT,
        HN.organic_nucleation_rate,
        (FT(0), FT(1e3), FT(1e3), FT(1e3), FT(300), FT(1), organ_nuc),
        850,
    )
    bench_press(
        FT,
        HN.organic_and_h2so4_nucleation_rate,
        (FT(2.6e6), FT(1), FT(1), FT(300), FT(1), mixed_nuc),
        120,
    )
end

TT.@testset "Performance Tests ($FT)" for FT in (Float64, Float32)
    benchmark_test(FT)
end
nothing
