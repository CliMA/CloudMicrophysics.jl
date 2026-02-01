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
    min_run_time,  # minimum allowed run time, in nanoseconds
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

    @info "Artifacts"
    bench_press(String, AFC.AIDA_ice_nucleation, ("in05_17_aida.edf"), 50_000, 30_000, 300)

    # 0-moment parameters
    p0m = CMP.Parameters0M(FT)

    # 1-moment parameters
    liquid = CMP.CloudLiquid(FT)
    ice = CMP.CloudIce(FT)
    rain = CMP.Rain(FT)
    snow = CMP.Snow(FT)
    ce = CMP.CollisionEff(FT)

    # 2-moment parameters
    override_file = joinpath(pkgdir(CM), "src", "parameters", "toml", "SB2006_limiters.toml")
    toml_dict = CP.create_toml_dict(FT; override_file)
    sb2006 = CMP.SB2006(toml_dict)
    sb2006_no_limiters = CMP.SB2006(toml_dict, false)

    # P3 parameters
    params_P3 = CMP.ParametersP3(FT)

    # Terminal velocity parameters
    blk1mvel = CMP.Blk1MVelType(FT)
    sb2006vel = CMP.SB2006VelType(FT)
    ch2022 = CMP.Chen2022VelType(FT)
    stokes_vel = CMP.StokesRegimeVelType(FT)

    # Aerosol parameters
    ap = CMP.AerosolActivationParameters(FT)

    # Ice nucleation parameters
    desert_dust = CMP.DesertDust(FT)
    kaolinite = CMP.Kaolinite(FT)
    ip = CMP.IceNucleationParameters(FT)
    H2SO4_prs = CMP.H2SO4SolutionParameters(FT)
    ip_frostenberg = CMP.Frostenberg2023(FT)

    # Aerosol nucleation parameters
    mixed_nuc = CMP.MixedNucleationParameters(FT)
    h2so4_nuc = CMP.H2S04NucleationParameters(FT)
    organ_nuc = CMP.OrganicNucleationParameters(FT)

    # Thermodynamics parameters and atmospheric state
    aps = CMP.AirProperties(FT)
    wtr = CMP.WaterProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    ρ_air = FT(1.2)
    T_air = FT(280)
    T_air_2 = FT(250)
    T_air_cold = FT(230)
    p_air = FT(1e5)
    w_air = FT(0.5)
    RH_2 = FT(1.5)
    S_ice = FT(1.2)
    dSi_dt = FT(0.05)

    q_liq = FT(5e-4)
    q_ice = FT(5e-4)
    q_rai = FT(1e-4)
    q_sno = FT(1e-4)
    q_tot = FT(1e-3)

    N_liq = FT(1e8)
    N_rai = FT(1e8)
    N_ice = FT(1e8)
    L_ice = ρ_air * q_ice

    e = FT(600)
    ρ_rim = FT(400)
    F_rim = FT(0.95)

    r_aer = FT(0.05 * 1e-6)
    σ_aer = FT(2)
    N_aer = FT(100.0 * 1e6)
    M_seasalt = FT(0.058443)
    κ_seasalt = FT(1.12)
    seasalt_mode = AM.Mode_κ(r_aer, σ_aer, N_aer, (FT(1),), (FT(1),), (M_seasalt,), (κ_seasalt,))
    aer_distr = AM.AerosolDistribution((seasalt_mode,))

    x_sulph = FT(0.1)
    Delta_a_w = FT(0.27)
    r_liq = FT(1e-6)
    V_liq = FT(4 / 3 * π * r_liq^3)
    Δt = FT(25)
    INPC = FT(1e5)

    @info "P3 Scheme"
    state = P3.get_state(params_P3; L_ice, N_ice, F_rim, ρ_rim)
    logλ = P3.get_distribution_logλ(state)
    bench_press(
        P3.P3State{FT, CMP.ParametersP3{FT, CMP.SlopePowerLaw{FT}}},
        (params, L_ice, N_ice, F_rim, ρ_rim) -> P3.get_state(params; L_ice, N_ice, F_rim, ρ_rim),
        (params_P3, L_ice, N_ice, F_rim, ρ_rim),
        25,
    )
    bench_press(FT, P3.get_distribution_logλ, (state,), 30_000)
    bench_press(FT, P3.get_distribution_logλ, (params_P3, L_ice, N_ice, F_rim, ρ_rim), 30_000)
    bench_press(FT, P3.ice_terminal_velocity_number_weighted, (ch2022, ρ_air, state, logλ), 120_000)
    bench_press(FT, P3.ice_terminal_velocity_mass_weighted, (ch2022, ρ_air, state, logλ), 135_000)
    bench_press(FT, P3.integrate, (x -> x^4, FT(0), FT(1)), 3_000)
    bench_press(FT, P3.D_m, (state, logλ), 3_000)

    @info "P3 Ice Nucleation"
    bench_press(
        @NamedTuple{dNdt::FT, dLdt::FT},
        P3.het_ice_nucleation,
        (kaolinite, tps, q_liq, N_liq, RH_2, T_air_2, ρ_air, Δt),
        200,
    )
    bench_press(
        @NamedTuple{dNdt::FT, dLdt::FT},
        P3.ice_melt,
        (ch2022, aps, tps, T_air, ρ_air, Δt, state, logλ),
        150_000,
    )
    bench_press(FT, CMI_het.P3_deposition_N_i, (ip.p3, T_air_cold), 230)
    bench_press(FT, CMI_het.P3_het_N_i, (ip.p3, T_air_cold, N_liq, V_liq, Δt), 230)

    @info "Cloud/Ice Terminal Velocity (Non-Eq)"
    bench_press(FT, CMN.terminal_velocity, (liquid, stokes_vel, ρ_air, q_liq), 350)
    bench_press(FT, CMN.terminal_velocity, (ice, ch2022.small_ice, ρ_air, q_ice), 400)

    @info "Precipitation Terminal Velocity (1M)"
    bench_press(FT, CM1.terminal_velocity, (rain, ch2022.rain, ρ_air, q_rai), 850)
    bench_press(FT, CM1.terminal_velocity, (snow, ch2022.large_ice, ρ_air, q_sno), 850)

    @info "Aerosol Activation"
    bench_press(FT, AA.total_N_activated, (ap, aer_distr, aps, tps, T_air, p_air, w_air, q_tot, FT(0), FT(0)), 1300)

    @info "Common Functions"
    bench_press(FT, CO.H2SO4_soln_saturation_vapor_pressure, (H2SO4_prs, x_sulph, T_air_cold), 50)
    bench_press(FT, CO.a_w_xT, (H2SO4_prs, tps, x_sulph, T_air_cold), 230)
    bench_press(FT, CO.a_w_eT, (tps, e, T_air_cold), 230)
    bench_press(FT, CO.a_w_ice, (tps, T_air_cold), 230)

    @info "Ice Nucleation"
    bench_press(FT, CMI_het.dust_activated_number_fraction, (desert_dust, ip.deposition, S_ice, T_air_2), 50)
    bench_press(FT, CMI_het.MohlerDepositionRate, (desert_dust, ip.deposition, S_ice, T_air_2, dSi_dt, N_aer), 80)
    bench_press(FT, CMI_het.deposition_J, (kaolinite, Delta_a_w), 230)
    bench_press(FT, CMI_het.ABIFM_J, (desert_dust, Delta_a_w), 230)
    bench_press(FT, CMI_het.INP_concentration_frequency, (ip_frostenberg, INPC, T_air_cold), 150)
    bench_press(FT, CMI_hom.homogeneous_J_cubic, (ip.homogeneous, Delta_a_w), 230)
    bench_press(FT, CMI_hom.homogeneous_J_linear, (ip.homogeneous, Delta_a_w), 230)

    @info "Non-equilibrium Microphysics"
    bench_press(FT, CMN.τ_relax, (liquid,), 15)
    bench_press(FT, CMN.conv_q_vap_to_q_lcl_icl, (ice, FT(2e-3), FT(1e-3)), 15)
    bench_press(
        FT,
        CMN.conv_q_vap_to_q_lcl_icl_MM2015,
        (liquid, tps, FT(0.00145), FT(0), FT(0), FT(0), FT(0), FT(0.8), FT(263)),
        70,
    )

    @info "0-Moment Scheme"
    bench_press(FT, CM0.remove_precipitation, (p0m, q_liq, q_ice), 12)

    @info "1-Moment Scheme"
    bench_press(FT, CM1.accretion, (liquid, rain, blk1mvel.rain, ce, q_liq, q_rai, ρ_air), 360)
    bench_press(FT, CMD.radar_reflectivity_1M, (rain, q_rai, ρ_air), 300)

    @info "2-Moment Scheme"
    for sb in [sb2006, sb2006_no_limiters]
        bench_press(
            @NamedTuple{au::CM2.LclRaiRates{FT}, sc::FT},
            CM2.autoconversion_and_cloud_liquid_self_collection,
            (sb, q_liq, q_rai, ρ_air, N_liq),
            300,
        )
        bench_press(@NamedTuple{sc::FT, br::FT}, CM2.rain_self_collection_and_breakup, (sb, q_rai, ρ_air, N_rai), 1200)
        bench_press(
            @NamedTuple{evap_rate_0::FT, evap_rate_1::FT},
            CM2.rain_evaporation,
            (sb, aps, tps, q_tot, q_liq, q_ice, q_rai, q_sno, ρ_air, N_rai, T_air),
            2000,
        )
        bench_press(Tuple{FT, FT}, CM2.rain_terminal_velocity, (sb, sb2006vel, q_rai, ρ_air, N_rai), 700)
        bench_press(Tuple{FT, FT}, CM2.rain_terminal_velocity, (sb, ch2022.rain, q_rai, ρ_air, N_rai), 2200)
        bench_press(FT, CMD.radar_reflectivity_2M, (sb, q_liq, q_rai, N_liq, N_rai, ρ_air), 2000)
        bench_press(FT, CMD.effective_radius_2M, (sb, q_liq, q_rai, N_liq, N_rai, ρ_air), 2000)

        @info "P3 Collisions"
        bench_press(@NamedTuple{∂ₜq_c::FT, ∂ₜq_r::FT, ∂ₜN_c::FT, ∂ₜN_r::FT, ∂ₜL_rim::FT, ∂ₜL_ice::FT, ∂ₜB_rim::FT},
            P3.bulk_liquid_ice_collision_sources,
            (
                params_P3, logλ, L_ice, N_ice, F_rim, ρ_rim,
                sb.pdf_c, sb.pdf_r, ρ_air * q_liq, N_liq, ρ_air * q_rai, N_rai,
                aps, tps, ch2022,
                ρ_air, T_air,
            ), 1e9)
    end
    bench_press(FT, CMD.effective_radius_Liu_Hallet_97, (wtr, ρ_air, q_liq, N_liq, q_rai, N_rai), 300)
    bench_press(FT, CM2.number_increase_for_mass_limit, (sb2006.numadj, FT(5e-6), q_rai, ρ_air, N_rai), 50)
    bench_press(FT, CM2.number_decrease_for_mass_limit, (sb2006.numadj, FT(2.6e-10), q_rai, ρ_air, N_rai), 50)

    @info "Homogeneous Nucleation"
    bench_press(
        @NamedTuple{binary_rate::FT, ternary_rate::FT},
        HN.h2so4_nucleation_rate,
        (FT(1e12), FT(1), FT(1), FT(208), h2so4_nuc),
        470,
    )
    bench_press(FT, HN.organic_nucleation_rate, (FT(0), FT(1e3), FT(1e3), FT(1e3), FT(300), FT(1), organ_nuc), 850)
    bench_press(FT, HN.organic_and_h2so4_nucleation_rate, (FT(2.6e6), FT(1), FT(1), FT(300), FT(1), mixed_nuc), 120)
end

TT.@testset "Performance Tests ($FT)" for FT in (Float64, Float32)
    benchmark_test(FT)
end
nothing
