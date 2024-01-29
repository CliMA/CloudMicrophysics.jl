import Test as TT
import BenchmarkTools as BT

import Thermodynamics as TD

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

@info "Performance Tests"

function bench_press(
    foo,
    args,
    min_run_time,
    min_memory = 0.0,
    min_allocs = 0.0,
)
    println("Testing ", "$foo")
    # Calling foo once before benchmarking
    # to make sure compile time is not included in the benchmark
    foo(args...)
    # Benchmark foo
    trail = BT.@benchmark $foo($args...)
    show(stdout, MIME("text/plain"), trail)
    println("\n")

    TT.@test BT.minimum(trail).time < min_run_time
    TT.@test trail.memory <= min_memory
    TT.@test trail.allocs <= min_allocs
end

function benchmark_test(FT)

    # 0-moment microphysics
    p0m = CMP.Parameters0M(FT)
    # 1-moment microphysics
    liquid = CMP.CloudLiquid(FT)
    rain = CMP.Rain(FT)
    ce = CMP.CollisionEff(FT)
    # 2-moment microphysics
    sb2006 = CMP.SB2006(FT)
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
    # aerosol nucleation parameters
    mixed_nuc = CMP.MixedNucleationParameters(FT)
    h2so4_nuc = CMP.H2S04NucleationParameters(FT)
    organ_nuc = CMP.OrganicNucleationParameters(FT)
    # air and thermodunamics parameters
    aps = CMP.AirProperties(FT)
    tps = TD.Parameters.ThermodynamicsParameters(FT)

    ρ_air = FT(1.2)
    T_air = FT(280)
    q_liq = FT(5e-4)
    q_ice = FT(5e-4)
    q_tot = FT(1e-3)
    q = TD.PhasePartition(q_tot)
    q_rai = FT(1e-4)
    q_sno = FT(1e-4)
    N_liq = FT(1e8)
    N_rai = FT(1e8)
    e = FT(600)

    ρ_r = FT(400.0)
    F_r = FT(0.95)

    T_air_2 = FT(250)
    T_air_cold = FT(230)
    S_ice = FT(1.2)

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
        1,
    )
    aer_distr = AM.AerosolDistribution((seasalt_mode,))

    x_sulph = FT(0.1)
    Delta_a_w = FT(0.27)

    # P3 scheme
    bench_press(P3.thresholds, (p3, ρ_r, F_r), 12e6, 2048, 80)

    # aerosol activation
    bench_press(
        AA.total_N_activated,
        (ap, aer_distr, aps, tps, T_air, p_air, w_air, q),
        1300,
    )

    # Common
    bench_press(
        CO.H2SO4_soln_saturation_vapor_pressure,
        (H2SO4_prs, x_sulph, T_air_cold),
        50,
    )
    bench_press(CO.a_w_xT, (H2SO4_prs, tps, x_sulph, T_air_cold), 230)
    bench_press(CO.a_w_eT, (tps, e, T_air_cold), 230)
    bench_press(CO.a_w_ice, (tps, T_air_cold), 230)

    # ice nucleation
    bench_press(
        CMI_het.dust_activated_number_fraction,
        (desert_dust, ip.deposition, S_ice, T_air_2),
        50,
    )
    bench_press(CMI_het.deposition_J, (kaolinite, Delta_a_w), 230)
    bench_press(CMI_het.ABIFM_J, (desert_dust, Delta_a_w), 230)
    bench_press(CMI_hom.homogeneous_J, (ip.homogeneous, Delta_a_w), 230)

    # non-equilibrium
    bench_press(CMN.τ_relax, (liquid,), 10)

    # 0-moment
    bench_press(CM0.remove_precipitation, (p0m, q), 10)

    # 1-moment
    bench_press(
        CM1.accretion,
        (liquid, rain, blk1mvel.rain, ce, q_liq, q_rai, ρ_air),
        350,
    )

    # 2-moment
    bench_press(
        CM2.autoconversion_and_liquid_self_collection,
        (sb2006, q_liq, q_rai, ρ_air, N_liq),
        250,
    )
    bench_press(
        CM2.rain_self_collection_and_breakup,
        (sb2006, q_rai, ρ_air, N_rai),
        820,
    )
    bench_press(
        CM2.rain_evaporation,
        (sb2006, aps, tps, q, q_rai, ρ_air, N_rai, T_air),
        2000,
    )
    bench_press(
        CM2.rain_terminal_velocity,
        (sb2006, sb2006vel, q_rai, ρ_air, N_rai),
        300,
    )
    bench_press(
        CM2.rain_terminal_velocity,
        (sb2006, ch2022.rain, q_rai, ρ_air, N_rai),
        1700,
    )
    # Homogeneous Nucleation
    bench_press(HN.h2so4_nucleation_rate, (1e12, 1.0, 1.0, 208, h2so4_nuc), 470)
    bench_press(
        HN.organic_nucleation_rate,
        (0.0, 1e3, 1e3, 1e3, 300, 1, organ_nuc),
        650,
    )
    bench_press(
        HN.organic_and_h2so4_nucleation_rate,
        (2.6e6, 1.0, 1.0, 300, 1, mixed_nuc),
        120,
    )
end

println("Testing Float64")
benchmark_test(Float64)

println("Testing Float32")
benchmark_test(Float32)
