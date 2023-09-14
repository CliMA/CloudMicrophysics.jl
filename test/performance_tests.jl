import Test as TT
import BenchmarkTools as BT

import CloudMicrophysics as CM
import Thermodynamics as TD
import CLIMAParameters as CP

import CloudMicrophysics.CommonTypes as CMT
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
import CloudMicrophysics.Parameters.MixedNucleationParameters
import CloudMicrophysics.Parameters.H2S04NucleationParameters
import CloudMicrophysics.Parameters.OrganicNucleationParameters
import CloudMicrophysics.P3Scheme as P3

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

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

    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)
    p3 = CMP.CloudMicrophysicsParametersP3(FT)
    liquid = CMT.LiquidType()
    rain = CMT.RainType()
    sb2006 = CMT.SB2006Type()
    sb2006vel = CMT.SB2006VelType()
    ch2022 = CMT.Chen2022Type()
    dust = CMT.DesertDustType()

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
        (prs, aer_distr, T_air, p_air, w_air, q),
        1300,
    )

    # Common
    bench_press(
        CO.H2SO4_soln_saturation_vapor_pressure,
        (prs, x_sulph, T_air_cold),
        50,
    )
    bench_press(CO.Delta_a_w, (prs, x_sulph, T_air_cold), 230)

    # ice nucleation
    bench_press(
        CMI_het.dust_activated_number_fraction,
        (prs, S_ice, T_air_2, dust),
        50,
    )
    bench_press(CMI_het.ABIFM_J, (prs, dust, Delta_a_w), 230)
    bench_press(CMI_hom.homogeneous_J, (prs, Delta_a_w), 230)

    # non-equilibrium
    bench_press(CMN.τ_relax, (prs, liquid), 10)

    # 0-moment
    bench_press(CM0.remove_precipitation, (prs, q), 10)

    # 1-moment
    bench_press(CM1.accretion, (prs, liquid, rain, q_liq, q_rai, ρ_air), 350)

    # 2-moment
    bench_press(
        CM2.autoconversion_and_liquid_self_collection,
        (prs, sb2006, q_liq, q_rai, ρ_air, N_liq),
        250,
    )
    bench_press(
        CM2.rain_self_collection_and_breakup,
        (prs, sb2006, q_rai, ρ_air, N_rai),
        820,
    )
    bench_press(
        CM2.rain_evaporation,
        (prs, sb2006, q, q_rai, ρ_air, N_rai, T_air),
        2000,
    )
    bench_press(
        CM2.rain_terminal_velocity,
        (prs, sb2006, sb2006vel, q_rai, ρ_air, N_rai),
        300,
    )
    bench_press(
        CM2.rain_terminal_velocity,
        (prs, sb2006, ch2022, q_rai, ρ_air, N_rai),
        1700,
    )
    # Homogeneous Nucleation
    bench_press(
        HN.h2so4_nucleation_rate,
        (1e12, 1.0, 1.0, 208, H2S04NucleationParameters(FT)),
        470,
    )
    bench_press(
        HN.organic_nucleation_rate,
        (0.0, 1e3, 1e3, 1e3, 300, 1, OrganicNucleationParameters(FT)),
        650,
    )
    bench_press(
        HN.organic_and_h2so4_nucleation_rate,
        (2.6e6, 1.0, 1.0, 300, 1, MixedNucleationParameters(FT)),
        120,
    )
end

println("Testing Float64")
benchmark_test(Float64)

println("Testing Float32")
benchmark_test(Float32)
