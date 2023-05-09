import Test as TT
import BenchmarkTools as BT

import CloudMicrophysics as CM
import Thermodynamics as TD
import CLIMAParameters as CP

const CMT = CM.CommonTypes
const AM = CM.AerosolModel
const AA = CM.AerosolActivation
const CMI = CM.HetIceNucleation
const CMN = CM.MicrophysicsNonEq
const CM0 = CM.Microphysics0M
const CM1 = CM.Microphysics1M
const CM2 = CM.Microphysics2M

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

function bench_press(foo, args, min_run_time)

    println("Testing ", "$foo")
    # Calling foo once before benchmarking
    # to make sure compile time is not included in the benchmark
    foo(args...)
    # Benchmark foo
    trail = BT.@benchmark $foo($args...)
    show(stdout, MIME("text/plain"), trail)
    println("\n")

    TT.@test BT.minimum(trail).time < min_run_time
    TT.@test trail.memory == 0
    TT.@test trail.allocs == 0
end

function benchmark_test(FT)

    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    liquid = CMT.LiquidType()
    rain = CMT.RainType()
    sb2006 = CMT.SB2006Type()
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

    T_air_2 = FT(250)
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

    # aerosol activation
    bench_press(
        AA.total_N_activated,
        (prs, aer_distr, T_air, p_air, w_air, q),
        1000,
    )

    # ice nucleation
    bench_press(CMI.dust_activated_number_fraction, (S_ice, T_air_2, dust), 30)

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
        800,
    )
    bench_press(
        CM2.rain_evaporation,
        (prs, sb2006, q, q_rai, ρ_air, N_rai, T_air),
        2000,
    )

end

benchmark_test(Float64)
benchmark_test(Float32)
