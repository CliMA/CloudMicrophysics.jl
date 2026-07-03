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
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.Nucleation as HN
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.CloudDiagnostics as CMD
import Profile

# Timing bounds below are scaled by a measured machine-speed factor.

# Fixed scalar workload, dominated by transcendental functions, for measuring
# the machine-speed factor.
function calibration_workload(x)
    FT = typeof(x)
    s = zero(FT)
    for _ in 1:256
        x = FT(0.5) + FT(0.25) * (one(FT) + sin(x))
        s += exp(-x) + log1p(x) + sqrt(x) + cbrt(x) + x^FT(1.5)
    end
    return s
end

# Minimum runtime of `calibration_workload(1.0)` on the reference machine, in nanoseconds.
const REFERENCE_CALIBRATION_NS = 12_200.0

# Upper limit on the machine-speed factor applied to timing bounds.
const MAX_SPEED_FACTOR = 20.0

# Machine-speed tolerance multiplier for the upper bound, measured once per suite.
const SPEED_FACTOR = Ref(1.0)

# Machine-speed multiplier for the lower bound, measured once per suite.
# The upper factor is clamped to [1, MAX_SPEED_FACTOR] so slower machines gain
# headroom. The lower bound instead tracks machines faster than the reference,
# where the raw factor drops below one; clamping to
# [1 / MAX_SPEED_FACTOR, MAX_SPEED_FACTOR] keeps it proportional to measured time
# so it does not trip on a fast machine.
const SPEED_FACTOR_LOWER = Ref(1.0)

function measure_speed_factor()
    trial = BT.@benchmark calibration_workload($(1.0)) samples = 100 evals = 100
    raw = BT.minimum(trial).time / REFERENCE_CALIBRATION_NS
    return clamp(raw, 1.0, MAX_SPEED_FACTOR)
end

function measure_speed_factor_lower()
    trial = BT.@benchmark calibration_workload($(1.0)) samples = 100 evals = 100
    raw = BT.minimum(trial).time / REFERENCE_CALIBRATION_NS
    return clamp(raw, 1.0 / MAX_SPEED_FACTOR, MAX_SPEED_FACTOR)
end

# Timing is two-sided: the upper bound catches regressions; the lower bound
# (min_expected_time ≈ one third of the fastest measured time, 0 for
# jitter-dominated sub-50 ns calls) flags a genuine speedup so the bounds get
# retuned.
function bench_press(
    type,
    foo,
    args,
    min_run_time,  # upper bound on run time, in nanoseconds
    min_memory = 0.0,
    min_allocs = 0.0;
    min_expected_time = 0.0,  # lower bound on run time, in nanoseconds; 0 disables
)
    @info "Testing $foo"

    # converts a string argument to tuple
    args isa String && (args = (args,))

    # Benchmark foo
    # Note: `$(splat(foo))($args)` is equivalent to `$foo($args...)`
    # but the latter might allocate during the benchmark call due to `...`
    trail = BT.@benchmark $(splat(foo))($args) samples = 100 evals = 100

    show(stdout, MIME("text/plain"), trail)
    println("\n")

    TT.@test min_expected_time * SPEED_FACTOR_LOWER[] < BT.minimum(trail).time < min_run_time * SPEED_FACTOR[]
    TT.@test trail.memory <= min_memory

    if !(trail.allocs <= min_allocs)
        # If allocations are above the threshold, print the allocations
        Profile.clear()
        Profile.Allocs.@profile sample_rate = 1 foo(args...)
        results = Profile.Allocs.fetch()
        sorted = sort(results.allocs, by = x -> x.size)
        if isempty(sorted)
            @info "No allocations detected using Profile.Allocs.@profile"
        else
            sorted_str = join(sorted, "\n\n")
            @error sorted_str
            # largest allocation
            trace = sorted[end].stacktrace
            @error join(trace, "\n")
        end
        Profile.clear()
        TT.@test trail.allocs <= min_allocs
    end

    # Test that foo is free from optimization failures
    # and unresolved method dispatches
    JET.@test_opt(foo(args...))

    # Test that the return type of foo is of correct type (Float32 vs Float64)
    TT.@test typeof(foo(args...)) <: type
end

function benchmark_test(FT)

    @info "Artifacts"
    bench_press(String, AFC.AIDA_ice_nucleation, ("in05_17_aida.edf"), 14_000, 30_000, 300; min_expected_time = 2_900)

    # 0-moment parameters
    p0m = CMP.Parameters0M(FT)

    # 1-moment parameters
    liquid = CMP.CloudLiquid(FT)
    ice = CMP.CloudIce(FT)
    rain = CMP.Rain(FT)
    snow = CMP.Snow(FT)

    # 2-moment parameters
    override_file = joinpath(pkgdir(CM), "src", "parameters", "toml", "SB2006_limiters.toml")
    toml_dict = CP.create_toml_dict(FT; override_file)
    sb2006 = CMP.SB2006(toml_dict)
    sb2006_no_limiters = CMP.SB2006(toml_dict; is_limited = false)

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

    # 1-moment microphysics unified params (both liquid autoconversion flavours)
    mp_1m = CMP.Microphysics1MParams(FT)
    mp_1m_2M = CMP.Microphysics1MParams(FT;
        rain_autoconversion = CMP.PrescribedNd(CP.create_toml_dict(FT)),
    )
    E_lcl_rai = mp_1m.options.cloud_liquid_rain_accretion.e

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
    INPC = FT(1e5)

    @info "P3 Scheme"
    state = P3.P3State(params_P3, L_ice, N_ice, F_rim, ρ_rim)
    logλ = P3.get_distribution_logλ(state)
    bench_press(
        P3.P3State,
        P3.P3State,
        (params_P3, L_ice, N_ice, F_rim, ρ_rim),
        200,
        min_expected_time = 35,
    )
    bench_press(FT, P3.get_distribution_logλ, (state,), 35_000; min_expected_time = 3_500)  # 10 (F64) / 8 (F32) FixedIterations BrentsMethod, zero warp divergence
    # The weighted-velocity integrals build a nested terminal-velocity closure
    # that escapes into `integrate`. With `ice_particle_terminal_velocity`
    # returning a single (concretely-typed) closure, this path is type-stable
    # on both 1.10 and 1.12 and allocates nothing — keep the default zero
    # allocation/memory budget so a future closure regression is caught here.
    bench_press(
        FT,
        P3.ice_terminal_velocity_number_weighted,
        (ch2022, ρ_air, state, logλ),
        30_000;
        min_expected_time = 4_700,
    )
    bench_press(
        FT,
        P3.ice_terminal_velocity_mass_weighted,
        (ch2022, ρ_air, state, logλ),
        35_000;
        min_expected_time = 4_800,
    )
    bench_press(FT, P3.integrate, (x -> x^4, FT(0), FT(1)), 4_500; min_expected_time = 400)
    bench_press(FT, P3.D_m, (state, logλ), 3_000; min_expected_time = 310)

    @info "P3 Ice Nucleation"
    bench_press(
        @NamedTuple{dNdt::FT, dLdt::FT},
        P3.het_ice_nucleation,
        (kaolinite, tps, q_liq, N_liq, RH_2, T_air_2, ρ_air),
        140,
        min_expected_time = 21,
    )
    bench_press(
        @NamedTuple{dNdt::FT, dLdt::FT},
        P3.ice_melt,
        (ch2022, aps, tps, T_air, ρ_air, state, logλ),
        35_000,
        min_expected_time = 5_100,
    )
    bench_press(FT, CMI_het.P3_deposition_N_i, (ip.p3, T_air_cold), 18)
    bench_press(FT, CMI_het.P3_het_N_i, (ip.p3, T_air_cold, N_liq, V_liq, FT(25)), 30)

    @info "Cloud/Ice Terminal Velocity (Non-Eq)"
    bench_press(FT, CMN.terminal_velocity, (liquid, stokes_vel, ρ_air, q_liq), 35)
    bench_press(FT, CMN.terminal_velocity, (ice, ch2022.small_ice, ρ_air, q_ice), 300; min_expected_time = 51)

    @info "Precipitation Terminal Velocity (1M)"
    bench_press(FT, CM1.terminal_velocity, (rain, ch2022.rain, ρ_air, q_rai), 600; min_expected_time = 110)
    bench_press(FT, CM1.terminal_velocity, (snow, ch2022.large_ice, ρ_air, q_sno), 600; min_expected_time = 100)

    @info "Aerosol Activation"
    bench_press(
        FT,
        AA.total_N_activated,
        (ap, aer_distr, aps, tps, T_air, p_air, w_air, q_tot, FT(0), FT(0)),
        500;
        min_expected_time = 94,
    )

    @info "Common Functions"
    bench_press(FT, CO.H2SO4_soln_saturation_vapor_pressure, (H2SO4_prs, x_sulph, T_air_cold), 20)
    bench_press(FT, CO.a_w_xT, (H2SO4_prs, tps, x_sulph, T_air_cold), 70)
    bench_press(FT, CO.a_w_eT, (tps, e, T_air_cold), 50)
    bench_press(FT, CO.a_w_ice, (tps, T_air_cold), 100)

    @info "Ice Nucleation"
    bench_press(FT, CMI_het.dust_activated_number_fraction, (desert_dust, ip.deposition, S_ice, T_air_2), 20)
    bench_press(FT, CMI_het.MohlerDepositionRate, (desert_dust, ip.deposition, S_ice, T_air_2, dSi_dt, N_aer), 14)
    bench_press(FT, CMI_het.deposition_J, (kaolinite, Delta_a_w), 50)
    bench_press(FT, CMI_het.ABIFM_J, (desert_dust, Delta_a_w), 50)
    bench_press(FT, CMI_het.INP_concentration_frequency, (ip_frostenberg, INPC, T_air_cold), 70)
    bench_press(FT, CMI_hom.homogeneous_J_cubic, (ip.homogeneous, Delta_a_w), 60)
    bench_press(FT, CMI_hom.homogeneous_J_linear, (ip.homogeneous, Delta_a_w), 50)

    @info "Non-equilibrium Microphysics"
    bench_press(FT, CMN.τ_relax, (ice, aps, ip_frostenberg, FT(1e-4), FT(250)), 90)

    mp_mock = (; cloud = (; liquid = liquid))
    micro_mock = (; q_tot = FT(0.00145), q_lcl = FT(0), q_icl = FT(0), q_rai = FT(0), q_sno = FT(0))
    thermo_mock = (; ρ = FT(0.8), T = FT(263))
    bench_press(FT, CMN.conv_q_vap_to_q_lcl,
        (CMP.CloudLiquidFormation(CP.create_toml_dict(FT)), mp_mock, tps, micro_mock, thermo_mock), 70)

    @info "0-Moment Scheme"
    bench_press(FT, CM0.remove_precipitation, (p0m, q_liq, q_ice), 14)
    bench_press(FT, CM0.∂remove_precipitation_∂q_tot, (p0m, q_liq, q_ice), 14)

    @info "1-Moment Scheme"
    micro_1m = (; q_tot, q_lcl = q_liq, q_icl = q_ice, q_rai, q_sno)
    thermo_1m = (; ρ = ρ_air, T = T_air)
    bench_press(
        FT, CM1.conv_q_lcl_to_q_rai,
        (mp_1m.options.rain_autoconversion, mp_1m, tps, micro_1m, thermo_1m),
        100,
    )
    bench_press(
        FT, CM1.conv_q_lcl_to_q_rai,
        (mp_1m_2M.options.rain_autoconversion, mp_1m_2M, tps, micro_1m, thermo_1m),
        16,
    )
    bench_press(
        FT,
        CM1.accretion,
        (liquid, rain, blk1mvel.rain, E_lcl_rai, q_liq, q_rai, ρ_air),
        140;
        min_expected_time = 21,
    )
    bench_press(
        FT, CM1.accretion, (mp_1m.options.cloud_liquid_rain_accretion, mp_1m, tps, micro_1m, thermo_1m), 300;
        min_expected_time = 53,
    )
    bench_press(
        @NamedTuple{S_accr::FT, S_melt::FT},
        CM1.accretion,
        (mp_1m.options.cloud_liquid_snow_accretion, mp_1m, tps, micro_1m, thermo_1m),
        300,
        min_expected_time = 51,
    )
    bench_press(
        FT, CM1.accretion, (mp_1m.options.cloud_ice_rain_accretion, mp_1m, tps, micro_1m, thermo_1m), 300;
        min_expected_time = 53,
    )
    bench_press(
        FT, CM1.accretion, (mp_1m.options.cloud_ice_snow_accretion, mp_1m, tps, micro_1m, thermo_1m), 300;
        min_expected_time = 53,
    )
    bench_press(
        @NamedTuple{S_rai_sno::FT, S_sno_rai::FT, S_melt::FT},
        CM1.accretion_snow_rain,
        (mp_1m.options.rain_snow_accretion, mp_1m, tps, micro_1m, thermo_1m),
        600,
        min_expected_time = 100,
    )
    bench_press(FT, CMD.radar_reflectivity_1M, (rain, q_rai, ρ_air), 120)

    @info "1-Moment Bulk Tendencies"
    CM1M = BMT.Microphysics1Moment()
    Δt = FT(1)
    bench_press(
        @NamedTuple{dq_lcl_dt::FT, dq_icl_dt::FT, dq_rai_dt::FT, dq_sno_dt::FT},
        BMT.bulk_microphysics_tendencies,
        (BMT.Instantaneous(), CM1M, mp_1m, tps, ρ_air, T_air, q_tot, q_liq, q_ice, q_rai, q_sno),
        1500,
        min_expected_time = 230,
    )
    bench_press(
        @NamedTuple{dq_lcl_dt::FT, dq_icl_dt::FT, dq_rai_dt::FT, dq_sno_dt::FT},
        BMT.bulk_microphysics_tendencies,
        (BMT.LinearizedAverage(), CM1M, mp_1m, tps, ρ_air, T_air, q_tot, q_liq, q_ice, q_rai, q_sno, Δt),
        1500,
        min_expected_time = 230,
    )
    bench_press(
        @NamedTuple{dq_lcl_dt::FT, dq_icl_dt::FT, dq_rai_dt::FT, dq_sno_dt::FT},
        BMT.bulk_microphysics_tendencies,
        (BMT.LinearizedAverage(), CM1M, mp_1m, tps, ρ_air, T_air, q_tot, q_liq, q_ice, q_rai, q_sno, Δt, 3),
        4000,
        min_expected_time = 710,
    )



    @info "2-Moment Scheme"
    for sb in [sb2006, sb2006_no_limiters]
        bench_press(
            @NamedTuple{au::CM2.LclRaiRates{FT}, sc::FT},
            CM2.autoconversion_and_cloud_liquid_self_collection,
            (sb, q_liq, q_rai, ρ_air, N_liq),
            140,
            min_expected_time = 22,
        )
        bench_press(
            @NamedTuple{sc::FT, br::FT}, CM2.rain_self_collection_and_breakup, (sb, q_rai, ρ_air, N_rai), 180;
            min_expected_time = 18,
        )
        bench_press(
            @NamedTuple{∂ₜρn_rai::FT, ∂ₜq_rai::FT},
            CM2.rain_evaporation,
            (sb, aps, tps, q_tot, q_liq, q_ice, q_rai, q_sno, ρ_air, N_rai, T_air),
            400,
            min_expected_time = 66,
        )
        bench_press(Tuple{FT, FT}, CM2.rain_terminal_velocity, (sb, sb2006vel, q_rai, ρ_air, N_rai), 120)
        bench_press(
            Tuple{FT, FT},
            CM2.rain_terminal_velocity,
            (sb, ch2022.rain, q_rai, ρ_air, N_rai),
            800;
            min_expected_time = 150,
        )
        bench_press(FT, CMD.radar_reflectivity_2M, (sb, q_liq, q_rai, N_liq, N_rai, ρ_air), 400; min_expected_time = 68)
        bench_press(FT, CMD.effective_radius_2M, (sb, q_liq, q_rai, N_liq, N_rai, ρ_air), 600; min_expected_time = 95)

        @info "P3 Collisions"
        # Julia <= 1.11 inference exceeds its depth budget on this collision
        # assembly, widening intermediates to Any (runtime dispatch + boxing;
        # correct but unoptimized), so the JET and allocation assertions fail
        # spuriously there. 1.12 resolves the chain: zero JET reports, zero
        # allocations. TODO: drop the gate once CI runs >= 1.12.
        if VERSION >= v"1.12"
            bench_press(@NamedTuple{∂ₜq_c::FT, ∂ₜq_r::FT, ∂ₜN_c::FT, ∂ₜN_r::FT, ∂ₜL_rim::FT, ∂ₜL_ice::FT, ∂ₜB_rim::FT},
                P3.bulk_liquid_ice_collision_sources,
                (
                    state, logλ,
                    sb.pdf_c, sb.pdf_r, ρ_air * q_liq, N_liq, ρ_air * q_rai, N_rai,
                    aps, tps, ch2022,
                    ρ_air, T_air,
                ), 17_000_000; min_expected_time = 2_700_000)
        end
    end
    bench_press(FT, CMD.effective_radius_Liu_Hallet_97, (wtr, ρ_air, q_liq, N_liq, q_rai, N_rai), 50)
    bench_press(FT, CM2.number_increase_for_mass_limit, (sb2006.numadj, FT(5e-6), q_rai, ρ_air, N_rai), 16)
    bench_press(FT, CM2.number_decrease_for_mass_limit, (sb2006.numadj, FT(2.6e-10), q_rai, ρ_air, N_rai), 14)

    @info "Homogeneous Nucleation"
    bench_press(
        @NamedTuple{binary_rate::FT, ternary_rate::FT},
        HN.h2so4_nucleation_rate,
        (FT(1e12), FT(1), FT(1), FT(208), h2so4_nuc),
        250,
        min_expected_time = 37,
    )
    bench_press(
        FT,
        HN.organic_nucleation_rate,
        (FT(0), FT(1e3), FT(1e3), FT(1e3), FT(300), FT(1), organ_nuc),
        160;
        min_expected_time = 25,
    )
    bench_press(FT, HN.organic_and_h2so4_nucleation_rate, (FT(2.6e6), FT(1), FT(1), FT(300), FT(1), mixed_nuc), 25)
end

SPEED_FACTOR[] = measure_speed_factor()
SPEED_FACTOR_LOWER[] = measure_speed_factor_lower()
@info "Machine-speed factor: $(SPEED_FACTOR[]) (lower bound: $(SPEED_FACTOR_LOWER[]))"
TT.@testset "Performance Tests ($FT)" for FT in (Float64, Float32)
    benchmark_test(FT)
end
nothing
