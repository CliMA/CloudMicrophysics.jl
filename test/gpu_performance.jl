import Test as TT
using KernelAbstractions
using ClimaComms
ClimaComms.@import_required_backends
import BenchmarkTools as BT
import SpecialFunctions as SF

# Needed for parameters
import ClimaParams as CP

# Modules to test
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI
const TD = TDI.TD
import Random
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.Utilities as UT

work_groups = 2

if ClimaComms.device() isa ClimaComms.CUDADevice
    using CUDA
    backend = CUDABackend()
    CUDA.allowscalar(false)
    ArrayType = CuArray
    @info "GPU Performance Tests running on CUDA GPU" backend ArrayType
else
    backend = CPU()
    ArrayType = Array
    @info "No CUDA GPU found. GPU Performance Tests running on CPU fallback via KernelAbstractions" backend ArrayType
end

@kernel inbounds = true function benchmark_1m_bulk_tendencies_kernel!(
    mp, tps, output, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    i = @index(Global, Linear)
    output[i] = BMT.bulk_microphysics_tendencies(
        BMT.Instantaneous(), BMT.Microphysics1Moment(), mp, tps,
        ρ[i], T[i], q_tot[i], q_lcl[i], q_icl[i], q_rai[i], q_sno[i],
    )
end

@kernel inbounds = true function benchmark_2m_bulk_tendencies_kernel!(
    mp, tps, output, ρ, T, q_tot, q_lcl, q_rai, n_lcl, n_rai,
)
    i = @index(Global, Linear)
    output[i] = BMT.bulk_microphysics_tendencies(
        BMT.Microphysics2Moment(), mp, tps,
        ρ[i], T[i], q_tot[i], q_lcl[i], n_lcl[i], q_rai[i], n_rai[i],
    )
end

@kernel inbounds = true function benchmark_p3_kernel!(
    p3_params, vel_params, output, L_ice, N_ice, F_rim, ρ_rim, ρₐ,
)
    i = @index(Global, Linear)
    state = P3.P3State(p3_params, L_ice[i], N_ice[i], F_rim[i], ρ_rim[i])
    logλ = P3.get_distribution_logλ(state)
    sc = P3.ice_self_collection(state, logλ, vel_params, ρₐ[i])
    output[i] = (; logλ, sc)
end

function setup_output(dims, FT, default_value = nothing)
    output = allocate(backend, FT, dims...)
    !isnothing(default_value) && fill!(output, FT(default_value))
    return (; output, ndrange = (dims[end],))
end

function constant_data(value; ndrange)
    FT = eltype(value)
    return fill!(allocate(backend, FT, ndrange), value)
end

function generate_atmospheric_states(N_states, FT, tps)
    # Construct realistic atmospheric states using DecayingTemperatureProfile
    z_arr = FT.(range(0, 15000, length = N_states))
    profile = TD.TemperatureProfiles.DecayingTemperatureProfile{FT}(tps, FT(300), FT(215))
    RH_arr = FT.(range(0.05, 1.05, length = N_states))
    rng = Random.MersenneTwister(1234)

    T_vec = Vector{FT}(undef, N_states)
    ρ_vec = Vector{FT}(undef, N_states)
    q_tot_vec = Vector{FT}(undef, N_states)
    q_lcl_vec = Vector{FT}(undef, N_states)
    q_icl_vec = Vector{FT}(undef, N_states)
    q_rai_vec = Vector{FT}(undef, N_states)
    q_sno_vec = Vector{FT}(undef, N_states)

    for i in 1:N_states
        T, p = profile(tps, z_arr[i])
        T_vec[i] = T
        RH = RH_arr[i]

        q_vap = TD.q_vap_from_RH(tps, p, T, RH, TD.Liquid())
        ρ = TD.air_density(tps, T, p, q_vap, zero(FT), zero(FT))
        ρ_vec[i] = ρ

        q_sat_liq = TD.q_vap_saturation(tps, T, ρ, TD.Liquid())
        q_sat_ice = TD.q_vap_saturation(tps, T, ρ, TD.Ice())

        # Construct liquid and ice condensates from supersaturation + random noise
        excess_liq = max(zero(FT), q_vap - q_sat_liq)
        excess_ice = max(zero(FT), q_vap - q_sat_ice)

        noise_lcl = FT(rand(rng, FT) * 1e-4)
        noise_icl = FT(rand(rng, FT) * 1e-4)
        noise_rai = FT(rand(rng, FT) * 1e-4)
        noise_sno = FT(rand(rng, FT) * 1e-4)

        q_lcl = excess_liq + noise_lcl
        q_icl = excess_ice + noise_icl

        q_lcl_vec[i] = q_lcl
        q_icl_vec[i] = q_icl
        q_rai_vec[i] = noise_rai
        q_sno_vec[i] = noise_sno

        q_tot_vec[i] = q_vap + q_lcl + q_icl + noise_rai + noise_sno
    end

    return (
        T = ArrayType(T_vec),
        ρ = ArrayType(ρ_vec),
        q_tot = ArrayType(q_tot_vec),
        q_lcl = ArrayType(q_lcl_vec),
        q_icl = ArrayType(q_icl_vec),
        q_rai = ArrayType(q_rai_vec),
        q_sno = ArrayType(q_sno_vec),
    )
end

function run_gpu_performance_benchmarks(FT)
    @info "--- Benchmarking with FT = $FT ---"

    # 1. Quality and Performance of Approximations (SF vs UT)
    TT.@testset "SpecialFunctions vs UT Approximations ($FT)" begin
        @info "Evaluating quality of gamma approximations..."
        for a in [FT(1.0), FT(2.0), FT(3.5)]
            for x in [FT(0.1), FT(1.0), FT(2.5), FT(5.0), FT(10.0)]
                P_sf, Q_sf = SF.gamma_inc(a, x)
                P_ut, Q_ut = UT.gamma_inc(a, x)
                TT.@test isapprox(P_sf, P_ut, atol = 1e-6)
                TT.@test isapprox(Q_sf, Q_ut, atol = 1e-6)
            end
        end

        for a in [FT(1.0), FT(2.0), FT(3.5)]
            for p in [FT(0.01), FT(0.1), FT(0.5), FT(0.9), FT(0.99)]
                q = 1 - p
                x_sf = SF.gamma_inc_inv(a, p, q)
                x_ut = UT.gamma_inc_inv(a, p, q)
                TT.@test isapprox(x_sf, x_ut, atol = 1e-5)
            end
        end

        a = FT(2.5);
        x = FT(3.0);
        p = FT(0.6);
        q = FT(0.4)
        b_sf = BT.@benchmark SF.gamma_inc($a, $x)
        b_ut = BT.@benchmark UT.gamma_inc($a, $x)
        @info "gamma_inc benchmark (SF vs UT):" SF=BT.minimum(b_sf) UT=BT.minimum(b_ut)

        b_inv_sf = BT.@benchmark SF.gamma_inc_inv($a, $p, $q)
        b_inv_ut = BT.@benchmark UT.gamma_inc_inv($a, $p, $q)
        @info "gamma_inc_inv benchmark (SF vs UT):" SF=BT.minimum(b_inv_sf) UT=BT.minimum(b_inv_ut)
    end

    # 2. Kernel Benchmarks (1-Moment BMT, 2-Moment BMT, P3)
    TT.@testset "Kernel Performance & Compile Time ($FT)" begin
        tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
        N_states = 1024
        states = generate_atmospheric_states(N_states, FT, tps)

        # 1-Moment Bulk Tendencies setup (realistic atmospheric states)
        mp_bulk = CMP.Microphysics1MParams(FT)
        DT_bulk = @NamedTuple{dq_lcl_dt::FT, dq_icl_dt::FT, dq_rai_dt::FT, dq_sno_dt::FT}
        (; output, ndrange) = setup_output(N_states, DT_bulk)

        T_arr = states.T
        ρ_arr = states.ρ
        q_tot_arr = states.q_tot
        q_lcl_arr = states.q_lcl
        q_icl_arr = states.q_icl
        q_rai_arr = states.q_rai
        q_sno_arr = states.q_sno

        kernel_1m_bulk! = benchmark_1m_bulk_tendencies_kernel!(backend, work_groups)
        t_compile_1m_bulk = @elapsed kernel_1m_bulk!(
            mp_bulk, tps, output, ρ_arr, T_arr, q_tot_arr, q_lcl_arr, q_icl_arr, q_rai_arr, q_sno_arr; ndrange,
        )
        KernelAbstractions.synchronize(backend)
        @info "1-Moment Bulk Tendencies Kernel first call (compile + run time): $(round(t_compile_1m_bulk, digits=4)) seconds"

        b_1m_bulk = BT.@benchmark (
            $kernel_1m_bulk!(
                $mp_bulk, $tps, $output, $ρ_arr, $T_arr, $q_tot_arr, $q_lcl_arr, $q_icl_arr, $q_rai_arr, $q_sno_arr;
                ndrange = $ndrange,
            );
            KernelAbstractions.synchronize($backend)
        )
        @info "1-Moment Bulk Tendencies Kernel runtime benchmark:" BT.minimum(b_1m_bulk)

        # 2-Moment Bulk Tendencies setup (warm rain only, realistic atmospheric states)
        mp_2m = CMP.Microphysics2MParams(FT)
        DT_2m_bulk = @NamedTuple{
            dq_lcl_dt::FT, dn_lcl_dt::FT, dq_rai_dt::FT, dn_rai_dt::FT,
            dq_ice_dt::FT, dq_rim_dt::FT, db_rim_dt::FT, dn_lcl_activation_dt::FT,
        }
        (; output, ndrange) = setup_output(N_states, DT_2m_bulk)

        # Realistic number concentrations for cloud and rain
        n_lcl_arr = constant_data(FT(1e8); ndrange)   # typical cloud droplet number
        n_rai_arr = constant_data(FT(1e4); ndrange)   # typical raindrop number

        kernel_2m_bulk! = benchmark_2m_bulk_tendencies_kernel!(backend, work_groups)
        t_compile_2m_bulk = @elapsed kernel_2m_bulk!(
            mp_2m, tps, output, ρ_arr, T_arr, q_tot_arr, q_lcl_arr, q_rai_arr, n_lcl_arr, n_rai_arr; ndrange,
        )
        KernelAbstractions.synchronize(backend)
        @info "2-Moment Bulk Tendencies Kernel first call (compile + run time): $(round(t_compile_2m_bulk, digits=4)) seconds"

        b_2m_bulk = BT.@benchmark (
            $kernel_2m_bulk!(
                $mp_2m, $tps, $output, $ρ_arr, $T_arr, $q_tot_arr, $q_lcl_arr, $q_rai_arr, $n_lcl_arr, $n_rai_arr;
                ndrange = $ndrange,
            );
            KernelAbstractions.synchronize($backend)
        )
        @info "2-Moment Bulk Tendencies Kernel runtime benchmark:" BT.minimum(b_2m_bulk)

        # P3 setup
        p3_params = CMP.ParametersP3(FT)
        Ch2022 = CMP.Chen2022VelType(FT)

        DT3 = @NamedTuple{logλ::FT, sc::@NamedTuple{dNdt::FT}}
        (; output, ndrange) = setup_output(N_states, DT3)
        L_ice = ArrayType(states.ρ .* states.q_icl)

        N_ice = constant_data(FT(1e8 * 1.2); ndrange)
        F_rim = constant_data(FT(0.95); ndrange)
        ρ_rim = constant_data(FT(400.0); ndrange)

        kernel_p3! = benchmark_p3_kernel!(backend, work_groups)
        t_compile_p3 = @elapsed kernel_p3!(p3_params, Ch2022, output, L_ice, N_ice, F_rim, ρ_rim, ρ_arr; ndrange)
        KernelAbstractions.synchronize(backend)
        @info "P3 Kernel first call (compile + run time): $(round(t_compile_p3, digits=4)) seconds"

        b_p3 = BT.@benchmark (
            $kernel_p3!($p3_params, $Ch2022, $output, $L_ice, $N_ice, $F_rim, $ρ_rim, $ρ_arr; ndrange = $ndrange);
            KernelAbstractions.synchronize($backend)
        )
        @info "P3 Kernel runtime benchmark:" BT.minimum(b_p3)
    end
end


for FT in (Float64, Float32)
    run_gpu_performance_benchmarks(FT)
end
