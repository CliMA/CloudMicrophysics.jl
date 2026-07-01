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

@kernel inbounds = true function benchmark_1m_accretion_kernel!(
    lcl, rain, icl, snow, e_lr, e_is, e_ls, e_ir, e_rs, coeff_disp, blk1mvel, output, ρ, qi, qs, ql, qr,
)
    i = @index(Global, Linear)
    liq_rai = CM1.accretion(lcl, rain, blk1mvel.rain, e_lr, ql[i], qr[i], ρ[i])
    ice_sno = CM1.accretion(icl, snow, blk1mvel.snow, e_is, qi[i], qs[i], ρ[i])
    liq_sno = CM1.accretion(lcl, snow, blk1mvel.snow, e_ls, ql[i], qs[i], ρ[i])
    ice_rai = CM1.accretion(icl, rain, blk1mvel.rain, e_ir, qi[i], qr[i], ρ[i])
    sno_rai = CM1.accretion_snow_rain(
        snow, rain, blk1mvel.snow, blk1mvel.rain, e_rs, coeff_disp, qs[i], qr[i], ρ[i],
    )
    output[i] = (; liq_rai, ice_sno, liq_sno, ice_rai, sno_rai)
end

@kernel inbounds = true function benchmark_2m_kernel!(
    KK2000, B1994, TC1980, LD2004, output, ql, qr, ρ, Nd,
)
    i = @index(Global, Linear)
    S_LD2004 = CM2.conv_q_lcl_to_q_rai(LD2004, ql[i], ρ[i], Nd[i])
    S_KK2000 = CM2.accretion(KK2000, ql[i], qr[i], ρ[i])
    S_B1994 = CM2.accretion(B1994, ql[i], qr[i], ρ[i])
    output[i] = (; S_LD2004, S_KK2000, S_B1994)
end

@kernel inbounds = true function benchmark_p3_kernel!(
    p3_params, vel_params, output, L_ice, N_ice, F_rim, ρ_rim, ρₐ, quad,
)
    i = @index(Global, Linear)
    state = P3.P3State(p3_params, L_ice[i], N_ice[i], F_rim[i], ρ_rim[i])
    logλ = P3.get_distribution_logλ(state)
    sc = P3.ice_self_collection(state, logλ, vel_params, ρₐ[i]; quad)
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

    # 2. Kernel Benchmarks (1-Moment, 2-Moment, P3)
    TT.@testset "Kernel Performance & Compile Time ($FT)" begin
        # 1-Moment setup
        lcl = CMP.CloudLiquid(FT)
        icl = CMP.CloudIce(FT)
        rain = CMP.Rain(FT)
        snow = CMP.Snow(FT)
        mp = CMP.Microphysics1MParams(FT)
        opts = mp.options
        e_lr = opts.cloud_liquid_rain_accretion.e
        e_is = opts.cloud_ice_snow_accretion.e
        e_ls = opts.cloud_liquid_snow_accretion.e
        e_ir = opts.cloud_ice_rain_accretion.e
        e_rs = opts.rain_snow_accretion.e
        coeff_disp = opts.rain_snow_accretion.coeff_disp
        blk1mvel = CMP.Blk1MVelType(FT)

        DT1 = @NamedTuple{liq_rai::FT, ice_sno::FT, liq_sno::FT, ice_rai::FT, sno_rai::FT}
        (; output, ndrange) = setup_output(1000, DT1)
        ρ = constant_data(FT(1.2); ndrange)
        qi = constant_data(FT(5e-4); ndrange)
        qs = constant_data(FT(5e-4); ndrange)
        ql = constant_data(FT(5e-4); ndrange)
        qr = constant_data(FT(5e-4); ndrange)

        kernel_1m! = benchmark_1m_accretion_kernel!(backend, work_groups)

        # Measure compile time (first execution)
        t_compile_1m = @elapsed kernel_1m!(
            lcl,
            rain,
            icl,
            snow,
            e_lr,
            e_is,
            e_ls,
            e_ir,
            e_rs,
            coeff_disp,
            blk1mvel,
            output,
            ρ,
            qi,
            qs,
            ql,
            qr;
            ndrange,
        )
        KernelAbstractions.synchronize(backend)
        @info "1-Moment Accretion Kernel first call (compile + run time): $(round(t_compile_1m, digits=4)) seconds"

        b_1m = BT.@benchmark (
            $kernel_1m!(
                $lcl,
                $rain,
                $icl,
                $snow,
                $e_lr,
                $e_is,
                $e_ls,
                $e_ir,
                $e_rs,
                $coeff_disp,
                $blk1mvel,
                $output,
                $ρ,
                $qi,
                $qs,
                $ql,
                $qr;
                ndrange = $ndrange,
            );
            KernelAbstractions.synchronize($backend)
        )
        @info "1-Moment Accretion Kernel runtime benchmark:" BT.minimum(b_1m)

        # 2-Moment setup
        KK2000 = CMP.KK2000(FT)
        B1994 = CMP.B1994(FT)
        TC1980 = CMP.TC1980(FT)
        LD2004 = CMP.LD2004(FT)

        DT2 = @NamedTuple{S_LD2004::FT, S_KK2000::FT, S_B1994::FT}
        (; output, ndrange) = setup_output(1000, DT2)
        Nd = constant_data(FT(1e8); ndrange)

        kernel_2m! = benchmark_2m_kernel!(backend, work_groups)
        t_compile_2m = @elapsed kernel_2m!(KK2000, B1994, TC1980, LD2004, output, ql, qr, ρ, Nd; ndrange)
        KernelAbstractions.synchronize(backend)
        @info "2-Moment Kernel first call (compile + run time): $(round(t_compile_2m, digits=4)) seconds"

        b_2m = BT.@benchmark (
            $kernel_2m!($KK2000, $B1994, $TC1980, $LD2004, $output, $ql, $qr, $ρ, $Nd; ndrange = $ndrange);
            KernelAbstractions.synchronize($backend)
        )
        @info "2-Moment Kernel runtime benchmark:" BT.minimum(b_2m)

        # P3 setup
        p3_params = CMP.ParametersP3(FT)
        Ch2022 = CMP.Chen2022VelType(FT)

        DT3 = @NamedTuple{logλ::FT, sc::@NamedTuple{dNdt::FT}}
        (; output, ndrange) = setup_output(1000, DT3)
        L_ice = constant_data(FT(5e-4 * 1.2); ndrange)

        N_ice = constant_data(FT(1e8 * 1.2); ndrange)
        F_rim = constant_data(FT(0.95); ndrange)
        ρ_rim = constant_data(FT(400.0); ndrange)

        _glq = P3.GaussLegendre(FT, 12)
        kernel_p3! = benchmark_p3_kernel!(backend, work_groups)
        t_compile_p3 = @elapsed kernel_p3!(p3_params, Ch2022, output, L_ice, N_ice, F_rim, ρ_rim, ρ, _glq; ndrange)
        KernelAbstractions.synchronize(backend)
        @info "P3 Kernel first call (compile + run time): $(round(t_compile_p3, digits=4)) seconds"

        b_p3 = BT.@benchmark (
            $kernel_p3!($p3_params, $Ch2022, $output, $L_ice, $N_ice, $F_rim, $ρ_rim, $ρ, $_glq; ndrange = $ndrange);
            KernelAbstractions.synchronize($backend)
        )
        @info "P3 Kernel runtime benchmark:" BT.minimum(b_p3)
    end
end


for FT in (Float64, Float32)
    run_gpu_performance_benchmarks(FT)
end
