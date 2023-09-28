using Test
using KernelAbstractions

import CUDAKernels as CK
using ClimaComms

import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CLIMAParameters as CP
import Thermodynamics as TD

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.HetIceNucleation as CMI_het
import CloudMicrophysics.HomIceNucleation as CMI_hom
import CloudMicrophysics.Common as CO
import CloudMicrophysics.CommonTypes as CMT
import CloudMicrophysics.Microphysics0M as CM0
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Nucleation as HN
import CloudMicrophysics.Parameters.H2S04NucleationParameters
import CloudMicrophysics.Parameters.MixedNucleationParameters
import CloudMicrophysics.Parameters.OrganicNucleationParameters
import CloudMicrophysics.P3Scheme as P3

const kaolinite = CMT.KaoliniteType()
const illite = CMT.IlliteType()

@info "GPU Tests"

if ClimaComms.device() isa ClimaComms.CUDADevice
    import CUDA
    ArrayType = CUDA.CuArray
    CUDA.allowscalar(false)
    device(::Type{T}) where {T <: CUDA.CuArray} = CK.CUDADevice()
else
    ArrayType = Array
    device(::Type{T}) where {T <: Array} = CK.CPU()
end
@show ArrayType

@kernel function test_aerosol_activation_kernel!(
    prs,
    output::AbstractArray{FT},
    r,
    stdev,
    N,
    ϵ,
    ϕ,
    M,
    ν,
    ρ,
    κ,
) where {FT}

    i = @index(Group, Linear)
    thermo_params = CMP.thermodynamics_params(prs)
    # atmospheric conditions (taken from aerosol activation tests)
    T = FT(294)      # air temperature K
    p = FT(1e5)    # air pressure Pa
    w = FT(0.5)         # vertical velocity m/s
    p_vs = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    q_vs = 1 / (1 - CMP.molmass_ratio(prs) * (p_vs - p) / p_vs)
    # water vapor specific humidity (saturated)
    q = TD.PhasePartition(q_vs)

    args = (T, p, w, q)

    @inbounds begin
        mode_B = AM.Mode_B(
            r[i],
            stdev[i],
            N[i],
            (FT(1.0),),
            (ϵ[i],),
            (ϕ[i],),
            (M[i],),
            (ν[i],),
            (ρ[i],),
            1,
        )
        mode_κ = AM.Mode_κ(
            r[i],
            stdev[i],
            N[i],
            (FT(1.0),),
            (FT(1.0),),
            (M[i],),
            (κ[i],),
            1,
        )

        arsl_dst_B = AM.AerosolDistribution((mode_B,))
        arsl_dst_κ = AM.AerosolDistribution((mode_κ,))

        output[1, i] = AA.mean_hygroscopicity_parameter(prs, arsl_dst_B)[1]
        output[2, i] = AA.mean_hygroscopicity_parameter(prs, arsl_dst_κ)[1]

        output[3, i] = AA.total_N_activated(prs, arsl_dst_B, args...)
        output[4, i] = AA.total_N_activated(prs, arsl_dst_κ, args...)

        output[5, i] = AA.total_M_activated(prs, arsl_dst_B, args...)
        output[6, i] = AA.total_M_activated(prs, arsl_dst_κ, args...)
    end
end

@kernel function test_0_moment_micro_kernel!(
    prs,
    output::AbstractArray{FT},
    liquid_frac,
    qc,
    qt,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        ql = qc[i] * liquid_frac[i]
        qi = (1 - liquid_frac[i]) * qc[i]
        q = TD.PhasePartition(FT(qt[i]), ql, qi)

        output[1, i] = CM0.remove_precipitation(prs, q)

        _τ_precip = CMP.τ_precip(prs)
        _qc_0 = CMP.qc_0(prs)

        output[2, i] = -max(0, ql + qi - _qc_0) / _τ_precip
    end
end

@kernel function test_1_moment_micro_accretion_kernel!(
    liquid,
    rain,
    ice,
    snow,
    ce,
    rain_blk1mvel,
    snow_blk1mvel,
    output::AbstractArray{FT},
    ρ,
    qt,
    qi,
    qs,
    ql,
    qr,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1, i] = CM1.accretion(ce, liquid, rain, ql[i], qr[i], ρ[i])
        output[2, i] = CM1.accretion(ce, ice, snow, qi[i], qs[i], ρ[i])
        output[3, i] = CM1.accretion(ce, liquid, snow, ql[i], qs[i], ρ[i])
        output[4, i] = CM1.accretion(ce, ice, rain, qi[i], qr[i], ρ[i])
        output[5, i] =
            CM1.accretion_rain_sink(rain, ice, ce, qi[i], qr[i], ρ[i])
        output[6, i] = CM1.accretion_snow_rain(
            ce,
            snow,
            rain,
            snow_blk1mvel,
            rain_blk1mvel,
            qs[i],
            qr[i],
            ρ[i],
        )
        output[7, i] = CM1.accretion_snow_rain(
            ce,
            rain,
            snow,
            rain_blk1mvel,
            snow_blk1mvel,
            qr[i],
            qs[i],
            ρ[i],
        )
    end
end

@kernel function test_1_moment_micro_snow_melt_kernel!(
    snow,
    air_props,
    thermo_params,
    output::AbstractArray{FT},
    ρ,
    T,
    qs,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] =
            CM1.snow_melt(snow, air_props, thermo_params, qs[i], ρ[i], T[i])
    end
end

@kernel function test_2_moment_acnv_kernel!(
    prs,
    acnvKK2000,
    acnvB1994,
    acnvTC1980,
    acnvLD2004,
    acnvVarT,
    output::AbstractArray{FT},
    ql,
    ρ,
    Nd,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1, i] =
            CM2.conv_q_liq_to_q_rai(prs, acnvVarT, ql[i], ρ[i], N_d = Nd[i])
        output[2, i] =
            CM2.conv_q_liq_to_q_rai(acnvLD2004, ql[i], ρ[i], N_d = Nd[i])
        output[3, i] =
            CM2.conv_q_liq_to_q_rai(acnvTC1980, ql[i], ρ[i], N_d = Nd[i])
        output[4, i] =
            CM2.conv_q_liq_to_q_rai(acnvB1994, ql[i], ρ[i], N_d = Nd[i])
        output[5, i] =
            CM2.conv_q_liq_to_q_rai(acnvKK2000, ql[i], ρ[i], N_d = Nd[i])
    end
end

@kernel function test_2_moment_accr_kernel!(
    accKK2000,
    accB1994,
    accTC1980,
    output::AbstractArray{FT},
    ql,
    qr,
    ρ,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1, i] = CM2.accretion(accKK2000, ql[i], qr[i], ρ[i])
        output[2, i] = CM2.accretion(accB1994, ql[i], qr[i], ρ[i])
        output[3, i] = CM2.accretion(accTC1980, ql[i], qr[i])
    end
end

@kernel function test_2_moment_SB2006_kernel!(
    rain,
    air_props,
    thermo_params,
    acnvSB2006,
    accrSB2006,
    scSB2006,
    brkSB2006,
    tvSB2006,
    evpSB2006,
    output::AbstractArray{FT},
    qt,
    ql,
    qr,
    Nl,
    Nr,
    ρ,
    T,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        q = TD.PhasePartition(FT(qt[i]), FT(ql[i]), FT(0))

        output[1, i] =
            CM2.autoconversion_and_liquid_self_collection(
                acnvSB2006,
                ql[i],
                qr[i],
                ρ[i],
                Nl[i],
            ).au.dq_liq_dt
        output[2, i] =
            CM2.autoconversion_and_liquid_self_collection(
                acnvSB2006,
                ql[i],
                qr[i],
                ρ[i],
                Nl[i],
            ).au.dN_liq_dt
        output[3, i] =
            CM2.autoconversion_and_liquid_self_collection(
                acnvSB2006,
                ql[i],
                qr[i],
                ρ[i],
                Nl[i],
            ).au.dq_rai_dt
        output[4, i] =
            CM2.autoconversion_and_liquid_self_collection(
                acnvSB2006,
                ql[i],
                qr[i],
                ρ[i],
                Nl[i],
            ).au.dN_rai_dt

        output[5, i] =
            CM2.autoconversion_and_liquid_self_collection(
                acnvSB2006,
                ql[i],
                qr[i],
                ρ[i],
                Nl[i],
            ).sc
        output[6, i] =
            CM2.accretion(accrSB2006, ql[i], qr[i], ρ[i], Nl[i]).dq_liq_dt
        output[7, i] =
            CM2.accretion(accrSB2006, ql[i], qr[i], ρ[i], Nl[i]).dN_liq_dt
        output[8, i] =
            CM2.accretion(accrSB2006, ql[i], qr[i], ρ[i], Nl[i]).dq_rai_dt
        output[9, i] =
            CM2.accretion(accrSB2006, ql[i], qr[i], ρ[i], Nl[i]).dN_rai_dt
        output[10, i] =
            CM2.rain_self_collection_and_breakup(
                scSB2006,
                brkSB2006,
                tvSB2006,
                qr[i],
                ρ[i],
                Nr[i],
            ).sc
        output[11, i] =
            CM2.rain_self_collection_and_breakup(
                scSB2006,
                brkSB2006,
                tvSB2006,
                qr[i],
                ρ[i],
                Nr[i],
            ).br
        output[12, i] =
            CM2.rain_terminal_velocity(rain, tvSB2006, qr[i], ρ[i], Nr[i])[1]
        output[13, i] =
            CM2.rain_terminal_velocity(rain, tvSB2006, qr[i], ρ[i], Nr[i])[2]
        #TODO - I think the incomplete gamma function doesnt work on the GPU
        #output[14, i] = CM2.rain_evaporation(
        #    evpSB2006,
        #    air_props,
        #    thermo_params,
        #    q,
        #    qr[i],
        #    ρ[i],
        #    Nr[i],
        #    T[i],
        #)[1]
        #output[15, i] = CM2.rain_evaporation(
        #    evpSB2006,
        #    air_props,
        #    thermo_params,
        #    q,
        #    qr[i],
        #    ρ[i],
        #    Nr[i],
        #    T[i],
        #)[2]
    end
end

@kernel function test_h2so4_nucleation_kernel!(
    prs,
    output::AbstractArray{FT},
    h2so4_conc,
    nh3_conc,
    negative_ion_conc,
    temp,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = sum(
            HN.h2so4_nucleation_rate(
                h2so4_conc[i],
                nh3_conc[i],
                negative_ion_conc[i],
                temp[i],
                prs,
            ),
        )
    end
end

@kernel function test_organic_nucleation_kernel!(
    prs,
    output::AbstractArray{FT},
    negative_ion_conc,
    monoterpene_conc,
    O3_conc,
    OH_conc,
    temp,
    condensation_sink,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = HN.organic_nucleation_rate(
            negative_ion_conc[i],
            monoterpene_conc[i],
            O3_conc[i],
            OH_conc[i],
            temp[i],
            condensation_sink[i],
            prs,
        )
    end
end

@kernel function test_organic_and_h2so4_nucleation_kernel!(
    prs,
    output::AbstractArray{FT},
    h2so4_conc,
    monoterpene_conc,
    OH_conc,
    temp,
    condensation_sink,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = HN.organic_and_h2so4_nucleation_rate(
            h2so4_conc[i],
            monoterpene_conc[i],
            OH_conc[i],
            temp[i],
            condensation_sink[i],
            prs,
        )
    end
end

@kernel function test_apparent_nucleation_rate_kernel!(
    output::AbstractArray{FT},
    output_diam,
    nucleation_rate,
    condensation_growth_rate,
    coag_sink,
    coag_sink_input_diam,
    input_diam,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = HN.apparent_nucleation_rate(
            output_diam[i],
            nucleation_rate[i],
            condensation_growth_rate[i],
            coag_sink[i],
            coag_sink_input_diam[i],
            input_diam[i],
        )
    end
end

@kernel function test_Common_H2SO4_soln_saturation_vapor_pressure_kernel!(
    prs,
    output::AbstractArray{FT},
    x_sulph,
    T,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] =
            CO.H2SO4_soln_saturation_vapor_pressure(prs, x_sulph[i], T[i])
    end
end

@kernel function test_Common_a_w_xT_kernel!(
    prs,
    output::AbstractArray{FT},
    x_sulph,
    T,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = CO.a_w_xT(prs, x_sulph[i], T[i])
    end
end

@kernel function test_Common_a_w_eT_kernel!(
    prs,
    output::AbstractArray{FT},
    e,
    T,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = CO.a_w_eT(prs, e[i], T[i])
    end
end

@kernel function test_Common_a_w_ice_kernel!(
    prs,
    output::AbstractArray{FT},
    T,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = CO.a_w_ice(prs, T[i])
    end
end

@kernel function test_IceNucleation_ABIFM_J_kernel!(
    prs,
    output::AbstractArray{FT},
    Delta_a_w,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1] = CMI_het.ABIFM_J(prs, kaolinite, Delta_a_w[1])
        output[2] = CMI_het.ABIFM_J(prs, illite, Delta_a_w[2])
    end
end

@kernel function test_IceNucleation_homogeneous_J_kernel!(
    prs,
    output::AbstractArray{FT},
    Delta_a_w,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1] = CMI_hom.homogeneous_J(prs, Delta_a_w[1])
    end
end

@kernel function test_h2so4_nucleation_kernel!(
    prs,
    output::AbstractArray{FT},
    h2so4_conc,
    nh3_conc,
    negative_ion_conc,
    temp,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = sum(
            HN.h2so4_nucleation_rate(
                h2so4_conc[i],
                nh3_conc[i],
                negative_ion_conc[i],
                temp[i],
                prs,
            ),
        )
    end
end

@kernel function test_organic_nucleation_kernel!(
    prs,
    output::AbstractArray{FT},
    negative_ion_conc,
    monoterpene_conc,
    O3_conc,
    OH_conc,
    temp,
    condensation_sink,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = HN.organic_nucleation_rate(
            negative_ion_conc[i],
            monoterpene_conc[i],
            O3_conc[i],
            OH_conc[i],
            temp[i],
            condensation_sink[i],
            prs,
        )
    end
end

@kernel function test_organic_and_h2so4_nucleation_kernel!(
    prs,
    output::AbstractArray{FT},
    h2so4_conc,
    monoterpene_conc,
    OH_conc,
    temp,
    condensation_sink,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = HN.organic_and_h2so4_nucleation_rate(
            h2so4_conc[i],
            monoterpene_conc[i],
            OH_conc[i],
            temp[i],
            condensation_sink[i],
            prs,
        )
    end
end

@kernel function test_apparent_nucleation_rate_kernel!(
    output::AbstractArray{FT},
    output_diam,
    nucleation_rate,
    condensation_growth_rate,
    coag_sink,
    coag_sink_input_diam,
    input_diam,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = HN.apparent_nucleation_rate(
            output_diam[i],
            nucleation_rate[i],
            condensation_growth_rate[i],
            coag_sink[i],
            coag_sink_input_diam[i],
            input_diam[i],
        )
    end
end

@kernel function test_P3_scheme_kernel!(
    prs,
    output::AbstractArray{FT},
    F_r,
    ρ_r,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1, i] = P3.thresholds(prs, ρ_r[i], F_r[i])[1]
        output[2, i] = P3.thresholds(prs, ρ_r[i], F_r[i])[2]
    end
end

function test_gpu(FT)

    make_prs(::Type{FT}) where {FT} = cloud_microphysics_parameters(
        CP.create_toml_dict(FT; dict_type = "alias"),
    )
    prs = make_prs(FT)
    air_props = CMT.AirProperties(FT)
    thermo_params = CMP.thermodynamics_params(prs)

    liquid = CMT.LiquidType()
    ice = CMT.IceType(FT)
    rain = CMT.RainType(FT)
    snow = CMT.SnowType(FT)

    ce = CMT.CollisionEfficiency(FT)
    p3_prs = CMP.CloudMicrophysicsParametersP3(FT)
    rain_blk1mvel = CMT.Blk1MVelType(FT, rain)
    snow_blk1mvel = CMT.Blk1MVelType(FT, snow)

    acnvSB2006 = CMT.AutoconversionSB2006(FT)
    accrSB2006 = CMT.AccretionSB2006(FT)
    scSB2006 = CMT.SelfCollectionSB2006(FT)
    brkSB2006 = CMT.BreakupSB2006(FT)
    tvSB2006 = CMT.TerminalVelocitySB2006(FT)
    evpSB2006 = CMT.EvaporationSB2006(FT)
    acnvKK2000 = CMT.AutoconversionKK2000(FT)
    acnvB1994 = CMT.AutoconversionB1994(FT)
    acnvTC1980 = CMT.AutoconversionTC1980(FT)
    acnvLD2004 = CMT.AutoconversionLD2004(FT)
    acnvVarT = CMT.VarTimeScaleAcnvType()
    accKK2000 = CMT.AccretionKK2000(FT)
    accB1994 = CMT.AccretionB1994(FT)
    accTC1980 = CMT.AccretionTC1980(FT)

    @testset "Aerosol activation kernels" begin
        data_length = 2
        output = ArrayType(Array{FT}(undef, 6, data_length))
        fill!(output, FT(-44.0))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        r = ArrayType([FT(0.243 * 1e-6), FT(1.5 * 1e-6)])
        stdev = ArrayType([FT(1.4), FT(2.1)])
        N = ArrayType([FT(100 * 1e6), FT(1 * 1e6)])
        ϵ = ArrayType([FT(1), FT(1)])
        ϕ = ArrayType([FT(1), FT(0.9)])
        M = ArrayType([FT(0.132), FT(0.058443)])
        ν = ArrayType([FT(3), FT(2)])
        ρ = ArrayType([FT(1770), FT(2170)])
        κ = ArrayType([FT(0.53), FT(1.12)])

        kernel! = test_aerosol_activation_kernel!(dev, work_groups)
        event = kernel!(
            make_prs(FT),
            output,
            r,
            stdev,
            N,
            ϵ,
            ϕ,
            M,
            ν,
            ρ,
            κ,
            ndrange = ndrange,
        )
        wait(dev, event)

        # test if all aerosol activation output is positive
        @test all(Array(output)[:, :] .>= FT(0))
        # test if higroscopicity parameter is the same for κ and B modes
        @test all(
            isapprox(Array(output)[1, :], Array(output)[2, :], rtol = 0.3),
        )
        # test if the number and mass activated are the same for κ and B modes
        @test all(
            isapprox(Array(output)[3, :], Array(output)[4, :], rtol = 1e-5),
        )
        @test all(
            isapprox(Array(output)[5, :], Array(output)[6, :], rtol = 1e-5),
        )
    end

    @testset "0-moment microphysics kernels" begin
        data_length = 3
        output = ArrayType(Array{FT}(undef, 2, data_length))
        fill!(output, FT(-44.0))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        liquid_frac = ArrayType([FT(0), FT(0.5), FT(1)])
        qt = ArrayType([FT(13e-3), FT(13e-3), FT(13e-3)])
        qc = ArrayType([FT(3e-3), FT(4e-3), FT(5e-3)])

        kernel! = test_0_moment_micro_kernel!(dev, work_groups)
        event = kernel!(
            make_prs(FT),
            output,
            liquid_frac,
            qc,
            qt,
            ndrange = ndrange,
        )
        wait(dev, event)

        # test 0-moment rain removal is callable and returns a reasonable value
        @test all(isequal(Array(output)[1, :], Array(output)[2, :]))
    end

    @testset "1-moment microphysics kernels" begin
        data_length = 2
        output = ArrayType(Array{FT}(undef, 7, data_length))
        fill!(output, FT(-44.0))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        ρ = ArrayType([FT(1.2), FT(1.2)])
        qt = ArrayType([FT(0), FT(20e-3)])
        qi = ArrayType([FT(0), FT(5e-4)])
        qs = ArrayType([FT(0), FT(5e-4)])
        ql = ArrayType([FT(0), FT(5e-4)])
        qr = ArrayType([FT(0), FT(5e-4)])

        kernel! = test_1_moment_micro_accretion_kernel!(dev, work_groups)
        event = kernel!(
            liquid,
            rain,
            ice,
            snow,
            ce,
            rain_blk1mvel,
            snow_blk1mvel,
            output,
            ρ,
            qt,
            qi,
            qs,
            ql,
            qr,
            ndrange = ndrange,
        )
        wait(dev, event)

        # test 1-moment accretion is callable and returns a reasonable value
        @test all(Array(output)[:, 1] .== FT(0))
        @test Array(output)[1, 2] ≈ FT(1.4150106417043544e-6)
        @test Array(output)[2, 2] ≈ FT(2.453070979562392e-7)
        @test Array(output)[3, 2] ≈ FT(2.453070979562392e-7)
        @test Array(output)[4, 2] ≈ FT(1.768763302130443e-6)
        @test Array(output)[5, 2] ≈ FT(3.085229094251214e-5)
        @test Array(output)[6, 2] ≈ FT(2.1705865794293408e-4)
        @test Array(output)[7, 2] ≈ FT(6.0118801860768854e-5)

        data_length = 3
        output = ArrayType(Array{FT}(undef, 1, data_length))
        fill!(output, FT(-44.0))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        ρ = ArrayType([FT(1.2), FT(1.2), FT(1.2)])
        T = ArrayType([FT(273.15 + 2), FT(273.15 + 2), FT(273.15 - 2)])
        qs = ArrayType([FT(1e-4), FT(0), FT(1e-4)])

        kernel! = test_1_moment_micro_snow_melt_kernel!(dev, work_groups)
        event = kernel!(
            snow,
            air_props,
            thermo_params,
            output,
            ρ,
            T,
            qs,
            ndrange = ndrange,
        )
        wait(dev, event)

        # test if 1-moment snow melt is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(9.518235437405256e-6)
        @test Array(output)[2] ≈ FT(0)
        @test Array(output)[3] ≈ FT(0)
    end

    @testset "2-moment microphysics kernels" begin

        data_length = 1
        output = ArrayType(Array{FT}(undef, 5, data_length))
        fill!(output, FT(-44))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        ql = ArrayType([FT(2e-3)])
        ρ = ArrayType([FT(1.2)])
        Nd = ArrayType([FT(1e8)])

        kernel! = test_2_moment_acnv_kernel!(dev, work_groups)
        event = kernel!(
            prs,
            acnvKK2000,
            acnvB1994,
            acnvTC1980,
            acnvLD2004,
            acnvVarT,
            output,
            ql,
            ρ,
            Nd,
            ndrange = ndrange,
        )
        wait(dev, event)

        @test Array(output)[1] ≈ FT(2e-6)
        @test Array(output)[2] ≈ FT(1.6963072465911614e-6)
        @test Array(output)[3] ≈ FT(3.5482867084128596e-6)
        @test Array(output)[4] ≈ FT(9.825462758968215e-7)
        @test Array(output)[5] ≈ FT(5.855332513368727e-8)

        data_length = 1
        output = ArrayType(Array{FT}(undef, 3, data_length))
        fill!(output, FT(-44))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        ql = ArrayType([FT(2e-3)])
        qr = ArrayType([FT(5e-4)])
        ρ = ArrayType([FT(1.2)])

        kernel! = test_2_moment_accr_kernel!(dev, work_groups)
        event = kernel!(
            accKK2000,
            accB1994,
            accTC1980,
            output,
            ql,
            qr,
            ρ,
            ndrange = ndrange,
        )
        wait(dev, event)

        @test isapprox(Array(output)[1], FT(6.6548664e-6), rtol = 1e-6)
        @test Array(output)[2] ≈ FT(7.2e-6)
        @test Array(output)[3] ≈ FT(4.7e-6)

        data_length = 1
        output = ArrayType(Array{FT}(undef, 13, data_length))
        fill!(output, FT(-44))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        T = ArrayType([FT(290)])
        qt = ArrayType([FT(7e-3)])
        ql = ArrayType([FT(2e-3)])
        qr = ArrayType([FT(5e-4)])
        ρ = ArrayType([FT(1.2)])
        Nl = ArrayType([FT(1e8)])
        Nr = ArrayType([FT(1e7)])

        kernel! = test_2_moment_SB2006_kernel!(dev, work_groups)
        event = kernel!(
            rain,
            air_props,
            thermo_params,
            acnvSB2006,
            accrSB2006,
            scSB2006,
            brkSB2006,
            tvSB2006,
            evpSB2006,
            output,
            qt,
            ql,
            qr,
            Nl,
            Nr,
            ρ,
            T,
            ndrange = ndrange,
        )
        wait(dev, event)

        @test isapprox(Array(output)[1], FT(-4.083606e-7), rtol = 1e-6)
        @test isapprox(Array(output)[2], FT(-3769.4827), rtol = 1e-6)
        @test isapprox(Array(output)[3], FT(4.083606e-7), rtol = 1e-6)
        @test isapprox(Array(output)[4], FT(1884.7413), rtol = 1e-6)
        @test isapprox(Array(output)[5], FT(-31040.115), rtol = 1e-6)
        @test isapprox(Array(output)[6], FT(-6.358926e-6), rtol = 1e-6)
        @test isapprox(Array(output)[7], FT(-317946.28), rtol = 1e-6)
        @test isapprox(Array(output)[8], FT(6.358926e-6), rtol = 1e-6)
        @test isapprox(Array(output)[9], FT(0.0), rtol = 1e-6)
        @test isapprox(Array(output)[10], FT(-21187.494), rtol = 1e-6)
        @test isapprox(Array(output)[11], FT(14154.027), rtol = 1e-6)
        @test isapprox(Array(output)[12], FT(0.9868878), rtol = 1e-6)
        @test isapprox(Array(output)[13], FT(4.517734), rtol = 1e-6)
        #TODO - doesn't work on the GPU
        #@test isapprox(Array(output)[14], FT(-241447.55), rtol = 1e-6)
        #@test isapprox(Array(output)[15], FT(-0.0034601581), rtol = 1e-6)
    end

    @testset "Common Kernels" begin
        data_length = 1
        output = ArrayType(Array{FT}(undef, 1, data_length))
        fill!(output, FT(-44))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        T = ArrayType([FT(230)])
        x_sulph = ArrayType([FT(0.1)])

        kernel! = test_Common_H2SO4_soln_saturation_vapor_pressure_kernel!(
            dev,
            work_groups,
        )
        event = kernel!(make_prs(FT), output, x_sulph, T, ndrange = ndrange)
        wait(dev, event)

        # test H2SO4_soln_saturation_vapor_pressure is callable and returns a reasonable value
        @test Array(output)[1] ≈ FT(12.685507586924)

        data_length = 1
        output = ArrayType(Array{FT}(undef, 1, data_length))
        fill!(output, FT(-44.0))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        T = ArrayType([FT(230)])
        x_sulph = ArrayType([FT(0.1)])

        kernel! = test_Common_a_w_xT_kernel!(dev, work_groups)
        event = kernel!(make_prs(FT), output, x_sulph, T, ndrange = ndrange)
        wait(dev, event)

        # test if a_w_xT is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(0.92824538441)

        data_length = 1
        output = ArrayType(Array{FT}(undef, 1, data_length))
        fill!(output, FT(-44.0))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        T = ArrayType([FT(282)])
        e = ArrayType([FT(1001)])

        kernel! = test_Common_a_w_eT_kernel!(dev, work_groups)
        event = kernel!(make_prs(FT), output, e, T, ndrange = ndrange)
        wait(dev, event)

        # test if a_w_eT is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(0.880978146)

        data_length = 1
        output = ArrayType(Array{FT}(undef, 1, data_length))
        fill!(output, FT(-44.0))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        T = ArrayType([FT(230)])

        kernel! = test_Common_a_w_ice_kernel!(dev, work_groups)
        event = kernel!(make_prs(FT), output, T, ndrange = ndrange)
        wait(dev, event)

        # test if a_w_ice is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(0.653191723)

    end

    @testset "Ice Nucleation kernels" begin
        data_length = 2
        output = ArrayType(Array{FT}(undef, 1, data_length))
        fill!(output, FT(-44.0))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        Delta_a_w = ArrayType([FT(0.16), FT(0.15)])

        kernel! = test_IceNucleation_ABIFM_J_kernel!(dev, work_groups)
        event = kernel!(make_prs(FT), output, Delta_a_w, ndrange = ndrange)
        wait(dev, event)

        # test if ABIFM_J is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(153.65772539109)
        @test Array(output)[2] ≈ FT(31.870032033791)

        data_length = 1
        output = ArrayType(Array{FT}(undef, 1, data_length))
        fill!(output, FT(-44))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        T = ArrayType([FT(220)])
        x_sulph = ArrayType([FT(0.15)])
        Delta_a_w = ArrayType([FT(0.2907389666103033)])

        kernel! = test_IceNucleation_homogeneous_J_kernel!(dev, work_groups)
        event = kernel!(make_prs(FT), output, Delta_a_w, ndrange = ndrange)
        wait(dev, event)

        # test homogeneous_J is callable and returns a reasonable value
        @test Array(output)[1] ≈ FT(2.66194650334444e12)
    end

    @testset "Homogeneous nucleation kernels" begin
        data_length = 2
        output = ArrayType(Array{FT}(undef, 1, data_length))
        fill!(output, FT(-44.0))
        ndrange = (data_length,)

        dev = device(ArrayType)
        work_groups = (1,)

        # h2so4 nucleation
        h2so4_conc = ArrayType([FT(1e12), FT(1e12)])
        nh3_conc = ArrayType([FT(1), FT(1)])
        negative_ion_conc = ArrayType([FT(1), FT(1)])
        temp = ArrayType([FT(208), FT(208)])

        kernel! = test_h2so4_nucleation_kernel!(dev, work_groups)
        event = kernel!(
            H2S04NucleationParameters(FT),
            output,
            h2so4_conc,
            nh3_conc,
            negative_ion_conc,
            temp,
            ndrange = ndrange,
        )
        wait(dev, event)

        @test all(Array(output) .> FT(0))

        # Organic nucleation
        fill!(output, FT(-44.0))

        negative_ion_conc = ArrayType([FT(0.0), FT(0.0)])
        monoterpene_conc = ArrayType([FT(1e24), FT(1e24)])
        O3_conc = ArrayType([FT(1e24), FT(1e24)])
        OH_conc = ArrayType([FT(1e24), FT(1e24)])
        temp = ArrayType([FT(300), FT(300)])
        condensation_sink = ArrayType([FT(1), FT(1)])

        kernel! = test_organic_nucleation_kernel!(dev, work_groups)
        event = kernel!(
            OrganicNucleationParameters(FT),
            output,
            negative_ion_conc,
            monoterpene_conc,
            O3_conc,
            OH_conc,
            temp,
            condensation_sink,
            ndrange = ndrange,
        )
        wait(dev, event)

        @test all(Array(output) .> FT(0))

        # Organic and h2so4 nucleation
        fill!(output, FT(-44.0))

        h2so4_conc = ArrayType([FT(2.6e6), FT(2.6e6)])
        monoterpene_conc = ArrayType([FT(1), FT(1)])
        OH_conc = ArrayType([FT(1), FT(1)])
        temp = ArrayType([FT(300), FT(300)])
        condensation_sink = ArrayType([FT(1), FT(1)])

        kernel! = test_organic_and_h2so4_nucleation_kernel!(dev, work_groups)
        event = kernel!(
            MixedNucleationParameters(FT),
            output,
            h2so4_conc,
            monoterpene_conc,
            OH_conc,
            temp,
            condensation_sink,
            ndrange = ndrange,
        )
        wait(dev, event)

        @test all(Array(output) .> FT(0))

        # Apparent nucleation rate
        fill!(output, FT(-44.0))

        output_diam = ArrayType([FT(1e-9), FT(1e-9)])
        nucleation_rate = ArrayType([FT(1e6), FT(1e6)])
        condensation_growth_rate = ArrayType([FT(1e6), FT(1e6)])
        coag_sink = ArrayType([FT(1e6), FT(1e6)])
        coag_sink_input_diam = ArrayType([FT(1e-9), FT(1e-9)])
        input_diam = ArrayType([FT(1e-9), FT(1e-9)])

        kernel! = test_apparent_nucleation_rate_kernel!(dev, work_groups)
        event = kernel!(
            output,
            output_diam,
            nucleation_rate,
            condensation_growth_rate,
            coag_sink,
            coag_sink_input_diam,
            input_diam,
            ndrange = ndrange,
        )
        wait(dev, event)

        @test all(Array(output) .> FT(0))
    end
    #=
    @testset "P3 scheme kernels" begin
        data_length = 2
        output = ArrayType(Array{FT}(undef, 2, data_length))
        fill!(output, FT(-44.0))

        ndrange = (data_length,)
        dev = device(ArrayType)
        work_groups = (1,)

        F_r = ArrayType([FT(0.5), FT(0.95)])
        ρ_r = ArrayType([FT(400), FT(800)])

        kernel! = test_P3_scheme_kernel!(dev, work_groups)
        event = kernel!(p3_prs, output, F_r, ρ_r, ndrange = ndrange)
        wait(dev, event)

        # test if all output is positive...
        @test all(Array(output) .> FT(0))
        #... and returns reasonable numbers
        @test isapprox(
            Array(output)[1, 1],
            FT(0.4946323381999426 * 1e-3),
            rtol = 1e-2,
        )
        @test isapprox(
            Array(output)[2, 1],
            FT(0.26151186272014415 * 1e-3),
            rtol = 1e-2,
        )
        @test isapprox(
            Array(output)[1, 2],
            FT(1.7400778369620664 * 1e-3),
            rtol = 1e-2,
        )
        @test isapprox(
            Array(output)[2, 2],
            FT(0.11516682512848 * 1e-3),
            rtol = 1e-2,
        )
    end
    =#
end

println("Testing Float64")
test_gpu(Float64)

println("Testing Float32")
test_gpu(Float32)
