using Test
using KernelAbstractions

import CUDAKernels as CK

import CloudMicrophysics as CM
import CLIMAParameters as CP
import Thermodynamics as TD

const CMP = CM.Parameters

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

const AM = CM.AerosolModel
const AA = CM.AerosolActivation
const CMI_het = CM.HetIceNucleation
const CMI_hom = CM.HomIceNucleation
const CO = CM.Common
const CMT = CM.CommonTypes
const CM0 = CM.Microphysics0M
const CM1 = CM.Microphysics1M
const HN = CM.Nucleation
const P3 = CM.P3Scheme

const liquid = CMT.LiquidType()
const ice = CMT.IceType()
const rain = CMT.RainType()
const snow = CMT.SnowType()
const kaolinite = CMT.KaoliniteType()
const illite = CMT.IlliteType()

@info "GPU Tests"

if get(ARGS, 1, "Array") == "CuArray"
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
    prs,
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
        output[1, i] = CM1.accretion(prs, liquid, rain, ql[i], qr[i], ρ[i])
        output[2, i] = CM1.accretion(prs, ice, snow, qi[i], qs[i], ρ[i])
        output[3, i] = CM1.accretion(prs, liquid, snow, ql[i], qs[i], ρ[i])
        output[4, i] = CM1.accretion(prs, ice, rain, qi[i], qr[i], ρ[i])
        output[5, i] = CM1.accretion_rain_sink(prs, qi[i], qr[i], ρ[i])
        output[6, i] =
            CM1.accretion_snow_rain(prs, snow, rain, qs[i], qr[i], ρ[i])
        output[7, i] =
            CM1.accretion_snow_rain(prs, rain, snow, qr[i], qs[i], ρ[i])
    end
end

@kernel function test_1_moment_micro_snow_melt_kernel!(
    prs,
    output::AbstractArray{FT},
    ρ,
    T,
    qs,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = CM1.snow_melt(prs, qs[i], ρ[i], T[i])
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
    output::AbstractArray{FT},
    x_sulph,
    T,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = CO.H2SO4_soln_saturation_vapor_pressure(x_sulph[i], T[i])
    end
end

@kernel function test_Common_Delta_a_w_kernel!(
    prs,
    output::AbstractArray{FT},
    x_sulph,
    T,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = CO.Delta_a_w(prs, x_sulph[i], T[i])
    end
end

@kernel function test_IceNucleation_ABIFM_J_kernel!(
    output::AbstractArray{FT},
    Delta_a_w,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1] = CMI_het.ABIFM_J(kaolinite, Delta_a_w[1])
        output[2] = CMI_het.ABIFM_J(illite, Delta_a_w[2])
    end
end

@kernel function test_IceNucleation_homogeneous_J_kernel!(
    output::AbstractArray{FT},
    Delta_a_w,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1] = CMI_hom.homogeneous_J(Delta_a_w[1])
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

@kernel function test_p3_thresholds_kernel!(
    output::AbstractArray{FT},
    ρ_r::FT,
    F_r::FT) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = P3.thresholds(
            ρ_r,
            F_r,
        )[i]
    end
end

function test_gpu(FT)

    make_prs(::Type{FT}) where {FT} = cloud_microphysics_parameters(
        CP.create_toml_dict(FT; dict_type = "alias"),
    )

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
            make_prs(FT),
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
        event = kernel!(make_prs(FT), output, ρ, T, qs, ndrange = ndrange)
        wait(dev, event)

        # test if 1-moment snow melt is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(9.518235437405256e-6)
        @test Array(output)[2] ≈ FT(0)
        @test Array(output)[3] ≈ FT(0)
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
        event = kernel!(output, x_sulph, T, ndrange = ndrange)
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

        kernel! = test_Common_Delta_a_w_kernel!(dev, work_groups)
        event = kernel!(make_prs(FT), output, x_sulph, T, ndrange = ndrange)
        wait(dev, event)

        # test if Delta_a_w is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(0.2750536615)
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
        event = kernel!(output, Delta_a_w, ndrange = ndrange)
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
        event = kernel!(output, Delta_a_w, ndrange = ndrange)
        wait(dev, event)

        # test homogeneous_J is callable and returns a reasonable value
        @test Array(output)[1] ≈ FT(2.66194650334444e12)
    end

    @testset "Homogeneous nucleation kernels" begin
        make_nuc_prs(::Type{FT}) where {FT} =
            nucleation_parameters(CP.create_toml_dict(FT; dict_type = "alias"))

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
            make_nuc_prs(FT),
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
            make_nuc_prs(FT),
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
            make_nuc_prs(FT),
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

    @testset "P3 Kernels" begin
        data_length = 4
        output = ArrayType(Array{FT}(undef, 1, data_length))
        fill!(output, FT(-44.0))

        dev = device(ArrayType)
        work_groups = (1,)
        ndrange = (data_length,)

        ρ_r = FT(400)
        F_r = FT(0.5)

        kernel! = test_p3_thresholds_kernel!(
            dev,
            work_groups,
        )
        event = kernel!(output, ρ_r, F_r, ndrange = ndrange)
        wait(dev, event)

        # test thresholds is callable and returns a reasonable value
        @test Array(output)[1] ≈ FT(0.0004943308543980244)
        @test Array(output)[2] ≈ FT(0.00026324133585593025)
        @test Array(output)[3] ≈ FT(306.6678474961574)
        @test Array(output)[4] ≈ FT(213.3356949923147)

    end
end

println("Testing Float64")
test_gpu(Float64)

println("Testing Float32")
test_gpu(Float32)
