using Test
using KernelAbstractions
using CUDAKernels

import CloudMicrophysics
import CLIMAParameters
import Thermodynamics

const CM = CloudMicrophysics
const CP = CLIMAParameters
const TD = Thermodynamics

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const prs = EarthParameterSet()

const AM = CloudMicrophysics.AerosolModel
const AA = CloudMicrophysics.AerosolActivation
const CM0 = CloudMicrophysics.Microphysics_0M
const CM1 = CloudMicrophysics.Microphysics_1M
const CP0 = CLIMAParameters.Atmos.Microphysics_0M

const liquid = CM.Microphysics_1M.LiquidType()
const ice = CM.Microphysics_1M.IceType()
const rain = CM.Microphysics_1M.RainType()
const snow = CM.Microphysics_1M.SnowType()

if get(ARGS, 1, "Array") == "CuArray"
    using CUDA
    ArrayType = CUDA.CuArray
    CUDA.allowscalar(false)
    device(::Type{T}) where {T <: CuArray} = CUDADevice()
else
    ArrayType = Array
    device(::Type{T}) where {T <: Array} = CPU()
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

    # atmospheric conditions (taken from aerosol activation tests)
    T::FT = 294.0       # air temperature K
    p::FT = 100000.0    # air pressure Pa
    w::FT = 0.5         # vertical velocity m/s
    p_vs::FT = TD.saturation_vapor_pressure(prs, T, TD.Liquid())
    q_vs::FT = 1 / (1 - CP.Planet.molmass_ratio(prs) * (p_vs - p) / p_vs)
    # water vapor specific humidity (saturated)
    q = TD.PhasePartition(q_vs, FT(0.0), FT(0.0))

    args = (T, p, w, q)

    @inbounds begin
        mode_B = AM.Mode_B(
            FT(r[i]),
            FT(stdev[i]),
            FT(N[i]),
            (FT(1.0),),
            (FT(ϵ[i]),),
            (FT(ϕ[i]),),
            (FT(M[i]),),
            (FT(ν[i]),),
            (FT(ρ[i]),),
            1,
        )
        mode_κ = AM.Mode_κ(
            FT(r[i]),
            FT(stdev[i]),
            FT(N[i]),
            (FT(1.0),),
            (FT(1.0),),
            (FT(M[i]),),
            (FT(κ[i]),),
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
        ql::FT = qc[i] * liquid_frac[i]
        qi::FT = (1 - liquid_frac[i]) * qc[i]
        q = TD.PhasePartition(FT(qt[i]), ql, qi)

        output[1, i] = CM0.remove_precipitation(prs, q)

        _τ_precip = CP0.τ_precip(prs)
        _qc_0 = CP0.qc_0(prs)

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

@testset "Aerosol activation kernels" begin
    FT = Float32
    data_length = 2
    output = ArrayType(Array{FT}(undef, 6, data_length))
    fill!(output, -44.0)

    dev = device(ArrayType)
    work_groups = (1,)
    ndrange = (data_length,)

    r = ArrayType([0.243 * 1e-6, 1.5 * 1e-6])
    stdev = ArrayType([1.4, 2.1])
    N = ArrayType([100.0 * 1e6, 1.0 * 1e6])
    ϵ = ArrayType([1.0, 1.0])
    ϕ = ArrayType([1.0, 0.9])
    M = ArrayType([0.132, 0.058443])
    ν = ArrayType([3.0, 2.0])
    ρ = ArrayType([1770.0, 2170.0])
    κ = ArrayType([0.53, 1.12])

    kernel! = test_aerosol_activation_kernel!(dev, work_groups)
    event =
        kernel!(prs, output, r, stdev, N, ϵ, ϕ, M, ν, ρ, κ, ndrange = ndrange)
    wait(dev, event)

    # test if all aerosol activation output is positive
    @test all(Array(output)[:, :] .>= FT(0))
    # test if higroscopicity parameter is the same for κ and B modes
    @test all(isapprox(Array(output)[1, :], Array(output)[2, :], rtol = 0.3))
    # test if the number and mass activated are the same for κ and B modes
    @test all(isapprox(Array(output)[3, :], Array(output)[4, :], rtol = 1e-5))
    @test all(isapprox(Array(output)[5, :], Array(output)[6, :], rtol = 1e-5))
end

@testset "0-moment microphysics kernels" begin
    FT = Float32
    data_length = 3
    output = ArrayType(Array{FT}(undef, 2, data_length))
    fill!(output, -44.0)

    dev = device(ArrayType)
    work_groups = (1,)
    ndrange = (data_length,)

    liquid_frac = ArrayType([0.0, 0.5, 1.0])
    qt = ArrayType([13e-3, 13e-3, 13e-3])
    qc = ArrayType([3e-3, 4e-3, 5e-3])

    kernel! = test_0_moment_micro_kernel!(dev, work_groups)
    event = kernel!(prs, output, liquid_frac, qc, qt, ndrange = ndrange)
    wait(dev, event)

    # test 0-moment rain removal is callable and returns a reasonable value
    @test all(isequal(Array(output)[1, :], Array(output)[2, :]))
end

@testset "1-moment microphysics kernels" begin
    FT = Float32
    data_length = 2
    output = ArrayType(Array{FT}(undef, 7, data_length))
    fill!(output, -44.0)

    dev = device(ArrayType)
    work_groups = (1,)
    ndrange = (data_length,)

    ρ = ArrayType([1.2, 1.2])
    qt = ArrayType([0.0, 20e-3])
    qi = ArrayType([0.0, 5e-4])
    qs = ArrayType([0.0, 5e-4])
    ql = ArrayType([0.0, 5e-4])
    qr = ArrayType([0.0, 5e-4])

    kernel! = test_1_moment_micro_accretion_kernel!(dev, work_groups)
    event = kernel!(prs, output, ρ, qt, qi, qs, ql, qr, ndrange = ndrange)
    wait(dev, event)

    # test 1-moment accretion is callable and returns a reasonable value
    @test all(Array(output)[:, 1] .== FT(0))
    @test Array(output)[1, 2] ≈ 1.4150106417043544e-6
    @test Array(output)[2, 2] ≈ 2.453070979562392e-7
    @test Array(output)[3, 2] ≈ 2.453070979562392e-7
    @test Array(output)[4, 2] ≈ 1.768763302130443e-6
    @test Array(output)[5, 2] ≈ 3.085229094251214e-5
    @test Array(output)[6, 2] ≈ 2.1705865794293408e-4
    @test Array(output)[7, 2] ≈ 6.0118801860768854e-5

    data_length = 3
    output = ArrayType(Array{FT}(undef, 1, data_length))
    fill!(output, -44.0)

    dev = device(ArrayType)
    work_groups = (1,)
    ndrange = (data_length,)

    ρ = ArrayType([1.2, 1.2, 1.2])
    T = ArrayType([273.15 + 2, 273.15 + 2, 273.15 - 2])
    qs = ArrayType([1e-4, 0.0, 1e-4])

    kernel! = test_1_moment_micro_snow_melt_kernel!(dev, work_groups)
    event = kernel!(prs, output, ρ, T, qs, ndrange = ndrange)
    wait(dev, event)

    # test if 1-moment snow melt is callable and returns reasonable values
    @test Array(output)[1] ≈ 9.518235437405256e-6
    @test Array(output)[2] ≈ 0.0
    @test Array(output)[3] ≈ 0.0
end
