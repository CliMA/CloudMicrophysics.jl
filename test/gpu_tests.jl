using Test
using KernelAbstractions
using ClimaComms

# Needed for parameters
import CloudMicrophysics.Parameters as CMP
import CLIMAParameters as CP
import Thermodynamics as TD

# Modules to test
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.HetIceNucleation as CMI_het
import CloudMicrophysics.HomIceNucleation as CMI_hom
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Microphysics0M as CM0
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Nucleation as MN
import CloudMicrophysics.P3Scheme as P3

const work_groups = (1,)

ClimaComms.device() isa ClimaComms.CUDADevice || error("No GPU found")

# Set up GPU
using CUDA
const backend = CUDABackend()
CUDA.allowscalar(false)
const ArrayType = CuArray

# For debugging on the CPU
#const backend = CPU()
#const ArrayType = Array

@info "GPU Tests" backend ArrayType

@kernel function aerosol_activation_kernel!(
    ap,
    aps,
    tps,
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
    T = FT(294)      # air temperature K
    p = FT(1e5)    # air pressure Pa
    w = FT(0.5)         # vertical velocity m/s
    p_vs = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
    q_vs = 1 / (1 - TD.Parameters.molmass_ratio(tps) * (p_vs - p) / p_vs)
    # water vapor specific humidity (saturated)
    q = TD.PhasePartition(q_vs)

    args = (aps, tps, T, p, w, q)

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

        output[1, i] = AA.mean_hygroscopicity_parameter(ap, arsl_dst_B)[1]
        output[2, i] = AA.mean_hygroscopicity_parameter(ap, arsl_dst_κ)[1]

        output[3, i] = AA.total_N_activated(ap, arsl_dst_B, args...)
        output[4, i] = AA.total_N_activated(ap, arsl_dst_κ, args...)

        output[5, i] = AA.total_M_activated(ap, arsl_dst_B, args...)
        output[6, i] = AA.total_M_activated(ap, arsl_dst_κ, args...)
    end
end

@kernel function test_0_moment_micro_kernel!(
    p0m,
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

        output[1, i] = CM0.remove_precipitation(p0m, q)
        output[2, i] = -max(0, ql + qi - p0m.qc_0) / p0m.τ_precip
    end
end

@kernel function test_1_moment_micro_accretion_kernel!(
    liquid,
    rain,
    ice,
    snow,
    ce,
    blk1mvel,
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
        output[1, i] =
            CM1.accretion(liquid, rain, blk1mvel.rain, ce, ql[i], qr[i], ρ[i])
        output[2, i] =
            CM1.accretion(ice, snow, blk1mvel.snow, ce, qi[i], qs[i], ρ[i])
        output[3, i] =
            CM1.accretion(liquid, snow, blk1mvel.snow, ce, ql[i], qs[i], ρ[i])
        output[4, i] =
            CM1.accretion(ice, rain, blk1mvel.rain, ce, qi[i], qr[i], ρ[i])
        output[5, i] = CM1.accretion_rain_sink(
            rain,
            ice,
            blk1mvel.rain,
            ce,
            qi[i],
            qr[i],
            ρ[i],
        )
        output[6, i] = CM1.accretion_snow_rain(
            snow,
            rain,
            blk1mvel.snow,
            blk1mvel.rain,
            ce,
            qs[i],
            qr[i],
            ρ[i],
        )
        output[7, i] = CM1.accretion_snow_rain(
            rain,
            snow,
            blk1mvel.rain,
            blk1mvel.snow,
            ce,
            qr[i],
            qs[i],
            ρ[i],
        )
    end
end

@kernel function test_1_moment_micro_snow_melt_kernel!(
    snow,
    blk1mvel,
    aps,
    tps,
    output::AbstractArray{FT},
    ρ,
    T,
    qs,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] =
            CM1.snow_melt(snow, blk1mvel.snow, aps, tps, qs[i], ρ[i], T[i])
    end
end

@kernel function test_2_moment_acnv_kernel!(
    KK2000,
    B1994,
    TC1980,
    LD2004,
    VarTSc,
    output::AbstractArray{FT},
    ql,
    ρ,
    Nd,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1, i] = CM2.conv_q_liq_to_q_rai(VarTSc, ql[i], ρ[i], N_d = Nd[i])
        output[2, i] = CM2.conv_q_liq_to_q_rai(LD2004, ql[i], ρ[i], N_d = Nd[i])
        output[3, i] = CM2.conv_q_liq_to_q_rai(TC1980, ql[i], ρ[i], N_d = Nd[i])
        output[4, i] = CM2.conv_q_liq_to_q_rai(B1994, ql[i], ρ[i], N_d = Nd[i])
        output[5, i] = CM2.conv_q_liq_to_q_rai(KK2000, ql[i], ρ[i], N_d = Nd[i])
    end
end

@kernel function test_2_moment_accr_kernel!(
    KK2000,
    B1994,
    TC1980,
    output::AbstractArray{FT},
    ql,
    qr,
    ρ,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1, i] = CM2.accretion(KK2000, ql[i], qr[i], ρ[i])
        output[2, i] = CM2.accretion(B1994, ql[i], qr[i], ρ[i])
        output[3, i] = CM2.accretion(TC1980, ql[i], qr[i])
    end
end

@kernel function test_2_moment_SB2006_kernel!(
    aps,
    tps,
    SB2006,
    SB2006Vel,
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
                SB2006,
                ql[i],
                qr[i],
                ρ[i],
                Nl[i],
            ).au.dq_liq_dt
        output[2, i] =
            CM2.autoconversion_and_liquid_self_collection(
                SB2006,
                ql[i],
                qr[i],
                ρ[i],
                Nl[i],
            ).au.dN_liq_dt
        output[3, i] =
            CM2.autoconversion_and_liquid_self_collection(
                SB2006,
                ql[i],
                qr[i],
                ρ[i],
                Nl[i],
            ).au.dq_rai_dt
        output[4, i] =
            CM2.autoconversion_and_liquid_self_collection(
                SB2006,
                ql[i],
                qr[i],
                ρ[i],
                Nl[i],
            ).au.dN_rai_dt

        output[5, i] =
            CM2.autoconversion_and_liquid_self_collection(
                SB2006,
                ql[i],
                qr[i],
                ρ[i],
                Nl[i],
            ).sc
        output[6, i] =
            CM2.accretion(SB2006, ql[i], qr[i], ρ[i], Nl[i]).dq_liq_dt
        output[7, i] =
            CM2.accretion(SB2006, ql[i], qr[i], ρ[i], Nl[i]).dN_liq_dt
        output[8, i] =
            CM2.accretion(SB2006, ql[i], qr[i], ρ[i], Nl[i]).dq_rai_dt
        output[9, i] =
            CM2.accretion(SB2006, ql[i], qr[i], ρ[i], Nl[i]).dN_rai_dt
        output[10, i] =
            CM2.rain_self_collection_and_breakup(SB2006, qr[i], ρ[i], Nr[i]).sc
        output[11, i] =
            CM2.rain_self_collection_and_breakup(SB2006, qr[i], ρ[i], Nr[i]).br
        output[12, i] =
            CM2.rain_terminal_velocity(SB2006, SB2006Vel, qr[i], ρ[i], Nr[i])[1]
        output[13, i] =
            CM2.rain_terminal_velocity(SB2006, SB2006Vel, qr[i], ρ[i], Nr[i])[2]
        output[14, i] =
            CM2.rain_evaporation(SB2006, aps, tps, q, qr[i], ρ[i], Nr[i], T[i])[1]
        output[15, i] =
            CM2.rain_evaporation(SB2006, aps, tps, q, qr[i], ρ[i], Nr[i], T[i])[2]
    end
end

@kernel function h2so4_nucleation_kernel!(
    h2so4_nuc,
    output::AbstractArray{FT},
    h2so4_conc,
    nh3_conc,
    negative_ion_conc,
    temp,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = sum(
            MN.h2so4_nucleation_rate(
                h2so4_conc[i],
                nh3_conc[i],
                negative_ion_conc[i],
                temp[i],
                h2so4_nuc,
            ),
        )
    end
end

@kernel function organic_nucleation_kernel!(
    org_nuc,
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
        output[i] = MN.organic_nucleation_rate(
            negative_ion_conc[i],
            monoterpene_conc[i],
            O3_conc[i],
            OH_conc[i],
            temp[i],
            condensation_sink[i],
            org_nuc,
        )
    end
end

@kernel function organic_and_h2so4_nucleation_kernel!(
    mix_nuc,
    output::AbstractArray{FT},
    h2so4_conc,
    monoterpene_conc,
    OH_conc,
    temp,
    condensation_sink,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = MN.organic_and_h2so4_nucleation_rate(
            h2so4_conc[i],
            monoterpene_conc[i],
            OH_conc[i],
            temp[i],
            condensation_sink[i],
            mix_nuc,
        )
    end
end

@kernel function apparent_nucleation_rate_kernel!(
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
        output[i] = MN.apparent_nucleation_rate(
            output_diam[i],
            nucleation_rate[i],
            condensation_growth_rate[i],
            coag_sink[i],
            coag_sink_input_diam[i],
            input_diam[i],
        )
    end
end

@kernel function Common_H2SO4_soln_saturation_vapor_pressure_kernel!(
    H2SO4_prs,
    output::AbstractArray{FT},
    x_sulph,
    T,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] =
            CO.H2SO4_soln_saturation_vapor_pressure(H2SO4_prs, x_sulph[i], T[i])
    end
end

@kernel function Common_a_w_xT_kernel!(
    H2SO4_prs,
    tps,
    output::AbstractArray{FT},
    x_sulph,
    T,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = CO.a_w_xT(H2SO4_prs, tps, x_sulph[i], T[i])
    end
end

@kernel function Common_a_w_eT_kernel!(
    tps,
    output::AbstractArray{FT},
    e,
    T,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = CO.a_w_eT(tps, e[i], T[i])
    end
end

@kernel function Common_a_w_ice_kernel!(
    tps,
    output::AbstractArray{FT},
    T,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[i] = CO.a_w_ice(tps, T[i])
    end
end

@kernel function IceNucleation_deposition_J_kernel!(
    output::AbstractArray{FT},
    kaolinite,
    feldspar,
    ferrihydrite,
    Delta_a_w,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1] = CMI_het.deposition_J(kaolinite, Delta_a_w[1])
        output[2] = CMI_het.deposition_J(feldspar, Delta_a_w[2])
        output[3] = CMI_het.deposition_J(ferrihydrite, Delta_a_w[3])
    end
end

@kernel function IceNucleation_ABIFM_J_kernel!(
    output::AbstractArray{FT},
    kaolinite,
    illite,
    Delta_a_w,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1] = CMI_het.ABIFM_J(kaolinite, Delta_a_w[1])
        output[2] = CMI_het.ABIFM_J(illite, Delta_a_w[2])
    end
end

@kernel function IceNucleation_homogeneous_J_kernel!(
    ip,
    output::AbstractArray{FT},
    Delta_a_w,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1] = CMI_hom.homogeneous_J(ip.homogeneous, Delta_a_w[1])
    end
end

@kernel function P3_scheme_kernel!(
    p3,
    output::AbstractArray{FT},
    F_r,
    ρ_r,
) where {FT}

    i = @index(Group, Linear)

    @inbounds begin
        output[1, i] = P3.thresholds(p3, ρ_r[i], F_r[i])[1]
        output[2, i] = P3.thresholds(p3, ρ_r[i], F_r[i])[2]
    end
end

"""
    setup_output(dims, FT)
Helper function for GPU tests. Allocates an array of type `FT` with dimensions
`dims`. The second element of `dims` is the `data_length`.
"""
function setup_output(dims, FT)
    output = allocate(backend, FT, dims...)
    return (; output, ndrange = (dims[2],))
end

function test_gpu(FT)

    # thermodynamics and air properties
    aps = CMP.AirProperties(FT)
    tps = TD.Parameters.ThermodynamicsParameters(FT)

    # aerosol activation
    ap = CMP.AerosolActivationParameters(FT)

    # 0-moment microphysocs
    p0m = CMP.Parameters0M(FT)

    # 1-momeny microphysics
    liquid = CMP.CloudLiquid(FT)
    ice = CMP.CloudIce(FT)
    rain = CMP.Rain(FT)
    snow = CMP.Snow(FT)
    ce = CMP.CollisionEff(FT)

    # Terminal velocity
    blk1mvel = CMP.Blk1MVelType(FT)
    SB2006Vel = CMP.SB2006VelType(FT)

    # 2-moment microphysics
    SB2006 = CMP.SB2006(FT)
    KK2000 = CMP.KK2000(FT)
    B1994 = CMP.B1994(FT)
    TC1980 = CMP.TC1980(FT)
    LD2004 = CMP.LD2004(FT)
    VarTSc = CMP.VarTimescaleAcnv(FT)

    # p3 microphysics
    p3 = CMP.ParametersP3(FT)

    # aerosol nucleation
    mix_nuc = CMP.MixedNucleationParameters(FT)
    h2so4_nuc = CMP.H2S04NucleationParameters(FT)
    org_nuc = CMP.OrganicNucleationParameters(FT)

    # ice nucleation
    H2SO4_prs = CMP.H2SO4SolutionParameters(FT)
    illite = CMP.Illite(FT)
    kaolinite = CMP.Kaolinite(FT)
    feldspar = CMP.Feldspar(FT)
    ferrihydrite = CMP.Ferrihydrite(FT)
    ip = CMP.IceNucleationParameters(FT)

    @testset "Aerosol activation kernels" begin
        dims = (6, 2)
        (; output, ndrange) = setup_output(dims, FT)

        r = ArrayType([FT(0.243 * 1e-6), FT(1.5 * 1e-6)])
        stdev = ArrayType([FT(1.4), FT(2.1)])
        N = ArrayType([FT(100 * 1e6), FT(1 * 1e6)])
        ϵ = ArrayType([FT(1), FT(1)])
        ϕ = ArrayType([FT(1), FT(0.9)])
        M = ArrayType([FT(0.132), FT(0.058443)])
        ν = ArrayType([FT(3), FT(2)])
        ρ = ArrayType([FT(1770), FT(2170)])
        κ = ArrayType([FT(0.53), FT(1.12)])

        kernel! = aerosol_activation_kernel!(backend, work_groups)
        kernel!(ap, aps, tps, output, r, stdev, N, ϵ, ϕ, M, ν, ρ, κ, ; ndrange)

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
        dims = (3, 2)
        (; output, ndrange) = setup_output(dims, FT)

        liquid_frac = ArrayType([FT(0), FT(0.5), FT(1)])
        qt = ArrayType([FT(13e-3), FT(13e-3), FT(13e-3)])
        qc = ArrayType([FT(3e-3), FT(4e-3), FT(5e-3)])

        kernel! = test_0_moment_micro_kernel!(backend, work_groups)
        kernel!(p0m, output, liquid_frac, qc, qt, ; ndrange)

        # test 0-moment rain removal is callable and returns a reasonable value
        @test all(isequal(Array(output)[1, :], Array(output)[2, :]))
    end

    @testset "1-moment microphysics kernels" begin
        dims = (7, 2)
        (; output, ndrange) = setup_output(dims, FT)

        ρ = ArrayType([FT(1.2), FT(1.2)])
        qt = ArrayType([FT(0), FT(20e-3)])
        qi = ArrayType([FT(0), FT(5e-4)])
        qs = ArrayType([FT(0), FT(5e-4)])
        ql = ArrayType([FT(0), FT(5e-4)])
        qr = ArrayType([FT(0), FT(5e-4)])

        kernel! = test_1_moment_micro_accretion_kernel!(backend, work_groups)
        kernel!(
            liquid,
            rain,
            ice,
            snow,
            ce,
            blk1mvel,
            output,
            ρ,
            qt,
            qi,
            qs,
            ql,
            qr;
            ndrange,
        )

        # test 1-moment accretion is callable and returns a reasonable value
        @test all(Array(output)[:, 1] .== FT(0))
        @test Array(output)[1, 2] ≈ FT(1.4150106417043544e-6)
        @test Array(output)[2, 2] ≈ FT(2.453070979562392e-7)
        @test Array(output)[3, 2] ≈ FT(2.453070979562392e-7)
        @test Array(output)[4, 2] ≈ FT(1.768763302130443e-6)
        @test Array(output)[5, 2] ≈ FT(3.085229094251214e-5)
        @test Array(output)[6, 2] ≈ FT(2.1705865794293408e-4)
        @test Array(output)[7, 2] ≈ FT(6.0118801860768854e-5)

        dims = (1, 3)
        (; output, ndrange) = setup_output(dims, FT)

        ρ = ArrayType([FT(1.2), FT(1.2), FT(1.2)])
        T = ArrayType([FT(273.15 + 2), FT(273.15 + 2), FT(273.15 - 2)])
        qs = ArrayType([FT(1e-4), FT(0), FT(1e-4)])

        kernel! = test_1_moment_micro_snow_melt_kernel!(backend, work_groups)
        kernel!(snow, blk1mvel, aps, tps, output, ρ, T, qs, ; ndrange)

        # test if 1-moment snow melt is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(9.518235437405256e-6)
        @test Array(output)[2] ≈ FT(0)
        @test Array(output)[3] ≈ FT(0)
    end

    @testset "2-moment microphysics kernels" begin

        dims = (5, 1)
        (; output, ndrange) = setup_output(dims, FT)

        ql = ArrayType([FT(2e-3)])
        ρ = ArrayType([FT(1.2)])
        Nd = ArrayType([FT(1e8)])

        kernel! = test_2_moment_acnv_kernel!(backend, work_groups)
        kernel!(
            KK2000,
            B1994,
            TC1980,
            LD2004,
            VarTSc,
            output,
            ql,
            ρ,
            Nd,
            ndrange = ndrange,
        )

        @test Array(output)[1] ≈ FT(2e-6)
        @test Array(output)[2] ≈ FT(1.6963072465911614e-6)
        @test Array(output)[3] ≈ FT(3.5482867084128596e-6)
        @test Array(output)[4] ≈ FT(9.825462758968215e-7)
        @test Array(output)[5] ≈ FT(5.855332513368727e-8)

        dims = (3, 1)
        (; output, ndrange) = setup_output(dims, FT)

        ql = ArrayType([FT(2e-3)])
        qr = ArrayType([FT(5e-4)])
        ρ = ArrayType([FT(1.2)])

        kernel! = test_2_moment_accr_kernel!(backend, work_groups)
        kernel!(KK2000, B1994, TC1980, output, ql, qr, ρ, ndrange = ndrange)

        @test isapprox(Array(output)[1], FT(6.6548664e-6), rtol = 1e-6)
        @test Array(output)[2] ≈ FT(7.2e-6)
        @test Array(output)[3] ≈ FT(4.7e-6)

        dims = (15, 1)
        (; output, ndrange) = setup_output(dims, FT)

        T = ArrayType([FT(290)])
        qt = ArrayType([FT(7e-3)])
        ql = ArrayType([FT(2e-3)])
        qr = ArrayType([FT(5e-4)])
        ρ = ArrayType([FT(1.2)])
        Nl = ArrayType([FT(1e8)])
        Nr = ArrayType([FT(1e7)])

        kernel! = test_2_moment_SB2006_kernel!(backend, work_groups)
        kernel!(
            aps,
            tps,
            SB2006,
            SB2006Vel,
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
        @test isapprox(Array(output)[14], FT(-243259.75126), rtol = 1e-6)
        @test isapprox(Array(output)[15], FT(-0.0034601581), rtol = 1e-6)
    end

    @testset "Common Kernels" begin
        dims = (1, 1)
        (; output, ndrange) = setup_output(dims, FT)

        T = ArrayType([FT(230)])
        x_sulph = ArrayType([FT(0.1)])

        kernel! = Common_H2SO4_soln_saturation_vapor_pressure_kernel!(
            backend,
            work_groups,
        )
        kernel!(H2SO4_prs, output, x_sulph, T; ndrange)

        # test H2SO4_soln_saturation_vapor_pressure is callable and returns a reasonable value
        @test Array(output)[1] ≈ FT(12.685507586924)

        dims = (1, 1)
        (; output, ndrange) = setup_output(dims, FT)

        T = ArrayType([FT(230)])
        x_sulph = ArrayType([FT(0.1)])

        kernel! = Common_a_w_xT_kernel!(backend, work_groups)
        kernel!(H2SO4_prs, tps, output, x_sulph, T; ndrange)

        # test if a_w_xT is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(0.92824538441)

        dims = (1, 1)
        (; output, ndrange) = setup_output(dims, FT)

        T = ArrayType([FT(282)])
        e = ArrayType([FT(1001)])

        kernel! = Common_a_w_eT_kernel!(backend, work_groups)
        kernel!(tps, output, e, T; ndrange)

        # test if a_w_eT is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(0.880978146)

        dims = (1, 1)
        (; output, ndrange) = setup_output(dims, FT)

        T = ArrayType([FT(230)])

        kernel! = Common_a_w_ice_kernel!(backend, work_groups)
        kernel!(tps, output, T; ndrange)

        # test if a_w_ice is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(0.653191723)
    end

    @testset "Ice Nucleation kernels" begin
        dims = (1, 3)
        (; output, ndrange) = setup_output(dims, FT)

        Delta_a_w = ArrayType([FT(0.16), FT(0.15), FT(0.15)])

        kernel! = IceNucleation_deposition_J_kernel!(backend, work_groups)
        kernel!(output, kaolinite, feldspar, ferrihydrite, Delta_a_w; ndrange)

        # test if deposition_J is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(1.5390757663075784e6)
        @test Array(output)[2] ≈ FT(5.693312205851678e6)
        @test Array(output)[3] ≈ FT(802555.3607426438)

        dims = (1, 2)
        (; output, ndrange) = setup_output(dims, FT)

        Delta_a_w = ArrayType([FT(0.16), FT(0.15)])

        kernel! = IceNucleation_ABIFM_J_kernel!(backend, work_groups)
        kernel!(output, kaolinite, illite, Delta_a_w; ndrange)

        # test if ABIFM_J is callable and returns reasonable values
        @test Array(output)[1] ≈ FT(153.65772539109)
        @test Array(output)[2] ≈ FT(31.870032033791)

        dims = (1, 1)
        (; output, ndrange) = setup_output(dims, FT)

        T = ArrayType([FT(220)])
        x_sulph = ArrayType([FT(0.15)])
        Delta_a_w = ArrayType([FT(0.2907389666103033)])

        kernel! = IceNucleation_homogeneous_J_kernel!(backend, work_groups)
        kernel!(ip, output, Delta_a_w; ndrange)

        # test homogeneous_J is callable and returns a reasonable value
        @test Array(output)[1] ≈ FT(2.66194650334444e12)
    end

    @testset "Modal nucleation kernels" begin
        dims = (1, 2)
        (; output, ndrange) = setup_output(dims, FT)

        # h2so4 nucleation
        h2so4_conc = ArrayType([FT(1e12), FT(1e12)])
        nh3_conc = ArrayType([FT(1), FT(1)])
        negative_ion_conc = ArrayType([FT(1), FT(1)])
        temp = ArrayType([FT(208), FT(208)])

        kernel! = h2so4_nucleation_kernel!(backend, work_groups)
        kernel!(
            h2so4_nuc,
            output,
            h2so4_conc,
            nh3_conc,
            negative_ion_conc,
            temp;
            ndrange,
        )

        @test all(Array(output) .> FT(0))

        # Organic nucleation

        negative_ion_conc = ArrayType([FT(0.0), FT(0.0)])
        monoterpene_conc = ArrayType([FT(1e24), FT(1e24)])
        O3_conc = ArrayType([FT(1e24), FT(1e24)])
        OH_conc = ArrayType([FT(1e24), FT(1e24)])
        temp = ArrayType([FT(300), FT(300)])
        condensation_sink = ArrayType([FT(1), FT(1)])

        kernel! = organic_nucleation_kernel!(backend, work_groups)
        kernel!(
            org_nuc,
            output,
            negative_ion_conc,
            monoterpene_conc,
            O3_conc,
            OH_conc,
            temp,
            condensation_sink;
            ndrange,
        )

        @test all(Array(output) .> FT(0))

        # Organic and h2so4 nucleation

        h2so4_conc = ArrayType([FT(2.6e6), FT(2.6e6)])
        monoterpene_conc = ArrayType([FT(1), FT(1)])
        OH_conc = ArrayType([FT(1), FT(1)])
        temp = ArrayType([FT(300), FT(300)])
        condensation_sink = ArrayType([FT(1), FT(1)])

        kernel! = organic_and_h2so4_nucleation_kernel!(backend, work_groups)
        kernel!(
            mix_nuc,
            output,
            h2so4_conc,
            monoterpene_conc,
            OH_conc,
            temp,
            condensation_sink;
            ndrange,
        )

        @test all(Array(output) .> FT(0))

        # Apparent nucleation rate

        output_diam = ArrayType([FT(1e-9), FT(1e-9)])
        nucleation_rate = ArrayType([FT(1e6), FT(1e6)])
        condensation_growth_rate = ArrayType([FT(1e6), FT(1e6)])
        coag_sink = ArrayType([FT(1e6), FT(1e6)])
        coag_sink_input_diam = ArrayType([FT(1e-9), FT(1e-9)])
        input_diam = ArrayType([FT(1e-9), FT(1e-9)])

        kernel! = apparent_nucleation_rate_kernel!(backend, work_groups)
        kernel!(
            output,
            output_diam,
            nucleation_rate,
            condensation_growth_rate,
            coag_sink,
            coag_sink_input_diam,
            input_diam;
            ndrange,
        )

        @test all(Array(output) .> FT(0))
    end
    @testset "P3 scheme kernels" begin
        dims = (2, 2)
        (; output, ndrange) = setup_output(dims, FT)

        F_r = ArrayType([FT(0.5), FT(0.95)])
        ρ_r = ArrayType([FT(400), FT(800)])

        kernel! = P3_scheme_kernel!(backend, work_groups)
        kernel!(p3, output, F_r, ρ_r; ndrange)

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
end

println("Testing Float64")
test_gpu(Float64)

println("Testing Float32")
test_gpu(Float32)
