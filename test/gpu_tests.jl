import Test as TT
using KernelAbstractions
using ClimaComms
ClimaComms.@import_required_backends

# Needed for parameters
import ClimaParams as CP

# Modules to test
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.HetIceNucleation as CMI_het
import CloudMicrophysics.HomIceNucleation as CMI_hom
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Microphysics0M as CM0
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.MicrophysicsNonEq as CMN
import CloudMicrophysics.Nucleation as MN
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT

work_groups = 2

ClimaComms.device() isa ClimaComms.CUDADevice || error("No GPU found")

# Set up GPU
using CUDA
backend = CUDABackend()
CUDA.allowscalar(false)
ArrayType = CuArray

# For debugging on the CPU
# backend = CPU()
# ArrayType = Array

@info "GPU Tests" backend ArrayType

@kernel inbounds = true function aerosol_activation_kernel!(
    ap, aps, tps, output, r, stdev, N, ϵ, ϕ, M, ν, ρ, κ,
)
    i = @index(Global, Linear)
    FT = eltype(ap)
    # atmospheric conditions (taken from aerosol activation tests)
    T = FT(294)      # air temperature K
    p = FT(1e5)    # air pressure Pa
    w = FT(0.5)         # vertical velocity m/s
    p_vs = TDI.saturation_vapor_pressure_over_liquid(tps, T)
    q_vs = 1 / (1 - 1 / TDI.Rd_over_Rv(tps) * (p_vs - p) / p_vs)
    q_liq = FT(0)
    q_ice = FT(0)

    args = (aps, tps, T, p, w, q_vs, q_liq, q_ice)

    mode_B = AM.Mode_B(
        r[i], stdev[i], N[i], (FT(1.0),), (ϵ[i],), (ϕ[i],), (M[i],), (ν[i],), (ρ[i],),
    )
    mode_κ = AM.Mode_κ(r[i], stdev[i], N[i], (FT(1.0),), (FT(1.0),), (M[i],), (κ[i],))

    arsl_dst_B = AM.AerosolDistribution((mode_B,))
    arsl_dst_κ = AM.AerosolDistribution((mode_κ,))

    mean_B = AA.mean_hygroscopicity_parameter(ap, arsl_dst_B)[1]
    mean_κ = AA.mean_hygroscopicity_parameter(ap, arsl_dst_κ)[1]

    N_act_B = AA.total_N_activated(ap, arsl_dst_B, args...)
    N_act_κ = AA.total_N_activated(ap, arsl_dst_κ, args...)

    M_act_B = AA.total_M_activated(ap, arsl_dst_B, args...)
    M_act_κ = AA.total_M_activated(ap, arsl_dst_κ, args...)

    output[i] = (; mean_B, N_act_B, M_act_B, mean_κ, N_act_κ, M_act_κ)
end

@kernel inbounds = true function test_noneq_micro_kernel!(
    lcl, icl, tps, output, ρ, T, qᵥ_sl, qᵢ, qᵢ_s,
)
    i = @index(Global, Linear)
    FT = eltype(tps)
    q_lcl = q_icl = q_rai = q_sno = FT(0) # set to zero in this test
    S_cond_MM2015 = CMN.conv_q_vap_to_q_lcl_MM2015(
        tps, qᵥ_sl[i], q_lcl, q_icl, q_rai, q_sno, ρ[i], T[i], lcl.τ_relax,
    )
    S_cond = CMN.conv_q_vap_to_q_lcl_icl(icl, qᵢ_s[i], qᵢ[i])
    output[i] = (; S_cond_MM2015, S_cond)
end

@kernel inbounds = true function test_chen2022_terminal_velocity_kernel!(
    lcl, icl, rain, snow, Ch2022, STVel, output, ρₐ, qₗ, qᵢ, qᵣ, qₛ,
)
    i = @index(Global, Linear)
    # CloudLiquid supports StokesRegimeVelType, not Chen2022VelTypeRain
    lcl_stokes = CMN.terminal_velocity(lcl, STVel, ρₐ[i], qₗ[i])
    lci_chen2022 = CMN.terminal_velocity(icl, Ch2022.small_ice, ρₐ[i], qᵢ[i])
    rai_chen2022 = CM1.terminal_velocity(rain, Ch2022.rain, ρₐ[i], qᵣ[i])
    sno_chen2022 = CM1.terminal_velocity(snow, Ch2022.large_ice, ρₐ[i], qₛ[i])
    output[i] = (; lcl_stokes, lci_chen2022, rai_chen2022, sno_chen2022)
end

@kernel inbounds = true function test_0_moment_micro_kernel!(p0m, liquid_frac, output, qc)
    i = @index(Global, Linear)
    ql = qc[i] * liquid_frac[i]
    qi = (1 - liquid_frac[i]) * qc[i]
    S_pr = CM0.remove_precipitation(p0m, ql, qi)
    S_pr_ref = -max(0, ql + qi - p0m.qc_0) / p0m.τ_precip
    output[i] = (; S_pr, S_pr_ref)
end

@kernel inbounds = true function test_1_moment_micro_acnv_kernel!(output, acnv1M, ql)
    i = @index(Global, Linear)
    output[i] = CM1.conv_q_lcl_to_q_rai(acnv1M, ql[i], false)
end

@kernel inbounds = true function test_1_moment_micro_accretion_kernel!(
    lcl, rain, icl, snow, ce, blk1mvel, output, ρ, qi, qs, ql, qr,
)
    i = @index(Global, Linear)
    liq_rai = CM1.accretion(lcl, rain, blk1mvel.rain, ce, ql[i], qr[i], ρ[i])
    ice_sno = CM1.accretion(icl, snow, blk1mvel.snow, ce, qi[i], qs[i], ρ[i])
    liq_sno = CM1.accretion(lcl, snow, blk1mvel.snow, ce, ql[i], qs[i], ρ[i])
    ice_rai = CM1.accretion(icl, rain, blk1mvel.rain, ce, qi[i], qr[i], ρ[i])
    rai_sink = CM1.accretion_rain_sink(rain, icl, blk1mvel.rain, ce, qi[i], qr[i], ρ[i])
    sno_rai = CM1.accretion_snow_rain(
        snow, rain, blk1mvel.snow, blk1mvel.rain, ce, qs[i], qr[i], ρ[i],
    )
    rai_sno = CM1.accretion_snow_rain(
        rain, snow, blk1mvel.rain, blk1mvel.snow, ce, qr[i], qs[i], ρ[i],
    )
    output[i] = (; liq_rai, ice_sno, liq_sno, ice_rai, rai_sink, sno_rai, rai_sno)
end

@kernel inbounds = true function test_1_moment_micro_snow_melt_kernel!(
    snow, blk1mvel, aps, tps, output, ρ, T, qs)
    i = @index(Global, Linear)
    output[i] = CM1.snow_melt(snow, blk1mvel.snow, aps, tps, qs[i], ρ[i], T[i])
end

@kernel inbounds = true function test_2_moment_acnv_kernel!(
    KK2000, B1994, TC1980, LD2004, VarTSc, output, ql, ρ, Nd,
)
    i = @index(Global, Linear)
    S_VarTSc = CM2.conv_q_lcl_to_q_rai(VarTSc, ql[i], ρ[i], Nd[i])
    S_LD2004 = CM2.conv_q_lcl_to_q_rai(LD2004, ql[i], ρ[i], Nd[i])
    S_TC1980 = CM2.conv_q_lcl_to_q_rai(TC1980, ql[i], ρ[i], Nd[i])
    S_B1994 = CM2.conv_q_lcl_to_q_rai(B1994, ql[i], ρ[i], Nd[i])
    S_KK2000 = CM2.conv_q_lcl_to_q_rai(KK2000, ql[i], ρ[i], Nd[i])
    output[i] = (; S_VarTSc, S_LD2004, S_TC1980, S_B1994, S_KK2000)
end

@kernel inbounds = true function test_2_moment_accr_kernel!(
    KK2000, B1994, TC1980, output, ql, qr, ρ,
)
    i = @index(Global, Linear)
    S_KK2000 = CM2.accretion(KK2000, ql[i], qr[i], ρ[i])
    S_B1994 = CM2.accretion(B1994, ql[i], qr[i], ρ[i])
    S_TC1980 = CM2.accretion(TC1980, ql[i], qr[i])
    output[i] = (; S_KK2000, S_B1994, S_TC1980)
end

function SB2006_2M_kernel(aps, tps, SB2006, SB2006Vel, qt, ql, qr, Nl, Nr, ρ, T)
    FT = eltype(aps)
    lcl_aconv_scoll =
        CM2.autoconversion_and_cloud_liquid_self_collection(SB2006, ql, qr, ρ, Nl)
    accr = CM2.accretion(SB2006, ql, qr, ρ, Nl)
    rain_scoll = CM2.rain_self_collection_and_breakup(SB2006, qr, ρ, Nr)
    rain_vel = CM2.rain_terminal_velocity(SB2006, SB2006Vel, qr, ρ, Nr)
    rain_evap = CM2.rain_evaporation(SB2006, aps, tps, qt, ql, FT(0), qr, FT(0), ρ, Nr, T)
    num_incr_mass_limit =
        CM2.number_increase_for_mass_limit(SB2006.numadj, SB2006.pdf_r.xr_max, qr, ρ, Nr)
    num_decr_mass_limit =
        CM2.number_decrease_for_mass_limit(SB2006.numadj, SB2006.pdf_r.xr_min, qr, ρ, Nr)
    return (;
        lcl_aconv_scoll, accr, rain_scoll, rain_vel,
        rain_evap, num_incr_mass_limit, num_decr_mass_limit,
    )
end

@kernel inbounds = true function test_2_moment_SB2006_kernel!(
    aps, tps, SB2006, SB2006Vel, output, qt, ql, qr, Nl, Nr, ρ, T,
)
    i = @index(Global, Linear)
    output[i] = SB2006_2M_kernel(
        aps, tps, SB2006, SB2006Vel, qt[i], ql[i], qr[i], Nl[i], Nr[i], ρ[i], T[i],
    )
end

@kernel inbounds = true function h2so4_nucleation_kernel!(
    h2so4_nuc, output, h2so4_conc, nh3_conc, negative_ion_conc, temp,
)
    i = @index(Global, Linear)
    output[i] = sum(MN.h2so4_nucleation_rate(
        h2so4_conc[i], nh3_conc[i], negative_ion_conc[i], temp[i], h2so4_nuc,
    ))
end

@kernel inbounds = true function organic_nucleation_kernel!(
    org_nuc, output, negative_ion_conc, monoterpene_conc, O3_conc, OH_conc,
    temp, condensation_sink,
)
    i = @index(Global, Linear)
    output[i] = MN.organic_nucleation_rate(
        negative_ion_conc[i], monoterpene_conc[i], O3_conc[i], OH_conc[i],
        temp[i], condensation_sink[i], org_nuc,
    )
end

@kernel inbounds = true function organic_and_h2so4_nucleation_kernel!(
    mix_nuc, output, h2so4_conc, monoterpene_conc, OH_conc, temp, condensation_sink,
)
    i = @index(Global, Linear)
    output[i] = MN.organic_and_h2so4_nucleation_rate(
        h2so4_conc[i], monoterpene_conc[i], OH_conc[i],
        temp[i], condensation_sink[i], mix_nuc,
    )
end

@kernel inbounds = true function apparent_nucleation_rate_kernel!(
    output, output_diam, nucleation_rate, condensation_growth_rate, coag_sink,
    coag_sink_input_diam, input_diam,
)
    i = @index(Global, Linear)
    output[i] = MN.apparent_nucleation_rate(
        output_diam[i], nucleation_rate[i], condensation_growth_rate[i], coag_sink[i],
        coag_sink_input_diam[i], input_diam[i],
    )
end

@kernel inbounds = true function Common_H2SO4_kernel!(H2SO4_prs, tps, output, x_sulph, T)
    i = @index(Global, Linear)
    H2SO4_soln_qv_sat = CO.H2SO4_soln_saturation_vapor_pressure(H2SO4_prs, x_sulph[i], T[i])
    H2SO4_soln_a_w = CO.a_w_xT(H2SO4_prs, tps, x_sulph[i], T[i])
    output[i] = (; H2SO4_soln_qv_sat, H2SO4_soln_a_w)
end

@kernel inbounds = true function Common_a_w_eT_kernel!(tps, output, e, T)
    i = @index(Global, Linear)
    output[i] = CO.a_w_eT(tps, e[i], T[i])
end

@kernel inbounds = true function Common_a_w_ice_kernel!(tps, output, T)
    i = @index(Global, Linear)
    output[i] = CO.a_w_ice(tps, T[i])
end

@kernel inbounds = true function IceNucleation_dust_activated_number_fraction_kernel!(
    desert_dust, arizona_test_dust, ip, output, S_i, T,
)
    i = @index(Global, Linear)
    desert_act_N_frac =
        CMI_het.dust_activated_number_fraction(desert_dust, ip, S_i[i], T[i])
    arizona_act_N_frac =
        CMI_het.dust_activated_number_fraction(arizona_test_dust, ip, S_i[i], T[i])
    output[i] = (; desert_act_N_frac, arizona_act_N_frac)
end

@kernel inbounds = true function IceNucleation_MohlerDepositionRate_kernel!(
    desert_dust, arizona_test_dust, ip, output, S_i, T, dSi_dt, N_aero,
)
    i = @index(Global, Linear)
    args = (ip, S_i[i], T[i], dSi_dt[i], N_aero[i])
    desert_dep_rate = CMI_het.MohlerDepositionRate(desert_dust, args...)
    arizona_dep_rate = CMI_het.MohlerDepositionRate(arizona_test_dust, args...)
    output[i] = (; desert_dep_rate, arizona_dep_rate)
end

@kernel inbounds = true function IceNucleation_deposition_J_kernel!(mineral, output, Delta_a_w)
    i = @index(Global, Linear)
    output[i] = CMI_het.deposition_J(mineral, Delta_a_w[i])
end

@kernel inbounds = true function IceNucleation_ABIFM_J_kernel!(mineral, output, Delta_a_w)
    i = @index(Global, Linear)
    output[i] = CMI_het.ABIFM_J(mineral, Delta_a_w[i])
end

@kernel inbounds = true function IceNucleation_P3_deposition_N_i_kernel!(ip, output, T)
    i = @index(Global, Linear)
    output[i] = CMI_het.P3_deposition_N_i(ip, T[i])
end

@kernel inbounds = true function IceNucleation_P3_het_N_i_kernel!(ip, output, T, N_l, V_l, Δt)
    i = @index(Global, Linear)
    output[i] = CMI_het.P3_het_N_i(ip, T[i], N_l[i], V_l[i], Δt[i])
end

@kernel inbounds = true function IceNucleation_INPC_frequency_kernel!(ip, output, INPC, T)
    i = @index(Global, Linear)
    output[i] = CMI_het.INP_concentration_frequency(ip, INPC[i], T[i])
end

@kernel inbounds = true function IceNucleation_homogeneous_J_kernel!(ip, output, Delta_a_w)
    i = @index(Global, Linear)
    J_cubic = CMI_hom.homogeneous_J_cubic(ip.homogeneous, Delta_a_w[i])
    J_linear = CMI_hom.homogeneous_J_linear(ip.homogeneous, Delta_a_w[i])
    output[i] = (; J_cubic, J_linear)
end

@kernel inbounds = true function IceNucleation_homogeneous_J_linear_kernel!(
    ip, output, Delta_a_w,
)
    i = @index(Global, Linear)
    output[i] = CMI_hom.homogeneous_J_linear(ip.homogeneous, Delta_a_w[i])
end

@kernel inbounds = true function test_bulk_tendencies_0m_kernel!(
    mp, tps, output, liquid_frac, qc, T,
)
    i = @index(Global, Linear)
    ql = qc[i] * liquid_frac[i]
    qi = (1 - liquid_frac[i]) * qc[i]
    CM0M = BMT.Microphysics0Moment()
    output[i] = BMT.bulk_microphysics_tendencies(CM0M, mp, tps, T[i], ql, qi)
end

@kernel inbounds = true function test_bulk_tendencies_0m_S0_kernel!(
    mp, tps, output, liquid_frac, qc, T, q_vap_sat,
)
    i = @index(Global, Linear)
    ql = qc[i] * liquid_frac[i]
    qi = (1 - liquid_frac[i]) * qc[i]
    CM0M = BMT.Microphysics0Moment()
    output[i] = BMT.bulk_microphysics_tendencies(CM0M, mp, tps, T[i], ql, qi, q_vap_sat[i])
end

@kernel inbounds = true function test_bulk_tendencies_1m_kernel!(
    mp, tps, output, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    i = @index(Global, Linear)
    CM1M = BMT.Microphysics1Moment()
    output[i] = BMT.bulk_microphysics_tendencies(
        CM1M, mp, tps, ρ[i], T[i], q_tot[i], q_lcl[i], q_icl[i], q_rai[i], q_sno[i],
    )
end

@kernel inbounds = true function test_bulk_tendencies_2m_warm_kernel!(
    mp, tps, output, ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai,
)
    i = @index(Global, Linear)
    CM2M = BMT.Microphysics2Moment()
    output[i] = BMT.bulk_microphysics_tendencies(
        CM2M, mp, tps, ρ[i], T[i], q_tot[i], q_lcl[i], n_lcl[i], q_rai[i], n_rai[i],
    )
end

@kernel inbounds = true function test_bulk_tendencies_2m_p3_kernel!(
    mp, tps, output, ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim,
)
    i = @index(Global, Linear)
    CM2M = BMT.Microphysics2Moment()
    output[i] = BMT.bulk_microphysics_tendencies(
        CM2M, mp, tps, ρ[i], T[i], q_tot[i], q_lcl[i], n_lcl[i], q_rai[i], n_rai[i],
        q_ice[i], n_ice[i], q_rim[i], b_rim[i],
    )
end

"""
    setup_output(dims, FT)
Helper function for GPU tests. Allocates an array of type `FT` with dimensions
`dims`. The last element of `dims` is the `data_length`.

Optionally, `default_value` can be provided to
uniformly initialize the output array values.
"""
function setup_output(dims, FT, default_value = nothing)
    output = allocate(backend, FT, dims...)
    !isnothing(default_value) && fill!(output, FT(default_value))
    return (; output, ndrange = (dims[end],))
end

"""
    constant_data(value; ndrange)
Helper function for GPU tests. Allocates a vector of type `FT` with length
`ndrange` and fills it with `value`.
"""
function constant_data(value; ndrange)
    FT = eltype(value)
    return fill!(allocate(backend, FT, ndrange), value)
end

function test_gpu(FT)

    # thermodynamics and air properties
    aps = CMP.AirProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    # aerosol activation
    ap = CMP.AerosolActivationParameters(FT)

    # 0-moment microphysocs
    p0m = CMP.Parameters0M(FT)

    # 1-momeny microphysics
    lcl = CMP.CloudLiquid(FT)
    icl = CMP.CloudIce(FT)
    rain = CMP.Rain(FT)
    snow = CMP.Snow(FT)
    ce = CMP.CollisionEff(FT)

    # Terminal velocity
    blk1mvel = CMP.Blk1MVelType(FT)
    SB2006Vel = CMP.SB2006VelType(FT)
    Ch2022 = CMP.Chen2022VelType(FT)
    STVel = CMP.StokesRegimeVelType(FT)

    # 2-moment microphysics
    SB2006 = CMP.SB2006(FT)
    SB2006_no_limiters = CMP.SB2006(FT, false)
    KK2000 = CMP.KK2000(FT)
    B1994 = CMP.B1994(FT)
    TC1980 = CMP.TC1980(FT)
    LD2004 = CMP.LD2004(FT)
    VarTSc = CMP.VarTimescaleAcnv(FT)

    # Bulk microphysics parameters
    mp_0m = CMP.Microphysics0MParams(FT)
    mp_1m = CMP.Microphysics1MParams(FT)
    mp_2m_warm = CMP.Microphysics2MParams(FT; with_ice = false)
    mp_2m_p3 = CMP.Microphysics2MParams(FT; with_ice = true)

    # aerosol nucleation
    mix_nuc = CMP.MixedNucleationParameters(FT)
    h2so4_nuc = CMP.H2S04NucleationParameters(FT)
    org_nuc = CMP.OrganicNucleationParameters(FT)

    # ice nucleation
    H2SO4_prs = CMP.H2SO4SolutionParameters(FT)
    desert_dust = CMP.DesertDust(FT)
    arizona_test_dust = CMP.ArizonaTestDust(FT)
    illite = CMP.Illite(FT)
    kaolinite = CMP.Kaolinite(FT)
    feldspar = CMP.Feldspar(FT)
    ferrihydrite = CMP.Ferrihydrite(FT)
    ip = CMP.IceNucleationParameters(FT)
    ip_frostenberg = CMP.Frostenberg2023(FT)

    TT.@testset "Aerosol activation kernels" begin
        r = ArrayType([FT(0.243 * 1e-6), FT(1.5 * 1e-6)])
        stdev = ArrayType([FT(1.4), FT(2.1)])
        N = ArrayType([FT(100 * 1e6), FT(1 * 1e6)])
        ϵ = ArrayType([FT(1), FT(1)])
        ϕ = ArrayType([FT(1), FT(0.9)])
        M = ArrayType([FT(0.132), FT(0.058443)])
        ν = ArrayType([FT(3), FT(2)])
        ρ = ArrayType([FT(1770), FT(2170)])
        κ = ArrayType([FT(0.53), FT(1.12)])
        DT = @NamedTuple{
            mean_B::FT, N_act_B::FT, M_act_B::FT,
            mean_κ::FT, N_act_κ::FT, M_act_κ::FT,
        }
        (; ndrange, output) = setup_output(2, DT)

        kernel! = aerosol_activation_kernel!(backend, work_groups)
        kernel!(ap, aps, tps, output, r, stdev, N, ϵ, ϕ, M, ν, ρ, κ; ndrange)
        for nt in Array(output)
            # test if all aerosol activation output is positive
            TT.@test isapprox(nt.mean_B, nt.mean_κ, rtol = 0.3)
            # test if higroscopicity parameter is the same for κ and B modes
            TT.@test isapprox(nt.N_act_B, nt.N_act_κ, rtol = 1e-5)
            # test if the number and mass activated are the same for κ and B modes
            TT.@test isapprox(nt.M_act_B, nt.M_act_κ, rtol = 1e-5)
        end
    end

    TT.@testset "non-equilibrium microphysics kernels" begin
        DT = @NamedTuple{S_cond_MM2015::FT, S_cond::FT}
        (; ndrange, output) = setup_output(1, DT)

        ρ = ArrayType([FT(0.8)])
        T = ArrayType([FT(263)])
        qᵥ_sl = ArrayType([FT(0.0035)])
        qᵢ = ArrayType([FT(0.003)])
        qᵢ_s = ArrayType([FT(0.002)])

        kernel! = test_noneq_micro_kernel!(backend, work_groups)
        kernel!(lcl, icl, tps, output, ρ, T, qᵥ_sl, qᵢ, qᵢ_s; ndrange)
        (; S_cond_MM2015, S_cond) = Array(output)[1]
        # test that nonequilibrium cloud formation is callable and returns a reasonable value
        TT.@test S_cond_MM2015 ≈ FT(3.76347635339803e-5)
        TT.@test S_cond ≈ FT(-1e-4)
    end

    TT.@testset "Chen 2022 terminal velocity kernels" begin
        DT = @NamedTuple{
            lcl_stokes::FT, lci_chen2022::FT,
            rai_chen2022::FT, sno_chen2022::FT,
        }
        (; ndrange, output) = setup_output(1, DT)

        ρ = ArrayType([FT(0.95)])
        qₗ = ArrayType([FT(0.004)])
        qᵢ = ArrayType([FT(0.003)])
        qᵣ = ArrayType([FT(0.002)])
        qₛ = ArrayType([FT(0.001)])

        kernel! = test_chen2022_terminal_velocity_kernel!(backend, work_groups)
        kernel!(lcl, icl, rain, snow, Ch2022, STVel, output, ρ, qₗ, qᵢ, qᵣ, qₛ; ndrange)
        v_term = Array(output)[1]

        # test that terminal velocity is callable and returns a reasonable value
        TT.@test v_term.lcl_stokes ≈ FT(0.021314907475574747)  # CloudLiquid with Stokes velocity
        TT.@test v_term.lci_chen2022 ≈ FT(0.01696129041896599)   # CloudIce with Chen2022 small ice
        TT.@test v_term.rai_chen2022 ≈ FT(6.9241079942767305)    # Rain with Chen2022 rain
        TT.@test v_term.sno_chen2022 ≈ FT(0.9514450529349796)    # Snow with Chen2022 large ice
    end

    TT.@testset "0-moment microphysics kernels" begin
        DT = @NamedTuple{S_pr::FT, S_pr_ref::FT}
        (; ndrange, output) = setup_output(10, DT)

        liquid_frac = ArrayType([FT(0), FT(0.5), FT(1)])
        qc = ArrayType([FT(3e-3), FT(4e-3), FT(5e-3)])

        kernel! = test_0_moment_micro_kernel!(backend, work_groups)
        kernel!(p0m, liquid_frac, output, qc; ndrange)
        # test that remove_precipitation matches the direct formula
        TT.@test all(map(nt -> nt.S_pr_ref == nt.S_pr, Array(output)))
    end

    TT.@testset "1-moment microphysics kernels" begin
        bad_value = -99999.99
        (; output, ndrange) = setup_output(2, FT, bad_value)
        ql = ArrayType([FT(1e-3), FT(5e-4)])

        kernel! = test_1_moment_micro_acnv_kernel!(backend, work_groups)
        kernel!(output, rain.acnv1M, ql; ndrange)
        out = Array(output)

        # Sanity checks for the GPU KernelAbstractions workflow
        # See https://github.com/CliMA/SurfaceFluxes.jl/issues/142
        TT.@test !any(isequal(out, FT(bad_value)))
        TT.@test out == FT[5e-7, 0]

        DT = @NamedTuple{
            liq_rai::FT, ice_sno::FT, liq_sno::FT, ice_rai::FT,
            rai_sink::FT, sno_rai::FT, rai_sno::FT,
        }
        (; output, ndrange) = setup_output(2, DT)

        ρ = ArrayType([FT(1.2), FT(1.2)])
        qi = ArrayType([FT(0), FT(5e-4)])
        qs = ArrayType([FT(0), FT(5e-4)])
        ql = ArrayType([FT(0), FT(5e-4)])
        qr = ArrayType([FT(0), FT(5e-4)])

        kernel! = test_1_moment_micro_accretion_kernel!(backend, work_groups)
        kernel!(lcl, rain, icl, snow, ce, blk1mvel, output, ρ, qi, qs, ql, qr; ndrange)
        out0, out1 = Array(output)

        # test 1-moment accretion is callable and returns a reasonable value
        TT.@test all(iszero, out0)
        TT.@test out1.liq_rai ≈ FT(1.4150106417043544e-6)
        TT.@test out1.ice_sno ≈ FT(2.453070979562392e-7)
        TT.@test out1.liq_sno ≈ FT(2.453070979562392e-7)
        TT.@test out1.ice_rai ≈ FT(1.768763302130443e-6)
        TT.@test out1.rai_sink ≈ FT(3.590060148920767e-5)
        TT.@test out1.sno_rai ≈ FT(2.466313958248222e-4)
        TT.@test out1.rai_sno ≈ FT(6.830957197816771e-5)

        (; output, ndrange) = setup_output(3, FT)

        T_freeze = FT(273.15)
        ρ = ArrayType([FT(1.2), FT(1.2), FT(1.2)])
        T = ArrayType([T_freeze + 2, T_freeze + 2, T_freeze - 2])
        qs = ArrayType([FT(1e-4), FT(0), FT(1e-4)])

        kernel! = test_1_moment_micro_snow_melt_kernel!(backend, work_groups)
        kernel!(snow, blk1mvel, aps, tps, output, ρ, T, qs; ndrange)
        out = Array(output)

        # test if 1-moment snow melt is callable and returns reasonable values
        TT.@test out ≈ FT[9.516553267013084e-6, 0, 0]
    end

    TT.@testset "2-moment microphysics kernels" begin

        TT.@testset "autoconversion" begin
            DT = @NamedTuple{S_VarTSc::FT, S_LD2004::FT, S_TC1980::FT, S_B1994::FT, S_KK2000::FT}
            (; output, ndrange) = setup_output(10, DT)

            ql = constant_data(FT(2e-3); ndrange)
            ρ = constant_data(FT(1.2); ndrange)
            Nd = constant_data(FT(1e8); ndrange)

            kernel! = test_2_moment_acnv_kernel!(backend, work_groups)
            kernel!(KK2000, B1994, TC1980, LD2004, VarTSc, output, ql, ρ, Nd; ndrange)
            out = Array(output)

            TT.@test allequal(out)
            (; S_VarTSc, S_LD2004, S_TC1980, S_B1994, S_KK2000) = out[1]
            TT.@test S_VarTSc ≈ FT(2e-6)
            TT.@test S_LD2004 ≈ FT(1.6963072465911614e-6)
            TT.@test S_TC1980 ≈ FT(3.5482867084128596e-6)
            TT.@test S_B1994 ≈ FT(9.825462758968215e-7)
            TT.@test S_KK2000 ≈ FT(5.855332513368727e-8)
        end

        TT.@testset "accretion" begin
            DT = @NamedTuple{S_KK2000::FT, S_B1994::FT, S_TC1980::FT}
            (; output, ndrange) = setup_output(10, DT)

            ql = constant_data(FT(2e-3); ndrange)
            qr = constant_data(FT(5e-4); ndrange)
            ρ = constant_data(FT(1.2); ndrange)

            kernel! = test_2_moment_accr_kernel!(backend, work_groups)
            kernel!(KK2000, B1994, TC1980, output, ql, qr, ρ; ndrange)
            out = Array(output)
            TT.@test allequal(out)
            (; S_KK2000, S_B1994, S_TC1980) = out[1]
            TT.@test S_KK2000 ≈ FT(6.6548664e-6) rtol = 1e-6
            TT.@test S_B1994 ≈ FT(7.2e-6)
            TT.@test S_TC1980 ≈ FT(4.7e-6)
        end

        for SB in [SB2006, SB2006_no_limiters]
            param_types = typeof.((aps, tps, SB, SB2006Vel))
            data = ntuple(i -> FT, 7)
            DT = Core.Compiler.return_type(SB2006_2M_kernel, Tuple{param_types..., data...})
            (; output, ndrange) = setup_output(10, DT)

            T = constant_data(FT(290); ndrange)
            qt = constant_data(FT(7e-3); ndrange)
            ql = constant_data(FT(2e-3); ndrange)
            qr = constant_data(FT(5e-4); ndrange)
            ρ = constant_data(FT(1.2); ndrange)
            Nl = constant_data(FT(1e8); ndrange)
            Nr = constant_data(FT(1e7); ndrange)

            kernel! = test_2_moment_SB2006_kernel!(backend, work_groups)
            kernel!(aps, tps, SB, SB2006Vel, output, qt, ql, qr, Nl, Nr, ρ, T; ndrange)
            out = Array(output)
            TT.@test allequal(out)
            tendencies = out[1]

            (; lcl_aconv_scoll, accr, rain_scoll, rain_vel, rain_evap, num_incr_mass_limit, num_decr_mass_limit) =
                tendencies

            TT.@test isapprox(lcl_aconv_scoll.au.dq_lcl_dt, FT(-5.742569998787898e-7), rtol = 1e-6)
            TT.@test isapprox(lcl_aconv_scoll.au.dN_lcl_dt, FT(-5300.833845034984), rtol = 1e-6)
            TT.@test isapprox(lcl_aconv_scoll.au.dq_rai_dt, FT(5.742569998787898e-7), rtol = 1e-6)
            TT.@test isapprox(lcl_aconv_scoll.au.dN_rai_dt, FT(2650.416922517492), rtol = 1e-6)
            TT.@test isapprox(lcl_aconv_scoll.sc, FT(-33859.96615496501), rtol = 1e-6)
            TT.@test isapprox(accr.dq_lcl_dt, FT(-6.358926e-6), rtol = 1e-6)
            TT.@test isapprox(accr.dN_lcl_dt, FT(-317946.28), rtol = 1e-6)
            TT.@test isapprox(accr.dq_rai_dt, FT(6.358926e-6), rtol = 1e-6)
            TT.@test isapprox(accr.dN_rai_dt, FT(0.0), rtol = 1e-6)
            if SB == SB2006
                TT.@test isapprox(rain_scoll.sc, FT(-21187.494), rtol = 1e-6)
                TT.@test isapprox(rain_scoll.br, FT(14154.027), rtol = 1e-6)
                TT.@test isapprox(rain_vel[1], FT(0.9868878), rtol = 1e-6)
                TT.@test isapprox(rain_vel[2], FT(4.517734), rtol = 1e-6)
                TT.@test isapprox(rain_evap.evap_rate_0, FT(-260791.30068415933), rtol = 1e-6)
                TT.@test isapprox(rain_evap.evap_rate_1, FT(-0.003709529301871412), rtol = 1e-6)
            end
            if SB == SB2006_no_limiters
                TT.@test isapprox(rain_scoll.sc, FT(-40447.855), rtol = 1e-6)
                TT.@test isapprox(rain_scoll.br, FT(0), rtol = 1e-6)
                TT.@test isapprox(rain_vel[1], FT(2.6429e-3), rtol = 1e-4)
                TT.@test isapprox(rain_vel[2], FT(0.1149338), rtol = 1e-5)
                TT.@test isapprox(rain_evap.evap_rate_0, FT(-56716.556198709244), rtol = 1e-6)
                TT.@test isapprox(rain_evap.evap_rate_1, FT(-0.00010034697555076008), rtol = 1e-6)
            end
            TT.@test isapprox(num_incr_mass_limit, FT(0), rtol = 1e-6)
            TT.@test isapprox(num_decr_mass_limit, FT(-7.692307e4), rtol = 1e-6)
        end
    end

    TT.@testset "Common Kernels" begin

        TT.@testset "H2SO4 kernels" begin
            DT = @NamedTuple{H2SO4_soln_qv_sat::FT, H2SO4_soln_a_w::FT}
            (; output, ndrange) = setup_output(10, DT)

            T = constant_data(FT(230); ndrange)
            x_sulph = constant_data(FT(0.1); ndrange)

            kernel! = Common_H2SO4_kernel!(backend, work_groups)
            kernel!(H2SO4_prs, tps, output, x_sulph, T; ndrange)
            out = Array(output)
            TT.@test allequal(out)
            (; H2SO4_soln_qv_sat, H2SO4_soln_a_w) = out[1]

            # test CO.H2SO4_soln_saturation_vapor_pressure
            TT.@test H2SO4_soln_qv_sat ≈ FT(12.685507586924)
            # test CO.a_w_xT
            TT.@test H2SO4_soln_a_w ≈ FT(0.928418590276476)
        end

        TT.@testset "CO.a_w_eT" begin
            (; output, ndrange) = setup_output(10, FT)

            T = constant_data(FT(282); ndrange)
            e = constant_data(FT(1001); ndrange)

            kernel! = Common_a_w_eT_kernel!(backend, work_groups)
            kernel!(tps, output, e, T; ndrange)
            out = Array(output)
            TT.@test allequal(out)
            a_w_eT = out[1]

            # test CO.a_w_eT
            TT.@test a_w_eT ≈ FT(0.880951366899518)
        end

        TT.@testset "CO.a_w_ice" begin
            (; output, ndrange) = setup_output(10, FT)

            T = constant_data(FT(230); ndrange)

            kernel! = Common_a_w_ice_kernel!(backend, work_groups)
            kernel!(tps, output, T; ndrange)
            out = Array(output)
            TT.@test allequal(out)
            a_w_ice = out[1]

            # test CO.a_w_ice
            TT.@test a_w_ice ≈ FT(0.6538439184585567)
        end
    end   # TT.@testset "Common Kernels"

    TT.@testset "Ice Nucleation kernels" begin
        TT.@testset "CMI_het.IN.dust_activated_number_fraction" begin
            DT = @NamedTuple{desert_act_N_frac::FT, arizona_act_N_frac::FT}
            (; output, ndrange) = setup_output(10, DT)

            T = constant_data(FT(240); ndrange)
            S_i = constant_data(FT(1.2); ndrange)

            kernel! = IceNucleation_dust_activated_number_fraction_kernel!(backend, work_groups)
            kernel!(desert_dust, arizona_test_dust, ip.deposition, output, S_i, T; ndrange)
            out = Array(output)
            TT.@test allequal(out)
            (; desert_act_N_frac, arizona_act_N_frac) = out[1]

            # test IN.dust_activated_number_fraction
            TT.@test desert_act_N_frac ≈ FT(0.0129835639)
            TT.@test arizona_act_N_frac ≈ FT(1.2233164999)
        end

        TT.@testset "CMI_het.IN.MohlerDepositionRate" begin
            DT = @NamedTuple{desert_dep_rate::FT, arizona_dep_rate::FT}
            (; output, ndrange) = setup_output(10, DT)

            T = constant_data(FT(240); ndrange)
            S_i = constant_data(FT(1.2); ndrange)
            dSi_dt = FT(0.03)
            N_aero = FT(3000)

            kernel! = IceNucleation_MohlerDepositionRate_kernel!(backend, work_groups)
            kernel!(desert_dust, arizona_test_dust, ip.deposition, output, S_i, T, dSi_dt, N_aero; ndrange)
            out = Array(output)
            TT.@test allequal(out)
            (; desert_dep_rate, arizona_dep_rate) = out[1]

            # test IN.MohlerDepositionRate
            TT.@test desert_dep_rate ≈ FT(38.6999999999)
            TT.@test arizona_dep_rate ≈ FT(423)
        end

        TT.@testset "CMI_het.IN.deposition_J" begin
            (; output, ndrange) = setup_output(10, FT)

            Δa_ws = FT[0.16, 0.15, 0.15]
            minerals = [kaolinite, feldspar, ferrihydrite]
            ref_dep_J = FT[1.5390757663075784e6, 5.693312205851678e6, 802555.3607426438]
            for (Δa_w, mineral, ref_J) in zip(Δa_ws, minerals, ref_dep_J)
                Delta_a_w = constant_data(Δa_w; ndrange)
                kernel! = IceNucleation_deposition_J_kernel!(backend, work_groups)
                kernel!(mineral, output, Delta_a_w; ndrange)
                out = Array(output)
                TT.@test allequal(out)
                dep_J = out[1]
                TT.@test dep_J ≈ ref_J
            end
        end

        TT.@testset "CMI_het.IN.ABIFM_J" begin
            (; output, ndrange) = setup_output(2, FT)

            Δa_ws = FT[0.16, 0.15]
            minerals = [kaolinite, illite]
            ref_J = FT[153.65772539109, 31.870032033791]
            for (Δa_w, mineral, ref_J) in zip(Δa_ws, minerals, ref_J)
                Delta_a_w = constant_data(Δa_w; ndrange)
                kernel! = IceNucleation_ABIFM_J_kernel!(backend, work_groups)
                kernel!(mineral, output, Delta_a_w; ndrange)
                out = Array(output)
                TT.@test allequal(out)
                J = out[1]
                TT.@test J ≈ ref_J
            end
        end

        TT.@testset "CMI_het.P3_deposition_N_i" begin
            (; output, ndrange) = setup_output(1, FT)

            ip_P3 = ip.p3
            T = constant_data(FT(240); ndrange)

            kernel! = IceNucleation_P3_deposition_N_i_kernel!(backend, work_groups)
            kernel!(ip_P3, output, T; ndrange)

            # test if P3_deposition_N_i is callable and returns reasonable values
            TT.@test Array(output)[1] ≈ FT(119018.93920746)
        end

        TT.@testset "CMI_het.P3_het_N_i" begin
            (; output, ndrange) = setup_output(1, FT)

            ip_P3 = ip.p3
            T = constant_data(FT(240); ndrange)
            N_l = constant_data(FT(2000); ndrange)
            V_l = constant_data(FT(3e-18); ndrange)
            Δt = constant_data(FT(0.1); ndrange)

            kernel! = IceNucleation_P3_het_N_i_kernel!(backend, work_groups)
            kernel!(ip_P3, output, T, N_l, V_l, Δt; ndrange)

            # test if P3_het_N_i is callable and returns reasonable values
            ref_val = if FT == Float64
                FT(0.0002736160475969029)
            else
                # loss of precision due to
                # `exp(-B * V_l_converted * Δt * exp(a * Tₛ))` -> 0.9999999 (Float32)
                # instead of 0.9999998631919762 (Float64).
                FT(0.00023841858f0)
            end
            TT.@test Array(output)[1] ≈ ref_val
        end

        TT.@testset "CMI_het.INP_concentration_frequency" begin
            (; output, ndrange) = setup_output(10, FT)

            T = constant_data(FT(233); ndrange)
            INPC = constant_data(FT(220000); ndrange)

            kernel! = IceNucleation_INPC_frequency_kernel!(backend, work_groups)
            kernel!(ip_frostenberg, output, INPC, T; ndrange)
            out = Array(output)

            # test INPC_frequency is callable and returns a reasonable value
            TT.@test allequal(out)
            TT.@test out[1] ≈ FT(0.26) rtol = 0.1
        end

        TT.@testset "CMI_het.homogeneous_J_[cubic/linear]" begin
            DT = @NamedTuple{J_cubic::FT, J_linear::FT}
            (; output, ndrange) = setup_output(10, DT)

            T = constant_data(FT(220); ndrange)
            Delta_a_w = constant_data(FT(0.2907389666103033); ndrange)

            kernel! = IceNucleation_homogeneous_J_kernel!(backend, work_groups)
            kernel!(ip, output, Delta_a_w; ndrange)
            out = Array(output)
            TT.@test allequal(out)
            (; J_cubic, J_linear) = out[1]

            # test homogeneous_J_cubic is callable and returns a reasonable value
            TT.@test J_cubic ≈ FT(2.66194650334444e12)
            # test homogeneous_J_linear is callable and returns a reasonable value
            TT.@test J_linear ≈ FT(7.156568123338207e11)
        end
    end  # TT.@testset "Ice nucleation kernels"

    TT.@testset "Modal nucleation kernels" begin
        (; output, ndrange) = setup_output(5, FT)

        TT.@testset "h2so4 nucleation" begin
            h2so4_conc = constant_data(FT(1e12); ndrange)
            nh3_conc = constant_data(FT(1); ndrange)
            negative_ion_conc = constant_data(FT(1); ndrange)
            temp = constant_data(FT(208); ndrange)

            kernel! = h2so4_nucleation_kernel!(backend, work_groups)
            kernel!(h2so4_nuc, output, h2so4_conc, nh3_conc, negative_ion_conc, temp; ndrange)

            TT.@test all(Array(output) .> FT(0))
        end

        TT.@testset "Organic nucleation" begin
            negative_ion_conc = constant_data(FT(0.0); ndrange)
            monoterpene_conc = constant_data(FT(1e24); ndrange)
            O3_conc = constant_data(FT(1e24); ndrange)
            OH_conc = constant_data(FT(1e24); ndrange)
            temp = constant_data(FT(300); ndrange)
            condensation_sink = constant_data(FT(1); ndrange)

            kernel! = organic_nucleation_kernel!(backend, work_groups)
            kernel!(
                org_nuc, output, negative_ion_conc, monoterpene_conc,
                O3_conc, OH_conc, temp, condensation_sink; ndrange,
            )

            TT.@test all(Array(output) .> FT(0))
        end

        TT.@testset "Organic and h2so4 nucleation" begin
            h2so4_conc = constant_data(FT(2.6e6); ndrange)
            monoterpene_conc = constant_data(FT(1); ndrange)
            OH_conc = constant_data(FT(1); ndrange)
            temp = constant_data(FT(300); ndrange)
            condensation_sink = constant_data(FT(1); ndrange)

            kernel! = organic_and_h2so4_nucleation_kernel!(backend, work_groups)
            kernel!(mix_nuc, output, h2so4_conc, monoterpene_conc, OH_conc, temp, condensation_sink; ndrange)

            TT.@test all(Array(output) .> FT(0))
        end

        TT.@testset "Apparent nucleation rate" begin
            output_diam = constant_data(FT(1e-9); ndrange)
            nucleation_rate = constant_data(FT(1e6); ndrange)
            condensation_growth_rate = constant_data(FT(1e6); ndrange)
            coag_sink = constant_data(FT(1e6); ndrange)
            coag_sink_input_diam = constant_data(FT(1e-9); ndrange)
            input_diam = constant_data(FT(1e-9); ndrange)

            kernel! = apparent_nucleation_rate_kernel!(backend, work_groups)
            kernel!(
                output, output_diam, nucleation_rate, condensation_growth_rate,
                coag_sink, coag_sink_input_diam, input_diam; ndrange,
            )
        end
    end  # TT.@testset "Modal nucleation kernels"

    TT.@testset "Bulk microphysics tendencies kernels" begin
        # 0M tests (returns dq_tot_dt, e_int_precip)
        ndrange = 10
        DT = @NamedTuple{dq_tot_dt::FT, e_int_precip::FT}
        (; output) = setup_output(ndrange, DT)
        liquid_frac = constant_data(FT(0.5); ndrange)
        qc = constant_data(FT(1e-3); ndrange)
        T = constant_data(FT(280.0); ndrange)

        kernel! = test_bulk_tendencies_0m_kernel!(backend, work_groups)
        TT.@testset "0M" begin
            kernel!(mp_0m, tps, output, liquid_frac, qc, T; ndrange)
            TT.@test allequal(Array(output))
            tendencies = Array(output)[1]
            TT.@test tendencies.dq_tot_dt ≤ 0  # Precipitation removal (dq_tot_dt) always negative or zero
            TT.@test isfinite(tendencies.e_int_precip)
        end

        # 0M S_0 tests (with q_vap_sat)
        (; output) = setup_output(ndrange, DT)
        q_vap_sat = constant_data(FT(0.01); ndrange)
        kernel_s0! = test_bulk_tendencies_0m_S0_kernel!(backend, work_groups)
        TT.@testset "0M S_0" begin
            kernel_s0!(mp_0m, tps, output, liquid_frac, qc, T, q_vap_sat; ndrange)
            TT.@test allequal(Array(output))
            tendencies = Array(output)[1]
            TT.@test tendencies.dq_tot_dt ≤ 0
            TT.@test isfinite(tendencies.e_int_precip)
        end


        # 1M tests
        DT = @NamedTuple{dq_lcl_dt::FT, dq_icl_dt::FT, dq_rai_dt::FT, dq_sno_dt::FT}
        (; output) = setup_output(ndrange, DT)
        ρ = constant_data(FT(1.0); ndrange)
        T = constant_data(FT(280.0); ndrange)
        q_tot = constant_data(FT(5e-3); ndrange)
        q_lcl = constant_data(FT(1e-3); ndrange)
        q_icl = constant_data(FT(0.5e-3); ndrange)
        q_rai = constant_data(FT(0.2e-3); ndrange)
        q_sno = constant_data(FT(0.1e-3); ndrange)

        kernel! = test_bulk_tendencies_1m_kernel!(backend, work_groups)
        TT.@testset "1M" begin
            kernel!(mp_1m, tps, output, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno; ndrange)
            TT.@test allequal(Array(output))
            tendencies = Array(output)[1]
            TT.@test all(isfinite, tendencies)
        end

        # 2M warm rain tests
        DT = @NamedTuple{
            dq_lcl_dt::FT, dn_lcl_dt::FT, dq_rai_dt::FT, dn_rai_dt::FT,
            dq_ice_dt::FT, dq_rim_dt::FT, db_rim_dt::FT,
        }
        (; output) = setup_output(ndrange, DT)
        n_lcl = constant_data(FT(1e8); ndrange)
        n_rai = constant_data(FT(1e6); ndrange)

        kernel! = test_bulk_tendencies_2m_warm_kernel!(backend, work_groups)
        TT.@testset "2M warm" begin
            kernel!(mp_2m_warm, tps, output, ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai; ndrange)
            TT.@test allequal(Array(output))
            tendencies = Array(output)[1]
            TT.@test all(isfinite, tendencies)
            TT.@test iszero(tendencies.dq_ice_dt)  # Ice tendency is zero for warm-only
        end

        # 2M+P3 tests
        q_ice = constant_data(FT(0.3e-3); ndrange)
        n_ice = constant_data(FT(1e5); ndrange)
        q_rim = constant_data(FT(0.1e-3); ndrange)
        b_rim = constant_data(FT(1e-10); ndrange)

        kernel! = test_bulk_tendencies_2m_p3_kernel!(backend, work_groups)
        TT.@testset "2M+P3" begin
            kernel!(mp_2m_p3, tps, output, ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim; ndrange)
            TT.@test allequal(Array(output))
            tendencies = Array(output)[1]
            TT.@test all(isfinite, tendencies)
            TT.@test !iszero(tendencies.dq_ice_dt)
        end
    end  # TT.@testset "Bulk microphysics tendencies kernels"
end  # function test_gpu(FT)

TT.@testset "GPU tests ($FT)" for FT in (Float64, Float32)
    test_gpu(FT)
end
nothing
