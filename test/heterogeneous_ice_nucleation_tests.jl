import Test as TT

import ClimaParams as CP

import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.HetIceNucleation as CMI_het

function test_heterogeneous_ice_nucleation(FT)

    # parameters for parameterizations
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    H2SO4_prs = CMP.H2SO4SolutionParameters(FT)
    ip = CMP.IceNucleationParameters(FT)
    ip_frostenberg = CMP.Frostenberg2023(FT)
    # more parameters for aerosol properties
    ATD = CMP.ArizonaTestDust(FT)
    desert_dust = CMP.DesertDust(FT)
    illite = CMP.Illite(FT)
    kaolinite = CMP.Kaolinite(FT)
    feldspar = CMP.Feldspar(FT)
    ferrihydrite = CMP.Ferrihydrite(FT)
    unsupported_sea_salt = CMP.Seasalt(FT)

    TT.@testset "dust_activation" begin

        T_warm = FT(250)
        T_cold = FT(210)
        Si_low = FT(1.01)
        Si_med = FT(1.2)
        Si_hgh = FT(1.34)
        Si_too_hgh = FT(1.5)
        dSi_dt = FT(0.05)
        dSi_dt_negative = FT(-0.3)
        N_aer = FT(3000)

        # Activate more in cold temperatures and higher supersaturations
        for dust in [ATD, desert_dust]
            TT.@test CMI_het.dust_activated_number_fraction(
                dust,
                ip.deposition,
                Si_hgh,
                T_warm,
            ) > CMI_het.dust_activated_number_fraction(
                dust,
                ip.deposition,
                Si_med,
                T_warm,
            )
            TT.@test CMI_het.dust_activated_number_fraction(
                dust,
                ip.deposition,
                Si_med,
                T_cold,
            ) > CMI_het.dust_activated_number_fraction(
                dust,
                ip.deposition,
                Si_med,
                T_warm,
            )
            TT.@test CMI_het.MohlerDepositionRate(
                dust,
                ip.deposition,
                Si_med,
                T_cold,
                dSi_dt,
                N_aer,
            ) > CMI_het.MohlerDepositionRate(
                dust,
                ip.deposition,
                Si_med,
                T_warm,
                dSi_dt,
                N_aer,
            )
        end

        # no activation if saturation exceeds allowed value
        for dust in [ATD, desert_dust]
            for T in [T_warm, T_cold]
                TT.@test_throws AssertionError("Si < ip.Sᵢ_max") CMI_het.dust_activated_number_fraction(
                    dust,
                    ip.deposition,
                    Si_too_hgh,
                    T,
                )
                TT.@test_throws AssertionError("Si < ip.Sᵢ_max") CMI_het.MohlerDepositionRate(
                    dust,
                    ip.deposition,
                    Si_too_hgh,
                    T,
                    dSi_dt,
                    N_aer,
                )
            end
        end

        # no activation if dSi_dt is negative
        for dust in [ATD, desert_dust]
            for T in [T_warm, T_cold]
                TT.@test CMI_het.MohlerDepositionRate(
                    dust,
                    ip.deposition,
                    Si_low,
                    T,
                    dSi_dt_negative,
                    N_aer,
                ) == FT(0)
            end
        end
    end

    TT.@testset "Deposition Nucleation J" begin

        T_warm_1 = FT(229.2)
        T_cold_1 = FT(228.8)
        x_sulph = FT(0.1)

        T_warm_2 = FT(285)
        T_cold_2 = FT(251)
        e_warm = FT(1088)
        e_cold = FT(544)

        # higher nucleation rate at colder temperatures
        for dust in [feldspar, ferrihydrite, kaolinite]
            TT.@test CMI_het.deposition_J(
                dust,
                CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_cold_1) -
                CO.a_w_ice(tps, T_cold_1),
            ) > CMI_het.deposition_J(
                dust,
                CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_warm_1) -
                CO.a_w_ice(tps, T_warm_1),
            )

            TT.@test CMI_het.deposition_J(
                dust,
                CO.a_w_eT(tps, e_cold, T_cold_2) - CO.a_w_ice(tps, T_cold_2),
            ) > CMI_het.deposition_J(
                dust,
                CO.a_w_eT(tps, e_warm, T_warm_2) - CO.a_w_ice(tps, T_warm_2),
            )
        end

        # if unsupported aerosol type, default to J = 0
        TT.@test CMI_het.deposition_J(
            unsupported_sea_salt,
            CO.a_w_eT(tps, e_cold, T_cold_2) - CO.a_w_ice(tps, T_cold_2),
        ) == 0
    end

    TT.@testset "P3 Deposition Nᵢ" begin

        T_warm = FT(235)
        T_cold = FT(234)

        T_too_cold = FT(232)

        # higher ice concentration at colder temperatures
        TT.@test CMI_het.P3_deposition_N_i(ip.p3, T_cold) >
                 CMI_het.P3_deposition_N_i(ip.p3, T_warm)

        # if colder than threshold T, use threshold T
        TT.@test CMI_het.P3_deposition_N_i(ip.p3, T_too_cold) ==
                 CMI_het.P3_deposition_N_i(ip.p3, ip.p3.T_dep_thres)
    end

    TT.@testset "ABIFM J" begin

        T_warm_1 = FT(229.2)
        T_cold_1 = FT(228.8)
        x_sulph = FT(0.1)

        T_warm_2 = FT(285)
        T_cold_2 = FT(251)
        e_warm = FT(1088)
        e_cold = FT(544)

        # higher nucleation rate at colder temperatures
        for dust in [illite, kaolinite, desert_dust]
            TT.@test CMI_het.ABIFM_J(
                dust,
                CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_cold_1) -
                CO.a_w_ice(tps, T_cold_1),
            ) > CMI_het.ABIFM_J(
                dust,
                CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_warm_1) -
                CO.a_w_ice(tps, T_warm_1),
            )

            TT.@test CMI_het.ABIFM_J(
                dust,
                CO.a_w_eT(tps, e_cold, T_cold_2) - CO.a_w_ice(tps, T_cold_2),
            ) > CMI_het.ABIFM_J(
                dust,
                CO.a_w_eT(tps, e_warm, T_warm_2) - CO.a_w_ice(tps, T_warm_2),
            )
        end

        # if unsupported aerosol type, default to J = 0
        TT.@test CMI_het.ABIFM_J(
            unsupported_sea_salt,
            CO.a_w_eT(tps, e_cold, T_cold_2) - CO.a_w_ice(tps, T_cold_2),
        ) == 0
    end

    TT.@testset "P3 Heterogeneous Nᵢ" begin

        T_warm = FT(235)
        T_cold = FT(234)
        N_liq = FT(2e5)
        r_l = FT(2e-5)
        V_l = FT(4 / 3 * FT(π) * r_l^3)
        Δt = FT(0.1)

        # higher ice concentration at colder temperatures
        TT.@test CMI_het.P3_het_N_i(ip.p3, T_cold, N_liq, V_l, Δt) >
                 CMI_het.P3_het_N_i(ip.p3, T_warm, N_liq, V_l, Δt)
    end

    TT.@testset "Frostenberg" begin

        temperatures = FT.([233, 257])
        INPCs = FT.([220000, 9])
        frequencies = FT.([0.26, 0.08])

        for (T, INPC, frequency) in zip(temperatures, INPCs, frequencies)
            TT.@test CMI_het.INP_concentration_frequency(
                ip_frostenberg,
                INPC,
                T,
            ) ≈ frequency rtol = 0.1
        end
    end
end

TT.@testset "Heterogeneous Ice Nucleation Tests ($FT)" for FT in (Float64, Float32)
    test_heterogeneous_ice_nucleation(FT)
end
nothing
