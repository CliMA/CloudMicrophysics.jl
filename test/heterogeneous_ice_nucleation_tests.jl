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
        N_lcl = FT(2e5)
        r_l = FT(2e-5)
        V_l = FT(4 / 3 * FT(π) * r_l^3)
        Δt = FT(0.1)

        # higher ice concentration at colder temperatures
        TT.@test CMI_het.P3_het_N_i(ip.p3, T_cold, N_lcl, V_l, Δt) >
                 CMI_het.P3_het_N_i(ip.p3, T_warm, N_lcl, V_l, Δt)
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

        # test T > T_freeze
        T_warm = ip_frostenberg.T_freeze + FT(1)
        TT.@test CMI_het.INP_concentration_frequency(
            ip_frostenberg,
            INPCs[1],
            T_warm,
        ) == FT(0)
    end

    TT.@testset "Frostenberg b coefficient" begin

        # T_celsius = -20
        T_cold = ip_frostenberg.T_freeze - FT(20)

        # b = 1 default
        TT.@test CMI_het.INP_concentration_mean(ip_frostenberg, T_cold) ≈
                 FT(6.238324625039508)

        # b shifts the mean by +9 log(b)
        ip_b = CMP.Frostenberg2023{FT}(;
            σ = ip_frostenberg.σ,
            b = FT(2),
            T_freeze = ip_frostenberg.T_freeze,
            τ = ip_frostenberg.τ,
            T_thr = ip_frostenberg.T_thr,
            S_thr = ip_frostenberg.S_thr,
        )
        TT.@test CMI_het.INP_concentration_mean(ip_b, T_cold) ≈
                 FT(12.476649250079015)
    end

    TT.@testset "F23 immersion limit rate" begin

        T_freeze = ip_frostenberg.T_freeze
        ρ = FT(1)

        # Above T_freeze: zero
        r_warm = CMI_het.immersion_limit_rate(ip_frostenberg, T_freeze + FT(0.1), ρ)
        TT.@test r_warm.∂ₜn_frz == FT(0)

        # Cold T: rate is positive and matches INPC/(ρ·τ)
        T_cold = T_freeze - FT(20)
        r_cold = CMI_het.immersion_limit_rate(ip_frostenberg, T_cold, ρ)
        INPC_expected = exp(CMI_het.INP_concentration_mean(ip_frostenberg, T_cold)) / ρ / ip_frostenberg.τ
        TT.@test r_cold.∂ₜn_frz ≈ INPC_expected rtol = sqrt(eps(FT))

        # Colder ⇒ larger rate
        r_colder = CMI_het.immersion_limit_rate(ip_frostenberg, T_freeze - FT(30), ρ)
        TT.@test r_colder.∂ₜn_frz > r_cold.∂ₜn_frz

        # A stochastic spectrum with log_shift = 1 ⇒ rate larger by exp(1); zero below T_min
        ip_shifted = CMP.StochasticINPSpectrum(;
            base = ip_frostenberg, log_shift = FT(1), T_min = T_freeze - FT(25), N_max = FT(Inf),
        )
        r_shifted = CMI_het.immersion_limit_rate(ip_shifted, T_cold, ρ)
        TT.@test r_shifted.∂ₜn_frz ≈ r_cold.∂ₜn_frz * exp(FT(1)) rtol = sqrt(eps(FT))
        TT.@test CMI_het.immersion_limit_rate(ip_shifted, T_freeze - FT(30), ρ).∂ₜn_frz == FT(0)

        # Type stability
        TT.@test eltype(r_cold.∂ₜn_frz) == FT
    end

    TT.@testset "Deposition nucleation toward an INP target spectrum" begin
        mp = CMP.Microphysics2MParams(FT; with_ice = true)
        f23 = mp.ice.ice_nucleation
        m_nuc = CO.nucleated_ice_mass(mp.ice.scheme)
        ρ = FT(0.8)
        function rate(mp, T, S_i, n_ice)
            q_vap = (1 + S_i) * TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
            micro = (; q_tot = q_vap, q_lcl = FT(0), q_rai = FT(0), q_ice = FT(0), n_ice)
            return CMI_het.deposition_rate(mp, tps, micro, (; ρ, T))
        end
        with_spectrum(inp) = CMP.Microphysics2MParams(FT; with_ice = true, ice_nucleation = inp)

        # Relaxation toward the median of the climatology over τ, with crystals of the nucleated mass
        T, S_i, n_ice = FT(250), FT(0.2), FT(10)
        r = rate(mp, T, S_i, n_ice)
        ∂ₜn_ref = (f23(T) / ρ - n_ice) / f23.τ
        TT.@test r.∂ₜn_frz > 0
        TT.@test r.∂ₜn_frz ≈ ∂ₜn_ref rtol = 8 * eps(FT)
        TT.@test r.∂ₜq_frz ≈ m_nuc * ∂ₜn_ref rtol = 8 * eps(FT)
        TT.@test eltype(r.∂ₜn_frz) == FT && eltype(r.∂ₜq_frz) == FT

        # Zero once the activated concentration reaches the target
        r_depleted = rate(mp, T, S_i, 2 * f23(T) / ρ)
        TT.@test r_depleted.∂ₜn_frz == 0 && r_depleted.∂ₜq_frz == 0

        # Active below T_thr at an ice supersaturation of at least S_thr
        for (T, S_i) in ((f23.T_thr, FT(0.2)), (FT(250), FT(0.04)), (FT(250), FT(-0.1)))
            r_inactive = rate(mp, T, S_i, FT(0))
            TT.@test r_inactive.∂ₜn_frz == 0 && r_inactive.∂ₜq_frz == 0
        end
        TT.@test CMI_het.is_active(f23, FT(250), f23.S_thr)

        # A stochastic spectrum shifts the log of the target, with a ceiling and a lowest temperature
        shifted =
            with_spectrum(CMP.StochasticINPSpectrum(; base = f23, log_shift = FT(1), T_min = FT(240), N_max = FT(Inf)))
        TT.@test rate(shifted, T, S_i, FT(0)).∂ₜn_frz ≈ exp(FT(1)) * rate(mp, T, S_i, FT(0)).∂ₜn_frz rtol = 8 * eps(FT)
        TT.@test rate(shifted, FT(239), S_i, FT(0)).∂ₜn_frz == 0
        capped = with_spectrum(CMP.StochasticINPSpectrum(; base = f23, log_shift = FT(0), T_min = FT(0), N_max = FT(1)))
        TT.@test rate(capped, T, S_i, FT(0)).∂ₜn_frz ≈ FT(1) / ρ / f23.τ rtol = 8 * eps(FT)

        # The Frostenberg 2023 mass rate is at most half the vapor excess per relaxation time
        dense =
            with_spectrum(CMP.StochasticINPSpectrum(; base = f23, log_shift = FT(30), T_min = FT(0), N_max = FT(Inf)))
        q_excess = S_i * TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        TT.@test rate(dense, T, S_i, FT(0)).∂ₜq_frz ≈ q_excess / (2 * f23.τ) rtol = 8 * eps(FT)
    end

    TT.@testset "Cloud-droplet immersion freezing (Bigg + cloud PSD)" begin

        toml_dict = CP.create_toml_dict(FT)
        rf = CMP.RainFreezing(toml_dict)
        pdf_c = CMP.CloudParticlePDF_SB2006(toml_dict)
        pdf_r = CMP.RainParticlePDF_SB2006_limited(toml_dict)
        T_freeze = TDI.TD.Parameters.T_freeze(tps)

        ρ = FT(1)
        N_lcl = FT(1e8)
        q_lcl = FT(5e-4)

        # Below the −4 °C gate ⇒ positive number and mass rates
        T_cold = T_freeze - FT(20)
        r_cold = CMI_het.liquid_freezing_rate(rf, pdf_c, tps, q_lcl, ρ, N_lcl, T_cold)
        TT.@test r_cold.∂ₜn_frz > FT(0)
        TT.@test r_cold.∂ₜq_frz > FT(0)

        # Colder ⇒ larger rate (Bigg's exp(a·ΔT) is monotonically increasing)
        r_colder = CMI_het.liquid_freezing_rate(
            rf, pdf_c, tps, q_lcl, ρ, N_lcl, T_freeze - FT(30),
        )
        TT.@test r_colder.∂ₜn_frz > r_cold.∂ₜn_frz
        TT.@test r_colder.∂ₜq_frz > r_cold.∂ₜq_frz

        # Above −4 °C gate ⇒ both rates zero
        r_warm = CMI_het.liquid_freezing_rate(
            rf, pdf_c, tps, q_lcl, ρ, N_lcl, T_freeze - FT(2),
        )
        TT.@test r_warm.∂ₜn_frz == FT(0)
        TT.@test r_warm.∂ₜq_frz == FT(0)

        # Zero N or q ⇒ both rates zero
        r_zero_N = CMI_het.liquid_freezing_rate(rf, pdf_c, tps, q_lcl, ρ, FT(0), T_cold)
        TT.@test r_zero_N.∂ₜn_frz == FT(0)
        TT.@test r_zero_N.∂ₜq_frz == FT(0)
        r_zero_q = CMI_het.liquid_freezing_rate(rf, pdf_c, tps, FT(0), ρ, N_lcl, T_cold)
        TT.@test r_zero_q.∂ₜn_frz == FT(0)
        TT.@test r_zero_q.∂ₜq_frz == FT(0)

        # Cloud and rain methods both run and return finite rates with the
        # same RainFreezing parameters — order-of-magnitude sanity only,
        # not equality (PSD shape and number/mass differ).
        r_cld = CMI_het.liquid_freezing_rate(rf, pdf_c, tps, q_lcl, ρ, N_lcl, T_cold)
        r_rai = CMI_het.liquid_freezing_rate(
            rf, pdf_r, tps, FT(1e-4), ρ, FT(1e3), T_cold,
        )
        TT.@test isfinite(r_cld.∂ₜn_frz) && r_cld.∂ₜn_frz > FT(0)
        TT.@test isfinite(r_rai.∂ₜn_frz) && r_rai.∂ₜn_frz > FT(0)

        # Type stability
        TT.@test eltype(r_cold.∂ₜn_frz) == FT
        TT.@test eltype(r_cold.∂ₜq_frz) == FT
    end
end

TT.@testset "Heterogeneous Ice Nucleation Tests ($FT)" for FT in (Float64, Float32)
    test_heterogeneous_ice_nucleation(FT)
end
nothing
