import Test as TT

import ClimaParams as CP

import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.HetIceNucleation as CMI_het
import CloudMicrophysics.HomIceNucleation as CMI_hom
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.DistributionTools as DT
import SpecialFunctions as SF

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

    TT.@testset "Deposition nucleation toward the Cooper spectrum" begin
        inp = CMP.ExponentialSupercoolingINP(FT)
        mp = CMP.Microphysics2MParams(FT; with_ice = true, ice_nucleation = inp)
        p3 = mp.ice.scheme
        # A cold cirrus state at about 215 hPa
        T, ρ = FT(215), FT(0.348)
        function rate(T, S_i, n_ice)
            q_vap = (1 + S_i) * TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
            micro = (; q_tot = q_vap, q_lcl = FT(0), q_rai = FT(0), q_ice = FT(0), n_ice)
            return CMI_het.deposition_rate(mp, tps, micro, (; ρ, T))
        end

        # The target is exponential in supercooling up to the ceiling
        TT.@test inp(T) == inp.N_max
        TT.@test inp(FT(250)) ≈ inp.a * exp(inp.b * (inp.T₀ - FT(250))) < inp.N_max

        # Relaxation at the rate at which deposition grows a crystal to the nucleated radius
        S_i = FT(0.3)
        k = CMI_het.delivery_rate(inp, mp, tps, T, S_i)
        G = CO.G_func_ice(mp.warm_rain.air_properties, tps, T)
        TT.@test k ≈ 2 * G * S_i / (p3.ρ_i * (p3.D_nuc / 2)^2) rtol = 4 * eps(FT)
        TT.@test FT(1) < 1 / CMI_het.delivery_rate(inp, mp, tps, T, inp.S_thr) < FT(100)
        TT.@test CMI_het.delivery_rate(inp, mp, tps, T, FT(-0.3)) == 0

        # Crystals of the nucleated mass, relaxing toward the target
        r = rate(T, S_i, FT(0))
        TT.@test r.∂ₜn_frz ≈ k * inp(T) / ρ rtol = 8 * eps(FT)
        TT.@test r.∂ₜq_frz == CO.nucleated_ice_mass(p3) * r.∂ₜn_frz
        TT.@test FT(0) < rate(T, S_i, inp(T) / ρ / 2).∂ₜn_frz < r.∂ₜn_frz
        TT.@test rate(T, S_i, 2 * inp(T) / ρ).∂ₜn_frz == 0

        # Active below T_thr at an ice supersaturation of at least S_thr
        TT.@test rate(inp.T_thr, S_i, FT(0)).∂ₜn_frz == 0
        TT.@test rate(T, FT(0.049), FT(0)).∂ₜn_frz == 0
        TT.@test CMI_het.is_active(inp, T, inp.S_thr)

        # No INP budget for immersion freezing
        TT.@test CMI_het.immersion_limit_rate(inp, T, ρ).∂ₜn_frz == FT(Inf)
    end
end

function test_drop_freezing(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    T_freeze = TDI.TD.Parameters.T_freeze(tps)
    # The Cooper spectrum has no immersion budget
    mp = CMP.Microphysics2MParams(FT; with_ice = true, ice_nucleation = CMP.ExponentialSupercoolingINP(FT))
    mp_f23 = CMP.Microphysics2MParams(FT; with_ice = true)
    ρ, qᵥ = FT(0.9), FT(3e-3)
    state(q_r, N_r, q_c, N_c; n_ice = FT(0)) = (;
        q_tot = qᵥ + q_c + q_r, q_lcl = q_c, n_lcl = N_c / ρ, q_rai = q_r, n_rai = N_r / ρ, q_ice = FT(0), n_ice,
    )
    rates(mp, micro, T) = (
        CMI_het.rain_freezing_rate(mp, tps, micro, (; ρ, T)),
        CMI_het.cloud_freezing_rate(mp, tps, micro, (; ρ, T)),
    )
    q_r, N_r, q_c, N_c = FT(1e-4), FT(1e3), FT(1e-4), FT(1e8)

    TT.@testset "Gauss-Laguerre rule" begin
        for α in (FT(0), FT(1), FT(1.5))
            quad = CMI_het.QUAD.GaussLaguerre(FT, 8; α)
            TT.@test sum(quad.weights) ≈ 1 rtol = 8 * eps(FT)
            # Exact for the moments of the normalized gamma weight up to order 2n - 1
            for k in (1, 5, 15)
                exact = SF.gamma(Float64(α) + k + 1) / SF.gamma(Float64(α) + 1)
                TT.@test sum(quad.weights .* quad.nodes .^ k) ≈ exact rtol = 100 * eps(FT)
            end
        end
        TT.@test mp.ice.cloud_freezing_quad.α == (mp.ice.cloud_pdf.νc + 1) / mp.ice.cloud_pdf.μc - 1
    end

    TT.@testset "freezing needs supercooled liquid" begin
        for (micro, T) in (
            (state(q_r, N_r, q_c, N_c), T_freeze),
            (state(FT(0), N_r, FT(0), N_c), T_freeze - 10),
            (state(q_r, FT(0), q_c, FT(0)), T_freeze - 10),
        )
            for r in rates(mp, micro, T)
                TT.@test r.∂ₜn_frz === FT(0) && r.∂ₜq_frz === FT(0)
            end
        end
        for r in rates(mp, state(q_r, N_r, q_c, N_c), prevfloat(T_freeze))
            TT.@test r.∂ₜn_frz > 0 && r.∂ₜq_frz > 0
        end
    end

    TT.@testset "freezing increases with supercooling" begin
        # Through immersion freezing, homogeneous freezing, and the supercooling at which the drop absorbs its
        # latent heat as it warms
        prev_r, prev_c = (FT(0), FT(0)), (FT(0), FT(0))
        for ΔT in FT(10) .^ range(FT(-3), log10(FT(90)); length = 40)
            r, c = rates(mp, state(q_r, N_r, q_c, N_c), T_freeze - ΔT)
            for v in (r.∂ₜn_frz, r.∂ₜq_frz, c.∂ₜn_frz, c.∂ₜq_frz)
                TT.@test isfinite(v) && v > 0
            end
            TT.@test r.∂ₜq_frz >= prev_r[2] * (1 - sqrt(eps(FT)))
            TT.@test c.∂ₜn_frz >= prev_c[1] * (1 - sqrt(eps(FT)))
            TT.@test c.∂ₜq_frz >= prev_c[2] * (1 - sqrt(eps(FT)))
            prev_r, prev_c = (r.∂ₜn_frz, r.∂ₜq_frz), (c.∂ₜn_frz, c.∂ₜq_frz)
        end
    end

    TT.@testset "rain number freezes at most at the mass fraction" begin
        micro = (; q_tot = FT(3e-3), q_lcl = FT(0), n_lcl = FT(0), q_rai = FT(1.89e-3), n_rai = FT(26400),
            q_ice = FT(0), n_ice = FT(0))
        for T in (T_freeze - 37, T_freeze - 36, T_freeze - 15, T_freeze - 5)
            r = CMI_het.rain_freezing_rate(mp, tps, micro, (; ρ = FT(0.8), T))
            k_q = r.∂ₜq_frz / micro.q_rai
            k_n = r.∂ₜn_frz / micro.n_rai
            TT.@test 0 < k_n <= k_q * (1 + 4 * eps(FT))
        end
        # Where nucleation limits freezing, the large drops freeze first, and the number rate is that of the size
        # distribution, about 1/20 of the mass fraction for an exponential distribution
        r = CMI_het.rain_freezing_rate(mp, tps, micro, (; ρ = FT(0.8), T = T_freeze - 5))
        TT.@test r.∂ₜn_frz / micro.n_rai < r.∂ₜq_frz / micro.q_rai / 10
    end

    TT.@testset "cloud droplet immersion freezing and its INP budget" begin
        # At 255 K the homogeneous coefficient is zero and nucleation limits the freezing of cloud droplets, so the
        # rates are the immersion moments of the size distribution
        T, ρ_c, q_lcl, N_lcl = FT(255), FT(0.9), FT(5e-4), FT(1e8)
        micro = (; q_tot = FT(2e-3) + q_lcl, q_lcl, n_lcl = N_lcl / ρ_c, q_rai = FT(0), n_rai = FT(0),
            q_ice = FT(0), n_ice = FT(10))
        pdf_c = mp.ice.cloud_pdf
        J = mp.ice.rain_freezing(T, T_freeze)
        (; λc, νcD, μcD) = CM2.pdf_cloud_parameters(pdf_c, q_lcl, ρ_c, N_lcl)
        M_D⁶ = DT.generalized_gamma_Mⁿ(νcD, μcD, λc, N_lcl / ρ_c, 6)
        unbounded = CMI_het.cloud_freezing_rate(mp, tps, micro, (; ρ = ρ_c, T))
        TT.@test unbounded.∂ₜn_frz ≈ J * q_lcl / pdf_c.ρw rtol = FT(1e-4)
        TT.@test unbounded.∂ₜq_frz ≈ J * pdf_c.ρw * (FT(π) / 6)^2 * M_D⁶ rtol = FT(1e-4)

        # The Frostenberg 2023 budget bounds the immersion number rate
        bounded = CMI_het.cloud_freezing_rate(mp_f23, tps, micro, (; ρ = ρ_c, T))
        ∂ₜn_lim = CMI_het.immersion_limit_rate(mp_f23.ice.ice_nucleation, T, ρ_c; n_active = micro.n_ice).∂ₜn_frz
        TT.@test FT(0) < ∂ₜn_lim < unbounded.∂ₜn_frz
        TT.@test bounded.∂ₜn_frz ≈ ∂ₜn_lim rtol = FT(1e-3)
    end

    TT.@testset "homogeneous freezing coefficient" begin
        hom = mp.ice.homogeneous
        for T in (FT(236), FT(237), FT(238))
            Δa_w = 1 - CO.a_w_ice(tps, T)
            TT.@test CMI_hom.homogeneous_freezing_rate_coefficient(hom, tps, T) ≈
                     CMI_hom.homogeneous_J_cubic(hom, Δa_w) rtol = 8 * eps(FT)
        end
        TT.@test CMI_hom.homogeneous_freezing_rate_coefficient(hom, tps, FT(260)) == 0
        J_edge = CMI_hom.homogeneous_J_cubic(hom, hom.Δa_w_max)
        TT.@test CMI_hom.homogeneous_freezing_rate_coefficient(hom, tps, FT(200)) ≈ J_edge rtol = 8 * eps(FT)
    end
end

TT.@testset "Heterogeneous Ice Nucleation Tests ($FT)" for FT in (Float64, Float32)
    test_heterogeneous_ice_nucleation(FT)
    test_drop_freezing(FT)
end
nothing
