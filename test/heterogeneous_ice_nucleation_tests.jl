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

    TT.@testset "F23 immersion limit rate" begin

        T_freeze = ip_frostenberg.T_freeze
        ρ = FT(1)
        τ = FT(300)

        # Above T_freeze: zero
        r_warm = CMI_het.f23_immersion_limit_rate(
            ip_frostenberg, T_freeze + FT(0.1), ρ; τ,
        )
        TT.@test r_warm.∂ₜn_frz == FT(0)

        # Cold T: rate is positive and matches INPC/(ρ·τ)
        T_cold = T_freeze - FT(20)
        r_cold = CMI_het.f23_immersion_limit_rate(ip_frostenberg, T_cold, ρ; τ)
        INPC_expected = exp(CMI_het.INP_concentration_mean(ip_frostenberg, T_cold)) / ρ / τ
        TT.@test r_cold.∂ₜn_frz ≈ INPC_expected rtol = sqrt(eps(FT))

        # Colder ⇒ larger rate
        r_colder = CMI_het.f23_immersion_limit_rate(
            ip_frostenberg, T_freeze - FT(30), ρ; τ,
        )
        TT.@test r_colder.∂ₜn_frz > r_cold.∂ₜn_frz

        # log_inpc_shift > 0 ⇒ larger rate
        r_shifted = CMI_het.f23_immersion_limit_rate(
            ip_frostenberg, T_cold, ρ; τ, inpc_log_shift = FT(1),
        )
        TT.@test r_shifted.∂ₜn_frz ≈ r_cold.∂ₜn_frz * exp(FT(1)) rtol = sqrt(eps(FT))

        # Type stability
        TT.@test eltype(r_cold.∂ₜn_frz) == FT
    end

    TT.@testset "F23 deposition rate" begin

        T_freeze = ip_frostenberg.T_freeze
        ρ = FT(1)
        τ_act = FT(300)
        n_ice = FT(0)
        ρ_i = FT(916.7)
        D_nuc = FT(10e-6)
        m_nuc = FT(π) / 6 * ρ_i * D_nuc^3
        # Strict-MM15 defaults (these are now the function's defaults; we
        # still pass them explicitly here to make the test setup
        # self-documenting).
        T_thresh_default = T_freeze - FT(15)
        S_i_thresh_default = FT(0.05)

        # Strict defaults are picked up automatically when not overridden.
        T_test = T_freeze - FT(20)   # below the -15 °C strict gate
        q_sat_test = TDI.saturation_vapor_specific_content_over_ice(tps, T_test, ρ)
        r_default = CMI_het.f23_deposition_rate(
            ip_frostenberg, tps, T_test, ρ, 2 * q_sat_test, FT(0), FT(0), n_ice;
            m_nuc, τ_act,
        )
        TT.@test r_default.∂ₜn_frz > FT(0)
        # And they're closed at T = -10 °C (above the -15 °C gate)
        T_warm_test = T_freeze - FT(10)
        q_sat_warm_test = TDI.saturation_vapor_specific_content_over_ice(tps, T_warm_test, ρ)
        r_default_closed = CMI_het.f23_deposition_rate(
            ip_frostenberg, tps, T_warm_test, ρ, 2 * q_sat_warm_test, FT(0), FT(0), n_ice;
            m_nuc, τ_act,
        )
        TT.@test r_default_closed.∂ₜn_frz == FT(0)

        # T_freeze - 20 °C, vapor strongly supersaturated wrt ice ⇒ both
        # gates open. The function takes q_tot/q_liq/q_ice and computes
        # q_vap internally, so we pass q_tot = q_vap_super, q_liq = 0,
        # q_ice = 0.
        T_cold = T_freeze - FT(20)
        q_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T_cold, ρ)
        q_vap_super = 2 * q_sat_ice          # S_i ≈ 1.0
        r_active = CMI_het.f23_deposition_rate(
            ip_frostenberg, tps, T_cold, ρ, q_vap_super, FT(0), FT(0), n_ice;
            m_nuc, T_thresh = T_thresh_default, S_i_thresh = S_i_thresh_default,
            τ_act,
        )
        TT.@test r_active.∂ₜn_frz > FT(0)
        TT.@test r_active.∂ₜq_frz > FT(0)
        # When the starter-mass term is the binding constraint (large
        # q_excess), ∂ₜq_frz = m_nuc · ∂ₜn_frz exactly:
        TT.@test r_active.∂ₜq_frz ≈ m_nuc * r_active.∂ₜn_frz   rtol = sqrt(eps(FT))

        # n_ice = INPC ⇒ depleted to zero ⇒ both n and q rates vanish
        INPC_at_T = exp(CMI_het.INP_concentration_mean(ip_frostenberg, T_cold)) / ρ
        r_depleted = CMI_het.f23_deposition_rate(
            ip_frostenberg, tps, T_cold, ρ, q_vap_super, FT(0), FT(0), INPC_at_T;
            m_nuc, T_thresh = T_thresh_default, S_i_thresh = S_i_thresh_default,
            τ_act,
        )
        TT.@test r_depleted.∂ₜn_frz == FT(0)
        TT.@test r_depleted.∂ₜq_frz == FT(0)

        # T above the -15 °C gate ⇒ zero. Use q_vap super-saturated AT
        # T_warm so we isolate the T-gate effect (otherwise S_i would
        # also be < threshold).
        T_warm = T_freeze - FT(10)
        q_sat_warm = TDI.saturation_vapor_specific_content_over_ice(tps, T_warm, ρ)
        q_vap_super_warm = 2 * q_sat_warm
        r_T_gate = CMI_het.f23_deposition_rate(
            ip_frostenberg, tps, T_warm, ρ, q_vap_super_warm, FT(0), FT(0), n_ice;
            m_nuc, T_thresh = T_thresh_default, S_i_thresh = S_i_thresh_default,
            τ_act,
        )
        TT.@test r_T_gate.∂ₜn_frz == FT(0)
        TT.@test r_T_gate.∂ₜq_frz == FT(0)

        # Subsaturated wrt ice ⇒ both rates zero (S_i gate closed, and
        # the vapor cap independently zeros ∂ₜq_frz).
        r_sub = CMI_het.f23_deposition_rate(
            ip_frostenberg, tps, T_cold, ρ, FT(0.5) * q_sat_ice, FT(0), FT(0), n_ice;
            m_nuc, T_thresh = T_thresh_default, S_i_thresh = S_i_thresh_default,
            τ_act,
        )
        TT.@test r_sub.∂ₜn_frz == FT(0)
        TT.@test r_sub.∂ₜq_frz == FT(0)

        # Tunable thresholds: warming the T_thresh opens the warm-side gate
        r_relaxed = CMI_het.f23_deposition_rate(
            ip_frostenberg, tps, T_warm, ρ, q_vap_super_warm, FT(0), FT(0), n_ice;
            m_nuc, T_thresh = T_freeze - FT(5), S_i_thresh = S_i_thresh_default,
            τ_act,
        )
        TT.@test r_relaxed.∂ₜn_frz > FT(0)
        TT.@test r_relaxed.∂ₜq_frz > FT(0)

        # Permissive thresholds (used by BMT) ⇒ always passes the gate
        r_permissive = CMI_het.f23_deposition_rate(
            ip_frostenberg, tps, T_warm, ρ, q_vap_super_warm, FT(0), FT(0), n_ice;
            m_nuc, T_thresh = FT(2000), S_i_thresh = FT(-2), τ_act,
        )
        TT.@test r_permissive.∂ₜn_frz > FT(0)
        TT.@test r_permissive.∂ₜq_frz > FT(0)

        # Vapor-excess cap: with an absurdly large m_nuc, the
        # m_nuc·∂ₜn_frz term dominates and the q_excess cap should bind.
        # Verify ∂ₜq_frz = q_excess/(2τ_act) exactly.
        T_super = T_freeze - FT(20)
        q_sat_super = TDI.saturation_vapor_specific_content_over_ice(tps, T_super, ρ)
        q_vap_tiny_excess = q_sat_super * FT(1.001)   # S_i = 0.001
        q_excess_tiny = q_vap_tiny_excess - q_sat_super
        r_capped = CMI_het.f23_deposition_rate(
            ip_frostenberg, tps, T_super, ρ, q_vap_tiny_excess, FT(0), FT(0), FT(0);
            m_nuc = FT(1e3),  # absurd m_nuc forces vapor cap to bind
            T_thresh = FT(2000), S_i_thresh = FT(-2), τ_act,
        )
        TT.@test r_capped.∂ₜq_frz ≈ q_excess_tiny / (2τ_act)   rtol = sqrt(eps(FT))

        # Type stability
        TT.@test eltype(r_active.∂ₜn_frz) == FT
        TT.@test eltype(r_active.∂ₜq_frz) == FT
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
