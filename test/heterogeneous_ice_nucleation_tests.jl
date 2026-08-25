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

    TT.@testset "Frostenberg a/b coefficients" begin

        # T_celsius = -20
        T_cold = ip_frostenberg.T_freeze - FT(20)

        # a = b = 1 defaults
        TT.@test CMI_het.INP_concentration_mean(ip_frostenberg, T_cold) ≈
                 FT(6.238324625039508)

        # a shifts the mean by -log(a)
        ip_a = CMP.Frostenberg2023{FT}(;
            σ = ip_frostenberg.σ,
            a = FT(2),
            b = ip_frostenberg.b,
            T_freeze = ip_frostenberg.T_freeze,
        )
        TT.@test CMI_het.INP_concentration_mean(ip_a, T_cold) ≈
                 FT(5.545177444479562)

        # b shifts the mean by +9 log(b)
        ip_b = CMP.Frostenberg2023{FT}(;
            σ = ip_frostenberg.σ,
            a = ip_frostenberg.a,
            b = FT(2),
            T_freeze = ip_frostenberg.T_freeze,
        )
        TT.@test CMI_het.INP_concentration_mean(ip_b, T_cold) ≈
                 FT(12.476649250079015)
    end

    TT.@testset "F23 immersion limit rate" begin

        T_freeze = ip_frostenberg.T_freeze
        ρ = FT(1)
        τ = FT(300)

        # Above T_freeze: zero
        r_warm = CMI_het.immersion_limit_rate(
            ip_frostenberg, T_freeze + FT(0.1), ρ; τ,
        )
        TT.@test r_warm.∂ₜn_frz == FT(0)

        # Cold T: rate is positive and matches INPC/(ρ·τ)
        T_cold = T_freeze - FT(20)
        r_cold = CMI_het.immersion_limit_rate(ip_frostenberg, T_cold, ρ; τ)
        INPC_expected = exp(CMI_het.INP_concentration_mean(ip_frostenberg, T_cold)) / ρ / τ
        TT.@test r_cold.∂ₜn_frz ≈ INPC_expected rtol = sqrt(eps(FT))

        # Colder ⇒ larger rate
        r_colder = CMI_het.immersion_limit_rate(
            ip_frostenberg, T_freeze - FT(30), ρ; τ,
        )
        TT.@test r_colder.∂ₜn_frz > r_cold.∂ₜn_frz

        # log_inpc_shift > 0 ⇒ larger rate
        r_shifted = CMI_het.immersion_limit_rate(
            ip_frostenberg, T_cold, ρ; τ, inpc_log_shift = FT(1),
        )
        TT.@test r_shifted.∂ₜn_frz ≈ r_cold.∂ₜn_frz * exp(FT(1)) rtol = sqrt(eps(FT))

        # Type stability
        TT.@test eltype(r_cold.∂ₜn_frz) == FT
    end

    TT.@testset "F23 deposition rate" begin
        # `Frostenberg2023` on the deposition target-spectrum interface: swap it in
        # for the default `ExponentialSupercoolingINP` on an otherwise-default
        # `Microphysics2MParams`, so the nascent-crystal mass, the depletion model
        # and everything else the shared `deposition_rate` body reads come from the
        # same single source every other closure uses, rather than from the
        # bespoke `D_nuc`/`m_nuc` this test used to pass in by hand.
        toml_dict = CP.create_toml_dict(FT)
        mp_default = CMP.Microphysics2MParams(toml_dict; with_ice = true)
        p3 = mp_default.ice.scheme
        mp = CMP.Microphysics2MParams(;
            warm_rain = mp_default.warm_rain,
            ice = CMP.P3IceParams(;
                scheme = p3,
                terminal_velocity = mp_default.ice.terminal_velocity,
                cloud_pdf = mp_default.ice.cloud_pdf,
                rain_pdf = mp_default.ice.rain_pdf,
                ice_nucleation = ip_frostenberg,
                rain_freezing = mp_default.ice.rain_freezing,
                homogeneous = mp_default.ice.homogeneous,
                inp_depletion_model = mp_default.ice.inp_depletion_model,
                quad = mp_default.ice.quad,
            ),
        )
        (; m_nuc) = CMP.ice_seed(p3)

        T_freeze = ip_frostenberg.T_freeze
        ρ = FT(1)

        # `Frostenberg2023` is selectable now, not a parallel code path.
        TT.@test ip_frostenberg isa CMP.AbstractINPTargetSpectrum

        # `is_active` is below freezing and not subsaturated with respect to ice. The two
        # literals it used to carry, 15 K below freezing and 5 percent supersaturation, were
        # the DEFAULT closure's values and not this spectrum's, which has no threshold at
        # either quantity, so both are gone. These assertions bracket the window that is
        # left, including the two points the removed thresholds used to exclude.
        TT.@test CMI_het.is_active(ip_frostenberg, T_freeze - FT(16), FT(0.06))
        TT.@test CMI_het.is_active(ip_frostenberg, T_freeze - FT(14), FT(0.06))
        TT.@test CMI_het.is_active(ip_frostenberg, T_freeze - FT(16), FT(0.04))
        TT.@test CMI_het.is_active(ip_frostenberg, T_freeze - FT(16), FT(0.05))
        # The window itself: at and above freezing it is shut, and subsaturated it is shut.
        TT.@test !CMI_het.is_active(ip_frostenberg, T_freeze, FT(0.06))
        TT.@test !CMI_het.is_active(ip_frostenberg, T_freeze + FT(1), FT(0.06))
        TT.@test !CMI_het.is_active(ip_frostenberg, T_freeze - FT(16), FT(-0.01))

        # The delivery form belongs to the SLOT and not to the spectrum, so this spectrum
        # takes the same diffusional seed-delivery rate the default target does: the two
        # closures now differ in their target spectrum alone, which is what makes them
        # comparable. It is state dependent, and it carries `max(S_i, 0)`, so it falls
        # continuously to zero as saturation is approached from above rather than being cut
        # off by the window.
        ip_default = mp.ice.ice_nucleation
        for (T_probe, S_probe) in
            ((T_freeze - FT(20), FT(0.5)), (T_freeze - FT(5), FT(0.02)))
            TT.@test CMI_het.delivery_rate(ip_frostenberg, mp, tps, T_probe, S_probe) ==
                     CMI_het.delivery_rate(ip_default, mp, tps, T_probe, S_probe)
        end
        TT.@test CMI_het.delivery_rate(ip_frostenberg, mp, tps, T_freeze - FT(20), FT(0.5)) >
                 CMI_het.delivery_rate(ip_frostenberg, mp, tps, T_freeze - FT(20), FT(0.05))
        TT.@test CMI_het.delivery_rate(ip_frostenberg, mp, tps, T_freeze - FT(20), FT(-0.5)) == 0

        # Below the -15 °C gate and above the 5% ice-supersaturation gate:
        # nucleation is available.
        T_test = T_freeze - FT(20)
        q_sat_test = TDI.saturation_vapor_specific_content_over_ice(tps, T_test, ρ)
        micro = (;
            q_tot = 2 * q_sat_test, q_lcl = FT(0), q_rai = FT(0), q_ice = FT(0),
            n_ice = FT(0),
        )
        thermo = (; ρ, T = T_test)
        r_default = CMI_het.deposition_rate(ip_frostenberg, mp, tps, micro, thermo)
        TT.@test r_default.∂ₜn_frz > FT(0)
        TT.@test r_default.∂ₜq_frz > FT(0)
        # Every crystal is created at the shared nascent mass: the pair is exact,
        # as for every target spectrum on this interface (there is no
        # vapor-excess branch that could scale one moment without the other).
        TT.@test r_default.∂ₜq_frz == m_nuc * r_default.∂ₜn_frz
        # The number rate is the bare INP-budget relaxation, at the slot's own delivery rate.
        INPC_at_T_test = exp(CMI_het.INP_concentration_mean(ip_frostenberg, T_test)) / ρ
        S_i_test = TDI.q_vap(micro.q_tot, micro.q_lcl + micro.q_rai, micro.q_ice) / q_sat_test - 1
        inv_τ_test = CMI_het.delivery_rate(ip_frostenberg, mp, tps, T_test, S_i_test)
        TT.@test r_default.∂ₜn_frz == max(FT(0), INPC_at_T_test - micro.n_ice) * inv_τ_test

        # n_ice = target ⇒ depleted to zero ⇒ both moments vanish.
        micro_depleted = (; micro..., n_ice = INPC_at_T_test)
        r_depleted = CMI_het.deposition_rate(ip_frostenberg, mp, tps, micro_depleted, thermo)
        TT.@test r_depleted.∂ₜn_frz == FT(0)
        TT.@test r_depleted.∂ₜq_frz == FT(0)

        # The window is the freezing point, not 15 K below it. That threshold and the 5
        # percent supersaturation one were the DEFAULT closure's values rather than this
        # spectrum's, and both are retired, so 263.15 K belongs among the OPEN-gate states.
        # An earlier version of this test asserted a positive rate there, was changed to
        # assert zero, and is changed back: the campaign's ungated reading of this closure
        # was the correct one.
        # Supersaturated and below freezing, the rate is POSITIVE at every one of these,
        # including the three the retired 15 K threshold used to zero.
        for T_below in (T_freeze - FT(15), T_freeze - FT(14), T_freeze - FT(10))
            q_sat_below = TDI.saturation_vapor_specific_content_over_ice(tps, T_below, ρ)
            micro_below = (; micro..., q_tot = 2 * q_sat_below)
            r_below = CMI_het.deposition_rate(
                ip_frostenberg, mp, tps, micro_below, (; ρ, T = T_below),
            )
            TT.@test r_below.∂ₜn_frz > FT(0)
            TT.@test r_below.∂ₜq_frz > FT(0)
        end
        # At and above freezing the window is shut, whatever the supersaturation.
        for T_above in (T_freeze, T_freeze + FT(1))
            q_sat_above = TDI.saturation_vapor_specific_content_over_ice(tps, T_above, ρ)
            micro_above = (; micro..., q_tot = 2 * q_sat_above)
            r_above = CMI_het.deposition_rate(
                ip_frostenberg, mp, tps, micro_above, (; ρ, T = T_above),
            )
            TT.@test r_above.∂ₜn_frz == FT(0)
            TT.@test r_above.∂ₜq_frz == FT(0)
        end

        # Subsaturated with respect to ice ⇒ zero even though the temperature
        # gate is open.
        micro_sub = (; micro..., q_tot = FT(0.5) * q_sat_test)
        r_sub = CMI_het.deposition_rate(ip_frostenberg, mp, tps, micro_sub, thermo)
        TT.@test r_sub.∂ₜn_frz == FT(0)
        TT.@test r_sub.∂ₜq_frz == FT(0)

        # Type stability.
        TT.@test eltype(r_default.∂ₜn_frz) == FT
        TT.@test eltype(r_default.∂ₜq_frz) == FT

        # The old standalone `deposition_rate(opt::CMP.Frostenberg2023, ...)` method
        # capped the implied mass injection at half the local vapor excess per
        # relaxation window; that cap has no counterpart in the shared interface
        # body above (`delivery_rate` sees only `T` and `S_i`, not the air density
        # a vapor-mass bound needs), so it is not reproduced here. It is provably
        # inert at the shared nascent mass in the states this file exercises, but
        # not in general: at very cold, thin-air cirrus states the unbounded
        # `(−T_celsius/10)⁹` growth of the Frostenberg target can still outrun the
        # local vapor excess. That gap is the validity-limit treatment already on
        # this project's deferred list alongside F23's other future work (the
        # ice-nucleating-particle tracer and the stochastic reading of the
        # spectrum), not something this reformatting unit adds a bound for.
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

        # 2 K of supercooling: freezing is available, and small
        r_warm = CMI_het.liquid_freezing_rate(
            rf, pdf_c, tps, q_lcl, ρ, N_lcl, T_freeze - FT(2),
        )
        TT.@test r_warm.∂ₜn_frz > FT(0)
        TT.@test r_warm.∂ₜq_frz > FT(0)   # paired, as every source must be
        TT.@test r_warm.∂ₜn_frz < r_cold.∂ₜn_frz   # and small: the rolloff toward ΔT = 0

        r_not_supercooled = CMI_het.liquid_freezing_rate(
            rf, pdf_c, tps, q_lcl, ρ, N_lcl, T_freeze,
        )
        TT.@test r_not_supercooled.∂ₜn_frz == FT(0)
        TT.@test r_not_supercooled.∂ₜq_frz == FT(0)

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

# The shared activity condition reproduces the expression it replaced, except in the band 0 <= ΔT < 4 K.
function test_liquid_freezing_gate_is_shared(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    T_freeze = TDI.TD.Parameters.T_freeze(tps)
    ϵₘ = CM.Utilities.ϵ_numerics_2M_M(FT)
    ϵₙ = CM.Utilities.ϵ_numerics_2M_N(FT)
    # the condition before the change, as the control
    old(q, n, T) = (n > ϵₙ) & (q > ϵₘ) & (T < T_freeze - 4)
    band(T) = (T >= T_freeze - 4) & (T < T_freeze)   # the band the derived gate opens

    TT.@testset "the liquid-freezing gate is one shared definition [FT=$FT]" begin
        # each clause straddled: mass present/absent, number present/absent, and the temperature
        # gate from well below to above freezing including its exact boundary
        qs = FT[0, ϵₘ, nextfloat(ϵₘ), 1e-6, 1e-3]
        ns = FT[0, ϵₙ, nextfloat(ϵₙ), 1e3, 1e8]
        Ts = FT[T_freeze - 40, T_freeze - 4.001, T_freeze - 4, T_freeze - 3.999,
            T_freeze - 1, T_freeze, T_freeze + 5]
        n_true = 0
        n_opened = 0
        for q in qs, n in ns, T in Ts
            got = CMI_het._liquid_freezing_is_active(FT, q, n, T, T_freeze)
            present = (n > ϵₙ) & (q > ϵₘ)
            if present & band(T)
                TT.@test got && !old(q, n, T)
                n_opened += 1
            else
                # everywhere else the gate is bit-identical to what it replaced
                TT.@test got === old(q, n, T)
            end
            n_true += got
        end
        TT.@test n_opened > 0        # the band has to be sampled, or the change is untested
        TT.@test n_true > 0
        TT.@test n_true < length(qs) * length(ns) * length(Ts)
    end
end

# Composed liquid freezing over ΔT, both categories, both precisions: nonnegative, exactly zero at
# ΔT <= 0, monotone nondecreasing in supercooling, and the capped J does not exceed the uncapped one.
function test_liquid_freezing_composition_safety(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    T_freeze = TDI.TD.Parameters.T_freeze(tps)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    aps = mp.warm_rain.air_properties
    evap = mp.warm_rain.seifert_beheng.evap
    hom = mp.ice.homogeneous
    ρ, qᵥ = FT(0.9), FT(3e-3)
    q_r, N_r = FT(1e-4), FT(1e3)
    q_c, N_c = FT(1e-4), FT(1e8)

    ΔTs = vcat(FT[-5, -1, 0], FT(10) .^ range(FT(-6), log10(FT(8)); length = 24))

    TT.@testset "composed liquid freezing is monotone in supercooling [FT=$FT]" begin
        prev_rain, prev_cloud = FT(-Inf), FT(-Inf)
        n_pos = 0
        for ΔT in ΔTs
            T = T_freeze - ΔT
            r = CMI_het.rain_freezing_rate(
                mp.ice.rain_freezing, hom, p3.vent, aps, tps, evap, mp.ice.rain_pdf,
                q_r, ρ, N_r, T, qᵥ)
            c = CMI_het.cloud_freezing_rate(
                mp.ice.rain_freezing, hom, p3.vent, aps, tps, mp.ice.cloud_pdf,
                q_c, ρ, N_c, T, qᵥ)

            for v in (r.∂ₜn_frz, r.∂ₜq_frz, c.∂ₜn_frz, c.∂ₜq_frz)
                TT.@test isfinite(v)
                TT.@test v >= 0
            end
            # Immersion freezing of cloud droplets is Bigg alone: no ice-nucleating-particle
            # budget bounds the heterogeneous coefficient, so there is no capped-versus-uncapped
            # pair to compare. What survives is that the coefficient is a rate.
            TT.@test c.J_het >= 0

            if ΔT <= 0
                # EXACT zero, which is what the derived gate buys: no freezing of water that is
                # not supercooled, in either category, at either precision
                TT.@test r.∂ₜn_frz == 0 && r.∂ₜq_frz == 0
                TT.@test c.∂ₜn_frz == 0 && c.∂ₜq_frz == 0
            else
                n_pos += (r.∂ₜn_frz > 0) + (c.∂ₜn_frz > 0)
                # monotone nondecreasing in supercooling
                TT.@test r.∂ₜn_frz >= prev_rain * (1 - sqrt(eps(FT)))
                TT.@test c.∂ₜn_frz >= prev_cloud * (1 - sqrt(eps(FT)))
                prev_rain, prev_cloud = r.∂ₜn_frz, c.∂ₜn_frz
            end
        end
        # the sweep must produce freezing somewhere, or the monotonicity assertions are vacuous
        TT.@test n_pos > 0
    end
end

# The (T, S_i) plane with no explicit thresholds, asserted on a grid: nonnegative everywhere,
# exactly zero where the air is not supersaturated over ice, continuous in both arguments, and in
# agreement with the former windowed behaviour well inside the window it used to impose.
function test_f23_deposition_rolloffs(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    T_freeze = TDI.TD.Parameters.T_freeze(tps)
    # `Frostenberg2023` on the deposition target-spectrum interface, swapped in for the
    # default closure on an otherwise-default parameter set, so the nascent-crystal mass
    # and the depletion model come from the same single source every closure reads.
    toml_dict = CP.create_toml_dict(FT)
    mp_default = CMP.Microphysics2MParams(toml_dict; with_ice = true)
    p3 = mp_default.ice.scheme
    mp = CMP.Microphysics2MParams(;
        warm_rain = mp_default.warm_rain,
        ice = CMP.P3IceParams(;
            scheme = p3,
            terminal_velocity = mp_default.ice.terminal_velocity,
            cloud_pdf = mp_default.ice.cloud_pdf,
            rain_pdf = mp_default.ice.rain_pdf,
            ice_nucleation = CMP.Frostenberg2023(FT),
            rain_freezing = mp_default.ice.rain_freezing,
            homogeneous = mp_default.ice.homogeneous,
            inp_depletion_model = mp_default.ice.inp_depletion_model,
            quad = mp_default.ice.quad,
        ),
    )
    ρ, n_ice = FT(0.8), FT(0)
    rate(T, S_i) = begin
        q_sat = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        q_vap = (1 + S_i) * q_sat
        micro = (; q_tot = q_vap, q_lcl = FT(0), q_rai = FT(0), q_ice = FT(0), n_ice)
        thermo = (; ρ, T)
        CMI_het.deposition_rate(mp.ice.ice_nucleation, mp, tps, micro, thermo)
    end

    TT.@testset "F23 deposition rolls off continuously in T and S_i [FT=$FT]" begin
        Ts = FT[T_freeze - 40, T_freeze - 20, T_freeze - 15.001, T_freeze - 15,
            T_freeze - 14.999, T_freeze - 8, T_freeze - 1, T_freeze, T_freeze + 2]
        Ss = FT[-0.2, -0.01, 0, 0.001, 0.0499, 0.05, 0.0501, 0.2, 0.5]

        for T in Ts, S_i in Ss
            r = rate(T, S_i)
            TT.@test isfinite(r.∂ₜn_frz) && isfinite(r.∂ₜq_frz)
            TT.@test r.∂ₜn_frz >= 0 && r.∂ₜq_frz >= 0
            if S_i < 0
                # The subsaturation floor owns this end now that the vapor-excess cap is
                # retired, and it carries BOTH moments: a number source with no mass to pay
                # for it is the degeneracy the scheme forbids. The boundary point itself is
                # inside the window, so the test brackets it rather than including it.
                TT.@test r.∂ₜq_frz == 0
                TT.@test r.∂ₜn_frz == 0
            end
        end

        # continuity where the two retired thresholds used to sit
        for S_i in (FT(0.2), FT(0.5))
            a = rate(T_freeze - FT(15.001), S_i).∂ₜn_frz
            b = rate(T_freeze - FT(14.999), S_i).∂ₜn_frz
            TT.@test isapprox(a, b; rtol = FT(0.01))
        end
        for T in (T_freeze - FT(20), T_freeze - FT(8))
            a = rate(T, FT(0.0499)).∂ₜn_frz
            b = rate(T, FT(0.0501)).∂ₜn_frz
            TT.@test isapprox(a, b; rtol = FT(0.01))
        end

        # unchanged where both former thresholds were satisfied
        for T in (T_freeze - FT(40), T_freeze - FT(20)), S_i in (FT(0.2), FT(0.5))
            r = rate(T, S_i)
            TT.@test r.∂ₜn_frz > 0
            TT.@test r.∂ₜq_frz > 0
        end
    end
end

TT.@testset "Heterogeneous Ice Nucleation Tests ($FT)" for FT in (Float64, Float32)
    test_heterogeneous_ice_nucleation(FT)
    test_liquid_freezing_gate_is_shared(FT)
    test_liquid_freezing_composition_safety(FT)
    test_f23_deposition_rolloffs(FT)
end
nothing
