import Test as TT

import ClimaParams

import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMNe

function test_microphysics_noneq(FT)

    ice = CMP.CloudIce(FT)
    liquid = CMP.CloudLiquid(FT)
    aps = CMP.AirProperties(FT)
    frs = CMP.Frostenberg2023(FT)
    fit = CMP.IceNumberTemperatureFit(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    Ch2022 = CMP.Chen2022VelType(FT)

    TT.@testset "τ_relax" begin
        # Default τ_relax values live in process_params
        mp = CMP.Microphysics1MParams(FT)
        TT.@test mp.process_params.cloud_liquid_formation.τ_relax ≈ FT(10)
        TT.@test mp.process_params.cloud_ice_formation.τ_relax ≈ FT(10)
    end

    TT.@testset "τ_relax Frostenberg" begin
        q_icl = FT(1e-4)
        T_cold = FT(250)
        ρ_frs = FT(0.8)
        τ_frs = CMNe.τ_relax(ice, aps, frs, q_icl, T_cold, ρ_frs)
        TT.@test τ_frs > FT(0)
        TT.@test isfinite(τ_frs)
        # Increasing ice content should decrease τ (larger crystals → faster relaxation)
        τ_frs_more = CMNe.τ_relax(ice, aps, frs, FT(10) * q_icl, T_cold, ρ_frs)
        TT.@test τ_frs_more < τ_frs
    end

    TT.@testset "τ_relax PrescribedIceNumber" begin
        q_icl = FT(1e-4)
        ρ_pin = FT(0.8)

        # Basic smoke test
        τ_pin = CMNe.τ_relax(ice, aps, q_icl, ρ_pin)
        TT.@test τ_pin > FT(0)
        TT.@test isfinite(τ_pin)

        # Increasing ice content should decrease τ (larger crystals → faster relaxation)
        τ_pin_more = CMNe.τ_relax(ice, aps, FT(10) * q_icl, ρ_pin)
        TT.@test τ_pin_more < τ_pin

        # Higher density (same q) → more ice mass per m³ → larger r
        # Since τ = 1/(4π D_v N_0 r) with N_0 fixed, larger r → smaller τ
        τ_dense = CMNe.τ_relax(ice, aps, q_icl, FT(1.2))
        TT.@test τ_dense < τ_pin
    end

    TT.@testset "ice_number_concentration temperature fit" begin
        T_freeze = fit.T_freeze
        N_0C = CMNe.ice_number_concentration(fit, T_freeze)
        TT.@test N_0C ≈ fit.N_ref * exp(fit.a)   # ~6e1 m⁻³ with the default fit
        TT.@test N_0C isa FT
        # held at the freezing-point value above freezing
        TT.@test CMNe.ice_number_concentration(fit, T_freeze + FT(5)) == N_0C
        # exponential increase toward colder temperatures: a factor exp(10 b) per 10 K
        N_m10 = CMNe.ice_number_concentration(fit, T_freeze - FT(10))
        N_m20 = CMNe.ice_number_concentration(fit, T_freeze - FT(20))
        TT.@test N_0C < N_m10 < N_m20
        TT.@test N_m10 / N_0C ≈ exp(10 * fit.b) rtol = FT(1e-5)
        TT.@test N_m20 / N_m10 ≈ exp(10 * fit.b) rtol = FT(1e-5)
        # cold-end cap
        TT.@test CMNe.ice_number_concentration(fit, FT(200)) == fit.N_max
    end

    TT.@testset "τ_relax TemperatureDependentIceNumber" begin
        q_icl = FT(1e-5)
        ρ = FT(0.8)
        T1, T2 = FT(263), FT(253)
        τ1 = CMNe.τ_relax(ice, aps, fit, q_icl, T1, ρ)
        τ2 = CMNe.τ_relax(ice, aps, fit, q_icl, T2, ρ)
        TT.@test τ1 > FT(0) && isfinite(τ1)
        TT.@test τ1 isa FT
        TT.@test τ2 < τ1   # more crystals at colder T → faster relaxation
        # above the 1 μm radius floor, τ ∝ N^(-2/3)
        N1 = CMNe.ice_number_concentration(fit, T1)
        N2 = CMNe.ice_number_concentration(fit, T2)
        TT.@test τ2 / τ1 ≈ (N2 / N1)^(-FT(2) / 3) rtol = FT(1e-4)
        # identical to the prescribed-number timescale evaluated with N_0 = N_ice(T)
        ice_N1 = CMP.CloudIce(; ice.pdf, ice.mass, ice.ρᵢ, ice.r_eff, N_0 = N1)
        TT.@test τ1 ≈ CMNe.τ_relax(ice_N1, aps, q_icl, ρ)
        # hours at -10 °C for a typical ice content (seconds for N_0 = 5e8 m⁻³)
        TT.@test FT(3600) < τ1 < FT(24 * 3600)
        TT.@test τ1 > FT(100) * CMNe.τ_relax(ice, aps, q_icl, ρ)
    end

    TT.@testset "TemperatureDependentIceNumber conv_q_vap_to_q_icl" begin
        ρ = FT(0.8)
        T = FT(273 - 10)
        qᵥ_si = TDI.p2q(tps, T, ρ, TDI.saturation_vapor_pressure_over_ice(tps, T))
        mp_fit = (; cloud = (; ice), air_properties = aps, process_params = (; cloud_ice_formation = fit))
        opt = CMP.TemperatureDependentIceNumber()

        #! format: off
        _conv(q_tot, q_icl, ρ, T) = CMNe.conv_q_vap_to_q_icl(
            opt, mp_fit, tps,
            (; q_tot, q_lcl = FT(0), q_icl, q_rai = FT(0), q_sno = FT(0)),
            (; ρ, T),
        )
        #! format: on

        # sign tests
        TT.@test _conv(FT(1.5) * qᵥ_si, FT(0), ρ, T) > FT(0)
        TT.@test _conv(FT(0.5) * qᵥ_si, FT(1e-4), ρ, T) < FT(0)
        TT.@test _conv(FT(0.5) * qᵥ_si, FT(0), ρ, T) == FT(0)   # nothing to sublimate
        TT.@test _conv(qᵥ_si, FT(0), ρ, T) ≈ FT(0)
        TT.@test _conv(FT(1.5) * qᵥ_si, FT(1e-5), ρ, T) isa FT

        # equals the constant-τ kernel evaluated at the fit's τ, which the accessor returns
        micro = (; q_tot = FT(1.5) * qᵥ_si, q_lcl = FT(0), q_icl = FT(1e-5), q_rai = FT(0), q_sno = FT(0))
        τ = CMNe.τ_vap_to_q_icl(opt, mp_fit, tps, micro, (; ρ, T))
        TT.@test τ == CMNe.τ_relax(ice, aps, fit, micro.q_icl, T, ρ)
        TT.@test _conv(micro.q_tot, micro.q_icl, ρ, T) ==
                 CMNe._conv_q_vap_to_q_icl_const(τ, tps, micro, (; ρ, T))

        # deposition per unit supersaturation is faster at colder temperatures (more crystals)
        T2 = FT(273 - 20)
        qᵥ_si2 = TDI.p2q(tps, T2, ρ, TDI.saturation_vapor_pressure_over_ice(tps, T2))
        rate_per_excess(T_, q_si_) = _conv(FT(1.5) * q_si_, FT(1e-5), ρ, T_) / (FT(0.5) * q_si_)
        TT.@test rate_per_excess(T2, qᵥ_si2) > rate_per_excess(T, qᵥ_si)

        # INP limiter: no deposition above freezing; sublimation still allowed
        T_warm = FT(280)
        qᵥ_si_w = TDI.p2q(tps, T_warm, ρ, TDI.saturation_vapor_pressure_over_ice(tps, T_warm))
        TT.@test _conv(FT(1.5) * qᵥ_si_w, FT(0), ρ, T_warm) == FT(0)
        TT.@test _conv(FT(0.5) * qᵥ_si_w, FT(1e-3), ρ, T_warm) < FT(0)
    end

    TT.@testset "CondEvap_DepSub" begin



        ρ = FT(0.8)
        T = FT(273 - 10)

        pᵥ_sl = TDI.saturation_vapor_pressure_over_liquid(tps, T)
        qᵥ_sl = TDI.p2q(tps, T, ρ, pᵥ_sl)

        pᵥ_si = TDI.saturation_vapor_pressure_over_ice(tps, T)
        qᵥ_si = TDI.p2q(tps, T, ρ, pᵥ_si)

        #! format: off
        # Test helpers (call the timescale kernels directly with an explicit τ)
        T_hom = FT(233)
        _conv_lcl(q_tot, q_lcl, q_icl, ρ, T) = CMNe._conv_q_vap_to_q_lcl_const(
            FT(10),
            tps,
            (; q_tot, q_lcl, q_icl, q_rai = FT(0), q_sno = FT(0)),
            (; ρ, T);
            T_hom,
        )

        _conv_icl(q_tot, q_lcl, q_icl, ρ, T) = CMNe._conv_q_vap_to_q_icl_const(
            FT(10),
            tps,
            (; q_tot, q_lcl, q_icl, q_rai = FT(0), q_sno = FT(0)),
            (; ρ, T)
        )

        _conv_icl_dep(q_tot, q_lcl, q_icl, ρ, T) = CMNe.conv_q_vap_to_q_icl(
            CMP.TemperatureDependent(),
            (; cloud = (; ice), air_properties = aps,
                process_params = (; cloud_ice_formation = (; τ_relax = FT(10), frostenberg = frs))),
            tps,
            (; q_tot, q_lcl, q_icl, q_rai = FT(0), q_sno = FT(0)),
            (; ρ, T)
        )

        # test sign
        TT.@test _conv_lcl(FT(0.5 * qᵥ_sl), FT(0), FT(0), ρ, T) == FT(0)
        TT.@test _conv_lcl(FT(1.5 * qᵥ_sl), FT(0), FT(0), ρ, T) > FT(0)
        TT.@test _conv_lcl(         qᵥ_sl,  FT(0), FT(0), ρ, T) ≈ FT(0)

        TT.@test _conv_icl(FT(0.5 * qᵥ_si), FT(0), FT(0), ρ, T) == FT(0)
        TT.@test _conv_icl(FT(1.5 * qᵥ_si), FT(0), FT(0), ρ, T) > FT(0)
        TT.@test _conv_icl(         qᵥ_si,  FT(0), FT(0), ρ, T) ≈ FT(0)

        # smoke test for values
        TT.@test _conv_lcl(FT(1.2 * qᵥ_sl), FT(0), FT(0), ρ, T) ≈ 3.763045798130144e-5  rtol = 1e-6
        TT.@test _conv_icl(FT(1.2 * qᵥ_si), FT(0), FT(0), ρ, T) ≈ 3.235984203087906e-5 rtol = 1e-6

        TT.@test _conv_lcl(FT(1.2 * qᵥ_sl), FT(0), FT(0), ρ, T) ≈ 3.7630474f-5 rtol = 1e-6
        TT.@test _conv_icl(FT(1.2 * qᵥ_si), FT(0), FT(0), ρ, T) ≈ 3.2359854f-5 rtol = 1e-6

        # ice grows faster than liquid (same saturation excess w.r.t. liquid)
        TT.@test _conv_lcl(FT(1.2 * qᵥ_sl), FT(0), FT(0), ρ, T) <
                 _conv_icl(FT(1.2 * qᵥ_sl), FT(0), FT(0), ρ, T)

        # --- INP limiter: above freezing, positive ice deposition should be zero ---
        T_warm = FT(280) # above freezing
        pᵥ_si_w = TDI.saturation_vapor_pressure_over_ice(tps, T_warm)
        qᵥ_si_w = TDI.p2q(tps, T_warm, ρ, pᵥ_si_w)
        pᵥ_sl_w = TDI.saturation_vapor_pressure_over_liquid(tps, T_warm)
        qᵥ_sl_w = TDI.p2q(tps, T_warm, ρ, pᵥ_sl_w)

        # Supersaturated w.r.t. ice above freezing → would deposit, but INP limiter zeros it
        TT.@test _conv_icl(FT(1.5 * qᵥ_si_w), FT(0), FT(0), ρ, T_warm) == FT(0)
        # Liquid condensation above freezing should be unaffected
        TT.@test _conv_lcl(FT(1.5 * qᵥ_sl_w), FT(0), FT(0), ρ, T_warm) > FT(0)
        # Ice sublimation (negative tendency) above freezing should NOT be limited
        TT.@test _conv_icl(FT(0.5 * qᵥ_si_w), FT(0), FT(0.001), ρ, T_warm) <= FT(0)

        # --- Homogeneous limiter: below T_hom, positive liquid condensation should be zero ---
        T_cold = FT(220) # well below T_hom ≈ 233 K
        pᵥ_sl_c = TDI.saturation_vapor_pressure_over_liquid(tps, T_cold)
        qᵥ_sl_c = TDI.p2q(tps, T_cold, ρ, pᵥ_sl_c)
        # Supersaturated w.r.t. liquid below T_hom → would condense, but limiter zeros it
        TT.@test _conv_lcl(FT(1.5 * qᵥ_sl_c), FT(0), FT(0), ρ, T_cold) == FT(0)
        # Liquid evaporation below T_hom should NOT be limited
        TT.@test _conv_lcl(FT(0.5 * qᵥ_sl_c), FT(0.001), FT(0), ρ, T_cold) < FT(0)

        # --- Asymmetric τ_dep ≠ τ_sub ---
        # Faster sublimation (smaller τ_relax) should give a larger magnitude sublimation rate
        # We test by calling the constant-timescale kernel with different τ values
        function _conv_icl_custom(τ, q_tot, q_lcl, q_icl, ρ, T)
            CMNe._conv_q_vap_to_q_icl_const(
                τ,
                tps,
                (; q_tot, q_lcl, q_icl, q_rai = FT(0), q_sno = FT(0)),
                (; ρ, T)
            )
        end
        sub_fast = _conv_icl_custom(FT(1), FT(0.5 * qᵥ_si), FT(0), FT(0.001), ρ, T)
        sub_slow = _conv_icl_custom(FT(100), FT(0.5 * qᵥ_si), FT(0), FT(0.001), ρ, T)
        TT.@test sub_fast < sub_slow  # faster sublimation → more negative

        # Deposition rate with TemperatureDependent should depend only on τ_dep, not τ_sub
        # (τ_sub is only used for sublimation, deposition uses the Frostenberg timescale)
        function _conv_icl_dep_custom(τ_sub, q_tot, q_lcl, q_icl, ρ, T)
            CMNe.conv_q_vap_to_q_icl(
                CMP.TemperatureDependent(),
                (; cloud = (; ice), air_properties = aps,
                    process_params = (; cloud_ice_formation = (; τ_relax = τ_sub, frostenberg = frs))),
                tps,
                (; q_tot, q_lcl, q_icl, q_rai = FT(0), q_sno = FT(0)),
                (; ρ, T)
            )
        end
        dep_a = _conv_icl_dep_custom(FT(1), FT(1.5 * qᵥ_si), FT(0), FT(0), ρ, T)
        dep_b = _conv_icl_dep_custom(FT(100), FT(1.5 * qᵥ_si), FT(0), FT(0), ρ, T)
        TT.@test dep_a ≈ dep_b  # τ_sub doesn't affect deposition

        #! format: on
    end

    TT.@testset "PrescribedIceNumber conv_q_vap_to_q_icl" begin
        ρ = FT(0.8)
        T = FT(273 - 10)

        pᵥ_si = TDI.saturation_vapor_pressure_over_ice(tps, T)
        qᵥ_si = TDI.p2q(tps, T, ρ, pᵥ_si)

        #! format: off
        _conv_pin(q_tot, q_lcl, q_icl, ρ, T) = CMNe.conv_q_vap_to_q_icl(
            CMP.PrescribedIceNumber(),
            (; cloud = (; ice), air_properties = aps, process_params = (;)),
            tps,
            (; q_tot, q_lcl, q_icl, q_rai = FT(0), q_sno = FT(0)),
            (; ρ, T)
        )
        #! format: on

        # Sign tests
        TT.@test _conv_pin(FT(0.5 * qᵥ_si), FT(0), FT(0), ρ, T) == FT(0)
        TT.@test _conv_pin(FT(1.5 * qᵥ_si), FT(0), FT(0), ρ, T) > FT(0)
        TT.@test _conv_pin(qᵥ_si, FT(0), FT(0), ρ, T) ≈ FT(0)

        # INP limiter: above freezing, positive ice deposition should be zero
        T_warm = FT(280)
        pᵥ_si_w = TDI.saturation_vapor_pressure_over_ice(tps, T_warm)
        qᵥ_si_w = TDI.p2q(tps, T_warm, ρ, pᵥ_si_w)
        TT.@test _conv_pin(FT(1.5 * qᵥ_si_w), FT(0), FT(0), ρ, T_warm) == FT(0)

        # Ice sublimation above freezing should NOT be limited
        TT.@test _conv_pin(FT(0.5 * qᵥ_si_w), FT(0), FT(0.001), ρ, T_warm) <= FT(0)

        # The tendency should differ from a constant-timescale result (τ depends on q_icl, ρ)
        rate_a = _conv_pin(FT(1.5 * qᵥ_si), FT(0), FT(1e-5), ρ, T)
        rate_b = _conv_pin(FT(1.5 * qᵥ_si), FT(0), FT(1e-3), ρ, T)
        TT.@test rate_a != rate_b  # timescale changes with q_icl
    end



    TT.@testset "Cloud condensate sedimentation - liquid" begin
        #setup
        ρ = FT(1.1)
        q_lcl = FT(1 * 1e-3)
        stokes_vel = CMP.StokesRegimeVelType(FT)

        #action
        vt_zero = CMNe.terminal_velocity(liquid, stokes_vel, ρ, FT(0))
        vt_liq = CMNe.terminal_velocity(liquid, stokes_vel, ρ, q_lcl)
        v_double = CMNe.terminal_velocity(liquid, stokes_vel, ρ, q_lcl * 2)

        #test
        TT.@test vt_zero == FT(0)
        TT.@test vt_liq > FT(0)

        # Test Stokes law scaling: v ∝ D² and D ∝ (q/N)^(1/3), so v ∝ q^(2/3)
        # Doubling q should scale velocity by 2^(2/3) ≈ 1.587
        expected_ratio = FT(2)^(FT(2) / 3)
        TT.@test v_double / vt_liq ≈ expected_ratio rtol = 1e-6
    end

    TT.@testset "Cloud condensate sedimentation - ice" begin
        #setup
        ρ = FT(0.75)
        q_icl = FT(0.5 * 1e-3)

        #action
        vt_zero = CMNe.terminal_velocity(ice, Ch2022.small_ice, ρ, FT(0))
        vt_ice = CMNe.terminal_velocity(ice, Ch2022.small_ice, ρ, q_icl)
        v_bigger = CMNe.terminal_velocity(ice, Ch2022.small_ice, ρ, q_icl * 2)

        #test
        TT.@test vt_zero == FT(0)
        TT.@test vt_ice > FT(0)  # Updated regression: 6/π factor added
        TT.@test v_bigger > vt_ice
    end
end

TT.@testset "Microphysics Non-Equilibrium Tests ($FT)" for FT in (Float64, Float32)
    test_microphysics_noneq(FT)
end
nothing
