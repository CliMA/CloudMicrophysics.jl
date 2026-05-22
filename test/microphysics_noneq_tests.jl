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
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    Ch2022 = CMP.Chen2022VelType(FT)

    TT.@testset "τ_relax" begin
        TT.@test CMNe.τ_relax(liquid) ≈ FT(10)
        TT.@test CMNe.τ_relax(ice) ≈ FT(10)
    end

    TT.@testset "τ_relax Frostenberg" begin
        q_icl = FT(1e-4)
        T_cold = FT(250)
        τ_frs = CMNe.τ_relax(ice, aps, frs, q_icl, T_cold)
        TT.@test τ_frs > FT(0)
        TT.@test isfinite(τ_frs)
        # Increasing ice content should decrease τ (larger crystals → faster relaxation)
        τ_frs_more = CMNe.τ_relax(ice, aps, frs, FT(10) * q_icl, T_cold)
        TT.@test τ_frs_more < τ_frs
    end

    TT.@testset "CondEvap_DepSub" begin

        τ_liq = CMNe.τ_relax(liquid)
        τ_ice = CMNe.τ_relax(ice)

        ρ = FT(0.8)
        T = FT(273 - 10)

        pᵥ_sl = TDI.saturation_vapor_pressure_over_liquid(tps, T)
        qᵥ_sl = TDI.p2q(tps, T, ρ, pᵥ_sl)

        pᵥ_si = TDI.saturation_vapor_pressure_over_ice(tps, T)
        qᵥ_si = TDI.p2q(tps, T, ρ, pᵥ_si)

        #! format: off
        # Test helpers
        _conv_lcl(q_tot, q_lcl, q_icl, ρ, T) = CMNe.conv_q_vap_to_q_lcl(
            CMP.ConstantTimescaleCloudLiquidFormation(),
            (; cloud = (; liquid = liquid)),
            tps,
            (; q_tot, q_lcl, q_icl, q_rai = FT(0), q_sno = FT(0)),
            (; ρ, T)
        )

        _conv_icl(q_tot, q_lcl, q_icl, ρ, T) = CMNe.conv_q_vap_to_q_icl(
            CMP.ConstantTimescaleCloudIceFormation(),
            (; cloud = (; ice = ice)),
            tps,
            (; q_tot, q_lcl, q_icl, q_rai = FT(0), q_sno = FT(0)),
            (; ρ, T)
        )

        _conv_icl_dep(q_tot, q_lcl, q_icl, ρ, T) = CMNe.conv_q_vap_to_q_icl(
            CMP.TemperatureDependentCloudIceFormation(),
            (; cloud = (; ice = ice), air_properties = aps, frostenberg2023 = frs),
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

        # --- Asymmetric τ_dep ≠ τ_sub ---
        # Faster sublimation (smaller τ_sub) should give a larger magnitude sublimation rate
        # We simulate this by changing the mp ice struct
        ice_base = CMP.CloudIce(FT)
        ice_slow = CMP.CloudIce(ice_base.pdf, ice_base.mass, ice_base.r0, ice_base.r_ice_snow, FT(100), ice_base.ρᵢ, ice_base.r_eff, ice_base.N_0)
        ice_fast = CMP.CloudIce(ice_base.pdf, ice_base.mass, ice_base.r0, ice_base.r_ice_snow, FT(1), ice_base.ρᵢ, ice_base.r_eff, ice_base.N_0)
        
        function _conv_icl_custom(ice_mock, q_tot, q_lcl, q_icl, ρ, T)
            CMNe.conv_q_vap_to_q_icl(
                CMP.ConstantTimescaleCloudIceFormation(),
                (; cloud = (; ice = ice_mock)),
                tps,
                (; q_tot, q_lcl, q_icl, q_rai = FT(0), q_sno = FT(0)),
                (; ρ, T)
            )
        end
        sub_fast = _conv_icl_custom(ice_fast, FT(0.5 * qᵥ_si), FT(0), FT(0.001), ρ, T)
        sub_slow = _conv_icl_custom(ice_slow, FT(0.5 * qᵥ_si), FT(0), FT(0.001), ρ, T)
        TT.@test sub_fast < sub_slow  # faster sublimation → more negative

        # Deposition rate with TemperatureDependentCloudIceFormation should depend only on τ_dep, not τ_sub
        function _conv_icl_dep_custom(ice_mock, q_tot, q_lcl, q_icl, ρ, T)
            CMNe.conv_q_vap_to_q_icl(
                CMP.TemperatureDependentCloudIceFormation(),
                (; cloud = (; ice = ice_mock), air_properties = aps, frostenberg2023 = frs),
                tps,
                (; q_tot, q_lcl, q_icl, q_rai = FT(0), q_sno = FT(0)),
                (; ρ, T)
            )
        end
        dep_a = _conv_icl_dep_custom(ice_fast, FT(1.5 * qᵥ_si), FT(0), FT(0), ρ, T)
        dep_b = _conv_icl_dep_custom(ice_slow, FT(1.5 * qᵥ_si), FT(0), FT(0), ρ, T)
        TT.@test dep_a ≈ dep_b  # τ_sub doesn't affect deposition

        #! format: on
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
