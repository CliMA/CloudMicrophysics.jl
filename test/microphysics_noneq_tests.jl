import Test as TT

import ClimaParams

import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMNe

function test_microphysics_noneq(FT)

    ice = CMP.CloudIce(FT)
    liquid = CMP.CloudLiquid(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    Ch2022 = CMP.Chen2022VelType(FT)

    TT.@testset "τ_relax" begin
        TT.@test CMNe.τ_relax(liquid) ≈ FT(10)
        TT.@test CMNe.τ_relax(ice) ≈ FT(10)
    end

    TT.@testset "CloudLiquidCondEvap" begin

        q_liq_sat = FT(5e-3)
        frac = [FT(0), FT(0.5), FT(1), FT(1.5)]

        _τ_cond_evap = CMNe.τ_relax(liquid)

        for fr in frac
            q_lcl = q_liq_sat * fr
            TT.@test CMNe.conv_q_vap_to_q_lcl_icl(liquid, q_liq_sat, q_lcl) ≈ (1 - fr) * q_liq_sat / _τ_cond_evap
        end
    end

    TT.@testset "CloudIceCondEvap" begin

        q_ice_sat = FT(2e-3)
        frac = [FT(0), FT(0.5), FT(1), FT(1.5)]

        _τ_cond_evap = CMNe.τ_relax(ice)

        for fr in frac
            q_icl = q_ice_sat * fr
            TT.@test CMNe.conv_q_vap_to_q_lcl_icl(ice, q_ice_sat, q_icl) ≈ (1 - fr) * q_ice_sat / _τ_cond_evap
        end
    end

    TT.@testset "CondEvap_DepSub_MM2015" begin

        ρ = FT(0.8)
        T = FT(273 - 10)

        pᵥ_sl = TDI.saturation_vapor_pressure_over_liquid(tps, T)
        qᵥ_sl = TDI.p2q(tps, T, ρ, pᵥ_sl)

        pᵥ_si = TDI.saturation_vapor_pressure_over_ice(tps, T)
        qᵥ_si = TDI.p2q(tps, T, ρ, pᵥ_si)

        #! format: off
        # test sign
        TT.@test first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps, FT(0.5 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T)) == FT(0)
        TT.@test first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps, FT(1.5 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T)) > FT(0)
        TT.@test first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps,          qᵥ_sl,  FT(0), FT(0), FT(0), FT(0), ρ, T)) ≈ FT(0)

        TT.@test first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice, tps, FT(0.5 * qᵥ_si), FT(0), FT(0), FT(0), FT(0), ρ, T)) == FT(0)
        TT.@test first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice, tps, FT(1.5 * qᵥ_si), FT(0), FT(0), FT(0), FT(0), ρ, T)) > FT(0)
        TT.@test first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice, tps,          qᵥ_si,  FT(0), FT(0), FT(0), FT(0), ρ, T)) ≈ FT(0)

        # smoke test for values
        TT.@test first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps, FT(1.2 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T)) ≈ 3.763045798130144e-5  rtol = 1e-6
        TT.@test first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice,    tps, FT(1.2 * qᵥ_si), FT(0), FT(0), FT(0), FT(0), ρ, T)) ≈ 3.235984203087906e-5 rtol = 1e-6

        TT.@test first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps, FT(1.2 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T)) ≈ 3.7630474f-5 rtol = 1e-6
        TT.@test first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice,    tps, FT(1.2 * qᵥ_si), FT(0), FT(0), FT(0), FT(0), ρ, T)) ≈ 3.2359854f-5 rtol = 1e-6

        # ice grows faster than liquid
        TT.@test first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps, FT(1.2 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T)) <
                 first(CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice,    tps, FT(1.2 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T))

        #! format: on
    end

    TT.@testset "CondEvap_DepSub_MM2015_derivatives" begin

        ρ = FT(0.8)
        T = FT(273 - 10)
        q_tot = FT(0.01)   # 10 g/kg total water
        q_lcl = FT(0.001)  # 1 g/kg cloud liquid
        q_icl = FT(0.0005) # 0.5 g/kg cloud ice
        q_rai = FT(0)
        q_sno = FT(0)
        τ = CMNe.τ_relax(liquid)

        # --- Liquid condensation derivative ---
        (S_liq, dS_liq_simple, dS_liq) = CMNe.conv_q_vap_to_q_lcl_icl_MM2015(
            liquid, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T,
        )

        # derivative should be negative (increasing q_lcl decreases tendency)
        TT.@test dS_liq < FT(0)
        TT.@test dS_liq_simple < FT(0)

        # simple derivative should be -1/τ
        TT.@test dS_liq_simple ≈ FT(-1) / τ

        # full derivative should be close to but not equal to simple
        TT.@test dS_liq ≈ dS_liq_simple rtol = 0.5
        TT.@test dS_liq != dS_liq_simple

        # --- Ice deposition derivative ---
        (S_ice, dS_ice_simple, dS_ice) = CMNe.conv_q_vap_to_q_lcl_icl_MM2015(
            ice, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T,
        )

        TT.@test dS_ice < FT(0)
        TT.@test dS_ice_simple ≈ FT(-1) / CMNe.τ_relax(ice)

        # --- Finite difference validation (Float64 only) ---
        # The analytical derivative is a total derivative: dS/dq
        # accounting for implicit temperature feedback: dT/dq = L/cₚ.
        # The FD must perturb both q and T simultaneously.
        # We use rtol=0.01 because the analytical formula neglects
        # higher-order dL/dT and dcₚ/dT terms.
        if FT == Float64
            ε = FT(1e-8)

            # Validate total dS/dq_lcl for liquid
            Lᵥ = TDI.Lᵥ(tps, T)
            cₚ = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)
            dT_dq_liq = Lᵥ / cₚ

            (S_plus, _, _) = CMNe.conv_q_vap_to_q_lcl_icl_MM2015(
                liquid, tps, q_tot, q_lcl + ε, q_icl, q_rai, q_sno, ρ, T + ε * dT_dq_liq,
            )
            (S_minus, _, _) = CMNe.conv_q_vap_to_q_lcl_icl_MM2015(
                liquid, tps, q_tot, q_lcl - ε, q_icl, q_rai, q_sno, ρ, T - ε * dT_dq_liq,
            )
            dS_fd_liq = (S_plus - S_minus) / (2ε)
            TT.@test dS_liq ≈ dS_fd_liq rtol = 0.01

            # Validate total dS/dq_icl for ice
            Lₛ = TDI.Lₛ(tps, T)
            dT_dq_ice = Lₛ / cₚ

            (S_plus_i, _, _) = CMNe.conv_q_vap_to_q_lcl_icl_MM2015(
                ice, tps, q_tot, q_lcl, q_icl + ε, q_rai, q_sno, ρ, T + ε * dT_dq_ice,
            )
            (S_minus_i, _, _) = CMNe.conv_q_vap_to_q_lcl_icl_MM2015(
                ice, tps, q_tot, q_lcl, q_icl - ε, q_rai, q_sno, ρ, T - ε * dT_dq_ice,
            )
            dS_fd_ice = (S_plus_i - S_minus_i) / (2ε)
            TT.@test dS_ice ≈ dS_fd_ice rtol = 0.01
        end

        # --- Limiter: when tendency < 0 and q <= 0, all outputs should be zero ---
        (S_z, dS_z_s, dS_z) = CMNe.conv_q_vap_to_q_lcl_icl_MM2015(
            liquid, tps, FT(0.5 * 0.002), FT(0), FT(0), FT(0), FT(0), ρ, T,
        )
        TT.@test S_z == FT(0)
        TT.@test dS_z == FT(0)
        TT.@test dS_z_s == FT(0)
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
