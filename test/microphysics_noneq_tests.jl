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
        TT.@test CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps, FT(0.5 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T) < FT(0)
        TT.@test CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps, FT(1.5 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T) > FT(0)
        TT.@test CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps,          qᵥ_sl,  FT(0), FT(0), FT(0), FT(0), ρ, T) ≈ FT(0)

        TT.@test CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice, tps, FT(0.5 * qᵥ_si), FT(0), FT(0), FT(0), FT(0), ρ, T) < FT(0)
        TT.@test CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice, tps, FT(1.5 * qᵥ_si), FT(0), FT(0), FT(0), FT(0), ρ, T) > FT(0)
        TT.@test CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice, tps,          qᵥ_si,  FT(0), FT(0), FT(0), FT(0), ρ, T) ≈ FT(0)

        # smoke test for values
        TT.@test CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps, FT(1.2 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T) ≈ 3.763045798130144e-5  rtol = 1e-6
        TT.@test CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice,    tps, FT(1.2 * qᵥ_si), FT(0), FT(0), FT(0), FT(0), ρ, T) ≈ 3.235984203087906e-5 rtol = 1e-6

        TT.@test CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps, FT(1.2 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T) ≈ 3.7630474f-5 rtol = 1e-6
        TT.@test CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice,    tps, FT(1.2 * qᵥ_si), FT(0), FT(0), FT(0), FT(0), ρ, T) ≈ 3.2359854f-5 rtol = 1e-6

        # ice grows faster than liquid
        TT.@test CMNe.conv_q_vap_to_q_lcl_icl_MM2015(liquid, tps, FT(1.2 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T) <
                 CMNe.conv_q_vap_to_q_lcl_icl_MM2015(ice,    tps, FT(1.2 * qᵥ_sl), FT(0), FT(0), FT(0), FT(0), ρ, T)

        #! format: on
    end

    TT.@testset "Cloud condensate sedimentation - liquid" begin
        #setup
        ρ = FT(1.1)
        q_lcl = FT(1 * 1e-3)

        #action
        vt_zero = CMNe.terminal_velocity(liquid, Ch2022.rain, ρ, FT(0))
        vt_liq = CMNe.terminal_velocity(liquid, Ch2022.rain, ρ, q_lcl)
        v_bigger = CMNe.terminal_velocity(liquid, Ch2022.rain, ρ, q_lcl * 2)

        #test
        TT.@test vt_zero == FT(0)
        TT.@test vt_liq ≈ 0.004249204292098458 rtol = sqrt(eps(FT))
        TT.@test v_bigger > vt_liq
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
        TT.@test vt_ice ≈ 0.003199067212454443 rtol = sqrt(eps(FT))
        TT.@test v_bigger > vt_ice
    end
end

TT.@testset "Microphysics Non-Equilibrium Tests ($FT)" for FT in (Float64, Float32)
    test_microphysics_noneq(FT)
end
nothing
