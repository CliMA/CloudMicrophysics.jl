import Test as TT

import ClimaParams
import Thermodynamics as TD

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMNe

@info "Non-equilibrium Microphysics Tests"
function test_microphysics_noneq(FT)

    ice = CMP.CloudIce(FT)
    liquid = CMP.CloudLiquid(FT)
    tps = TD.Parameters.ThermodynamicsParameters(FT)

    TT.@testset "τ_relax" begin
        TT.@test CMNe.τ_relax(liquid) ≈ FT(10)
        TT.@test CMNe.τ_relax(ice) ≈ FT(10)
    end

    TT.@testset "CloudLiquidCondEvap" begin

        q_liq_sat = FT(5e-3)
        frac = [FT(0), FT(0.5), FT(1), FT(1.5)]

        _τ_cond_evap = CMNe.τ_relax(liquid)

        for fr in frac
            q_liq = q_liq_sat * fr

            TT.@test CMNe.conv_q_vap_to_q_liq_ice(
                liquid,
                TD.PhasePartition(FT(0), q_liq_sat, FT(0)),
                TD.PhasePartition(FT(0), q_liq, FT(0)),
            ) ≈ (1 - fr) * q_liq_sat / _τ_cond_evap
        end
    end

    TT.@testset "CloudIceCondEvap" begin

        q_ice_sat = FT(2e-3)
        frac = [FT(0), FT(0.5), FT(1), FT(1.5)]

        _τ_cond_evap = CMNe.τ_relax(ice)

        for fr in frac
            q_ice = q_ice_sat * fr

            TT.@test CMNe.conv_q_vap_to_q_liq_ice(
                ice,
                TD.PhasePartition(FT(0), FT(0), q_ice_sat),
                TD.PhasePartition(FT(0), FT(0), q_ice),
            ) ≈ (1 - fr) * q_ice_sat / _τ_cond_evap
        end
    end

    TT.@testset "CondEvap_DepSub_MM2015" begin

        ρ = FT(0.8)
        T = FT(273 - 10)

        pᵥ_sl = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
        qᵥ_sl = TD.q_vap_saturation_from_density(tps, T, ρ, pᵥ_sl)

        pᵥ_si = TD.saturation_vapor_pressure(tps, T, TD.Ice())
        qᵥ_si = TD.q_vap_saturation_from_density(tps, T, ρ, pᵥ_si)

        qᵥ_s = TD.PhasePartition(FT(0), qᵥ_sl, qᵥ_si)

        qₚ(qᵥ) = TD.PhasePartition(FT(qᵥ))

        #! format: off
        # test sign
        TT.@test CMNe.conv_q_vap_to_q_liq_ice_MM2015(liquid, tps, qᵥ_s, qₚ(FT(0.5 * qᵥ_sl)), T) < FT(0)
        TT.@test CMNe.conv_q_vap_to_q_liq_ice_MM2015(liquid, tps, qᵥ_s, qₚ(FT(1.5 * qᵥ_sl)), T) > FT(0)
        TT.@test CMNe.conv_q_vap_to_q_liq_ice_MM2015(liquid, tps, qᵥ_s, qₚ(      qᵥ_sl), T) ≈ FT(0)

        TT.@test CMNe.conv_q_vap_to_q_liq_ice_MM2015(ice, tps, qᵥ_s, qₚ(FT(0.5 * qᵥ_si)), T) < FT(0)
        TT.@test CMNe.conv_q_vap_to_q_liq_ice_MM2015(ice, tps, qᵥ_s, qₚ(FT(1.5 * qᵥ_si)), T) > FT(0)
        TT.@test CMNe.conv_q_vap_to_q_liq_ice_MM2015(ice, tps, qᵥ_s, qₚ(      qᵥ_si), T) ≈ FT(0)

        # smoke test for values
        TT.@test CMNe.conv_q_vap_to_q_liq_ice_MM2015(liquid, tps, qᵥ_s, qₚ(FT(1.2 * qᵥ_sl)), T) ≈ 3.7631f-5 rtol = 1e-6
        TT.@test CMNe.conv_q_vap_to_q_liq_ice_MM2015(ice,    tps, qᵥ_s, qₚ(FT(1.2 * qᵥ_si)), T) ≈ 3.2356777f-5 rtol = 1e-6

        # ice grows faster than liquid
        TT.@test CMNe.conv_q_vap_to_q_liq_ice_MM2015(liquid, tps, qᵥ_s, qₚ(FT(1.2 * qᵥ_sl)), T) <
                 CMNe.conv_q_vap_to_q_liq_ice_MM2015(ice,    tps, qᵥ_s, qₚ(FT(1.2 * qᵥ_sl)), T)

        #! format: on
    end
end

test_microphysics_noneq(Float32)
test_microphysics_noneq(Float64)
