import Test as TT

import Thermodynamics as TD

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMNe

@info "Non-equilibrium Microphysics Tests"
function test_microphysics_noneq(FT)

    ice = CMP.CloudIce(FT)
    liquid = CMP.CloudLiquid(FT)

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
end

test_microphysics_noneq(Float32)
test_microphysics_noneq(Float64)
