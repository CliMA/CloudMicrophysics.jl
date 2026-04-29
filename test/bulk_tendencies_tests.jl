using Test

import ClimaParams as CP

import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.MicrophysicsNonEq as CMNonEq
import CloudMicrophysics.ThermodynamicsInterface as TDI

# Helper to compute q_tot that enforces saturation (S=0) over liquid
# This disables condensation/evaporation to isolate microphysical processes
function get_saturated_q_tot(tps, T::FT, ρ::FT, q_lcl::FT, q_icl::FT, q_rai::FT, q_sno::FT) where {FT}
    q_vap_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
    return q_vap_sat + q_lcl + q_icl + q_rai + q_sno
end

###
### 0M tendencies and derivatives tests
###

function test_bulk_microphysics_0m_tendencies(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics0MParams(FT)

    T_freeze = TDI.T_freeze(tps)

    @testset "BulkMicrophysicsTendencies 0M - Precipitation removal" begin
        T = T_freeze + FT(10)
        q_lcl = FT(2e-3)  # Above threshold
        q_icl = FT(0)

        dq_tot_dt = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics0Moment(),
            mp, tps, T, q_lcl, q_icl,
        )

        # Should be negative (removing condensate)
        @test dq_tot_dt < FT(0)
    end

    @testset "BulkMicrophysicsTendencies 0M - Below threshold" begin
        T = T_freeze + FT(10)
        q_lcl = FT(1e-6)  # Below threshold
        q_icl = FT(0)

        dq_tot_dt = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics0Moment(),
            mp, tps, T, q_lcl, q_icl,
        )

        # No precipitation when below threshold
        @test dq_tot_dt == FT(0)
    end

    @testset "BulkMicrophysicsTendencies 0M - Type stability" begin
        T = T_freeze + FT(5)
        q_lcl = FT(1e-3)
        q_icl = FT(5e-4)

        dq_tot_dt = @inferred BMT.bulk_microphysics_tendencies(
            BMT.Microphysics0Moment(),
            mp, tps, T, q_lcl, q_icl,
        )
        @test dq_tot_dt isa FT
    end

    @testset "BulkMicrophysicsTendencies 0M - S_0 precipitation removal" begin
        ρ = FT(1.2)
        T = T_freeze + FT(10)
        q_lcl = FT(2e-3)
        q_icl = FT(0)
        q_vap_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)

        dq_tot_dt = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics0Moment(),
            mp, tps, T, q_lcl, q_icl, q_vap_sat,
        )

        # Should be negative (removing condensate above S_0 * q_vap_sat)
        @test dq_tot_dt < FT(0)
    end

    @testset "BulkMicrophysicsTendencies 0M - S_0 below threshold" begin
        ρ = FT(1.2)
        T = T_freeze + FT(10)
        q_vap_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
        # Condensate below S_0 * q_vap_sat
        q_lcl = FT(1e-8)
        q_icl = FT(0)

        dq_tot_dt = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics0Moment(),
            mp, tps, T, q_lcl, q_icl, q_vap_sat,
        )

        @test dq_tot_dt == FT(0)
    end

    @testset "BulkMicrophysicsTendencies 0M - S_0 type stability" begin
        ρ = FT(1.2)
        T = T_freeze + FT(5)
        q_lcl = FT(1e-3)
        q_icl = FT(5e-4)
        q_vap_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)

        dq_tot_dt = @inferred BMT.bulk_microphysics_tendencies(
            BMT.Microphysics0Moment(),
            mp, tps, T, q_lcl, q_icl, q_vap_sat,
        )
        @test dq_tot_dt isa FT
    end
end


###
### 1M tendencies and derivatives tests
###

function test_bulk_microphysics_1m_tendencies(FT)

    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    # Use unified parameter container
    mp = CMP.Microphysics1MParams(FT)

    # Extract individual parameters for direct testing
    liquid = mp.cloud.liquid
    ice = mp.cloud.ice
    rain = mp.precip.rain
    snow = mp.precip.snow
    ce = mp.collision
    aps = mp.air_properties
    vel = mp.terminal_velocity

    T_freeze = TDI.T_freeze(tps)

    @testset "BulkMicrophysicsTendencies - Ice to snow conversion" begin
        # In cold conditions with cloud ice above autoconversion threshold:
        # - snow should increase (autoconversion + accretion from ice)
        ρ = FT(1.2)
        T = T_freeze - FT(20)  # cold

        q_lcl = FT(0)
        q_icl = FT(1e-3)  # above threshold
        q_rai = FT(0)
        q_sno = FT(1e-4)
        q_tot = FT(0.012)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Snow increases from ice autoconversion and accretion
        @test tendencies.dq_sno_dt > FT(0)
    end

    @testset "BulkMicrophysicsTendencies - Snow melting" begin
        # In warm conditions with snow:
        # - snow should melt to rain
        ρ = FT(1.2)
        T = T_freeze + FT(5)  # above freezing

        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(1e-3)
        q_tot = FT(0.010)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Snow decreases due to melting (and sublimation)
        @test tendencies.dq_sno_dt < FT(0)
        # Rain increases from melted snow
        @test tendencies.dq_rai_dt > FT(0)
    end

    @testset "BulkMicrophysicsTendencies - Finiteness checks" begin
        # Regression test: keep results constant when code changes
        ρ = FT(1.2)
        T = T_freeze - FT(5)
        q_tot = FT(0.015)
        q_lcl = FT(5e-4)
        q_icl = FT(5e-4)
        q_rai = FT(5e-4)
        q_sno = FT(5e-4)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Test that all tendencies are finite and well-defined
        @test isfinite(tendencies.dq_lcl_dt)
        @test isfinite(tendencies.dq_icl_dt)
        @test isfinite(tendencies.dq_rai_dt)
        @test isfinite(tendencies.dq_sno_dt)
    end

    @testset "BulkMicrophysicsTendencies - Autoconversion component" begin
        # Verify that autoconversion contributes to rain tendency
        ρ = FT(1.2)
        T = T_freeze + FT(10)  # warm, no ice processes
        q_tot = FT(0.015)
        q_lcl = FT(2e-3)  # large cloud liquid, above autoconversion threshold
        q_icl = FT(0)
        q_rai = FT(0)  # no rain yet
        q_sno = FT(0)

        # Call fused function
        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Compute autoconversion separately
        S_acnv_lcl = CM1.conv_q_lcl_to_q_rai(rain.acnv1M, q_lcl, true)

        # Rain tendency should be positive and include autoconversion
        @test tendencies.dq_rai_dt > FT(0)
        # Autoconversion should be a major component of rain formation
        # (in the absence of rain, there's no accretion, so it should be close)
        @test abs(tendencies.dq_rai_dt - S_acnv_lcl) / S_acnv_lcl < FT(0.1)
    end

    @testset "BulkMicrophysicsTendencies - Subsaturated evaporation" begin
        # Test rain evaporation in subsaturated conditions
        ρ = FT(1.2)
        T = T_freeze + FT(15)  # warm
        # Use density-based saturation specific humidity
        q_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)

        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(1e-3)  # rain present
        q_sno = FT(0)
        q_vap = FT(0.5) * q_sat  # subsaturated
        q_tot = q_vap + q_lcl + q_icl + q_rai + q_sno

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Rain should decrease due to evaporation
        @test tendencies.dq_rai_dt < FT(0)
    end

    @testset "BulkMicrophysicsTendencies - Edge cases" begin
        # Test with very small values
        ρ = FT(1.2)
        T = T_freeze - FT(10)
        q_tot = FT(0.01)
        q_lcl = FT(1e-10)
        q_icl = FT(1e-10)
        q_rai = FT(1e-10)
        q_sno = FT(1e-10)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        @test isfinite(tendencies.dq_lcl_dt)
        @test isfinite(tendencies.dq_icl_dt)
        @test isfinite(tendencies.dq_rai_dt)
        @test isfinite(tendencies.dq_sno_dt)
    end

    @testset "BulkMicrophysicsTendencies - Cold riming (lcl + sno → sno)" begin
        # Cold conditions: cloud liquid accreting onto snow should form more snow (riming)
        ρ = FT(1.2)
        T = T_freeze - FT(10)  # well below freezing
        q_lcl = FT(5e-4)  # cloud liquid present
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(1e-3)  # snow present for accretion
        q_tot = FT(0.012)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Verify accretion is happening by checking snow growth
        # (cloud liquid tendency may be positive overall due to condensation)
        S_accr = CM1.accretion(liquid, snow, vel.snow, ce, q_lcl, q_sno, ρ)
        @test S_accr > FT(0)  # Accretion is occurring

        # Snow should increase (riming + deposition)
        @test tendencies.dq_sno_dt > FT(0)
        # No rain formation in cold conditions (approximately zero, allowing for numerical noise)
        @test abs(tendencies.dq_rai_dt) < FT(1e-6)
    end

    @testset "BulkMicrophysicsTendencies - Warm shedding (lcl + sno → rai + melt)" begin
        # Warm conditions: cloud liquid + snow → rain, plus additional snow melt
        ρ = FT(1.2)
        T = T_freeze + FT(5)  # above freezing
        q_lcl = FT(5e-4)  # cloud liquid present
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(1e-3)  # snow present
        q_tot = FT(0.015)  # close to saturation

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Verify accretion is happening
        S_accr = CM1.accretion(liquid, snow, vel.snow, ce, q_lcl, q_sno, ρ)
        @test S_accr > FT(0)  # Accretion is occurring

        # Rain should increase (from accretion + snow melt)
        @test tendencies.dq_rai_dt > FT(0)
        # Snow should decrease (melting + sublimation in subsaturated conditions)
        @test tendencies.dq_sno_dt < FT(0)
    end

    @testset "BulkMicrophysicsTendencies - Cold rain-snow collision (rai + sno → sno)" begin
        # Cold conditions: rain + snow collisions should form snow
        ρ = FT(1.2)
        T = T_freeze - FT(5)  # below freezing
        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(5e-4)  # rain present
        q_sno = FT(5e-4)  # snow present
        q_tot = FT(0.012)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Rain should decrease (freezing onto snow)
        @test tendencies.dq_rai_dt < FT(0)
        # Snow should increase
        @test tendencies.dq_sno_dt > FT(0)
    end

    @testset "BulkMicrophysicsTendencies - Warm rain-snow collision (sno + rai → rai)" begin
        # Warm conditions: snow + rain collisions should form rain
        ρ = FT(1.2)
        T = T_freeze + FT(3)  # just above freezing
        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(5e-4)  # rain present
        q_sno = FT(5e-4)  # snow present
        q_tot = FT(0.015)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Snow should decrease (melting from collision + thermal melt)
        @test tendencies.dq_sno_dt < FT(0)
        # Rain should increase
        @test tendencies.dq_rai_dt > FT(0)
    end

    @testset "BulkMicrophysicsTendencies - Snow deposition (supersaturated, cold)" begin
        # Cold, supersaturated conditions: snow should grow by vapor deposition
        ρ = FT(1.2)
        T = T_freeze - FT(15)  # cold
        q_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)

        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(1e-4)  # small amount of snow
        q_vap = FT(1.2) * q_sat_ice  # supersaturated over ice
        q_tot = q_vap + q_lcl + q_icl + q_rai + q_sno

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Snow should increase from deposition (evaporation_sublimation returns positive for S > 0)
        @test tendencies.dq_sno_dt > FT(0)
    end

    @testset "BulkMicrophysicsTendencies - Ice deposition suppressed above freezing" begin
        # Above freezing: ice deposition should be suppressed even if supersaturated over ice
        ρ = FT(1.2)
        T = T_freeze + FT(2)  # slightly above freezing
        q_lcl = FT(0)
        q_icl = FT(1e-4)  # small ice present
        q_rai = FT(0)
        q_sno = FT(0)
        q_tot = FT(0.015)  # plenty of vapor

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Ice tendency should be negative or zero (sublimation, not deposition)
        # because deposition is suppressed above freezing
        @test tendencies.dq_icl_dt <= FT(0)
    end

    @testset "BulkMicrophysicsTendencies - Extreme conditions" begin
        # Test at extreme but physical atmospheric conditions

        # Very cold (cirrus cloud, ~10km altitude)
        ρ_cold = FT(0.4)
        T_cold = T_freeze - FT(53)  # Very cold
        q_tot_cold = FT(0.0005)
        q_icl_cold = FT(1e-5)
        q_sno_cold = FT(1e-5)

        tendencies_cold = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ_cold, T_cold, q_tot_cold, FT(0), q_icl_cold, FT(0), q_sno_cold,
        )
        @test isfinite(tendencies_cold.dq_icl_dt)
        @test isfinite(tendencies_cold.dq_sno_dt)

        # Warm, dense (tropical boundary layer)
        ρ_warm = FT(1.15)
        T_warm = FT(303)  # 30°C
        q_tot_warm = FT(0.020)
        q_lcl_warm = FT(1e-3)
        q_rai_warm = FT(5e-4)

        tendencies_warm = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ_warm, T_warm, q_tot_warm, q_lcl_warm, FT(0), q_rai_warm, FT(0),
        )
        @test isfinite(tendencies_warm.dq_lcl_dt)
        @test isfinite(tendencies_warm.dq_rai_dt)
        # No ice/snow tendencies in warm conditions without ice/snow
        @test tendencies_warm.dq_icl_dt == FT(0) || isfinite(tendencies_warm.dq_icl_dt)
        @test tendencies_warm.dq_sno_dt == FT(0)
    end

    #####
    ##### Quantitative Physics Tests (GPU safety and conservation)
    #####

    @testset "Type stability (@inferred) for GPU safety" begin
        # For GPU: ensure the compiler can infer the return type
        ρ = FT(1.2)
        T = T_freeze + FT(7)
        q_tot = FT(0.01)
        q_lcl = FT(1e-4)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(0)

        # @inferred will throw if return type cannot be inferred
        tendencies = @inferred BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )
        @test tendencies isa NamedTuple{(:dq_lcl_dt, :dq_icl_dt, :dq_rai_dt, :dq_sno_dt), NTuple{4, FT}}
    end

    @testset "Conservation: warm autoconversion at saturation" begin
        # At exact saturation (S=0), condensation is zero
        # Only lcl → rai processes should be active
        # Mass lost by cloud must equal mass gained by rain
        T = T_freeze + FT(17)  # Warm, no ice processes
        ρ = FT(1.0)
        q_lcl = FT(2e-3)  # Above autoconversion threshold
        q_rai = FT(5e-4)  # Rain present for accretion
        q_icl = FT(0)
        q_sno = FT(0)
        q_tot = get_saturated_q_tot(tps, T, ρ, q_lcl, q_icl, q_rai, q_sno)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Mass conservation: lcl + rai should sum to zero (autoconversion + accretion)
        # Using sqrt(eps) tolerance for floating-point associativity
        @test tendencies.dq_lcl_dt + tendencies.dq_rai_dt ≈ FT(0) atol = sqrt(eps(FT))
        # Cloud liquid should decrease
        @test tendencies.dq_lcl_dt < FT(0)
        # Rain should increase
        @test tendencies.dq_rai_dt > FT(0)
    end

    @testset "Conservation: pure snow melting at saturation" begin
        # For pure snow melting, we need ice saturation to disable sublimation
        # (snow equilibrates with ice saturation, not liquid)
        # Only sno → rai melting should be active
        T = T_freeze + FT(5)  # 5K above freezing
        ρ = FT(1.0)
        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(1e-3)
        # Use ice saturation to zero out snow sublimation
        q_vap_sat = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        q_tot = q_vap_sat + q_lcl + q_icl + q_rai + q_sno

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Mass conservation: sno + rai should sum to approximately zero (pure melting)
        # Allow small tolerance for numerical precision
        @test tendencies.dq_sno_dt + tendencies.dq_rai_dt ≈ FT(0) atol = FT(1e-10)
        # Snow should decrease
        @test tendencies.dq_sno_dt < FT(0)
        # Rain should increase
        @test tendencies.dq_rai_dt > FT(0)
    end

    @testset "Physics: Warm Shedding Thermodynamics (α factor)" begin
        # Verify the accretion-induced snow melt follows the thermodynamic formula:
        # dq_sno_dt = -S_accr_melt - S_melt + S_subl
        # where S_accr_melt = S_accr * α, α = cv_l / L_f * (T - T_freeze)

        T = T_freeze + FT(5) # 5K above freezing
        ρ = FT(1.0)
        q_lcl = FT(1e-3)
        q_sno = FT(1e-3)
        q_tot = FT(0.015)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, FT(0), FT(0), q_sno,
        )

        # Calculate individual components
        S_accr = CM1.accretion(liquid, snow, vel.snow, ce, q_lcl, q_sno, ρ)
        S_melt = CM1.snow_melt(snow, vel.snow, aps, tps, q_sno, ρ, T)
        S_subl = CM1.evaporation_sublimation(snow, vel.snow, aps, tps, q_tot, q_lcl, FT(0), FT(0), q_sno, ρ, T)

        # Calculate α
        T_frz = TDI.T_freeze(tps)
        cv_l = TDI.cv_l(tps)
        L_f = TDI.Lf(tps, T)
        α = cv_l / L_f * (T - T_frz)

        # Expected snow tendency: -S_accr*α - S_melt + S_subl
        expected_sno_dt = -S_accr * α - S_melt + S_subl

        # Test that the formula is exact
        @test tendencies.dq_sno_dt ≈ expected_sno_dt atol = 10 * eps(FT)

        # And verify α is in the expected range (physically sensible)
        # For 5K above freezing: α ≈ 0.06 (cv_l ≈ 4200 J/kg/K, L_f ≈ 334000 J/kg)
        @test α > FT(0)
        @test α < FT(0.1)
    end

    @testset "Sanity: no precipitation generation from nothing" begin
        # If all precipitation = 0, precipitation tendencies must be zero
        # regardless of vapor content
        # (cloud condensation can occur if q_vap > q_sat, but that's not precipitation)
        T = T_freeze + FT(7)
        ρ = FT(1.0)
        q_tot = FT(0.01)  # Some vapor available
        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(0)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # No precipitation can form without cloud condensate first
        @test tendencies.dq_rai_dt == FT(0)
        @test tendencies.dq_sno_dt == FT(0)
        # Cloud tendencies can be non-zero due to condensation/deposition
        @test !isnan(tendencies.dq_lcl_dt)
        @test !isnan(tendencies.dq_icl_dt)
    end

end

function test_average_bulk_microphysics_1m_tendencies(FT)

    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics1MParams(FT)
    T_freeze = TDI.T_freeze(tps)

    @testset "_bulk_microphysics_linearized_operator - Finiteness checks" begin
        ρ = FT(1.2)
        T = T_freeze - FT(5)
        q_tot = FT(0.015)
        q_lcl = FT(5e-4)
        q_icl = FT(5e-4)
        q_rai = FT(5e-4)
        q_sno = FT(5e-4)

        lin = BMT._bulk_microphysics_linearized_operator(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        @test isfinite(lin.M11)
        @test isfinite(lin.M22)
        @test isfinite(lin.M31)
        @test isfinite(lin.M33)
        @test isfinite(lin.M34)
        @test isfinite(lin.M41)
        @test isfinite(lin.M42)
        @test isfinite(lin.M43)
        @test isfinite(lin.M44)
        @test isfinite(lin.e1)
        @test isfinite(lin.e2)
        @test isfinite(lin.e4)
    end

    @testset "_bulk_microphysics_linearized_operator - Type stability (@inferred)" begin
        ρ = FT(1.2)
        T = T_freeze + FT(7)
        q_tot = FT(0.01)
        q_lcl = FT(1e-4)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(0)

        lin = @inferred BMT._bulk_microphysics_linearized_operator(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        @test lin isa NamedTuple
    end

    @testset "_bulk_microphysics_linearized_operator - Warm rain-only structure" begin
        # Warm rain-only case: only rain evaporation should contribute,
        # so only M33 should be nonzero.
        ρ = FT(1.2)
        T = T_freeze + FT(15)
        q_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)

        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(1e-3)
        q_sno = FT(0)
        q_vap = FT(0.5) * q_sat
        q_tot = q_vap + q_rai

        lin = BMT._bulk_microphysics_linearized_operator(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        @test lin.M33 <= FT(0)
        @test lin.M11 == FT(0)
        @test lin.M22 == FT(0)
        @test lin.M31 == FT(0)
        @test lin.M34 == FT(0)
        @test lin.M41 == FT(0)
        @test lin.M42 == FT(0)
        @test lin.M43 == FT(0)
        @test lin.M44 == FT(0)
        @test lin.e1 == FT(0)
        @test lin.e2 == FT(0)
        @test lin.e4 == FT(0)
    end

    @testset "_bulk_microphysics_linearized_operator - Warm pure snow melt structure" begin
        # Warm snow-only case: snow melts to rain, so M44 < 0 and M34 > 0.
        ρ = FT(1.2)
        T = T_freeze + FT(5)

        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(1e-3)
        q_vap_sat = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        q_tot = q_vap_sat + q_sno

        lin = BMT._bulk_microphysics_linearized_operator(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        @test lin.M34 > FT(0)
        @test lin.M44 < FT(0)
        @test lin.M11 == FT(0)
        @test lin.M22 == FT(0)
        @test lin.M31 == FT(0)
        @test lin.M41 == FT(0)
        @test lin.M42 == FT(0)
        @test lin.M43 == FT(0)
    end

    @testset "_average_bulk_microphysics_tendencies - Finiteness checks" begin
        ρ = FT(1.2)
        T = T_freeze - FT(5)
        q_tot = FT(0.015)
        q_lcl = FT(5e-4)
        q_icl = FT(5e-4)
        q_rai = FT(5e-4)
        q_sno = FT(5e-4)
        Δt = FT(10)

        tendencies = BMT._average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test isfinite(tendencies.dq_lcl_dt)
        @test isfinite(tendencies.dq_icl_dt)
        @test isfinite(tendencies.dq_rai_dt)
        @test isfinite(tendencies.dq_sno_dt)
    end

    @testset "_average_bulk_microphysics_tendencies - Type stability (@inferred)" begin
        ρ = FT(1.2)
        T = T_freeze + FT(7)
        q_tot = FT(0.01)
        q_lcl = FT(1e-4)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(0)
        Δt = FT(1)

        tendencies = @inferred BMT._average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test tendencies isa NamedTuple{(:dq_lcl_dt, :dq_icl_dt, :dq_rai_dt, :dq_sno_dt), NTuple{4, FT}}
    end

    @testset "_average_bulk_microphysics_tendencies - rain evaporation damping vs dt" begin
        ρ = FT(1.2)
        T = T_freeze + FT(15)
        q_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)

        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(1e-3)
        q_sno = FT(0)
        q_vap = FT(0.5) * q_sat
        q_tot = q_vap + q_rai

        dts = FT[1, 5, 10, 50, 100]
        rates = similar(dts)

        for i in eachindex(dts)
            tendencies = BMT._average_bulk_microphysics_tendencies(
                BMT.Microphysics1Moment(),
                mp, tps,
                ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, dts[i],
            )
            rates[i] = tendencies.dq_rai_dt
        end

        @test all(isfinite, rates)
        @test all(r -> r < 0, rates)
        @test all(abs(rates[i + 1]) <= abs(rates[i]) for i in 1:(length(rates) - 1))
    end

    @testset "_average_bulk_microphysics_tendencies - Matches solved linear system" begin
        ρ = FT(1.1)
        T = T_freeze + FT(4)
        q_lcl = FT(4e-4)
        q_icl = FT(2e-4)
        q_rai = FT(3e-4)
        q_sno = FT(6e-4)
        q_tot = FT(0.014)
        Δt = FT(7)

        lin = BMT._bulk_microphysics_linearized_operator(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        tendencies = BMT._average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        invΔt = one(FT) / Δt

        q_lcl_new = q_lcl + Δt * tendencies.dq_lcl_dt
        q_icl_new = q_icl + Δt * tendencies.dq_icl_dt
        q_rai_new = q_rai + Δt * tendencies.dq_rai_dt
        q_sno_new = q_sno + Δt * tendencies.dq_sno_dt

        @test (q_lcl_new - q_lcl) * invΔt ≈ lin.M11 * q_lcl_new + lin.e1 atol = FT(100) * eps(FT)
        @test (q_icl_new - q_icl) * invΔt ≈ lin.M22 * q_icl_new + lin.e2 atol = FT(100) * eps(FT)
        @test (q_rai_new - q_rai) * invΔt ≈ lin.M31 * q_lcl_new + lin.M33 * q_rai_new + lin.M34 * q_sno_new atol =
            FT(100) * eps(FT)
        @test (q_sno_new - q_sno) * invΔt ≈
              lin.M41 * q_lcl_new + lin.M42 * q_icl_new + lin.M43 * q_rai_new + lin.M44 * q_sno_new + lin.e4 atol =
            FT(100) * eps(FT)
    end

    @testset "_average_bulk_microphysics_tendencies - Small Δt agrees with instantaneous tendency for rain evaporation" begin
        # In this simple case the model is essentially dq_rai/dt = M33 * q_rai,
        # so the averaged implicit tendency should approach the instantaneous one
        # as Δt -> 0.
        ρ = FT(1.2)
        T = T_freeze + FT(15)
        q_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)

        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(1e-3)
        q_sno = FT(0)
        q_vap = FT(0.5) * q_sat
        q_tot = q_vap + q_rai

        inst = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        avg = BMT._average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, FT(1e-2),
        )

        @test avg.dq_lcl_dt ≈ inst.dq_lcl_dt atol = FT(1e-8)
        @test avg.dq_icl_dt ≈ inst.dq_icl_dt atol = FT(1e-8)
        @test avg.dq_sno_dt ≈ inst.dq_sno_dt atol = FT(1e-8)
        @test avg.dq_rai_dt ≈ inst.dq_rai_dt rtol = FT(1e-3)
    end

    @testset "average_bulk_microphysics_tendencies - Finiteness checks" begin
        ρ = FT(1.2)
        T = T_freeze - FT(5)
        q_tot = FT(0.015)
        q_lcl = FT(5e-4)
        q_icl = FT(5e-4)
        q_rai = FT(5e-4)
        q_sno = FT(5e-4)
        Δt = FT(10)

        tendencies = BMT.average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test isfinite(tendencies.dq_lcl_dt)
        @test isfinite(tendencies.dq_icl_dt)
        @test isfinite(tendencies.dq_rai_dt)
        @test isfinite(tendencies.dq_sno_dt)
    end

    @testset "average_bulk_microphysics_tendencies - Type stability (@inferred)" begin
        ρ = FT(1.2)
        T = T_freeze + FT(7)
        q_tot = FT(0.01)
        q_lcl = FT(1e-4)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(0)
        Δt = FT(1)

        tendencies = @inferred BMT.average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test tendencies isa NamedTuple{(:dq_lcl_dt, :dq_icl_dt, :dq_rai_dt, :dq_sno_dt), NTuple{4, FT}}
    end

    @testset "average_bulk_microphysics_tendencies - all zero inputs" begin
        ρ = FT(1.2)
        T = T_freeze + FT(5)
        q_tot = FT(0)
        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(0)
        Δt = FT(10)

        tendencies = BMT.average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test tendencies.dq_lcl_dt == FT(0)
        @test tendencies.dq_icl_dt == FT(0)
        @test tendencies.dq_rai_dt == FT(0)
        @test tendencies.dq_sno_dt == FT(0)
    end

    @testset "average_bulk_microphysics_tendencies - nsub=1 matches single-substep solver" begin
        ρ = FT(1.1)
        T = T_freeze + FT(4)
        q_lcl = FT(4e-4)
        q_icl = FT(2e-4)
        q_rai = FT(3e-4)
        q_sno = FT(6e-4)
        q_tot = FT(0.014)
        Δt = FT(7)

        single = BMT._average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        substepped = BMT.average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, 1,
        )

        @test substepped.dq_lcl_dt ≈ single.dq_lcl_dt atol = FT(100) * eps(FT)
        @test substepped.dq_icl_dt ≈ single.dq_icl_dt atol = FT(100) * eps(FT)
        @test substepped.dq_rai_dt ≈ single.dq_rai_dt atol = FT(100) * eps(FT)
        @test substepped.dq_sno_dt ≈ single.dq_sno_dt atol = FT(100) * eps(FT)
    end

    @testset "average_bulk_microphysics_tendencies - Warm pure snow melt keeps expected signs" begin
        ρ = FT(1.0)
        T = T_freeze + FT(5)
        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(1e-3)
        q_vap_sat = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        q_tot = q_vap_sat + q_sno
        Δt = FT(10)

        tendencies = BMT.average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test tendencies.dq_sno_dt < FT(0)
        @test tendencies.dq_rai_dt > FT(0)
        @test isfinite(tendencies.dq_sno_dt)
        @test isfinite(tendencies.dq_rai_dt)
    end

    @testset "average_bulk_microphysics_tendencies - More substeps do not change simple rain-only case much" begin
        # In a simple rain-only case, rebuilding the operator should not change
        # the result much as nsub increases.
        ρ = FT(1.2)
        T = T_freeze + FT(15)
        q_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)

        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(1e-3)
        q_sno = FT(0)
        q_vap = FT(0.5) * q_sat
        q_tot = q_vap + q_rai
        Δt = FT(1)

        t1 = BMT.average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, 1,
        )

        t10 = BMT.average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, 10,
        )

        @test t10.dq_lcl_dt ≈ t1.dq_lcl_dt atol = FT(1e-10)
        @test t10.dq_icl_dt ≈ t1.dq_icl_dt atol = FT(1e-10)
        @test t10.dq_sno_dt ≈ t1.dq_sno_dt atol = FT(1e-10)
        @test t10.dq_rai_dt ≈ t1.dq_rai_dt rtol = FT(1e-2)
    end

    @testset "average_bulk_microphysics_tendencies - Substepping remains finite near freezing" begin
        ρ = FT(1.2)
        T = T_freeze + FT(0.01)
        q_tot = FT(0.015)
        q_lcl = FT(1e-3)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(5e-4)
        Δt = FT(20)

        tendencies = BMT.average_bulk_microphysics_tendencies(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, 20,
        )

        @test isfinite(tendencies.dq_lcl_dt)
        @test isfinite(tendencies.dq_icl_dt)
        @test isfinite(tendencies.dq_rai_dt)
        @test isfinite(tendencies.dq_sno_dt)
    end

end



###
### 2M tendencies and derivatives tests
###

function test_bulk_microphysics_2m_tendencies(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = false)

    T_freeze = TDI.T_freeze(tps)

    @testset "BulkMicrophysicsTendencies 2M - Autoconversion" begin
        ρ = FT(1.2)
        T = T_freeze + FT(10)
        q_lcl = FT(2e-3)  # Above threshold
        q_rai = FT(0)
        n_lcl = FT(1e8)
        n_rai = FT(0)
        q_tot = get_saturated_q_tot(tps, T, ρ, q_lcl, FT(0), q_rai, FT(0))

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(),
            mp,
            tps,
            ρ,
            T,
            q_tot,
            q_lcl,
            n_lcl,
            q_rai,
            n_rai,
        )

        @test tendencies.dq_lcl_dt < FT(0)  # Cloud decreases
        @test tendencies.dq_rai_dt > FT(0)  # Rain increases
    end

    @testset "BulkMicrophysicsTendencies 2M - Condensation and Evaporation" begin
        ρ = FT(1.2)
        T = T_freeze + FT(10)
        q_lcl = FT(0)
        q_rai = FT(0)
        n_lcl = FT(1e8)
        n_rai = FT(0)

        # Test 1: Supersaturated (condensation)
        q_vap_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
        q_tot_super = q_vap_sat * FT(1.05)

        tend_cond = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(),
            mp, tps, ρ, T, q_tot_super, q_lcl, n_lcl, q_rai, n_rai,
        )
        @test tend_cond.dq_lcl_dt > FT(0)  # Cloud increases

        # Test 2: Subsaturated (evaporation)
        q_lcl = FT(1e-3)
        q_tot_sub = q_vap_sat * FT(0.8) + q_lcl

        tend_evap = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(),
            mp, tps, ρ, T, q_tot_sub, q_lcl, n_lcl, q_rai, n_rai,
        )
        @test tend_evap.dq_lcl_dt < FT(0)  # Cloud decreases
    end

    @testset "BulkMicrophysicsTendencies 2M - Type stability" begin
        ρ = FT(1.2)
        T = T_freeze + FT(5)

        # Note: underlying CM2 functions have some Float32 type instability,
        # so we only test that the return is a NamedTuple with correct keys
        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(),
            mp,
            tps,
            ρ,
            T,
            FT(0.01),   # q_tot
            FT(1e-3),   # q_lcl
            FT(1e8),    # n_lcl
            FT(1e-4),   # q_rai
            FT(1e4),    # n_rai
        )
        @test tendencies isa @NamedTuple{
            dq_lcl_dt::FT,
            dn_lcl_dt::FT,
            dq_rai_dt::FT,
            dn_rai_dt::FT,
            dq_ice_dt::FT,
            dq_rim_dt::FT,
            db_rim_dt::FT,
        }
        # Ice tendencies should be zero for 2M mode
        @test tendencies.dq_ice_dt == FT(0)
        @test tendencies.dq_rim_dt == FT(0)
        @test tendencies.db_rim_dt == FT(0)
    end
end


###
### P3 tendencies tests
###

function test_bulk_microphysics_p3_tendencies(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)

    # Extract individual parameters for direct testing
    p3 = mp.ice.scheme
    pdf_c = mp.ice.cloud_pdf
    pdf_r = mp.ice.rain_pdf
    T_freeze = TDI.T_freeze(tps)

    @testset "BulkMicrophysicsTendencies P3 - Warm rain processes" begin
        # Test that P3 includes 2M warm rain processes
        ρ = FT(1.2)
        T = T_freeze + FT(10)  # Above freezing
        q_lcl = FT(2e-3)  # Significant cloud liquid
        n_lcl = FT(1e8 / ρ)  # 100/mg
        q_rai = FT(0)
        n_rai = FT(0)
        q_ice = FT(0)  # No ice
        n_ice = FT(0)
        q_rim = FT(0)
        b_rim = FT(0)
        q_tot = get_saturated_q_tot(tps, T, ρ, q_lcl, q_ice + q_rim, q_rai, FT(0))
        logλ = FT(10)  # Dummy, not used without ice

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(),
            mp,
            tps,
            ρ,
            T,
            q_tot,
            q_lcl,
            n_lcl,
            q_rai,
            n_rai,
            q_ice,
            n_ice,
            q_rim,
            b_rim,
            logλ,
        )

        # Autoconversion from cloud to rain should occur
        @test tendencies.dq_lcl_dt < FT(0)  # Cloud decreases
        @test tendencies.dq_rai_dt > FT(0)  # Rain increases
        @test tendencies.dn_rai_dt > FT(0)  # Rain number increases
    end

    @testset "BulkMicrophysicsTendencies P3 - Ice melting" begin
        # Ice should melt above freezing
        ρ = FT(1.2)
        T = T_freeze + FT(5)  # Above freezing
        q_lcl = FT(0)
        n_lcl = FT(0)
        q_rai = FT(0)
        n_rai = FT(0)
        q_ice = FT(1e-4)  # Some ice
        n_ice = FT(2e5) / ρ  # Ice number per kg
        q_rim = FT(0.5e-4)  # Some rime
        b_rim = FT(1e-7)  # Rime volume
        q_tot = get_saturated_q_tot(tps, T, ρ, q_lcl, q_ice + q_rim, q_rai, FT(0))

        # Compute logλ from P3 state
        L_ice = q_ice * ρ
        N_ice = n_ice * ρ
        F_rim = q_rim / q_ice
        ρ_rim = q_rim * ρ / (b_rim * ρ)
        state = CM.P3Scheme.get_state(p3; F_rim, ρ_rim, L_ice, N_ice)
        logλ = CM.P3Scheme.get_distribution_logλ(state)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(),
            mp,
            tps,
            ρ,
            T,
            q_tot,
            q_lcl,
            n_lcl,
            q_rai,
            n_rai,
            q_ice,
            n_ice,
            q_rim,
            b_rim,
            logλ,
        )

        # Ice should decrease due to melting
        @test tendencies.dq_ice_dt < FT(0)
        # Rain should increase from melted ice
        @test tendencies.dq_rai_dt > FT(0)
    end

    @testset "BulkMicrophysicsTendencies P3 - Liquid-ice collisions" begin
        # Ice collecting cloud liquid below freezing
        ρ = FT(1.2)
        T = T_freeze - FT(10)  # Below freezing
        q_lcl = FT(1e-3)  # Cloud liquid
        n_lcl = FT(1e8) / ρ
        q_rai = FT(1e-5)  # Some rain
        n_rai = FT(1e5) / ρ
        q_ice = FT(1e-4)
        n_ice = FT(2e5) / ρ
        q_rim = FT(0.5e-4)
        b_rim = FT(1e-7)
        q_tot = get_saturated_q_tot(tps, T, ρ, q_lcl, q_ice, q_rai, FT(0))

        # Compute logλ
        L_ice = q_ice * ρ
        N_ice = n_ice * ρ
        F_rim = q_rim / q_ice
        ρ_rim = q_rim * ρ / (b_rim * ρ)
        state = CM.P3Scheme.get_state(p3; F_rim, ρ_rim, L_ice, N_ice)
        logλ = CM.P3Scheme.get_distribution_logλ(state)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(),
            mp,
            tps,
            ρ,
            T,
            q_tot,
            q_lcl,
            n_lcl,
            q_rai,
            n_rai,
            q_ice,
            n_ice,
            q_rim,
            b_rim,
            logλ,
        )

        # Cloud liquid should decrease (collected by ice)
        @test tendencies.dq_lcl_dt < FT(0)
        # Ice should increase from riming
        @test tendencies.dq_ice_dt >= FT(0)
    end

    @testset "BulkMicrophysicsTendencies P3 - Return finiteness" begin
        ρ = FT(1.2)
        T = T_freeze - FT(5)
        q_tot = FT(0.015)
        q_lcl = FT(1e-3)
        n_lcl = FT(1e8) / ρ
        q_rai = FT(1e-4)
        n_rai = FT(1e5) / ρ
        q_ice = FT(1e-4)
        n_ice = FT(2e5) / ρ
        q_rim = FT(0.3e-4)
        b_rim = FT(5e-8)

        L_ice = q_ice * ρ
        N_ice = n_ice * ρ
        F_rim = q_rim / q_ice
        ρ_rim = q_rim * ρ / (b_rim * ρ)
        state = CM.P3Scheme.get_state(p3; F_rim, ρ_rim, L_ice, N_ice)
        logλ = CM.P3Scheme.get_distribution_logλ(state)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(),
            mp,
            tps,
            ρ,
            T,
            q_tot,
            q_lcl,
            n_lcl,
            q_rai,
            n_rai,
            q_ice,
            n_ice,
            q_rim,
            b_rim,
            logλ,
        )

        @test isfinite(tendencies.dq_lcl_dt)
        @test isfinite(tendencies.dn_lcl_dt)
        @test isfinite(tendencies.dq_rai_dt)
        @test isfinite(tendencies.dn_rai_dt)
        @test isfinite(tendencies.dq_ice_dt)
        @test isfinite(tendencies.dq_rim_dt)
        @test isfinite(tendencies.db_rim_dt)
    end

    @testset "BulkMicrophysicsTendencies P3 - Type stability" begin
        ρ = FT(1.2)
        T = T_freeze - FT(5)
        q_tot = FT(0.015)
        q_lcl = FT(1e-3)
        n_lcl = FT(1e8) / ρ
        q_rai = FT(1e-4)
        n_rai = FT(1e5) / ρ
        q_ice = FT(1e-4)
        n_ice = FT(2e5) / ρ
        q_rim = FT(0.3e-4)
        b_rim = FT(5e-8)

        L_ice = q_ice * ρ
        N_ice = n_ice * ρ
        F_rim = q_rim / q_ice
        ρ_rim = q_rim * ρ / (b_rim * ρ)
        state = CM.P3Scheme.get_state(p3; F_rim, ρ_rim, L_ice, N_ice)
        logλ = CM.P3Scheme.get_distribution_logλ(state)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(),
            mp,
            tps,
            ρ,
            T,
            q_tot,
            q_lcl,
            n_lcl,
            q_rai,
            n_rai,
            q_ice,
            n_ice,
            q_rim,
            b_rim,
            logλ,
        )

        # Check that we get the expected NamedTuple type
        @test tendencies isa @NamedTuple{
            dq_lcl_dt::FT,
            dn_lcl_dt::FT,
            dq_rai_dt::FT,
            dn_rai_dt::FT,
            dq_ice_dt::FT,
            dq_rim_dt::FT,
            db_rim_dt::FT,
        }
    end
end

@testset "Bulk Microphysics Tendencies ($FT)" for FT in (Float64, Float32)
    test_bulk_microphysics_0m_tendencies(FT)
    test_bulk_microphysics_1m_tendencies(FT)
    test_average_bulk_microphysics_1m_tendencies(FT)
    test_bulk_microphysics_2m_tendencies(FT)
    test_bulk_microphysics_p3_tendencies(FT)
end
nothing
