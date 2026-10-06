using Test

import ClimaParams as CP

import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP

"""
    stiff_prescribed_ice_params(FT)

`Microphysics1MParams` with `PrescribedIceNumber` and a dense prescribed ice number (5e8 m⁻³),
the value the stiff vapor <-> ice relaxation tests below were written for. Pinned here so the
tests do not depend on the ClimaParams default (1e5 m⁻³ since ClimaParams 1.1.12).
"""
function stiff_prescribed_ice_params(FT; max_latent_heating_rate = Inf, N_0 = 5.0e8)
    override = Dict(
        "cloud_ice_sedimentation_number_concentration" => Dict("value" => N_0, "type" => "float"),
        # the solver tests probe the relaxation itself; the latent-heating limiter has its own tests
        "microphysics_max_latent_heating_rate" => Dict("value" => max_latent_heating_rate, "type" => "float"),
    )
    td = CP.create_toml_dict(FT; override_file = override)
    return CMP.Microphysics1MParams(td; cloud_ice_formation = CMP.PrescribedIceNumber())
end

"""
1-moment parameters with the latent-heating limiter disabled, for the solver tests that compare
against unlimited references; the limiter has its own tests.
"""
function params_1m_no_limiter(FT)
    override = Dict("microphysics_max_latent_heating_rate" => Dict("value" => Inf, "type" => "float"))
    return CMP.Microphysics1MParams(CP.create_toml_dict(FT; override_file = override))
end
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
    E_lcl_sno = mp.process_params.cloud_liquid_snow_accretion.e
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Compute autoconversion separately using the same option-dispatched function
        micro_acnv = (; q_tot, q_lcl, q_icl, q_rai, q_sno)
        thermo_acnv = (; ρ, T, w = FT(0))
        S_acnv_lcl = CM1.conv_q_lcl_to_q_rai(mp.processes.rain_autoconversion, mp, tps, micro_acnv, thermo_acnv)

        # Rain tendency should be positive and include autoconversion
        @test tendencies.dq_rai_dt > FT(0)
        # Autoconversion should be a major component of rain formation
        # (in the absence of rain, there's no accretion, so it should be close)
        @test abs(tendencies.dq_rai_dt - S_acnv_lcl) / S_acnv_lcl < FT(0.1)
    end

    @testset "BulkMicrophysicsTendencies - Autoconversion component (PrescribedNd)" begin
        # Same check as above but using the variable-timescale 2M autoconversion option
        # (PrescribedNd), which uses the prescribed {τ,α,Nc}.
        mp_2m = CMP.Microphysics1MParams(FT;
            rain_autoconversion = CMP.PrescribedNd(),
        )

        ρ = FT(1.2)
        T = T_freeze + FT(10)  # warm, no ice processes
        q_tot = FT(0.015)
        q_lcl = FT(2e-3)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(0)

        tendencies_2m = BMT.bulk_microphysics_tendencies(
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp_2m, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Compute autoconversion separately using the option-dispatched function
        micro_acnv = (; q_tot, q_lcl, q_icl, q_rai, q_sno)
        thermo_acnv = (; ρ, T, w = FT(0))
        S_acnv_2m = CM1.conv_q_lcl_to_q_rai(mp_2m.processes.rain_autoconversion, mp_2m, tps, micro_acnv, thermo_acnv)

        # Rain tendency should be positive and dominated by autoconversion
        @test tendencies_2m.dq_rai_dt > FT(0)
        @test abs(tendencies_2m.dq_rai_dt - S_acnv_2m) / S_acnv_2m < FT(0.1)
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Verify accretion is happening by checking snow growth
        # (cloud liquid tendency may be positive overall due to condensation)
        S_accr = CM1.accretion(liquid, snow, vel.snow, E_lcl_sno, q_lcl, q_sno, ρ)
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Verify accretion is happening
        S_accr = CM1.accretion(liquid, snow, vel.snow, E_lcl_sno, q_lcl, q_sno, ρ)
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # Snow should increase from deposition (conv_q_sno_to_q_vap returns positive for S > 0)
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ_cold, T_cold, FT(0), q_tot_cold, FT(0), q_icl_cold, FT(0), q_sno_cold,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ_warm, T_warm, FT(0), q_tot_warm, q_lcl_warm, FT(0), q_rai_warm, FT(0),
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, FT(0), FT(0), q_sno,
        )

        # Calculate individual components
        S_accr = CM1.accretion(liquid, snow, vel.snow, E_lcl_sno, q_lcl, q_sno, ρ)
        micro_s = (; q_tot, q_lcl, q_icl = FT(0), q_rai = FT(0), q_sno)
        thermo_s = (; ρ, T)
        S_melt = CM1.conv_q_sno_to_q_rai(CMP.SnowMelt(), mp, tps, micro_s, thermo_s)
        S_subl = CM1.conv_q_sno_to_q_vap(CMP.DepositionAndSublimation(), mp, tps, micro_s, thermo_s)

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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        # No precipitation can form without cloud condensate first
        @test tendencies.dq_rai_dt == FT(0)
        @test tendencies.dq_sno_dt == FT(0)
        # Cloud tendencies can be non-zero due to condensation/deposition
        @test !isnan(tendencies.dq_lcl_dt)
        @test !isnan(tendencies.dq_icl_dt)
    end

end

function test_linearized_bulk_microphysics_1m_tendencies(FT)

    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = params_1m_no_limiter(FT)
    T_freeze = TDI.T_freeze(tps)

    @testset "LinearizedAverage - stale call shape without w throws" begin
        mp_stale = CMP.Microphysics1MParams(FT)
        ρ = FT(1.2)
        T = FT(280)
        q_tot = FT(1e-3)
        q_lcl = FT(2e-4)
        Δt = FT(60)
        # pre-#770 argument order (no w): Δt would receive an Integer -> must not run with shifted arguments
        @test_throws MethodError BMT.bulk_microphysics_tendencies(
            BMT.LinearizedAverage(), BMT.Microphysics1Moment(), mp_stale, tps,
            ρ, T, q_tot, q_lcl, FT(0), FT(0), FT(0), Δt, 3,
        )
    end

    @testset "_linearize (via _microphysics_source_terms) - Finiteness checks" begin
        ρ = FT(1.2)
        T = T_freeze - FT(5)
        q_tot = FT(0.015)
        q_lcl = FT(5e-4)
        q_icl = FT(5e-4)
        q_rai = FT(5e-4)
        q_sno = FT(5e-4)

        q_min = TDI.TD.Parameters.q_min(tps)

        src = BMT._microphysics_source_terms(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        lin = BMT._linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min, FT(60))

        @test isfinite(lin.M11)
        @test isfinite(lin.M12)
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

    @testset "_linearize (via _microphysics_source_terms) - Type stability (@inferred)" begin
        ρ = FT(1.2)
        T = T_freeze + FT(7)
        q_tot = FT(0.01)
        q_lcl = FT(1e-4)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(0)

        q_min = TDI.TD.Parameters.q_min(tps)

        src = @inferred BMT._microphysics_source_terms(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        lin = BMT._linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min, FT(60))

        @test lin isa NamedTuple
    end

    @testset "_linearize (via _microphysics_source_terms) - Warm rain-only structure" begin
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

        q_min = TDI.TD.Parameters.q_min(tps)

        src = BMT._microphysics_source_terms(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        lin = BMT._linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min, FT(60))

        @test lin.M33 <= FT(0)
        @test lin.M11 == FT(0)
        @test lin.M12 == FT(0)
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

    @testset "_linearize (via _microphysics_source_terms) - Warm pure snow melt structure" begin
        # Warm snow-only case: snow melts to rain, so M44 < 0 and M34 > 0.
        ρ = FT(1.2)
        T = T_freeze + FT(5)

        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(1e-3)
        q_vap_sat = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        q_tot = q_vap_sat + q_sno

        q_min = TDI.TD.Parameters.q_min(tps)

        src = BMT._microphysics_source_terms(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        lin = BMT._linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min, FT(60))

        @test lin.M34 > FT(0)
        @test lin.M44 < FT(0)
        @test lin.M11 == FT(0)
        @test lin.M12 == FT(0)
        @test lin.M22 == FT(0)
        @test lin.M31 == FT(0)
        @test lin.M41 == FT(0)
        @test lin.M42 == FT(0)
        @test lin.M43 == FT(0)
    end

    @testset "_relaxation_transfer limits" begin
        S, Δt = FT(2e-6), FT(60)
        # disabled process: τ = Inf with S = 0 gives exactly zero, not 0 ⋅ Inf
        @test BMT._relaxation_transfer(FT(0), FT(Inf), Δt) == FT(0)
        @test BMT._relaxation_transfer(S, FT(Inf), Δt) ≈ S * Δt rtol = FT(1e-5)
        # Δt ≪ τ: the instantaneous rate is recovered
        @test BMT._relaxation_transfer(S, FT(1e6), Δt) ≈ S * Δt rtol = FT(1e-4)
        # intermediate: exact exponential relaxation
        @test BMT._relaxation_transfer(S, FT(300), Δt) ≈ S * FT(300) * (1 - exp(FT(-0.2)))
        # Δt ≫ τ: saturates at the equilibrium amount S τ
        @test BMT._relaxation_transfer(S, FT(1), Δt) ≈ S * FT(1) rtol = FT(1e-6)
        # odd in S; type-stable
        @test BMT._relaxation_transfer(-S, FT(300), Δt) == -BMT._relaxation_transfer(S, FT(300), Δt)
        @test BMT._relaxation_transfer(S, FT(Inf), Δt) isa FT
    end

    @testset "LinearizedAverage - stiff sublimation with competing sinks: closes the deficit, q ≥ 0, no overshoot" begin
        # A 10 % subsaturated ice cloud, ice above the autoconversion threshold, snow present:
        # with N_0 = 5e8 m⁻³ the sublimation timescale is ~2 s against a 60 s substep. The
        # implicit sink must remove about the vapor deficit (not all the ice, as the plain
        # S/q decay does, which then overshoots to ~13 % supersaturation) while the other
        # sinks act, and the ice must stay non-negative.
        mp_presc = stiff_prescribed_ice_params(FT)
        Ls = TDI.TD.Parameters.LH_s0(tps)
        cp = TDI.TD.Parameters.cp_d(tps)
        Δt = FT(60)
        ρ = FT(0.7)
        T = FT(255)
        q_icl = FT(3e-4)
        q_sno = FT(2e-4)
        q_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        q_tot = FT(0.9) * q_sat_ice + q_icl + q_sno
        r = BMT.bulk_microphysics_tendencies(
            BMT.LinearizedAverage(), BMT.Microphysics1Moment(),
            mp_presc, tps, ρ, T, FT(0), q_tot, FT(0), q_icl, FT(0), q_sno, Δt, 1,
        )
        q_icl_new = q_icl + r.dq_icl_dt * Δt
        q_sno_new = q_sno + r.dq_sno_dt * Δt
        T_new = T + Ls / cp * (r.dq_icl_dt + r.dq_sno_dt) * Δt
        q_v_new = q_tot - q_icl_new - q_sno_new - (r.dq_lcl_dt + r.dq_rai_dt) * Δt
        RHi_new = q_v_new / TDI.saturation_vapor_specific_content_over_ice(tps, T_new, ρ)
        @test all(isfinite, (q_icl_new, q_sno_new, T_new))
        @test q_icl_new >= FT(0)
        @test q_icl_new < q_icl
        @test RHi_new >= FT(0.9)
        # cloud-ice sublimation closes the deficit exactly; the small excess over
        # saturation comes from snow sublimation (constant-rate, pre-existing treatment).
        # A plain S/q decay of the cloud ice gives RHi_new ≈ 1.13 here.
        @test RHi_new <= FT(1.05)
        # deficit larger than the condensate: coefficients finite, ice non-negative
        q_tot2 = FT(0.4) * q_sat_ice + FT(5e-5)
        r2 = BMT.bulk_microphysics_tendencies(
            BMT.LinearizedAverage(), BMT.Microphysics1Moment(),
            mp_presc, tps, ρ, T, FT(0), q_tot2, FT(0), FT(5e-5), FT(0), FT(0), Δt, 1,
        )
        @test isfinite(r2.dq_icl_dt)
        @test FT(5e-5) + r2.dq_icl_dt * Δt >= FT(0)
    end

    @testset "LinearizedAverage - TemperatureDependentIceNumber: slow deposition in mixed-phase conditions" begin
        # With N_ice(T) ~ 8e2 m⁻³ at -10 °C the deposition timescale is hours, so a
        # liquid-saturated supercooled cloud loses only a small fraction of its liquid per
        # step; with the constant N_0 = 5e8 m⁻³ the vapor excess over ice saturation is
        # deposited within the step and the liquid evaporates toward it (WBF glaciation).
        # The net cloud-ice tendency is not tested: with deposition this slow, the ice sinks
        # (autoconversion to snow) dominate it.
        opt = CMP.TemperatureDependentIceNumber()
        mp_fit = CMP.Microphysics1MParams(FT; cloud_ice_formation = opt)
        mp_presc = stiff_prescribed_ice_params(FT)
        @test mp_fit.process_params.cloud_ice_formation isa CMP.IceNumberTemperatureFit
        Δt = FT(180)
        ρ = FT(0.8)
        T = FT(263)
        q_lcl = FT(2e-4)
        q_icl = FT(1e-5)
        q_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
        q_tot = q_sat_liq + q_lcl + q_icl   # liquid-saturated, hence supersaturated over ice
        micro = (; q_tot, q_lcl, q_icl, q_rai = FT(0), q_sno = FT(0))
        thermo = (; ρ, T)
        τ = CMNonEq.τ_vap_to_q_icl(opt, mp_fit, tps, micro, thermo)
        @test FT(3600) < τ < FT(48 * 3600)
        # instantaneous deposition rates
        S_dep = CMNonEq.conv_q_vap_to_q_icl(opt, mp_fit, tps, micro, thermo)
        S_dep_presc = CMNonEq.conv_q_vap_to_q_icl(CMP.PrescribedIceNumber(), mp_presc, tps, micro, thermo)
        @test S_dep > FT(0)
        @test S_dep * Δt < FT(0.05) * q_lcl      # a few percent of the liquid mass per step at most
        @test S_dep_presc > FT(100) * S_dep      # the constant N_0 deposits orders of magnitude faster
        # coupled substep solver: finite tendencies, liquid survives the step
        r = BMT.bulk_microphysics_tendencies(
            BMT.LinearizedAverage(), BMT.Microphysics1Moment(),
            mp_fit, tps, ρ, T, FT(0), q_tot, q_lcl, q_icl, FT(0), FT(0), Δt, 3,
        )
        @test all(isfinite, (r.dq_lcl_dt, r.dq_icl_dt, r.dq_rai_dt, r.dq_sno_dt))
        q_lcl_new = q_lcl + r.dq_lcl_dt * Δt
        @test q_lcl_new > FT(0.9) * q_lcl
        # the same state with the default N_0 = 5e8 m⁻³ loses much more liquid in one step
        r_presc = BMT.bulk_microphysics_tendencies(
            BMT.LinearizedAverage(), BMT.Microphysics1Moment(),
            mp_presc, tps, ρ, T, FT(0), q_tot, q_lcl, q_icl, FT(0), FT(0), Δt, 3,
        )
        q_lcl_new_presc = q_lcl + r_presc.dq_lcl_dt * Δt
        @test q_lcl_new_presc < q_lcl_new
        @test q_lcl - q_lcl_new_presc > FT(10) * max(q_lcl - q_lcl_new, FT(0))
    end

    @testset "LinearizedAverage - stiff vapor↔ice relaxation (PrescribedIceNumber) does not cross saturation" begin
        # With the default N_0 = 5e8 m⁻³ the deposition/sublimation timescale is 1-7 s,
        # far below a 60 s substep. Using the time-averaged relaxation over the substep must
        # bring the parcel toward ice saturation without crossing it; treating the
        # instantaneous rate as a constant source produced a period-2
        # deposition/sublimation flip-flop with ±1 K temperature swings in supercooled clouds.
        mp_presc = stiff_prescribed_ice_params(FT)
        Lv = TDI.TD.Parameters.LH_v0(tps)
        Ls = TDI.TD.Parameters.LH_s0(tps)
        cp = TDI.TD.Parameters.cp_d(tps)
        Δt = FT(60)
        ρ = FT(0.79)
        cases = (
            (FT(265), FT(1.10), FT(0), FT(0)),        # deposition leg, no liquid
            (FT(265), FT(0.90), FT(0), FT(3e-4)),     # sublimation leg, plenty of ice
            (FT(250), FT(1.15), FT(2e-4), FT(1e-6)),  # mixed phase (liquid evaporates while ice grows)
        )
        for (T, RHi, q_lcl, q_icl) in cases
            q_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
            q_tot = RHi * q_sat_ice + q_lcl + q_icl
            τ = CMNonEq.τ_relax(mp_presc.cloud.ice, mp_presc.air_properties, q_icl, ρ)
            @test Δt / τ > FT(5)  # the regime under test is genuinely stiff
            r = BMT.bulk_microphysics_tendencies(
                BMT.LinearizedAverage(), BMT.Microphysics1Moment(),
                mp_presc, tps, ρ, T, FT(0), q_tot, q_lcl, q_icl, FT(0), FT(0), Δt, 1,
            )
            q_lcl_new = q_lcl + r.dq_lcl_dt * Δt
            q_icl_new = q_icl + r.dq_icl_dt * Δt
            q_rai_new = r.dq_rai_dt * Δt
            q_sno_new = r.dq_sno_dt * Δt
            T_new = T + (Lv / cp * (r.dq_lcl_dt + r.dq_rai_dt) + Ls / cp * (r.dq_icl_dt + r.dq_sno_dt)) * Δt
            q_v_new = q_tot - q_lcl_new - q_icl_new - q_rai_new - q_sno_new
            RHi_new = q_v_new / TDI.saturation_vapor_specific_content_over_ice(tps, T_new, ρ)
            @test all(isfinite, (q_lcl_new, q_icl_new, q_rai_new, q_sno_new, T_new))
            @test q_icl_new >= FT(0)
            @test q_lcl_new >= FT(0)
            # The parcel must approach ice saturation from its initial side and not
            # cross it (2% tolerance for the first-order Γ linearization of the
            # latent-heat feedback). Other processes (e.g. ice → snow autoconversion
            # in the sublimation leg) may legitimately stop it short of saturation.
            if RHi > 1
                @test RHi_new >= FT(0.98)
                q_lcl == 0 && @test RHi_new <= RHi
            else
                @test RHi_new <= FT(1.02)
                @test RHi_new >= RHi
            end
        end
    end

    @testset "_linearized_implicit_step - Finiteness checks" begin
        ρ = FT(1.2)
        T = T_freeze - FT(5)
        q_tot = FT(0.015)
        q_lcl = FT(5e-4)
        q_icl = FT(5e-4)
        q_rai = FT(5e-4)
        q_sno = FT(5e-4)
        Δt = FT(10)

        tendencies = BMT._linearized_implicit_step(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test isfinite(tendencies.dq_lcl_dt)
        @test isfinite(tendencies.dq_icl_dt)
        @test isfinite(tendencies.dq_rai_dt)
        @test isfinite(tendencies.dq_sno_dt)
    end

    @testset "_linearized_implicit_step - Type stability (@inferred)" begin
        ρ = FT(1.2)
        T = T_freeze + FT(7)
        q_tot = FT(0.01)
        q_lcl = FT(1e-4)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(0)
        Δt = FT(1)

        tendencies = @inferred BMT._linearized_implicit_step(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test tendencies isa NamedTuple{(:dq_lcl_dt, :dq_icl_dt, :dq_rai_dt, :dq_sno_dt), NTuple{4, FT}}
    end

    @testset "_linearized_implicit_step - rain evaporation damping vs dt" begin
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
            tendencies = BMT._linearized_implicit_step(
                BMT.Microphysics1Moment(),
                mp, tps,
                ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, dts[i],
            )
            rates[i] = tendencies.dq_rai_dt
        end

        @test all(isfinite, rates)
        @test all(r -> r < 0, rates)
        @test all(abs(rates[i + 1]) <= abs(rates[i]) for i in 1:(length(rates) - 1))
    end

    @testset "_linearized_implicit_step - Matches solved linear system" begin
        ρ = FT(1.1)
        T = T_freeze + FT(4)
        q_lcl = FT(4e-4)
        q_icl = FT(2e-4)
        q_rai = FT(3e-4)
        q_sno = FT(6e-4)
        q_tot = FT(0.014)
        Δt = FT(7)

        q_min = TDI.TD.Parameters.q_min(tps)

        src = BMT._microphysics_source_terms(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        lin = BMT._linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min, Δt)

        tendencies = BMT._linearized_implicit_step(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        invΔt = one(FT) / Δt

        q_lcl_new = q_lcl + Δt * tendencies.dq_lcl_dt
        q_icl_new = q_icl + Δt * tendencies.dq_icl_dt
        q_rai_new = q_rai + Δt * tendencies.dq_rai_dt
        q_sno_new = q_sno + Δt * tendencies.dq_sno_dt

        @test (q_lcl_new - q_lcl) * invΔt ≈ lin.M11 * q_lcl_new + lin.M12 * q_icl_new + lin.e1 atol = FT(100) * eps(FT)
        @test (q_icl_new - q_icl) * invΔt ≈ lin.M22 * q_icl_new + lin.e2 atol = FT(100) * eps(FT)
        @test (q_rai_new - q_rai) * invΔt ≈ lin.M31 * q_lcl_new + lin.M33 * q_rai_new + lin.M34 * q_sno_new atol =
            FT(100) * eps(FT)
        @test (q_sno_new - q_sno) * invΔt ≈
              lin.M41 * q_lcl_new + lin.M42 * q_icl_new + lin.M43 * q_rai_new + lin.M44 * q_sno_new + lin.e4 atol =
            FT(100) * eps(FT)
    end

    @testset "_linearized_implicit_step - Small Δt agrees with instantaneous tendency for rain evaporation" begin
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
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        avg = BMT._linearized_implicit_step(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, FT(1e-2),
        )

        @test avg.dq_lcl_dt ≈ inst.dq_lcl_dt atol = FT(1e-8)
        @test avg.dq_icl_dt ≈ inst.dq_icl_dt atol = FT(1e-8)
        @test avg.dq_sno_dt ≈ inst.dq_sno_dt atol = FT(1e-8)
        @test avg.dq_rai_dt ≈ inst.dq_rai_dt rtol = FT(1e-3)
    end

    @testset "LinearizedAverage small Δt agrees with Instantaneous (all species, warm)" begin
        # With all species active the linearized tendency should approach the
        # instantaneous one as Δt → 0.  This cross-checks _aggregate_tendencies
        # against _linearize: any sign or routing mismatch would show up here.
        # Note: Δt must be small for linearization accuracy but not so small
        # that Float32 suffers catastrophic cancellation in (q_new - q_old)/Δt.
        ρ = FT(1.2)
        T = T_freeze + FT(5)  # warm regime
        q_lcl = FT(5e-4)
        q_icl = FT(2e-4)
        q_rai = FT(3e-4)
        q_sno = FT(3e-4)
        q_tot = FT(0.012)

        inst = BMT.bulk_microphysics_tendencies(
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        lin = BMT.bulk_microphysics_tendencies(
            BMT.LinearizedAverage(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, FT(1e-2),
        )

        @test lin.dq_lcl_dt ≈ inst.dq_lcl_dt rtol = FT(5e-2)
        @test lin.dq_icl_dt ≈ inst.dq_icl_dt rtol = FT(5e-2)
        @test lin.dq_rai_dt ≈ inst.dq_rai_dt rtol = FT(5e-2)
        @test lin.dq_sno_dt ≈ inst.dq_sno_dt rtol = FT(5e-2)
    end

    @testset "LinearizedAverage small Δt agrees with Instantaneous (all species, cold)" begin
        ρ = FT(1.2)
        T = T_freeze - FT(10)  # cold regime
        q_lcl = FT(3e-4)
        q_icl = FT(5e-4)
        q_rai = FT(2e-4)
        q_sno = FT(4e-4)
        q_tot = FT(0.012)

        inst = BMT.bulk_microphysics_tendencies(
            BMT.Instantaneous(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno,
        )

        lin = BMT.bulk_microphysics_tendencies(
            BMT.LinearizedAverage(), BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, FT(1e-2),
        )

        @test lin.dq_lcl_dt ≈ inst.dq_lcl_dt rtol = FT(5e-2)
        @test lin.dq_icl_dt ≈ inst.dq_icl_dt rtol = FT(5e-2)
        @test lin.dq_rai_dt ≈ inst.dq_rai_dt rtol = FT(5e-2)
        @test lin.dq_sno_dt ≈ inst.dq_sno_dt rtol = FT(5e-2)
    end


    @testset "bulk_microphysics_tendencies(LinearizedAverage()) - Finiteness checks" begin
        ρ = FT(1.2)
        T = T_freeze - FT(5)
        q_tot = FT(0.015)
        q_lcl = FT(5e-4)
        q_icl = FT(5e-4)
        q_rai = FT(5e-4)
        q_sno = FT(5e-4)
        Δt = FT(10)

        tendencies = BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(),
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test isfinite(tendencies.dq_lcl_dt)
        @test isfinite(tendencies.dq_icl_dt)
        @test isfinite(tendencies.dq_rai_dt)
        @test isfinite(tendencies.dq_sno_dt)
    end

    @testset "bulk_microphysics_tendencies(LinearizedAverage()) - Type stability (@inferred)" begin
        ρ = FT(1.2)
        T = T_freeze + FT(7)
        q_tot = FT(0.01)
        q_lcl = FT(1e-4)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(0)
        Δt = FT(1)

        tendencies = @inferred BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(),
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test tendencies isa NamedTuple{(:dq_lcl_dt, :dq_icl_dt, :dq_rai_dt, :dq_sno_dt), NTuple{4, FT}}
    end

    @testset "bulk_microphysics_tendencies(LinearizedAverage()) - all zero inputs" begin
        ρ = FT(1.2)
        T = T_freeze + FT(5)
        q_tot = FT(0)
        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(0)
        Δt = FT(10)

        tendencies = BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(),
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test tendencies.dq_lcl_dt == FT(0)
        @test tendencies.dq_icl_dt == FT(0)
        @test tendencies.dq_rai_dt == FT(0)
        @test tendencies.dq_sno_dt == FT(0)
    end

    @testset "bulk_microphysics_tendencies(LinearizedAverage()) - nsub=1 matches single-substep solver" begin
        ρ = FT(1.1)
        T = T_freeze + FT(4)
        q_lcl = FT(4e-4)
        q_icl = FT(2e-4)
        q_rai = FT(3e-4)
        q_sno = FT(6e-4)
        q_tot = FT(0.014)
        Δt = FT(7)

        single = BMT._linearized_implicit_step(
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        substepped = BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(),
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, 1,
        )

        @test substepped.dq_lcl_dt ≈ single.dq_lcl_dt atol = FT(100) * eps(FT)
        @test substepped.dq_icl_dt ≈ single.dq_icl_dt atol = FT(100) * eps(FT)
        @test substepped.dq_rai_dt ≈ single.dq_rai_dt atol = FT(100) * eps(FT)
        @test substepped.dq_sno_dt ≈ single.dq_sno_dt atol = FT(100) * eps(FT)
    end

    @testset "bulk_microphysics_tendencies(LinearizedAverage()) - Warm pure snow melt keeps expected signs" begin
        ρ = FT(1.0)
        T = T_freeze + FT(5)
        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(1e-3)
        q_vap_sat = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        q_tot = q_vap_sat + q_sno
        Δt = FT(10)

        tendencies = BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(),
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
        )

        @test tendencies.dq_sno_dt < FT(0)
        @test tendencies.dq_rai_dt > FT(0)
        @test isfinite(tendencies.dq_sno_dt)
        @test isfinite(tendencies.dq_rai_dt)
    end

    @testset "bulk_microphysics_tendencies(LinearizedAverage()) - More substeps do not change simple rain-only case much" begin
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

        t1 = BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(),
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, 1,
        )

        t10 = BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(),
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, 10,
        )

        @test t10.dq_lcl_dt ≈ t1.dq_lcl_dt atol = FT(1e-10)
        @test t10.dq_icl_dt ≈ t1.dq_icl_dt atol = FT(1e-10)
        @test t10.dq_sno_dt ≈ t1.dq_sno_dt atol = FT(1e-10)
        @test t10.dq_rai_dt ≈ t1.dq_rai_dt rtol = FT(1e-2)
    end

    @testset "bulk_microphysics_tendencies(LinearizedAverage()) - Substepping remains finite near freezing" begin
        ρ = FT(1.2)
        T = T_freeze + FT(0.01)
        q_tot = FT(0.015)
        q_lcl = FT(1e-3)
        q_icl = FT(0)
        q_rai = FT(0)
        q_sno = FT(5e-4)
        Δt = FT(20)

        tendencies = BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(),
            BMT.Microphysics1Moment(),
            mp, tps,
            ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, 20,
        )

        @test isfinite(tendencies.dq_lcl_dt)
        @test isfinite(tendencies.dq_icl_dt)
        @test isfinite(tendencies.dq_rai_dt)
        @test isfinite(tendencies.dq_sno_dt)
    end

    @testset "_donor_limiter_scale: realized transfer of the scaled decays is f times the unscaled one" begin
        Δt = FT(60)
        q, Dpc, Dcol, f = FT(1e-3), FT(0.05), FT(0.01), FT(0.3)
        s = BMT._donor_limiter_scale(f, Dpc, Dcol, Δt)
        realized(sc) = q * sc * Dpc * Δt / (1 + (Dcol + sc * Dpc) * Δt)
        @test realized(s) ≈ f * realized(one(FT)) rtol = FT(1e-5)
        @test BMT._donor_limiter_scale(one(FT), Dpc, Dcol, Δt) == one(FT)
        @test BMT._donor_limiter_scale(f, FT(0), Dcol, Δt) == f
        @test (@inferred BMT._donor_limiter_scale(f, Dpc, Dcol, Δt)) isa FT
    end

    @testset "LinearizedAverage - latent heating limiter bounds every substep and conserves water" begin
        # Glaciating updraft state of the AMIP blow-up (PrescribedIceNumber, N_0 = 5e8); with the bound at
        # 0.005 K/s (0.2 K per 40 s substep) the
        # heating of each substep must not exceed the bound (up to the implicit-decay approximation),
        # all phase-change transfers scale together, and total water is unchanged. With the bound
        # disabled (Inf) the limiter factor is 1.
        mp_lim = stiff_prescribed_ice_params(FT; max_latent_heating_rate = 0.005)
        mp_off = stiff_prescribed_ice_params(FT)
        @test mp_lim.max_latent_heating_rate == FT(0.005)
        @test isinf(mp_off.max_latent_heating_rate)
        ρ = FT(0.6884);
        T = FT(266.55)
        q_tot = FT(8.129e-3);
        q_lcl = FT(6.102e-4);
        q_icl = FT(8.345e-4);
        q_rai = FT(1.408e-5);
        q_sno = FT(1.634e-3)
        Δts = FT(40)
        Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
        Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)
        heating(r) = (Lv_over_cp * (r.dq_lcl_dt + r.dq_rai_dt) + Ls_over_cp * (r.dq_icl_dt + r.dq_sno_dt)) * Δts
        r_lim = BMT._linearized_implicit_step(
            BMT.Microphysics1Moment(),
            mp_lim,
            tps,
            ρ,
            T,
            FT(0),
            q_tot,
            q_lcl,
            q_icl,
            q_rai,
            q_sno,
            Δts,
        )
        r_off = BMT._linearized_implicit_step(
            BMT.Microphysics1Moment(),
            mp_off,
            tps,
            ρ,
            T,
            FT(0),
            q_tot,
            q_lcl,
            q_icl,
            q_rai,
            q_sno,
            Δts,
        )
        @test heating(r_off) > FT(0.005) * Δts          # the limiter is needed for this state
        @test heating(r_lim) <= FT(0.005) * Δts * FT(1.05)
        @test heating(r_lim) > FT(0.5) * FT(0.005) * Δts # and it does not switch the processes off
        # the limiter factor itself
        args = (ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δts)
        f_lim = BMT._linearized_implicit_step_factors(BMT.Microphysics1Moment(), mp_lim, tps, args...).f_lim
        f_off = BMT._linearized_implicit_step_factors(BMT.Microphysics1Moment(), mp_off, tps, args...).f_lim
        @test FT(0) < f_lim < FT(1)
        @test f_off == FT(1)
        # the realized heating sits on the bound (both solves share the linearization; the scale factors are exact
        # for each donor alone, so the bound is met up to the coupling between donors)
        @test heating(r_lim) ≈ FT(0.005) * Δts rtol = FT(0.05)
        # gentle state: limiter inactive
        T_w = T_freeze + FT(10);
        q_sl = TDI.saturation_vapor_specific_content_over_liquid(tps, T_w, FT(1.1))
        gentle = (FT(1.1), T_w, FT(0), q_sl * FT(1.001) + FT(1e-4), FT(1e-4), FT(0), FT(0), FT(0), Δts)
        @test BMT._linearized_implicit_step_factors(BMT.Microphysics1Moment(), mp_lim, tps, gentle...).f_lim == FT(1)
        # water conservation: the tendencies redistribute water among the species (vapor implied)
        for r in (r_lim, r_off)
            @test all(isfinite, (r.dq_lcl_dt, r.dq_icl_dt, r.dq_rai_dt, r.dq_sno_dt))
            @test q_lcl + r.dq_lcl_dt * Δts >= FT(0) && q_icl + r.dq_icl_dt * Δts >= FT(0)
            @test q_rai + r.dq_rai_dt * Δts >= FT(0) && q_sno + r.dq_sno_dt * Δts >= FT(0)
        end
    end

    @testset "LinearizedAverage - latent heating limiter in a melting layer (fusion-only cooling)" begin
        # 276 K, saturated, snow 3 g/kg and rain 1 g/kg: melting cools ~1 K in a 60 s substep. With the bound
        # 0.002 K/s the realized cooling must sit on the bound (-0.12 K), the melted water must still land in
        # rain, and the vapor check must stay inactive (no vapor sources).
        mp_lim = stiff_prescribed_ice_params(FT; max_latent_heating_rate = 0.002)
        mp_off = stiff_prescribed_ice_params(FT)
        ρ = FT(1);
        T = FT(276)
        q_sl = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
        q_lcl = FT(0);
        q_icl = FT(0);
        q_rai = FT(1e-3);
        q_sno = FT(3e-3)
        q_tot = q_sl + q_rai + q_sno
        Δt = FT(60)
        Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
        Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)
        heating(r) = (Lv_over_cp * (r.dq_lcl_dt + r.dq_rai_dt) + Ls_over_cp * (r.dq_icl_dt + r.dq_sno_dt)) * Δt
        args = (ρ, T, FT(0), q_tot, q_lcl, q_icl, q_rai, q_sno, Δt)
        st_lim = BMT._linearized_implicit_step_factors(BMT.Microphysics1Moment(), mp_lim, tps, args...)
        st_off = BMT._linearized_implicit_step_factors(BMT.Microphysics1Moment(), mp_off, tps, args...)
        r_lim = st_lim.rates;
        r_off = st_off.rates
        @test heating(r_off) < -FT(0.5)                        # the bound is needed
        @test st_off.f_lim == FT(1) && st_off.α_cap == FT(1)
        @test FT(0) < st_lim.f_lim < FT(1)
        @test st_lim.α_cap == FT(1)
        @test heating(r_lim) ≈ -FT(0.002) * Δt rtol = FT(0.02)  # realized cooling on the bound
        @test r_lim.dq_sno_dt < FT(0) && r_lim.dq_rai_dt > FT(0)
        @test r_lim.dq_rai_dt ≈ -r_lim.dq_sno_dt rtol = FT(0.02)  # melted snow lands in rain
        @test q_sno + r_lim.dq_sno_dt * Δt >= FT(0)
    end

    @testset "LinearizedAverage - default limiter (2 K/min) is inactive for an ordinary updraft" begin
        mp_def = stiff_prescribed_ice_params(FT; max_latent_heating_rate = 2 / 60)
        Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
        Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)
        # warm updraft, 2 % supersaturated over liquid, one 60 s substep
        ρ = FT(1);
        T = FT(285)
        q_sl = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
        q_lcl = FT(3e-4);
        q_rai = FT(1e-4)
        q_tot = FT(1.02) * q_sl + q_lcl + q_rai
        Δt = FT(60)
        heating(r, Δt) = (Lv_over_cp * (r.dq_lcl_dt + r.dq_rai_dt) + Ls_over_cp * (r.dq_icl_dt + r.dq_sno_dt)) * Δt
        args = (ρ, T, FT(0), q_tot, q_lcl, FT(0), q_rai, FT(0), Δt)
        st = BMT._linearized_implicit_step_factors(BMT.Microphysics1Moment(), mp_def, tps, args...)
        @test st.f_lim == FT(1)
        r = st.rates
        ref = BMT.bulk_microphysics_tendencies(
            BMT.LinearizedAverage(),
            BMT.Microphysics1Moment(),
            mp_def,
            tps,
            args...,
            60,
        )
        @test FT(0) < heating(r, Δt) < FT(2 / 60) * Δt
        @test heating(r, Δt) ≈ heating(ref, Δt) rtol = FT(0.02)
        T1 = T + heating(r, Δt)
        q_v = q_tot - (q_lcl + r.dq_lcl_dt * Δt) - (q_rai + r.dq_rai_dt * Δt)
        @test q_v ≈ TDI.saturation_vapor_specific_content_over_liquid(tps, T1, ρ) rtol = FT(2e-3)  # Γ-consistent
        @test r.dq_lcl_dt > FT(0) && r.dq_rai_dt > FT(0)       # cloud liquid condenses; rain grows by accretion (its 1M vapor exchange is evaporation only)
    end

    @testset "LinearizedAverage - negative condensate inputs are clamped; limiter parameter is validated" begin
        mp = stiff_prescribed_ice_params(FT)
        ρ = FT(0.8);
        T = FT(250)
        q_tot = FT(1.1) * TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        r_neg = BMT._linearized_implicit_step(
            BMT.Microphysics1Moment(), mp, tps, ρ, T, FT(0), q_tot, -FT(1e-6), -FT(1e-7), -FT(1e-8), -FT(1e-9),
            FT(60),
        )
        r_zero = BMT._linearized_implicit_step(
            BMT.Microphysics1Moment(), mp, tps, ρ, T, FT(0), q_tot, FT(0), FT(0), FT(0), FT(0), FT(60),
        )
        @test all(isfinite, values(r_neg))
        @test all(k -> r_neg[k] == r_zero[k], keys(r_zero))   # no spurious source or sink from a negative pool
        @test r_neg.dq_icl_dt > FT(0)                          # deposition onto (nucleating) ice proceeds
        for bad in (0.0, -1.0, NaN)
            @test_throws ArgumentError stiff_prescribed_ice_params(FT; max_latent_heating_rate = bad)
        end
        @test_throws ArgumentError CMP._validated_max_latent_heating_rate(FT(0))
        @test CMP._validated_max_latent_heating_rate(FT(Inf)) == FT(Inf)
    end

    @testset "LinearizedAverage - latent heating limiter holds when rain and snow feed each other" begin
        # 274 K with cloud ice, rain and snow: rain freezes on the ice (rai → sno) while the snow melts
        # (sno → rai), two large opposing fusion transfers through pools that refill each other. The
        # per-donor scaling alone can then miss the bound; the uniform fallback must enforce it.
        mp_lim = stiff_prescribed_ice_params(FT; max_latent_heating_rate = 0.002)
        ρ = FT(1);
        T = FT(274)
        q_sl = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
        q = (FT(1e-4), FT(5e-4), FT(1e-4), FT(1e-3))
        q_tot = FT(1.03) * q_sl + sum(q)
        Δts = FT(40)
        Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
        Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)
        heating(r) = (Lv_over_cp * (r.dq_lcl_dt + r.dq_rai_dt) + Ls_over_cp * (r.dq_icl_dt + r.dq_sno_dt)) * Δts
        Tloc = T
        engaged = false
        for _ in 1:3
            st = BMT._linearized_implicit_step_factors(
                BMT.Microphysics1Moment(),
                mp_lim,
                tps,
                ρ,
                Tloc,
                FT(0),
                q_tot,
                q...,
                Δts,
            )
            r = st.rates
            @test abs(heating(r)) <= FT(0.002) * Δts * (one(FT) + FT(1e-3))
            engaged |= st.f_lim < one(FT)
            q = (
                q[1] + r.dq_lcl_dt * Δts,
                q[2] + r.dq_icl_dt * Δts,
                q[3] + r.dq_rai_dt * Δts,
                q[4] + r.dq_sno_dt * Δts,
            )
            @test all(>=(-eps(FT)), q)
            Tloc += heating(r)
        end
        @test engaged                                          # the bound was needed in at least one substep
    end

    @testset "LinearizedAverage - vapor check floor follows the saturation curve for a large substep heating" begin
        # RH_ice 3 at 250 K with stiff ice: one 120 s substep heats ~3.8 K; the floor of the vapor check is
        # linear in ΔT, so the vapor ends a few % below the new ice saturation (≤ 6 %; the default 2 K/min
        # limiter keeps the substep heating, and hence this error, much smaller in production).
        mp = stiff_prescribed_ice_params(FT)
        ρ = FT(0.8);
        T = FT(250)
        q_si = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        q_icl = FT(1e-4)
        q_tot = FT(3) * q_si + q_icl
        Δt = FT(120)
        Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
        Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)
        st = BMT._linearized_implicit_step_factors(
            BMT.Microphysics1Moment(),
            mp,
            tps,
            ρ,
            T,
            FT(0),
            q_tot,
            FT(0),
            q_icl,
            FT(0),
            FT(0),
            Δt,
        )
        r = st.rates
        ΔT = (Lv_over_cp * (r.dq_lcl_dt + r.dq_rai_dt) + Ls_over_cp * (r.dq_icl_dt + r.dq_sno_dt)) * Δt
        q_v = q_tot - (q_icl + r.dq_icl_dt * Δt) - (r.dq_lcl_dt + r.dq_rai_dt + r.dq_sno_dt) * Δt
        @test ΔT > FT(3)
        @test st.α_cap < one(FT)                               # the check is what stops the deposition
        q_si1 = TDI.saturation_vapor_specific_content_over_ice(tps, T + ΔT, ρ)
        @test FT(0.94) * q_si1 <= q_v <= q_si1 * (one(FT) + FT(1e-3))
    end

    @testset "LinearizedAverage - hard positivity of every tracer under stiff processes" begin
        # stiff processes on tiny pools from an almost dry or very moist atmosphere, huge substeps:
        # condensate pools stay ≥ 0 (implicit decays) and the vapor stays ≥ 0 (uniform guard)
        for Δt in (FT(1), FT(120), FT(3600))
            mp = stiff_prescribed_ice_params(FT)
            for (T, ρ, q) in (
                (FT(250), FT(0.8), (FT(1e-3), FT(1e-6), FT(0), FT(0))),        # liquid in ice-subsaturated dry air
                (FT(230), FT(0.5), (FT(3e-3), FT(1e-4), FT(1e-4), FT(1e-4))),  # homogeneous freezing + deposition
                (FT(274), FT(1.0), (FT(1e-4), FT(5e-4), FT(1e-4), FT(3e-3))),  # melting churn above freezing
            )
                for q_v in (FT(0), FT(1e-7), FT(3) * TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ))
                    q_tot = q_v + sum(q)
                    r = BMT._linearized_implicit_step(BMT.Microphysics1Moment(), mp, tps, ρ, T, FT(0), q_tot, q..., Δt)
                    q_new = q .+ (r.dq_lcl_dt, r.dq_icl_dt, r.dq_rai_dt, r.dq_sno_dt) .* Δt
                    @test all(isfinite, q_new)
                    @test all(>=(-eps(FT)), q_new)
                    @test q_tot - sum(q_new) >= -FT(4) * eps(FT) * q_tot     # vapor never negative
                end
            end
        end
    end

    @testset "_linearized_implicit_step_factors - Type stability (@inferred)" begin
        mp = stiff_prescribed_ice_params(FT)
        st = @inferred BMT._linearized_implicit_step_factors(
            BMT.Microphysics1Moment(), mp, tps, FT(0.8), FT(250), FT(0), FT(1e-3), FT(1e-4), FT(1e-5), FT(1e-5),
            FT(1e-5), FT(60),
        )
        @test st.α_cap isa FT && st.f_lim isa FT && st.g_uniform isa FT
        @test all(v -> v isa FT, values(st.rates))
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
            dn_lcl_activation_dt::FT,
        }
        # Ice tendencies should be zero for 2M mode
        @test tendencies.dq_ice_dt == FT(0)
        @test tendencies.dq_rim_dt == FT(0)
        @test tendencies.db_rim_dt == FT(0)
        # Default NoActivation scheme yields zero activation tendency
        @test tendencies.dn_lcl_activation_dt == FT(0)
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
        state = CM.P3Scheme.P3State(p3, L_ice, N_ice, F_rim, ρ_rim)
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
        state = CM.P3Scheme.P3State(p3, L_ice, N_ice, F_rim, ρ_rim)
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
        state = CM.P3Scheme.P3State(p3, L_ice, N_ice, F_rim, ρ_rim)
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
        state = CM.P3Scheme.P3State(p3, L_ice, N_ice, F_rim, ρ_rim)
        logλ = CM.P3Scheme.get_distribution_logλ(state)

        tendencies = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(), mp, tps, ρ, T, q_tot,
            q_lcl, n_lcl, q_rai, n_rai,
            q_ice, n_ice, q_rim, b_rim, logλ,
        )

        # Check that we get the expected NamedTuple type
        @test tendencies isa @NamedTuple{
            dq_lcl_dt::FT, dn_lcl_dt::FT, dq_rai_dt::FT, dn_rai_dt::FT,
            dq_ice_dt::FT, dn_ice_dt::FT, dq_rim_dt::FT, db_rim_dt::FT,
            dn_lcl_activation_dt::FT,
        }
    end
end

@testset "Bulk Microphysics Tendencies ($FT)" for FT in (Float64, Float32)
    test_bulk_microphysics_0m_tendencies(FT)
    test_bulk_microphysics_1m_tendencies(FT)
    test_linearized_bulk_microphysics_1m_tendencies(FT)
    test_bulk_microphysics_2m_tendencies(FT)
    test_bulk_microphysics_p3_tendencies(FT)
end
nothing
