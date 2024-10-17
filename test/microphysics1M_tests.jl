import Test as TT

import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaParams as CP

import Thermodynamics.Parameters as TDP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CMC
import CloudMicrophysics.Microphysics1M as CM1

@info "Microphysics 1M Tests"

function test_microphysics1M(FT)

    tps = TD.Parameters.ThermodynamicsParameters(FT)
    aps = CMP.AirProperties(FT)

    rain = CMP.Rain(FT)
    snow = CMP.Snow(FT)
    ice = CMP.CloudIce(FT)
    liquid = CMP.CloudLiquid(FT)

    ce = CMP.CollisionEff(FT)

    Ch2022 = CMP.Chen2022VelType(FT)
    blk1mvel = CMP.Blk1MVelType(FT)

    TT.@testset "1M_microphysics RainFallSpeed" begin
        # eq. 5d in [Grabowski1996](@cite)
        function terminal_velocity_empir(
            q_rai::FT,
            q_tot::FT,
            ρ::FT,
            ρ_air_ground::FT,
        ) where {FT <: Real}
            rr = q_rai / (1 - q_tot)
            vel =
                FT(14.34) * ρ_air_ground^FT(0.5) * ρ^-FT(0.3654) * rr^FT(0.1346)
            return vel
        end

        # some example values
        q_rain_range = range(FT(1e-8), stop = FT(5e-3), length = 10)
        ρ_air, q_tot, ρ_air_ground = FT(1.2), FT(20 * 1e-3), FT(1.22)

        for q_rai in q_rain_range
            TT.@test CM1.terminal_velocity(rain, blk1mvel.rain, ρ_air, q_rai) ≈
                     terminal_velocity_empir(q_rai, q_tot, ρ_air, ρ_air_ground) atol =
                0.2 * terminal_velocity_empir(q_rai, q_tot, ρ_air, ρ_air_ground)
        end
    end

    TT.@testset "1M_microphysics RadarReflectivity" begin

        # some example values
        ρ_air = FT(1)
        q_rai = FT(0.18e-3)

        TT.@test CM1.radar_reflectivity(rain, q_rai, ρ_air) ≈ FT(12.17) atol =
            0.2

        q_rai = FT(0.89e-4)

        TT.@test CM1.radar_reflectivity(rain, q_rai, ρ_air) ≈ FT(6.68) atol =
            0.2

    end

    TT.@testset "1M_microphysics - Chen 2022 rain terminal velocity" begin
        #setup
        ρ = FT(1.2)
        q_rai = FT(5e-4)
        N_rai = FT(1e4)

        #action
        vt_rai = CM1.terminal_velocity(rain, Ch2022.rain, ρ, q_rai)
        v_bigger = CM1.terminal_velocity(rain, Ch2022.rain, ρ, q_rai * 2)

        #test
        TT.@test vt_rai ≈ 3.0721397260869705 rtol = 2e-6
        TT.@test CM1.terminal_velocity(rain, Ch2022.rain, ρ, FT(0))[1] ≈ 0 atol =
            eps(FT)
        TT.@test v_bigger > vt_rai
    end

    TT.@testset "1M_microphysics - Chen 2022 snow terminal velocity" begin
        #setup
        ρ = FT(1.1)
        q_sno = FT(5e-4)
        N_rai = FT(1e4)

        #action
        vt_sno = CM1.terminal_velocity(snow, Ch2022.snow_ice, ρ, q_sno)
        v_bigger = CM1.terminal_velocity(snow, Ch2022.snow_ice, ρ, q_sno * 2)

        #test
        TT.@test vt_sno ≈ 2.0622134820974636 rtol = 2e-6
        TT.@test CM1.terminal_velocity(snow, Ch2022.snow_ice, ρ, FT(0)) ≈ 0 atol =
            eps(FT)
        TT.@test v_bigger > vt_sno
    end


    TT.@testset "RainAutoconversion" begin

        q_liq_threshold = rain.acnv1M.q_threshold
        τ_acnv_rai = rain.acnv1M.τ

        q_liq_small = FT(0.5) * q_liq_threshold
        TT.@test CM1.conv_q_liq_to_q_rai(rain.acnv1M, q_liq_small) == FT(0)

        TT.@test CM1.conv_q_liq_to_q_rai(rain.acnv1M, q_liq_small, true) ≈
                 FT(0.0) atol = 0.15 * q_liq_threshold / τ_acnv_rai

        q_liq_big = FT(1.5) * q_liq_threshold
        TT.@test CM1.conv_q_liq_to_q_rai(rain.acnv1M, q_liq_big) ≈
                 FT(0.5) * q_liq_threshold / τ_acnv_rai

        TT.@test CM1.conv_q_liq_to_q_rai(rain.acnv1M, q_liq_big, true) ≈
                 FT(0.5) * q_liq_threshold / τ_acnv_rai atol =
            FT(0.15) * q_liq_threshold / τ_acnv_rai

    end

    TT.@testset "SnowAutoconversionNoSupersat" begin

        q_ice_threshold = snow.acnv1M.q_threshold
        τ_acnv_sno = snow.acnv1M.τ

        q_ice_small = FT(0.5) * q_ice_threshold
        TT.@test CM1.conv_q_ice_to_q_sno_no_supersat(
            snow.acnv1M,
            q_ice_small,
        ) == FT(0)

        TT.@test CM1.conv_q_ice_to_q_sno_no_supersat(
            snow.acnv1M,
            q_ice_small,
            true,
        ) ≈ FT(0.0) atol = FT(0.15) * q_ice_threshold / τ_acnv_sno

        q_ice_big = FT(1.5) * q_ice_threshold
        TT.@test CM1.conv_q_ice_to_q_sno_no_supersat(snow.acnv1M, q_ice_big) ≈
                 FT(0.5) * q_ice_threshold / τ_acnv_sno

        TT.@test CM1.conv_q_ice_to_q_sno_no_supersat(
            snow.acnv1M,
            q_ice_big,
            true,
        ) ≈ FT(0.5) * q_ice_threshold / τ_acnv_sno atol =
            FT(0.15) * q_ice_threshold / τ_acnv_sno

        TT.@test CM1.conv_q_ice_to_q_sno_no_supersat(snow.acnv1M, q_ice_big) ==
                 CM1.conv_q_ice_to_q_sno_no_supersat(
            snow.acnv1M,
            q_ice_big,
            false,
        )

    end

    TT.@testset "SnowAutoconversion" begin

        ρ = FT(1.0)

        # above freezing temperatures -> no snow
        q = TD.PhasePartition(FT(15e-3), FT(2e-3), FT(1e-3))
        T = FT(273.15 + 30)
        TT.@test CM1.conv_q_ice_to_q_sno(ice, aps, tps, q, ρ, T) == FT(0)

        # no ice -> no snow
        q = TD.PhasePartition(FT(15e-3), FT(2e-3), FT(0))
        T = FT(273.15 - 30)
        TT.@test CM1.conv_q_ice_to_q_sno(ice, aps, tps, q, ρ, T) == FT(0)

        # no supersaturation -> no snow
        T = FT(273.15 - 5)
        q_sat_ice = TD.q_vap_saturation_generic(tps, T, ρ, TD.Ice())
        q = TD.PhasePartition(q_sat_ice, FT(2e-3), FT(3e-3))
        TT.@test CM1.conv_q_ice_to_q_sno(ice, aps, tps, q, ρ, T) == FT(0)

        # TODO - coudnt find a plot of what it should be from the original paper
        # just chacking if the number stays the same
        T = FT(273.15 - 10)
        q_vap = FT(1.02) * TD.q_vap_saturation_generic(tps, T, ρ, TD.Ice())
        q_liq = FT(0)
        q_ice = FT(0.03) * q_vap
        q = TD.PhasePartition(q_vap + q_liq + q_ice, q_liq, q_ice)

        TT.@test CM1.conv_q_ice_to_q_sno(ice, aps, tps, q, ρ, T) ≈
                 FT(1.8512022335645584e-9)
    end

    TT.@testset "RainLiquidAccretion" begin

        # eq. 5b in [Grabowski1996](@cite)
        function accretion_empir(
            q_rai::FT,
            q_liq::FT,
            q_tot::FT,
        ) where {FT <: Real}
            rr = q_rai / (FT(1) - q_tot)
            rl = q_liq / (FT(1) - q_tot)
            return FT(2.2) * rl * rr^FT(7 / 8)
        end

        # some example values
        q_rain_range = range(FT(1e-8), stop = FT(5e-3), length = 10)
        ρ_air, q_liq, q_tot = FT(1.2), FT(5e-4), FT(20e-3)

        for q_rai in q_rain_range
            if q_rai > eps(FT)
                TT.@test CM1.accretion(
                    liquid,
                    rain,
                    blk1mvel.rain,
                    ce,
                    q_liq,
                    q_rai,
                    ρ_air,
                ) ≈ accretion_empir(q_rai, q_liq, q_tot) atol =
                    (0.1 * accretion_empir(q_rai, q_liq, q_tot))
            else
                TT.@test CM1.accretion(
                    liquid,
                    rain,
                    blk1mvel.rain,
                    ce,
                    q_liq,
                    q_rai,
                    ρ_air,
                ) ≈ FT(0) atol = FT(0.2)
            end
        end
    end

    TT.@testset "Accretion" begin
        # TODO - coudnt find a plot of what it should be from the original paper
        # just chacking if the number stays the same

        # some example values
        ρ = FT(1.2)
        q_tot = FT(20e-3)
        q_ice = FT(5e-4)
        q_sno = FT(5e-4)
        q_liq = FT(5e-4)
        q_rai = FT(5e-4)

        TT.@test CM1.accretion(
            liquid,
            rain,
            blk1mvel.rain,
            ce,
            q_liq,
            q_rai,
            ρ,
        ) ≈ FT(1.4150106417043544e-6)
        TT.@test CM1.accretion(ice, snow, blk1mvel.snow, ce, q_ice, q_sno, ρ) ≈
                 FT(2.453070979562392e-7)
        TT.@test CM1.accretion(
            liquid,
            snow,
            blk1mvel.snow,
            ce,
            q_liq,
            q_sno,
            ρ,
        ) ≈ FT(2.453070979562392e-7)
        TT.@test CM1.accretion(ice, rain, blk1mvel.rain, ce, q_ice, q_rai, ρ) ≈
                 FT(1.768763302130443e-6)

        TT.@test CM1.accretion_rain_sink(
            rain,
            ice,
            blk1mvel.rain,
            ce,
            q_ice,
            q_rai,
            ρ,
        ) ≈ FT(3.085229094251214e-5)

        TT.@test CM1.accretion_snow_rain(
            snow,
            rain,
            blk1mvel.snow,
            blk1mvel.rain,
            ce,
            q_sno,
            q_rai,
            ρ,
        ) ≈ FT(2.1705865794293408e-4)
        TT.@test CM1.accretion_snow_rain(
            rain,
            snow,
            blk1mvel.rain,
            blk1mvel.snow,
            ce,
            q_rai,
            q_sno,
            ρ,
        ) ≈ FT(6.0118801860768854e-5)
    end

    TT.@testset "RainEvaporation" begin

        # eq. 5c in [Grabowski1996](@cite)
        function rain_evap_empir(
            tps::TDP.ThermodynamicsParameters,
            q_rai::FT,
            q::TD.PhasePartition,
            T::FT,
            p::FT,
            ρ::FT,
        ) where {FT <: Real}

            q_sat = TD.q_vap_saturation_generic(tps, T, ρ, TD.Liquid())
            q_vap = q.tot - q.liq
            rr = q_rai / (1 - q.tot)
            rv_sat = q_sat / (1 - q.tot)
            S = q_vap / q_sat - 1

            ag, bg = FT(5.4 * 1e2), FT(2.55 * 1e5)
            G = FT(1) / (ag + bg / p / rv_sat) / ρ

            av, bv = FT(1.6), FT(124.9)
            F =
                av * (ρ / FT(1e3))^FT(0.525) * rr^FT(0.525) +
                bv * (ρ / FT(1e3))^FT(0.7296) * rr^FT(0.7296)

            return 1 / (1 - q.tot) * S * F * G
        end

        # example values
        T, p = FT(273.15 + 15), FT(90000)
        ϵ = 1 / TD.Parameters.molmass_ratio(tps)
        p_sat = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))
        q_rain_range = range(FT(1e-8), stop = FT(5e-3), length = 10)
        q_tot = FT(15e-3)
        q_vap = FT(0.15) * q_sat
        q_ice = FT(0)
        q_liq = q_tot - q_vap - q_ice
        q = TD.PhasePartition(q_tot, q_liq, q_ice)
        R = TD.gas_constant_air(tps, q)
        ρ = p / R / T

        for q_rai in q_rain_range
            if q_rai > eps(FT)
                TT.@test CM1.evaporation_sublimation(
                    rain,
                    blk1mvel.rain,
                    aps,
                    tps,
                    q,
                    q_rai,
                    ρ,
                    T,
                ) ≈ rain_evap_empir(tps, q_rai, q, T, p, ρ) atol =
                    -0.5 * rain_evap_empir(tps, q_rai, q, T, p, ρ)
            else
                TT.@test CM1.evaporation_sublimation(
                    rain,
                    blk1mvel.rain,
                    aps,
                    tps,
                    q,
                    q_rai,
                    ρ,
                    T,
                ) ≈ FT(0) atol = FT(0.2)
            end
        end

        # no condensational growth for rain
        T, p = FT(273.15 + 15), FT(90000)
        ϵ = 1 / TD.Parameters.molmass_ratio(tps)
        p_sat = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))
        q_rai = FT(1e-4)
        q_tot = FT(15e-3)
        q_vap = FT(1.15) * q_sat
        q_ice = FT(0)
        q_liq = q_tot - q_vap - q_ice
        q = TD.PhasePartition(q_tot, q_liq, q_ice)
        R = TD.gas_constant_air(tps, q)
        ρ = p / R / T

        TT.@test CM1.evaporation_sublimation(
            rain,
            blk1mvel.rain,
            aps,
            tps,
            q,
            q_rai,
            ρ,
            T,
        ) ≈ FT(0)

    end

    TT.@testset "SnowSublimation" begin
        # TODO - coudnt find a plot of what it should be from the original paper
        # just chacking if the number stays the same

        cnt = 0
        ref_val = [
            FT(-1.9756907119482267e-7),
            FT(1.9751292385808357e-7),
            FT(-1.6641552112891826e-7),
            FT(1.663814937710236e-7),
        ]
        # some example values
        for T in [FT(273.15 + 2), FT(273.15 - 2)]
            p = FT(90000)
            ϵ = 1 / TD.Parameters.molmass_ratio(tps)
            p_sat = TD.saturation_vapor_pressure(tps, T, TD.Ice())
            q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))

            for eps in [FT(0.95), FT(1.05)]
                cnt += 1

                q_tot = eps * q_sat
                q_ice = FT(0)
                q_liq = FT(0)
                q = TD.PhasePartition(q_tot, q_liq, q_ice)

                q_sno = FT(1e-4)

                R = TD.gas_constant_air(tps, q)
                ρ = p / R / T

                TT.@test CM1.evaporation_sublimation(
                    snow,
                    blk1mvel.snow,
                    aps,
                    tps,
                    q,
                    q_sno,
                    ρ,
                    T,
                ) ≈ ref_val[cnt]

            end
        end
    end

    TT.@testset "SnowMelt" begin

        # TODO - find a good reference to compare with
        T = FT(273.15 + 2)
        ρ = FT(1.2)
        q_sno = FT(1e-4)
        TT.@test CM1.snow_melt(snow, blk1mvel.snow, aps, tps, q_sno, ρ, T) ≈
                 FT(9.518235437405256e-6)

        # no snow -> no snow melt
        T = FT(273.15 + 2)
        ρ = FT(1.2)
        q_sno = FT(0)
        TT.@test CM1.snow_melt(snow, blk1mvel.snow, aps, tps, q_sno, ρ, T) ≈
                 FT(0)

        # T < T_freeze -> no snow melt
        T = FT(273.15 - 2)
        ρ = FT(1.2)
        q_sno = FT(1e-4)
        TT.@test CM1.snow_melt(snow, blk1mvel.snow, aps, tps, q_sno, ρ, T) ≈
                 FT(0)

    end
end

println("Testing Float64")
test_microphysics1M(Float64)

println("Testing Float32")
test_microphysics1M(Float32)
