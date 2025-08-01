import Test as TT

import ClimaParams as CP

import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CMC
import CloudMicrophysics.Microphysics1M as CM1

function test_microphysics1M(FT)

    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
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
        vt_sno = CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, q_sno)
        v_bigger = CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, q_sno * 2)

        #test
        TT.@test vt_sno ≈ 0.5151154754853068 rtol = 3e-6
        TT.@test CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, FT(0)) ≈ 0 atol =
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
        qᵣ = FT(1e-4)
        qₛ = FT(1e-4)
        T₀ = FT(273.15)

        # above freezing temperatures -> no snow
        qᵥ = FT(15e-3)
        qₗ = FT(2e-3)
        qᵢ = FT(1e-3)
        qₜ = qᵥ + qₗ + qᵢ + qᵣ + qₛ
        T = T₀ + FT(30)
        TT.@test CM1.conv_q_ice_to_q_sno(ice, aps, tps, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T) == FT(0)

        # no ice -> no snow
        qᵢ = FT(0)
        qₜ = qᵥ + qₗ + qᵢ + qᵣ + qₛ
        T = T₀ - FT(30)
        TT.@test CM1.conv_q_ice_to_q_sno(ice, aps, tps, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T) == FT(0)

        # no supersaturation -> no snow
        T = T₀ - FT(5)
        qᵢ = FT(3e-3)
        q_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        qₜ = q_sat_ice
        TT.@test CM1.conv_q_ice_to_q_sno(ice, aps, tps, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T) ≈ FT(0)

        # Coudnt find a plot of what it should be from the original paper
        # just checking if the number stays the same
        T = T₀ - FT(10)
        qᵥ = FT(1.02) * TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        qₗ = FT(0)
        qᵢ = FT(0.03) * qᵥ
        qₜ = qᵥ + qₗ + qᵢ + qᵣ + qₛ
        ref = FT(2.5405159487552133e-9)
        TT.@test CM1.conv_q_ice_to_q_sno(ice, aps, tps, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T) ≈ ref
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
        ) ≈ FT(3.590060148920766e-5)

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
            tps::TDI.TD.Parameters.ThermodynamicsParameters,
            q_rai::FT,
            q_tot::FT,
            q_lcl::FT,
            q_ice::FT,
            T::FT,
            p::FT,
            ρ::FT,
        ) where {FT <: Real}

            q_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
            q_vap = q_tot - q_lcl - q_ice - q_rai
            rr = q_rai / (1 - q_tot)
            rv_sat = q_sat / (1 - q_tot)
            S = q_vap / q_sat - 1

            ag, bg = FT(5.4 * 1e2), FT(2.55 * 1e5)
            G = FT(1) / (ag + bg / p / rv_sat) / ρ

            av, bv = FT(1.6), FT(124.9)
            F =
                av * (ρ / FT(1e3))^FT(0.525) * rr^FT(0.525) +
                bv * (ρ / FT(1e3))^FT(0.7296) * rr^FT(0.7296)

            return 1 / (1 - q_tot) * S * F * G
        end

        # example values
        T, p = FT(273.15 + 15), FT(90000)
        ϵ = TDI.Rd_over_Rv(tps)
        p_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))
        q_rain_range = range(FT(1e-8), stop = FT(5e-3), length = 10)
        q_tot = FT(15e-3)
        q_vap = FT(0.15) * q_sat
        for q_rai in q_rain_range
            q_lcl = q_tot - q_vap - q_rai
            R = TDI.Rₘ(tps, q_tot, q_lcl + q_rai, FT(0))
            ρ = p / R / T

            tmp1 = CM1.evaporation_sublimation(rain, blk1mvel.rain, aps, tps, q_tot, q_lcl, FT(0), q_rai, FT(0), ρ, T)
            tmp2 = rain_evap_empir(tps, q_rai, q_tot, q_lcl, FT(0), T, p, ρ)

            @assert isapprox(tmp1, tmp2; atol = 1e-6)
        end

        # no condensational growth for rain
        T, p = FT(273.15 + 15), FT(90000)
        ϵ = TDI.Rd_over_Rv(tps)
        p_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))
        q_rai = FT(1e-4)
        q_tot = FT(15e-3)
        q_vap = FT(1.15) * q_sat
        q_ice = FT(0)
        q_liq = q_tot - q_vap - q_ice
        q_sno = FT(0) # TODO - add non-zero case
        R = TDI.Rₘ(tps, q_tot, q_liq, q_ice)
        ρ = p / R / T

        TT.@test CM1.evaporation_sublimation(
            rain,
            blk1mvel.rain,
            aps,
            tps,
            q_tot,
            q_liq,
            q_ice,
            q_rai,
            q_sno,
            ρ,
            T,
        ) ≈ FT(0)

    end

    TT.@testset "SnowSublimation" begin
        # TODO - coudnt find a plot of what it should be from the original paper
        # just checking if the number stays the same

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
            ϵ = TDI.Rd_over_Rv(tps)
            p_sat = TDI.saturation_vapor_pressure_over_ice(tps, T)
            q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))

            for eps in [FT(0.95), FT(1.05)]
                cnt += 1

                q_sno = FT(1e-4)
                q_rai = FT(0) # TODO -add non-zero case
                q_tot = eps * q_sat + q_sno
                q_ice = FT(0)
                q_liq = FT(0)

                R = TDI.Rₘ(tps, q_tot, q_liq + q_rai, q_ice + q_sno)
                ρ = p / R / T

                tmp1 = CM1.evaporation_sublimation(
                    snow,
                    blk1mvel.snow,
                    aps,
                    tps,
                    q_tot,
                    q_liq,
                    q_ice,
                    q_rai,
                    q_sno,
                    ρ,
                    T,
                )
                TT.@test tmp1 ≈ ref_val[cnt] rtol = 1e-2
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

TT.@testset "Microphysics 1M Tests ($FT)" for FT in (Float64, Float32)
    test_microphysics1M(FT)
end
nothing
