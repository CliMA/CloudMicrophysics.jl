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

        #action
        vt_rai = CM1.terminal_velocity(rain, Ch2022.rain, ρ, q_rai)
        v_bigger = CM1.terminal_velocity(rain, Ch2022.rain, ρ, q_rai * 2)

        #test
        TT.@test vt_rai ≈ 5.25213637238494 rtol = 1e-5
        TT.@test CM1.terminal_velocity(rain, Ch2022.rain, ρ, FT(0))[1] ≈ 0 atol =
            eps(FT)
        TT.@test v_bigger > vt_rai
    end

    TT.@testset "1M_microphysics - Chen 2022 snow terminal velocity" begin
        #setup
        ρ = FT(1.1)
        q_sno = FT(5e-4)

        #action
        vt_sno = CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, q_sno)
        v_bigger = CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, q_sno * 2)

        #test
        TT.@test vt_sno ≈ 0.8573952434834717 rtol = 3e-6
        TT.@test CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, FT(0)) ≈ 0 atol =
            eps(FT)
        TT.@test v_bigger > vt_sno
    end

    TT.@testset "1M_microphysics - Chen 2022 snow terminal velocity with snow_shape" begin
        #setup
        ρ = FT(1.1)
        q_sno = FT(5e-4)

        # action - test both Oblate and Prolate shapes
        vt_oblate = CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, q_sno, CM1.Oblate())
        vt_prolate = CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, q_sno, CM1.Prolate())
        v_bigger_oblate = CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, q_sno * 2, CM1.Oblate())
        v_bigger_prolate = CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, q_sno * 2, CM1.Prolate())

        # test - regression values
        TT.@test vt_oblate > 0
        TT.@test vt_prolate > 0
        # zero input returns zero
        TT.@test CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, FT(0), CM1.Oblate()) ≈ 0 atol = eps(FT)
        TT.@test CM1.terminal_velocity(snow, Ch2022.large_ice, ρ, FT(0), CM1.Prolate()) ≈ 0 atol = eps(FT)
        # monotonicity: more snow -> higher velocity
        TT.@test v_bigger_oblate > vt_oblate
        TT.@test v_bigger_prolate > vt_prolate
        # both shapes give reasonable velocities (within an order of magnitude of each other)
        TT.@test 0.1 < vt_oblate / vt_prolate < 10
    end

    TT.@testset "1M_microphysics - 1M snow terminal velocity" begin
        # NaN check for edge case near zero
        TT.@test !isnan(CM1.terminal_velocity(snow, blk1mvel.snow, FT(0.2439843), FT(3.0f-45)))

        # zero input returns zero
        TT.@test CM1.terminal_velocity(snow, blk1mvel.snow, FT(1.2), FT(0)) ≈ 0 atol = eps(FT)

        # monotonicity: more snow -> higher velocity
        ρ = FT(1.2)
        q_sno = FT(5e-4)
        vt_sno = CM1.terminal_velocity(snow, blk1mvel.snow, ρ, q_sno)
        v_bigger = CM1.terminal_velocity(snow, blk1mvel.snow, ρ, q_sno * 2)
        TT.@test vt_sno > 0
        TT.@test v_bigger > vt_sno
    end

    TT.@testset "lambda_inverse" begin
        # Direct tests for the Marshall-Palmer rate parameter
        ρ = FT(1.2)

        # zero specific content returns minimum value (r0 * 1e-5)
        λ_inv_zero_rain = CM1.lambda_inverse(rain.pdf, rain.mass, FT(0), ρ)
        λ_inv_zero_snow = CM1.lambda_inverse(snow.pdf, snow.mass, FT(0), ρ)
        TT.@test λ_inv_zero_rain ≈ rain.mass.r0 * FT(1e-5)
        TT.@test λ_inv_zero_snow ≈ snow.mass.r0 * FT(1e-5)

        # λ⁻¹ increases with specific content (larger mean particle size)
        q_small = FT(1e-5)
        q_large = FT(1e-3)
        λ_inv_small = CM1.lambda_inverse(rain.pdf, rain.mass, q_small, ρ)
        λ_inv_large = CM1.lambda_inverse(rain.pdf, rain.mass, q_large, ρ)
        TT.@test λ_inv_large > λ_inv_small

        # λ⁻¹ for snow also increases with content
        λ_inv_snow_small = CM1.lambda_inverse(snow.pdf, snow.mass, q_small, ρ)
        λ_inv_snow_large = CM1.lambda_inverse(snow.pdf, snow.mass, q_large, ρ)
        TT.@test λ_inv_snow_large > λ_inv_snow_small

        # reasonable order of magnitude: λ⁻¹ ~ 100 μm to 1 mm for typical rain
        q_rai = FT(1e-4)
        λ_inv = CM1.lambda_inverse(rain.pdf, rain.mass, q_rai, ρ)
        TT.@test FT(1e-5) < λ_inv < FT(1e-2)  # 10 μm to 10 mm
    end

    TT.@testset "MixedPhaseEvaporation" begin
        # Test evaporation with non-zero snow content (previously untested)
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        T, p = T_freeze + FT(10), FT(90000)
        ϵ = TDI.Rd_over_Rv(tps)
        p_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))

        q_rai = FT(1e-4)
        q_sno = FT(1e-4)  # Non-zero snow in warm conditions
        q_tot = FT(15e-3)
        q_vap = FT(0.7) * q_sat  # subsaturated
        q_ice = FT(0)
        q_liq = q_tot - q_vap - q_ice - q_rai - q_sno
        R = TDI.Rₘ(tps, q_tot, q_liq + q_rai, q_ice + q_sno)
        ρ = p / R / T

        # Rain evaporates, rate should be negative
        evap_rate = CM1.evaporation_sublimation(
            rain, blk1mvel.rain, aps, tps,
            q_tot, q_liq, q_ice, q_rai, q_sno, ρ, T,
        )
        TT.@test evap_rate < 0
    end

    TT.@testset "MixedPhaseSublimation" begin
        # Test sublimation with non-zero rain content (previously untested)
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        T, p = T_freeze - FT(10), FT(90000)
        ϵ = TDI.Rd_over_Rv(tps)
        p_sat = TDI.saturation_vapor_pressure_over_ice(tps, T)
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))

        q_sno = FT(1e-4)
        q_rai = FT(1e-4)  # Non-zero rain in cold conditions
        q_vap = FT(0.9) * q_sat  # slightly subsaturated over ice
        q_tot = q_vap + q_sno + q_rai
        q_ice = FT(0)
        q_liq = FT(0)
        R = TDI.Rₘ(tps, q_tot, q_liq + q_rai, q_ice + q_sno)
        ρ = p / R / T

        # Snow sublimation rate should be negative (losing mass)
        subl_rate = CM1.evaporation_sublimation(
            snow, blk1mvel.snow, aps, tps,
            q_tot, q_liq, q_ice, q_rai, q_sno, ρ, T,
        )
        TT.@test subl_rate < 0
    end

    TT.@testset "RainAutoconversion" begin

        q_lcl_threshold = rain.acnv1M.q_threshold
        τ_acnv_rai = rain.acnv1M.τ

        q_lcl_small = FT(0.5) * q_lcl_threshold
        TT.@test CM1.conv_q_lcl_to_q_rai(rain.acnv1M, q_lcl_small) == FT(0)

        TT.@test CM1.conv_q_lcl_to_q_rai(rain.acnv1M, q_lcl_small, true) ≈
                 FT(0.0) atol = 0.15 * q_lcl_threshold / τ_acnv_rai

        q_lcl_big = FT(1.5) * q_lcl_threshold
        TT.@test CM1.conv_q_lcl_to_q_rai(rain.acnv1M, q_lcl_big) ≈
                 FT(0.5) * q_lcl_threshold / τ_acnv_rai

        TT.@test CM1.conv_q_lcl_to_q_rai(rain.acnv1M, q_lcl_big, true) ≈
                 FT(0.5) * q_lcl_threshold / τ_acnv_rai atol =
            FT(0.15) * q_lcl_threshold / τ_acnv_rai

    end

    TT.@testset "SnowAutoconversionNoSupersat" begin

        q_icl_threshold = snow.acnv1M.q_threshold
        τ_acnv_sno = snow.acnv1M.τ

        q_icl_small = FT(0.5) * q_icl_threshold
        TT.@test CM1.conv_q_icl_to_q_sno_no_supersat(
            snow.acnv1M,
            q_icl_small,
        ) == FT(0)

        TT.@test CM1.conv_q_icl_to_q_sno_no_supersat(
            snow.acnv1M,
            q_icl_small,
            true,
        ) ≈ FT(0.0) atol = FT(0.15) * q_icl_threshold / τ_acnv_sno

        q_icl_big = FT(1.5) * q_icl_threshold
        TT.@test CM1.conv_q_icl_to_q_sno_no_supersat(snow.acnv1M, q_icl_big) ≈
                 FT(0.5) * q_icl_threshold / τ_acnv_sno

        TT.@test CM1.conv_q_icl_to_q_sno_no_supersat(
            snow.acnv1M,
            q_icl_big,
            true,
        ) ≈ FT(0.5) * q_icl_threshold / τ_acnv_sno atol =
            FT(0.15) * q_icl_threshold / τ_acnv_sno

        TT.@test CM1.conv_q_icl_to_q_sno_no_supersat(snow.acnv1M, q_icl_big) ==
                 CM1.conv_q_icl_to_q_sno_no_supersat(
            snow.acnv1M,
            q_icl_big,
            false,
        )

    end

    TT.@testset "SnowAutoconversion" begin

        ρ = FT(1.0)
        qᵣ = FT(1e-4)
        qₛ = FT(1e-4)
        T_freeze = TDI.TD.Parameters.T_freeze(tps)

        # above freezing temperatures -> no snow
        qᵥ = FT(15e-3)
        qₗ = FT(2e-3)
        qᵢ = FT(1e-3)
        qₜ = qᵥ + qₗ + qᵢ + qᵣ + qₛ
        T = T_freeze + FT(30)
        TT.@test CM1.conv_q_icl_to_q_sno(ice, aps, tps, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T) == FT(0)

        # no cloud ice -> no snow
        qᵢ = FT(0)
        qₜ = qᵥ + qₗ + qᵢ + qᵣ + qₛ
        T = T_freeze - FT(30)
        TT.@test CM1.conv_q_icl_to_q_sno(ice, aps, tps, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T) == FT(0)

        # no supersaturation -> no snow
        T = T_freeze - FT(5)
        qᵢ = FT(3e-3)
        q_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        qₜ = q_sat_ice
        TT.@test CM1.conv_q_icl_to_q_sno(ice, aps, tps, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T) ≈ FT(0)

        # Regression test: keep result constant when code changes
        T = T_freeze - FT(10)
        qᵥ = FT(1.02) * TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        qₗ = FT(0)
        qᵢ = FT(0.03) * qᵥ
        qₜ = qᵥ + qₗ + qᵢ + qᵣ + qₛ
        ref = FT(2.5408135723057333e-9)
        TT.@test CM1.conv_q_icl_to_q_sno(ice, aps, tps, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T) ≈ ref
    end

    TT.@testset "RainLiquidAccretion" begin
        # Compare against eq. 5b in [Grabowski1996](@cite)
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
        # Regression test: keep results constant when code changes

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
        ) ≈ FT(2.466313958248222e-4) # Includes velocity dispersion correction
        TT.@test CM1.accretion_snow_rain(
            rain,
            snow,
            blk1mvel.rain,
            blk1mvel.snow,
            ce,
            q_rai,
            q_sno,
            ρ,
        ) ≈ FT(6.830957197816771e-5) # Includes velocity dispersion correction
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
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        T, p = T_freeze + FT(15), FT(90000)
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
            dtmp1 = CM1.∂evaporation_sublimation_∂q_precip(
                rain,
                blk1mvel.rain,
                aps,
                tps,
                q_tot,
                q_lcl,
                FT(0),
                q_rai,
                FT(0),
                ρ,
                T,
            )
            tmp2 = rain_evap_empir(tps, q_rai, q_tot, q_lcl, FT(0), T, p, ρ)

            TT.@test tmp1 ≈ tmp2 atol = 1e-6
            TT.@test dtmp1 ≈ tmp1 / q_rai
        end

        # no condensational growth for rain
        T, p = T_freeze + FT(15), FT(90000)
        ϵ = TDI.Rd_over_Rv(tps)
        p_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))
        q_rai = FT(1e-4)
        q_tot = FT(15e-3)
        q_vap = FT(1.15) * q_sat
        q_ice = FT(0)
        q_liq = q_tot - q_vap - q_ice
        q_sno = FT(0)
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
        # Regression test: keep results constant when code changes

        cnt = 0
        ref_val = [
            FT(-1.9756907119482267e-7),
            FT(1.9751292385808357e-7),
            FT(-1.6641552112891826e-7),
            FT(1.663814937710236e-7),
        ]
        # some example values
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        for T in [T_freeze + FT(2), T_freeze - FT(2)]
            p = FT(90000)
            ϵ = TDI.Rd_over_Rv(tps)
            p_sat = TDI.saturation_vapor_pressure_over_ice(tps, T)
            q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))

            for eps in [FT(0.95), FT(1.05)]
                cnt += 1

                q_sno = FT(1e-4)
                q_rai = FT(0)
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
                dtmp1 = CM1.∂evaporation_sublimation_∂q_precip(
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
                TT.@test dtmp1 ≈ tmp1 / q_sno
            end
        end
    end

    TT.@testset "SnowMelt" begin

        # Regression test: keep result constant when code changes
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        T = T_freeze + FT(2)
        ρ = FT(1.2)
        q_sno = FT(1e-4)
        (melt_rate) = CM1.snow_melt(snow, blk1mvel.snow, aps, tps, q_sno, ρ, T)
        melt_deriv = CM1.∂snow_melt_∂q_sno(snow, blk1mvel.snow, aps, tps, q_sno, ρ, T)
        TT.@test melt_rate ≈ FT(9.516553267013085e-6)
        TT.@test melt_deriv ≈ melt_rate / q_sno

        # no snow -> no snow melt
        T = T_freeze + FT(2)
        ρ = FT(1.2)
        q_sno = FT(0)
        TT.@test CM1.snow_melt(snow, blk1mvel.snow, aps, tps, q_sno, ρ, T) ≈
                 FT(0)

        # T < T_freeze -> no snow melt
        T = T_freeze - FT(2)
        ρ = FT(1.2)
        q_sno = FT(1e-4)
        TT.@test CM1.snow_melt(snow, blk1mvel.snow, aps, tps, q_sno, ρ, T) ≈
                 FT(0)

    end

    TT.@testset "RainEvaporationDerivative_FiniteDiff" begin
        # Verify ∂(evap)/∂q_rai via central finite differences at constant e_tot, ρ, q_tot.
        # The analytical derivative uses rate/q_rai, which neglects the n0(q)‐dependence
        # in the Marshall–Palmer distribution and temperature feedback through latent heat.
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        T, p = T_freeze + FT(15), FT(90000)
        ϵ = TDI.Rd_over_Rv(tps)
        p_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))
        q_rai = FT(1e-3)
        q_tot = FT(15e-3)
        q_vap = FT(0.15) * q_sat  # subsaturated
        q_lcl = q_tot - q_vap - q_rai
        q_icl = FT(0)
        q_sno = FT(0)
        R = TDI.Rₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)
        ρ = p / R / T

        rate = CM1.evaporation_sublimation(
            rain, blk1mvel.rain, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T,
        )
        drate = CM1.∂evaporation_sublimation_∂q_precip(
            rain, blk1mvel.rain, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T,
        )
        TT.@test rate < 0  # evaporation
        TT.@test drate ≈ rate / q_rai

        # Central finite difference at constant e_tot, ρ, q_tot:
        # perturb q_rai, recover T from energy conservation
        q_liq = q_lcl + q_rai
        q_ice = q_icl + q_sno
        e_int = TDI.TD.internal_energy(tps, T, q_tot, q_liq, q_ice)

        Δq = FT(1e-8)
        q_liq_p = q_lcl + (q_rai + Δq)
        q_liq_m = q_lcl + (q_rai - Δq)
        T_p = TDI.TD.air_temperature(tps, e_int, q_tot, q_liq_p, q_ice)
        T_m = TDI.TD.air_temperature(tps, e_int, q_tot, q_liq_m, q_ice)
        rate_p = CM1.evaporation_sublimation(
            rain, blk1mvel.rain, aps, tps, q_tot, q_lcl, q_icl, q_rai + Δq, q_sno, ρ, T_p,
        )
        rate_m = CM1.evaporation_sublimation(
            rain, blk1mvel.rain, aps, tps, q_tot, q_lcl, q_icl, q_rai - Δq, q_sno, ρ, T_m,
        )
        fd_deriv = (rate_p - rate_m) / (2Δq)
        TT.@test sign(drate) == sign(fd_deriv)
        TT.@test drate ≈ fd_deriv rtol = FT(0.2)
    end

    TT.@testset "SnowSublimationDerivative_FiniteDiff" begin
        # Verify ∂(subl)/∂q_sno via central finite differences at constant e_tot, ρ, q_tot.
        # The analytical derivative uses rate/q_sno, which neglects the n0(q_sno)‐dependence
        # in the Marshall–Palmer distribution (n0 ~ q_sno^ν for snow) and temperature
        # feedback through latent heat. This is a cruder approximation for snow than
        # for rain because snow n0 depends on q_sno.
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        T, p = T_freeze - FT(10), FT(90000)
        ϵ = TDI.Rd_over_Rv(tps)
        p_sat = TDI.saturation_vapor_pressure_over_ice(tps, T)
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))
        q_sno = FT(1e-4)
        q_tot = FT(0.95) * q_sat + q_sno
        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        R = TDI.Rₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)
        ρ = p / R / T

        rate = CM1.evaporation_sublimation(
            snow, blk1mvel.snow, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T,
        )
        drate = CM1.∂evaporation_sublimation_∂q_precip(
            snow, blk1mvel.snow, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T,
        )
        TT.@test rate < 0  # sublimation
        TT.@test drate ≈ rate / q_sno

        # Central finite difference at constant e_tot, ρ, q_tot:
        # perturb q_sno, recover T from energy conservation
        q_liq = q_lcl + q_rai
        q_ice = q_icl + q_sno
        e_int = TDI.TD.internal_energy(tps, T, q_tot, q_liq, q_ice)

        Δq = FT(1e-8)
        q_ice_p = q_icl + (q_sno + Δq)
        q_ice_m = q_icl + (q_sno - Δq)
        T_p = TDI.TD.air_temperature(tps, e_int, q_tot, q_liq, q_ice_p)
        T_m = TDI.TD.air_temperature(tps, e_int, q_tot, q_liq, q_ice_m)
        rate_p = CM1.evaporation_sublimation(
            snow, blk1mvel.snow, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno + Δq, ρ, T_p,
        )
        rate_m = CM1.evaporation_sublimation(
            snow, blk1mvel.snow, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno - Δq, ρ, T_m,
        )
        fd_deriv = (rate_p - rate_m) / (2Δq)
        TT.@test sign(drate) == sign(fd_deriv)
        TT.@test drate ≈ fd_deriv rtol = FT(0.7)
    end

    TT.@testset "SnowMeltDerivative_FiniteDiff" begin
        # Verify ∂(melt)/∂q_sno via central finite differences at constant e_tot, ρ, q_tot.
        # The analytical derivative uses rate/q_sno.  Snow melt has a weaker n0(q_sno)
        # dependence than sublimation, so the approximation is more accurate (~3% error).
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        T = T_freeze + FT(4)
        ρ = FT(1.2)
        q_sno = FT(1e-4)
        q_lcl = FT(0)
        q_icl = FT(0)
        q_rai = FT(0)
        q_tot = FT(10e-3)

        rate = CM1.snow_melt(snow, blk1mvel.snow, aps, tps, q_sno, ρ, T)
        drate = CM1.∂snow_melt_∂q_sno(snow, blk1mvel.snow, aps, tps, q_sno, ρ, T)
        TT.@test rate > 0  # melting
        TT.@test drate ≈ rate / q_sno

        # Central finite difference at constant e_tot, ρ, q_tot:
        # perturb q_sno, recover T from energy conservation
        q_liq = q_lcl + q_rai
        q_ice = q_icl + q_sno
        e_int = TDI.TD.internal_energy(tps, T, q_tot, q_liq, q_ice)

        Δq = FT(1e-8)
        q_ice_p = q_icl + (q_sno + Δq)
        q_ice_m = q_icl + (q_sno - Δq)
        T_p = TDI.TD.air_temperature(tps, e_int, q_tot, q_liq, q_ice_p)
        T_m = TDI.TD.air_temperature(tps, e_int, q_tot, q_liq, q_ice_m)
        rate_p = CM1.snow_melt(snow, blk1mvel.snow, aps, tps, q_sno + Δq, ρ, T_p)
        rate_m = CM1.snow_melt(snow, blk1mvel.snow, aps, tps, q_sno - Δq, ρ, T_m)
        fd_deriv = (rate_p - rate_m) / (2Δq)
        TT.@test sign(drate) == sign(fd_deriv)
        TT.@test drate ≈ fd_deriv rtol = FT(0.05)
    end
end

TT.@testset "Microphysics 1M Tests ($FT)" for FT in (Float64, Float32)
    test_microphysics1M(FT)
end
nothing
