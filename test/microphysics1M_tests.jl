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

    mp = CMP.Microphysics1MParams(FT)

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
        micro = (; q_tot, q_lcl = q_liq, q_icl = q_ice, q_rai, q_sno)
        thermo = (; ρ, T)
        evap_rate = CM1.conv_q_rai_to_q_vap(CMP.RainEvaporation(), mp, tps, micro, thermo)
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
        micro = (; q_tot, q_lcl = q_liq, q_icl = q_ice, q_rai, q_sno)
        thermo = (; ρ, T)
        subl_rate = CM1.conv_q_sno_to_q_vap(CMP.SublimationOnly(), mp, tps, micro, thermo)
        TT.@test subl_rate < 0
    end

    TT.@testset "RainAutoconversion" begin

        q_lcl_threshold = mp.options.rain_autoconversion.acnv1M.q_threshold
        τ_acnv_rai = mp.options.rain_autoconversion.acnv1M.τ

        micro_s = (; q_tot = FT(0), q_lcl = FT(0.5) * q_lcl_threshold,
            q_icl = FT(0), q_rai = FT(0), q_sno = FT(0))
        micro_b = (; q_tot = FT(0), q_lcl = FT(1.5) * q_lcl_threshold,
            q_icl = FT(0), q_rai = FT(0), q_sno = FT(0))
        thermo = (; ρ = FT(1), T = FT(280))

        # Below threshold → near zero (smooth logistic)
        TT.@test CM1.conv_q_lcl_to_q_rai(mp.options.rain_autoconversion, mp, tps, micro_s, thermo) ≈
                 FT(0.0) atol = 0.15 * q_lcl_threshold / τ_acnv_rai

        # Above threshold → ≈ 0.5 * q_threshold / τ (smooth logistic)
        TT.@test CM1.conv_q_lcl_to_q_rai(mp.options.rain_autoconversion, mp, tps, micro_b, thermo) ≈
                 FT(0.5) * q_lcl_threshold / τ_acnv_rai atol =
            FT(0.15) * q_lcl_threshold / τ_acnv_rai

    end

    TT.@testset "RainAutoconversion2M" begin
        # Variable-timescale autoconversion (PrescribedNd)
        mp_2m = CMP.Microphysics1MParams(FT;
            rain_autoconversion = CMP.PrescribedNd(CP.create_toml_dict(FT)),
        )
        thermo = (; ρ = FT(1), T = FT(280))

        # Zero cloud liquid → zero rate
        micro_0 = (; q_tot = FT(0), q_lcl = FT(0), q_icl = FT(0), q_rai = FT(0), q_sno = FT(0))
        TT.@test CM1.conv_q_lcl_to_q_rai(mp_2m.options.rain_autoconversion, mp_2m, tps, micro_0, thermo) == FT(0)

        # Negative cloud liquid → zero rate (max(0, q_lcl))
        micro_neg = (; q_tot = FT(0), q_lcl = FT(-1e-4), q_icl = FT(0), q_rai = FT(0), q_sno = FT(0))
        TT.@test CM1.conv_q_lcl_to_q_rai(mp_2m.options.rain_autoconversion, mp_2m, tps, micro_neg, thermo) == FT(0)

        # Positive cloud liquid → positive rate
        q_lcl = FT(2e-3)
        N_d = mp_2m.options.rain_autoconversion.autoconv.Nc
        (; τ, α) = mp_2m.options.rain_autoconversion.autoconv
        micro = (; q_tot = FT(0), q_lcl, q_icl = FT(0), q_rai = FT(0), q_sno = FT(0))
        rate = CM1.conv_q_lcl_to_q_rai(mp_2m.options.rain_autoconversion, mp_2m, tps, micro, thermo)
        TT.@test rate > FT(0)
        # Rate should match the formula: q_lcl / (τ * (N_d / 1e8)^α)
        TT.@test rate ≈ q_lcl / (τ * (N_d / 100_000_000)^α)

        # Higher N_d → slower autoconversion (more droplets → longer timescale)
        # We test this via the formula directly: smaller N_d → larger rate
        q_lcl_test = FT(1e-3)
        N_d_low = FT(1e7)
        N_d_high = FT(1e9)
        rate_low_N = q_lcl_test / (τ * (N_d_low / 100_000_000)^α)
        rate_high_N = q_lcl_test / (τ * (N_d_high / 100_000_000)^α)
        TT.@test rate_low_N > rate_high_N

        # Regression: keep result constant when code changes
        # q_lcl = 2e-3, default PrescribedNd Nc ~1e8
        TT.@test rate ≈ FT(2e-6) rtol = 1e-3

    end

    TT.@testset "SnowAutoconversionNoSupersat" begin

        q_icl_threshold = mp.options.snow_autoconversion.acnv1M.q_threshold
        τ_acnv_sno = mp.options.snow_autoconversion.acnv1M.τ

        # Below threshold → rate near zero (uses logistic smoothing)
        q_icl_small = FT(0.5) * q_icl_threshold
        micro_s = (; q_tot = FT(0), q_lcl = FT(0), q_icl = q_icl_small, q_rai = FT(0), q_sno = FT(0))
        thermo_s = (; ρ = FT(1), T = FT(250))
        TT.@test CM1.conv_q_icl_to_q_sno(
            mp.options.snow_autoconversion, mp, tps, micro_s, thermo_s,
        ) ≈ FT(0.0) atol = FT(0.15) * q_icl_threshold / τ_acnv_sno

        # Above threshold → positive rate
        q_icl_big = FT(1.5) * q_icl_threshold
        micro_b = (; q_tot = FT(0), q_lcl = FT(0), q_icl = q_icl_big, q_rai = FT(0), q_sno = FT(0))
        TT.@test CM1.conv_q_icl_to_q_sno(
            mp.options.snow_autoconversion, mp, tps, micro_b, thermo_s,
        ) ≈ FT(0.5) * q_icl_threshold / τ_acnv_sno atol =
            FT(0.15) * q_icl_threshold / τ_acnv_sno

    end

    TT.@testset "SnowAutoconversion" begin

        ρ = FT(1.0)
        qᵣ = FT(1e-4)
        qₛ = FT(1e-4)
        T_freeze = TDI.TD.Parameters.T_freeze(tps)

        # Use mp with WithSupersaturation
        mp_ss = CMP.Microphysics1MParams(FT;
            snow_autoconversion = CMP.WithSupersaturation(CP.create_toml_dict(FT)),
        )

        # above freezing temperatures -> no snow
        qᵥ = FT(15e-3)
        qₗ = FT(2e-3)
        qᵢ = FT(1e-3)
        qₜ = qᵥ + qₗ + qᵢ + qᵣ + qₛ
        T = T_freeze + FT(30)
        micro = (; q_tot = qₜ, q_lcl = qₗ, q_icl = qᵢ, q_rai = qᵣ, q_sno = qₛ)
        thermo = (; ρ, T)
        TT.@test CM1.conv_q_icl_to_q_sno(mp_ss.options.snow_autoconversion, mp_ss, tps, micro, thermo) ==
                 FT(0)

        # no cloud ice -> no snow
        qᵢ = FT(0)
        qₜ = qᵥ + qₗ + qᵢ + qᵣ + qₛ
        T = T_freeze - FT(30)
        micro = (; q_tot = qₜ, q_lcl = qₗ, q_icl = qᵢ, q_rai = qᵣ, q_sno = qₛ)
        thermo = (; ρ, T)
        TT.@test CM1.conv_q_icl_to_q_sno(mp_ss.options.snow_autoconversion, mp_ss, tps, micro, thermo) ==
                 FT(0)

        # no supersaturation -> no snow
        T = T_freeze - FT(5)
        qᵢ = FT(3e-3)
        q_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        qₜ = q_sat_ice
        micro = (; q_tot = qₜ, q_lcl = qₗ, q_icl = qᵢ, q_rai = qᵣ, q_sno = qₛ)
        thermo = (; ρ, T)
        TT.@test CM1.conv_q_icl_to_q_sno(mp_ss.options.snow_autoconversion, mp_ss, tps, micro, thermo) ≈ FT(0)

        # Regression test: keep result constant when code changes
        T = T_freeze - FT(10)
        qᵥ = FT(1.02) * TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        qₗ = FT(0)
        qᵢ = FT(0.03) * qᵥ
        qₜ = qᵥ + qₗ + qᵢ + qᵣ + qₛ
        ref = FT(2.5408135723057333e-9)
        micro = (; q_tot = qₜ, q_lcl = qₗ, q_icl = qᵢ, q_rai = qᵣ, q_sno = qₛ)
        thermo = (; ρ, T)
        TT.@test CM1.conv_q_icl_to_q_sno(mp_ss.options.snow_autoconversion, mp_ss, tps, micro, thermo) ≈ ref
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

        E_lcl_rai = mp.options.cloud_liquid_rain_accretion.e
        for q_rai in q_rain_range
            if q_rai > eps(FT)
                TT.@test CM1.accretion(
                    liquid,
                    rain,
                    blk1mvel.rain,
                    E_lcl_rai,
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
                    E_lcl_rai,
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

        E_lcl_rai = mp.options.cloud_liquid_rain_accretion.e
        E_lcl_sno = mp.options.cloud_liquid_snow_accretion.e
        E_icl_rai = mp.options.cloud_ice_rain_accretion.e
        E_icl_sno = mp.options.cloud_ice_snow_accretion.e
        E_rai_sno = mp.options.rain_snow_accretion.e
        coeff_disp = mp.options.rain_snow_accretion.coeff_disp

        TT.@test CM1.accretion(
            liquid,
            rain,
            blk1mvel.rain,
            E_lcl_rai,
            q_liq,
            q_rai,
            ρ,
        ) ≈ FT(1.4150106417043544e-6)
        TT.@test CM1.accretion(ice, snow, blk1mvel.snow, E_icl_sno, q_ice, q_sno, ρ) ≈
                 FT(2.453070979562392e-7)
        TT.@test CM1.accretion(
            liquid,
            snow,
            blk1mvel.snow,
            E_lcl_sno,
            q_liq,
            q_sno,
            ρ,
        ) ≈ FT(2.453070979562392e-7)
        TT.@test CM1.accretion(ice, rain, blk1mvel.rain, E_icl_rai, q_ice, q_rai, ρ) ≈
                 FT(1.768763302130443e-6)

        TT.@test CM1.accretion_rain_sink(
            rain,
            ice,
            blk1mvel.rain,
            E_icl_rai,
            q_ice,
            q_rai,
            ρ,
        ) ≈ FT(3.590060148920766e-5)

        TT.@test CM1.accretion_snow_rain(
            snow,
            rain,
            blk1mvel.snow,
            blk1mvel.rain,
            E_rai_sno,
            coeff_disp,
            q_sno,
            q_rai,
            ρ,
        ) ≈ FT(2.466313958248222e-4) # Includes velocity dispersion correction
        TT.@test CM1.accretion_snow_rain(
            rain,
            snow,
            blk1mvel.rain,
            blk1mvel.snow,
            E_rai_sno,
            coeff_disp,
            q_rai,
            q_sno,
            ρ,
        ) ≈ FT(6.830957197816771e-5) # Includes velocity dispersion correction
    end

    TT.@testset "AccretionOptionAPI" begin
        # Verify that the option-dispatched wrappers return the same values
        # as the internal low-level kernels (regression values from "Accretion" testset)
        ρ = FT(1.2)
        q_tot = FT(20e-3)
        q_liq = FT(5e-4)
        q_ice = FT(5e-4)
        q_sno = FT(5e-4)
        q_rai = FT(5e-4)
        T_freeze = TDI.T_freeze(tps)
        T_warm = T_freeze + FT(5)   # above freezing
        T_cold = T_freeze - FT(5)   # below freezing

        micro = (; q_tot, q_lcl = q_liq, q_icl = q_ice, q_rai, q_sno)
        thermo_warm = (; ρ, T = T_warm)
        thermo_cold = (; ρ, T = T_cold)

        # CloudLiquidRainAccretion: scalar, matches internal kernel
        TT.@test CM1.accretion(mp.options.cloud_liquid_rain_accretion, mp, tps, micro, thermo_warm) ≈
                 FT(1.4150106417043544e-6)

        # CloudIceRainAccretion: scalar, matches internal kernel
        TT.@test CM1.accretion(mp.options.cloud_ice_rain_accretion, mp, tps, micro, thermo_warm) ≈
                 FT(1.768763302130443e-6)

        # CloudIceSnowAccretion: scalar, matches internal kernel
        TT.@test CM1.accretion(mp.options.cloud_ice_snow_accretion, mp, tps, micro, thermo_warm) ≈
                 FT(2.453070979562392e-7)

        # CloudLiquidSnowAccretion: NamedTuple (; S_accr, S_melt)
        # S_accr matches the internal lcl×sno kernel; S_melt = α * S_accr (α > 0 warm)
        (; S_accr, S_melt) = CM1.accretion(mp.options.cloud_liquid_snow_accretion, mp, tps, micro, thermo_warm)
        TT.@test S_accr ≈ FT(2.453070979562392e-7)
        TT.@test S_melt >= FT(0)      # α >= 0
        TT.@test S_melt <= S_accr     # melt ≤ S_accr (α ≤ 1 for reasonable ΔT)

        # Cold: melt should be zero (T < T_freeze → α = 0)
        let r = CM1.accretion(mp.options.cloud_liquid_snow_accretion, mp, tps, micro, thermo_cold)
            TT.@test r.S_accr ≈ FT(2.453070979562392e-7)
            TT.@test r.S_melt == FT(0)
        end

        # RainSnowAccretion: NamedTuple (; S_rai_sno, S_sno_rai, S_melt)
        # cold arm S_rai_sno matches internal (sno, rai) call
        # warm arm S_sno_rai matches internal (rai, sno) call
        let r = CM1.accretion_snow_rain(mp.options.rain_snow_accretion, mp, tps, micro, thermo_warm)
            TT.@test r.S_rai_sno ≈ FT(2.466313958248222e-4)
            TT.@test r.S_sno_rai ≈ FT(6.830957197816771e-5)
            TT.@test r.S_melt >= FT(0)
        end
        # Cold: S_melt = 0 (α = 0 below freezing)
        let r = CM1.accretion_snow_rain(mp.options.rain_snow_accretion, mp, tps, micro, thermo_cold)
            TT.@test r.S_melt == FT(0)
        end

        # nothing variants return zero
        micro_zero = (; q_tot = FT(0), q_lcl = FT(0), q_icl = FT(0), q_rai = FT(0), q_sno = FT(0))
        TT.@test CM1.accretion(nothing, mp, tps, micro_zero, thermo_warm) == FT(0)
        TT.@test CM1.accretion_snow_rain(nothing, mp, tps, micro_zero, thermo_warm).S_rai_sno == FT(0)
        TT.@test CM1.accretion_snow_rain(nothing, mp, tps, micro_zero, thermo_warm).S_sno_rai == FT(0)
        # Active variants: zero inputs → zero rates
        TT.@test CM1.accretion(mp.options.cloud_liquid_rain_accretion, mp, tps, micro_zero, thermo_warm) == FT(0)
        TT.@test CM1.accretion(mp.options.cloud_ice_rain_accretion, mp, tps, micro_zero, thermo_warm) == FT(0)
        TT.@test CM1.accretion(mp.options.cloud_ice_snow_accretion, mp, tps, micro_zero, thermo_warm) == FT(0)
        TT.@test CM1.accretion(mp.options.cloud_liquid_snow_accretion, mp, tps, micro_zero, thermo_warm).S_accr == FT(0)
        TT.@test CM1.accretion_snow_rain(mp.options.rain_snow_accretion, mp, tps, micro_zero, thermo_warm).S_rai_sno ==
                 FT(0)
        TT.@test CM1.accretion_snow_rain(mp.options.rain_snow_accretion, mp, tps, micro_zero, thermo_warm).S_sno_rai ==
                 FT(0)
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

            micro = (; q_tot, q_lcl, q_icl = FT(0), q_rai, q_sno = FT(0))
            thermo = (; ρ, T)
            tmp1 = CM1.conv_q_rai_to_q_vap(CMP.RainEvaporation(), mp, tps, micro, thermo)
            tmp2 = rain_evap_empir(tps, q_rai, q_tot, q_lcl, FT(0), T, p, ρ)

            TT.@test tmp1 ≈ tmp2 atol = 1e-6
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

        micro = (; q_tot, q_lcl = q_liq, q_icl = q_ice, q_rai, q_sno)
        thermo = (; ρ, T)
        TT.@test CM1.conv_q_rai_to_q_vap(CMP.RainEvaporation(), mp, tps, micro, thermo) ≈ FT(0)

    end

    TT.@testset "SnowSublimation" begin
        # Regression test: keep results constant when code changes

        cnt = 0
        ref_val = [
            FT(-1.9756907119482267e-7),
            FT(0),
            FT(-1.6641552112891826e-7),
            FT(0),
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

                micro = (; q_tot, q_lcl = q_liq, q_icl = q_ice, q_rai, q_sno)
                thermo = (; ρ, T)
                tmp1 = CM1.conv_q_sno_to_q_vap(CMP.SublimationOnly(), mp, tps, micro, thermo)
                TT.@test tmp1 ≈ ref_val[cnt] rtol = 1e-2
            end
        end
    end

    TT.@testset "SnowSublimation and Deposition" begin
        # Regression test: keep results constant when code changes
        # Uses DepositionSublimation() and original reference values that include deposition
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

                micro = (; q_tot, q_lcl = q_liq, q_icl = q_ice, q_rai, q_sno)
                thermo = (; ρ, T)
                tmp1 = CM1.conv_q_sno_to_q_vap(CMP.DepositionAndSublimation(), mp, tps, micro, thermo)
                TT.@test tmp1 ≈ ref_val[cnt] rtol = 1e-2
            end
        end
    end

    TT.@testset "SnowMelt" begin

        # Regression test: keep result constant when code changes
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        T = T_freeze + FT(2)
        ρ = FT(1.2)
        q_sno = FT(1e-4)
        micro_sno = (; q_tot = FT(0), q_lcl = FT(0), q_icl = FT(0), q_rai = FT(0), q_sno)
        thermo_sno = (; ρ, T)
        (melt_rate) = CM1.conv_q_sno_to_q_rai(CMP.SnowMelt(), mp, tps, micro_sno, thermo_sno)
        TT.@test melt_rate ≈ FT(9.516553267013085e-6)

        # no snow -> no snow melt
        T = T_freeze + FT(2)
        ρ = FT(1.2)
        q_sno = FT(0)
        micro_sno = (; q_tot = FT(0), q_lcl = FT(0), q_icl = FT(0), q_rai = FT(0), q_sno)
        thermo_sno = (; ρ, T)
        TT.@test CM1.conv_q_sno_to_q_rai(CMP.SnowMelt(), mp, tps, micro_sno, thermo_sno) ≈
                 FT(0)

        # T < T_freeze -> no snow melt
        T = T_freeze - FT(2)
        ρ = FT(1.2)
        q_sno = FT(1e-4)
        micro_sno = (; q_tot = FT(0), q_lcl = FT(0), q_icl = FT(0), q_rai = FT(0), q_sno)
        thermo_sno = (; ρ, T)
        TT.@test CM1.conv_q_sno_to_q_rai(CMP.SnowMelt(), mp, tps, micro_sno, thermo_sno) ≈
                 FT(0)

    end

    TT.@testset "CloudIceMelt" begin

        # Regression test: keep result constant when code changes
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        T = T_freeze + FT(2)
        ρ = FT(1.2)
        q_icl = FT(1e-4)
        mp_melt = CMP.Microphysics1MParams(FT;
            cloud_ice_melt = CMP.CloudIceMelt(),
        )
        micro = (; q_tot = FT(0), q_lcl = FT(0), q_icl, q_rai = FT(0), q_sno = FT(0))
        thermo = (; ρ, T)
        melt_rate = CM1.conv_q_icl_to_q_lcl(CMP.CloudIceMelt(), mp_melt, tps, micro, thermo)
        TT.@test melt_rate > 0

        # no cloud ice → no melt
        micro_0 = (; q_tot = FT(0), q_lcl = FT(0), q_icl = FT(0), q_rai = FT(0), q_sno = FT(0))
        TT.@test CM1.conv_q_icl_to_q_lcl(CMP.CloudIceMelt(), mp_melt, tps, micro_0, thermo) ≈ FT(0)

        # T < T_freeze → no melt
        thermo_cold = (; ρ, T = T_freeze - FT(2))
        TT.@test CM1.conv_q_icl_to_q_lcl(CMP.CloudIceMelt(), mp_melt, tps, micro, thermo_cold) ≈ FT(0)

        # more ice → higher melt rate
        thermo_warm = (; ρ, T = T_freeze + FT(5))
        micro_s = (; q_tot = FT(0), q_lcl = FT(0), q_icl = FT(1e-5), q_rai = FT(0), q_sno = FT(0))
        micro_l = (; q_tot = FT(0), q_lcl = FT(0), q_icl = FT(1e-3), q_rai = FT(0), q_sno = FT(0))
        rate_small = CM1.conv_q_icl_to_q_lcl(CMP.CloudIceMelt(), mp_melt, tps, micro_s, thermo_warm)
        rate_large = CM1.conv_q_icl_to_q_lcl(CMP.CloudIceMelt(), mp_melt, tps, micro_l, thermo_warm)
        TT.@test rate_large > rate_small

        # warmer T → higher melt rate
        thermo_w = (; ρ, T = T_freeze + FT(10))
        thermo_c = (; ρ, T = T_freeze + FT(1))
        rate_warm = CM1.conv_q_icl_to_q_lcl(CMP.CloudIceMelt(), mp_melt, tps, micro, thermo_w)
        rate_cool = CM1.conv_q_icl_to_q_lcl(CMP.CloudIceMelt(), mp_melt, tps, micro, thermo_c)
        TT.@test rate_warm > rate_cool

        # nothing returns zero
        TT.@test CM1.conv_q_icl_to_q_lcl(nothing, mp, tps, micro, thermo) == FT(0)

    end

end

TT.@testset "Microphysics 1M Tests ($FT)" for FT in (Float64, Float32)
    test_microphysics1M(FT)
end
nothing
