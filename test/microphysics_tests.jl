import Test

import Thermodynamics
import CloudMicrophysics
import CLIMAParameters
const CP = CLIMAParameters

const TT = Test
const TD = Thermodynamics
const CMT = CloudMicrophysics.CommonTypes
const CMNe = CloudMicrophysics.MicrophysicsNonEq
const CM0 = CloudMicrophysics.Microphysics0M
const CM1 = CloudMicrophysics.Microphysics1M
const CM2 = CloudMicrophysics.Microphysics2M

include(joinpath(pkgdir(CloudMicrophysics), "test", "create_parameters.jl"))
FT = Float64
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
const prs = cloud_microphysics_parameters(toml_dict)

const liquid = CMT.LiquidType()
const ice = CMT.IceType()
const rain = CMT.RainType()
const snow = CMT.SnowType()

TT.@testset "τ_relax" begin

    TT.@test CMNe.τ_relax(prs, liquid) ≈ 10
    TT.@test CMNe.τ_relax(prs, ice) ≈ 10

end

TT.@testset "0M_microphysics" begin

    _τ_precip = CMP.τ_precip(prs)
    _qc_0 = CMP.qc_0(prs)
    _S_0 = CMP.S_0(prs)

    q_vap_sat = 10e-3
    qc = 3e-3
    q_tot = 13e-3
    frac = [0.0, 0.5, 1.0]

    # no rain if no cloud
    q = TD.PhasePartition(q_tot, 0.0, 0.0)
    TT.@test CM0.remove_precipitation(prs, q) ≈ 0
    TT.@test CM0.remove_precipitation(prs, q, q_vap_sat) ≈ 0

    # rain based on qc threshold
    for lf in frac
        q_liq = qc * lf
        q_ice = (1 - lf) * qc

        q = TD.PhasePartition(q_tot, q_liq, q_ice)

        TT.@test CM0.remove_precipitation(prs, q) ≈
                 -max(0, q_liq + q_ice - _qc_0) / _τ_precip
    end

    # rain based on supersaturation threshold
    for lf in frac
        q_liq = qc * lf
        q_ice = (1 - lf) * qc

        q = TD.PhasePartition(q_tot, q_liq, q_ice)

        TT.@test CM0.remove_precipitation(prs, q, q_vap_sat) ≈
                 -max(0, q_liq + q_ice - _S_0 * q_vap_sat) / _τ_precip
    end
end

TT.@testset "RainFallSpeed" begin
    # eq. 5d in [Grabowski1996](@cite)
    function terminal_velocity_empir(
        q_rai::FT,
        q_tot::FT,
        ρ::FT,
        ρ_air_ground::FT,
    ) where {FT <: Real}
        rr = q_rai / (1 - q_tot)
        vel = FT(14.34) * ρ_air_ground^FT(0.5) * ρ^-FT(0.3654) * rr^FT(0.1346)
        return vel
    end

    # some example values
    q_rain_range = range(1e-8, stop = 5e-3, length = 10)
    ρ_air, q_tot, ρ_air_ground = 1.2, 20 * 1e-3, 1.22

    for q_rai in q_rain_range
        TT.@test CM1.terminal_velocity(prs, rain, ρ_air, q_rai) ≈
                 terminal_velocity_empir(q_rai, q_tot, ρ_air, ρ_air_ground) atol =
            0.2 * terminal_velocity_empir(q_rai, q_tot, ρ_air, ρ_air_ground)

    end
end

TT.@testset "CloudLiquidCondEvap" begin

    q_liq_sat = 5e-3
    frac = [0.0, 0.5, 1.0, 1.5]

    _τ_cond_evap = CMNe.τ_relax(prs, liquid)

    for fr in frac
        q_liq = q_liq_sat * fr

        TT.@test CMNe.conv_q_vap_to_q_liq_ice(
            prs,
            liquid,
            TD.PhasePartition(0.0, q_liq_sat, 0.0),
            TD.PhasePartition(0.0, q_liq, 0.0),
        ) ≈ (1 - fr) * q_liq_sat / _τ_cond_evap
    end
end

TT.@testset "CloudIceCondEvap" begin

    q_ice_sat = 2e-3
    frac = [0.0, 0.5, 1.0, 1.5]

    _τ_cond_evap = CMNe.τ_relax(prs, ice)

    for fr in frac
        q_ice = q_ice_sat * fr

        TT.@test CMNe.conv_q_vap_to_q_liq_ice(
            prs,
            ice,
            TD.PhasePartition(0.0, 0.0, q_ice_sat),
            TD.PhasePartition(0.0, 0.0, q_ice),
        ) ≈ (1 - fr) * q_ice_sat / _τ_cond_evap
    end
end

TT.@testset "RainAutoconversion" begin

    _q_liq_threshold = CMP.q_liq_threshold(prs)
    _τ_acnv_rai = CMP.τ_acnv_rai(prs)

    q_liq_small = 0.5 * _q_liq_threshold
    TT.@test CM1.conv_q_liq_to_q_rai(prs, q_liq_small) == 0.0

    TT.@test CM1.conv_q_liq_to_q_rai(
        prs,
        q_liq_small,
        smooth_transition = true,
    ) ≈ 0.0 atol = 0.15 * _q_liq_threshold / _τ_acnv_rai

    q_liq_big = 1.5 * _q_liq_threshold
    TT.@test CM1.conv_q_liq_to_q_rai(prs, q_liq_big) ==
             0.5 * _q_liq_threshold / _τ_acnv_rai

    TT.@test CM1.conv_q_liq_to_q_rai(prs, q_liq_big, smooth_transition = true) ≈
             0.5 * _q_liq_threshold / _τ_acnv_rai atol =
        0.15 * _q_liq_threshold / _τ_acnv_rai

end

TT.@testset "SnowAutoconversionNoSupersat" begin

    _q_ice_threshold = CMP.q_ice_threshold(prs)
    _τ_acnv_sno = CMP.τ_acnv_sno(prs)

    q_ice_small = 0.5 * _q_ice_threshold
    TT.@test CM1.conv_q_ice_to_q_sno_no_supersat(prs, q_ice_small) == 0.0

    TT.@test CM1.conv_q_ice_to_q_sno_no_supersat(
        prs,
        q_ice_small,
        smooth_transition = true,
    ) ≈ 0.0 atol = 0.15 * _q_ice_threshold / _τ_acnv_sno

    q_ice_big = 1.5 * _q_ice_threshold
    TT.@test CM1.conv_q_ice_to_q_sno_no_supersat(prs, q_ice_big) ≈
             0.5 * _q_ice_threshold / _τ_acnv_sno

    TT.@test CM1.conv_q_ice_to_q_sno_no_supersat(
        prs,
        q_ice_big,
        smooth_transition = true,
    ) ≈ 0.5 * _q_ice_threshold / _τ_acnv_sno atol =
        0.15 * _q_ice_threshold / _τ_acnv_sno
end

TT.@testset "SnowAutoconversion" begin

    thermo_params = CMP.thermodynamics_params(prs)
    ρ = 1.0

    # above freezing temperatures -> no snow
    q = TD.PhasePartition(15e-3, 2e-3, 1e-3)
    T = 273.15 + 30
    TT.@test CM1.conv_q_ice_to_q_sno(prs, q, ρ, T) == 0.0

    # no ice -> no snow
    q = TD.PhasePartition(15e-3, 2e-3, 0.0)
    T = 273.15 - 30
    TT.@test CM1.conv_q_ice_to_q_sno(prs, q, ρ, T) == 0.0

    # no supersaturation -> no snow
    T = 273.15 - 5
    q_sat_ice = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
    q = TD.PhasePartition(q_sat_ice, 2e-3, 3e-3)
    TT.@test CM1.conv_q_ice_to_q_sno(prs, q, ρ, T) == 0.0

    # TODO - coudnt find a plot of what it should be from the original paper
    # just chacking if the number stays the same
    T = 273.15 - 10
    q_vap = 1.02 * TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
    q_liq = 0.0
    q_ice = 0.03 * q_vap
    q = TD.PhasePartition(q_vap + q_liq + q_ice, q_liq, q_ice)
    TT.@test CM1.conv_q_ice_to_q_sno(prs, q, ρ, T) ≈ 1.8512022335645584e-9

end

TT.@testset "RainLiquidAccretion" begin

    # eq. 5b in [Grabowski1996](@cite)
    function accretion_empir(q_rai::FT, q_liq::FT, q_tot::FT) where {FT <: Real}
        rr = q_rai / (FT(1) - q_tot)
        rl = q_liq / (FT(1) - q_tot)
        return FT(2.2) * rl * rr^FT(7 / 8)
    end

    # some example values
    q_rain_range = range(1e-8, stop = 5e-3, length = 10)
    ρ_air, q_liq, q_tot = 1.2, 5e-4, 20e-3

    for q_rai in q_rain_range
        TT.@test CM1.accretion(prs, liquid, rain, q_liq, q_rai, ρ_air) ≈
                 accretion_empir(q_rai, q_liq, q_tot) atol =
            (0.1 * accretion_empir(q_rai, q_liq, q_tot))
    end
end

TT.@testset "Accretion" begin
    # TODO - coudnt find a plot of what it should be from the original paper
    # just chacking if the number stays the same

    # some example values
    ρ = 1.2
    q_tot = 20e-3
    q_ice = 5e-4
    q_sno = 5e-4
    q_liq = 5e-4
    q_rai = 5e-4

    TT.@test CM1.accretion(prs, liquid, rain, q_liq, q_rai, ρ) ≈
             1.4150106417043544e-6
    TT.@test CM1.accretion(prs, ice, snow, q_ice, q_sno, ρ) ≈
             2.453070979562392e-7
    TT.@test CM1.accretion(prs, liquid, snow, q_liq, q_sno, ρ) ≈
             2.453070979562392e-7
    TT.@test CM1.accretion(prs, ice, rain, q_ice, q_rai, ρ) ≈
             1.768763302130443e-6

    TT.@test CM1.accretion_rain_sink(prs, q_ice, q_rai, ρ) ≈
             3.085229094251214e-5

    TT.@test CM1.accretion_snow_rain(prs, snow, rain, q_sno, q_rai, ρ) ≈
             2.1705865794293408e-4
    TT.@test CM1.accretion_snow_rain(prs, rain, snow, q_rai, q_sno, ρ) ≈
             6.0118801860768854e-5
end

TT.@testset "RainEvaporation" begin

    # eq. 5c in [Grabowski1996](@cite)
    function rain_evap_empir(
        prs::CMP.AbstractCloudMicrophysicsParameters,
        q_rai::FT,
        q::TD.PhasePartition,
        T::FT,
        p::FT,
        ρ::FT,
    ) where {FT <: Real}

        thermo_params = CMP.thermodynamics_params(prs)
        q_sat = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
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

    thermo_params = CMP.thermodynamics_params(prs)
    # example values
    T, p = 273.15 + 15, 90000.0
    ϵ = 1.0 / CMP.molmass_ratio(prs)
    p_sat = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
    q_rain_range = range(1e-8, stop = 5e-3, length = 10)
    q_tot = 15e-3
    q_vap = 0.15 * q_sat
    q_ice = 0.0
    q_liq = q_tot - q_vap - q_ice
    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    R = TD.gas_constant_air(thermo_params, q)
    ρ = p / R / T

    for q_rai in q_rain_range
        TT.@test CM1.evaporation_sublimation(prs, rain, q, q_rai, ρ, T) ≈
                 rain_evap_empir(prs, q_rai, q, T, p, ρ) atol =
            -0.5 * rain_evap_empir(prs, q_rai, q, T, p, ρ)
    end

    # no condensational growth for rain
    T, p = 273.15 + 15, 90000.0
    ϵ = 1.0 / CMP.molmass_ratio(prs)
    p_sat = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
    q_rai = 1e-4
    q_tot = 15e-3
    q_vap = 1.15 * q_sat
    q_ice = 0.0
    q_liq = q_tot - q_vap - q_ice
    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    R = TD.gas_constant_air(thermo_params, q)
    ρ = p / R / T

    TT.@test CM1.evaporation_sublimation(prs, rain, q, q_rai, ρ, T) ≈ 0.0

end

TT.@testset "SnowSublimation" begin
    # TODO - coudnt find a plot of what it should be from the original paper
    # just chacking if the number stays the same

    thermo_params = CMP.thermodynamics_params(prs)
    cnt = 0
    ref_val = [
        -1.9756907119482267e-7,
        1.9751292385808357e-7,
        -1.6641552112891826e-7,
        1.663814937710236e-7,
    ]
    # some example values
    for T in [273.15 + 2, 273.15 - 2]
        p = 90000.0
        ϵ = 1.0 / CMP.molmass_ratio(prs)
        p_sat = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))

        for eps in [0.95, 1.05]
            cnt += 1

            q_tot = eps * q_sat
            q_ice = 0.0
            q_liq = 0.0
            q = TD.PhasePartition(q_tot, q_liq, q_ice)

            q_sno = 1e-4

            R = TD.gas_constant_air(thermo_params, q)
            ρ = p / R / T

            TT.@test CM1.evaporation_sublimation(prs, snow, q, q_sno, ρ, T) ≈
                     ref_val[cnt]

        end
    end
end

TT.@testset "SnowMelt" begin

    # TODO - find a good reference to compare with
    T = 273.15 + 2
    ρ = 1.2
    q_sno = 1e-4
    TT.@test CM1.snow_melt(prs, q_sno, ρ, T) ≈ 9.518235437405256e-6

    # no snow -> no snow melt
    T = 273.15 + 2
    ρ = 1.2
    q_sno = 0.0
    TT.@test CM1.snow_melt(prs, q_sno, ρ, T) ≈ 0

    # T < T_freeze -> no snow melt
    T = 273.15 - 2
    ρ = 1.2
    q_sno = 1e-4
    TT.@test CM1.snow_melt(prs, q_sno, ρ, T) ≈ 0

end

TT.@testset "2M_microphysics - unit tests" begin

    ρ = 1.0

    # no reference data available - checking if callable and not NaN
    q_liq = 0.5e-3
    q_rai = 1e-6
    TT.@test CM2.accretion_KK2000(prs, q_liq, q_rai, ρ) != NaN
    TT.@test CM2.accretion_B1994(prs, q_liq, q_rai, ρ) != NaN
    TT.@test CM2.accretion_TC1980(prs, q_liq, q_rai) != NaN

    # output should be zero if either q_liq or q_rai are zero
    q_liq = 0.0
    q_rai = 1e-6
    TT.@test CM2.conv_q_liq_to_q_rai_KK2000(prs, q_liq, ρ) == 0.0
    TT.@test CM2.conv_q_liq_to_q_rai_B1994(prs, q_liq, ρ) == 0.0
    TT.@test CM2.conv_q_liq_to_q_rai_TC1980(prs, q_liq, ρ) == 0.0
    TT.@test CM2.conv_q_liq_to_q_rai_LD2004(prs, q_liq, ρ) == 0.0
    TT.@test CM2.accretion_KK2000(prs, q_liq, q_rai, ρ) == 0.0
    TT.@test CM2.accretion_B1994(prs, q_liq, q_rai, ρ) == 0.0
    TT.@test CM2.accretion_TC1980(prs, q_liq, q_rai) == 0.0
    q_liq = 0.5e-3
    q_rai = 0.0
    TT.@test CM2.accretion_KK2000(prs, q_liq, q_rai, ρ) == 0.0
    TT.@test CM2.accretion_B1994(prs, q_liq, q_rai, ρ) == 0.0
    TT.@test CM2.accretion_TC1980(prs, q_liq, q_rai) == 0.0

    # far from threshold points, autoconversion with and without smooth transition should 
    # be approximately equal
    q_liq = 0.5e-3
    TT.@test CM2.conv_q_liq_to_q_rai_B1994(
        prs,
        q_liq,
        ρ,
        smooth_transition = true,
    ) ≈ CM2.conv_q_liq_to_q_rai_B1994(
        prs,
        q_liq,
        ρ,
        smooth_transition = false,
    ) rtol = 0.2
    TT.@test CM2.conv_q_liq_to_q_rai_TC1980(
        prs,
        q_liq,
        ρ,
        smooth_transition = true,
    ) ≈ CM2.conv_q_liq_to_q_rai_TC1980(
        prs,
        q_liq,
        ρ,
        smooth_transition = false,
    ) rtol = 0.2
    TT.@test CM2.conv_q_liq_to_q_rai_LD2004(
        prs,
        q_liq,
        ρ,
        smooth_transition = true,
    ) ≈ CM2.conv_q_liq_to_q_rai_LD2004(
        prs,
        q_liq,
        ρ,
        smooth_transition = false,
    ) rtol = 0.2

end

TT.@testset "2M_microphysics - compare with Wood_2005" begin

    ρ = 1.0
    q_liq = 0.5e-3

    # compare with Wood 2005 Fig 1 panel a
    function compare(f, input, output; eps = 0.1)
        TT.@test f(prs, input * 1e-3, ρ) ≈ output atol = eps * output
    end
    compare(
        CM2.conv_q_liq_to_q_rai_KK2000,
        0.03138461538461537,
        2.636846054348105e-12,
    )
    compare(
        CM2.conv_q_liq_to_q_rai_KK2000,
        0.8738461538461537,
        9.491665962977648e-9,
    )
    compare(
        CM2.conv_q_liq_to_q_rai_B1994,
        0.13999999999999999,
        4.584323122458155e-12,
        eps = 1,
    )
    compare(
        CM2.conv_q_liq_to_q_rai_B1994,
        0.9000000000000006,
        5.4940586176564715e-8,
        eps = 1,
    )
    compare(
        CM2.conv_q_liq_to_q_rai_TC1980,
        0.2700000000000001,
        3.2768635256661366e-8,
    )
    compare(
        CM2.conv_q_liq_to_q_rai_TC1980,
        0.9000000000000006,
        5.340418612468997e-7,
    )
    compare(
        CM2.conv_q_liq_to_q_rai_LD2004,
        0.3700000000000002,
        8.697439193234471e-9,
    )
    compare(
        CM2.conv_q_liq_to_q_rai_LD2004,
        0.9000000000000006,
        1.1325570516983242e-7,
    )

    # compare with Wood 2005 Fig 1 panel b
    function compare_Nd(f, input, output; eps = 0.1)
        TT.@test f(prs, q_liq, N_d = input * 1e6, ρ) ≈ output atol =
            eps * output
    end
    compare_Nd(
        CM2.conv_q_liq_to_q_rai_KK2000,
        16.13564081404141,
        6.457285532394289e-8,
    )
    compare_Nd(
        CM2.conv_q_liq_to_q_rai_KK2000,
        652.093931356625,
        8.604011482409198e-11,
    )
    compare_Nd(
        CM2.conv_q_liq_to_q_rai_B1994,
        14.47851799831075,
        4.2829062386778675e-7,
    )
    compare_Nd(
        CM2.conv_q_liq_to_q_rai_B1994,
        693.0425211336465,
        6.076294746898778e-12,
    )
    compare_Nd(
        CM2.conv_q_liq_to_q_rai_TC1980,
        13.658073017575544,
        2.7110779872658386e-7,
    )
    compare_Nd(
        CM2.conv_q_liq_to_q_rai_TC1980,
        205.0970632305975,
        1.0928660431622176e-7,
    )
    compare_Nd(
        CM2.conv_q_liq_to_q_rai_LD2004,
        15.122629721719655,
        1.1647783461546477e-7,
    )
    compare_Nd(
        CM2.conv_q_liq_to_q_rai_LD2004,
        149.01220754857331,
        1.3917890403908125e-8,
        eps = 1,
    )

end
