import Plots

import Thermodynamics
import CloudMicrophysics
import ClimaParams
import CloudMicrophysics.Common as CO

const PL = Plots
const CM1 = CloudMicrophysics.Microphysics1M
const TD = Thermodynamics
const CMP = CloudMicrophysics.Parameters

FT = Float64

const tps = Thermodynamics.Parameters.ThermodynamicsParameters(FT)
const aps = CMP.AirProperties(FT)
const liquid = CMP.CloudLiquid(FT)
const ice = CMP.CloudIce(FT)
const rain = CMP.Rain(FT)
const snow = CMP.Snow(FT)
const Chen2022 = CMP.Chen2022VelType(FT)
const Blk1MVel = CMP.Blk1MVelType(FT)
const ce = CMP.CollisionEff(FT)

# eq. 5b in [Grabowski1996](@cite)
function accretion_empirical(q_rai::DT, q_liq::DT, q_tot::DT) where {DT <: Real}
    rr = q_rai / (DT(1) - q_tot)
    rl = q_liq / (DT(1) - q_tot)
    return DT(2.2) * rl * rr^DT(7 / 8)
end

# eq. 5c in [Grabowski1996](@cite)
function rain_evap_empirical(
    q_rai::DT,
    q::TD.PhasePartition,
    T::DT,
    p::DT,
    ρ::DT,
) where {DT <: Real}

    ts_neq = TD.PhaseNonEquil_ρTq(tps, ρ, T, q)
    q_sat = TD.q_vap_saturation(tps, ts_neq)
    q_vap = q.tot - q.liq
    rr = q_rai / (DT(1) - q.tot)
    rv_sat = q_sat / (DT(1) - q.tot)
    S = q_vap / q_sat - DT(1)

    ag, bg = 5.4 * 1e2, 2.55 * 1e5
    G = DT(1) / (ag + bg / p / rv_sat) / ρ

    av, bv = 1.6, 124.9
    F =
        av * (ρ / DT(1e3))^DT(0.525) * rr^DT(0.525) +
        bv * (ρ / DT(1e3))^DT(0.7296) * rr^DT(0.7296)

    return DT(1) / (DT(1) - q.tot) * S * F * G
end

# example values
q_liq_range = range(1e-8, stop = 5e-3, length = 100)
q_ice_range = range(1e-8, stop = 5e-3, length = 100)
q_rain_range = range(1e-8, stop = 5e-3, length = 100)
q_snow_range = range(1e-8, stop = 5e-3, length = 100)
ρ_air, ρ_air_ground = 1.2, 1.22
q_liq, q_ice, q_tot = 5e-4, 5e-4, 20e-3

T = 273.15
PL.plot(
    q_liq_range * 1e3,
    [CM1.conv_q_liq_to_q_rai(rain.acnv1M, q_liq) for q_liq in q_liq_range],
    linewidth = 3,
    xlabel = "q_liq or q_ice [g/kg]",
    ylabel = "autoconversion rate [1/s]",
    label = "Rain",
)
PL.plot!(
    q_ice_range * 1e3,
    [
        CM1.conv_q_ice_to_q_sno(
            ice,
            aps,
            tps,
            TD.PhasePartition(q_tot, 0.0, q_ice),
            ρ_air,
            T - 5,
        ) for q_ice in q_ice_range
    ],
    linewidth = 3,
    label = "Snow T=-5C",
)
PL.plot!(
    q_ice_range * 1e3,
    [
        CM1.conv_q_ice_to_q_sno(
            ice,
            aps,
            tps,
            TD.PhasePartition(q_tot, 0.0, q_ice),
            ρ_air,
            T - 10,
        ) for q_ice in q_ice_range
    ],
    linewidth = 3,
    label = "Snow T=-10C",
)
PL.plot!(
    q_ice_range * 1e3,
    [
        CM1.conv_q_ice_to_q_sno(
            ice,
            aps,
            tps,
            TD.PhasePartition(q_tot, 0.0, q_ice),
            ρ_air,
            T - 15,
        ) for q_ice in q_ice_range
    ],
    linewidth = 3,
    label = "Snow T=-15C",
)
PL.savefig("autoconversion_rate.svg") # hide

PL.plot(
    q_rain_range * 1e3,
    [
        CM1.accretion(liquid, rain, Blk1MVel.rain, ce, q_liq, q_rai, ρ_air) for
        q_rai in q_rain_range
    ],
    linewidth = 3,
    xlabel = "q_rain or q_snow [g/kg]",
    ylabel = "accretion rate [1/s]",
    label = "Liq+Rain-CLIMA",
)
PL.plot!(
    q_rain_range * 1e3,
    [
        CM1.accretion(ice, rain, Blk1MVel.rain, ce, q_ice, q_rai, ρ_air) for
        q_rai in q_rain_range
    ],
    linewidth = 3,
    label = "Ice+Rain-CLIMA",
)
PL.plot!(
    q_snow_range * 1e3,
    [
        CM1.accretion(liquid, snow, Blk1MVel.snow, ce, q_liq, q_sno, ρ_air) for
        q_sno in q_snow_range
    ],
    linewidth = 3,
    label = "Liq+Snow-CLIMA",
)
PL.plot!(
    q_snow_range * 1e3,
    [
        CM1.accretion(ice, snow, Blk1MVel.snow, ce, q_ice, q_sno, ρ_air) for
        q_sno in q_snow_range
    ],
    linewidth = 4,
    linestyle = :dash,
    label = "Ice+Snow-CLIMA",
)
PL.plot!(
    q_rain_range * 1e3,
    [accretion_empirical(q_rai, q_liq, q_tot) for q_rai in q_rain_range],
    linewidth = 3,
    label = "Liq+Rain-Empirical",
)
PL.savefig("accretion_rate.svg") # hide

q_ice = 1e-6
PL.plot(
    q_rain_range * 1e3,
    [
        CM1.accretion_rain_sink(
            rain,
            ice,
            Blk1MVel.rain,
            ce,
            q_ice,
            q_rai,
            ρ_air,
        ) for q_rai in q_rain_range
    ],
    linewidth = 3,
    xlabel = "q_rain or q_snow [g/kg]",
    ylabel = "accretion rain sink rate [1/s]",
    label = "q_ice=1e-6",
)
q_ice = 1e-5
PL.plot!(
    q_rain_range * 1e3,
    [
        CM1.accretion_rain_sink(
            rain,
            ice,
            Blk1MVel.rain,
            ce,
            q_ice,
            q_rai,
            ρ_air,
        ) for q_rai in q_rain_range
    ],
    linewidth = 3,
    xlabel = "q_rain or q_snow [g/kg]",
    ylabel = "accretion rain sink rate [1/s]",
    label = "q_ice=1e-5",
)
q_ice = 1e-4
PL.plot!(
    q_rain_range * 1e3,
    [
        CM1.accretion_rain_sink(
            rain,
            ice,
            Blk1MVel.rain,
            ce,
            q_ice,
            q_rai,
            ρ_air,
        ) for q_rai in q_rain_range
    ],
    linewidth = 3,
    xlabel = "q_rain or q_snow [g/kg]",
    ylabel = "accretion rain sink rate [1/s]",
    label = "q_ice=1e-4",
)
PL.savefig("accretion_rain_sink_rate.svg") # hide

q_sno = 1e-6
PL.plot(
    q_rain_range * 1e3,
    [
        CM1.accretion_snow_rain(
            rain,
            snow,
            Blk1MVel.rain,
            Blk1MVel.snow,
            ce,
            q_rai,
            q_sno,
            ρ_air,
        ) for q_rai in q_rain_range
    ],
    linewidth = 3,
    xlabel = "q_rain [g/kg]",
    ylabel = "snow-rain accretion rate [1/s] T>0",
    label = "q_snow=1e-6",
)
q_sno = 1e-5
PL.plot!(
    q_rain_range * 1e3,
    [
        CM1.accretion_snow_rain(
            rain,
            snow,
            Blk1MVel.rain,
            Blk1MVel.snow,
            ce,
            q_rai,
            q_sno,
            ρ_air,
        ) for q_rai in q_rain_range
    ],
    linewidth = 3,
    label = "q_snow=1e-5",
)
q_sno = 1e-4
PL.plot!(
    q_rain_range * 1e3,
    [
        CM1.accretion_snow_rain(
            rain,
            snow,
            Blk1MVel.rain,
            Blk1MVel.snow,
            ce,
            q_rai,
            q_sno,
            ρ_air,
        ) for q_rai in q_rain_range
    ],
    linewidth = 3,
    label = "q_snow=1e-4",
)
PL.savefig("accretion_snow_rain_above_freeze.svg") # hide

q_rai = 1e-6
PL.plot(
    q_snow_range * 1e3,
    [
        CM1.accretion_snow_rain(
            snow,
            rain,
            Blk1MVel.snow,
            Blk1MVel.rain,
            ce,
            q_sno,
            q_rai,
            ρ_air,
        ) for q_sno in q_snow_range
    ],
    linewidth = 3,
    xlabel = "q_snow [g/kg]",
    ylabel = "snow-rain accretion rate [1/s] T<0",
    label = "q_rain=1e-6",
)
q_rai = 1e-5
PL.plot!(
    q_snow_range * 1e3,
    [
        CM1.accretion_snow_rain(
            snow,
            rain,
            Blk1MVel.snow,
            Blk1MVel.rain,
            ce,
            q_sno,
            q_rai,
            ρ_air,
        ) for q_sno in q_snow_range
    ],
    linewidth = 3,
    label = "q_rain=1e-5",
)
q_rai = 1e-4
PL.plot!(
    q_snow_range * 1e3,
    [
        CM1.accretion_snow_rain(
            snow,
            rain,
            Blk1MVel.snow,
            Blk1MVel.rain,
            ce,
            q_sno,
            q_rai,
            ρ_air,
        ) for q_sno in q_snow_range
    ],
    linewidth = 3,
    label = "q_rain=1e-4",
)
PL.savefig("accretion_snow_rain_below_freeze.svg") # hide

# example values
T, p = 273.15 + 15, 90000.0
ϵ = 1.0 / TD.Parameters.molmass_ratio(tps)
p_sat = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
q_rain_range = range(1e-8, stop = 5e-3, length = 100)
q_tot = 15e-3
q_vap = 0.15 * q_sat
q_ice = 0.0
q_liq = q_tot - q_vap - q_ice
q = TD.PhasePartition(q_tot, q_liq, q_ice)
R = TD.gas_constant_air(tps, q)
ρ = p / R / T

PL.plot(
    q_rain_range * 1e3,
    [
        CM1.evaporation_sublimation(
            rain,
            Blk1MVel.rain,
            aps,
            tps,
            q,
            q_rai,
            ρ,
            T,
        ) for q_rai in q_rain_range
    ],
    xlabel = "q_rain [g/kg]",
    linewidth = 3,
    ylabel = "rain evaporation rate [1/s]",
    label = "ClimateMachine",
)
PL.plot!(
    q_rain_range * 1e3,
    [rain_evap_empirical(q_rai, q, T, p, ρ) for q_rai in q_rain_range],
    linewidth = 3,
    label = "empirical",
)
PL.savefig("rain_evaporation_rate.svg") # hide

# example values
T, p = 273.15 - 15, 90000.0
ϵ = 1.0 / TD.Parameters.molmass_ratio(tps)
p_sat = TD.saturation_vapor_pressure(tps, T, TD.Ice())
q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
q_snow_range = range(1e-8, stop = 5e-3, length = 100)
q_tot = 15e-3
q_vap = 0.15 * q_sat
q_liq = 0.0
q_ice = q_tot - q_vap - q_liq
q = TD.PhasePartition(q_tot, q_liq, q_ice)
R = TD.gas_constant_air(tps, q)
ρ = p / R / T

PL.plot(
    q_snow_range * 1e3,
    [
        CM1.evaporation_sublimation(
            snow,
            Blk1MVel.snow,
            aps,
            tps,
            q,
            q_sno,
            ρ,
            T,
        ) for q_sno in q_snow_range
    ],
    xlabel = "q_snow [g/kg]",
    linewidth = 3,
    ylabel = "snow deposition sublimation rate [1/s]",
    label = "T<0",
)

T, p = 273.15 + 15, 90000.0
ϵ = 1.0 / TD.Parameters.molmass_ratio(tps)
p_sat = TD.saturation_vapor_pressure(tps, T, TD.Ice())
q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
q_snow_range = range(1e-8, stop = 5e-3, length = 100)
q_tot = 15e-3
q_vap = 0.15 * q_sat
q_liq = 0.0
q_ice = q_tot - q_vap - q_liq
q = TD.PhasePartition(q_tot, q_liq, q_ice)
R = TD.gas_constant_air(tps, q)
ρ = p / R / T

PL.plot!(
    q_snow_range * 1e3,
    [
        CM1.evaporation_sublimation(
            snow,
            Blk1MVel.snow,
            aps,
            tps,
            q,
            q_sno,
            ρ,
            T,
        ) for q_sno in q_snow_range
    ],
    xlabel = "q_snow [g/kg]",
    linewidth = 3,
    ylabel = "snow deposition sublimation rate [1/s]",
    label = "T>0",
)
PL.savefig("snow_sublimation_deposition_rate.svg") # hide

T = 273.15
PL.plot(
    q_snow_range * 1e3,
    [
        CM1.snow_melt(snow, Blk1MVel.snow, aps, tps, q_sno, ρ, T + 2) for
        q_sno in q_snow_range
    ],
    xlabel = "q_snow [g/kg]",
    linewidth = 3,
    ylabel = "snow melt rate [1/s]",
    label = "T=2C",
)
PL.plot!(
    q_snow_range * 1e3,
    [
        CM1.snow_melt(snow, Blk1MVel.snow, aps, tps, q_sno, ρ, T + 4) for
        q_sno in q_snow_range
    ],
    xlabel = "q_snow [g/kg]",
    linewidth = 3,
    label = "T=4C",
)
PL.plot!(
    q_snow_range * 1e3,
    [
        CM1.snow_melt(snow, Blk1MVel.snow, aps, tps, q_sno, ρ, T + 6) for
        q_sno in q_snow_range
    ],
    xlabel = "q_snow [g/kg]",
    linewidth = 3,
    label = "T=6C",
)
PL.savefig("snow_melt_rate.svg") # hide
