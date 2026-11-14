import CairoMakie as MK

import CloudMicrophysics
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.ThermodynamicsInterface as TDI

FT = Float64

tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
aps = CMP.AirProperties(FT)
liquid = CMP.CloudLiquid(FT)
ice = CMP.CloudIce(FT)
rain = CMP.Rain(FT)
snow = CMP.Snow(FT)
Chen2022 = CMP.Chen2022VelType(FT)
Blk1MVel = CMP.Blk1MVelType(FT)
ce = CMP.CollisionEff(FT)

# eq. 5b in [Grabowski1996](@cite)
function accretion_empirical(q_rai::DT, q_liq::DT, q_tot::DT) where {DT <: Real}
    rr = q_rai / (DT(1) - q_tot)
    rl = q_liq / (DT(1) - q_tot)
    return DT(2.2) * rl * rr^DT(7 / 8)
end

# eq. 5c in [Grabowski1996](@cite)
function rain_evap_empirical(tps, q_rai, q_tot, q_liq, T, p, ρ)
    DT = eltype(q_rai)

    p_v_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)
    q_sat = TDI.p2q(tps, T, ρ, p_v_sat)

    q_vap = q_tot - q_liq
    rr = q_rai / (DT(1) - q_tot)
    rv_sat = q_sat / (DT(1) - q_tot)
    S = q_vap / q_sat - 1

    ag, bg = 5.4 * 1e2, 2.55 * 1e5
    G = 1 / (ag + bg / p / rv_sat) / ρ

    av, bv = 1.6, 124.9
    F =
        av * (ρ / DT(1e3))^DT(0.525) * rr^DT(0.525) +
        bv * (ρ / DT(1e3))^DT(0.7296) * rr^DT(0.7296)

    return 1 / (1 - q_tot) * S * F * G
end

# example values
q_min, q_max = 1e-8, 5e-3
q_liq_range = range(q_min, stop = q_max, length = 100)
q_ice_range = range(q_min, stop = q_max, length = 100)
q_rain_range = range(q_min, stop = q_max, length = 100)
q_snow_range = range(q_min, stop = q_max, length = 100)
ρ_air, ρ_air_ground = 1.2, 1.22
q_liq, q_ice, q_tot = 5e-4, 5e-4, 20e-3
q_rai = 1e-3
q_sno = 1e-4
limits = (0, q_max * 1e3, 0, nothing)

MK.set_theme!(MK.theme_minimal())

# autoconversion rate figure
T = 273.15
q_range = TDI.PP.(q_tot, 0.0, q_ice_range)
fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_liq or q_ice [g/kg]", ylabel = "autoconversion rate [1/s]", limits)
q_ice_to_q_sno_rate = T -> CM1.conv_q_ice_to_q_sno.(ice, aps, tps, q_range, q_rai, q_sno, ρ_air, T)
MK.lines!(q_liq_range * 1e3, CM1.conv_q_liq_to_q_rai.(rain.acnv1M, q_liq_range), label = "Rain")
MK.lines!(q_ice_range * 1e3, q_ice_to_q_sno_rate(T - 5), label = "Snow T = −5°C")
MK.lines!(q_ice_range * 1e3, q_ice_to_q_sno_rate(T - 10), label = "Snow T = −10°C")
MK.lines!(q_ice_range * 1e3, q_ice_to_q_sno_rate(T - 15), label = "Snow T = −15°C")
MK.axislegend(ax; position = :lt)
MK.save("autoconversion_rate.svg", fig) # hide

# accretion rate figure
fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_rain or q_snow [g/kg]", ylabel = "accretion rate [1/s]", limits)
MK.lines!(
    q_rain_range * 1e3, CM1.accretion.(liquid, rain, Blk1MVel.rain, ce, q_liq, q_rain_range, ρ_air),
    label = "Liq+Rain-CliMA",
)
MK.lines!(
    q_rain_range * 1e3, CM1.accretion.(ice, rain, Blk1MVel.rain, ce, q_ice, q_rain_range, ρ_air),
    label = "Ice+Rain-CliMA",
)
MK.lines!(
    q_snow_range * 1e3, CM1.accretion.(liquid, snow, Blk1MVel.snow, ce, q_liq, q_snow_range, ρ_air),
    label = "Liq+Snow-CliMA",
)
MK.lines!(
    q_snow_range * 1e3, CM1.accretion.(ice, snow, Blk1MVel.snow, ce, q_ice, q_snow_range, ρ_air),
    label = "Ice+Snow-CliMA", linewidth = 4, linestyle = :dash,
)
MK.lines!(
    q_rain_range * 1e3, accretion_empirical.(q_rain_range, q_liq, q_tot),
    label = "Liq+Rain-Empirical",
)
MK.axislegend(ax; position = :lt)
MK.save("accretion_rate.svg", fig) # hide

# accretion rain sink rate figure
fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_rain or q_snow [g/kg]", ylabel = "accretion rain sink rate [1/s]", limits)
_accr_rain_sink(q_ice) = CM1.accretion_rain_sink.(rain, ice, Blk1MVel.rain, ce, q_ice, q_rain_range, ρ_air) # hide
MK.lines!(q_rain_range * 1e3, _accr_rain_sink(1e-6), label = "q_ice = 1e-6")
MK.lines!(q_rain_range * 1e3, _accr_rain_sink(1e-5), label = "q_ice = 1e-5")
MK.lines!(q_rain_range * 1e3, _accr_rain_sink(1e-4), label = "q_ice = 1e-4")
MK.axislegend(ax; position = :lt)
MK.save("accretion_rain_sink_rate.svg", fig) # hide

# accretion snow-rain rate figure
fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_rain [g/kg]", ylabel = "snow-rain accretion rate [1/s] T>0", limits)
_accr_snow_rain(q_sno) =
    CM1.accretion_snow_rain.(rain, snow, Blk1MVel.rain, Blk1MVel.snow, ce, q_rain_range, q_sno, ρ_air)
MK.lines!(q_rain_range * 1e3, _accr_snow_rain(1e-6), label = "q_snow = 1e-6")
MK.lines!(q_rain_range * 1e3, _accr_snow_rain(1e-5), label = "q_snow = 1e-5")
MK.lines!(q_rain_range * 1e3, _accr_snow_rain(1e-4), label = "q_snow = 1e-4")
MK.axislegend(ax; position = :lt)
MK.save("accretion_snow_rain_above_freeze.svg", fig) # hide

# accretion snow-rain rate figure
fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_snow [g/kg]", ylabel = "snow-rain accretion rate [1/s] T<0", limits)
_accr_snow_rain2(q_sno) =
    CM1.accretion_snow_rain.(snow, rain, Blk1MVel.snow, Blk1MVel.rain, ce, q_snow_range, q_sno, ρ_air)
MK.lines!(q_snow_range * 1e3, _accr_snow_rain2(1e-6), label = "q_rain = 1e-6")
MK.lines!(q_snow_range * 1e3, _accr_snow_rain2(1e-5), label = "q_rain = 1e-5")
MK.lines!(q_snow_range * 1e3, _accr_snow_rain2(1e-4), label = "q_rain = 1e-4")
MK.axislegend(ax; position = :lt)
MK.save("accretion_snow_rain_below_freeze.svg", fig) # hide

# example values
T, p = 273.15 + 15, 90000.0
ϵ = TDI.Rd_over_Rv(tps)
p_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)
q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
q_rain_range = range(1e-8, stop = 5e-3, length = 100)
q_tot = 15e-3
q_vap = 0.15 * q_sat
q_ice = 0.0
q_liq = q_tot - q_vap - q_ice
q_sno = 0.0
q = TDI.PP(q_tot, q_liq, q_ice)
R = TDI.Rₘ(tps, q_tot, q_liq, q_ice)
ρ = p / R / T

fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_rain [g/kg]", ylabel = "rain evaporation rate [1/s]")
MK.xlims!(ax, 0, q_max * 1e3)
MK.lines!(
    q_rain_range * 1e3, CM1.evaporation_sublimation.(rain, Blk1MVel.rain, aps, tps, Ref(q), q_rain_range, q_sno, ρ, T),
    label = "ClimateMachine",
)
MK.lines!(q_rain_range * 1e3, rain_evap_empirical.(tps, q_rain_range, q_tot, q_liq, T, p, ρ), label = "empirical")
MK.axislegend(ax; position = :rt)
MK.save("rain_evaporation_rate.svg", fig) # hide

# snow sublimation rate figure
fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_snow [g/kg]", ylabel = "snow deposition sublimation rate [1/s]")
MK.xlims!(ax, 0, q_max * 1e3)

# example values 1
T, p = 273.15 - 15, 90000.0
ϵ = TDI.Rd_over_Rv(tps)
p_sat = TDI.saturation_vapor_pressure_over_ice(tps, T)
q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
q_snow_range = range(1e-8, stop = 5e-3, length = 100)
q_tot = 15e-3
q_vap = 0.15 * q_sat
q_liq = 0.0
q_ice = q_tot - q_vap - q_liq
q_rai = 0.0
q = TDI.PP(q_tot, q_liq, q_ice)
R = TDI.Rₘ(tps, q_tot, q_liq, q_ice)
ρ = p / R / T
rate = CM1.evaporation_sublimation.(snow, Blk1MVel.snow, aps, tps, Ref(q), q_rai, q_snow_range, ρ, T)
MK.lines!(q_snow_range * 1e3, rate, label = "T < 0°C")

# example values 2
T, p = 273.15 + 15, 90000.0
ϵ = TDI.Rd_over_Rv(tps)
p_sat = TDI.saturation_vapor_pressure_over_ice(tps, T)
q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
q_snow_range = range(1e-8, stop = 5e-3, length = 100)
q_tot = 15e-3
q_vap = 0.15 * q_sat
q_liq = 0.0
q_ice = q_tot - q_vap - q_liq
q_rai = 0.0
q = TDI.PP(q_tot, q_liq, q_ice)
R = TDI.Rₘ(tps, q_tot, q_liq, q_ice)
ρ = p / R / T
rate = CM1.evaporation_sublimation.(snow, Blk1MVel.snow, aps, tps, Ref(q), q_rai, q_snow_range, ρ, T)
MK.lines!(q_snow_range * 1e3, rate, label = "T > 0°C")

MK.axislegend(ax; position = :rt)
MK.save("snow_sublimation_deposition_rate.svg", fig) # hide

# snow melt rate figure
fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_snow [g/kg]", ylabel = "snow melt rate [1/s]", limits)
T = 273.15
_snow_melt(ΔT) = CM1.snow_melt.(snow, Blk1MVel.snow, aps, tps, q_snow_range, ρ, T + ΔT)
MK.lines!(q_snow_range * 1e3, _snow_melt(2), label = "T = 2°C")
MK.lines!(q_snow_range * 1e3, _snow_melt(4), label = "T = 4°C")
MK.lines!(q_snow_range * 1e3, _snow_melt(6), label = "T = 6°C")
MK.axislegend(ax; position = :lt)
MK.save("snow_melt_rate.svg", fig) # hide
