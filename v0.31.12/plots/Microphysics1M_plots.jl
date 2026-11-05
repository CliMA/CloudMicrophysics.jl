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
function accretion_empirical(q_rai::DT, q_lcl::DT, q_tot::DT) where {DT <: Real}
    rr = q_rai / (DT(1) - q_tot)
    rl = q_lcl / (DT(1) - q_tot)
    return DT(2.2) * rl * rr^DT(7 / 8)
end

# eq. 5c in [Grabowski1996](@cite)
function rain_evap_empirical(tps, q_rai, q_tot, q_lcl, T, p, ρ)
    DT = eltype(q_rai)

    p_v_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)
    q_sat = TDI.p2q(tps, T, ρ, p_v_sat)

    q_vap = q_tot - q_lcl
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
q_lcl_range = range(q_min, stop = q_max, length = 100)
q_icl_range = range(q_min, stop = q_max, length = 100)
q_rain_range = range(q_min, stop = q_max, length = 100)
q_snow_range = range(q_min, stop = q_max, length = 100)
ρ_air, ρ_air_ground = 1.2, 1.22
q_lcl, q_icl, q_tot = 5e-4, 5e-4, 20e-3
q_rai = 1e-3
q_sno = 1e-4
limits = (0, q_max * 1e3, 0, nothing)

MK.set_theme!(MK.theme_minimal())

# autoconversion rate figure
T = 273.15
fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_lcl or q_icl [g/kg]", ylabel = "autoconversion rate [1/s]", limits)
q_icl_to_q_sno_rate = T -> CM1.conv_q_icl_to_q_sno.(ice, aps, tps, q_tot, FT(0), q_icl_range, q_rai, q_sno, ρ_air, T)
MK.lines!(q_lcl_range * 1e3, CM1.conv_q_lcl_to_q_rai.(rain.acnv1M, q_lcl_range), label = "Rain")
MK.lines!(q_icl_range * 1e3, q_icl_to_q_sno_rate(T - 5), label = "Snow T = −5°C")
MK.lines!(q_icl_range * 1e3, q_icl_to_q_sno_rate(T - 10), label = "Snow T = −10°C")
MK.lines!(q_icl_range * 1e3, q_icl_to_q_sno_rate(T - 15), label = "Snow T = −15°C")
MK.axislegend(ax; position = :lt)
MK.save("autoconversion_rate.svg", fig) # hide

# accretion rate figure
fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_rain or q_snow [g/kg]", ylabel = "accretion rate [1/s]", limits)
MK.lines!(
    q_rain_range * 1e3, CM1.accretion.(liquid, rain, Blk1MVel.rain, ce, q_lcl, q_rain_range, ρ_air),
    label = "Liq+Rain-CliMA",
)
MK.lines!(
    q_rain_range * 1e3, CM1.accretion.(ice, rain, Blk1MVel.rain, ce, q_icl, q_rain_range, ρ_air),
    label = "Ice+Rain-CliMA",
)
MK.lines!(
    q_snow_range * 1e3, CM1.accretion.(liquid, snow, Blk1MVel.snow, ce, q_lcl, q_snow_range, ρ_air),
    label = "Liq+Snow-CliMA",
)
MK.lines!(
    q_snow_range * 1e3, CM1.accretion.(ice, snow, Blk1MVel.snow, ce, q_icl, q_snow_range, ρ_air),
    label = "Ice+Snow-CliMA", linewidth = 4, linestyle = :dash,
)
MK.lines!(
    q_rain_range * 1e3, accretion_empirical.(q_rain_range, q_lcl, q_tot),
    label = "Liq+Rain-Empirical",
)
MK.axislegend(ax; position = :lt)
MK.save("accretion_rate.svg", fig) # hide

# accretion rain sink rate figure
fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_rain or q_snow [g/kg]", ylabel = "accretion rain sink rate [1/s]", limits)
_accr_rain_sink(q_icl) = CM1.accretion_rain_sink.(rain, ice, Blk1MVel.rain, ce, q_icl, q_rain_range, ρ_air) # hide
MK.lines!(q_rain_range * 1e3, _accr_rain_sink(1e-6), label = "q_icl = 1e-6")
MK.lines!(q_rain_range * 1e3, _accr_rain_sink(1e-5), label = "q_icl = 1e-5")
MK.lines!(q_rain_range * 1e3, _accr_rain_sink(1e-4), label = "q_icl = 1e-4")
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
q_icl = 0.0
q_lcl = q_tot - q_vap - q_icl
q_sno = 0.0
R = TDI.Rₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)
ρ = p / R / T

fig = MK.Figure()
ax = MK.Axis(fig[1, 1]; xlabel = "q_rain [g/kg]", ylabel = "rain evaporation rate [1/s]")
MK.xlims!(ax, 0, q_max * 1e3)
MK.lines!(
    q_rain_range * 1e3,
    CM1.evaporation_sublimation.(
        rain,
        Blk1MVel.rain,
        aps,
        tps,
        q_tot,
        q_lcl .- q_rain_range,
        q_icl,
        q_rain_range,
        q_sno,
        ρ,
        T,
    ),
    label = "ClimateMachine",
)
MK.lines!(q_rain_range * 1e3, rain_evap_empirical.(tps, q_rain_range, q_tot, q_lcl, T, p, ρ), label = "empirical")
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
q_lcl = 0.0
q_icl = q_tot - q_vap - q_lcl
q_rai = 0.0
R = TDI.Rₘ(tps, q_tot, q_lcl + q_rai, q_icl)
ρ = p / R / T
rate =
    CM1.evaporation_sublimation.(
        snow,
        Blk1MVel.snow,
        aps,
        tps,
        q_tot,
        q_lcl,
        q_icl .- q_snow_range,
        q_rai,
        q_snow_range,
        ρ,
        T,
    )
MK.lines!(q_snow_range * 1e3, rate, label = "T < 0°C")

# example values 2
T, p = 273.15 + 15, 90000.0
ϵ = TDI.Rd_over_Rv(tps)
p_sat = TDI.saturation_vapor_pressure_over_ice(tps, T)
q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
q_snow_range = range(1e-8, stop = 5e-3, length = 100)
q_tot = 15e-3
q_vap = 0.15 * q_sat
q_lcl = 0.0
q_icl = q_tot - q_vap - q_lcl
q_rai = 0.0
R = TDI.Rₘ(tps, q_tot, q_lcl + q_rai, q_icl)
ρ = p / R / T
rate =
    CM1.evaporation_sublimation.(
        snow,
        Blk1MVel.snow,
        aps,
        tps,
        q_tot,
        q_lcl,
        q_icl .- q_snow_range,
        q_rai,
        q_snow_range,
        ρ,
        T,
    )
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

# -------------------------------------------------------------------
# Derivative plots: tendency and ∂tendency/∂q for rain evap, snow subl, snow melt
# -------------------------------------------------------------------

T_freeze = TDI.TD.Parameters.T_freeze(tps)

# ---------- Rain evaporation vs temperature ----------
fig = MK.Figure(size = (900, 450))
ax1 = MK.Axis(fig[1, 1];
    xlabel = "Temperature (K)", ylabel = "Evaporation rate (1/s)",
    title = "Rain evaporation vs T\n(q_rai=1e-3, q_tot=15e-3 g/kg, 15% RH, ρ from p=90kPa)",
)
ax2 = MK.Axis(fig[1, 2];
    xlabel = "Temperature (K)", ylabel = "∂(rate)/∂q_rai (1/s)",
    title = "Derivative of rain evaporation",
)

T_range = collect(range(250.0, stop = 310.0, length = 121))
for q_rai_val in [FT(5e-4), FT(1e-3), FT(2e-3)]
    rates = zeros(FT, length(T_range))
    drates = zeros(FT, length(T_range))
    for (i, T) in enumerate(T_range)
        p = FT(90000)
        ϵ = TDI.Rd_over_Rv(tps)
        p_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))
        q_tot = FT(15e-3)
        q_vap = FT(0.15) * q_sat
        q_lcl = max(q_tot - q_vap - q_rai_val, FT(0))
        R = TDI.Rₘ(tps, q_tot, q_lcl + q_rai_val, FT(0))
        ρ_loc = p / R / T
        r = CM1.evaporation_sublimation(
            rain, Blk1MVel.rain, aps, tps, q_tot, q_lcl, FT(0), q_rai_val, FT(0), ρ_loc, T,
        )
        dr = CM1.∂evaporation_sublimation_∂q_precip(
            rain, Blk1MVel.rain, aps, tps, q_tot, q_lcl, FT(0), q_rai_val, FT(0), ρ_loc, T,
        )
        rates[i] = r
        drates[i] = dr
    end
    lbl = "q_rai=$(q_rai_val*1e3) g/kg"
    MK.lines!(ax1, T_range, rates; label = lbl)
    MK.lines!(ax2, T_range, drates; label = lbl)
end
MK.hlines!(ax1, [0]; color = :black, linewidth = 0.5)
MK.axislegend(ax1; position = :lb, framevisible = false)
MK.axislegend(ax2; position = :lb, framevisible = false)
MK.save("rain_evap_deriv_vs_T.svg", fig) # hide

# ---------- Rain evaporation vs q_rai ----------
fig = MK.Figure(size = (900, 450))
ax1 = MK.Axis(fig[1, 1];
    xlabel = "q_rai (g/kg)", ylabel = "Evaporation rate (1/s)",
    title = "Rain evaporation vs q_rai",
)
ax2 = MK.Axis(fig[1, 2];
    xlabel = "q_rai (g/kg)", ylabel = "∂(rate)/∂q_rai (1/s)",
    title = "Derivative of rain evaporation",
)

q_rai_range = collect(range(FT(1e-6), stop = FT(5e-3), length = 100))
for T in [FT(275), FT(290), FT(305)]
    p = FT(90000)
    ϵ = TDI.Rd_over_Rv(tps)
    p_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)
    q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))
    q_tot = FT(15e-3)
    q_vap = FT(0.15) * q_sat
    rates = zeros(FT, length(q_rai_range))
    drates = zeros(FT, length(q_rai_range))
    for (i, q_rai_val) in enumerate(q_rai_range)
        q_lcl = max(q_tot - q_vap - q_rai_val, FT(0))
        R = TDI.Rₘ(tps, q_tot, q_lcl + q_rai_val, FT(0))
        ρ_loc = p / R / T
        r = CM1.evaporation_sublimation(
            rain, Blk1MVel.rain, aps, tps, q_tot, q_lcl, FT(0), q_rai_val, FT(0), ρ_loc, T,
        )
        dr = CM1.∂evaporation_sublimation_∂q_precip(
            rain, Blk1MVel.rain, aps, tps, q_tot, q_lcl, FT(0), q_rai_val, FT(0), ρ_loc, T,
        )
        rates[i] = r
        drates[i] = dr
    end
    lbl = "T=$(Int(T))K"
    MK.lines!(ax1, q_rai_range * 1e3, rates; label = lbl)
    MK.lines!(ax2, q_rai_range * 1e3, drates; label = lbl)
end
MK.hlines!(ax1, [0]; color = :black, linewidth = 0.5)
MK.axislegend(ax1; position = :lb, framevisible = false)
MK.axislegend(ax2; position = :lt, framevisible = false)
MK.save("rain_evap_deriv_vs_qrai.svg", fig) # hide

# ---------- Snow sublimation/deposition vs temperature ----------
fig = MK.Figure(size = (900, 450))
ax1 = MK.Axis(fig[1, 1];
    xlabel = "Temperature (K)", ylabel = "Subl/Dep rate (1/s)",
    title = "Snow sublimation/deposition vs T\n(q_sno=1e-4, subsaturated over ice)",
)
ax2 = MK.Axis(fig[1, 2];
    xlabel = "Temperature (K)", ylabel = "∂(rate)/∂q_sno (1/s)",
    title = "Derivative of snow subl/dep",
)

T_range_cold = collect(range(220.0, stop = 280.0, length = 121))
for (q_sno_val, rh_ice) in [(FT(5e-5), FT(0.8)), (FT(1e-4), FT(0.8)), (FT(1e-4), FT(1.05))]
    rates = zeros(FT, length(T_range_cold))
    drates = zeros(FT, length(T_range_cold))
    for (i, T) in enumerate(T_range_cold)
        p = FT(90000)
        ϵ = TDI.Rd_over_Rv(tps)
        p_sat = TDI.saturation_vapor_pressure_over_ice(tps, T)
        q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1))
        q_tot = rh_ice * q_sat + q_sno_val
        R = TDI.Rₘ(tps, q_tot, FT(0), q_sno_val)
        ρ_loc = p / R / T
        r = CM1.evaporation_sublimation(
            snow, Blk1MVel.snow, aps, tps, q_tot, FT(0), FT(0), FT(0), q_sno_val, ρ_loc, T,
        )
        dr = CM1.∂evaporation_sublimation_∂q_precip(
            snow, Blk1MVel.snow, aps, tps, q_tot, FT(0), FT(0), FT(0), q_sno_val, ρ_loc, T,
        )
        rates[i] = r
        drates[i] = dr
    end
    rh_pct = Int(round(rh_ice * 100))
    lbl = "q_sno=$(q_sno_val*1e3)g/kg, RH_ice=$(rh_pct)%"
    MK.lines!(ax1, T_range_cold, rates; label = lbl)
    MK.lines!(ax2, T_range_cold, drates; label = lbl)
end
MK.hlines!(ax1, [0]; color = :black, linewidth = 0.5)
MK.axislegend(ax1; position = :lb, framevisible = false)
MK.axislegend(ax2; position = :lb, framevisible = false)
MK.save("snow_subl_deriv_vs_T.svg", fig) # hide

# ---------- Snow sublimation vs q_sno ----------
fig = MK.Figure(size = (900, 450))
ax1 = MK.Axis(fig[1, 1];
    xlabel = "q_sno (g/kg)", ylabel = "Sublimation rate (1/s)",
    title = "Snow sublimation vs q_sno\n(subsaturated: RH_ice = 80%)",
)
ax2 = MK.Axis(fig[1, 2];
    xlabel = "q_sno (g/kg)", ylabel = "∂(rate)/∂q_sno (1/s)",
    title = "Derivative of snow sublimation",
)

q_sno_range = collect(range(FT(1e-6), stop = FT(5e-3), length = 100))
for T in [FT(240), FT(255), FT(270)]
    p = FT(90000)
    ϵ = TDI.Rd_over_Rv(tps)
    p_sat_ice = TDI.saturation_vapor_pressure_over_ice(tps, T)
    q_sat_ice = ϵ * p_sat_ice / (p + p_sat_ice * (ϵ - 1))
    rates = zeros(FT, length(q_sno_range))
    drates = zeros(FT, length(q_sno_range))
    for (i, q_sno_val) in enumerate(q_sno_range)
        q_tot = FT(0.8) * q_sat_ice + q_sno_val
        R = TDI.Rₘ(tps, q_tot, FT(0), q_sno_val)
        ρ_loc = p / R / T
        r = CM1.evaporation_sublimation(
            snow, Blk1MVel.snow, aps, tps, q_tot, FT(0), FT(0), FT(0), q_sno_val, ρ_loc, T,
        )
        dr = CM1.∂evaporation_sublimation_∂q_precip(
            snow, Blk1MVel.snow, aps, tps, q_tot, FT(0), FT(0), FT(0), q_sno_val, ρ_loc, T,
        )
        rates[i] = r
        drates[i] = dr
    end
    lbl = "T=$(Int(T))K"
    MK.lines!(ax1, q_sno_range * 1e3, rates; label = lbl)
    MK.lines!(ax2, q_sno_range * 1e3, drates; label = lbl)
end
MK.hlines!(ax1, [0]; color = :black, linewidth = 0.5)
MK.axislegend(ax1; position = :lb, framevisible = false)
MK.axislegend(ax2; position = :lt, framevisible = false)
MK.save("snow_subl_deriv_vs_qsno.svg", fig) # hide

# ---------- Snow melt vs temperature ----------
fig = MK.Figure(size = (900, 450))
ax1 = MK.Axis(fig[1, 1];
    xlabel = "Temperature (K)", ylabel = "Melt rate (1/s)",
    title = "Snow melt rate vs T",
)
ax2 = MK.Axis(fig[1, 2];
    xlabel = "Temperature (K)", ylabel = "∂(rate)/∂q_sno (1/s)",
    title = "Derivative of snow melt",
)

T_range_warm = collect(range(T_freeze - 2, stop = T_freeze + 10, length = 121))
for q_sno_val in [FT(5e-5), FT(1e-4), FT(5e-4)]
    rates = zeros(FT, length(T_range_warm))
    drates = zeros(FT, length(T_range_warm))
    ρ_loc = FT(1.2)
    for (i, T) in enumerate(T_range_warm)
        r = CM1.snow_melt(snow, Blk1MVel.snow, aps, tps, q_sno_val, ρ_loc, T)
        dr = CM1.∂snow_melt_∂q_sno(snow, Blk1MVel.snow, aps, tps, q_sno_val, ρ_loc, T)
        rates[i] = r
        drates[i] = dr
    end
    lbl = "q_sno=$(q_sno_val*1e3) g/kg"
    MK.lines!(ax1, T_range_warm, rates; label = lbl)
    MK.lines!(ax2, T_range_warm, drates; label = lbl)
end
MK.hlines!(ax1, [0]; color = :black, linewidth = 0.5)
MK.axislegend(ax1; position = :lt, framevisible = false)
MK.axislegend(ax2; position = :lt, framevisible = false)
MK.save("snow_melt_deriv_vs_T.svg", fig) # hide

# ---------- Snow melt vs q_sno ----------
fig = MK.Figure(size = (900, 450))
ax1 = MK.Axis(fig[1, 1];
    xlabel = "q_sno (g/kg)", ylabel = "Melt rate (1/s)",
    title = "Snow melt rate vs q_sno",
)
ax2 = MK.Axis(fig[1, 2];
    xlabel = "q_sno (g/kg)", ylabel = "∂(rate)/∂q_sno (1/s)",
    title = "Derivative of snow melt",
)

q_sno_range2 = collect(range(FT(1e-6), stop = FT(5e-3), length = 100))
ρ_loc = FT(1.2)
for ΔT in [FT(2), FT(4), FT(6)]
    T = T_freeze + ΔT
    rates = zeros(FT, length(q_sno_range2))
    drates = zeros(FT, length(q_sno_range2))
    for (i, q_sno_val) in enumerate(q_sno_range2)
        r = CM1.snow_melt(snow, Blk1MVel.snow, aps, tps, q_sno_val, ρ_loc, T)
        dr = CM1.∂snow_melt_∂q_sno(snow, Blk1MVel.snow, aps, tps, q_sno_val, ρ_loc, T)
        rates[i] = r
        drates[i] = dr
    end
    lbl = "T = T_freeze + $(Int(ΔT))K"
    MK.lines!(ax1, q_sno_range2 * 1e3, rates; label = lbl)
    MK.lines!(ax2, q_sno_range2 * 1e3, drates; label = lbl)
end
MK.hlines!(ax1, [0]; color = :black, linewidth = 0.5)
MK.axislegend(ax1; position = :lt, framevisible = false)
MK.axislegend(ax2; position = :lt, framevisible = false)
MK.save("snow_melt_deriv_vs_qsno.svg", fig) # hide
