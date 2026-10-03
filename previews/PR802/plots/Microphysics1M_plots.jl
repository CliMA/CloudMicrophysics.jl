import CairoMakie as MK

import ClimaParams
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.ThermodynamicsInterface as TDI

FT = Float64

tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
mp = CMP.Microphysics1MParams(FT)

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

# Example values
q_min, q_max = 1e-8, 5e-3
q_lcl_range = range(q_min, stop = q_max, length = 100)
q_icl_range = range(q_min, stop = q_max, length = 100)
q_rain_range = range(q_min, stop = q_max, length = 100)
q_snow_range = range(q_min, stop = q_max, length = 100)
ρ_air = 1.2
q_lcl, q_icl, q_tot = 5e-4, 5e-4, 20e-3
q_rai = 1e-3
q_sno = 1e-4
T = 273.15
limits = (0, q_max * 1e3, 0, nothing)

# Autoconversion rate figure
mp_ss = CMP.Microphysics1MParams(FT; snow_autoconversion = CMP.WithSupersaturation())
autoconversion_rain = [
    CM1.conv_q_lcl_to_q_rai(
        mp.processes.rain_autoconversion, mp, tps,
        (; q_tot, q_lcl = q, q_icl, q_rai, q_sno),
        (; ρ = ρ_air, T, w = FT(0)),
    ) for q in q_lcl_range
]
autoconversion_snow(T) = [
    CM1.conv_q_icl_to_q_sno(
        mp_ss.processes.snow_autoconversion, mp_ss, tps,
        (; q_tot, q_lcl = FT(0), q_icl = q, q_rai, q_sno),
        (; ρ = ρ_air, T),
    ) for q in q_icl_range
]
MK.with_theme(MK.theme_minimal()) do
    fig = MK.Figure()
    ax = MK.Axis(fig[1, 1]; xlabel = "q_lcl or q_icl [g/kg]", ylabel = "autoconversion rate [1/s]", limits)
    MK.lines!(ax, q_lcl_range * 1e3, autoconversion_rain; label = "Rain")
    MK.lines!(ax, q_icl_range * 1e3, autoconversion_snow(T - 5); label = "Snow T = −5°C")
    MK.lines!(ax, q_icl_range * 1e3, autoconversion_snow(T - 10); label = "Snow T = −10°C")
    MK.lines!(ax, q_icl_range * 1e3, autoconversion_snow(T - 15); label = "Snow T = −15°C")
    MK.axislegend(ax; position = :lt)
    MK.save("autoconversion_rate.svg", fig)
end

# Accretion rate figure
accretion_rain(process) = [
    CM1.accretion(process, mp, tps, (; q_tot, q_lcl, q_icl, q_rai = q, q_sno), (; ρ = ρ_air, T))
    for q in q_rain_range
]
accretion_snow(process) = [
    CM1.accretion(process, mp, tps, (; q_tot, q_lcl, q_icl, q_rai, q_sno = q), (; ρ = ρ_air, T))
    for q in q_snow_range
]
accretion_liq_rain = accretion_rain(mp.processes.cloud_liquid_rain_accretion)
accretion_ice_rain = accretion_rain(mp.processes.cloud_ice_rain_accretion)
accretion_liq_snow = getproperty.(accretion_snow(mp.processes.cloud_liquid_snow_accretion), :S_accr)
accretion_ice_snow = accretion_snow(mp.processes.cloud_ice_snow_accretion)
MK.with_theme(MK.theme_minimal()) do
    fig = MK.Figure()
    ax = MK.Axis(fig[1, 1]; xlabel = "q_rain or q_snow [g/kg]", ylabel = "accretion rate [1/s]", limits)
    MK.lines!(ax, q_rain_range * 1e3, accretion_liq_rain; label = "Liq+Rain-CliMA")
    MK.lines!(ax, q_rain_range * 1e3, accretion_ice_rain; label = "Ice+Rain-CliMA")
    MK.lines!(ax, q_snow_range * 1e3, accretion_liq_snow; label = "Liq+Snow-CliMA")
    MK.lines!(ax, q_snow_range * 1e3, accretion_ice_snow; label = "Ice+Snow-CliMA", linewidth = 4, linestyle = :dash)
    MK.lines!(ax, q_rain_range * 1e3, accretion_empirical.(q_rain_range, q_lcl, q_tot); label = "Liq+Rain-Empirical")
    MK.axislegend(ax; position = :lt)
    MK.save("accretion_rate.svg", fig)
end

# Accretion rain sink rate figure
accretion_rain_sink(q_icl) = [
    CM1.accretion_rain_sink(
        mp.processes.cloud_ice_rain_accretion, mp, tps,
        (; q_tot, q_lcl, q_icl, q_rai = q, q_sno),
        (; ρ = ρ_air, T),
    ) for q in q_rain_range
]
MK.with_theme(MK.theme_minimal()) do
    fig = MK.Figure()
    ax = MK.Axis(fig[1, 1]; xlabel = "q_rain or q_snow [g/kg]", ylabel = "accretion rain sink rate [1/s]", limits)
    MK.lines!(ax, q_rain_range * 1e3, accretion_rain_sink(1e-6); label = "q_icl = 1e-6")
    MK.lines!(ax, q_rain_range * 1e3, accretion_rain_sink(1e-5); label = "q_icl = 1e-5")
    MK.lines!(ax, q_rain_range * 1e3, accretion_rain_sink(1e-4); label = "q_icl = 1e-4")
    MK.axislegend(ax; position = :lt)
    MK.save("accretion_rain_sink_rate.svg", fig)
end

# Snow-rain accretion rate figures, above freezing (snow → rain) and below freezing (rain → snow)
accretion_snow_rain(q_rai, q_sno) = CM1.accretion_snow_rain(
    mp.processes.rain_snow_accretion, mp, tps,
    (; q_tot, q_lcl, q_icl, q_rai, q_sno),
    (; ρ = ρ_air, T),
)
accretion_snow_rain_warm(q_sno) = [accretion_snow_rain(q, q_sno).S_sno_rai for q in q_rain_range]
accretion_snow_rain_cold(q_rai) = [accretion_snow_rain(q_rai, q).S_rai_sno for q in q_snow_range]
MK.with_theme(MK.theme_minimal()) do
    fig = MK.Figure()
    ax = MK.Axis(fig[1, 1]; xlabel = "q_rain [g/kg]", ylabel = "snow-rain accretion rate [1/s] T>0", limits)
    MK.lines!(ax, q_rain_range * 1e3, accretion_snow_rain_warm(1e-6); label = "q_snow = 1e-6")
    MK.lines!(ax, q_rain_range * 1e3, accretion_snow_rain_warm(1e-5); label = "q_snow = 1e-5")
    MK.lines!(ax, q_rain_range * 1e3, accretion_snow_rain_warm(1e-4); label = "q_snow = 1e-4")
    MK.axislegend(ax; position = :lt)
    MK.save("accretion_snow_rain_above_freeze.svg", fig)
end
MK.with_theme(MK.theme_minimal()) do
    fig = MK.Figure()
    ax = MK.Axis(fig[1, 1]; xlabel = "q_snow [g/kg]", ylabel = "snow-rain accretion rate [1/s] T<0", limits)
    MK.lines!(ax, q_snow_range * 1e3, accretion_snow_rain_cold(1e-6); label = "q_rain = 1e-6")
    MK.lines!(ax, q_snow_range * 1e3, accretion_snow_rain_cold(1e-5); label = "q_rain = 1e-5")
    MK.lines!(ax, q_snow_range * 1e3, accretion_snow_rain_cold(1e-4); label = "q_rain = 1e-4")
    MK.axislegend(ax; position = :lt)
    MK.save("accretion_snow_rain_below_freeze.svg", fig)
end

# Rain evaporation rate figure
ϵ = TDI.Rd_over_Rv(tps)
T_evap, p_evap = 273.15 + 15, 90000.0
p_sat_evap = TDI.saturation_vapor_pressure_over_liquid(tps, T_evap)
q_sat_evap = ϵ * p_sat_evap / (p_evap + p_sat_evap * (ϵ - 1.0))
q_tot_evap = 15e-3
q_lcl_evap = q_tot_evap - 0.15 * q_sat_evap
ρ_evap = p_evap / TDI.Rₘ(tps, q_tot_evap, q_lcl_evap + q_rai, 0.0) / T_evap
evaporation = [
    CM1.conv_q_rai_to_q_vap(
        mp.processes.rain_condensation_evaporation, mp, tps,
        (; q_tot = q_tot_evap, q_lcl = q_lcl_evap - q, q_icl = 0.0, q_rai = q, q_sno = 0.0),
        (; ρ = ρ_evap, T = T_evap),
    ) for q in q_rain_range
]
evaporation_empirical = rain_evap_empirical.(tps, q_rain_range, q_tot_evap, q_lcl_evap, T_evap, p_evap, ρ_evap)
MK.with_theme(MK.theme_minimal()) do
    fig = MK.Figure()
    ax = MK.Axis(
        fig[1, 1];
        xlabel = "q_rain [g/kg]",
        ylabel = "rain evaporation rate [1/s]",
        limits = (0, q_max * 1e3, nothing, nothing),
    )
    MK.lines!(ax, q_rain_range * 1e3, evaporation; label = "ClimateMachine")
    MK.lines!(ax, q_rain_range * 1e3, evaporation_empirical; label = "empirical")
    MK.axislegend(ax; position = :rt)
    MK.save("rain_evaporation_rate.svg", fig)
end

# Snow deposition and sublimation rate figure, for air at 15% of saturation over ice
function snow_air(T; p = 90000.0, q_tot = 15e-3)
    p_sat = TDI.saturation_vapor_pressure_over_ice(tps, T)
    q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
    q_icl = q_tot - 0.15 * q_sat
    ρ = p / TDI.Rₘ(tps, q_tot, 0.0, q_icl) / T
    return (; T, q_tot, q_icl, ρ)
end
snow_sublimation(air) = [
    CM1.conv_q_sno_to_q_vap(
        mp.processes.snow_deposition_sublimation, mp, tps,
        (; air.q_tot, q_lcl = 0.0, q_icl = air.q_icl - q, q_rai = 0.0, q_sno = q),
        (; air.ρ, air.T),
    ) for q in q_snow_range
]
air_cold, air_warm = snow_air(273.15 - 15), snow_air(273.15 + 15)
MK.with_theme(MK.theme_minimal()) do
    fig = MK.Figure()
    ax = MK.Axis(
        fig[1, 1];
        xlabel = "q_snow [g/kg]",
        ylabel = "snow deposition sublimation rate [1/s]",
        limits = (0, q_max * 1e3, nothing, nothing),
    )
    MK.lines!(ax, q_snow_range * 1e3, snow_sublimation(air_cold); label = "T < 0°C")
    MK.lines!(ax, q_snow_range * 1e3, snow_sublimation(air_warm); label = "T > 0°C")
    MK.axislegend(ax; position = :rt)
    MK.save("snow_sublimation_deposition_rate.svg", fig)
end

# Snow melt rate figure
snow_melt(ΔT) = [
    CM1.conv_q_sno_to_q_rai(
        mp.processes.snow_melt, mp, tps,
        (; q_tot = FT(0), q_lcl = FT(0), q_icl = FT(0), q_rai = FT(0), q_sno = q),
        (; air_warm.ρ, T = T + ΔT),
    ) for q in q_snow_range
]
MK.with_theme(MK.theme_minimal()) do
    fig = MK.Figure()
    ax = MK.Axis(fig[1, 1]; xlabel = "q_snow [g/kg]", ylabel = "snow melt rate [1/s]", limits)
    MK.lines!(ax, q_snow_range * 1e3, snow_melt(2); label = "T = 2°C")
    MK.lines!(ax, q_snow_range * 1e3, snow_melt(4); label = "T = 4°C")
    MK.lines!(ax, q_snow_range * 1e3, snow_melt(6); label = "T = 6°C")
    MK.axislegend(ax; position = :lt)
    MK.save("snow_melt_rate.svg", fig)
end

# Heterogeneous freezing rate figures, versus q_lcl and versus T
mp_with_N_0(N_0) = CMP.Microphysics1MParams(;
    processes = mp.processes,
    process_params = mp.process_params,
    cloud = CMP.CloudPhaseParams1M(;
        liquid = CMP.CloudLiquid(; ρw = mp.cloud.liquid.ρw, r_eff = mp.cloud.liquid.r_eff, N_0),
        ice = mp.cloud.ice,
    ),
    precip = mp.precip,
    air_properties = mp.air_properties,
    terminal_velocity = mp.terminal_velocity,
)
q_tot_het = 15e-3
heterogeneous_freezing(N_0, q_lcl, T) = CM1.conv_q_lcl_to_q_icl(
    CMP.Heterogeneous(), mp_with_N_0(N_0), tps,
    (; q_tot = q_tot_het, q_lcl, q_icl = FT(0), q_rai = FT(0), q_sno = FT(0)),
    (; ρ = ρ_air, T),
)
N_0_values = [(FT(1e7), "N₀ = 10⁷ m⁻³"), (FT(1e8), "N₀ = 10⁸ m⁻³"), (FT(5e8), "N₀ = 5×10⁸ m⁻³")]
T_het = FT(273.15 - 15)
MK.with_theme(MK.theme_minimal()) do
    fig = MK.Figure()
    ax = MK.Axis(fig[1, 1]; xlabel = "q_lcl [g/kg]", ylabel = "het. freezing rate [1/s]", limits)
    for (N_0, label) in N_0_values
        MK.lines!(ax, q_lcl_range * 1e3, heterogeneous_freezing.(N_0, q_lcl_range, T_het); label)
    end
    MK.axislegend(ax; position = :lt)
    MK.save("het_freezing_rate.svg", fig)
end
T_range = range(FT(273 - 40), stop = FT(273), length = 100)
MK.with_theme(MK.theme_minimal()) do
    fig = MK.Figure()
    ax = MK.Axis(fig[1, 1]; xlabel = "T [K]", ylabel = "het. freezing rate [1/s]")
    for (N_0, label) in N_0_values
        MK.lines!(ax, T_range, heterogeneous_freezing.(N_0, FT(1e-4), T_range); label)
    end
    MK.axislegend(ax; position = :lt)
    MK.save("het_freezing_rate_vs_T.svg", fig)
end

# Velocity-dependent Kessler autoconversion: τ(w), q_threshold(w) and rate versus w,
# with illustrative quiescent-regime values that differ from the convective ones
toml_vd = ClimaParams.create_toml_dict(FT;
    override_file = Dict(
        "rain_autoconversion_timescale_stratiform" => Dict("value" => 14400.0, "type" => "float"),
        "cloud_liquid_water_specific_humidity_autoconversion_threshold_stratiform" =>
            Dict("value" => 1e-3, "type" => "float"),
    ),
)
mp_vd = CMP.Microphysics1MParams(toml_vd)
w_range = range(FT(-6), stop = FT(6), length = 200)
τ_vals = [CM1.rain_autoconversion_timescale(mp_vd.processes.rain_autoconversion, mp_vd, w) / 3600 for w in w_range]
qt_vals = [CM1.rain_autoconversion_threshold(mp_vd.processes.rain_autoconversion, mp_vd, w) * 1000 for w in w_range]
autoconversion_vs_w(q) = [
    CM1.conv_q_lcl_to_q_rai(
        mp_vd.processes.rain_autoconversion, mp_vd, tps,
        (; q_tot = FT(0), q_lcl = q, q_icl = FT(0), q_rai = FT(0), q_sno = FT(0)),
        (; ρ = FT(1.2), T = FT(280), w),
    ) for w in w_range
]
MK.with_theme(MK.theme_minimal()) do
    fig = MK.Figure(size = (1200, 400))
    ax1 = MK.Axis(fig[1, 1]; xlabel = "w [m/s]", ylabel = "τ(w) [hours]", title = "Effective autoconversion timescale")
    MK.lines!(ax1, collect(w_range), τ_vals)
    ax2 = MK.Axis(
        fig[1, 2];
        xlabel = "w [m/s]",
        ylabel = "q_threshold(w) [g/kg]",
        title = "Effective autoconversion threshold",
    )
    MK.lines!(ax2, collect(w_range), qt_vals)
    ax3 = MK.Axis(fig[1, 3]; xlabel = "w [m/s]", ylabel = "autoconversion rate [1/s]", title = "Autoconversion rate")
    for (q, label) in [(FT(5e-4), "q_lcl = 0.5 g/kg"), (FT(1e-3), "q_lcl = 1.0 g/kg"), (FT(2e-3), "q_lcl = 2.0 g/kg")]
        MK.lines!(ax3, collect(w_range), autoconversion_vs_w(q); label)
    end
    MK.axislegend(ax3; position = :ct)
    MK.save("velocity_dependent_autoconversion.svg", fig)
end
nothing
