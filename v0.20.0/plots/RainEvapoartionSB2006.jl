import CairoMakie as MK

import SpecialFunctions as SF

import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaParams as CP

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Microphysics2M as CM2

FT = Float64

const tps = TD.Parameters.ThermodynamicsParameters(FT)
const aps = CMP.AirProperties(FT)
const SB2006 = CMP.SB2006(FT)

function rain_evaporation_CPU(SB2006, aps, tps, q, q_rai, ρ, N_rai, T)

    evap_rate_0 = FT(0)
    evap_rate_1 = FT(0)
    S = TD.supersaturation(tps, q, ρ, T, TD.Liquid())

    if (q_rai > FT(0) && S < FT(0))

        (; ν_air, D_vapor) = aps
        (; av, bv, α, β, ρ0) = SB2006.evap
        ρw = SB2006.pdf.ρw
        x_star = SB2006.pdf.xr_min
        G = CO.G_func(aps, tps, T, TD.Liquid())

        xr = CM2.raindrops_limited_vars(SB2006.pdf, q_rai, ρ, N_rai).xr
        Dr = (FT(6) / FT(π) / ρw)^FT(1 / 3) * xr^FT(1 / 3)

        t_star = (FT(6) * x_star / xr)^FT(1 / 3)
        a_vent_0 = av * FT(SF.gamma(-1, t_star)) / FT(6)^FT(-2 / 3)
        b_vent_0 =
            bv * FT(SF.gamma((-1 / 2) + (3 / 2) * β, t_star)) /
            FT(6)^FT(β / 2 - 1 / 2)

        a_vent_1 = av * SF.gamma(FT(2)) / FT(6)^FT(1 / 3)
        b_vent_1 =
            bv * SF.gamma(FT(5 / 2) + FT(3 / 2) * β) / FT(6)^FT(β / 2 + 1 / 2)

        N_Re = α * xr^β * sqrt(ρ0 / ρ) * Dr / ν_air
        Fv0 = a_vent_0 + b_vent_0 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)
        Fv1 = a_vent_1 + b_vent_1 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)

        evap_rate_0 = min(FT(0), FT(2) * FT(π) * G * S * N_rai * Dr * Fv0 / xr)
        evap_rate_1 = min(FT(0), FT(2) * FT(π) * G * S * N_rai * Dr * Fv1 / ρ)
    end

    return (; evap_rate_0, evap_rate_1)
end

qᵥ = FT(1e-2)
q = TD.PhasePartition(qᵥ)
qᵣ = FT(1e-4)
ρ = FT(1.1)
Nᵣ = FT(1e8)
T = FT(300)

qᵣ_range = range(1e-9, stop = 1e-4, length = 1000)
Nᵣ_range = range(1e6, stop = 1e9, length = 1000)
T_range = range(273.15, stop = 273.15 + 50, length = 1000)

#! format: off
evap_qᵣ_0 = [rain_evaporation_CPU(SB2006, aps, tps, q, _qᵣ, ρ, Nᵣ, T).evap_rate_0 for _qᵣ in qᵣ_range]
evap_Nᵣ_0 = [rain_evaporation_CPU(SB2006, aps, tps, q, qᵣ, ρ, _Nᵣ, T).evap_rate_0 for _Nᵣ in Nᵣ_range]
evap_T_0 = [rain_evaporation_CPU(SB2006, aps, tps, q, qᵣ, ρ, Nᵣ, _T).evap_rate_0 for _T in T_range]

evap_qᵣ_0n = [CM2.rain_evaporation(SB2006, aps, tps, q, _qᵣ, ρ, Nᵣ, T).evap_rate_0 for _qᵣ in qᵣ_range]
evap_Nᵣ_0n = [CM2.rain_evaporation(SB2006, aps, tps, q, qᵣ, ρ, _Nᵣ, T).evap_rate_0 for _Nᵣ in Nᵣ_range]
evap_T_0n = [CM2.rain_evaporation(SB2006, aps, tps, q, qᵣ, ρ, Nᵣ, _T).evap_rate_0 for _T in T_range]

evap_qᵣ_3 = [rain_evaporation_CPU(SB2006, aps, tps, q, _qᵣ, ρ, Nᵣ, T).evap_rate_1 for _qᵣ in qᵣ_range]
evap_Nᵣ_3 = [rain_evaporation_CPU(SB2006, aps, tps, q, qᵣ, ρ, _Nᵣ, T).evap_rate_1 for _Nᵣ in Nᵣ_range]
evap_T_3 = [rain_evaporation_CPU(SB2006, aps, tps, q, qᵣ, ρ, Nᵣ, _T).evap_rate_1 for _T in T_range]

evap_qᵣ_3n = [CM2.rain_evaporation(SB2006, aps, tps, q, _qᵣ, ρ, Nᵣ, T).evap_rate_1 for _qᵣ in qᵣ_range]
evap_Nᵣ_3n = [CM2.rain_evaporation(SB2006, aps, tps, q, qᵣ, ρ, _Nᵣ, T).evap_rate_1 for _Nᵣ in Nᵣ_range]
evap_T_3n = [CM2.rain_evaporation(SB2006, aps, tps, q, qᵣ, ρ, Nᵣ, _T).evap_rate_1 for _T in T_range]

fig = MK.Figure(resolution = (800, 600))

ax1 = MK.Axis(fig[1, 1], xlabel = "q_rain [g/kg]", ylabel = "evap rate [1/cm3/s]")
ax2 = MK.Axis(fig[2, 1], xlabel = "N_rain [1/cm3]", ylabel = "evap rate [1/cm3/s]")
ax3 = MK.Axis(fig[3, 1], xlabel = "T [K]", ylabel = "evap rate [1/cm3/s]")

ax4 = MK.Axis(fig[1, 2], xlabel = "q_rain [g/kg]", ylabel = "evap rate [g/kg/s]")
ax5 = MK.Axis(fig[2, 2], xlabel = "N_rain [1/cm3]", ylabel = "evap rate [g/kg/s]")
ax6 = MK.Axis(fig[3, 2], xlabel = "T [K]", ylabel = "evap rate [g/kg/s]")

MK.lines!(ax1, qᵣ_range .* 1e3, evap_qᵣ_0 .* 1e-6, color = :blue, label="exact")
MK.lines!(ax1, qᵣ_range .* 1e3, evap_qᵣ_0n .* 1e-6, color = :orange, label="approx")
MK.lines!(ax2, Nᵣ_range .* 1e-6, evap_Nᵣ_0 .* 1e-6, color =:blue)
MK.lines!(ax2, Nᵣ_range .* 1e-6, evap_Nᵣ_0n .* 1e-6, color =:orange)
MK.lines!(ax3, T_range, evap_T_0 .* 1e-6, color =:blue)
MK.lines!(ax3, T_range, evap_T_0n .* 1e-6, color =:orange)
MK.lines!(ax4, qᵣ_range .* 1e3, evap_qᵣ_3 .* 1e3, color =:blue)
MK.lines!(ax4, qᵣ_range .* 1e3, evap_qᵣ_3n .* 1e3, color =:orange)
MK.lines!(ax5, Nᵣ_range .* 1e-6, evap_Nᵣ_3 .* 1e3, color =:blue)
MK.lines!(ax5, Nᵣ_range .* 1e-6, evap_Nᵣ_3n .* 1e3, color =:orange)
MK.lines!(ax6, T_range, evap_T_3 .* 1e3, color =:blue)
MK.lines!(ax6, T_range, evap_T_3n .* 1e3, color =:orange)
#! format: on

MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 12,
    orientation = :horizontal,
    nbanks = 2,
    position = :rt,
)

MK.save("SB2006_rain_evaporation.svg", fig)
