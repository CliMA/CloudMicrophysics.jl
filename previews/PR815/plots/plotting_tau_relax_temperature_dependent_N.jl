import CairoMakie as MK
CairoMakie = MK

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMNe

FT = Float64

# Parameters
aps = CMP.AirProperties(FT)
ice = CMP.CloudIce(FT)
fit = CMP.IceNumberTemperatureFit(FT)
ρ = FT(0.8)   # representative mid-troposphere air density

# Temperature range: -60 °C to +5 °C
T_range = range(fit.T_freeze - 60, fit.T_freeze + 5, length = 300)
T_celsius = T_range .- fit.T_freeze

# Ice number concentration N_ice(T)
N_ice = [CMNe.ice_number_concentration(fit, T) for T in T_range]

# Relaxation timescale for several cloud ice specific contents, compared with the
# prescribed N₀ = 5×10⁸ m⁻³ at the same q_icl
q_icl_values = [FT(1e-7), FT(1e-6), FT(1e-5), FT(1e-4)]
q_icl_labels = ["10⁻⁷", "10⁻⁶", "10⁻⁵", "10⁻⁴"]
τ_fit = [CMNe.τ_relax(ice, aps, fit, q_icl, T, ρ) for T in T_range, q_icl in q_icl_values]
τ_N0 = [CMNe.τ_relax(ice, aps, q_icl, ρ) for q_icl in q_icl_values]

# Colors — distinct, colorblind-friendly
colors = [
    MK.RGBf(0.12, 0.47, 0.71),  # blue
    MK.RGBf(1.0, 0.50, 0.05),   # orange
    MK.RGBf(0.17, 0.63, 0.17),  # green
    MK.RGBf(0.84, 0.15, 0.16),  # red
]
linestyles = [:solid, :dash, :dot, :dashdot]

fig = MK.Figure(size = (800, 900))
ax1 = MK.Axis(
    fig[1, 1],
    xlabel = "Temperature [°C]",
    ylabel = "N_ice [m⁻³]",
    title = "Temperature-dependent cloud ice number concentration",
    yscale = log10,
)
MK.lines!(ax1, T_celsius, N_ice, color = colors[1], linewidth = 2.5, label = "N_ice(T)")
MK.hlines!(ax1, [ice.N_0], color = :gray40, linestyle = :dash, linewidth = 2, label = "prescribed N₀ = $(ice.N_0) m⁻³")
MK.axislegend(ax1, position = :rt)

ax2 = MK.Axis(
    fig[2, 1],
    xlabel = "Temperature [°C]",
    ylabel = "Relaxation timescale τ [s]",
    title = "Ice relaxation timescale — temperature-dependent N_ice (ρ = $(ρ) kg/m³)\nsolid/dashed: N_ice(T); thin dotted: prescribed N₀ at the same q_icl",
    yscale = log10,
)
for (j, label) in enumerate(q_icl_labels)
    MK.lines!(
        ax2, T_celsius, τ_fit[:, j],
        label = "q_icl = $(label) kg/kg",
        color = colors[j],
        linewidth = 2.5,
        linestyle = linestyles[j],
    )
    MK.hlines!(ax2, [τ_N0[j]], color = colors[j], linestyle = :dot, linewidth = 1.2)
end
MK.hlines!(ax2, [FT(3600)], color = :gray40, linestyle = :dash, linewidth = 1.5, label = "1 hour")
MK.axislegend(ax2, position = :lt)

MK.save("tau_relax_temperature_dependent_N.svg", fig)
println("Plot saved to tau_relax_temperature_dependent_N.svg")
