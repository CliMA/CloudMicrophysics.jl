import CairoMakie as MK
CairoMakie = MK

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMNe

FT = Float64

# Parameters
aps = CMP.AirProperties(FT)
ρ = FT(0.8)   # representative mid-troposphere air density

# Cloud ice density (from defaults)
ice_default = CMP.CloudIce(FT)
ρᵢ = ice_default.ρᵢ

# Build CloudIce structs with different N_0 values
# (reuse the default struct's other fields, only vary N_0)
N_0_values = [FT(1e7), FT(5e7), FT(1e8), FT(5e8)]   # 1/m³
N_0_labels = ["10⁷", "5×10⁷", "10⁸", "5×10⁸"]

# q_icl range: 10⁻⁸ to 10⁻² kg/kg (log-spaced)
q_icl_range = 10 .^ range(-8, -2, length = 300)

# Compute τ for each (N_0, q_icl) pair
τ_data = Matrix{FT}(undef, length(q_icl_range), length(N_0_values))
for (j, N_0) in enumerate(N_0_values)
    ice_j = CMP.CloudIce(;
        pdf = ice_default.pdf,
        mass = ice_default.mass,
        ρᵢ = ice_default.ρᵢ,
        r_eff = ice_default.r_eff,
        N_0 = N_0,
    )
    for (i, q_icl) in enumerate(q_icl_range)
        τ_data[i, j] = CMNe.τ_relax(ice_j, aps, q_icl, ρ)
    end
end

# Colors — distinct, colorblind-friendly
colors = [
    MK.RGBf(0.12, 0.47, 0.71),  # blue
    MK.RGBf(1.0, 0.50, 0.05),   # orange
    MK.RGBf(0.17, 0.63, 0.17),  # green
    MK.RGBf(0.84, 0.15, 0.16),  # red
]

linestyles = [:solid, :dash, :dot, :dashdot]

# Plot
fig = MK.Figure(size = (800, 550))
ax = MK.Axis(
    fig[1, 1],
    xlabel = "Cloud ice specific content q_icl [kg/kg]",
    ylabel = "Relaxation timescale τ [s]",
    title = "Ice relaxation timescale — prescribed N₀\n(ρ = $(ρ) kg/m³)",
    xscale = log10,
    yscale = log10,
)

for (j, (N_0, label)) in enumerate(zip(N_0_values, N_0_labels))
    MK.lines!(
        ax, q_icl_range, τ_data[:, j],
        label = "N₀ = $(label) m⁻³",
        color = colors[j],
        linewidth = 2.5,
        linestyle = linestyles[j],
    )
end

MK.axislegend(ax, position = :rt)
MK.save("tau_relax_prescribed_N0.svg", fig)

println("Plot saved to tau_relax_prescribed_N0.svg")
