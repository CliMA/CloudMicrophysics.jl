import CairoMakie as MK
CairoMakie = MK

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMNe
import CloudMicrophysics.ThermodynamicsInterface as TDI

FT = Float64

# Parameters
ice = CMP.CloudIce(FT)
aps = CMP.AirProperties(FT)
frs = CMP.Frostenberg2023(FT)
tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

# Temperature range: -40°C to 20°C
T_freeze = FT(273.15)
T_range = range(T_freeze - 40, T_freeze + 20, length = 500)
T_celsius = T_range .- T_freeze

# q_icl values to plot
q_icl_values = [FT(0), FT(1e-8), FT(1e-6), FT(1e-4), FT(1e-2)]
q_icl_labels = ["0", "10⁻⁸", "10⁻⁶", "10⁻⁴", "10⁻²"]

# Air density (representative value for mid-troposphere)
ρ = FT(0.8)

import ClimaParams as CP

# Sublimation timescale (constant, from ClimaParams defaults)
params_1m = CMP.Microphysics1MParams(CP.create_toml_dict(FT))
τ_sub = params_1m.process_params.cloud_ice_formation.τ_relax

# Compute the effective timescale used in conv_q_vap_to_q_icl:
#   T < T_freeze (deposition):  τ_dep × Γᵢ  (Frostenberg τ with thermodynamic correction)
#   T ≥ T_freeze (sublimation): τ_sub × Γᵢ  (constant τ_sub with thermodynamic correction)
τ_data = Matrix{FT}(undef, length(T_range), length(q_icl_values))
for (j, q_icl) in enumerate(q_icl_values)
    q_tot = FT(0.01)  # representative total water content
    q_lcl = FT(0)
    q_rai = FT(0)
    q_sno = FT(0)
    for (i, T) in enumerate(T_range)
        Rᵥ = TDI.Rᵥ(tps)
        Lₛ = TDI.Lₛ(tps, T)
        cₚ_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)
        qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
        dqsi_dT = CMNe.dqcld_dT(qᵥ_sat_ice, Lₛ, Rᵥ, T)
        Γᵢ = CMNe.gamma_helper(Lₛ, cₚ_air, dqsi_dT)

        if T < T_freeze
            # Deposition branch: τ_dep × Γᵢ
            τ_dep = CMNe.τ_relax(ice, aps, frs, q_icl, T)
            τ_data[i, j] = τ_dep * Γᵢ
        else
            # Sublimation branch: τ_sub × Γᵢ
            τ_data[i, j] = τ_sub * Γᵢ
        end
    end
end

# Replace Inf with NaN for plotting (log scale can't handle Inf)
τ_data[.!isfinite.(τ_data)] .= NaN

# Colors — distinct, colorblind-friendly
colors = [
    MK.RGBf(0.12, 0.47, 0.71),  # blue
    MK.RGBf(1.0, 0.50, 0.05),  # orange
    MK.RGBf(0.17, 0.63, 0.17),  # green
    MK.RGBf(0.84, 0.15, 0.16),  # red
    MK.RGBf(0.58, 0.40, 0.74),  # purple
]

linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]

# Plot
fig = MK.Figure(size = (800, 550))
ax = MK.Axis(
    fig[1, 1],
    xlabel = "Temperature [°C]",
    ylabel = "Effective timescale [s]",
    title = "Ice relaxation timescale — Frostenberg et al. (2023)\nDeposition: τ_dep×Γᵢ (T < 0°C) | Sublimation: τ_sub×Γᵢ (T ≥ 0°C)",
    yscale = log10,
)

for (j, (q_icl, label)) in enumerate(zip(q_icl_values, q_icl_labels))
    MK.lines!(
        ax, T_celsius, τ_data[:, j],
        label = "q_icl = $(label) kg/kg",
        color = colors[j],
        linewidth = 2.5,
        linestyle = linestyles[j],
    )
end

# Add vertical line at T_freeze
MK.vlines!(ax, 0.0, color = :gray, linestyle = :dash, linewidth = 1)
MK.text!(ax, 0.5, 1e10, text = "T_freeze", fontsize = 12, color = :gray)

MK.axislegend(ax, position = :rt)
MK.save("tau_relax_frostenberg.svg", fig)

println("Plot saved to tau_relax_frostenberg.svg")
