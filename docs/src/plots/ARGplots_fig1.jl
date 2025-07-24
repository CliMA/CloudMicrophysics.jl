import CairoMakie as MK

import CloudMicrophysics as CM
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ClimaParams

FT = Float64
tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
aip = CMP.AirProperties(FT)
ap = CMP.AerosolActivationParameters(FT)

# Atmospheric conditions
T = 294.0         # air temperature
p = 1000.0 * 1e2   # air pressure
w = 0.5           # vertical velocity

# We need the phase partition here only so that we can compute the
# moist R_m and cp_m in aerosol activation module.
# We are assuming here saturated conditions and no liquid water or ice.
# This is consistent with the assumptions of the aerosol activation scheme.
p_vs = TDI.saturation_vapor_pressure_over_liquid(tps, T)
q_vs = 1 / (1 - 1 / TDI.Rd_over_Rv(tps) * (p_vs - p) / p_vs)

# Abdul-Razzak and Ghan 2000 Figure 1 mode 1
# https://doi.org/10.1029/1999JD901161
r_dry = 0.05 * 1e-6 # um
stdev = 2.0         # -
N_1 = 100.0 * 1e6   # 1/m3

# Sulfate - universal parameters
sulfate = CMP.Sulfate(FT)

mass_fractions_1 = (1.0,)
paper_mode_1_B = AM.Mode_B(
    r_dry, stdev, N_1, mass_fractions_1,
    (sulfate.ϵ,), (sulfate.ϕ,), (sulfate.M,), (sulfate.ν,), (sulfate.ρ,),
)

N_2_range = range(0, stop = 5000 * 1e6, length = 100)
N_act_frac_B = similar(N_2_range)

for (it, N_2) in enumerate(N_2_range)
    mass_fractions_2 = (1.0,)
    paper_mode_2_B = AM.Mode_B(
        r_dry, stdev, N_2, mass_fractions_2,
        (sulfate.ϵ,), (sulfate.ϕ,), (sulfate.M,), (sulfate.ν,), (sulfate.ρ,),
    )
    AD_B = AM.AerosolDistribution((paper_mode_1_B, paper_mode_2_B))
    N_act_frac_B[it] = AA.N_activated_per_mode(ap, AD_B, aip, tps, T, p, w, q_vs, FT(0), FT(0))[1] / N_1
end

# data read from Fig 1 in Abdul-Razzak and Ghan 2000
# using https://automeris.io/WebPlotDigitizer/
include(joinpath(pkgdir(CM), "docs", "src", "plots", "ARGdata.jl"))

# Create figure and axis
m_to_cm = FT(100)
fig = MK.with_theme(MK.theme_minimal(), linewidth = 2, fontsize = 14) do
    fig = MK.Figure()
    ax = MK.Axis(fig[1, 1];
        xlabel = "Mode 2 aerosol number concentration [1/cm³]", ylabel = "Mode 1 number fraction activated",
        limits = (extrema(N_2_range) ./ m_to_cm^3, (0, 1)),
    )

    # Plot the computed data
    MK.lines!(ax, N_2_range / m_to_cm^3, N_act_frac_B; label = "CliMA-B")
    # Plot the observed data points
    MK.scatter!(ax, Fig1_x_obs, Fig1_y_obs; label = "paper observations", color = :black, markersize = 15)
    # Plot the parameterization line
    MK.lines!(ax, Fig1_x_param, Fig1_y_param; label = "paper parameterization", color = :black)
    # Add legend
    MK.axislegend(ax, position = :rt)
    # Save the figure
    MK.save("Abdul-Razzak_and_Ghan_fig_1.svg", fig)

    fig
end
