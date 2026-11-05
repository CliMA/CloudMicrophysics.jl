import CairoMakie as MK

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

const CMO = CM.Common
const CMI = CM.HomIceNucleation
const CMP = CM.Parameters

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))
FT = Float64
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
const prs = cloud_microphysics_parameters(toml_dict)
thermo_params = CMP.thermodynamics_params(prs)

# Initializing
T_range = range(229.0, stop = 234.5, length = 50)  # air temperature
x_sulph = Vector{FT}([0.03, 0.04, 0.06])           # wt% sulphuric acid in droplets

# Solving for Δa and J values
Δa1 = [CMO.Delta_a_w(prs, x_sulph[1], T) for T in T_range]
Δa2 = [CMO.Delta_a_w(prs, x_sulph[2], T) for T in T_range]
Δa3 = [CMO.Delta_a_w(prs, x_sulph[3], T) for T in T_range]

J1 = @. CMI.homogeneous_J(prs, Δa1)
J2 = @. CMI.homogeneous_J(prs, Δa2)
J3 = @. CMI.homogeneous_J(prs, Δa3)

log10J_1 = @. log10(J1)
log10J_2 = @. log10(J2)
log10J_3 = @. log10(J3)

Δa_range = range(0.27, stop = 0.32, length = 50)
J_given_Δa = @. CMI.homogeneous_J(prs, Δa_range)

#! format: off
# Spichtinger et al 2023
# https://acp.copernicus.org/articles/23/2035/2023/acp-23-2035-2023.pdf
# Values are from Koop et al. (2000) line in Figure 1
Spichtinger2023_temp = [
    230.0139, 231.0179, 231.9860, 233.00797, 234.01195, 235.01594, 236.00199, 237.16733,
]
Spichtinger2023_log10J = [
    24.33735, 22.6506, 21.0843, 19.5181, 18.072289, 16.626506, 15.24096, 13.493975,
]

# Baumgartner et al 2022
# https://acp.copernicus.org/articles/22/65/2022/acp-22-65-2022.pdf
# Figure 2a
Baum_Delta_a = [0.26, 0.27, 0.28, 0.29, 0.3, 0.32, 0.33, 0.339]
Baum_J = [4.25e-5, 0.306, 454.09, 2.06e5, 6.31e7, 1.8e12, 4.45e14, 1.69e17]

# Plotting
fig = MK.Figure(resolution = (800, 500))
ax1 = MK.Axis(
    fig[1, 1],
    ylabel = "log10(J) with J in SI units",
    xlabel = "Temperature [K]",
    title = "CliMA vs Spichtinger2023",
    xticklabelsize = 14.0f0,
    xlabelsize = 14,
    ylabelsize = 14,
    limits = ((228.0, 240.0), nothing),
)
ax2 = MK.Axis(
    fig[1, 2],
    xlabel = "Δa_w [-]",
    ylabel = "J [cm^-3 s^-1]",
    title = "CliMA vs Baumgartner2022",
    yscale = log10,
    xticklabelsize = 14.0f0,
    xlabelsize = 14,
    ylabelsize = 14,
)

MK.lines!(ax1, Spichtinger2023_temp, Spichtinger2023_log10J, label = "Spichtinger 2023 x = ?")
MK.lines!(ax1, T_range, log10J_1, label = "CliMA x = 0.03")
MK.lines!(ax1, T_range, log10J_2, label = "CliMA x = 0.04")
MK.lines!(ax1, T_range, log10J_3, label = "CliMA x = 0.06")

MK.lines!(ax2, Baum_Delta_a, Baum_J, label = "Baumgartner 2022")
MK.lines!(ax2, Δa_range, J_given_Δa .* 1e-6, label = "CliMA")
#! format: on

MK.axislegend(ax1, position = :rt, labelsize = 13.0f0)
MK.axislegend(ax2, position = :rb, labelsize = 14.0f0)

MK.save("HomFreezingPlots.svg", fig)
