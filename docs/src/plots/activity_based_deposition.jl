import CairoMakie as MK

import CloudMicrophysics as CM
import CloudMicrophysics.Common as CO
import CloudMicrophysics.HetIceNucleation as IN
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI

FT = Float32
tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
feldspar = CMP.Feldspar(FT)
ferrihydrite = CMP.Ferrihydrite(FT)
kaolinite = CMP.Kaolinite(FT)

# Initializing
m_to_cm = 100
Δa_w = range(FT(0), stop = FT(0.32), length = 50)    # difference in solution and ice water activity
J_feldspar = @. IN.deposition_J(feldspar, Δa_w)      # J, in SI units (m⁻² s⁻¹)
J_ferrihydrite = @. IN.deposition_J(ferrihydrite, Δa_w)
J_kaolinite = @. IN.deposition_J(kaolinite, Δa_w)

# Alpert et al 2022 Figure 6
# https://doi.org/10.1039/d1ea00077b
# China et al 2017 Supplementary Figure S5
# https://doi.org/10.1002/2016JD025817

# data read from Fig 6 in Alpert et al 2022 and
# China et al 2017 figure S5
# using https://automeris.io/WebPlotDigitizer/

Alpert2022_Feldspar_Δa = [0.019459, 0.256216]
Alpert2022_Feldspar_J = exp10.([1.039735, 4.165563])

Alpert2022_Ferrihydrite_Δa = [0.0989189, 0.336486]
Alpert2022_Ferrihydrite_J = exp10.([1.2781457, 4.21854])

China2017_Δa = [
    0.029478, 0.033209, 0.034328, 0.043657, 0.045149, 0.079851, 0.084701, 0.10896, 0.11082, 0.12687,
    0.13545, 0.1597, 0.1709, 0.20448, 0.2403, 0.25597, 0.27575]
China2017_J =
    exp10.([
        -4.4108, -5.0923, -1.4199, -0.65547, 1.1709, 1.3786, -0.56696, 1.9061, 2.067, 2.2509,
        2.5335, -0.010108, 2.5005, 3.3701, 4.0395, 4.4841, 5.3108,
    ])
China2017_J_perror =
    exp10.([
        -2.0831, -2.7646, 0.94809, 1.6724, 3.5388, 3.7064, 1.7407, 4.2339, 4.3949, 5.2609,
        5.5436, 2.3579, 5.4905, 5.6377, 5.6248, 5.8687, 6.555,
    ])
China2017_J_merror =
    exp10.([
        -5.7353, -6.3766, -2.684, -1.9597, -0.093296, 0.074277, -1.8713, 0.60178, 0.76277, 0.80605,
        1.0888, -1.2743, 1.0155, 1.9253, 2.6348, 3.0794, 3.9463,
    ])

# Create figure and axis
blue, orange, green, pink, cyan, red, yellow = MK.Makie.wong_colors()
fig = MK.Figure()
ax = MK.Axis(fig[1, 1], xlabel = "Δa_w [unitless]", ylabel = "J [cm⁻² s⁻¹]"; yscale = MK.log10)

# Plot CM.jl model lines
MK.lines!(ax, Δa_w, J_feldspar / m_to_cm^2; label = "CM.jl Feldspar", linestyle = :dash, color = cyan)
MK.lines!(ax, Δa_w, J_ferrihydrite / m_to_cm^2; label = "CM.jl Ferrihydrite", linestyle = :dash, color = pink)
MK.lines!(ax, Δa_w, J_kaolinite / m_to_cm^2; label = "CM.jl Kaolinite", linestyle = :dash, color = green)

# Plot Alpert 2022 data
MK.lines!(ax, Alpert2022_Feldspar_Δa, Alpert2022_Feldspar_J, color = blue, label = "Alpert2022 Feldspar")
MK.lines!(ax, Alpert2022_Ferrihydrite_Δa, Alpert2022_Ferrihydrite_J; color = red, label = "Alpert2022 Ferrihydrite")

# Plot China 2017 scatter data
MK.rangebars!(ax, China2017_Δa, China2017_J_merror, China2017_J_perror; whiskerwidth = 5, color = yellow)
MK.scatter!(ax, China2017_Δa, China2017_J;
    color = yellow, markersize = 4, strokewidth = 0.5, label = "China2017 Kaolinite",
)

# Add legend
MK.axislegend(ax, position = :rb)

MK.save("water_activity_depo_nuc.svg", fig)
