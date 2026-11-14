import Plots as PL

import Thermodynamics as TD
import CloudMicrophysics as CM
import CloudMicrophysics.Common as CO
import CloudMicrophysics.HetIceNucleation as IN
import CloudMicrophysics.Parameters as CMP

FT = Float32
const tps = TD.Parameters.ThermodynamicsParameters(FT)
const feldspar = CMP.Feldspar(FT)
const ferrihydrite = CMP.Ferrihydrite(FT)
const kaolinite = CMP.Kaolinite(FT)

# Initializing
Δa_w = range(FT(0), stop = FT(0.32), length = 50)    # difference in solution and ice water activity
J_feldspar = @. IN.deposition_J(feldspar, Δa_w)      # J in SI units
log10J_converted_feldspar = @. log10(J_feldspar * 1e-4)       # converts J into cm^-2 s^-1 and takes log
J_ferrihydrite = @. IN.deposition_J(ferrihydrite, Δa_w)
log10J_converted_ferrihydrite = @. log10(J_ferrihydrite * 1e-4)
J_kaolinite = @. IN.deposition_J(kaolinite, Δa_w)
log10J_converted_kaolinite = @. log10(J_kaolinite * 1e-4)

# Alpert et al 2022 Figure 6
# https://doi.org/10.1039/d1ea00077b
# China et al 2017 Supplementary Figure S5
# https://doi.org/10.1002/2016JD025817

# data read from Fig 6 in Alpert et al 2022 and
# China et al 2017 figure S5
# using https://automeris.io/WebPlotDigitizer/

#! format: off
Alpert2022_Feldspar_Delta_a = [0.019459, 0.256216]
Alpert2022_Feldspar_log10J = [1.039735, 4.165563]

Alpert2022_Ferrihydrite_Delta_a = [0.0989189, 0.336486]
Alpert2022_Ferrihydrite_log10J = [1.2781457, 4.21854]

China2017_Delta_a = [0.01918, 0.02398, 0.02518, 0.03537, 0.07314, 0.07794, 0.10252, 0.10492, 0.1217, 0.1307, 0.15588, 0.16787, 0.20264, 0.23981, 0.256, 0.27638]
China2017_log10J = [-4.36923, -5.07692, -1.38462, -0.64615, 1.2, 1.35385, -0.58462, 1.90769, 2.06154, 2.24615, 2.52308, 0, 2.46154, 3.32308, 4, 4.43077, 5.26154]
#! format: on

PL.plot(
    Δa_w,
    log10J_converted_feldspar,
    label = "CM.jl Feldspar",
    xlabel = "Delta a_w [unitless]",
    ylabel = "log10(J) [cm^-2 s^-1]",
    linestyle = :dash,
    linecolor = :cyan,
)
PL.plot!(
    Δa_w,
    log10J_converted_ferrihydrite,
    label = "CM.jl Ferrihydrite",
    linestyle = :dash,
    linecolor = :pink,
)
PL.plot!(
    Δa_w,
    log10J_converted_kaolinite,
    label = "CM.jl Kaolinite",
    linestyle = :dash,
    linecolor = :orange,
)
PL.plot!(
    Alpert2022_Feldspar_Delta_a,
    Alpert2022_Feldspar_log10J,
    linecolor = :blue,
    label = "Alpert2022 Feldspar",
)
PL.plot!(
    Alpert2022_Ferrihydrite_Delta_a,
    Alpert2022_Ferrihydrite_log10J,
    linecolor = :red,
    label = "Alpert2022 Ferrihydrite",
)
PL.scatter!(
    China2017_Delta_a,
    China2017_log10J,
    markercolor = :yellow,
    markersize = 2,
    label = "China2017 Kaolinite",
)

PL.savefig("water_activity_depo_nuc.svg")
