import Plots as PL

import CloudMicrophysics as CM
import CloudMicrophysics.HetIceNucleation as IN
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP

FT = Float32
ip = CMP.Frostenberg2023(FT)

T_range_celsius = range(-40, stop = -2, length = 500)  # air temperature [°C]
T_range_kelvin = T_range_celsius .+ 273.15  # air temperature [K]
INPC_range = 10.0 .^ (range(-5, stop = 7, length = 500)) #ice nucleating particle concentration

frequency = [
    IN.INP_concentration_frequency(ip, INPC, T) > 0.015 ?
    IN.INP_concentration_frequency(ip, INPC, T) : missing for
    INPC in INPC_range, T in T_range_kelvin
]
mu = [exp(log(-(ip.b * T)^9 * 10^(-9))) for T in T_range_celsius]


PL.contourf(
    T_range_kelvin,
    INPC_range,
    frequency,
    xlabel = "T [K]",
    ylabel = "INPC [m⁻³]",
    colorbar_title = "frequency",
    yaxis = :log,
    color = :lighttest,
    gridlinewidth = 3,
    ylims = (1e-5, 1e7),
    lw = 0,
)

PL.plot!(
    T_range_kelvin,
    mu,
    label = "median INPC",
    legend = :topright,
    color = :darkred,
)

PL.plot!(
    repeat([257], 50),
    10.0 .^ (range(-2, stop = 6, length = 50)),
    label = "T = 257 K",
    color = :darkgreen,
    linestyle = :dash,
)
PL.savefig("Frostenberg_fig1.svg")


#plotting the distribution for T=257 K

T = FT(257.0) # [K]
INPC_range = FT(10) .^ range(FT(-1), 4, length = 100)
frequency = [IN.INP_concentration_frequency(ip, INPC, T) for INPC in INPC_range]

PL.plot(
    INPC_range,
    frequency,
    xlabel = "INPC [m⁻³]",
    ylabel = "frequency",
    xaxis = :log,
    color = :darkgreen,
    label = "T = 257 K",
    legend = :topright,
)

PL.savefig("Frostenberg_fig1_T16.svg")
