import CairoMakie as MK

import CloudMicrophysics as CM
import CloudMicrophysics.HetIceNucleation as IN
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP

FT = Float32
ip = CMP.Frostenberg2023(FT)

T_range = range(233, stop = 271, length = 500) # air temperature [K]
INPC_range = 10.0 .^ (range(-5, stop = 7, length = 500)) #ice nucleating particle concentration

frequency = [IN.INP_concentration_frequency(ip, INPC, T) for T in T_range, INPC in INPC_range]
mu = [exp(IN.INP_concentration_mean(ip, T)) for T in T_range]

fig = MK.Figure(size = (700, 500))
ax = MK.Axis(fig[1, 1]; xlabel = "T [K]", ylabel = "INPC [m⁻³]", yscale = log10)
MK.limits!(ax, extrema(T_range)..., 1e-5, 1e7)
# frequencies below 0.015 are left blank; the filled contours are rasterized,
# since their many polygons would otherwise make a large SVG
levels = range(0.015, maximum(frequency), length = 16)
cf = MK.contourf!(ax, T_range, INPC_range, frequency; levels, colormap = :lighttest, rasterize = 1)
MK.lines!(ax, T_range, mu; color = :darkred, label = "median INPC")
MK.lines!(
    ax,
    fill(257, 50),
    10.0 .^ (range(-2, stop = 6, length = 50));
    color = :darkgreen,
    linestyle = :dash,
    label = "T = 257 K",
)
MK.axislegend(ax; position = :rt)
MK.Colorbar(fig[1, 2], cf; label = "frequency")
MK.save("Frostenberg_fig1.svg", fig)

#plotting the distribution for T=257 K

T = FT(257.0) # [K]
INPC_range = FT(10) .^ range(FT(-1), 4, length = 100)
frequency = [IN.INP_concentration_frequency(ip, INPC, T) for INPC in INPC_range]

fig = MK.Figure(size = (700, 500))
ax = MK.Axis(fig[1, 1]; xlabel = "INPC [m⁻³]", ylabel = "frequency", xscale = log10)
MK.lines!(ax, INPC_range, frequency; color = :darkgreen, label = "T = 257 K")
MK.axislegend(ax; position = :rt)
MK.save("Frostenberg_fig1_T16.svg", fig)
nothing
