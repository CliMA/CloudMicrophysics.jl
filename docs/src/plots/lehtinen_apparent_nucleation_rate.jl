import CairoMakie as MK

import CloudMicrophysics.Nucleation as Nucleation

output_diams = 1.7:0.125:10
coag_sinks = 0.5:0.125:10
coag_sink_input_diam = 0.49
input_diam = 1.7
cond_growth_rate = 3
rates = map(zip(coag_sinks, output_diams)) do (cs, od)
    Nucleation.apparent_nucleation_rate(od, 1, cond_growth_rate, cs, coag_sink_input_diam, input_diam)
end

fig = MK.Figure(size = (700, 500))
ax = MK.Axis(fig[1, 1]; xlabel = "Diameter (nm)", ylabel = "Apparent Nucleation Rate Reduction")
MK.lines!(ax, collect(output_diams), rates; linewidth = 2)
MK.save("apparent_nucleation.svg", fig)
nothing
