using Plots

import CloudMicrophysics.Nucleation as Nucleation

output_diams = 1.7:0.125:10
coag_sinks = 0.5:0.125:10
coag_sink_input_diam = 0.49
input_diam = 1.7
cond_growth_rate = 3
rates = []
for (cs, od) in zip(coag_sinks, output_diams)
    rate = Nucleation.apparent_nucleation_rate(
        od,
        1,
        cond_growth_rate,
        cs,
        coag_sink_input_diam,
        input_diam,
    )
    append!(rates, rate)
end

plot(
    collect(output_diams),
    rates,
    # yaxis = :log,
    xlabel = "Diameter (nm)",
    ylabel = "Apparent Nucleation Rate Reduction",
    legend = false,
)
savefig("apparent_nucleation.svg");
