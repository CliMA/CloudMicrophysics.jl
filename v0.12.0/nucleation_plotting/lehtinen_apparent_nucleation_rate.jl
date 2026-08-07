import CLIMAParameters as CP
using Plots

include("../../src/Nucleation.jl")
using .Nucleation

FT = Float64
toml_dict = CP.create_toml_dict(FT)
param_names = [
    "u_b_n",
    "v_b_n",
    "w_b_n",
    "u_b_i",
    "v_b_i",
    "w_b_i",
    "u_t_n",
    "v_t_n",
    "w_t_n",
    "u_t_i",
    "v_t_i",
    "w_t_i",
    "p_t_n",
    "p_A_n",
    "a_n",
    "p_t_i",
    "p_A_i",
    "a_i",
    "p_b_n",
    "p_b_i",
    "p_t_n",
    "p_t_i",
]
params = CP.get_parameter_values!(toml_dict, param_names)
params = (; params...)

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
