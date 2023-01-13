using Test, Plots, CLIMAParameters

include("../src/Nucleation.jl")
using .Nucleation


parameter_file = "temp_params.toml"
local_exp_file = joinpath(@__DIR__, parameter_file)
FT = Float64
toml_dict = CLIMAParameters.create_toml_dict(FT; override_file = local_exp_file)
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
params = CLIMAParameters.get_parameter_values!(toml_dict, param_names)
params = (; params...)

function plot_pure_h2so4_nucleation_rate(
    h2so4_conc,
    nh3_conc,
    negative_ion_conc,
    temp,
    params,
)

    rates = map(
        x -> Nucleation.cloud_h2so4_nucleation_rate(
            x,
            nh3_conc,
            negative_ion_conc,
            temp,
            params,
        )[1],
        h2so4_conc,
    )
    title = "$temp K"
    Plots.plot(
        title = title,
        h2so4_conc,
        rates,
        xaxis = :log,
        yaxis = :log,
        lw = 3,
        ylabel = "Nucleation rate cm⁻³ s⁻¹",
        label = "Parameterization",
    )
    Plots.svg("$(temp)")
    return rates

end

h2so4_conc = [1e5, 3e5, 5e5, 7e5, 9e5, 1e6, 3e6, 5e6, 7e6, 9e6, 1e7, 3e7, 5e7]
nh3_conc = 0
negative_ion_conc = 0
temp = 208
plot_pure_h2so4_nucleation_rate(
    h2so4_conc,
    nh3_conc,
    negative_ion_conc,
    temp,
    params,
)
