import CloudMicrophysics as CM
import CLIMAParameters as CP
import Thermodynamics as TD

"""
    cloud_microphysics_parameters(toml_dict)

Construct a CloudMicrophysics parameter struct, using
a `toml_dict`, from the CliMA's centralized parameter file.

```julia
include(joinpath(pkgdir(CloudMicrophysics), "test", "create_parameters.jl"))
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
param_set = cloud_microphysics_parameters(toml_dict)
```
"""
function cloud_microphysics_parameters(toml_dict)
    FT = CP.float_type(toml_dict)
    aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
    thermo_params = TD.Parameters.ThermodynamicsParameters{FT}(; pairs...)
    TP = typeof(thermo_params)

    aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
    aliases = setdiff(aliases, ["thermo_params"])
    pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
    return CM.Parameters.CloudMicrophysicsParameters{FT, TP}(;
        pairs...,
        thermo_params,
    )
end

function nucleation_parameters(toml_dict)
    nucleation_param_names = (
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
        "a_1",
        "a_2",
        "a_3",
        "a_4",
        "a_5",
        "k_H2SO4org",
        "Y_MTO3",
        "Y_MTOH",
        "k_MTO3",
        "k_MTOH",
        "exp_MTO3",
        "exp_MTOH",
    )
    params = CP.get_parameter_values!(toml_dict, nucleation_param_names)
    return (; params...)
end
