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
