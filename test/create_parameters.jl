import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
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

    aliases = string.(fieldnames(CMP.ModalNucleationParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
    modal_nucleation_params = CMP.ModalNucleationParameters{FT}(; pairs...)
    MNP = typeof(modal_nucleation_params)

    aliases = string.(fieldnames(CMP.CloudMicrophysicsParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")

    return CMP.CloudMicrophysicsParameters{FT, TP, MNP}(;
        pairs...,
        thermo_params,
        modal_nucleation_params,
    )
end
