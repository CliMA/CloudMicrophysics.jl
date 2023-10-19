import Thermodynamics as TD
const TP = TD.Parameters

import CLIMAParameters as CP

"""
    ThermodynamicsParameters

A wrapper function for creating Thermodynamics.jl parameters
"""
function ThermodynamicsParameters(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(
        FT;
        dict_type = "alias",
    ),
) where {FT}
    aliases = string.(fieldnames(TP.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")

    return TP.ThermodynamicsParameters{FT}(; param_pairs...)
end
