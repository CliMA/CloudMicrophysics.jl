export Parameters0M

"""
    Parameters0M{FT}

Parameters for zero-moment bulk microphysics scheme

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Parameters0M{FT} <: ParametersType
    "precipitation timescale [s]"
    τ_precip::FT
    "condensate specific content precipitation threshold [-]"
    qc_0::FT
    "supersaturation precipitation threshold [-]"
    S_0::FT
end

function Parameters0M(td::CP.ParamDict)
    name_map = (;
        :precipitation_timescale => :τ_precip,
        :specific_humidity_precipitation_threshold => :qc_0,
        :supersaturation_precipitation_threshold => :S_0,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Parameters0M(; parameters...)
end
