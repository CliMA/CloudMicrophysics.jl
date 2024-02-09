export Parameters0M

"""
    Parameters0M{FT}

Parameters for zero-moment bulk microphysics scheme

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Parameters0M{FT} <: ParametersType{FT}
    "precipitation timescale [s]"
    τ_precip::FT
    "specific humidity precipitation threshold [-]"
    qc_0::FT
    "supersaturation precipitation threshold [-]"
    S_0::FT
end

Parameters0M(::Type{FT}) where {FT <: AbstractFloat} =
    Parameters0M(CP.create_toml_dict(FT))

function Parameters0M(td::CP.AbstractTOMLDict)
    name_map = (;
        :precipitation_timescale => :τ_precip,
        :specific_humidity_precipitation_threshold => :qc_0,
        :supersaturation_precipitation_threshold => :S_0,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return Parameters0M{FT}(; parameters...)
end
