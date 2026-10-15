export Parameters0M, Evaporation0M

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

"""
    Evaporation0M{FT}

Parameters for the evaporation (and sublimation) of the zero-moment
precipitation flux in sub-saturated air, following
[Kessler1969](@cite) as used in the IFS ([Tiedtke1993](@cite), [IFS2023](@cite)).

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Evaporation0M{FT} <: ParametersType
    "evaporation coefficient (`0` disables evaporation) [1/s]"
    k_E::FT
    "critical relative humidity above which precipitation does not evaporate [-]"
    RH_c::FT
    "fraction of the grid box covered by precipitation [-]"
    a_p::FT
    "precipitation flux scale [kg/m2/s]"
    α₂::FT
    "precipitation flux exponent [-]"
    α₃::FT
end

function Evaporation0M(td::CP.ParamDict)
    name_map = (;
        :precipitation_evaporation_coefficient => :k_E,
        :precipitation_evaporation_rh_crit => :RH_c,
        :precipitation_evaporation_area_fraction => :a_p,
        :precipitation_evaporation_flux_scale => :α₂,
        :precipitation_evaporation_exponent => :α₃,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Evaporation0M(; parameters...)
end
