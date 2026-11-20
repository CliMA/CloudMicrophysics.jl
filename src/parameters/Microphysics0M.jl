export Parameters0M, Evaporation0M

"""
    Parameters0M{FT}

Parameters for zero-moment bulk microphysics scheme

# Fields
$(DocStringExtensions.FIELDS)
"""
@params struct Parameters0M{FT} <: ParametersType
    "precipitation timescale [s]"
    τ_precip::FT
    "condensate specific content precipitation threshold [-]"
    qc_0::FT
    "supersaturation precipitation threshold [-]"
    S_0::FT
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
