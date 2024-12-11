"""
    ParametersType{FT}

The top-level super-type for all cloud microphysics parameters
"""
abstract type ParametersType{FT} end
Base.eltype(::ParametersType{FT}) where {FT} = FT
Base.broadcastable(o::ParametersType) = Ref(o)

"""
    AerosolType{FT}

The top-level super-type for all aerosol properties
"""
abstract type AerosolType{FT} <: ParametersType{FT} end

"""
    AerosolDistributionType

The top-level super-type for all aerosol distribution types
"""
abstract type AerosolDistributionType end

"""
    CloudCondensateType{FT}

The top-level super-type for cloud condensate types (liquid and ice)
"""
abstract type CloudCondensateType{FT} <: ParametersType{FT} end

"""
    PrecipitationType{FT}

The top-level super-type for precipitation types (rain and snow)
"""
abstract type PrecipitationType{FT} <: ParametersType{FT} end

"""
    TerminalVelocityType{FT}

The top-level super-type for terminal velocity parameterizations
"""
abstract type TerminalVelocityType{FT} <: ParametersType{FT} end

"""
    Precipitation2MType

The top-level super-type for 2-moment precipitation parameterizations
"""
abstract type Precipitation2MType{FT} <: ParametersType{FT} end
