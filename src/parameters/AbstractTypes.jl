"""
    ParametersType

The top-level super-type for all cloud microphysics parameters
"""
abstract type ParametersType end


# Temporary fallback until we stop checking eltype from the parameters
Base.eltype(p::ParametersType) = fieldcount(typeof(p)) > 0 ? eltype(getfield(p, 1)) : Any
Base.broadcastable(x::ParametersType) = tuple(x)
Base.show(io::IO, x::ParametersType) =
    ShowMethods.parseable_show_with_fields_no_type_header(io, x; with_module_prefix = false)
Base.show(io::IO, mime::MIME"text/plain", x::ParametersType) =
    ShowMethods.show_type_and_fields(io, mime, x)

"""
    AerosolType

The top-level super-type for all aerosol properties
"""
abstract type AerosolType <: ParametersType end

"""
    AerosolDistributionType

The top-level super-type for all aerosol distribution types
"""
abstract type AerosolDistributionType end

"""
    CloudCondensateType

The top-level super-type for cloud condensate types (liquid and ice)
"""
abstract type CloudCondensateType <: ParametersType end

"""
    PrecipitationType

The top-level super-type for precipitation types (rain and snow)
"""
abstract type PrecipitationType <: ParametersType end

"""
    TerminalVelocityType

The top-level super-type for terminal velocity parameterizations
"""
abstract type TerminalVelocityType <: ParametersType end

"""
    Precipitation2MType

The top-level super-type for 2-moment precipitation parameterizations
"""
abstract type Precipitation2MType <: ParametersType end
