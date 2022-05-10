"""

Additional type hierarchy to dispatch over for some microphysics parameters

"""
module CommonTypes

export AbstractCloudType
export AbstractPrecipType
export LiquidType
export IceType
export RainType
export SnowType
export AbstractAerosolDistribution

"""
    AbstractAerosolDistribution

The top-level super-type for all aerosol distribution types.
"""
abstract type AbstractAerosolDistribution{T} end

"""
    AbstractCloudType

The top-level super-type for cloud liquid water and cloud ice types
"""
abstract type AbstractCloudType end

"""
    AbstractPrecipType

The top-level super-type for precipitation types (rain and snow)
"""
abstract type AbstractPrecipType end

"""
    LiquidType

The type for cloud liquid water condensate
"""
struct LiquidType <: AbstractCloudType end

"""
    IceType

The type for cloud ice condensate
"""
struct IceType <: AbstractCloudType end

"""
    RainType

The type for rain
"""
struct RainType <: AbstractPrecipType end

"""
    SnowType

The type for snow
"""
struct SnowType <: AbstractPrecipType end

end #module CommoniTypes.jl
