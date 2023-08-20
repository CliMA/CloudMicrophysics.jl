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
export AbstractAerosolActivation
export AbstractParameterizedAerosolActivation
export ARG2000Type
export Abstract2MPrecipType
export KK2000Type
export B1994Type
export TC1980Type
export LD2004Type
export SB2006Type
export Chen2022Type
export AbstractAerosolType
export ArizonaTestDustType
export DesertDustType

"""
    AbstractAerosolDistribution

The top-level super-type for all aerosol distribution types.
"""
abstract type AbstractAerosolDistribution{T} end

"""
    AbstractAerosolActivation

The top-level super-type for all aerosol activation schemes.
"""
abstract type AbstractAerosolActivation end

"""
    AbstractParameterizedAerosolActivation

The type for all aerosol activation schemes that use a parameterization based
    on maximum supersaturation (as opposed to using an ML emulator)
"""
abstract type AbstractParameterizedAerosolActivation end

"""
    ARG2000Type

The type for the aerosol activation scheme formulated by Abdul-Razzak and Ghan (2000)
"""
struct ARG2000Type <: AbstractParameterizedAerosolActivation end

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

"""
    Abstract2MPrecipType

The top-level super-type for 2-moment precipitation formation types
"""
abstract type Abstract2MPrecipType end

"""
    KK2000Type

The type for 2-moment precipitation formation by Khairoutdinov and Kogan (2000)
"""
struct KK2000Type <: Abstract2MPrecipType end

"""
    B1994Type

The type for 2-moment precipitation formation by Beheng (1994)
"""
struct B1994Type <: Abstract2MPrecipType end

"""
    TC1980Type

The type for 2-moment precipitation formation by Tripoli and Cotton (1980)
"""
struct TC1980Type <: Abstract2MPrecipType end

"""
    LD2004Type

The type for 2-moment precipitation formation by Liu and Daum (2004)
"""
struct LD2004Type <: Abstract2MPrecipType end

"""
    SB2006Type

The type for 2-moment precipitation formation by Seifert and Beheng (2006)
"""
struct SB2006Type <: Abstract2MPrecipType end

"""
    Chen2022Type

The type for 2-moment precipitation terminal velocity by Chen et al 2022
"""
struct Chen2022Type <: Abstract2MPrecipType end

"""
    AbstractAerosolType

The top-level super-type for all aerosol types.
"""
abstract type AbstractAerosolType end

"""
    ArizonaTestDustType

The type for Arizona test dust for deposition activated fraction
"""
struct ArizonaTestDustType <: AbstractAerosolType end

"""
    DesertDustType

The type for desert dust for deposition activated fraction
"""
struct DesertDustType <: AbstractAerosolType end

"""
    KaoliniteType

The type for Kaolinite for ABIFM nucleation rate
"""
struct KaoliniteType <: AbstractAerosolType end

"""
    IlliteType

The type for illite for ABIFM nucleation rate
"""
struct IlliteType <: AbstractAerosolType end

Base.broadcastable(x::LiquidType) = tuple(x)
Base.broadcastable(x::IceType) = tuple(x)
Base.broadcastable(x::RainType) = tuple(x)
Base.broadcastable(x::SnowType) = tuple(x)

end #module CommoniTypes.jl
