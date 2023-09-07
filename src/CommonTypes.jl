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
export Abstract2MPrecipType
export KK2000Type
export B1994Type
export TC1980Type
export LD2004Type
export SB2006Type
export AbstractTerminalVelocityType
export Blk1MVelType
export SB2006VelType
export Chen2022Type
export AbstractAerosolType
export ArizonaTestDustType
export DesertDustType
export KaoliniteType
export IlliteType

"""
    AbstractAerosolDistribution

The top-level super-type for all aerosol distribution types.
"""
abstract type AbstractAerosolDistribution{T} end
Base.broadcastable(x::AbstractAerosolDistribution) = tuple(x)

"""
    AbstractCloudType

The top-level super-type for cloud liquid water and cloud ice types
"""
abstract type AbstractCloudType end
Base.broadcastable(x::AbstractCloudType) = tuple(x)

"""
    AbstractPrecipType

The top-level super-type for precipitation types (rain and snow)
"""
abstract type AbstractPrecipType end
Base.broadcastable(x::AbstractPrecipType) = tuple(x)

"""
    LiquidType

The type for cloud liquid water condensate
"""
struct LiquidType <: AbstractCloudType end
Base.broadcastable(x::LiquidType) = tuple(x)

"""
    IceType

The type for cloud ice condensate
"""
struct IceType <: AbstractCloudType end
Base.broadcastable(x::IceType) = tuple(x)

"""
    RainType

The type for rain
"""
struct RainType <: AbstractPrecipType end
Base.broadcastable(x::RainType) = tuple(x)

"""
    SnowType

The type for snow
"""
struct SnowType <: AbstractPrecipType end
Base.broadcastable(x::SnowType) = tuple(x)

"""
    Abstract2MPrecipType

The top-level super-type for 2-moment precipitation formation types
"""
abstract type Abstract2MPrecipType end
Base.broadcastable(x::Abstract2MPrecipType) = tuple(x)

"""
    KK2000Type

The type for 2-moment precipitation formation by Khairoutdinov and Kogan (2000)
"""
struct KK2000Type <: Abstract2MPrecipType end
Base.broadcastable(x::KK2000Type) = tuple(x)

"""
    B1994Type

The type for 2-moment precipitation formation by Beheng (1994)
"""
struct B1994Type <: Abstract2MPrecipType end
Base.broadcastable(x::B1994Type) = tuple(x)

"""
    TC1980Type

The type for 2-moment precipitation formation by Tripoli and Cotton (1980)
"""
struct TC1980Type <: Abstract2MPrecipType end
Base.broadcastable(x::TC1980Type) = tuple(x)

"""
    LD2004Type

The type for 2-moment precipitation formation by Liu and Daum (2004)
"""
struct LD2004Type <: Abstract2MPrecipType end
Base.broadcastable(x::LD2004Type) = tuple(x)

"""
    VarTimeScaleAcnvType

The type for 2-moment precipitation formation based on the 1-moment parameterization
"""
struct VarTimeScaleAcnvType <: Abstract2MPrecipType end
Base.broadcastable(x::VarTimeScaleAcnvType) = tuple(x)

"""
    SB2006Type

The type for 2-moment precipitation formation by Seifert and Beheng (2006)
"""
struct SB2006Type <: Abstract2MPrecipType end
Base.broadcastable(x::SB2006Type) = tuple(x)

"""
    AbstractTerminalVelocityType

The top-level super-type for terminal velocity parameterizations
"""
abstract type AbstractTerminalVelocityType end
Base.broadcastable(x::AbstractTerminalVelocityType) = tuple(x)

"""
    Blk1MVelType

The type for precipitation terminal velocity from the simple 1-moment scheme
"""
struct Blk1MVelType <: AbstractTerminalVelocityType end
Base.broadcastable(x::Blk1MVelType) = tuple(x)

"""
    SB2006VelType

The type for precipitation terminal velocity from Seifert and Beheng (2006)
"""
struct SB2006VelType <: AbstractTerminalVelocityType end
Base.broadcastable(x::SB2006VelType) = tuple(x)

"""
    Chen2022Type

The type for precipitation terminal velocity from Chen et. al. 2022
"""
struct Chen2022Type <: AbstractTerminalVelocityType end
Base.broadcastable(x::Chen2022Type) = tuple(x)

"""
    AbstractAerosolType

The top-level super-type for all aerosol types.
"""
abstract type AbstractAerosolType end
Base.broadcastable(x::AbstractAerosolType) = tuple(x)

"""
    ArizonaTestDustType

The type for Arizona test dust for deposition activated fraction
"""
struct ArizonaTestDustType <: AbstractAerosolType end
Base.broadcastable(x::ArizonaTestDustType) = tuple(x)

"""
    DesertDustType

The type for desert dust for deposition activated fraction
"""
struct DesertDustType <: AbstractAerosolType end
Base.broadcastable(x::DesertDustType) = tuple(x)

"""
    KaoliniteType

The type for Kaolinite for ABIFM nucleation rate
"""
struct KaoliniteType <: AbstractAerosolType end
Base.broadcastable(x::KaoliniteType) = tuple(x)

"""
    IlliteType

The type for illite for ABIFM nucleation rate
"""
struct IlliteType <: AbstractAerosolType end
Base.broadcastable(x::IlliteType) = tuple(x)

end #module CommoniTypes.jl
