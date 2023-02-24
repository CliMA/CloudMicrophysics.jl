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

export AbstractAerosolType
export ArizonaTestDustType
export DesertDustType

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

end #module CommoniTypes.jl
