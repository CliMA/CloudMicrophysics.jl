"""
    A container for information on aerosol size distribution
    and chemical properties.

    The size distribution is a sum of lognormal internally mixed modes.
    The chemical composition can be expressed using kappa parameter
    or hygroscopicity parameter B.
"""
module AerosolModel

export Mode_B
export Mode_κ

export AbstractAerosolDistribution

"""
    AbstractAerosolDistribution

The top-level super-type for all aerosol distribution types.
"""
abstract type AbstractAerosolDistribution{T} end

"""
    Mode_B

Represents the sizes and chemical composition
of aerosol particles in one size distribution mode.
The mode is assumed to be made up of internally mixed components
and follow a lognormal size distribution.
The chemical composition of aerosol particles in this mode
is described using the parameters from Abdul-Razzak and Ghan 2000.
"""
struct Mode_B{T, FT}
    "geometric mean dry radius"
    r_dry::FT
    "geometric standard deviation"
    stdev::FT
    "total number concentration"
    N::FT
    "tuple of mass mixing ratios for all components in this mode"
    mass_mix_ratio::T
    "tuple of mass fractions of soluble material for all components in this mode"
    soluble_mass_frac::T
    "tuple of osmotic coefficients for all components in this mode"
    osmotic_coeff::T
    "tuple of molar masses for all components in this mode"
    molar_mass::T
    "tuple of number of ions the salt dissociates into for all components in this mode"
    dissoc::T
    "tuple of aerosol densities for all components in this mode"
    aerosol_density::T
    "number of components in the mode"
    n_components::Int64
end

"""
    Mode_κ

Represents the sizes and chemical composition
of aerosol particles in one size distribution mode.
The mode is assumed to be made up of internally mixed components
and follow a lognormal size distribution.
The chemical composition of aerosol particles in this mode
is described using the parameters from Petters and Kreidenweis 2007.
"""
struct Mode_κ{T, FT}
    "geometric mean dry radius"
    r_dry::FT
    "geometric standard deviation"
    stdev::FT
    "total number concentration"
    N::FT
    "tuple of volume mixing ratios for all components in this mode"
    vol_mix_ratio::T
    "tuple of mass mixing ratios for all components in this mode"
    mass_mix_ratio::T
    "tuple of molar masses for all components in this mode"
    molar_mass::T
    "tuple of kappa-kohler values for all components in this mode"
    kappa::T
    "number of components in the mode"
    n_components::Int64
end

"""
    AerosolDistribution

Represents the aerosol size distribution as a tuple with different modes.
All modes have to either be of type Mode_B (Abdul-Razzak and Ghan 2000)
or of type Mode_κ (Petters and Kreidenweis 2007).

# Constructors

    AerosolDistribution(Modes::T)
"""
struct AerosolDistribution{T} <: AbstractAerosolDistribution{T}

    "tuple with all aerosol size distribution modes"
    Modes::T

    function AerosolDistribution(Modes::NTuple{N, T}) where {N, T}
        return new{typeof(Modes)}(Modes)
    end

end

end
