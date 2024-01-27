"""
    A container for information on aerosol size distribution
    and chemical properties.

    The size distribution is a sum of lognormal internally mixed modes.
    The chemical composition can be expressed using kappa parameter
    or hygroscopicity parameter B.
"""
module AerosolModel

import ..Parameters as CMP

export Mode_B
export Mode_κ

"""
    Mode_B

Represents the sizes and chemical composition
of aerosol particles in one size distribution mode.
The mode is assumed to be made up of internally mixed components
and follow a lognormal size distribution.
The chemical composition of aerosol particles in this mode
is described using the parameters from Abdul-Razzak and Ghan 2000.
"""
struct Mode_B{NCOMP, T, FT}
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
end

function Mode_B(
    r_dry::FT,
    stdev::FT,
    N::FT,
    mass_mix_ratio::T,
    soluble_mass_frac::T,
    osmotic_coeff::T,
    molar_mass::T,
    dissoc::T,
    aerosol_density::T,
    NCOMP::Int,
) where {T, FT}
    return Mode_B{NCOMP, T, FT}(
        r_dry,
        stdev,
        N,
        mass_mix_ratio,
        soluble_mass_frac,
        osmotic_coeff,
        molar_mass,
        dissoc,
        aerosol_density,
    )
end

""" number of components in the mode """
n_components(::Mode_B{NCOMP}) where {NCOMP} = NCOMP

"""
    Mode_κ

Represents the sizes and chemical composition
of aerosol particles in one size distribution mode.
The mode is assumed to be made up of internally mixed components
and follow a lognormal size distribution.
The chemical composition of aerosol particles in this mode
is described using the parameters from Petters and Kreidenweis 2007.
"""
struct Mode_κ{NCOMP, T, FT}
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
end

function Mode_κ(
    r_dry::FT,
    stdev::FT,
    N::FT,
    vol_mix_ratio::T,
    mass_mix_ratio::T,
    molar_mass::T,
    kappa::T,
    NCOMP::Int,
) where {T, FT}
    return Mode_κ{NCOMP, T, FT}(
        r_dry,
        stdev,
        N,
        vol_mix_ratio,
        mass_mix_ratio,
        molar_mass,
        kappa,
    )
end

""" number of components in the mode """
n_components(::Mode_κ{NCOMP}) where {NCOMP} = NCOMP

"""
    Mode_INP

Represents the sizes and types of particles in one size distribution mode.
The mode is assumed to be made up of internally mixed components
and follow a lognormal size distribution.
"""
struct Mode_INP{NCOMP, FT}
    "geometric mean dry radius"
    r_dry::FT
    "geometric standard deviation"
    stdev::FT
    "total number concentration"
    N::FT
    "aerosol type"
    aerosol_type::CMP.AerosolType{FT}
end

function Mode_INP(
    r_dry::FT,
    stdev::FT,
    N::FT,
    aerosol_type::CMP.AerosolType{FT},
    NCOMP::Int,
) where {FT}
    return Mode_INP{NCOMP, FT}(
        r_dry,
        stdev,
        N,
        aerosol_type::CMP.AerosolType{FT},
    )
end

""" number of components in the mode """
n_components(::Mode_INP{NCOMP}) where {NCOMP} = NCOMP

"""
    AerosolDistribution

Represents the aerosol size distribution as a tuple with different modes.
All modes have to either be of type Mode_B (Abdul-Razzak and Ghan 2000)
or of type Mode_κ (Petters and Kreidenweis 2007).

# Constructors

    AerosolDistribution(Modes::T)
"""
struct AerosolDistribution{T} <: CMP.AerosolDistributionType

    "tuple with all aerosol size distribution modes"
    Modes::T

end
function AerosolDistribution(Modes::NTuple{N, T}) where {N, T}
    return AerosolDistribution{typeof(Modes)}(Modes)
end
Base.broadcastable(x::AerosolDistribution) = tuple(x)
n_modes(::AerosolDistribution{NTuple{N, T}}) where {N, T} = N

end
