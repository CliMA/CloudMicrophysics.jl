"""
    AerosolModel

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
end

""" number of components in the mode """
n_components(::Mode_B{T}) where {T <: Tuple} = fieldcount(T)
n_components(::Mode_B{T}) where {T <: Real} = 1

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
end

""" number of components in the mode """
n_components(::Mode_κ{T}) where {T <: Tuple} = fieldcount(T)
n_components(::Mode_κ{T}) where {T <: Real} = 1

"""
    AerosolDistribution

Represents the aerosol size distribution as a tuple with different modes.
All modes have to either be of type Mode_B (Abdul-Razzak and Ghan 2000)
or of type Mode_κ (Petters and Kreidenweis 2007).

# Constructors

    AerosolDistribution(modes::T)
    AerosolDistribution(modes::Union{Mode_κ, Mode_B}...)
"""
struct AerosolDistribution{T} <: CMP.AerosolDistributionType

    "tuple with all aerosol size distribution modes"
    modes::T

end
AerosolDistribution(modes::Union{Mode_κ, Mode_B}...) =
    AerosolDistribution{typeof(modes)}(modes)

Base.broadcastable(x::AerosolDistribution) = tuple(x)
n_modes(d::AerosolDistribution) = length(d.modes)

"""
    aerosol_distribution(pa::CMP.PrescribedAerosol)

The two-mode [`AerosolDistribution`](@ref) described by a
[`CMP.PrescribedAerosol`](@ref) parameter set: one single-component `Mode_κ` per mode.

Built at the point of use rather than stored, because the parameter module is loaded before this
one and cannot hold a `Mode_κ`. Both modes are `isbits` and share a concrete type, which is what
`mean_hygroscopicity_parameter` dispatches on, so the construction compiles away inside a GPU
kernel.

The volume and mass mixing ratios are one and the molar mass zero, as they are for any
single-component mode: the κ-Köhler path reads only `vol_mix_ratio` and `kappa`. A zero molar
mass makes `M_activated_per_mode` return zero, so the mass-activation entry points are not
available on a distribution built this way.
"""
function aerosol_distribution(pa::CMP.PrescribedAerosol{FT}) where {FT}
    one_component = (one(FT),)
    no_molar_mass = (zero(FT),)
    accum = Mode_κ(
        pa.r_dry_accum, pa.stdev_accum, pa.N_accum,
        one_component, one_component, no_molar_mass, (pa.κ_accum,),
    )
    coarse = Mode_κ(
        pa.r_dry_coarse, pa.stdev_coarse, pa.N_coarse,
        one_component, one_component, no_molar_mass, (pa.κ_coarse,),
    )
    return AerosolDistribution(accum, coarse)
end

end
