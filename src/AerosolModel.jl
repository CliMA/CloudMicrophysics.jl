"""
    A container for information on aerosol size distribution
    and chemical properties.

    The size distribution is a sum of lognormal internally mixed modes.
"""
module AerosolModel

export Mode
export AerosolDistribution

"""
    Mode

Represents the sizes and chemical composition
of aerosol particles in one size distribution mode.
The mode is assumed to be made up of internally mixed components
and follow a lognormal size distribution.
"""
struct Mode{T, FT}
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
    AerosolDistribution

Represents the aerosol size distribution as a tuple with all modes.

# Constructors

    AerosolDistribution(Modes::T)
"""
struct AerosolDistribution{T}

    "tuple with all aerosol size distribution modes"
    Modes::T

    function AerosolDistribution(Modes::T) where {T}
        return new{T}(Modes)
    end

end

end
