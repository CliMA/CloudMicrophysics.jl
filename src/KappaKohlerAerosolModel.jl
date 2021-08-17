"""
    A container for information on aerosol size distribution
    and chemical properties.

    The size distribution is a sum of lognormal internally mixed modes.
"""
#module KappaKohlerAerosolModel

#export KappaKohlerMode
#export KappaKohlerAerosolDistribution

"""
    Mode

Represents the sizes and chemical composition
of aerosol particles in one size distribution mode.
The mode is assumed to be made up of internally mixed components
and follow a lognormal size distribution.
"""
struct KappaKohlerMode{T}
    "geometric mean dry radius"
    r_dry::Real
    "geometric standard deviation"
    stdev::Real
    "total number concentration"
    N::Real
    "tuple of mass mixing ratios for all components in this mode"
    mass_mix_ratio::T
    "tuple of molar masses for all components in this mode"
    molar_mass::T
    "tuple of kappa-kohler values for all components in this mode"
    kappa::Real
    "number of components in the mode"
    n_components::Int64
end


"""
    AerosolDistribution

Represents the aerosol size distribution as a tuple with all modes.

# Constructors

    AerosolDistribution(Modes::T)
"""
struct KappaKohlerAerosolDistribution{T}

    "tuple with all aerosol size distribution modes"
   KK_Modes::T

    function KappaKohlerAerosolDistribution(Modes::T) where {T}
        return new{T}(Modes)
    end

end

#end
