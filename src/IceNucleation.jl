"""
Parameterization for heterogenous ice nucleation.
"""
module HetIceNucleation

import ..CommonTypes
const CT = CommonTypes

export dust_activated_number_fraction

"""
    dust_activated_number_fraction(Si, T, dust_type)

 - `Si` - ice saturation ratio
 - `T` - air temperature [K]
 - `dust_type` - a type for different dustst.

Returns the number fraction of mineral dust particles acting as
deposition nuclei (n ice nuclei / n dust particles).
From Mohler et al 2006, empirical fit to data.
"""
function dust_activated_number_fraction(
    Si::FT,
    T::FT,
    dust_type::CT.ArizonaTestDustType,
) where {FT <: Real}

    # TODO - make S0 and a part of CLIMAParameters
    S0::FT = 1.05 # threshold ice saturatio ratio
    a::FT = T > FT(220) ? FT(4.7) : FT(9.2)

    return max(0, exp(a * (Si - S0)) - 1)
end
function dust_activated_number_fraction(
    Si::FT,
    T::FT,
    dust_type::CT.DesertDustType,
) where {FT <: Real}

    # TODO - make S0 and a part of CLIMAParameters
    S0::FT = T > FT(220) ? FT(1.14) : FT(1.05)
    a::FT = T > FT(220) ? FT(0.5) : FT(2)

    return max(0, exp(a * (Si - S0)) - 1)
end

end # end module
