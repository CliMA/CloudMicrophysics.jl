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
From Mohler et al 2006 Table 2 (averages from different measurements
excluding those where a was not measured)
"""
function dust_activated_number_fraction(
    Si::FT,
    T::FT,
    dust_type::CT.ArizonaTestDustType,
) where {FT <: Real}

    if Si > FT(1.35)
        @warn "Supersaturation exceedes the allowed value."
        @warn "No dust particles will be activated"
        return FT(0)
    else
        S0::FT = T > FT(220) ? FT(1.03) : FT(1.07)
        a::FT = T > FT(220) ? FT(4.7) : FT(9.2)
        return max(0, exp(a * (Si - S0)) - 1)
    end
end
function dust_activated_number_fraction(
    Si::FT,
    T::FT,
    dust_type::CT.DesertDustType,
) where {FT <: Real}

    if Si > FT(1.35)
        @warn "Supersaturation exceedes the allowed value."
        @warn "No dust particles will be activated"
        return FT(0)
    else
        S0::FT = T > FT(220) ? FT(1.17) : FT(1.03)
        a::FT = T > FT(220) ? FT(0.43) : FT(2.35)
        return max(0, exp(a * (Si - S0)) - 1)
    end
end

end # end module
