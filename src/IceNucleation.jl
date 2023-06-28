"""
Parameterization for heterogenous ice nucleation.
"""
module HetIceNucleation

import ..CommonTypes as CT
import ..Parameters as CMP

const APS = CMP.AbstractCloudMicrophysicsParameters

export dust_activated_number_fraction

S0_warm(prs::APS, ::CT.ArizonaTestDustType) = CMP.S0_warm_ATD_Mohler2006(prs)
S0_cold(prs::APS, ::CT.ArizonaTestDustType) = CMP.S0_cold_ATD_Mohler2006(prs)
a_warm(prs::APS, ::CT.ArizonaTestDustType) = CMP.a_warm_ATD_Mohler2006(prs)
a_cold(prs::APS, ::CT.ArizonaTestDustType) = CMP.a_cold_ATD_Mohler2006(prs)

S0_warm(prs::APS, ::CT.DesertDustType) = CMP.S0_warm_DD_Mohler2006(prs)
S0_cold(prs::APS, ::CT.DesertDustType) = CMP.S0_cold_DD_Mohler2006(prs)
a_warm(prs::APS, ::CT.DesertDustType) = CMP.a_warm_DD_Mohler2006(prs)
a_cold(prs::APS, ::CT.DesertDustType) = CMP.a_cold_DD_Mohler2006(prs)

"""
    dust_activated_number_fraction(prs, Si, T, dust_type)

 - `prs` - set with model parameters
 - `Si` - ice saturation ratio
 - `T` - air temperature [K]
 - `dust_type` - a type for different dustst.

Returns the number fraction of mineral dust particles acting as
deposition nuclei (n ice nuclei / n dust particles).
From Mohler et al 2006 Table 2 (averages from different measurements
excluding those where a was not measured)
"""
function dust_activated_number_fraction(
    prs::APS,
    Si::FT,
    T::FT,
    dust_type::CT.AbstractAerosolType,
) where {FT <: Real}

    Si_max::FT = CMP.Si_max_Mohler2006(prs)
    T_thr::FT = CMP.T_thr_Mohler2006(prs)

    if Si > Si_max
        @warn "Supersaturation exceedes the allowed value."
        @warn "No dust particles will be activated"
        return FT(0)
    else
        S0::FT = T > T_thr ? S0_warm(prs, dust_type) : S0_cold(prs, dust_type)
        a::FT = T > T_thr ? a_warm(prs, dust_type) : a_cold(prs, dust_type)
        return max(0, exp(a * (Si - S0)) - 1)
    end
end

end # end module
