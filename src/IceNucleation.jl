"""
Parameterization for heterogenous ice nucleation.
"""
module HetIceNucleation

import ..CommonTypes as CT
import ..Parameters as CMP
import Thermodynamics as TD

const APS = CMP.AbstractCloudMicrophysicsParameters

export dust_activated_number_fraction
export ABIFM_J

S0_warm(prs::APS, ::CT.ArizonaTestDustType) = CMP.S0_warm_ATD_Mohler2006(prs)
S0_cold(prs::APS, ::CT.ArizonaTestDustType) = CMP.S0_cold_ATD_Mohler2006(prs)
a_warm(prs::APS, ::CT.ArizonaTestDustType) = CMP.a_warm_ATD_Mohler2006(prs)
a_cold(prs::APS, ::CT.ArizonaTestDustType) = CMP.a_cold_ATD_Mohler2006(prs)

S0_warm(prs::APS, ::CT.DesertDustType) = CMP.S0_warm_DD_Mohler2006(prs)
S0_cold(prs::APS, ::CT.DesertDustType) = CMP.S0_cold_DD_Mohler2006(prs)
a_warm(prs::APS, ::CT.DesertDustType) = CMP.a_warm_DD_Mohler2006(prs)
a_cold(prs::APS, ::CT.DesertDustType) = CMP.a_cold_DD_Mohler2006(prs)

J_het_m(prs::APS, ::CT.DesertDustType) = CMP.J_ABIFM_m_KA2013_DesertDust(prs)
J_het_c(prs::APS, ::CT.DesertDustType) = CMP.J_ABIFM_c_KA2013_DesertDust(prs)
J_het_m(prs::APS, ::CT.KaoliniteType) = CMP.J_ABIFM_m_KA2013_Kaolinite(prs)
J_het_c(prs::APS, ::CT.KaoliniteType) = CMP.J_ABIFM_c_KA2013_Kaolinite(prs)
J_het_m(prs::APS, ::CT.IlliteType) = CMP.J_ABIFM_m_KA2013_Illite(prs)
J_het_c(prs::APS, ::CT.IlliteType) = CMP.J_ABIFM_c_KA2013_Illite(prs)

"""
    dust_activated_number_fraction(prs, Si, T, dust_type)

 - `prs` - set with model parameters
 - `Si` - ice saturation ratio
 - `T` - air temperature [K]
 - `dust_type` - a type for different dusts.

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

"""
    ABIFM_J(prs, dust_type, Δa_w)

 - `prs` - set with free parameters
 - `dust_type` - aerosol type
 - `Δa_w` - change in water activity [unitless].

Returns the immersion freezing nucleation rate coefficient, `J`, in m^-2 s^-1
for sulphuric acid solutions.
For other solutions, p_sol should be adjusted accordingly.
Delta_a_w can be found using the Delta_a_w function in Common.jl.
The free parameters `m` and `c` are taken from Knopf & Alpert 2013
see DOI: 10.1039/C3FD00035D
"""
function ABIFM_J(
    prs::APS,
    dust_type::CT.AbstractAerosolType,
    Δa_w::FT,
) where {FT <: Real}

    m::FT = J_het_m(prs, dust_type)
    c::FT = J_het_c(prs, dust_type)

    logJ::FT = m * Δa_w + c

    return max(FT(0), FT(10)^logJ * FT(1e4)) # converts cm^-2 s^-1 to m^-2 s^-1
end

end # end module

"""
Parameterization for homogeneous ice nucleation
"""
module HomIceNucleation

import ..CommonTypes as CT
import ..Parameters as CMP
import Thermodynamics as TD

const APS = CMP.AbstractCloudMicrophysicsParameters

export homogeneous_J

"""
    homogeneous_J(prs, Δa_w)

 - `prs` - a set with free parameters
 - `Δa_w` - change in water activity [-].

Returns the homogeneous freezing nucleation rate coefficient,
`J`, in m^-3 s^-1 for sulphuric acid solutions.
Parameterization based on Koop 2000, DOI: 10.1038/35020537.
Delta_a_w can be found using the Delta_a_w function in Common.jl.
"""
function homogeneous_J(prs::APS, Δa_w::FT) where {FT <: Real}

    Δa_w_min::FT = CMP.Koop2000_min_delta_aw(prs)
    Δa_w_max::FT = CMP.Koop2000_max_delta_aw(prs)
    c1::FT = CMP.Koop2000_J_hom_c1(prs)
    c2::FT = CMP.Koop2000_J_hom_c2(prs)
    c3::FT = CMP.Koop2000_J_hom_c3(prs)
    c4::FT = CMP.Koop2000_J_hom_c4(prs)

    @assert Δa_w > Δa_w_min
    @assert Δa_w < Δa_w_max

    logJ::FT = c1 + c2 * Δa_w - c3 * Δa_w^2 + c4 * Δa_w^3

    return FT(10)^(logJ) * 1e6
end

end # end module
