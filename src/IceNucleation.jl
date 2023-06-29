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

J_het_coeff_m(::CT.DesertDustType) = 22.62
J_het_coeff_c(::CT.DesertDustType) = -1.35

J_het_coeff_m(::CT.KaoliniteType) = 54.58834
J_het_coeff_c(::CT.KaoliniteType) = -10.54758

J_het_coeff_m(::CT.IlliteType) = 54.48075
J_het_coeff_c(::CT.IlliteType) = -10.66873

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
    ABIFM_J(dust_type, Δa_w)

 - `dust_type` - choosing aerosol type
 - `Δa_w` - change in water activity [unitless].

Returns the immersion freezing nucleation rate coefficient, `J`, in m^-2 s^-1 for sulphuric acid containing solutions. 
For other solutions, p_sol should be adjusted accordingly. Delta_a_w can be found using the Delta_a_w function in Common.jl. 
`m` and `c` constants are taken from Knopf & Alpert 2013.
"""
function ABIFM_J(dust_type::CT.AbstractAerosolType, Δa_w::FT) where {FT <: Real}

    m = J_het_coeff_m(dust_type)
    c = J_het_coeff_c(dust_type)

    logJ = m * Δa_w + c

    return max(0, 10^logJ * 100^2) # converts cm^-2 s^-1 to m^-2 s^-1
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
    homogeneous_J(Δa_w)

 - `Δa_w` - change in water activity

Returns the homogeneous freezing nucleation rate coefficient, `J`, in m^-3 s^-1 for sulphuric
acid containing solutions. Parameterization based off Koop 2000.
Delta_a_w can be found using the Delta_a_w function in Common.jl.
"""
function homogeneous_J(Δa_w::FT) where {FT <: Real}

    @assert Δa_w > 0.26
    @assert Δa_w < 0.34

    logJ = -906.7 + 8502 * Δa_w - 26924 * Δa_w^2 + 29180 * Δa_w^3
    J = 10^(logJ)

    return J * 1e6
end

end # end module
