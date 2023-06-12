"""
Parameterization for heterogenous ice nucleation.
"""
module HetIceNucleation

import ..CommonTypes as CT
import ..Parameters as CMP
import Thermodynamics as TD

const APS = CMP.AbstractCloudMicrophysicsParameters

export dust_activated_number_fraction
export H2SO4_soln_saturation_vapor_pressure
export ABIFM_Delta_a_w
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
    H2SO4_soln_saturation_vapor_pressure(x, T)

 - `x` - wt percent sulphuric acid [unitless] 
 - `T` - air temperature [K].

Returns the saturation vapor pressure above a sulphuric acid solution droplet in Pa.
`x` is, for example, 0.1 if droplets are 10 percent sulphuric acid by weight
"""
function H2SO4_soln_saturation_vapor_pressure(x::FT, T::FT) where {FT <: Real}

    @assert T < FT(235)
    @assert T > FT(185)

    w_h = 1.4408 * x
    p_sol =
        exp(
            23.306 - 5.3465 * x + 12 * x * w_h - 8.19 * x * w_h^2 +
            (-5814 + 928.9 * x - 1876.7 * x * w_h) / T,
        ) * 100 # * 100 converts mbar --> Pa
    return p_sol
end

"""
    ABIFM_Delta_a_w(prs, x, T)

 - `prs` - set with model parameters
 - `x` - wt percent sulphuric acid [unitless]
 - `T` - air temperature [K].

Returns the change in water activity when droplet undergoes immersion freezing.
`x` is, for example, 0.1 if droplets are 10 percent sulphuric acid by weight.
To be used in conjuction with ABIFM_J to obtain nucleation rate coefficient
"""
function ABIFM_Delta_a_w(prs::APS, x::FT, T::FT) where {FT <: Real}

    thermo_params = CMP.thermodynamics_params(prs)

    p_sol = H2SO4_soln_saturation_vapor_pressure(x, T)
    p_sat = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    p_ice = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())

    a_w = p_sol / p_sat
    a_w_ice = p_ice / p_sat
    Delta_a_w = a_w - a_w_ice

    return min(Delta_a_w, FT(1))
end

"""
    ABIFM_J(dust_type, Delta_a_w)

 - `dust_type` - choosing aerosol type
 - `Delta_a_w` - change in water activity [unitless].

Returns the immersion freezing nucleation rate coefficient, `J`, in m^-2 s^-1 for sulphuric acid containing solutions. 
For other solutions, p_sol should be adjusted accordingly. Delta_a_w can be found using ABIFM_Delta_a_w. 
`m` and `c` constants are taken from Knopf & Alpert 2013.
"""
function ABIFM_J(
    dust_type::CT.AbstractAerosolType,
    Delta_a_w::FT,
) where {FT <: Real}

    m = J_het_coeff_m(dust_type)
    c = J_het_coeff_c(dust_type)

    logJ = m * Delta_a_w + c

    return max(0, 10^logJ * 100^2) # converts cm^-2 s^-1 to m^-2 s^-1
end

end # end module
