"""
Parameterization for heterogenous ice nucleation.
"""
module HetIceNucleation

import ..Parameters as CMP
import Thermodynamics as TD

export dust_activated_number_fraction
export deposition_J
export ABIFM_J

"""
    dust_activated_number_fraction(dust, ip, Si, T)

 - `dust` - a struct with dust parameters
 - `ip` - a struct with ice nucleation parameters
 - `Si` - ice saturation ratio
 - `T` - air temperature [K]

Returns the number fraction of mineral dust particles acting as
deposition nuclei (n ice nuclei / n dust particles).
From Mohler et al 2006 Table 2 (averages from different measurements
excluding those where a was not measured)
"""
function dust_activated_number_fraction(
    dust::Union{CMP.DesertDust, CMP.ArizonaTestDust},
    ip::CMP.Mohler2006,
    Si::FT,
    T::FT,
) where {FT}

    if Si > ip.Sᵢ_max
        @warn "Supersaturation exceeds the allowed value."
        @warn "No dust particles will be activated"
        return FT(0)
    else
        S₀::FT = T > ip.T_thr ? dust.S₀_warm : dust.S₀_cold
        a::FT = T > ip.T_thr ? dust.a_warm : dust.a_cold
        return max(0, exp(a * (Si - S₀)) - 1)
    end
end

"""
    deposition_J(dust, Δa_w)

 - `dust` - a struct with dust parameters
 - `Δa_w` - change in water activity [unitless].

Returns the deposition nucleation rate coefficient, `J`, in m^-2 s^-1
for different minerals in liquid droplets.
The free parameters `m` and `c` are derived from China et al (2017)
see DOI: 10.1002/2016JD025817
"""
function deposition_J(
    dust::Union{CMP.Ferrihydrite, CMP.Feldspar, CMP.Kaolinite},
    Δa_w::FT,
) where {FT}

    logJ::FT = dust.deposition_m * Δa_w + dust.deposition_c

    return max(FT(0), FT(10)^logJ * FT(1e4)) # converts cm^-2 s^-1 to m^-2 s^-1
end

"""
    ABIFM_J(dust, Δa_w)

 - `dust` - a struct with dust parameters
 - `Δa_w` - change in water activity [unitless].

Returns the immersion freezing nucleation rate coefficient, `J`, in m^-2 s^-1
for different minerals in liquid droplets.
The free parameters `m` and `c` are taken from Knopf & Alpert 2013
see DOI: 10.1039/C3FD00035D
"""
function ABIFM_J(
    dust::Union{CMP.DesertDust, CMP.Illite, CMP.Kaolinite},
    Δa_w::FT,
) where {FT}

    logJ::FT = dust.ABIFM_m * Δa_w + dust.ABIFM_c

    return max(FT(0), FT(10)^logJ * FT(1e4)) # converts cm^-2 s^-1 to m^-2 s^-1
end

end # end module

"""
Parameterization for homogeneous ice nucleation
"""
module HomIceNucleation

import ..Parameters as CMP
import Thermodynamics as TD

export homogeneous_J

"""
    homogeneous_J(ip, Δa_w)

 - `ip` - a struct with ice nucleation parameters
 - `Δa_w` - change in water activity [-].

Returns the homogeneous freezing nucleation rate coefficient,
`J`, in m^-3 s^-1 for sulphuric acid solutions.
Parameterization based on Koop 2000, DOI: 10.1038/35020537.
"""
function homogeneous_J(ip::CMP.Koop2000, Δa_w::FT) where {FT}

    @assert Δa_w > ip.Δa_w_min
    @assert Δa_w < ip.Δa_w_max

    logJ::FT = ip.c₁ + ip.c₂ * Δa_w - ip.c₃ * Δa_w^2 + ip.c₄ * Δa_w^3

    return FT(10)^(logJ) * 1e6
end

end # end module
