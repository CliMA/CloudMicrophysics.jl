"""
Parameterization for heterogenous ice nucleation.
"""
module HetIceNucleation

import ..Parameters as CMP
import Thermodynamics as TD

export dust_activated_number_fraction
export MohlerDepositionRate
export deposition_J
export ABIFM_J
export P3_deposition_N_i
export P3_het_N_i
export INP_concentration_frequency
export INP_concentration_mean

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

    @assert Si < ip.Sᵢ_max

    S₀::FT = T > ip.T_thr ? dust.S₀_warm : dust.S₀_cold
    a::FT = T > ip.T_thr ? dust.a_warm : dust.a_cold
    return max(0, exp(a * (Si - S₀)) - 1)
end

"""
    MohlerDepositionRate(dust, ip, S_i, T, dSi_dt, N_aer)

 - `dust` - a struct with dust parameters
 - `ip` - a struct with ice nucleation parameters
 - `Si` - ice saturation
 - `T` - ambient temperature
 - `dSi_dt` - change in ice saturation over time
 - `N_aer` - number of unactivated aerosols

Returns the ice nucleation rate from deposition.
From Mohler et al 2006 equation 5.
"""
function MohlerDepositionRate(
    dust::Union{CMP.DesertDust, CMP.ArizonaTestDust},
    ip::CMP.Mohler2006,
    Si::FT,
    T::FT,
    dSi_dt::FT,
    N_aer::FT,
) where {FT}

    @assert Si < ip.Sᵢ_max

    a::FT = T > ip.T_thr ? dust.a_warm : dust.a_cold
    return max(0, N_aer * a * dSi_dt)
end

"""
    deposition_J(dust, Δa_w)

 - `dust` - a struct with dust parameters
 - `Δa_w` - change in water activity [unitless].

Returns the deposition nucleation rate coefficient, `J`, in m^-2 s^-1
for different minerals in liquid droplets.
The free parameters `m` and `c` are derived from China et al (2017)
see DOI: 10.1002/2016JD025817
Returns zero for unsupported aerosol types.
"""
function deposition_J(
    dust::Union{
        CMP.Ferrihydrite,
        CMP.Feldspar,
        CMP.Kaolinite,
        CMP.Illite,
        CMP.ArizonaTestDust,
    },
    Δa_w::FT,
) where {FT}

    logJ::FT = dust.deposition_m * Δa_w + dust.deposition_c

    return max(FT(0), FT(10)^logJ * FT(1e4)) # converts cm^-2 s^-1 to m^-2 s^-1
end
function deposition_J(dust::CMP.AerosolType, Δa_w::FT) where {FT}
    #println("Aerosol type not supported for ABDINM.")
    return FT(0)
end

"""
    ABIFM_J(dust, Δa_w)

 - `dust` - a struct with dust parameters
 - `Δa_w` - change in water activity [unitless].

Returns the immersion freezing nucleation rate coefficient, `J`, in m^-2 s^-1
for different minerals in liquid droplets.
The free parameters `m` and `c` are taken from Knopf & Alpert 2013
see DOI: 10.1039/C3FD00035D
Returns zero for unsupported aerosol types.
"""
function ABIFM_J(
    dust::Union{CMP.DesertDust, CMP.Illite, CMP.Kaolinite, CMP.ArizonaTestDust},
    Δa_w::FT,
) where {FT}

    logJ::FT = dust.ABIFM_m * Δa_w + dust.ABIFM_c

    return max(FT(0), FT(10)^logJ * FT(1e4)) # converts cm^-2 s^-1 to m^-2 s^-1
end
function ABIFM_J(dust::CMP.AerosolType, Δa_w::FT) where {FT}
    #println("Aerosol type not supported for ABIFM.")
    return FT(0)
end

"""
    P3_deposition_N_i(ip, T)

 - `ip` - a struct with ice nucleation parameters,
 - `T` - air temperature [K].

Returns the number of ice nucleated via deposition nucleation with units of m^-3.
From Thompson et al 2004 eqn 2 as used in Morrison & Milbrandt 2015.
"""
function P3_deposition_N_i(ip::CMP.MorrisonMilbrandt2014, T::FT) where {FT}

    T₀ = ip.T₀                  # 0°C
    T_thres = ip.T_dep_thres    # cutoff temperature

    Nᵢ = ifelse(
        T < T_thres,
        max(FT(0), FT(ip.c₁ * exp(ip.c₂ * (T₀ - T_thres)))),
        max(FT(0), FT(ip.c₁ * exp(ip.c₂ * (T₀ - T)))),
    )
    return Nᵢ * 1e3  # converts L^-1 to m^-3
end

"""
    P3_het_N_i(ip, T, N_l, B, V_l, a, Δt)

 - `ip` - a struct with ice nucleation parameters,
 - `T` - air temperature [K],
 - `N_l` - number of droplets [m^-3],
 - `B` - water-type dependent parameter [cm^-3 s^-1],
 - `V_l` - volume of droplets to be frozen [m^3],
 - `a` - empirical parameter [C^-1],
 - `Δt` - timestep.

Returns the number of ice nucleated within Δt via heterogeneous freezing
with units of m^-3. From Pruppacher & Klett 1997 eqn (9-51) as used in
Morrison & Milbrandt 2015.
"""
function P3_het_N_i(
    ip::CMP.MorrisonMilbrandt2014,
    T::FT,
    N_l::FT,
    V_l::FT,
    Δt::FT,
) where {FT}

    a = ip.het_a                # (celcius)^-1
    B = ip.het_B                # cm^-3 s^-1 for rain water
    V_l_converted = V_l * 1e6   # converted from m^3 to cm^3
    Tₛ = ip.T₀ - T

    return N_l * (1 - exp(-B * V_l_converted * Δt * exp(a * Tₛ)))
end

"""
    INP_concentration_frequency(params,INPC,T)

 - `params` - a struct with INPC(T) distribution parameters
 - `INPC` - concentration of ice nucleating particles [m^-3]
 - `T` - air temperature [K]

Returns the relative frequency of a given INP concentration,
depending on the temperature.
Based on Frostenberg et al., 2023. See DOI: 10.5194/acp-23-10883-2023
"""
function INP_concentration_frequency(
    params::CMP.Frostenberg2023,
    INPC::FT,
    T::FT,
) where {FT}

    μ = INP_concentration_mean(T)

    return 1 / (sqrt(2 * FT(π)) * params.σ) *
           exp(-(log(INPC) - μ)^2 / (2 * params.σ^2))
end

"""
    INP_concentration_mean(T)

 - `T` - air temperature [K]

Returns the logarithm of mean INP concentration (in m^-3), depending on the temperature.
Based on the function μ(T) in Frostenberg et al., 2023. See DOI: 10.5194/acp-23-10883-2023
"""
function INP_concentration_mean(T::FT) where {FT}

    T_celsius = T - FT(273.15)

    return log(-(T_celsius)^9 * 10^(-9))
end

end # end module

"""
Parameterization for homogeneous ice nucleation
"""
module HomIceNucleation

import ..Parameters as CMP
import Thermodynamics as TD

export homogeneous_J_cubic
export homogeneous_J_linear

"""
    homogeneous_J_cubic(ip, Δa_w)

 - `ip` - a struct with ice nucleation parameters
 - `Δa_w` - change in water activity [-].

Returns the homogeneous freezing nucleation rate coefficient,
`J`, in m^-3 s^-1 for sulphuric acid solutions.
Parameterization based on Koop 2000, DOI: 10.1038/35020537.
"""
function homogeneous_J_cubic(ip::CMP.Koop2000, Δa_w::FT) where {FT}

    @assert Δa_w >= ip.Δa_w_min
    @assert Δa_w <= ip.Δa_w_max

    logJ::FT = ip.c₁ + ip.c₂ * Δa_w - ip.c₃ * Δa_w^2 + ip.c₄ * Δa_w^3

    return FT(10)^(logJ) * 1e6
end

"""
    homogeneous_J_linear(ip, Δa_w)

 - `ip` - a struct with ice nucleation parameters
 - `Δa_w` - change in water activity [-].

Returns the homogeneous freezing nucleation rate coefficient,
`J`, in m^-3 s^-1 for sulphuric acid solutions.
Parameterization derived from a linear fit of the Koop 2000 parameterization, DOI: 10.1038/35020537.
"""
function homogeneous_J_linear(ip::CMP.Koop2000, Δa_w::FT) where {FT}

    logJ::FT = ip.linear_c₂ * Δa_w + ip.linear_c₁

    return FT(10)^(logJ) * 1e6
end

end # end module
