"""
    Two-moment bulk microphysics autoconversion and accretion rates, which includes:
    - autoconversion of cloud liquid water into rain
    - accretion rates of rain
    from Khairoutdinov and Kogan (2000), Beheng (1994), Tripoli and Cotton (1980), 
    Liu and Daum (2004).

"""
module Microphysics2M

import SpecialFunctions
const SF = SpecialFunctions

import Thermodynamics
const TD = Thermodynamics

import CLIMAParameters
const CP = CLIMAParameters
const APS = CP.AbstractParameterSet

import ..CommonTypes
const CT = CommonTypes

import ..InternalClimaParams
const ICP = InternalClimaParams

export conv_q_liq_to_q_rai_KK2000
export conv_q_liq_to_q_rai_B1994
export conv_q_liq_to_q_rai_TC1980
export heaviside

"""
    heaviside(x)

 - `x` - some variable or expression

Returns 1 if x > 0 and 0 if x ≤ 0
"""
function heaviside(x::FT) where {FT <: Real}
    if x > 0
        return 1
    else 
        return 0
    end
end

# autoconversion rates

"""
    conv_q_liq_to_q_rai_KK2000(param_set, q_liq, ρ, N_d)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - rain water humidity
 - `ρ` - rain water density (DOUBLE CHECK THIS)
 - `N_d` - cloud droplet concentration

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Khairoutdinov and Kogan (2000).
"""
function conv_q_liq_to_q_rai_KK2000(param_set::APS, q_liq::FT, ρ::FT, N_d::FT) where {FT <: Real}

    A::FT = 7.42e13
    a::FT = 2.47
    b::FT = -1.79
    c::FT = -1.47

    return A*q_liq^a*N_d^b*ρ^c
end

"""
    conv_q_liq_to_q_rai_B1994(param_set, q_liq, N_d)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - rain water humidity
 - `N_d` - cloud droplet concentration

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Beheng (1994).
"""
function conv_q_liq_to_q_rai_B1994(param_set::APS, q_liq::FT, N_d::FT) where {FT <: Real}

    C::FT = 3e34
    a::FT = -1.7
    b::FT = 4.7
    c::FT = -3.3
    if N_d::FT >= 2e8
        d = 3.9
    else
        d = 9.9
    end
        
    return C*d^a*q_liq^b*N_d^c
end

"""
    conv_q_liq_to_q_rai_TC1980(param_set, q_liq, N_d)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - rain water humidity
 - `N_d` - cloud droplet concentration

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Tripoli and Cotton (1980).
"""
function conv_q_liq_to_q_rai_TC1980(param_set::APS, q_liq::FT, N_d::FT) where {FT <: Real}

    D = 3268
    a = 7/3
    b = -1/3
    ρ_w = 1e3 # kg/m^3
    r_cm = 7e-6 # m
    _q_liq_threshold = 4/3*π*ρ_w*N_d*r_cm^3 # kg
    H = heaviside(q_liq - _q_liq_threshold)
    
    return D*q_liq^a*N_d^b*H
end

"""
    conv_q_liq_to_q_rai_LD2004(param_set, q_liq, N_d)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - rain water humidity
 - `N_d` - cloud droplet concentration

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Liu and Daum (2004).
"""
function conv_q_liq_to_q_rai_LD2004(param_set::APS, q_liq::FT, N_d::FT) where {FT <: Real}

    ρ_w::FT = 1e3 # kg/m^3
    r_0::FT = 1e-3 # m
    m_e::FT = 3
    Δ_m::FT = 0
    χ_m::FT = 1
    m_0::FT = 4/3*pi*ρ_w*r_0^3 # kg
    n_0::FT = 16e6 # 1/m^4
    ρ::FT = 1.2 # kg/m^3
    λ = (SF.gamma(m_e+Δ_m+1)*χ_m*m_0*n_0)/(q_liq*ρ*r_0^(m_e+Δ_m))^(1/(m_e+Δ_m+1))

    r_v::FT = 6^(1/3)/λ * 1e6
    β_6::FT = ((r_v+3)/r_v)^(1/3)
    E::FT = 1.08e10*β_6^6
    R_6::FT = β_6*r_v
    R_6C::FT = 2.25/(q_liq^(1/6)*R_6^(1/2))
    a::FT = 3
    b::FT = -1

    return E*q_liq^a*N_d^b*heaviside(R_6 - R_6C)
end

# accretion rates

"""
    accretion_KK2000(param_set, q_liq, q_rai, ρ)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - rain water humidity
 - `q_rai` - rain specific humidity
 - `ρ` - rain water density

 Returns the accretion rate of rain, parametrized
 following Khairoutdinov and Kogan (2000).
"""
function accretion_KK2000(
    param_set::APS,
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
) where {FT <: Real}

    A::FT = 67.0
    a::FT = 1.15
    b::FT = -1.3
    return A*(q_liq*q_rai)^a*ρ^b
end

"""
    accretion_B1994(param_set, q_liq, q_rai)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - rain water humidity
 - `q_rai` - rain specific humidity

 Returns the accretion rate of rain, parametrized
 following Beheng (1994).
"""
function accretion_B1994(
    param_set::APS,
    q_liq::FT,
    q_rai::FT, 
) where {FT <: Real}

    B::FT = 6.0
    return B*q_liq*q_rai
end

"""
    accretion_KK2000(param_set, q_liq, q_rai)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - rain water humidity
 - `q_rai` - rain specific humidity

 Returns the accretion rate of rain, parametrized
 following Tripoli and Cotton (1980).
"""
function accretion_TC1980(
    param_set::APS,
    q_liq::FT,
    q_rai::FT, 
) where {FT <: Real}

    C::FT = 4.7
    return C*q_liq*q_rai
end

end