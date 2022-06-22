"""
    Two-moment bulk microphysics autoconversion and accretion rates, which includes:

"""
module Microphysics2M

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

function heaviside(x::FT) where {FT <: Real}
    if x > 0
        return 1
    else 
        return 0
    end
end

function conv_q_liq_to_q_rai_KK2000(param_set::APS, q_liq::FT, ρ::FT; N_d::FT = 1e8) where {FT <: Real}

    A::FT = 7.42e13
    a::FT = 2.47
    b::FT = -1.79
    c::FT = -1.47

    return A*q_liq^a*N_d^b*ρ^c
end

function conv_q_liq_to_q_rai_B1994(param_set::APS, q_liq::FT, N_d::FT = 1e8) where {FT <: Real}

    C::FT = 3e34
    a::FT = -1.7
    b::FT = 4.7
    c::FT = -3.3
    if N_d >= 2e8
        d = 3.9
    else
        d = 9.9
    end
        
    return C*d^a*q_liq^b*N_d^c
end

function conv_q_liq_to_q_rai_TC1980(param_set::APS, q_liq::FT, N_d::FT = 1e8) where {FT <: Real}

    D = 3268
    a = 7/3
    b = -1/3
    ρ_w = 1e3 # kg/m^3
    r_cm = 7e-6 # m
    _q_liq_threshold = 4/3*π*ρ_w*N_d*r_cm^3 # kg
    H = heaviside(q_liq - _q_liq_threshold)
    
    return D*q_liq^a*N_d^b*H
end

function conv_q_liq_to_q_rai_LD2004(param_set::APS, q_liq::FT, N_d::FT = 1e8) where {FT <: Real}

end

end