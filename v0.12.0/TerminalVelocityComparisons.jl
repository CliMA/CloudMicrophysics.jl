import Plots as PL

FT = Float64

import CloudMicrophysics as CM
import CLIMAParameters as CP

const CMT = CM.CommonTypes
const CM1 = CM.Microphysics1M
const CM2 = CM.Microphysics2M
const CMP = CM.Parameters

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
const param_set = cloud_microphysics_parameters(toml_dict)

const rain = CMT.RainType()
const SB2006 = CMT.SB2006Type()
const Chen2022 = CMT.Chen2022Type()

"""
    rain_terminal_velocity_individual_C(ρ, scheme, D_r)

 - `ρ`` - air density
 - `scheme` - type for parameterization
 - `D_r` - diameter of the raindrops

 Returns the fall velocity of a raindrop using multiple gamma-function terms
 to increase accuracy for `scheme == Chen2022Type`
"""
function rain_terminal_velocity_individual_Chen(
    param_set,
    ρ::FT,
    scheme::CMT.Chen2022Type,
    D_r::FT, #in m
) where {FT <: Real}

    # coefficients from Table B1 from Chen et. al. 2022
    ρ0::FT = CMP.q_coeff_rain_Ch2022(param_set)
    a1::FT = CMP.a1_coeff_rain_Ch2022(param_set)
    a2::FT = CMP.a2_coeff_rain_Ch2022(param_set)
    a3::FT = CMP.a3_coeff_rain_Ch2022(param_set)
    a3_pow::FT = CMP.a3_pow_coeff_rain_Ch2022(param_set)
    b1::FT = CMP.b1_coeff_rain_Ch2022(param_set)
    b2::FT = CMP.b2_coeff_rain_Ch2022(param_set)
    b3::FT = CMP.b3_coeff_rain_Ch2022(param_set)
    b_ρ::FT = CMP.b_rho_coeff_rain_Ch2022(param_set)
    c1::FT = CMP.c1_coeff_rain_Ch2022(param_set)
    c2::FT = CMP.c2_coeff_rain_Ch2022(param_set)
    c3::FT = CMP.c3_coeff_rain_Ch2022(param_set)

    q = exp(ρ0 * ρ)
    ai = (a1 * q, a2 * q, a3 * q * ρ^a3_pow)
    bi = (b1 - b_ρ * ρ, b2 - b_ρ * ρ, b3 - b_ρ * ρ)
    ci = (c1, c2, c3)

    D_r = D_r * 1000 #D_r is in mm in the paper --> multiply D_r by 1000

    v = 0
    for i in 1:3
        v += (ai[i] * (D_r)^(bi[i]) * exp(-D_r * ci[i]))
    end
    return v
end

function rain_terminal_velocity_individual_SB(
    param_set,
    ρ::FT,
    D_r::FT,
) where {FT <: Real}

    a_r::FT = CMP.aR_tv_SB2006(param_set)
    b_r::FT = CMP.bR_tv_SB2006(param_set)
    c_r::FT = CMP.cR_tv_SB2006(param_set)
    ρ_air_ground::FT = CMP.ρ0_SB2006(param_set)

    v = (ρ_air_ground / ρ)^(1 / 2) * (a_r - b_r * exp(-c_r * D_r))
    return v
end

q_liq_range = range(1e-8, stop = 1e-3, length = 1000)
q_rai_range = range(1e-8, stop = 1e-3, length = 1000)
D_r_range = range(50e-6, stop = 9e-3, length = 1000)
N_d_range = range(1e7, stop = 1e9, length = 1000)
q_liq = 5e-4
# q_rai = 5e-4
ρ_air = 1.0 # kg m^-3
N_rai = 1e4 #this value matters for the individual terminal velocity

PL.plot(
    q_rai_range * 1e3,
    [
        CM2.rain_terminal_velocity(param_set, SB2006, q_rai, ρ_air, N_rai)[1]
        for q_rai in q_rai_range
    ],
    linewidth = 3,
    xlabel = "q_rain [g/kg]",
    ylabel = "terminal velocity [m/s]",
    label = "Rain-SB2006 [ND]",
)
PL.plot!(
    q_rai_range * 1e3,
    [
        CM2.rain_terminal_velocity(param_set, Chen2022, q_rai, ρ_air, N_rai)[1]
        for q_rai in q_rai_range
    ],
    linewidth = 3,
    xlabel = "q_rain [g/kg]",
    ylabel = "terminal velocity [m/s]",
    label = "Rain-Chen2022 [ND]",
)
PL.plot!(
    q_rai_range * 1e3,
    [
        CM2.rain_terminal_velocity(param_set, SB2006, q_rai, ρ_air, N_rai)[2]
        for q_rai in q_rai_range
    ],
    linewidth = 3,
    xlabel = "q_rain [g/kg]",
    ylabel = "terminal velocity [m/s]",
    label = "Rain-SB2006 [M]",
)
PL.plot!(
    q_rai_range * 1e3,
    [
        CM2.rain_terminal_velocity(param_set, Chen2022, q_rai, ρ_air, N_rai)[2]
        for q_rai in q_rai_range
    ],
    linewidth = 3,
    xlabel = "q_rain [g/kg]",
    ylabel = "terminal velocity [m/s]",
    label = "Rain-Chen2022 [M]",
)
PL.plot!(
    q_rai_range * 1e3,
    [
        CM1.terminal_velocity(param_set, rain, ρ_air, q_rai) for
        q_rai in q_rai_range
    ],
    linewidth = 3,
    xlabel = "q_rain [g/kg]",
    ylabel = "terminal velocity [m/s]",
    label = "Rain-Ogura-1-moment",
)

PL.title!("Terminal Velocity vs. Specific Humidity of Rain", titlefontsize = 9)
PL.savefig("terminal_velocity_bulk_comparisons.svg")

PL.plot(
    D_r_range,
    [
        rain_terminal_velocity_individual_SB(param_set, ρ_air, D_r) for
        D_r in D_r_range
    ],
    linewidth = 3,
    xlabel = "D [m]",
    ylabel = "terminal velocity [m/s]",
    label = "Rain-SB2006",
)
PL.plot!(
    D_r_range,
    [
        rain_terminal_velocity_individual_Chen(param_set, ρ_air, Chen2022, D_r)
        for D_r in D_r_range
    ],
    linewidth = 3,
    xlabel = "D [m]",
    ylabel = "terminal velocity [m/s]",
    label = "Rain-Chen2022",
)
PL.title!("Terminal Velocity vs. Diameter", titlefontsize = 13)
PL.savefig("terminal_velocity_individual_raindrop.svg")
