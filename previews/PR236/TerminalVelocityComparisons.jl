import Plots as PL

import CloudMicrophysics as CM
import CLIMAParameters as CP

FT = Float64

import CloudMicrophysics.CommonTypes as CMT
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CMO

const rain = CMT.RainType(FT)
const liquid = CMT.LiquidType()
const ice = CMT.IceType(FT)
const snow = CMT.SnowType(FT)
const tv_SB2006 = CMT.TerminalVelocitySB2006(FT)
const SB2006 = CMT.SB2006Type()
const Chen2022 = CMT.Chen2022Type(FT)
const SB2006Vel = CMT.SB2006VelType()
const Blk1MVel = CMT.Blk1MVelType(FT, rain)
const APS = CMP.AbstractCloudMicrophysicsParameters

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
const param_set = cloud_microphysics_parameters(toml_dict)

"""
    rain_terminal_velocity_individual_Chen(param_set, ρ, D_r)

 - `param_set` - set with free parameters
 - `ρ` - air density
 - `D_r` - diameter of the raindrops

Returns the fall velocity of a raindrop from Chen et al 2022
"""
function rain_terminal_velocity_individual_Chen(
    precip::CMT.RainType,
    velo_scheme::CMT.Chen2022Type,
    ρ::FT,
    D_r::FT, #in m
) where {FT <: Real}
    ai, bi, ci = CMO.Chen2022_vel_coeffs(precip, velo_scheme, ρ)
    D_r = D_r * 1000 #D_r is in mm in the paper --> multiply D_r by 1000

    v = 0
    for i in 1:3
        v += (ai[i] * (D_r)^(bi[i]) * exp(-D_r * ci[i]))
    end
    return v
end

"""
    rain_terminal_velocity_individual_SB(param_set, ρ, D_r)

 - `param_set` - set with free parameters
 - `ρ` - air density
 - `D_r` - diameter of the raindrops

Returns the fall velocity of a raindrop from Seifert and Beheng 2006
"""
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

"""
    terminal_velocity_individual_1M(prs, precip, ρ, D_r)

 - `prs` - set with free parameters
 - `precip` - precipitation type (rain or snow)
 - `ρ` - air density
 - `D_r` - particle diameter

Returns the fall velocity of a raindrop or snow from 1-moment scheme
"""
function terminal_velocity_individual_1M(
    precip::CMT.AbstractPrecipType,
    velo_scheme::CMT.Blk1MVelType,
    ρ::FT,
    D_r::FT,
) where {FT <: Real}
    (; r0) = precip
    (; χv, ve, Δv) = velo_scheme
    v0 = CM1.get_v0(precip, ρ)
    vt = χv * v0 * (D_r / (2 * r0))^(Δv + ve)
    return vt
end

"""
    ice_terminal_velocity_individual_Chen(param_set, ρ, D_r)

 - `param_set` - set with free parameters
 - `ρ` - air density
 - `D_r` - diameter of the raindrops

Returns the fall velocity of an ice particle from Chen et al 2022
"""
function ice_terminal_velocity_individual_Chen(
    precip::CMT.IceType,
    velo_scheme::CMT.Chen2022Type,
    ρ::FT,
    D_r::FT, #in m
) where {FT <: Real}

    ρ_i = precip.ρ
    (; As, Bs, Cs, Es, Fs, Gs) = velo_scheme.snowice

    ai = [Es * ρ^As, Fs * ρ^As]
    bi = [Bs + ρ * Cs, Bs + ρ * Cs]
    ci = [0, Gs]

    D_r = D_r * 1000 #D_r is in mm in the paper --> multiply D_r by 1000
    v = 0
    for i in 1:2
        v += (ai[i] * (D_r)^(bi[i]) * exp(-D_r * ci[i]))
    end
    return v
end

"""
    snow_terminal_velocity_individual_Chen(prs, ρ, D_r)

 - `prs` - set with free parameters
 - `ρ` - air density
 - `D_r` - particle diameter

Returns the fall velocity of snow from Chen et al 2022
"""
function snow_terminal_velocity_individual_Chen(
    precip::CMT.SnowType,
    velo_scheme::CMT.Chen2022Type,
    ρ::FT,
    D_r::FT,
) where {FT <: Real}
    (; r0, m0, me, Δm, χm, a0, ae, Δa, χa) = precip
    D_r = D_r * 1000
    ρ_i::FT = precip.ρ
    m0_comb = m0 * χm
    a0_comb = a0 * χa
    me_comb = me + Δm
    ae_comb = ae + Δa
    α = FT(-1 / 3)

    (; As, Bs, Cs, Es, Fs, Gs) = velo_scheme.snowice

    ai = [Es * ρ^As, Fs * ρ^As]
    bi = [Bs + ρ * Cs, Bs + ρ * Cs]
    ci = [0, Gs]

    aspect_ratio =
        ((16 * (a0_comb)^3 * (ρ_i)^2) / (9 * π * (m0_comb)^2)) *
        (D_r / (2 * r0 * 1000))^(3 * ae_comb - 2 * me_comb)
    aspect_ratio = aspect_ratio^(α)
    vt = 0
    for i in 1:2
        vt = vt + aspect_ratio * ai[i] * ((D_r))^(bi[i]) * exp(-1 * ci[i] * D_r)
    end
    return vt
end

function aspect_ratio_snow_1M(snow::CMT.SnowType, D_r::FT) where {FT <: Real}
    (; m0, me, Δm, a0, ae, Δa, χm, χa, r0) = snow
    ρ_i = snow.ρ

    m0_comb = m0 * χm
    a0_comb = a0 * χa
    me_comb = me + Δm
    ae_comb = ae + Δa

    aspect_ratio =
        ((16 * (a0_comb)^3 * (ρ_i)^2) / (9 * π * (m0_comb)^2)) *
        (D_r / (2 * r0))^(3 * ae_comb - 2 * me_comb)
    return aspect_ratio
end

ρ_air, ρ_air_ground = 1.2, 1.22
q_liq, q_ice, q_tot = 5e-4, 5e-4, 20e-3
N_rai = 1e7 #this value matters for the individual terminal velocity
q_rain_range = range(1e-8, stop = 5e-3, length = 100)
D_r_range = range(1e-6, stop = 5e-3, length = 1000)
N_d_range = range(1e6, stop = 1e9, length = 1000)
D_r_ar_range = range(1e-6, stop = 6.25e-4, length = 1000)

#! format: off

SB_rain_bN = [CM2.rain_terminal_velocity(rain, tv_SB2006, q_rai, ρ_air, N_rai)[1] for q_rai in q_rain_range]
Ch_rain_bN = [CM2.rain_terminal_velocity(rain, Chen2022, tv_SB2006, q_rai, ρ_air, N_rai)[1] for q_rai in q_rain_range]
SB_rain_bM = [CM2.rain_terminal_velocity(rain, tv_SB2006, q_rai, ρ_air, N_rai)[2] for q_rai in q_rain_range]
Ch_rain_bM = [CM2.rain_terminal_velocity(rain, Chen2022, tv_SB2006, q_rai, ρ_air, N_rai)[2] for q_rai in q_rain_range]
M1_rain_bM = [CM1.terminal_velocity(rain, Blk1MVel, ρ_air, q_rai) for q_rai in q_rain_range]

SB_rain = [rain_terminal_velocity_individual_SB(param_set, ρ_air, D_r) for D_r in D_r_range]
Ch_rain = [rain_terminal_velocity_individual_Chen(rain, Chen2022, ρ_air, D_r) for D_r in D_r_range]
M1_rain = [terminal_velocity_individual_1M(rain, Blk1MVel, ρ_air, D_r) for D_r in D_r_range]
Ch_snow = [snow_terminal_velocity_individual_Chen(snow, Chen2022, ρ_air, D_r) for D_r in D_r_range]
M1_snow = [terminal_velocity_individual_1M(snow, Blk1MVel, ρ_air, D_r) for  D_r in D_r_range]
Ch_ice = [ice_terminal_velocity_individual_Chen(ice, Chen2022, ρ_air, D_r) for  D_r in D_r_range]

Aspect_Ratio = [aspect_ratio_snow_1M(snow, D_r) for  D_r in D_r_ar_range]

# 2 Moment Scheme Figures

# single drop comparison
p1 = PL.plot(D_r_range *1e3,  SB_rain, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", label  = "Rain-SB2006", color = :blue)
p1 = PL.plot!(D_r_range * 1e3, Ch_rain, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", label = "Rain-Chen2022", color = :red)
p1 = PL.plot!(D_r_range * 1e3, M1_rain, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", label = "Rain-1M", color = :green)
# group velocity comparison
p2 = PL.plot(q_rain_range * 1e3,  SB_rain_bN, linewidth = 3, xlabel = "q_rain [g/kg]", ylabel = "terminal velocity [m/s]", label = "Rain-SB2006 [ND]", linestyle = :dot, color = :blue)
p2 = PL.plot!(q_rain_range * 1e3, Ch_rain_bN, linewidth = 3, xlabel = "q_rain [g/kg]", ylabel = "terminal velocity [m/s]", label = "Rain-Chen2022 [ND]", linestyle = :dot, color = :red)
p2 = PL.plot!(q_rain_range * 1e3, SB_rain_bM, linewidth = 3, xlabel = "q_rain [g/kg]", ylabel = "terminal velocity [m/s]", label = "Rain-SB2006 [M]", color = :blue)
p2 = PL.plot!(q_rain_range * 1e3, Ch_rain_bM, linewidth = 3, xlabel = "q_rain [g/kg]", ylabel = "terminal velocity [m/s]", label = "Rain-Chen2022 [M]", color = :red)
p2 = PL.plot!(q_rain_range * 1e3, M1_rain_bM, linewidth = 3, xlabel = "q_rain [g/kg]", ylabel = "terminal velocity [m/s]", label = "Rain-default-1M [M]", color = :green)
# save plot
PL.plot(p1, p2, layout = (1, 2), size = (800, 500), dpi = 300)
PL.savefig("2M_terminal_velocity_comparisons.svg")

# 1 Moment Scheme Figures

# individual particle comparison
p1 = PL.plot(D_r_range * 1e3,  Ch_rain, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", label = "Rain-Chen",       color = :blue,  linestyle = :dot)
p1 = PL.plot!(D_r_range * 1e3, M1_rain, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", label = "Rain-default-1M", color = :blue)
p1 = PL.plot!(D_r_range * 1e3, Ch_snow, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", label = "Snow-Chen",       color = :red,   linestyle = :dot)
p1 = PL.plot!(D_r_range * 1e3, M1_snow, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", label = "Snow-default-1M", color = :red)
p1 = PL.plot!(D_r_range * 1e3, Ch_ice, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]",  label = "Ice-Chen",        color = :pink,  linestyle = :dot)
# snow aspect ratio plot
p2 = PL.plot(D_r_ar_range * 1e3,  Aspect_Ratio, linewidth = 3, xlabel = "D [mm]", ylabel = "aspect ratio", label = "Aspect Ratio", color = :red)
# save plot
PL.plot(p1, p2, layout = (1, 2), size = (800, 500), dpi = 300)
PL.savefig("1M_individual_terminal_velocity_comparisons.svg")

#! format: on
