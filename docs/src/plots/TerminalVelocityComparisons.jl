import Plots as PL

import CLIMAParameters as CP

FT = Float64

import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CMO

const rain = CMP.Rain(FT)
const liquid = CMP.CloudLiquid(FT)
const ice = CMP.CloudIce(FT)
const snow = CMP.Snow(FT)

const SB2006 = CMP.SB2006(FT)

const Chen2022 = CMP.Chen2022VelType(FT)
const SB2006Vel = CMP.SB2006VelType(FT)
const Blk1MVel = CMP.Blk1MVelType(FT)

"""
    rain_terminal_velocity_individual_Chen(velo_scheme ρ, D_r)

 - `velo_scheme` - structs with free parameters
 - `ρ` - air density
 - `D_r` - diameter of the raindrops

Returns the fall velocity of a raindrop from Chen et al 2022
"""
function rain_terminal_velocity_individual_Chen(
    velo_scheme::CMP.Chen2022VelTypeRain,
    ρ::FT,
    D_r::FT, #in m
) where {FT <: Real}
    ai, bi, ci = CMO.Chen2022_vel_coeffs(velo_scheme, ρ)
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
    (; ρ0, aR, bR, cR)::CMP.SB2006VelType,
    ρ::FT,
    D_r::FT,
) where {FT <: Real}
    return v = (ρ0 / ρ)^(1 / 2) * (aR - bR * exp(-cR * D_r))
end

"""
    terminal_velocity_individual_1M(velo_scheme, ρ, D_r)

 - `velo_scheme` - set with free parameters
 - `ρ` - air density
 - `D_r` - particle diameter

Returns the fall velocity of a raindrop or snow from 1-moment scheme
"""
function terminal_velocity_individual_1M(
    velo_scheme::Union{CMP.Blk1MVelTypeRain, CMP.Blk1MVelTypeSnow},
    ρ::FT,
    D_r::FT,
) where {FT <: Real}
    (; χv, ve, Δv, r0) = velo_scheme
    v0 = CM1.get_v0(velo_scheme, ρ)
    vt = χv * v0 * (D_r / (2 * r0))^(Δv + ve)
    return vt
end

"""
    ice_terminal_velocity_individual_Chen(velo_scheme, ρ, D_r)

 - `velo_scheme` - set with free parameters
 - `ρ` - air density
 - `D_r` - diameter of the raindrops

Returns the fall velocity of an ice particle from Chen et al 2022
"""
function ice_terminal_velocity_individual_Chen(
    velo_scheme::CMP.Chen2022VelTypeSnowIce,
    ρ::FT,
    D_r::FT, #in m
) where {FT <: Real}

    (; As, Bs, Cs, Es, Fs, Gs) = velo_scheme

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
    snow_terminal_velocity_individual_Chen(precip, velo_scheme, ρ, D_r)

 - `precip`, `velo_scheme` - structs with free parameters
 - `ρ` - air density
 - `D_r` - particle diameter

Returns the fall velocity of snow from Chen et al 2022
"""
function snow_terminal_velocity_individual_Chen(
    precip::CMP.Snow,
    velo_scheme::CMP.Chen2022VelTypeSnowIce,
    ρ::FT,
    D_r::FT,
) where {FT <: Real}
    (; r0, m0, me, Δm, χm) = precip.mass
    (; a0, ae, Δa, χa) = precip.area

    D_r = D_r * 1000
    m0_comb = m0 * χm
    a0_comb = a0 * χa
    me_comb = me + Δm
    ae_comb = ae + Δa
    α = FT(-1 / 3)

    (; As, Bs, Cs, Es, Fs, Gs, ρᵢ) = velo_scheme
    ai = [Es * ρ^As, Fs * ρ^As]
    bi = [Bs + ρ * Cs, Bs + ρ * Cs]
    ci = [0, Gs]

    aspect_ratio =
        ((16 * (a0_comb)^3 * (ρᵢ)^2) / (9 * π * (m0_comb)^2)) *
        (D_r / (2 * r0 * 1000))^(3 * ae_comb - 2 * me_comb)
    aspect_ratio = aspect_ratio^(α)
    vt = 0
    for i in 1:2
        vt = vt + aspect_ratio * ai[i] * ((D_r))^(bi[i]) * exp(-1 * ci[i] * D_r)
    end
    return vt
end

function aspect_ratio_snow_1M(snow::CMP.Snow, D_r::FT) where {FT <: Real}
    (; r0, m0, me, χm, Δm) = snow.mass
    (; a0, ae, Δa, χa) = snow.area
    ρᵢ = snow.ρᵢ
    m0_comb = m0 * χm
    a0_comb = a0 * χa
    me_comb = me + Δm
    ae_comb = ae + Δa

    aspect_ratio =
        ((16 * (a0_comb)^3 * (ρᵢ)^2) / (9 * π * (m0_comb)^2)) *
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

SB_rain_bN = [CM2.rain_terminal_velocity(SB2006, SB2006Vel, q_rai, ρ_air, N_rai)[1] for q_rai in q_rain_range]
Ch_rain_bN = [CM2.rain_terminal_velocity(SB2006, Chen2022.rain, q_rai, ρ_air, N_rai)[1] for q_rai in q_rain_range]
SB_rain_bM = [CM2.rain_terminal_velocity(SB2006, SB2006Vel, q_rai, ρ_air, N_rai)[2] for q_rai in q_rain_range]
Ch_rain_bM = [CM2.rain_terminal_velocity(SB2006, Chen2022.rain, q_rai, ρ_air, N_rai)[2] for q_rai in q_rain_range]
M1_rain_bM = [CM1.terminal_velocity(rain, Blk1MVel.rain, ρ_air, q_rai) for q_rai in q_rain_range]

SB_rain = [rain_terminal_velocity_individual_SB(SB2006Vel, ρ_air, D_r) for D_r in D_r_range]
Ch_rain = [rain_terminal_velocity_individual_Chen(Chen2022.rain, ρ_air, D_r) for D_r in D_r_range]
M1_rain = [terminal_velocity_individual_1M(Blk1MVel.rain, ρ_air, D_r) for D_r in D_r_range]
Ch_snow = [snow_terminal_velocity_individual_Chen(snow, Chen2022.snow_ice, ρ_air, D_r) for D_r in D_r_range]
M1_snow = [terminal_velocity_individual_1M(Blk1MVel.snow, ρ_air, D_r) for  D_r in D_r_range]
Ch_ice = [ice_terminal_velocity_individual_Chen(Chen2022.snow_ice, ρ_air, D_r) for  D_r in D_r_range]

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
