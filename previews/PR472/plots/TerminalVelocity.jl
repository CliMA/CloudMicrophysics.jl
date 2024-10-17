import Plots as PL
using Measures

import ClimaParams as CP

FT = Float64

import CloudMicrophysics.MicrophysicsNonEq as CMNe
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CMO

const rain = CMP.Rain(FT)
const liquid = CMP.CloudLiquid(FT)
const ice = CMP.CloudIce(FT)
const snow = CMP.Snow(FT)

const SB2006 = CMP.SB2006(FT)
const SB2006_no_lim = CMP.SB2006(FT, false)

const Chen2022 = CMP.Chen2022VelType(FT)
const SB2006Vel = CMP.SB2006VelType(FT)
const Blk1MVel = CMP.Blk1MVelType(FT)

const oblate = CM1.Oblate()
const prolate = CM1.Prolate()

function aspect_ratio_snow_1M_oblate(snow::CMP.Snow, D::FT) where {FT <: Real}
    (; r0, m0, me, χm, Δm) = snow.mass
    (; a0, ae, Δa, χa) = snow.area
    ρᵢ = snow.ρᵢ

    aᵢ = χa * a0 * (D / 2 / r0)^(ae + Δa)
    mᵢ = χm * m0 * (D / 2 / r0)^(me + Δm)

    return 3 * sqrt(FT(π)) * mᵢ / (4 * ρᵢ * aᵢ^(3/2))
end
function aspect_ratio_snow_1M_prolate(snow::CMP.Snow, D::FT) where {FT <: Real}
    (; r0, m0, me, χm, Δm) = snow.mass
    (; a0, ae, Δa, χa) = snow.area
    ρᵢ = snow.ρᵢ

    aᵢ = χa * a0 * (D / 2 / r0)^(ae + Δa)
    mᵢ = χm * m0 * (D / 2 / r0)^(me + Δm)

    return 16 * ρᵢ^2 * aᵢ^3 / (9 * FT(π) * mᵢ^2)
end

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
    ai, bi, ci = CMO.Chen2022_vel_coeffs_B1(velo_scheme, ρ)

    v = 0
    for i in 1:3
        v += (ai[i] * (D_r)^(bi[i]) * exp(-D_r * ci[i]))
    end
    return v
end

function ice_terminal_velocity_individual_Chen(
    velo_scheme::CMP.Chen2022VelTypeSnowIce,
    ρ::FT,
    D_r::FT, #in m
) where {FT <: Real}
    ai, bi, ci = CMO.Chen2022_vel_coeffs_B2(velo_scheme, ρ)

    v = 0
    for i in 1:2
        v += (ai[i] * (D_r)^(bi[i]) * exp(-D_r * ci[i]))
    end
    return v
end

function snow_terminal_velocity_individual_Chen(
    snow::CMP.Snow,
    velo_scheme::CMP.Chen2022VelTypeSnowIce,
    ρ::FT,
    D_r::FT, #in m
) where {FT <: Real}
    ai, bi, ci = CMO.Chen2022_vel_coeffs_B4(velo_scheme, ρ)

    ϕ = FT(0.1)

    v = 0
    for i in 1:2
        v += ai[i] * (D_r)^(bi[i]) * exp(-D_r * ci[i])
    end
    v_term = ϕ^(1 / 3) * v
    return max(FT(0), v_term)
end
function snow_terminal_velocity_individual_Chen_oblate(
    snow::CMP.Snow,
    velo_scheme::CMP.Chen2022VelTypeSnowIce,
    ρ::FT,
    D_r::FT, #in m
) where {FT <: Real}
    ai, bi, ci = CMO.Chen2022_vel_coeffs_B4(velo_scheme, ρ)

    ϕ = aspect_ratio_snow_1M_oblate(snow, D_r)

    v = 0
    for i in 1:2
        v += ai[i] * (D_r)^(bi[i]) * exp(-D_r * ci[i])
    end
    v_term = ϕ^(1 / 3) * v
    return max(FT(0), v_term)
end
function snow_terminal_velocity_individual_Chen_prolate(
    snow::CMP.Snow,
    velo_scheme::CMP.Chen2022VelTypeSnowIce,
    ρ::FT,
    D_r::FT, #in m
) where {FT <: Real}
    ai, bi, ci = CMO.Chen2022_vel_coeffs_B4(velo_scheme, ρ)

    ϕ = aspect_ratio_snow_1M_prolate(snow, D_r)

    v = 0
    for i in 1:2
        v += ai[i] * (D_r)^(bi[i]) * exp(-D_r * ci[i])
    end
    v_term = ϕ^(-1 / 6) * v
    return max(FT(0), v_term)
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

ρ_air = 1.2
D_r_range = range(1e-6, stop = 6e-3, length = 1000)
D_r_range_small = range(0.1 * 1e-6, stop = 100 * 1e-6, length = 1000)
q_range = range(0, stop = 5 * 1e-3, length = 100)

#! format: off
# velocity values for cloud particle sizes
SB_rain_small = [rain_terminal_velocity_individual_SB(SB2006Vel, ρ_air, D_r)                 for D_r in D_r_range_small]
M1_rain_small = [terminal_velocity_individual_1M(Blk1MVel.rain, ρ_air, D_r)                  for D_r in D_r_range_small]
M1_snow_small = [terminal_velocity_individual_1M(Blk1MVel.snow, ρ_air, D_r)                  for D_r in D_r_range_small]
Ch_liq_small =  [rain_terminal_velocity_individual_Chen(Chen2022.rain, ρ_air, D_r)           for D_r in D_r_range_small]
Ch_ice_small =  [ice_terminal_velocity_individual_Chen(Chen2022.snow_ice, ρ_air, D_r)        for D_r in D_r_range_small]
Ch_rain_small = [rain_terminal_velocity_individual_Chen(Chen2022.rain, ρ_air, D_r)           for D_r in D_r_range_small]
Ch_snow_small = [snow_terminal_velocity_individual_Chen(snow, Chen2022.snow_ice, ρ_air, D_r) for D_r in D_r_range_small]
Ch_snow_small_oblate = [snow_terminal_velocity_individual_Chen_oblate(snow, Chen2022.snow_ice, ρ_air, D_r) for D_r in D_r_range_small]
Ch_snow_small_prolate = [snow_terminal_velocity_individual_Chen_prolate(snow, Chen2022.snow_ice, ρ_air, D_r) for D_r in D_r_range_small]
# velocity values for precip particle sizes
SB_rain = [rain_terminal_velocity_individual_SB(SB2006Vel, ρ_air, D_r)                 for D_r in D_r_range]
M1_rain = [terminal_velocity_individual_1M(Blk1MVel.rain, ρ_air, D_r)                  for D_r in D_r_range]
M1_snow = [terminal_velocity_individual_1M(Blk1MVel.snow, ρ_air, D_r)                  for D_r in D_r_range]
Ch_liq =  [rain_terminal_velocity_individual_Chen(Chen2022.rain, ρ_air, D_r)           for D_r in D_r_range]
Ch_ice =  [ice_terminal_velocity_individual_Chen(Chen2022.snow_ice, ρ_air, D_r)        for D_r in D_r_range]
Ch_rain = [rain_terminal_velocity_individual_Chen(Chen2022.rain, ρ_air, D_r)           for D_r in D_r_range]
Ch_snow  = [snow_terminal_velocity_individual_Chen(snow, Chen2022.snow_ice, ρ_air, D_r) for D_r in D_r_range]
Ch_snow_oblate  = [snow_terminal_velocity_individual_Chen_oblate(snow, Chen2022.snow_ice, ρ_air, D_r) for D_r in D_r_range]
Ch_snow_prolate = [snow_terminal_velocity_individual_Chen_prolate(snow, Chen2022.snow_ice, ρ_air, D_r) for D_r in D_r_range]
# obs data
D_Gunn_Kinzer = [0.0, 0.078, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8] .* 1e-3
u_Gunn_Kinzer = [0.0, 18.0, 27, 72, 117, 162, 206, 247, 287, 327, 367, 403, 464, 517, 565, 609, 649, 690, 727, 757, 782, 806, 826, 844, 860, 872, 883, 892, 898, 903, 907, 909, 912, 914, 916, 917] ./ 100
D_Gunn_Kinzer_small = [0.0, 0.078, 0.1] .* 1e-3
u_Gunn_Kinzer_small = [0.0, 18.0, 27] ./ 100

@info(aspect_ratio_snow_1M_oblate(snow, 317 * 1.0e-6))   # = 0.9999378568038546
@info(aspect_ratio_snow_1M_prolate(snow, 317 * 1.0e-6))  # = 1.0001242979785816

# aspect ratio plot
Aspect_Ratio_oblate = [aspect_ratio_snow_1M_oblate(snow, D_r) for  D_r in D_r_range]
Aspect_Ratio_prolate = [aspect_ratio_snow_1M_prolate(snow, D_r) for  D_r in D_r_range]
# group velocity values
bM1_rain = [CM1.terminal_velocity(rain, Blk1MVel.rain, ρ_air, q) for q in q_range]
bM1_snow = [CM1.terminal_velocity(snow, Blk1MVel.snow, ρ_air, q) for q in q_range]
bCh_rain = [CM1.terminal_velocity(rain, Chen2022.rain,     ρ_air, q) for q in q_range]
bCh_snow = [CM1.terminal_velocity(snow, Chen2022.snow_ice, ρ_air, q) for q in q_range]
bCh_snow_oblate =  [CM1.terminal_velocity(snow, Chen2022.snow_ice, ρ_air, q, oblate) for q in q_range]
bCh_snow_prolate = [CM1.terminal_velocity(snow, Chen2022.snow_ice, ρ_air, q, prolate) for q in q_range]
bCh_liq = [CMNe.terminal_velocity(liquid, Chen2022.rain,     ρ_air, q) for q in q_range]
bCh_ice = [CMNe.terminal_velocity(ice,    Chen2022.snow_ice, ρ_air, q) for q in q_range]

# individual particle comparison - cloud size range
p1 = PL.scatter(D_Gunn_Kinzer_small * 1e6, u_Gunn_Kinzer_small * 1e2, ms = 3, color = 4, label = "Gunn and Kinzer (1949)")
p1 = PL.plot!(D_r_range_small  * 1e6, M1_rain_small * 1e2, linewidth = 3, xlabel = "D [um]", ylabel = "terminal velocity [cm/s]", label = "Rain 1M",   color = :skyblue1, margin=10mm)
p1 = PL.plot!(D_r_range_small  * 1e6, M1_snow_small * 1e2, linewidth = 3, xlabel = "D [um]", ylabel = "terminal velocity [cm/s]", label = "Snow 1M",   color = :plum)
p1 = PL.plot!(D_r_range_small  * 1e6, SB_rain_small * 1e2, linewidth = 3, xlabel = "D [um]", ylabel = "terminal velocity [cm/s]", label = "Rain SB",   color = :cadetblue)
p1 = PL.plot!(D_r_range_small  * 1e6, Ch_rain_small * 1e2, linewidth = 3, xlabel = "D [um]", ylabel = "terminal velocity [cm/s]", label = "Rain Chen", color = :blue)
p1 = PL.plot!(D_r_range_small  * 1e6, Ch_snow_small         * 1e2, linewidth = 3, xlabel = "D [um]", ylabel = "terminal velocity [cm/s]", label = "Snow Chen", color = :darkviolet)
p1 = PL.plot!(D_r_range_small  * 1e6, Ch_snow_small_oblate  * 1e2, linewidth = 3, xlabel = "D [um]", ylabel = "terminal velocity [cm/s]", label = "Snow Chen", color = :darkviolet, style = :dash)
p1 = PL.plot!(D_r_range_small  * 1e6, Ch_snow_small_prolate * 1e2, linewidth = 3, xlabel = "D [um]", ylabel = "terminal velocity [cm/s]", label = "Snow Chen", color = :darkviolet, style = :dot)
p1 = PL.plot!(D_r_range_small  * 1e6, Ch_ice_small * 1e2,  linewidth = 3, xlabel = "D [um]", ylabel = "terminal velocity [cm/s]", label = "Ice Chen",  color = :orange)
# individual particle comparison - precip size range
p2 = PL.scatter(D_Gunn_Kinzer * 1e3, u_Gunn_Kinzer, ms = 3, color = 4, legend=false)
p2 = PL.plot!(D_r_range * 1e3, M1_rain, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", legend=false, color = :skyblue1)
p2 = PL.plot!(D_r_range * 1e3, M1_snow, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", legend=false, color = :plum)
p2 = PL.plot!(D_r_range * 1e3, SB_rain, linewidth = 3, xlabel = "D [um]", ylabel = "terminal velocity [m/s]", legend=false, color = :cadetblue)
p2 = PL.plot!(D_r_range * 1e3, Ch_rain, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", legend=false, color = :blue)
p2 = PL.plot!(D_r_range * 1e3, Ch_snow,         linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", legend=false, color = :darkviolet)
p2 = PL.plot!(D_r_range * 1e3, Ch_snow_oblate,  linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", legend=false, color = :darkviolet, style = :dash)
p2 = PL.plot!(D_r_range * 1e3, Ch_snow_prolate, linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", legend=false, color = :darkviolet, style = :dot)
p2 = PL.plot!(D_r_range * 1e3, Ch_ice,  linewidth = 3, xlabel = "D [mm]", ylabel = "terminal velocity [m/s]", legend=false, color = :orange)
# snow aspect ratio
p3 = PL.plot(D_r_range  * 1e3,  Aspect_Ratio_oblate,  linewidth = 3, xlabel = "D [mm]", ylabel = "aspect ratio", label = "Snow aspect ratio oblate",  color = :darkviolet, ylim = (0, 100), style = :dash)
p3 = PL.plot!(D_r_range * 1e3,  Aspect_Ratio_prolate, linewidth = 3, xlabel = "D [mm]", ylabel = "aspect ratio", label = "Snow aspect ratio prolate", color = :darkviolet, style = :dot)
# group velocity comparison
p4 = PL.plot(q_range * 1e3,  bCh_liq * 100,  linewidth = 3, xlabel = "q [g/kg]", ylabel = "terminal velocity [cm/s]", color = :green,      label="Liq Chen", margin=10mm)
p4 = PL.plot!(q_range * 1e3, bCh_ice * 100, linewidth = 3, xlabel = "q [g/kg]", ylabel = "terminal velocity [cm/s]", color = :orange,     label="Ice Chen")

p5 = PL.plot(q_range  * 1e3, bM1_rain, linewidth = 3, xlabel = "q [g/kg]", ylabel = "terminal velocity [m/s]", color = :skyblue1,   label="Rain 1M")
p5 = PL.plot!(q_range * 1e3, bM1_snow, linewidth = 3, xlabel = "q [g/kg]", ylabel = "terminal velocity [m/s]", color = :plum,       label="Snow 1M")
p5 = PL.plot!(q_range * 1e3, bCh_rain, linewidth = 3, xlabel = "q [g/kg]", ylabel = "terminal velocity [m/s]", color = :blue,       label="Rain Chen")
p5 = PL.plot!(q_range * 1e3, bCh_snow,         linewidth = 3, xlabel = "q [g/kg]", ylabel = "terminal velocity [m/s]", color = :darkviolet, label="Snow Chen")
p5 = PL.plot!(q_range * 1e3, bCh_snow_oblate,  linewidth = 3, xlabel = "q [g/kg]", ylabel = "terminal velocity [m/s]", color = :darkviolet, label="Snow Chen", style = :dash)
p5 = PL.plot!(q_range * 1e3, bCh_snow_prolate, linewidth = 3, xlabel = "q [g/kg]", ylabel = "terminal velocity [m/s]", color = :darkviolet, label="Snow Chen", style = :dot)
p5 = PL.plot!(q_range * 1e3, bCh_liq,  linewidth = 3, xlabel = "q [g/kg]", ylabel = "terminal velocity [m/s]", color = :green,      label="Liq Chen")
p5 = PL.plot!(q_range * 1e3, bCh_ice,  linewidth = 3, xlabel = "q [g/kg]", ylabel = "terminal velocity [m/s]", color = :orange,     label="Ice Chen")

PL.plot(p1, p2, p3, p4, p5, layout = (2, 3), size = (1200, 750), dpi=400)
PL.savefig("1M_individual_terminal_velocity_comparisons.svg")

#! format: on
