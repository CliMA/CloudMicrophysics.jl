import CairoMakie as MK

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
const SB2006_no_lim = CMP.SB2006(FT; is_limited = false)

const Chen2022 = CMP.Chen2022VelType(FT)
const STVel = CMP.StokesRegimeVelType(FT)
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

    return 3 * sqrt(FT(π)) * mᵢ / (4 * ρᵢ * aᵢ^(3 / 2))
end
function aspect_ratio_snow_1M_prolate(snow::CMP.Snow, D::FT) where {FT <: Real}
    (; r0, m0, me, χm, Δm) = snow.mass
    (; a0, ae, Δa, χa) = snow.area
    ρᵢ = snow.ρᵢ

    aᵢ = χa * a0 * (D / 2 / r0)^(ae + Δa)
    mᵢ = χm * m0 * (D / 2 / r0)^(me + Δm)

    return 16 * ρᵢ^2 * aᵢ^3 / (9 * FT(π) * mᵢ^2)
end

function snow_terminal_velocity_individual_Chen(
    snow::CMP.Snow,
    velo_scheme::CMP.Chen2022VelTypeLargeIce,
    ρₐ::FT,
    D::FT, #in m
) where {FT <: Real}
    vₜ = CMO.particle_terminal_velocity(velo_scheme, ρₐ, snow.ρᵢ)
    (; ϕ, κ) = snow.aspr
    v_term = ϕ^κ * vₜ(D)
    return max(FT(0), v_term)
end
function snow_terminal_velocity_individual_Chen_oblate(
    snow::CMP.Snow,
    velo_scheme::CMP.Chen2022VelTypeLargeIce,
    ρₐ::FT,
    D_r::FT, #in m
) where {FT <: Real}
    vₜ = CMO.particle_terminal_velocity(velo_scheme, ρₐ, snow.ρᵢ)
    ϕ = aspect_ratio_snow_1M_oblate(snow, D_r)
    v_term = ϕ^(1 // 3) * vₜ(D_r)
    return max(FT(0), v_term)
end
function snow_terminal_velocity_individual_Chen_prolate(
    snow::CMP.Snow,
    velo_scheme::CMP.Chen2022VelTypeLargeIce,
    ρₐ::FT,
    D_r::FT, #in m
) where {FT <: Real}
    vₜ = CMO.particle_terminal_velocity(velo_scheme, ρₐ, snow.ρᵢ)
    ϕ = aspect_ratio_snow_1M_prolate(snow, D_r)
    v_term = ϕ^(-1 // 6) * vₜ(D_r)
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
N_lcl = 500e6

#! format: off
# velocity values for cloud particle sizes
v_term_rain = CMO.particle_terminal_velocity(Chen2022.rain, ρ_air)
v_term_small_ice = CMO.particle_terminal_velocity(Chen2022.small_ice, ρ_air, ice.ρᵢ)
v_term_stokes = CMO.particle_terminal_velocity(STVel, ρ_air)
ST_cloud_small = v_term_stokes.(D_r_range_small)
SB_rain_small = [rain_terminal_velocity_individual_SB(SB2006Vel, ρ_air, D_r)                 for D_r in D_r_range_small]
M1_rain_small = [terminal_velocity_individual_1M(Blk1MVel.rain, ρ_air, D_r)                  for D_r in D_r_range_small]
M1_snow_small = [terminal_velocity_individual_1M(Blk1MVel.snow, ρ_air, D_r)                  for D_r in D_r_range_small]
Ch_lcl_small = v_term_rain.(D_r_range_small)
Ch_icl_small = v_term_small_ice.(D_r_range_small)
Ch_rain_small = v_term_rain.(D_r_range_small)
Ch_snow_small = [snow_terminal_velocity_individual_Chen(snow, Chen2022.large_ice, ρ_air, D_r) for D_r in D_r_range_small]
Ch_snow_small_oblate = [snow_terminal_velocity_individual_Chen_oblate(snow, Chen2022.large_ice, ρ_air, D_r) for D_r in D_r_range_small]
Ch_snow_small_prolate = [snow_terminal_velocity_individual_Chen_prolate(snow, Chen2022.large_ice, ρ_air, D_r) for D_r in D_r_range_small]
# velocity values for precip particle sizes
ST_cloud = v_term_stokes.(D_r_range)
SB_rain = [rain_terminal_velocity_individual_SB(SB2006Vel, ρ_air, D_r)                 for D_r in D_r_range]
M1_rain = [terminal_velocity_individual_1M(Blk1MVel.rain, ρ_air, D_r)                  for D_r in D_r_range]
M1_snow = [terminal_velocity_individual_1M(Blk1MVel.snow, ρ_air, D_r)                  for D_r in D_r_range]
Ch_lcl = v_term_rain.(D_r_range)
Ch_icl = v_term_small_ice.(D_r_range)
Ch_rain = v_term_rain.(D_r_range)
Ch_snow = [snow_terminal_velocity_individual_Chen(snow, Chen2022.large_ice, ρ_air, D_r) for D_r in D_r_range]
Ch_snow_oblate  = [snow_terminal_velocity_individual_Chen_oblate(snow, Chen2022.large_ice, ρ_air, D_r) for D_r in D_r_range]
Ch_snow_prolate = [snow_terminal_velocity_individual_Chen_prolate(snow, Chen2022.large_ice, ρ_air, D_r) for D_r in D_r_range]
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
bCh_snow = [CM1.terminal_velocity(snow, Chen2022.large_ice, ρ_air, q) for q in q_range]
bCh_snow_oblate =  [CM1.terminal_velocity(snow, Chen2022.large_ice, ρ_air, q, oblate) for q in q_range]
bCh_snow_prolate = [CM1.terminal_velocity(snow, Chen2022.large_ice, ρ_air, q, prolate) for q in q_range]
bSt_N_lcl = [CM2.cloud_terminal_velocity(SB2006.pdf_c, STVel, q, ρ_air, N_lcl)[1] for q in q_range]
bSt_lcl = [CM2.cloud_terminal_velocity(SB2006.pdf_c, STVel, q, ρ_air, N_lcl)[2] for q in q_range]
bCh_lcl = [CMNe.terminal_velocity(liquid, STVel, ρ_air, q) for q in q_range]
bCh_icl = [CMNe.terminal_velocity(ice,    Chen2022.small_ice, ρ_air, q) for q in q_range]

#! format: on

fig = MK.Figure(size = (1200, 750))
ax1 = MK.Axis(
    fig[1, 1];
    xlabel = "D [μm]",
    ylabel = "terminal velocity [cm/s]",
    title = "Individual particles, cloud sizes",
)
ax2 = MK.Axis(
    fig[1, 2];
    xlabel = "D [mm]",
    ylabel = "terminal velocity [m/s]",
    title = "Individual particles, precipitation sizes",
)
ax3 = MK.Axis(fig[1, 3]; xlabel = "D [mm]", ylabel = "aspect ratio", title = "Snow aspect ratio")
ax4 = MK.Axis(
    fig[2, 1];
    xlabel = "q [g/kg]",
    ylabel = "terminal velocity [cm/s]",
    title = "Group velocity, cloud condensate",
)
ax5 =
    MK.Axis(fig[2, 2]; xlabel = "q [g/kg]", ylabel = "terminal velocity [m/s]", title = "Group velocity, precipitation")

# individual particle velocities; each entry is (label, cloud-size values, precipitation-size values, line attributes)
individual = [
    ("Rain 1M", M1_rain_small, M1_rain, (; color = :skyblue1)),
    ("Snow 1M", M1_snow_small, M1_snow, (; color = :plum)),
    ("Stokes", ST_cloud_small, ST_cloud, (; color = :cadetblue, linestyle = :dash)),
    ("Rain SB", SB_rain_small, SB_rain, (; color = :cadetblue)),
    ("Rain Chen", Ch_rain_small, Ch_rain, (; color = :blue)),
    ("Snow Chen", Ch_snow_small, Ch_snow, (; color = :darkviolet)),
    ("Snow Chen, oblate", Ch_snow_small_oblate, Ch_snow_oblate, (; color = :darkviolet, linestyle = :dash)),
    ("Snow Chen, prolate", Ch_snow_small_prolate, Ch_snow_prolate, (; color = :darkviolet, linestyle = :dot)),
    ("Ice Chen", Ch_icl_small, Ch_icl, (; color = :orange)),
]
MK.scatter!(ax1, D_Gunn_Kinzer_small * 1e6, u_Gunn_Kinzer_small * 1e2; color = :black, label = "Gunn and Kinzer (1949)")
MK.scatter!(ax2, D_Gunn_Kinzer * 1e3, u_Gunn_Kinzer; color = :black, markersize = 7)
for (label, v_small, v, attrs) in individual
    MK.lines!(ax1, D_r_range_small * 1e6, v_small * 1e2; linewidth = 3, label, attrs...)
    MK.lines!(ax2, D_r_range * 1e3, v; linewidth = 3, attrs...)
end
MK.ylims!(ax2, -1, 11)
# the empty grid cell holds the legend shared by the two individual particle panels
MK.Legend(fig[2, 3], ax1, "Individual particles"; framevisible = false, tellwidth = false)

MK.lines!(
    ax3,
    D_r_range * 1e3,
    Aspect_Ratio_oblate;
    linewidth = 3,
    color = :darkviolet,
    linestyle = :dash,
    label = "oblate",
)
MK.lines!(
    ax3,
    D_r_range * 1e3,
    Aspect_Ratio_prolate;
    linewidth = 3,
    color = :darkviolet,
    linestyle = :dot,
    label = "prolate",
)
MK.ylims!(ax3, 0, 100)
MK.axislegend(ax3; position = :rt, framevisible = false)

MK.lines!(ax4, q_range * 1e3, bSt_lcl * 100; linewidth = 3, color = :cadetblue, label = "Liq Stokes")
MK.lines!(
    ax4,
    q_range * 1e3,
    bSt_N_lcl * 100;
    linewidth = 3,
    color = :cadetblue,
    linestyle = :dash,
    label = "Liq Num Stokes",
)
MK.lines!(ax4, q_range * 1e3, bCh_lcl * 100; linewidth = 3, color = :orange, label = "Liq Chen")
MK.lines!(ax4, q_range * 1e3, bCh_icl * 100; linewidth = 3, color = :blue, label = "Ice Chen")
MK.axislegend(ax4; position = (0.95, 0.6), framevisible = false, rowgap = 0, labelsize = 13)

MK.lines!(ax5, q_range * 1e3, bM1_rain; linewidth = 3, color = :skyblue1, label = "Rain 1M")
MK.lines!(ax5, q_range * 1e3, bM1_snow; linewidth = 3, color = :plum, label = "Snow 1M")
MK.lines!(ax5, q_range * 1e3, bCh_rain; linewidth = 3, color = :blue, label = "Rain Chen")
MK.lines!(ax5, q_range * 1e3, bCh_snow; linewidth = 3, color = :darkviolet, label = "Snow Chen")
MK.lines!(
    ax5,
    q_range * 1e3,
    bCh_snow_oblate;
    linewidth = 3,
    color = :darkviolet,
    linestyle = :dash,
    label = "Snow Chen, oblate",
)
MK.lines!(
    ax5,
    q_range * 1e3,
    bCh_snow_prolate;
    linewidth = 3,
    color = :darkviolet,
    linestyle = :dot,
    label = "Snow Chen, prolate",
)
MK.lines!(ax5, q_range * 1e3, bCh_lcl; linewidth = 3, color = :green, label = "Liq Chen")
MK.lines!(ax5, q_range * 1e3, bCh_icl; linewidth = 3, color = :orange, label = "Ice Chen")
MK.axislegend(ax5; position = (0.95, 0.6), framevisible = false, rowgap = 0, labelsize = 13)

MK.save("1M_individual_terminal_velocity_comparisons.svg", fig)
nothing
