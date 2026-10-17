import Plots as PL

import CloudMicrophysics as CM
import CLIMAParameters as CP

FT = Float64

import CloudMicrophysics.CommonTypes as CMT
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CMO

const rain = CMT.RainType()
const liquid = CMT.LiquidType()
const ice = CMT.IceType()
const snow = CMT.SnowType()
const SB2006 = CMT.SB2006Type()
const Chen2022 = CMT.Chen2022Type()
const SB2006Vel = CMT.SB2006VelType()
const Blk1MVel = CMT.Blk1MVelType()
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
    param_set,
    ρ::FT,
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
    prs::APS,
    precip::CMT.AbstractPrecipType,
    ρ::FT,
    D_r::FT,
) where {FT <: Real}
    _χv::FT = CM1.χv(prs, precip)
    _v0::FT = CM1.v0(prs, ρ, precip)
    _ve::FT = CM1.ve(prs, precip)
    _r0::FT = CM1.r0(prs, precip)
    _Δv::FT = CM1.Δv(prs, precip)
    vt = _χv * _v0 * (D_r / (2 * _r0))^(_Δv + _ve)
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
    prs::APS,
    ρ::FT,
    D_r::FT, #in m
) where {FT <: Real}

    ρ_i::FT = CMP.ρ_cloud_ice(prs)

    _As, _Bs, _Cs, _Es, _Fs, _Gs = CMO.Chen2022_snow_ice_coeffs(prs, ρ_i)

    ai = [_Es * ρ^_As, _Fs * ρ^_As]
    bi = [_Bs + ρ * _Cs, _Bs + ρ * _Cs]
    ci = [0, _Gs]

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
    prs::APS,
    ρ::FT,
    D_r::FT,
) where {FT <: Real}
    D_r = D_r * 1000
    _m0::FT = CM1.m0(prs, snow)
    _me::FT = CM1.me(prs, snow)
    _Δm::FT = CM1.Δm(prs, snow)
    _a0::FT = CM1.a0(prs, snow)
    _ae::FT = CM1.ae(prs, snow)
    _Δa::FT = CM1.Δa(prs, snow)
    _χm::FT = CM1.χm(prs, snow)
    _χa::FT = CM1.χa(prs, snow)
    _r0::FT = CM1.r0(prs, snow)
    ρ_i::FT = CMP.ρ_cloud_ice(prs)
    m0_comb::FT = _m0 * _χm
    a0_comb::FT = _a0 * _χa
    me_comb::FT = _me + _Δm
    ae_comb::FT = _ae + _Δa
    α = FT(-1 / 3)

    _As, _Bs, _Cs, _Es, _Fs, _Gs = CMO.Chen2022_snow_ice_coeffs(prs, ρ_i)

    ai = [_Es * ρ^_As, _Fs * ρ^_As]
    bi = [_Bs + ρ * _Cs, _Bs + ρ * _Cs]
    ci = [0, _Gs]

    aspect_ratio =
        ((16 * (a0_comb)^3 * (ρ_i)^2) / (9 * π * (m0_comb)^2)) *
        (D_r / (2 * _r0 * 1000))^(3 * ae_comb - 2 * me_comb)
    aspect_ratio = aspect_ratio^(α)
    vt = 0
    for i in 1:2
        vt = vt + aspect_ratio * ai[i] * ((D_r))^(bi[i]) * exp(-1 * ci[i] * D_r)
    end
    return vt
end

function aspect_ratio_snow_1M(prs::APS, D_r::FT) where {FT <: Real}
    _m0::FT = CM1.m0(prs, snow)
    _me::FT = CM1.me(prs, snow)
    _Δm::FT = CM1.Δm(prs, snow)
    _a0::FT = CM1.a0(prs, snow)
    _ae::FT = CM1.ae(prs, snow)
    _Δa::FT = CM1.Δa(prs, snow)
    _χm::FT = CM1.χm(prs, snow)
    _χa::FT = CM1.χa(prs, snow)
    _r0::FT = CM1.r0(prs, snow)
    ρ_i::FT = CMP.ρ_cloud_ice(prs)

    m0_comb::FT = _m0 * _χm
    a0_comb::FT = _a0 * _χa
    me_comb::FT = _me + _Δm
    ae_comb::FT = _ae + _Δa

    aspect_ratio =
        ((16 * (a0_comb)^3 * (ρ_i)^2) / (9 * π * (m0_comb)^2)) *
        (D_r / (2 * _r0))^(3 * ae_comb - 2 * me_comb)
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

SB_rain_bN = [CM2.rain_terminal_velocity(param_set, SB2006, SB2006Vel,  q_rai, ρ_air, N_rai)[1] for q_rai in q_rain_range]
Ch_rain_bN = [CM2.rain_terminal_velocity(param_set, SB2006, Chen2022, q_rai, ρ_air, N_rai)[1] for q_rai in q_rain_range]
SB_rain_bM = [CM2.rain_terminal_velocity(param_set, SB2006, SB2006Vel,  q_rai, ρ_air, N_rai)[2] for q_rai in q_rain_range]
Ch_rain_bM = [CM2.rain_terminal_velocity(param_set, SB2006, Chen2022, q_rai, ρ_air, N_rai)[2] for q_rai in q_rain_range]
M1_rain_bM = [CM1.terminal_velocity(param_set, rain, Blk1MVel, ρ_air, q_rai) for q_rai in q_rain_range]

SB_rain = [rain_terminal_velocity_individual_SB(param_set, ρ_air, D_r) for D_r in D_r_range]
Ch_rain = [rain_terminal_velocity_individual_Chen(param_set, ρ_air, D_r) for D_r in D_r_range]
M1_rain = [terminal_velocity_individual_1M(param_set, rain, ρ_air, D_r) for D_r in D_r_range]
Ch_snow = [snow_terminal_velocity_individual_Chen(param_set, ρ_air, D_r) for D_r in D_r_range]
M1_snow = [terminal_velocity_individual_1M(param_set, snow, ρ_air, D_r) for  D_r in D_r_range]
Ch_ice = [ice_terminal_velocity_individual_Chen(param_set, ρ_air, D_r) for  D_r in D_r_range]

Aspect_Ratio = [aspect_ratio_snow_1M(param_set, D_r) for  D_r in D_r_ar_range]

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
