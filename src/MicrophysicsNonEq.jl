"""
Non-equilibrium bulk microphysics scheme for cloud condensate formation.

Implements relaxation-to-equilibrium approach for:
- Condensation and evaporation of cloud liquid water
- Deposition and sublimation of cloud ice

See also: `Microphysics1M` for precipitating hydrometeor processes.
"""
module MicrophysicsNonEq

import ..Parameters as CMP
import ..ThermodynamicsInterface as TDI
import ..Common as CO
import ..Utilities as UT
import ..HetIceNucleation as IN

export τ_relax
export conv_q_vap_to_q_lcl
export conv_q_vap_to_q_icl
export INP_limiter
export dqcld_dT
export gamma_helper

"""
    τ_relax(ice, air_properties, frostenberg, q_icl, T)

Computes the deposition relaxation timescale through
the Frostenberg et al., (2023) parameterization.
See DOI: 10.5194/acp-23-10883-2023
"""
@inline function τ_relax(
    (; ρᵢ)::CMP.CloudIce, (; D_vapor)::CMP.AirProperties,
    ip::CMP.Frostenberg2023, q_icl, T,
)
    FT = eltype(ρᵢ)
    # Get the estimated number of INPs
    N_icl = exp(IN.INP_concentration_mean(ip, T))

    # Compute the radius assuming spherical particles and
    # mono-modal distribution
    r = N_icl > UT.ϵ_numerics(FT) ? cbrt((3 * q_icl) / (4 * FT(π) * N_icl * ρᵢ)) : zero(FT)
    r0 = FT(1e-6)
    r_safe = max(r, r0)

    # Compute the relaxation timescale
    τ = (4 * FT(π) * D_vapor * N_icl * r_safe)^(-1)
    return τ
end

"""
    INP_limiter(tendency, tps, T)

Returns `true` when ice deposition should be suppressed:
positive tendency (deposition) at T > T_freeze (no INPs available).
"""
@inline function INP_limiter(tendency, tps, T)
    return T > TDI.T_freeze(tps) && tendency > zero(tendency)
end

"""
    dqcld_dT(qᵥ_sat, L, Rᵥ, T)

Computes the derivative of the saturation specific humidity with respect to
temperature for a given phase of water.

# Arguments
- `qᵥ_sat` - saturation specific humidity [kg/kg]
- `L` - latent heat [J/kg]
- `Rᵥ` - gas constant for water vapor [J/kg/K]
- `T` - temperature [K]
"""
@inline function dqcld_dT(qᵥ_sat, L, Rᵥ, T)
    return qᵥ_sat * (L / (Rᵥ * T^2) - 1 / T)
end

"""
    gamma_helper(L, cₚ_air, dqcld_dT)

Computes the thermodynamic adjustment factor Γ.

# Arguments
- `L` - latent heat [J/kg]
- `cₚ_air` - specific heat capacity of air [J/kg/K]
- `dqcld_dT` - derivative of saturation specific humidity with respect to temperature [kg/kg/K]
"""
@inline function gamma_helper(L, cₚ_air, dqcld_dT)
    return 1 + (L / cₚ_air) * dqcld_dT
end

"""
    conv_q_vap_to_q_lcl(opt::CloudLiquidFormation, mp, tps, micro, thermo)
    conv_q_vap_to_q_lcl(::Nothing, mp, tps, micro, thermo)

Computes cloud liquid tendency from condensation and evaporation using the formulation from
Morrison & Grabowski (2008), https://doi.org/10.1175/2007JAS2374.1, and
Morrison & Milbrandt (2015), https://doi.org/10.1175/JAS-D-14-0065.1.

# Arguments
- `opt`: `CloudLiquidFormation(...)` or `nothing` (disabled)
- `mp`: 1-moment microphysics parameters
- `tps`: thermodynamics parameters
- `micro`: microphysics state `(; q_tot, q_lcl, q_icl, q_rai, q_sno)`
- `thermo`: thermodynamic state `(; ρ, T)`

# Returns
- Cloud condensate tendency [kg/kg/s]
"""
@inline conv_q_vap_to_q_lcl(::Nothing, mp, tps, micro, thermo) = zero(thermo.T)
@inline conv_q_vap_to_q_lcl(::Nothing, mp, tps, micro, thermo, p_vs_liq) = zero(thermo.T)
# Convenience method: compute the saturation vapor pressure over liquid, then
# delegate to the variant taking a precomputed `p_vs_liq`.
@inline conv_q_vap_to_q_lcl(opt::CMP.CloudLiquidFormation, mp, tps::TDI.PS, micro, thermo) =
    conv_q_vap_to_q_lcl(
        opt, mp, tps, micro, thermo,
        TDI.saturation_vapor_pressure_over_liquid(tps, thermo.T),
    )
function conv_q_vap_to_q_lcl(
    opt::CMP.CloudLiquidFormation, mp, tps::TDI.PS, micro, thermo, p_vs_liq,
)
    (; q_tot, q_lcl, q_icl, q_rai, q_sno) = micro
    (; ρ, T) = thermo
    τ = opt.τ_relax

    Rᵥ = TDI.Rᵥ(tps)
    Lᵥ = TDI.Lᵥ(tps, T)
    cₚ_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)

    qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_icl + q_sno)
    qᵥ_sat_liq = TDI.saturation_vapor_specific_content(tps, T, ρ, p_vs_liq)

    dqsl_dT = dqcld_dT(qᵥ_sat_liq, Lᵥ, Rᵥ, T)
    Γₗ = gamma_helper(Lᵥ, cₚ_air, dqsl_dT)

    sat_excess = qᵥ - qᵥ_sat_liq
    timescale = τ * Γₗ

    # compute the tendency
    return ifelse(
        sat_excess < 0,
        -min(-sat_excess, max(0, q_lcl)) / timescale,
        sat_excess / timescale,
    )
end

"""
    conv_q_vap_to_q_icl(opt::ConstantTimescale, mp, tps, micro, thermo)
    conv_q_vap_to_q_icl(opt::TemperatureDependent, mp, tps, micro, thermo)
    conv_q_vap_to_q_icl(::Nothing, mp, tps, micro, thermo)

Computes cloud ice tendency from deposition and sublimation using the formulation from
Morrison & Grabowski (2008), https://doi.org/10.1175/2007JAS2374.1, and
Morrison & Milbrandt (2015), https://doi.org/10.1175/JAS-D-14-0065.1.

# Arguments
- `opt`: `ConstantTimescale(...)`, `TemperatureDependent(...)`, or `nothing` (disabled)
- `mp`: 1-moment microphysics parameters
- `tps`: thermodynamics parameters
- `micro`: microphysics state `(; q_tot, q_lcl, q_icl, q_rai, q_sno)`
- `thermo`: thermodynamic state `(; ρ, T)`

# Returns
- Cloud condensate tendency [kg/kg/s]
"""
@inline conv_q_vap_to_q_icl(::Nothing, mp, tps, micro, thermo) = zero(thermo.T)
@inline conv_q_vap_to_q_icl(::Nothing, mp, tps, micro, thermo, p_vs_ice) = zero(thermo.T)
# Convenience method: compute the saturation vapor pressure over ice, then
# delegate to the variant taking a precomputed `p_vs_ice`.
@inline conv_q_vap_to_q_icl(opt::CMP.ConstantTimescale, mp, tps::TDI.PS, micro, thermo) =
    conv_q_vap_to_q_icl(
        opt, mp, tps, micro, thermo,
        TDI.saturation_vapor_pressure_over_ice(tps, thermo.T),
    )
function conv_q_vap_to_q_icl(
    opt::CMP.ConstantTimescale, mp, tps::TDI.PS, micro, thermo, p_vs_ice,
)
    (; q_tot, q_lcl, q_icl, q_rai, q_sno) = micro
    (; ρ, T) = thermo
    τ = opt.τ_relax

    Rᵥ = TDI.Rᵥ(tps)
    Lₛ = TDI.Lₛ(tps, T)
    cₚ_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)

    qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_icl + q_sno)
    qᵥ_sat_ice = TDI.saturation_vapor_specific_content(tps, T, ρ, p_vs_ice)

    dqsi_dT = dqcld_dT(qᵥ_sat_ice, Lₛ, Rᵥ, T)
    Γᵢ = gamma_helper(Lₛ, cₚ_air, dqsi_dT)

    sat_excess = qᵥ - qᵥ_sat_ice
    timescale = τ * Γᵢ

    # compute the tendency
    tendency = ifelse(
        sat_excess < 0,
        -min(-sat_excess, max(0, q_icl)) / timescale,
        sat_excess / timescale,
    )
    limiter = INP_limiter(tendency, tps, T)
    return ifelse(limiter, zero(tendency), tendency)
end
# `TemperatureDependent` ice formation has no precomputed-`p_vs` specialization;
# fall back to recomputing the saturation pressure internally (off the hot path).
@inline conv_q_vap_to_q_icl(
    opt::CMP.TemperatureDependent, mp, tps::TDI.PS, micro, thermo, p_vs_ice,
) = conv_q_vap_to_q_icl(opt, mp, tps, micro, thermo)
function conv_q_vap_to_q_icl(
    opt::CMP.TemperatureDependent, mp, tps::TDI.PS, micro, thermo,
)
    (; q_tot, q_lcl, q_icl, q_rai, q_sno) = micro
    (; ρ, T) = thermo
    τ_sub = opt.τ_relax
    τ_dep = τ_relax(mp.cloud.ice, mp.air_properties, opt.frostenberg, q_icl, T)

    Rᵥ = TDI.Rᵥ(tps)
    Lₛ = TDI.Lₛ(tps, T)
    cₚ_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)

    qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_icl + q_sno)
    qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)

    dqsi_dT = dqcld_dT(qᵥ_sat_ice, Lₛ, Rᵥ, T)
    Γᵢ = gamma_helper(Lₛ, cₚ_air, dqsi_dT)

    sat_excess = qᵥ - qᵥ_sat_ice
    sublimation_timescale = τ_sub * Γᵢ
    deposition_timescale = τ_dep * Γᵢ

    # compute the tendency
    tendency = ifelse(
        sat_excess < 0,
        -min(-sat_excess, max(0, q_icl)) / sublimation_timescale,
        sat_excess / deposition_timescale,
    )
    limiter = INP_limiter(tendency, tps, T)
    return ifelse(limiter, zero(tendency), tendency)
end

### -------------------------- ###
### 1-moment terminal velocity ###
### -------------------------- ###

"""
    terminal_velocity(sediment::CloudLiquid, vel::StokesRegimeVelType, ρₐ, q)
    terminal_velocity(sediment::CloudIce, vel::Chen2022VelTypeSmallIce, ρₐ, q)

Computes mass-weighted average terminal velocity for cloud droplets or ice crystals
assuming a monodisperse size distribution.

- **Cloud Liquid**: Uses Stokes Law (v ∝ D²), valid for small droplets (Re < 1).
- **Cloud Ice**: Uses Chen et al. (2022) parameterization,
  [DOI: 10.1016/j.atmosres.2022.106171](https://doi.org/10.1016/j.atmosres.2022.106171)

# Arguments
- `sediment` - cloud liquid or ice parameters struct (provides density and N_0)
- `vel` - velocity parameters (StokesRegimeVelType for liquid, Chen2022VelTypeSmallIce for ice)
- `ρₐ` - air density [kg/m³]
- `q` - cloud liquid water or ice specific content [kg/kg]

# Returns
- Mass-weighted terminal velocity [m/s]
"""
function terminal_velocity(
    (; ρw, N_0)::CMP.CloudLiquid{FT},
    vel::CMP.StokesRegimeVelType{FT},
    ρₐ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    if q > UT.ϵ_numerics(FT)
        # Stokes law: v(D) = C * D^2, valid for D < ~80 μm (Re < 1)
        v_term = CO.particle_terminal_velocity(vel, ρₐ)
        # Mean volume diameter from assumed number concentration
        D = cbrt(FT(6 / π) * ρₐ * q / N_0 / ρw)
        fall_w = v_term(D)
    end
    return fall_w
end

function terminal_velocity(
    (; ρᵢ, N_0)::CMP.CloudIce{FT},
    vel::CMP.Chen2022VelTypeSmallIce{FT},
    ρₐ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    if q > UT.ϵ_numerics(FT)
        v_term = CO.particle_terminal_velocity(vel, ρₐ, ρᵢ)
        # Mean volume diameter from assumed number concentration
        D = cbrt(FT(6 / π) * ρₐ * q / N_0 / ρᵢ)
        fall_w = max(FT(0), v_term(D))
    end
    return fall_w
end

end #module MicrophysicsNonEq.jl
