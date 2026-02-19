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
export conv_q_vap_to_q_lcl_icl
export conv_q_vap_to_q_lcl_MM2015
export conv_q_vap_to_q_icl_MM2015
export conv_q_vap_to_q_lcl_icl_MM2015

"""
    τ_relax(liquid::CloudLiquid)
    τ_relax(ice::CloudIce)

Returns the relaxation timescale for phase change processes.

# Arguments
- `liquid::CloudLiquid` or `ice::CloudIce` - cloud condensate parameters struct

# Returns
- Relaxation timescale [s] for condensation/evaporation (liquid) or
  deposition/sublimation (ice)
"""
@inline τ_relax(p::CMP.CloudLiquid) = p.τ_relax
@inline τ_relax(p::CMP.CloudIce) = p.τ_relax

""" Calculate relaxation timescales from some number concentration N """

""" Includes options to calculate relaxation timescales from cloud droplet
number concentration for both liquid and ice. Also includes
an option to calculate ice relaxation timescale through
the Frostenberg et al., (2023). See DOI: 10.5194/acp-23-10883-2023
parameterization that approximates ice droplet number from temperature. """

function τ_N(q_c, N_c, ρ_c, D_vapor)
    r = ((3 * q_c) / (4 * pi * N_c * ρ_c))^(1 / 3)
    τ = (4 * pi * D_vapor * N_c * r)^(-1)
    return τ
end

function τ_Frostenberg(
    (; ρᵢ)::CMP.CloudIce,
    (; D_vapor)::CMP.AirProperties,
    ip::CMP.Frostenberg2023,
    q_icl,
    T,
)
    N_ice = exp(IN.INP_concentration_mean(ip, T))
    τᵢ = τ_N(q_icl, N_ice, ρᵢ, D_vapor)
    return τᵢ
end

"""
    conv_q_vap_to_q_lcl_icl(params::CloudLiquid, q_sat_liq, q_lcl)
    conv_q_vap_to_q_lcl_icl(params::CloudIce, q_sat_ice, q_icl)

Computes the tendency of cloud condensate specific content using a
simple relaxation-to-equilibrium formulation with a constant timescale.

# Arguments
- `params` - cloud liquid or ice parameters struct containing `τ_relax`
- `q_sat_liq` or `q_sat_ice` - saturation specific humidity [kg/kg]
- `q_lcl` or `q_icl` - cloud liquid water or ice specific content [kg/kg]

# Returns
- Cloud condensate tendency in kg/kg/s, positive for condensation/deposition,
  negative for evaporation/sublimation

The tendency is computed as:
```math
\\frac{dq}{dt} = \\frac{q_{sat} - q}{\\tau_{relax}}
```
"""
@inline function conv_q_vap_to_q_lcl_icl((; τ_relax)::CMP.CloudLiquid, q_sat_liq, q_lcl)
    return (q_sat_liq - q_lcl) / τ_relax
end

@inline function conv_q_vap_to_q_lcl_icl((; τ_relax)::CMP.CloudIce, q_sat_ice, q_icl)
    return (q_sat_ice - q_icl) / τ_relax
end

"""
    conv_q_vap_to_q_lcl_MM2015(tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T, τ_relax)
    conv_q_vap_to_q_icl_MM2015(tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T, τ_relax)

Computes cloud condensate tendency using the formulation from
Morrison & Grabowski (2008), https://doi.org/10.1175/2007JAS2374.1, and
Morrison & Milbrandt (2015), https://doi.org/10.1175/JAS-D-14-0065.1.

This formulation includes a thermodynamic adjustment factor Γ that
accounts for latent heat release modifying the saturation state.

# Arguments
- `tps` - thermodynamics parameters struct
- `q_tot` - total water specific content [kg/kg]
- `q_lcl` - cloud liquid water specific content [kg/kg]
- `q_icl` - cloud ice specific content [kg/kg]
- `q_rai` - rain specific content [kg/kg]
- `q_sno` - snow specific content [kg/kg]
- `ρ` - air density [kg/m³]
- `T` - air temperature [K]
- `τ_relax` - relaxation timescale [s]

# Returns
- Cloud condensate tendency [kg/kg/s]

# Notes
This function does NOT apply limiters for small or negative specific humidities.
Users should apply appropriate bounds checking when integrating in a model.
"""
function conv_q_vap_to_q_lcl_MM2015(
    tps::TDI.PS,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    ρ,
    T,
    τ_relax,
)
    Rᵥ = TDI.Rᵥ(tps)
    Lᵥ = TDI.Lᵥ(tps, T)
    cₚ_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)

    qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_icl + q_sno)

    pᵥ_sat_liq = TDI.saturation_vapor_pressure_over_liquid(tps, T)
    qᵥ_sat_liq = TDI.p2q(tps, T, ρ, pᵥ_sat_liq)

    dqsldT = qᵥ_sat_liq * (Lᵥ / (Rᵥ * T^2) - 1 / T)
    Γₗ = 1 + (Lᵥ / cₚ_air) * dqsldT

    tendency = (qᵥ - qᵥ_sat_liq) / (τ_relax * Γₗ)
    return ifelse(tendency < 0 && q_lcl <= 0, zero(tendency), tendency)
end

# just to be able to still use old functionality in ClimaAtmos -- to be deleted
function conv_q_vap_to_q_lcl_icl_MM2015(
    liquid::CMP.CloudLiquid,
    tps::TDI.PS,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    ρ,
    T,
)
    return conv_q_vap_to_q_lcl_MM2015(
        tps,
        q_tot,
        q_lcl,
        q_icl,
        q_rai,
        q_sno,
        ρ,
        T,
        liquid.τ_relax,
    )
end

function conv_q_vap_to_q_icl_MM2015(
    tps::TDI.PS,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    ρ,
    T,
    τ_relax,
)
    Rᵥ = TDI.Rᵥ(tps)
    Lₛ = TDI.Lₛ(tps, T)
    cₚ_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)

    qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_icl + q_sno)

    pᵥ_sat_ice = TDI.saturation_vapor_pressure_over_ice(tps, T)
    qᵥ_sat_ice = TDI.p2q(tps, T, ρ, pᵥ_sat_ice)

    dqsidT = qᵥ_sat_ice * (Lₛ / (Rᵥ * T^2) - 1 / T)
    Γᵢ = 1 + (Lₛ / cₚ_air) * dqsidT

    tendency = (qᵥ - qᵥ_sat_ice) / (τ_relax * Γᵢ)
    return ifelse(tendency < 0 && q_icl <= 0, zero(tendency), tendency)
end

# just to be able to still use old functionality -- to be deleted
function conv_q_vap_to_q_lcl_icl_MM2015(
    ice::CMP.CloudIce,
    tps::TDI.PS,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    ρ,
    T,
)
    return conv_q_vap_to_q_icl_MM2015(
        tps,
        q_tot,
        q_lcl,
        q_icl,
        q_rai,
        q_sno,
        ρ,
        T,
        ice.τ_relax,
    )
end

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
