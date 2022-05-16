"""
    Non-equilibrium bulk microphysics scheme, which includes:

  - condensation and evaporation of cloud liquid water and
    deposition and sublimation of cloud ice (relaxation to equilibrium)
"""
module MicrophysicsNonEq

import Thermodynamics
const TD = Thermodynamics

import CloudMicrophysics.CloudMicrophysicsParameters
import CloudMicrophysics.NonEqMoistureParameters

import ..CommonTypes
const CT = CommonTypes

export τ_relax
export conv_q_vap_to_q_liq_ice

"""
    τ_relax(param_set, liquid)
    τ_relax(param_set, ice)

 - `param_set` - abstract set with Earth parameters
 - `liquid` - a type for cloud liquid water
 - `ice` - a type for cloud ice

Returns the relaxation timescale for condensation and evaporation of
cloud liquid water or the relaxation timescale for sublimation and
deposition of cloud ice.
"""
τ_relax(param_set::NonEqMoistureParameters, liquid::CT.LiquidType) = param_set.τ_cond_evap
τ_relax(param_set::NonEqMoistureParameters, ice::CT.IceType) = param_set.τ_sub_dep

"""
    conv_q_vap_to_q_liq_ice(param_set, liquid, q_sat, q)
    conv_q_vap_to_q_liq_ice(param_set, ice, q_sat, q)

 - `param_set` - abstract set with Earth parameters
 - `liquid` - a type for cloud water
 - `ice` - a type for cloud ice
 - `q_sat` - PhasePartition at equilibrium
 - `q` - current PhasePartition

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The tendency is obtained assuming a relaxation to equilibrium with
a constant timescale.
"""
function conv_q_vap_to_q_liq_ice(
    param_set::NonEqMoistureParameters,
    liquid::CT.LiquidType,
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    τ_cond_evap = τ_relax(param_set, liquid)

    return (q_sat.liq - q.liq) / τ_cond_evap
end
function conv_q_vap_to_q_liq_ice(
    param_set::NonEqMoistureParameters,
    ice::CT.IceType,
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    τ_sub_dep::FT = τ_relax(param_set, ice)

    return (q_sat.ice - q.ice) / τ_sub_dep
end

end #module MicrophysicsNonEq.jl
