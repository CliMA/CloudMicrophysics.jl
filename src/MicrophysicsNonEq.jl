"""
    Non-equilibrium bulk microphysics scheme, which includes:

  - condensation and evaporation of cloud liquid water and
    deposition and sublimation of cloud ice (relaxation to equilibrium)
"""
module MicrophysicsNonEq

import Thermodynamics as TD

import ..CommonTypes as CT
import ..Parameters as CMP

const APS = CMP.AbstractCloudMicrophysicsParameters

export τ_relax
export conv_q_vap_to_q_liq_ice

"""
    τ_relax(prs, liquid)
    τ_relax(prs, ice)

 - `prs` - abstract set with Earth parameters
 - `liquid` or `ice` - a type for cloud liquid water or ice

Returns the relaxation timescale for condensation and evaporation of
cloud liquid water or the relaxation timescale for sublimation and
deposition of cloud ice.
"""
τ_relax(prs::APS, ::CT.LiquidType) = CMP.τ_cond_evap(prs)
τ_relax(prs::APS, ::CT.IceType) = CMP.τ_sub_dep(prs)

"""
    conv_q_vap_to_q_liq_ice(prs, liquid, q_sat, q)
    conv_q_vap_to_q_liq_ice(prs, ice, q_sat, q)

 - `prs` - abstract set with Earth parameters
 - `liquid` or `ice` - a type for cloud water or ice
 - `q_sat` - PhasePartition at equilibrium
 - `q` - current PhasePartition

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The tendency is obtained assuming a relaxation to equilibrium with
a constant timescale.
"""
function conv_q_vap_to_q_liq_ice(
    prs::APS,
    liquid::CT.LiquidType,
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    _τ_cond_evap::FT = τ_relax(prs, liquid)

    return (q_sat.liq - q.liq) / _τ_cond_evap
end
function conv_q_vap_to_q_liq_ice(
    prs::APS,
    ice::CT.IceType,
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    _τ_sub_dep::FT = τ_relax(prs, ice)

    return (q_sat.ice - q.ice) / _τ_sub_dep
end

end #module MicrophysicsNonEq.jl
