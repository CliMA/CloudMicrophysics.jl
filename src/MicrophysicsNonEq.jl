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
    τ_relax(liquid)
    τ_relax(ice)

 - `liquid` or `ice` - a type for cloud liquid water or ice

Returns the relaxation timescale for condensation and evaporation of
cloud liquid water or the relaxation timescale for sublimation and
deposition of cloud ice.
"""
τ_relax(p::CT.LiquidType) = p.τ_relax
τ_relax(p::CT.IceType) = p.τ_relax

"""
    conv_q_vap_to_q_liq_ice(liquid, q_sat, q)
    conv_q_vap_to_q_liq_ice(ice, q_sat, q)

 - `liquid` or `ice` - a type for cloud water or ice
 - `q_sat` - PhasePartition at equilibrium
 - `q` - current PhasePartition

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The tendency is obtained assuming a relaxation to equilibrium with
a constant timescale.
"""
function conv_q_vap_to_q_liq_ice(
    liquid::CT.LiquidType,
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
) where {FT <: Real}
    return (q_sat.liq - q.liq) / liquid.τ_relax
end
function conv_q_vap_to_q_liq_ice(
    ice::CT.IceType,
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
) where {FT <: Real}
    return (q_sat.ice - q.ice) / ice.τ_relax
end

end #module MicrophysicsNonEq.jl
