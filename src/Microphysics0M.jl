"""
    Microphysics0M

Zero-moment bulk microphysics scheme that instantly removes
moisture above certain threshold.
This is equivalent to instanteneous conversion of cloud condensate
into precipitation and precipitation fallout with infinite
terminal velocity.
"""
module Microphysics0M

import Thermodynamics
const TD = Thermodynamics

import ..Parameters
const CMP = Parameters
const APS = CMP.AbstractCloudMicrophysicsParameters

export remove_precipitation

"""
    remove_precipitation(param_set::APS, q; q_vap_sat)

 - `param_set` - abstract parameter set
 - `q` - current PhasePartition
 - `q_vap_sat` - water vapor specific humidity at saturation

Returns the `q_tot` tendency due to the removal of precipitation.
The tendency is obtained assuming a relaxation with a constant timescale
to a state with precipitable water removed.
The threshold for when to remove `q_tot` is defined either by the
condensate specific humidity or supersaturation.
The thresholds and the relaxation timescale are defined in
CLIMAParameters.
"""
function remove_precipitation(
    param_set::APS,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    _τ_precip::FT = CMP.τ_precip(param_set)
    _qc_0::FT = CMP.qc_0(param_set)

    return -max(0, (q.liq + q.ice - _qc_0)) / _τ_precip
end
function remove_precipitation(
    param_set::APS,
    q::TD.PhasePartition{FT},
    q_vap_sat::FT,
) where {FT <: Real}

    _τ_precip::FT = CMP.τ_precip(param_set)
    _S_0::FT = CMP.S_0(param_set)

    return -max(0, (q.liq + q.ice - _S_0 * q_vap_sat)) / _τ_precip
end

end #module Microphysics0M.jl
