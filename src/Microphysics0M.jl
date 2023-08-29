"""
    Microphysics0M

Zero-moment bulk microphysics scheme that instantly removes
moisture above certain threshold.
This is equivalent to instanteneous conversion of cloud condensate
into precipitation and precipitation fallout with infinite
terminal velocity.
"""
module Microphysics0M

import Thermodynamics as TD

import ..Parameters as CMP
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
remove_precipitation((; τ_precip, qc_0), q::TD.PhasePartition) =
    -max(0, (q.liq + q.ice - qc_0)) / τ_precip

remove_precipitation((; τ_precip, S_0), q::TD.PhasePartition, q_vap_sat) =
    -max(0, (q.liq + q.ice - S_0 * q_vap_sat)) / τ_precip

end #module Microphysics0M.jl
