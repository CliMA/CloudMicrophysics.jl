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
import CloudMicrophysics.Parameters as CMP

export remove_precipitation

"""
    remove_precipitation(params_0M, q_liq, q_ice; q_vap_sat)
    remove_precipitation(params_0M, q; q_vap_sat)

 - `params_0M` - a struct with 0-moment parameters
 - `q` - current Thermodynamics.PhasePartition or `q_liq` and `q_ice` specific contents
 - `q_vap_sat` - specific humidity at saturation

Returns the `q_tot` tendency due to the removal of precipitation.
The tendency is obtained assuming a relaxation with a constant timescale
to a state with precipitable water removed.
The threshold for when to remove `q_tot` is defined either by the
condensate specific content or supersaturation.
The thresholds and the relaxation timescale are defined `Parameters0M` struct.
"""
remove_precipitation(
    (; τ_precip, qc_0)::CMP.Parameters0M,
    q_liq,
    q_ice,
) = -max(0, (q_liq + q_ice - qc_0)) / τ_precip
remove_precipitation(
    params::CMP.Parameters0M,
    q::TD.PhasePartition,
) = remove_precipitation(params, q.liq, q.ice)

remove_precipitation(
    (; τ_precip, S_0)::CMP.Parameters0M,
    q_liq,
    q_ice,
    q_vap_sat,
) = -max(0, (q_liq + q_ice - S_0 * q_vap_sat)) / τ_precip
remove_precipitation(
    params::CMP.Parameters0M,
    q::TD.PhasePartition,
    q_vap_sat,
) = remove_precipitation(params, q.liq, q.ice, q_vap_sat)

end #module Microphysics0M.jl
