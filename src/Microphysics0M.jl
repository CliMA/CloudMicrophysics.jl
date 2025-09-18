"""
    Microphysics0M

Zero-moment bulk microphysics scheme that instantly removes
moisture above certain threshold.
This is equivalent to instanteneous conversion of cloud condensate
into precipitation and precipitation fallout with infinite
terminal velocity.
"""
module Microphysics0M

import CloudMicrophysics.Parameters as CMP

# Needed for the wrapper for calling the 0-moment remove_precipitation
# with TD.PhasePartition in ClimaAtmos
import CloudMicrophysics.ThermodynamicsInterface as TDI

export remove_precipitation

"""
    remove_precipitation(params_0M, q_lcl, q_icl; q_vap_sat)
    remove_precipitation(params_0M, q; q_vap_sat)

 - `params_0M` - a struct with 0-moment parameters
 - `q` - Thermodynamics PhasePartition struct (deprecated)
 - `q_lcl` and `q_icl` cloud liquid water and cloud ice specific contents
 - `q_vap_sat` - specific humidity at saturation

Returns the tota water (`q_tot`) tendency due to the removal of precipitation.
The tendency is obtained assuming a relaxation with a constant timescale
to a state with precipitable water removed.
The threshold for when to remove `q_tot` is defined either by the
condensate specific content or supersaturation.
The thresholds and the relaxation timescale are defined `Parameters0M` struct.
"""
remove_precipitation(
    (; τ_precip, qc_0)::CMP.Parameters0M,
    q_lcl,
    q_icl,
) = -max(0, (q_lcl + q_icl - qc_0)) / τ_precip
remove_precipitation(
    (; τ_precip, S_0)::CMP.Parameters0M,
    q_lcl,
    q_icl,
    q_vap_sat,
) = -max(0, (q_lcl + q_icl - S_0 * q_vap_sat)) / τ_precip

###
### Wrappers for calling with TD.PhasePartition
###
### For now leaving the PhasePartition wrapper because I'm not sure how to get
### rid of equilibrium thermo state in the Atmos model.
###
remove_precipitation(params::CMP.Parameters0M, q::TDI.TD.PhasePartition) =
    remove_precipitation(params, q.liq, q.ice)
remove_precipitation(params::CMP.Parameters0M, q::TDI.TD.PhasePartition, q_vap_sat) =
    remove_precipitation(params, q.liq, q.ice, q_vap_sat)

end #module Microphysics0M.jl
