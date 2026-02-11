"""
    Microphysics0M

Zero-moment bulk microphysics scheme that removes cloud condensate
above a threshold, equivalent to instantaneous conversion to precipitation
with infinite terminal velocity.
"""
module Microphysics0M

import CloudMicrophysics.Parameters as CMP

import CloudMicrophysics.ThermodynamicsInterface as TDI

export remove_precipitation

"""
    remove_precipitation(params_0M::Parameters0M, q_lcl, q_icl)
    remove_precipitation(params_0M::Parameters0M, q_lcl, q_icl, q_vap_sat)

Compute the total water tendency due to precipitation removal.

The tendency assumes relaxation with constant timescale to a state
with condensate above a threshold removed. The threshold is defined
by either fixed condensate specific humidity (`qc_0`) or fixed
supersaturation (`S_0`), along with the relaxation timescale,
specified in the `Parameters0M` struct.

# Arguments
- `params_0M`: 0-moment microphysics parameters (contains `τ_precip` and either `qc_0` or `S_0`)
- `q_lcl`: cloud liquid water specific humidity (kg/kg)
- `q_icl`: cloud ice specific humidity (kg/kg)
- `q_vap_sat`: (second method only) saturation specific humidity (kg/kg)

# Returns
- Total water (`q_tot`) tendency (kg/kg/s)
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

end #module Microphysics0M.jl
