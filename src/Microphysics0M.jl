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
export ∂remove_precipitation_∂q_tot

"""
    remove_precipitation(params_0M::Parameters0M, q_lcl, q_icl)
    remove_precipitation(params_0M::Parameters0M, q_lcl, q_icl, q_tot, q_vap_sat)

Returns the total water tendency due to precipitation removal (kg/kg/s).

The tendency assumes relaxation with constant timescale to a state 
with total water above a threshold removed. The threshold is defined 
by either fixed condensate specific humidity (`qc_0`) or fixed 
supersaturation (`S_0`), along with the relaxation timescale, 
specified in the `Parameters0M` struct.

In the second method water is removed whenever `q_tot` exceeds 
`(1 + S_0) * q_vap_sat`, i.e. whenever the relative humidity that the 
total water would imply is above the threshold `1 + S_0`.

# Arguments
- `params_0M`: 0-moment microphysics parameters (contains `τ_precip` and either `qc_0` or `S_0`)
- `q_lcl`: cloud liquid water specific humidity (kg/kg); unused in the second method
- `q_icl`: cloud ice specific humidity (kg/kg); unused in the second method
- `q_tot`: (second method only) total water specific humidity (kg/kg)
- `q_vap_sat`: (second method only) saturation specific humidity (kg/kg)
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
    q_tot,
    q_vap_sat,
) = -max(0, (q_tot - (1 + S_0) * q_vap_sat)) / τ_precip

"""
    ∂remove_precipitation_∂q_tot(params_0M::Parameters0M, q_lcl, q_icl)
    ∂remove_precipitation_∂q_tot(params_0M::Parameters0M, q_lcl, q_icl, q_tot, q_vap_sat)

Returns the derivative of `remove_precipitation` with respect to `q_tot`.

In the 0-moment scheme `q_tot` is the only prognostic variable and
  `q_lcl` and `q_icl` are diagnosed.
The derivative is −1/τ_precip whenever precipitation is actively being removed, and 0 otherwise.

# Arguments
- `params_0M`: 0-moment microphysics parameters (contains `τ_precip` and either `qc_0` or `S_0`)
- `q_lcl`: cloud liquid water specific humidity (kg/kg); unused in the second method
- `q_icl`: cloud ice specific humidity (kg/kg); unused in the second method
- `q_tot`: (second method only) total water specific humidity (kg/kg)
- `q_vap_sat`: (second method only) saturation specific humidity (kg/kg)
"""
∂remove_precipitation_∂q_tot(
    (; τ_precip, qc_0)::CMP.Parameters0M,
    q_lcl,
    q_icl,
) = ifelse(q_lcl + q_icl > qc_0, -1 / τ_precip, zero(q_lcl))

∂remove_precipitation_∂q_tot(
    (; τ_precip, S_0)::CMP.Parameters0M,
    q_lcl,
    q_icl,
    q_tot,
    q_vap_sat,
) = ifelse(q_tot > (1 + S_0) * q_vap_sat, -1 / τ_precip, zero(q_tot))

end #module Microphysics0M.jl
