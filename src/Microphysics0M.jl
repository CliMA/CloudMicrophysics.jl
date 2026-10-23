"""
    Microphysics0M

Zero-moment bulk microphysics scheme that removes cloud condensate 
above a threshold, equivalent to instantaneous conversion to precipitation 
with infinite terminal velocity. The resulting precipitation flux can
evaporate in sub-saturated air below cloud.
"""
module Microphysics0M

import CloudMicrophysics.Parameters as CMP

import CloudMicrophysics.ThermodynamicsInterface as TDI

export remove_precipitation
export ∂remove_precipitation_∂q_tot
export evaporate_precipitation

"""
    remove_precipitation(params_0M::Parameters0M, q_lcl, q_icl)
    remove_precipitation(params_0M::Parameters0M, q_lcl, q_icl, q_vap_sat)

Returns the total water tendency due to precipitation removal (kg/kg/s).

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

"""
    ∂remove_precipitation_∂q_tot(params_0M::Parameters0M, q_lcl, q_icl)
    ∂remove_precipitation_∂q_tot(params_0M::Parameters0M, q_lcl, q_icl, q_vap_sat)

Returns the derivative of `remove_precipitation` with respect to `q_tot`.

In the 0-moment scheme `q_tot` is the only prognostic variable and
  `q_lcl` and `q_icl` are diagnosed.
The derivative is −1/τ_precip whenever precipitation is actively being removed, and 0 otherwise.

# Arguments
- `params_0M`: 0-moment microphysics parameters (contains `τ_precip` and either `qc_0` or `S_0`)
- `q_lcl`: cloud liquid water specific humidity (kg/kg)
- `q_icl`: cloud ice specific humidity (kg/kg)
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
    q_vap_sat,
) = ifelse(q_lcl + q_icl > S_0 * q_vap_sat, -1 / τ_precip, zero(q_lcl))

"""
    evaporate_precipitation(params::Evaporation0M, deficit, p_over_ps, F)

Returns the evaporation (or sublimation) rate of the precipitation flux
in sub-saturated air (kg/kg/s, positive for a vapor source),
following [Kessler1969](@cite) as used in the IFS ([Tiedtke1993](@cite)):

```math
E = a_p k_E \\max(0, \\mathrm{deficit})
    \\left[\\frac{\\sqrt{p/p_s} \\, F / a_p}{α_2}\\right]^{α_3}
```

The rate is not limited by the available precipitation or by the vapor
needed to saturate the air; time-stepping schemes should apply those limits.

# Arguments
- `params`: 0-moment evaporation parameters (contains `k_E`, `a_p`, `α₂`, `α₃`)
- `deficit`: humidity deficit `RH_c q_vap_sat - q_vap` (kg/kg), with `RH_c` from
  `params` and `q_vap_sat` over liquid for rain and over ice for snow
- `p_over_ps`: ratio of air pressure to surface pressure (≥ 0) (-)
- `F`: grid-mean precipitation flux, positive downward (≥ 0) (kg/m2/s)
"""
evaporate_precipitation(
    (; k_E, a_p, α₂, α₃)::CMP.Evaporation0M,
    deficit,
    p_over_ps,
    F,
) = a_p * k_E * max(0, deficit) * (sqrt(p_over_ps) * F / (a_p * α₂))^α₃

end #module Microphysics0M.jl
