"""
    Zero-moment bulk microphysics scheme that instantly removes
    moisture above certain threshold.
    This is equivalent to instanteneous conversion of cloud condensate
    into precipitation and precipitation fallout with infinite
    terminal velocity.

"""
module Microphysics_0M

import Thermodynamics

const TD = Thermodynamics

import CloudMicrophysics.Microphysics_0M_Parameters


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
    param_set::Microphysics_0M_Parameters,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    τ_precip = param_set.τ_precip
    qc_0 = param_set.qc_0

    return -max(0, (q.liq + q.ice - qc_0)) / τ_precip
end
function remove_precipitation(
    param_set::Microphysics_0M_Parameters,
    q::TD.PhasePartition{FT},
    q_vap_sat::FT,
) where {FT <: Real}

    τ_precip = param_set.τ_precip
    S_0 = param_set.S_0

    return -max(0, (q.liq + q.ice - S_0 * q_vap_sat)) / τ_precip
end



end #module Microphysics_0M.jl
