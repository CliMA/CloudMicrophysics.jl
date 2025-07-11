module EmulatorModelsExt

import MLJ
import DataFrames as DF

import ClimaParams as CP
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI

const TDP = TDI.TD.Parameters

"""
    N_activated_per_mode(machine, ap, ad, aip, tps, T, p, w, qₜ, qₗ, qᵢ)

  - `machine` - ML model
  - `ap`  - a struct with aerosol activation parameters
  - `ad`  - aerosol distribution struct
  - `aip` - a struct with air parameters
  - `tps` - a struct with thermodynamics parameters
  - `T`   - air temperature
  - `p`   - air pressure
  - `w`   - vertical velocity
  - `qₜ`  - total water specific content
  - `qₗ`  - liquid water specific content
  - `qᵢ`  - ice specific content

Returns the number of activated aerosol particles
in each aerosol size distribution mode by using a trained emulator.
"""
function AA.N_activated_per_mode(
    machine::MLJ.Machine,
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    T::FT,
    p::FT,
    w::FT,
    qₜ::FT,
    qₗ::FT,
    qᵢ::FT,
) where {FT <: Real}
    hygro = AA.mean_hygroscopicity_parameter(ap, ad)
    return ntuple(Val(AM.n_modes(ad))) do i
        # Model predicts activation of the first mode. So, swap each mode
        # with the first mode repeatedly to predict all activations.
        modes_perm = collect(1:AM.n_modes(ad))
        modes_perm[[1, i]] = modes_perm[[i, 1]]
        per_mode_data = [
            (;
                Symbol("mode_$(j)_N") => ad.modes[modes_perm[j]].N,
                Symbol("mode_$(j)_mean") => ad.modes[modes_perm[j]].r_dry,
                Symbol("mode_$(j)_stdev") => ad.modes[modes_perm[j]].stdev,
                Symbol("mode_$(j)_kappa") => hygro[modes_perm[j]],
            ) for j in 1:AM.n_modes(ad)
        ]
        additional_data = (;
            :velocity => w,
            :initial_temperature => T,
            :initial_pressure => p,
        )
        X = DF.DataFrame([merge(reduce(merge, per_mode_data), additional_data)])
        max(FT(0), min(FT(1), MLJ.predict(machine, X)[1])) * ad.modes[i].N
    end
end

"""
    total_N_activated(machine, ap, ad, aip, tps, T, p, w, qₜ, qₗ, qᵢ)

  - `machine` - ML model
  - `ap` - a struct with aerosol activation parameters
  - `ad` - aerosol distribution struct
  - `aip` - a struct with air properties
  - `tps` - a struct with thermodynamics parameters
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `qₜ`  - total water specific content
  - `qₗ`  - liquid water specific content
  - `qᵢ`  - ice specific content

Returns the total number of activated aerosol particles by using a trained emulator.
"""
function AA.total_N_activated(
    machine::MLJ.Machine,
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    T::FT,
    p::FT,
    w::FT,
    qₜ::FT,
    qₗ::FT,
    qᵢ::FT,
) where {FT}
    return sum(AA.N_activated_per_mode(machine, ap, ad, aip, tps, T, p, w, qₜ, qₗ, qᵢ))
end

"""
    AerosolActivationParameters(ekp_params)

    - `ekp_params` - parameters from the trained Ensemble Kalman Process
Returns a calibrated set of aerosol activation parameters
"""
function CMP.AerosolActivationParameters(ekp_params::AbstractArray{FT}) where {FT <: Real}
    default_param_set = CMP.AerosolActivationParameters(FT)
    (f1, f2, g1, g2, p1, p2) = FT.(ekp_params)
    cur_values = (;
        (
            name => getfield(default_param_set, name) for
            name in fieldnames(typeof(default_param_set))
        )...
    )
    overridden_values = merge(cur_values, (; f1, f2, g1, g2, p1, p2))
    return CMP.AerosolActivationParameters{FT}(; overridden_values...)
end

end # module
