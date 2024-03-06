module EmulatorModelsExt

import MLJ
import DataFrames as DF

import Thermodynamics as TD
import Thermodynamics.Parameters as TDP
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.Parameters as CMP

function AA.N_activated_per_mode(
    machine::MLJ.Machine,
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}
    hygro = AA.mean_hygroscopicity_parameter(ap, ad)
    return ntuple(Val(AM.n_modes(ad))) do i
        # Model predicts activation of the first mode. So, swap each mode
        # with the first mode repeatedly to predict all activations.
        modes_perm = collect(1:AM.n_modes(ad))
        modes_perm[[1, i]] = modes_perm[[i, 1]]
        per_mode_data = [
            (;
                Symbol("mode_$(j)_N") => ad.Modes[modes_perm[j]].N,
                Symbol("mode_$(j)_mean") => ad.Modes[modes_perm[j]].r_dry,
                Symbol("mode_$(j)_stdev") => ad.Modes[modes_perm[j]].stdev,
                Symbol("mode_$(j)_kappa") => hygro[modes_perm[j]],
            ) for j in 1:AM.n_modes(ad)
        ]
        additional_data = (;
            :velocity => w,
            :initial_temperature => T,
            :initial_pressure => p,
        )
        X = DF.DataFrame([merge(reduce(merge, per_mode_data), additional_data)])
        max(FT(0), min(FT(1), MLJ.predict(machine, X)[1])) * ad.Modes[i].N
    end
end

end # module
