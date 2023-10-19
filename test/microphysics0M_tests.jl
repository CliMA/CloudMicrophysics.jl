import Test as TT

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics0M as CM0

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

@info "0-moment microphysics tests"

function test_microphysics0M(FT)

    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    TT.@testset "0M_microphysics" begin
        params_0m = CMP.CloudMicrophysicsParameters0M(FT)
        (; τ_precip, qc_0, S_0) = params_0m

        q_vap_sat = FT(10e-3)
        qc = FT(3e-3)
        q_tot = FT(13e-3)
        frac = [FT(0), FT(0.5), FT(1.0)]

        # no rain if no cloud
        q = TD.PhasePartition(q_tot)
        TT.@test CM0.remove_precipitation(params_0m, q) ≈ FT(0)
        TT.@test CM0.remove_precipitation(params_0m, q, q_vap_sat) ≈ FT(0)

        # rain based on qc threshold
        for lf in frac
            q_liq = qc * lf
            q_ice = (1 - lf) * qc

            q = TD.PhasePartition(q_tot, q_liq, q_ice)

            TT.@test CM0.remove_precipitation(prs, q) ≈
                     -max(0, q_liq + q_ice - qc_0) / τ_precip
        end

        # rain based on supersaturation threshold
        for lf in frac
            q_liq = qc * lf
            q_ice = (1 - lf) * qc

            q = TD.PhasePartition(q_tot, q_liq, q_ice)

            TT.@test CM0.remove_precipitation(prs, q, q_vap_sat) ≈
                     -max(0, q_liq + q_ice - S_0 * q_vap_sat) / τ_precip
        end
    end
end

println("Testing Float64")
test_microphysics0M(Float64)

println("Testing Float32")
test_microphysics0M(Float32)
