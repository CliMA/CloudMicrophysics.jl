import Test as TT

import ClimaParams

import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics0M as CM0

function test_microphysics0M(FT)

    TT.@testset "0M_microphysics" begin
        p0m = CMP.Parameters0M(FT)
        (; τ_precip, qc_0, S_0) = p0m

        q_vap_sat = FT(10e-3)
        qc = FT(3e-3)
        q_tot = FT(13e-3)
        frac = [FT(0), FT(0.5), FT(1.0)]

        # no rain if no cloud
        TT.@test CM0.remove_precipitation(p0m, FT(0), FT(0)) ≈ FT(0)
        TT.@test CM0.remove_precipitation(p0m, FT(0), FT(0), q_vap_sat) ≈ FT(0)


        # rain based on qc threshold
        for lf in frac
            q_lcl = qc * lf
            q_icl = (1 - lf) * qc

            TT.@test CM0.remove_precipitation(p0m, q_lcl, q_icl) ≈
                     -max(0, q_lcl + q_icl - qc_0) / τ_precip
        end

        # rain based on supersaturation threshold
        for lf in frac
            q_lcl = qc * lf
            q_icl = (1 - lf) * qc

            TT.@test CM0.remove_precipitation(p0m, q_lcl, q_icl, q_vap_sat) ≈
                     -max(0, q_lcl + q_icl - S_0 * q_vap_sat) / τ_precip
        end
    end
end

TT.@testset "Microphysics 0M Tests ($FT)" for FT in (Float64, Float32)
    test_microphysics0M(FT)
end
nothing
