import Test as TT

import ClimaParams as CP

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
        TT.@test CM0.∂remove_precipitation_∂q_tot(p0m, FT(0), FT(0)) ≈ FT(0)
        TT.@test CM0.∂remove_precipitation_∂q_tot(p0m, FT(0), FT(0), q_vap_sat) ≈ FT(0)

        # rain based on qc threshold
        for lf in frac
            q_lcl = qc * lf
            q_icl = (1 - lf) * qc

            TT.@test CM0.remove_precipitation(p0m, q_lcl, q_icl) ≈
                     -max(0, q_lcl + q_icl - qc_0) / τ_precip
            TT.@test CM0.∂remove_precipitation_∂q_tot(p0m, q_lcl, q_icl) ≈ -1 / τ_precip

        end

        # rain based on supersaturation threshold
        for lf in frac
            q_lcl = qc * lf
            q_icl = (1 - lf) * qc

            TT.@test CM0.remove_precipitation(p0m, q_lcl, q_icl, q_vap_sat) ≈
                     -max(0, q_lcl + q_icl - S_0 * q_vap_sat) / τ_precip

            TT.@test CM0.∂remove_precipitation_∂q_tot(p0m, q_lcl, q_icl, q_vap_sat) ≈
                     -1 / τ_precip
        end
    end
end

function test_microphysics0M_params(FT)

    TT.@testset "0M_parameters" begin
        toml_dict = CP.create_toml_dict(FT)
        mp = CMP.Microphysics0MParams(toml_dict)

        TT.@test mp.precip isa CMP.Parameters0M{FT}
        TT.@test mp.precip_updraft isa CMP.Parameters0M{FT}
        TT.@test mp.evaporation isa CMP.Evaporation0M{FT}
        TT.@test mp.precip == CMP.Parameters0M(toml_dict)
        # default ratio of 1: updrafts use the grid-mean timescale
        TT.@test mp.precip_updraft.τ_precip == mp.precip.τ_precip
        TT.@test mp.precip_updraft.qc_0 == mp.precip.qc_0
        TT.@test mp.precip_updraft.S_0 == mp.precip.S_0
        TT.@test mp.evaporation == CMP.Evaporation0M(toml_dict)
        TT.@test CMP.Microphysics0MParams(FT) == mp
        TT.@test eltype(mp) == FT

        # updraft timescale is the ratio times the grid-mean timescale
        override = CP.create_toml_dict(
            FT;
            override_file = Dict(
                "precipitation_timescale_updraft_ratio" =>
                    Dict("value" => 0.2, "type" => "float"),
                "precipitation_evaporation_coefficient" =>
                    Dict("value" => 1.36e-4, "type" => "float"),
            ),
        )
        mp_o = CMP.Microphysics0MParams(override)
        TT.@test mp_o.precip_updraft.τ_precip ≈ FT(0.2) * mp.precip.τ_precip
        TT.@test mp_o.precip_updraft.qc_0 == mp.precip.qc_0
        TT.@test mp_o.precip_updraft.S_0 == mp.precip.S_0
        TT.@test mp_o.precip.τ_precip == mp.precip.τ_precip
        TT.@test mp_o.evaporation.k_E == FT(1.36e-4)

        # keyword construction without updraft or evaporation parameters
        mp_kw = CMP.Microphysics0MParams(; precip = mp.precip)
        TT.@test mp_kw.precip_updraft === mp.precip
        TT.@test isnothing(mp_kw.evaporation)
    end
end

function test_microphysics0M_evaporation(FT)

    TT.@testset "0M_evaporation" begin
        pe = CMP.Evaporation0M(FT)
        pe_on = CMP.Evaporation0M(;
            k_E = FT(5.44e-4),
            RH_c = pe.RH_c,
            a_p = pe.a_p,
            α₂ = pe.α₂,
            α₃ = pe.α₃,
        )
        (; k_E, a_p, α₂, α₃) = pe_on

        deficit = FT(2e-3)
        p_over_ps = FT(0.8)
        F = FT(1e-3)

        # off by default
        TT.@test pe.k_E == FT(0)
        TT.@test CM0.evaporate_precipitation(pe, deficit, p_over_ps, F) == FT(0)

        # no evaporation without a deficit or a flux
        TT.@test CM0.evaporate_precipitation(pe_on, FT(0), p_over_ps, F) == FT(0)
        TT.@test CM0.evaporate_precipitation(pe_on, -deficit, p_over_ps, F) == FT(0)
        TT.@test CM0.evaporate_precipitation(pe_on, deficit, p_over_ps, FT(0)) == FT(0)

        # rate value
        E = CM0.evaporate_precipitation(pe_on, deficit, p_over_ps, F)
        TT.@test E ≈ a_p * k_E * deficit * (sqrt(p_over_ps) * F / (a_p * α₂))^α₃
        TT.@test E > FT(0)
        TT.@test E isa FT

        # scaling: linear in k_E and deficit, F^α₃, (p/p_s)^(α₃/2)
        TT.@test CM0.evaporate_precipitation(pe_on, 2 * deficit, p_over_ps, F) ≈ 2 * E
        TT.@test CM0.evaporate_precipitation(pe_on, deficit, p_over_ps, 2 * F) ≈
                 2^α₃ * E
        TT.@test CM0.evaporate_precipitation(pe_on, deficit, p_over_ps / 4, F) ≈
                 (FT(1) / 2)^α₃ * E

        # type stability
        TT.@inferred CM0.evaporate_precipitation(pe_on, deficit, p_over_ps, F)
        TT.@test @allocated(CM0.evaporate_precipitation(pe_on, deficit, p_over_ps, F)) == 0
    end
end

TT.@testset "Microphysics 0M Tests ($FT)" for FT in (Float64, Float32)
    test_microphysics0M(FT)
    test_microphysics0M_params(FT)
    test_microphysics0M_evaporation(FT)
end
nothing
