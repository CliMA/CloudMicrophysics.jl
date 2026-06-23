using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ForwardDiff as FD

function test_activation_schemes(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    args = (tps, FT(1.05), FT(288), FT(0.012), FT(2e-3), FT(0), FT(5e7), FT(0), FT(0))
    splice(a; q_lcl = a[5], n_lcl = a[7]) = (a[1:4]..., q_lcl, a[6], n_lcl, a[8:9]...)

    @testset "NoActivation is the null source ($FT)" begin
        @test BMT.activation_source(CMP.NoActivation(), args...) === zero(FT)
    end

    @testset "DiagnosticNc relaxation ($FT)" begin
        scheme = CMP.DiagnosticNc(; N_c = FT(1e8))
        # cloudy: relax up toward the target
        s = BMT.activation_source(scheme, args...)
        @test s ≈ (FT(1e8) - FT(5e7)) / FT(60)
        # above target: relax down
        s_dn = BMT.activation_source(scheme, splice(args; n_lcl = FT(2e8))...)
        @test s_dn ≈ (FT(1e8) - FT(2e8)) / FT(60)
        # cloud-free: target zero (decay of leftover droplets)
        s0 = BMT.activation_source(scheme, splice(args; q_lcl = FT(0))...)
        @test s0 ≈ -FT(5e7) / FT(60)
        # custom timescale and threshold
        s2 = BMT.activation_source(
            CMP.DiagnosticNc(; N_c = FT(1e8), q_thresh = FT(1e-2), τ_relax = FT(30)), args...,
        )
        @test s2 ≈ -FT(5e7) / FT(30)  # q_lcl below the raised threshold
    end

    @testset "activation in the fused 2M tendency ($FT)" begin
        mp_on = CMP.Microphysics2MParams(
            FT; activation_scheme = CMP.DiagnosticNc(; N_c = FT(1e8)),
        )
        mp_off = CMP.Microphysics2MParams(FT)
        @test mp_off.warm_rain.activation_scheme isa CMP.NoActivation
        t_on = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(), mp_on, tps,
            FT(1.05), FT(288), FT(0.012), FT(2e-3), FT(5e7), FT(1e-4), FT(4e4),
        )
        t_off = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(), mp_off, tps,
            FT(1.05), FT(288), FT(0.012), FT(2e-3), FT(5e7), FT(1e-4), FT(4e4),
        )
        @test t_off.dn_lcl_activation_dt === zero(FT)
        @test t_on.dn_lcl_activation_dt ≈ (FT(1e8) - FT(5e7)) / FT(60)
        # the activation source is folded into the total droplet tendency
        @test t_on.dn_lcl_dt - t_off.dn_lcl_dt ≈ t_on.dn_lcl_activation_dt
    end

    @testset "activation source is concretely typed under mixed args ($FT)" begin
        D = FD.Dual{Nothing, FT, 8}
        scheme = CMP.DiagnosticNc(; N_c = FT(1e8))
        for sig in (
            (FT, FT, FT, D, FT, D, FT, FT),  # Dual q_lcl, n_lcl
            (FT, FT, FT, FT, FT, D, FT, FT), # Dual n_lcl only
            (FT, FT, FT, D, FT, FT, FT, FT), # Dual q_lcl only
        )
            rts = Base.return_types(
                BMT.activation_source, Tuple{typeof(scheme), typeof(tps), FT, sig...},
            )
            @test all(isconcretetype, rts)
        end
    end
end

test_activation_schemes(Float64)
test_activation_schemes(Float32)
