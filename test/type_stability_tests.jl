using Test
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.ThermodynamicsInterface as TDI



function run_type_stability_tests()

    @testset "BulkMicrophysicsTendencies - Broadcast Type Stability" begin

        for FT in (Float32, Float64)
            println("Testing type stability for $FT...")

            tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
            N = 10

            # --- 0-Moment ---
            mp0 = CMP.Microphysics0MParams(FT)
            T_0M = fill(FT(280), N)
            q_lcl_0M = fill(FT(1e-3), N)
            q_icl_0M = fill(FT(1e-4), N)

            tendencies_0M =
                BMT.bulk_microphysics_tendencies.(
                    Ref(BMT.Microphysics0Moment()),
                    Ref(mp0),
                    Ref(tps),
                    T_0M,
                    q_lcl_0M,
                    q_icl_0M,
                )

            @test tendencies_0M isa Vector
            @test eltype(tendencies_0M) <: NamedTuple
            # strict check
            val0 = tendencies_0M[1]
            for k in keys(val0)
                @test getproperty(val0, k) isa FT
            end

            # --- 0-Moment (S_0 mode) ---
            ρ_0M = fill(FT(1.2), N)
            q_vap_sat_0M = fill(FT(0.01), N)
            tendencies_0M_S0 =
                BMT.bulk_microphysics_tendencies.(
                    Ref(BMT.Microphysics0Moment()),
                    Ref(mp0),
                    Ref(tps),
                    T_0M,
                    q_lcl_0M,
                    q_icl_0M,
                    q_vap_sat_0M,
                )

            @test tendencies_0M_S0 isa Vector
            @test eltype(tendencies_0M_S0) <: NamedTuple
            val0s = tendencies_0M_S0[1]
            for k in keys(val0s)
                @test getproperty(val0s, k) isa FT
            end

            # --- 1-Moment ---
            mp1 = CMP.Microphysics1MParams(FT)
            ρ = fill(FT(1.2), N)
            T = fill(FT(280), N)
            q_tot = fill(FT(0.012), N)
            q_lcl = fill(FT(1e-3), N)
            q_icl = fill(FT(1e-4), N)
            q_rai = fill(FT(1e-4), N)
            q_sno = fill(FT(1e-4), N)

            tendencies_1M =
                BMT.bulk_microphysics_tendencies.(
                    Ref(BMT.Microphysics1Moment()),
                    Ref(mp1),
                    Ref(tps),
                    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
                )

            @test tendencies_1M isa Vector
            val1 = tendencies_1M[1]
            for k in keys(val1)
                @test getproperty(val1, k) isa FT
            end

            # --- 2-Moment (Warm Rain) ---
            mp2_warm = CMP.Microphysics2MParams(FT; with_ice = false)

            n_lcl = fill(FT(1e8), N)
            n_rai = fill(FT(1e4), N)

            tendencies_2M_warm =
                BMT.bulk_microphysics_tendencies.(
                    Ref(BMT.Microphysics2Moment()),
                    Ref(mp2_warm),
                    Ref(tps),
                    ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai,
                )

            @test tendencies_2M_warm isa Vector
            val2w = tendencies_2M_warm[1]
            for k in keys(val2w)
                @test getproperty(val2w, k) isa FT
            end

            # --- 2-Moment (Warm + Ice P3) ---
            mp2_p3 = CMP.Microphysics2MParams(FT; with_ice = true)
            q_ice = fill(FT(1e-4), N)
            n_ice = fill(FT(1e5), N)
            q_rim = fill(FT(1e-5), N)
            b_rim = fill(FT(1e-7), N)
            logλ = fill(FT(0), N) # Not used if ice=0, but we test ice path

            tendencies_2M_p3 =
                BMT.bulk_microphysics_tendencies.(
                    Ref(BMT.Microphysics2Moment()),
                    Ref(mp2_p3),
                    Ref(tps),
                    ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai,
                    q_ice, n_ice, q_rim, b_rim, logλ,
                )

            @test tendencies_2M_p3 isa Vector
            val2p3 = tendencies_2M_p3[1]
            for k in keys(val2p3)
                @test getproperty(val2p3, k) isa FT
            end

        end
    end
end

run_type_stability_tests()
