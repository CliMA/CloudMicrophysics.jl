import Test as TT
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP

@info "P3 Scheme Tests"

function test_p3_thresholds(FT)

    p3 = CMP.ParametersP3(FT)

    TT.@testset "thresholds (nonlinear solver function)" begin

        # initialize test values:
        ρ_r = FT(400)
        F_r = FT(0.8)
        ρ_r_good = (FT(200), FT(400), FT(800)) # representative ρ_r values
        F_r_good = (FT(0.5), FT(0.8), FT(0.95)) # representative F_r values

        # test asserts
        for _ρ_r in (FT(0), FT(-1))
            TT.@test_throws AssertionError("ρ_r > FT(0)") P3.thresholds(
                p3,
                _ρ_r,
                F_r,
            )
        end
        for _ρ_r in (FT(0), FT(-1))
            TT.@test_throws AssertionError("ρ_r <= p3.ρ_l") P3.thresholds(
                p3,
                FT(1200),
                F_r,
            )
        end

        for _F_r in (FT(-1 * eps(FT)), FT(-1))
            TT.@test_throws AssertionError("F_r >= FT(0)") P3.thresholds(
                p3,
                ρ_r,
                _F_r,
            )
        end
        for _F_r in (FT(1), FT(1.5))
            TT.@test_throws AssertionError("F_r < FT(1)") P3.thresholds(
                p3,
                ρ_r,
                _F_r,
            )
        end

        # Test if the P3 scheme solution satisifies the conditions
        # from eqs. 14-17 in Morrison and Milbrandt 2015
        for F_r in F_r_good
            for ρ_r in ρ_r_good
                sol = P3.thresholds(p3, ρ_r, F_r)
                atol = 5e3 * eps(FT)
                TT.@test sol.D_cr ≈ P3.D_cr_helper(p3, F_r, sol[3]) atol = atol
                TT.@test sol.D_gr ≈ P3.D_gr_helper(p3, sol[3]) atol = atol
                TT.@test sol.ρ_g ≈ P3.ρ_g_helper(ρ_r, F_r, sol[4]) atol = atol
                TT.@test sol.ρ_d ≈ P3.ρ_d_helper(p3, sol[1], sol[2]) atol = atol
            end
        end

        # Check that the P3 scheme solution matches the published values
        function diff(ρ_r, F_r, el, gold, rtol = 2e-2)
            TT.@test P3.thresholds(p3, ρ_r, F_r)[el] * 1e3 ≈ gold rtol = rtol
        end
        # D_cr and D_gr vs Fig. 1a Morrison and Milbrandt 2015
        D_cr_fig_1a_ref = [FT(0.4946323381999426), FT(1.0170979628696817)]
        D_gr_fig_1a_ref = [FT(0.26151186272014415), FT(0.23392868352755775)]
        for val in [1, 2]
            diff(ρ_r_good[2], F_r_good[val], 1, D_cr_fig_1a_ref[val])
            diff(ρ_r_good[2], F_r_good[val], 2, D_gr_fig_1a_ref[val])
        end
        # D_cr and D_gr vs Fig. 1b Morrison and Milbrandt 2015
        #! format: off
        D_cr_fig_1b_ref = [FT(6.152144691917768), FT(3.2718818175768405), FT(1.7400778369620664)]
        D_gr_fig_1b_ref = [FT(0.39875043123651077), FT(0.2147085163169669), FT(0.11516682512848)]
        #! format: on
        for val in [1, 2, 3]
            diff(ρ_r_good[val], F_r_good[3], 1, D_cr_fig_1b_ref[val])
            diff(ρ_r_good[val], F_r_good[3], 2, D_gr_fig_1b_ref[val])
        end
    end
end

function test_p3_mass(FT)

    p3 = CMP.ParametersP3(FT)

    TT.@testset "p3 mass_ functions" begin

        # Initialize test values
        Ds = (FT(1e-5), FT(1e-4), FT(1e-3))
        ρs = (FT(400), FT(600), FT(800))
        F_rs = (FT(0.0), FT(0.5), FT(0.8))

        # Test to see that the mass_ functions give positive, not NaN values
        for D in Ds
            for ρ in ρs
                for F_r in F_rs
                    TT.@test P3.mass_s(D, ρ) >= 0
                    TT.@test P3.mass_nl(p3, D) >= 0
                    TT.@test P3.mass_r(p3, D, F_r) >= 0
                end
            end
        end
    end

    # Test to see that p3_mass gives correct mass 
    TT.@testset "p3_mass() accurate values" begin

        # Initialize test values
        Ds = (FT(1e-5), FT(1e-4), FT(1e-3))
        ρs = (FT(400), FT(600), FT(800))
        F_rs = (FT(0.0), FT(0.5), FT(0.8))
        eps = FT(1e-3)

        for ρ in ρs
            for F_r in F_rs
                D_th = P3.D_th_helper(p3)
                D1 = D_th / 2
                th = P3.thresholds(p3, ρ, F_r)

                if (F_r > 0)
                    D2 = (th.D_gr + D_th) / 2
                    D3 = (th.D_cr + th.D_gr) / 2
                    D4 = th.D_cr + eps

                    TT.@test P3.p3_mass(p3, D1, F_r, th) ==
                             P3.mass_s(D1, p3.ρ_i)
                    TT.@test P3.p3_mass(p3, D2, F_r, th) == P3.mass_nl(p3, D2)
                    TT.@test P3.p3_mass(p3, D3, F_r, th) ==
                             P3.mass_s(D3, th.ρ_g)
                    TT.@test P3.p3_mass(p3, D4, F_r, th) ==
                             P3.mass_r(p3, D4, F_r)
                else
                    D2 = D1 + eps
                    TT.@test P3.p3_mass(p3, D1, F_r, th) ==
                             P3.mass_s(D1, p3.ρ_i)
                    TT.@test P3.p3_mass(p3, D2, F_r, th) == P3.mass_nl(p3, D2)
                end

            end
        end
    end

end

println("Testing Float32")
test_p3_thresholds(Float32)
test_p3_mass(Float32)

println("Testing Float64")
test_p3_thresholds(Float64)
test_p3_mass(Float64)
