import Test as TT
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP

@info "P3 Scheme Tests"

function test_p3_thresholds(FT)

    include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    # bulk density of ice:
    ρ_i::FT = CMP.ρ_cloud_ice(prs)
    # mass power law coefficient and exponent:
    α_va::FT = P3.α_va(prs)
    β_va = CMP.β_va_BF1995(prs)
    # threshold particle dimension
    D_th::FT = P3.D_th(prs, FT)

    TT.@testset "thresholds (nonlinear solver function)" begin

        # initialize test values:
        F_r_bad = [FT(0.0), FT(-1.0), FT(1.0), FT(1.5)] # unreasonable ("bad") values
        ρ_r_bad = [FT(0.0), FT(-1.0), FT(1200.0)] # unreasonable ("bad") values
        ρ_r_good = (FT(200), FT(400), FT(800)) # representative ρ_r values
        F_r_good = (FT(0.5), FT(0.8), FT(0.95)) # representative F_r values

        # If no rime present:
        for ρ_r in ρ_r_good
            TT.@test_throws DomainError(
                F_r_bad[1],
                "D_cr, D_gr, ρ_g, ρ_d are not physically relevant when no rime is present.",
            ) P3.thresholds(ρ_r, F_r_bad[1])
        end

        for F_r in F_r_good
            TT.@test_throws DomainError(
                ρ_r_bad[1],
                "D_cr, D_gr, ρ_g, ρ_d are not physically relevant when no rime is present.",
            ) P3.thresholds(ρ_r_bad[1], F_r)
        end

        # If unreasonably large values:
        for ρ_r in ρ_r_good
            for F_r in F_r_bad[3:4]
                TT.@test_throws DomainError(
                    F_r,
                    "The rime mass fraction F_r is not physically defined for values greater than or equal to 1 because some fraction of the total mass must always consist of the mass of the unrimed portion of the particle.",
                ) P3.thresholds(ρ_r, F_r)
            end
        end

        for F_r in F_r_good
            TT.@test_throws DomainError(
                ρ_r_bad[3],
                "Predicted rime density ρ_r, being a density of bulk ice, cannot exceed the density of water (997 kg m^-3).",
            ) P3.thresholds(ρ_r_bad[3], F_r)
        end

        # If negative values:
        for ρ_r in ρ_r_good
            TT.@test_throws DomainError(
                F_r_bad[2],
                "Rime mass fraction F_r cannot be negative.",
            ) P3.thresholds(ρ_r, F_r_bad[2])
        end

        for F_r in F_r_good
            TT.@test_throws DomainError(
                ρ_r_bad[2],
                "Predicted rime density ρ_r cannot be negative.",
            ) P3.thresholds(ρ_r_bad[2], F_r)
        end

        # Is the result consistent with the expressions for D_cr, D_gr, ρ_g, ρ_d?
        # Define function:
        function f(u, p)
            return [
                (u[1]) -
                (
                    (1 / (1 - p[2])) * ((6 * α_va) / (FT(π) * u[3]))
                )^(1 / (3 - β_va)),
                (u[2]) - (((6 * α_va) / (FT(π) * (u[3])))^(1 / (3 - β_va))),
                (u[3]) - (p[1] * p[2]) - ((1 - p[2]) * (u[4])),
                (u[4]) - (
                    ((6 * α_va) * ((u[1]^(β_va - 2)) - ((u[2])^(β_va - 2)))) /
                    (FT(π) * (β_va - 2) * (max((u[1]) - (u[2]), 1e-16)))
                ),
            ]
        end

        # Test for all "good" values if passing the output back
        # into the function gives 0, with tolerance 1.5e-6
        for F_r in F_r_good
            for ρ_r in ρ_r_good
                p = [ρ_r, F_r]
                vals = P3.thresholds(ρ_r, F_r)
                output = f(vals, p)
                for result in output
                    TT.@test abs(result) ≈ FT(0) atol = 1.1e3 * eps(FT)
                end
            end
        end

        # Define diff function for comparisons
        function diff(test, gold, delta = 2e-2)
            TT.@test test ≈ gold rtol = delta
        end

        # Check that D_cr returned by P3.thresholds() matches the value
        # displayed in Fig. 1a of Morrison and Milbrandt 2015 within 1% error:
        # MM2015 values against which we test are obtained with use of
        # WebPlotDigitizer (https://automeris.io/WebPlotDigitizer/)

        diff(
            FT(P3.thresholds(ρ_r_good[2], F_r_good[1])[1] * 1e3),
            FT(0.4946323381999426),
        )
        diff(
            FT(P3.thresholds(ρ_r_good[2], F_r_good[2])[1] * 1e3),
            FT(1.0170979628696817),
        )

        # same for D_gr:
        diff(
            FT(P3.thresholds(ρ_r_good[2], F_r_good[1])[2] * 1e3),
            FT(0.26151186272014415),
        )
        diff(
            FT(P3.thresholds(ρ_r_good[2], F_r_good[2])[2] * 1e3),
            FT(0.23392868352755775),
        )

        # Similarly, check that D_cr, D_gr returned by P3.thresholds()
        # matches the value displayed in Fig. 1b of MM2015 within 1% error
        # D_cr:
        diff(
            FT(P3.thresholds(ρ_r_good[1], F_r_good[3])[1] * 1e3),
            FT(6.152144691917768),
        )
        diff(
            FT(P3.thresholds(ρ_r_good[2], F_r_good[3])[1] * 1e3),
            FT(3.2718818175768405),
        )
        diff(
            FT(P3.thresholds(ρ_r_good[3], F_r_good[3])[1] * 1e3),
            FT(1.7400778369620664),
        )

        # D_gr
        diff(
            FT(P3.thresholds(ρ_r_good[1], F_r_good[3])[2] * 1e3),
            FT(0.39875043123651077),
        )
        diff(
            FT(P3.thresholds(ρ_r_good[2], F_r_good[3])[2] * 1e3),
            FT(0.2147085163169669),
        )
        diff(
            FT(P3.thresholds(ρ_r_good[3], F_r_good[3])[2] * 1e3),
            FT(0.11516682512848),
        )

    end
end

println("Testing Float32")
test_p3_thresholds(Float32)

println("Testing Float64")
test_p3_thresholds(Float64)
