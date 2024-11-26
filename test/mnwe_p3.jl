import Test as TT
import QuadGK as QGK

import ClimaParams
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP

function mnwe(FT)

    p3 = CMP.ParametersP3(FT)

    N_s = (FT(1e6), FT(1e7), FT(1e8), FT(1e9))
    L_s = (FT(1e-4), FT(1e-3), FT(1e-2))
    ρᵣ_s = (FT(200), FT(400), FT(600), FT(800))
    Fᵣ_s = (FT(0), FT(0.5), FT(0.8), FT(0.95))
    Fₗ_s = (FT(0), FT(0.33), FT(0.67), FT(1))

    TT.@testset "Numerical integrals sanity checks for N " begin
        for N in N_s
            for L in L_s
                for ρᵣ in ρᵣ_s
                    for Fᵣ in Fᵣ_s
                        for Fₗ in Fₗ_s
                            @info("Testing for: ", N, L, ρᵣ, Fᵣ, Fₗ)

                            # Get thresholds
                            th = P3.thresholds(p3, ρᵣ, Fᵣ)
                            @info("Thresholds: ", th)

                            # Get shape parameters
                            λ, N₀ = P3.distribution_parameter_solver(
                                p3, L, N, ρᵣ, Fᵣ, Fₗ
                            )
                            @info("Shape parameters: ", λ, N₀)

                            # Get intergal bounds
                            # P3.get_ice_bound(p3, λ, N, eps(FT))
                            #@info("Integral upper bound: ", bound)

                            bound = FT(1e-2)
                            # Number and mass sanity check
                            f(d) = P3.N′ice(p3, d, λ, N₀)
                            m(d) = P3.p3_mass(p3, d, Fᵣ, Fₗ, th)

                            qgk_N, = QGK.quadgk(d -> f(d), FT(0), bound)
                            qgk_L, = QGK.quadgk(d -> f(d) * m(d), FT(0), bound)

                            @info("Sanity check: ", N, qgk_N, L, qgk_L)
                            TT.@test N ≈ qgk_N rtol = 0.05
                            TT.@test L ≈ qgk_L rtol = 0.05
                        end
                    end
                end
            end
        end
    end
end

for FT in [Float64,]# Float32]
    mnwe(FT)
end
