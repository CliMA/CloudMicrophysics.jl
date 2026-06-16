using Test: @testset, @test
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import ClimaParams as CP
import ForwardDiff as FD

# Direct evaluation of the analytical solution for ρ_d. Evaluated in BigFloat it is the
# high-precision reference; a Float64 reference is not accurate enough, because the direct
# expression itself loses about two percent at F_rim = 1e-7.
function get_ρ_d_direct(::Type{FT}, F_rim, ρ_rim, β_va) where {FT}
    k = (one(FT) - F_rim)^(-one(FT) / (FT(3) - β_va))
    den = (β_va - FT(2)) * (k - one(FT)) / ((one(FT) - F_rim) * k - one(FT)) - (one(FT) - F_rim)
    return ρ_rim * F_rim / den
end

@testset "P3 get_ρ_d Float32 stability (small rime fraction)" begin
    mass = CMP.ParametersP3(Float32).mass
    βva = mass.β_va
    ρ_g_ref(F_rim, ρ_rim) =
        P3.get_ρ_g(F_rim, ρ_rim, get_ρ_d_direct(BigFloat, F_rim, ρ_rim, BigFloat(βva)))
    for F_rim in (1.0f-7, 1.0f-6, 1.0f-5, 1.0f-4, 1.0f-3, 1.0f-2, 1.0f-1, 4.0f-1, 9.0f-1),
        ρ_rim in (1.0f0, 5.0f1, 4.0f2, 9.16f2)

        ρ_d = P3.get_ρ_d(mass, F_rim, ρ_rim)
        ρ_g = P3.get_ρ_g(F_rim, ρ_rim, ρ_d)
        @test ρ_g > 0           # The direct implementation produces negative numbers for some inputs.
        @test isfinite(ρ_d)
        @test ρ_g ≈ Float32(ρ_g_ref(BigFloat(F_rim), BigFloat(ρ_rim))) rtol = 1.0f-5
    end
    # The value and its derivative with respect to F_rim stay finite under ForwardDiff in Float32.
    g(x) = (ρ_d = P3.get_ρ_d(mass, x, 4.0f2); P3.get_ρ_g(x, 4.0f2, ρ_d))
    @test isfinite(FD.derivative(g, 1.0f-4))
end
