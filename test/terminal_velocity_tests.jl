import Test as TT
import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.TerminalVelocity as TV

@info "Terminal Velocity Tests"

function test_chen_velocities(FT)

    p3 = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)

    ρ_a = FT(1.2)

    TT.@testset "Chen 2022 - Rain" begin
        Ds = range(FT(1e-6), stop = FT(1e-5), length = 5)
        expected = [0.002508, 0.009156, 0.01632, 0.02377, 0.03144]
        for i in axes(Ds, 1)
            vel = TV.velocity_chen(Ds[i], Chen2022.rain, ρ_a)
            TT.@test vel >= 0
            TT.@test vel ≈ expected[i] rtol = 1e-3
        end
    end

    TT.@testset "Chen 2022 - Ice" begin
        F_r = FT(0.5)
        ρ_r = FT(500)
        th = P3.thresholds(p3, ρ_r, F_r)
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.08109, 0.3838, 0.5917, 0.8637, 1.148]
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = TV.velocity_chen(
                D,
                Chen2022.snow_ice,
                ρ_a,
                P3.p3_mass(p3, D, F_r, th),
                P3.p3_area(p3, D, F_r, th),
                p3.ρ_i,
            )
            TT.@test vel >= 0
            TT.@test vel ≈ expected[i] rtol = 1e-3
        end

    end
end

println("Testing Float32")
test_chen_velocities(Float32)

println("Testing Float64")
test_chen_velocities(Float64)
