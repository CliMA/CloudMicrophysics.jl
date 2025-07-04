import Test as TT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics2M as CM2
import ClimaParams as CP

function test_ventilation_factor(FT)
    TT.@testset "Ventilation factor / P3 terminal velocity smoke test" begin
        params = CMP.ParametersP3(FT)
        vel = CMP.Chen2022VelType(FT)
        aps = CMP.AirProperties(FT)
        F_rim = FT(0.5)  # Riming fraction [-]
        ρ_rim = FT(500)  # Riming density [kg/m³]
        L_ice = FT(0.22) # Ice mass concentration [kg/m³] (not used in this test)
        N_ice = FT(1e6)  # Ice number concentration [1/m³] (not used in this test)
        ρₐ = FT(1.2)     # Air density [kg/m³]
        state = P3.get_state(params; F_rim, ρ_rim, L_ice, N_ice)
        vent = state.params.vent
        v_term = P3.ice_particle_terminal_velocity(vel, ρₐ, state)
        vent_factor = CO.ventilation_factor(vent, aps, v_term)
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        calc_vents = vent_factor.(Ds)
        smoke_vents = [0.91818553, 1.3191912, 1.7451854, 2.1598392, 2.5553002]
        for (calc_vent, smoke_vent) in zip(calc_vents, smoke_vents)
            TT.@test calc_vent ≈ smoke_vent rtol = 1e-6
        end
    end
end

TT.@testset "Ventilation Tests ($FT)" for FT in (Float64, Float32)
    test_ventilation_factor(FT)
end
nothing
