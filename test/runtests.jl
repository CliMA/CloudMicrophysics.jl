import Test as TT
import Preferences
import ClimaParams

TT.@testset "All tests" begin
    include("aerosol_activation_tests.jl")
    include("artifact_calling_tests.jl")
    include("heterogeneous_ice_nucleation_tests.jl")
    include("homogeneous_ice_nucleation_tests.jl")
    include("microphysics0M_tests.jl")
    include("microphysics1M_tests.jl")
    include("bulk_tendencies_tests.jl")
    include("bulk_tendencies_quadrature_tests.jl")
    include("type_stability_tests.jl")
    include("microphysics2M_tests.jl")
    include("microphysics_noneq_tests.jl")
    include("cloud_diagnostics.jl")
    include("common_functions_tests.jl")
    include("common_types_tests.jl")
    include("parameter_name_maps_tests.jl")
    include("nucleation_unit_tests.jl")
    include("precipitation_susceptibility_tests.jl")
    include("p3_tests.jl")
    include("p3_tables_tests.jl")
    include("p3_rho_d_stability.jl")
    include("aqua.jl")
    include("performance_tests.jl")
    include("ventilation_tests.jl")
    include("ad_compat_tests.jl")
    include("return_type_tests.jl")
    include("DistributionTools_tests.jl")
    include("utilities_tests.jl")
    include("gamma_inc_tests.jl")
    include("unrolled_logsumexp.jl")
end

# The calibration tests construct many parameter sets, each of which is a new
# type that has to be compiled when the `CLIMA_PARAMS_MODE` preference of
# ClimaParams is "type", so they run in a separate process in "bits" mode.
TT.@testset "Calibration tests" begin
    Preferences.set_preferences!(ClimaParams, "CLIMA_PARAMS_MODE" => "bits"; force = true)
    for file in ("aerosol_activation_calibration.jl", "ice_nucleation_calibration.jl")
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(joinpath(@__DIR__, file))`
        TT.@test success(run(ignorestatus(cmd)))
    end
    Preferences.delete_preferences!(ClimaParams, "CLIMA_PARAMS_MODE"; force = true)
end
nothing
