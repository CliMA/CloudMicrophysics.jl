import Test as TT

TT.@testset "All tests" begin
    include("aerosol_activation_tests.jl")
    include("artifact_calling_tests.jl")
    include("heterogeneous_ice_nucleation_tests.jl")
    include("homogeneous_ice_nucleation_tests.jl")
    include("microphysics0M_tests.jl")
    include("microphysics1M_tests.jl")
    include("microphysics2M_tests.jl")
    include("microphysics_noneq_tests.jl")
    include("cloud_diagnostics.jl")
    include("common_functions_tests.jl")
    include("common_types_tests.jl")
    include("nucleation_unit_tests.jl")
    include("precipitation_susceptibility_tests.jl")
    include("p3_tests.jl")
    include("aqua.jl")
    include("performance_tests.jl")
    include("aerosol_activation_calibration.jl")
    include("ice_nucleation_calibration.jl")
    include("ventilation_tests.jl")
    include("DistributionTools_tests.jl")
end
nothing
