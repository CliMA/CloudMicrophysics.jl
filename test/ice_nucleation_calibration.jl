import CloudMicrophysics as CM
import Test as TT

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
include(
    joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "perfect_model_J.jl"),
)

function test_J_calibration(FT, IN_mode)
    output = calibrate_J_parameters(FT, IN_mode)
    calibrated_parameters = [output[1], output[2]]
    coeff_true = [output[3], output[4]]

    TT.@test calibrated_parameters[1] ≈ coeff_true[1] rtol = FT(0.15)
    TT.@test calibrated_parameters[2] ≈ coeff_true[2] rtol = FT(0.15)
end

@info "Ice Nucleation Test"

TT.@testset "Perfect model calibration on ABDINM" begin
    test_J_calibration(Float64, "ABDINM")
end

TT.@testset "Perfect model calibration on ABIFM" begin
    test_J_calibration(Float64, "ABIFM")
end

TT.@testset "Perfect model calibration on ABHOM" begin
    test_J_calibration(Float64, "ABHOM")
end
