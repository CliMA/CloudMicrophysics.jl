import CloudMicrophysics as CM
import Test as TT

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "calibration.jl"))

function test_J_calibration(FT, IN_mode)
    params = perf_model_params(FT, IN_mode)
    IC = perf_model_IC(FT, IN_mode)

    pseudo_data = perf_model_pseudo_data(FT, IN_mode, params, IC)
    Γ = pseudo_data[2]
    y_truth = pseudo_data[1]
    coeff_true = pseudo_data[3]

    output = calibrate_J_parameters(
        FT,
        IN_mode,
        params,
        IC,
        y_truth,
        Γ,
        perfect_model = true,
    )
    calibrated_parameters = [output[1], output[2]]

    TT.@test calibrated_parameters[1] ≈ coeff_true[1] rtol = FT(0.15)
    TT.@test calibrated_parameters[2] ≈ coeff_true[2] rtol = FT(0.17)
end

@info "Ice Nucleation Calibration Test"

TT.@testset "Perfect model calibration on ABDINM" begin
    test_J_calibration(Float64, "ABDINM")
end

TT.@testset "Perfect model calibration on ABIFM" begin
    test_J_calibration(Float64, "ABIFM")
end

TT.@testset "Perfect model calibration on ABHOM" begin
    test_J_calibration(Float64, "ABHOM")
end
