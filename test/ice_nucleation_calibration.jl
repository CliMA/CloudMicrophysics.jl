import CloudMicrophysics as CM
import Test as TT

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "calibration.jl"))
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "calibration_setup.jl"))

function test_J_calibration(FT, IN_mode)
    params = perf_model_params(FT, IN_mode)
    IC = perf_model_IC(FT, IN_mode)
    end_sim = 25

    pseudo_data = perf_model_pseudo_data(FT, IN_mode, [params], [IC], end_sim)
    Γ = pseudo_data[2]
    y_truth = pseudo_data[1]
    coeff_true = pseudo_data[3]

    EKI_output = calibrate_J_parameters_EKI(
        FT,
        IN_mode,
        [params],
        [IC],
        y_truth,
        end_sim,
        Γ,
        perfect_model = true,
    )
    UKI_output = calibrate_J_parameters_UKI(
        FT,
        IN_mode,
        [params],
        [IC],
        y_truth,
        end_sim,
        Γ,
        perfect_model = true,
    )

    EKI_calibrated_parameters = EKI_output[1]
    UKI_calibrated_parameters = UKI_output[1]
    EKI_calibrated_soln = run_model(
        [params],
        IN_mode,
        EKI_calibrated_parameters,
        FT,
        [IC],
        end_sim,
    )
    UKI_calibrated_soln = run_model(
        [params],
        IN_mode,
        UKI_calibrated_parameters,
        FT,
        [IC],
        end_sim,
    )
    true_soln = run_model([params], IN_mode, coeff_true, FT, [IC], end_sim)

    # test that end ICNC are similar
    TT.@testset "EKI Perfect Model Calibrations on AIDA" begin
        TT.@test_skip EKI_calibrated_soln[9, end] ≈ true_soln[9, end] rtol = FT(0.3)
    end
    TT.@testset "UKI Perfect Model Calibrations on AIDA" begin
        TT.@test_skip UKI_calibrated_soln[9, end] ≈ true_soln[9, end] rtol = FT(0.3)
    end
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
