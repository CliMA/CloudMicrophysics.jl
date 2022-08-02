import Test

import CloudMicrophysics
const CO = CloudMicrophysics.Common


TT.@testset "logistic_function unit tests" begin

    TT.@test CO.logistic_function(-1.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function(0.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function(1.0, 1.0, 2.0) == 0.5
    TT.@test CO.logistic_function(2.0, 1.0, 2.0) ≈ 0.9525 atol = 1e-4

    TT.@test CO.logistic_function(1.0, 0.0, 2.0) == 1.0
    TT.@test CO.logistic_function(0.0, 0.0, 2.0) == 0.0
    TT.@test_throws AssertionError("x_0 >= 0") CO.logistic_function(
        1.0,
        -1.0,
        2.0,
    )
    TT.@test_throws AssertionError("k > 0") CO.logistic_function(1.0, 1.0, 0.0)

end

TT.@testset "logistic_function_integral unit tests" begin

    TT.@test CO.logistic_function_integral(-1.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function_integral(0.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function_integral(1.0, 1.0, 2.0) ≈ 0.3115 atol = 1e-4
    TT.@test CO.logistic_function_integral(3.0, 1.0, 2.0) ≈ 2.0 atol = 1e-2

    TT.@test CO.logistic_function_integral(1.0, 0.0, 2.0) == 1.0
    TT.@test CO.logistic_function_integral(0.0, 0.0, 2.0) == 0.0
    TT.@test_throws AssertionError("x_0 >= 0") CO.logistic_function_integral(
        1.0,
        -1.0,
        2.0,
    )
    TT.@test_throws AssertionError("k > 0") CO.logistic_function_integral(
        1.0,
        1.0,
        0.0,
    )

end
