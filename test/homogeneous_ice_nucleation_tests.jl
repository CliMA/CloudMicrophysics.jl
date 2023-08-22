import Test as TT

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

const CO = CM.Common
const CMH = CM.HomIceNucleation

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

@info "Homogeneous Ice Nucleation Tests"

function test_homogeneous_J(FT)
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    TT.@testset "Homogeneous J" begin

        T_warm = FT(229.2)
        T_cold = FT(228.8)
        x_sulph = FT(0.1)
        Δa_w_too_small = FT(0.25)
        Δa_w_too_large = FT(0.35)


        # higher nucleation rate at colder temperatures
        TT.@test CMH.homogeneous_J(prs, CO.Delta_a_w(prs, x_sulph, T_cold)) >
                 CMH.homogeneous_J(prs, CO.Delta_a_w(prs, x_sulph, T_warm))

        # If Δa_w out of range
        TT.@test_throws AssertionError("Δa_w > Δa_w_min") CMH.homogeneous_J(
            prs,
            Δa_w_too_small,
        )
        TT.@test_throws AssertionError("Δa_w < Δa_w_max") CMH.homogeneous_J(
            prs,
            Δa_w_too_large,
        )
    end
end

println("Testing Float64")
test_homogeneous_J(Float64)

println("Testing Float32")
test_homogeneous_J(Float32)
