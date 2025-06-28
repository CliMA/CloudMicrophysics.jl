import Test as TT

import ClimaParams as CP

import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.HomIceNucleation as CMH

function test_homogeneous_J_cubic(FT)

    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    H2SO4_prs = CMP.H2SO4SolutionParameters(FT)
    ip = CMP.IceNucleationParameters(FT)

    TT.@testset "Homogeneous J (cubic parameterization)" begin

        T_warm = FT(229.2)
        T_cold = FT(228.8)
        x_sulph = FT(0.1)
        Δa_w_too_small = FT(0.25)
        Δa_w_too_large = FT(0.35)

        # higher nucleation rate at colder temperatures
        TT.@test CMH.homogeneous_J_cubic(
            ip.homogeneous,
            CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_cold) -
            CO.a_w_ice(tps, T_cold),
        ) > CMH.homogeneous_J_cubic(
            ip.homogeneous,
            CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_warm) -
            CO.a_w_ice(tps, T_warm),
        )

        # If Δa_w out of range
        TT.@test_throws AssertionError("Δa_w >= ip.Δa_w_min") CMH.homogeneous_J_cubic(
            ip.homogeneous,
            Δa_w_too_small,
        )
        TT.@test_throws AssertionError("Δa_w <= ip.Δa_w_max") CMH.homogeneous_J_cubic(
            ip.homogeneous,
            Δa_w_too_large,
        )
    end
end

function test_homogeneous_J_linear(FT)

    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    H2SO4_prs = CMP.H2SO4SolutionParameters(FT)
    ip = CMP.IceNucleationParameters(FT)

    TT.@testset "Homogeneous J (linear parameterization)" begin

        T_warm = FT(229.2)
        T_cold = FT(228.8)
        x_sulph = FT(0.1)

        # higher nucleation rate at colder temperatures
        TT.@test CMH.homogeneous_J_linear(
            ip.homogeneous,
            CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_cold) -
            CO.a_w_ice(tps, T_cold),
        ) > CMH.homogeneous_J_linear(
            ip.homogeneous,
            CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_warm) -
            CO.a_w_ice(tps, T_warm),
        )

    end
end


TT.@testset "Homogeneous Ice Nucleation Tests ($FT)" for FT in (Float64, Float32)
    test_homogeneous_J_cubic(FT)
    test_homogeneous_J_linear(FT)
end
nothing
