import Test as TT

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

const CMT = CM.CommonTypes
const CMI = CM.HetIceNucleation
const ArizonaTestDust = CMT.ArizonaTestDustType()
const DesertDust = CMT.DesertDustType()

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

@info "Heterogeneous Ice Nucleation Tests"

function test_dust_activation(FT)
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    TT.@testset "dust_activation" begin

        T_warm = FT(250)
        T_cold = FT(210)
        Si_low = FT(1.01)
        Si_med = FT(1.2)
        Si_hgh = FT(1.34)
        Si_too_hgh = FT(1.5)

        # No activation below critical supersaturation
        for dust in [ArizonaTestDust, DesertDust]
            for T in [T_warm, T_cold]
                TT.@test CMI.dust_activated_number_fraction(
                    prs,
                    Si_low,
                    T,
                    dust,
                ) == FT(0)
            end
        end

        # Activate more in cold temperatures and higher supersaturations
        for dust in [ArizonaTestDust, DesertDust]
            TT.@test CMI.dust_activated_number_fraction(
                prs,
                Si_hgh,
                T_warm,
                dust,
            ) > CMI.dust_activated_number_fraction(
                prs,
                Si_med,
                T_warm,
                dust,
            )
            TT.@test CMI.dust_activated_number_fraction(
                prs,
                Si_med,
                T_cold,
                dust,
            ) > CMI.dust_activated_number_fraction(
                prs,
                Si_med,
                T_warm,
                dust,
            )
        end

        for dust in [ArizonaTestDust, DesertDust]
            for T in [T_warm, T_cold]
                TT.@test CMI.dust_activated_number_fraction(
                    prs,
                    Si_too_hgh,
                    T,
                    dust,
                ) == FT(0)
            end
        end
    end
end

function test_H2SO4_soln_saturation_vapor_pressure(FT)
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    TT.@testset "H2SO4 solution saturated vapor pressure" begin

        T_warm = FT(225.0)
        T_cold = FT(200.0)
        T_too_warm = FT(240)
        T_too_cold = FT(180)
        x_sulph = FT(0.1)

        # If T out of range
        TT.@test_throws AssertionError("T < FT(235)") CMI.H2SO4_soln_saturation_vapor_pressure(
            x_sulph,
            T_too_warm,
        )
        TT.@test_throws AssertionError("T > FT(185)") CMI.H2SO4_soln_saturation_vapor_pressure(
            x_sulph,
            T_too_cold,
        )

        # p_sol should be higher at warmer temperatures
        TT.@test CMI.H2SO4_soln_saturation_vapor_pressure(x_sulph, T_warm) >
                 CMI.H2SO4_soln_saturation_vapor_pressure(x_sulph, T_cold)
    end
end

function test_ABIFM_Delta_a_w(FT)
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    TT.@testset "ABIFM Delta_a_w" begin

        T_warm = FT(229.2)
        T_cold = FT(228.8)
        x_sulph = FT(0.1)

        # Delta_a_w never greater than 1
        for T in [T_warm, T_cold]
            TT.@test CMI.ABIFM_Delta_a_w(prs, x_sulph, T) <= FT(1)
        end
    end
end

function test_ABIFM_J(FT)
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    TT.@testset "ABIFM J" begin

        T_warm = FT(229.2)
        T_cold = FT(228.8)
        x_sulph = FT(0.1)

        # higher nucleation rate at colder temperatures
        for dust in
            [CMT.IlliteType(), CMT.KaoliniteType(), CMT.DesertDustType()]
            TT.@test CMI.ABIFM_J(
                dust,
                CMI.ABIFM_Delta_a_w(prs, x_sulph, T_cold),
            ) > CMI.ABIFM_J(
                dust,
                CMI.ABIFM_Delta_a_w(prs, x_sulph, T_warm),
            )
        end
    end
end



println("Testing Float64")
test_dust_activation(Float64)


println("Testing Float32")
test_dust_activation(Float32)


println("Testing Float64")
test_H2SO4_soln_saturation_vapor_pressure(Float64)


println("Testing Float32")
test_H2SO4_soln_saturation_vapor_pressure(Float32)


println("Testing Float64")
test_ABIFM_Delta_a_w(Float64)


println("Testing Float32")
test_ABIFM_Delta_a_w(Float32)


println("Testing Float64")
test_ABIFM_J(Float64)


println("Testing Float32")
test_ABIFM_J(Float32)
