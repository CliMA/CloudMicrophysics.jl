import Test as TT

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

import CloudMicrophysics.CommonTypes as CMT
import CloudMicrophysics.Common as CO
import CloudMicrophysics.HetIceNucleation as CMI_het
const ArizonaTestDust = CMT.ArizonaTestDustType()
const DesertDust = CMT.DesertDustType()

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

@info "Heterogeneous Ice Nucleation Tests"

function test_heterogeneous_ice_nucleation(FT)

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
                TT.@test CMI_het.dust_activated_number_fraction(
                    prs,
                    Si_low,
                    T,
                    dust,
                ) == FT(0)
            end
        end

        # Activate more in cold temperatures and higher supersaturations
        for dust in [ArizonaTestDust, DesertDust]
            TT.@test CMI_het.dust_activated_number_fraction(
                prs,
                Si_hgh,
                T_warm,
                dust,
            ) > CMI_het.dust_activated_number_fraction(
                prs,
                Si_med,
                T_warm,
                dust,
            )
            TT.@test CMI_het.dust_activated_number_fraction(
                prs,
                Si_med,
                T_cold,
                dust,
            ) > CMI_het.dust_activated_number_fraction(
                prs,
                Si_med,
                T_warm,
                dust,
            )
        end

        for dust in [ArizonaTestDust, DesertDust]
            for T in [T_warm, T_cold]
                TT.@test CMI_het.dust_activated_number_fraction(
                    prs,
                    Si_too_hgh,
                    T,
                    dust,
                ) == FT(0)
            end
        end
    end

    TT.@testset "ABIFM J" begin

        T_warm_1 = FT(229.2)
        T_cold_1 = FT(228.8)
        x_sulph = FT(0.1)

        T_warm_2 = FT(285)
        T_cold_2 = FT(251)
        e_warm = FT(1088)
        e_cold = FT(544)

        # higher nucleation rate at colder temperatures
        for dust in
            [CMT.IlliteType(), CMT.KaoliniteType(), CMT.DesertDustType()]
            TT.@test CMI_het.ABIFM_J(
                prs,
                dust,
                CO.a_w_xT(prs, x_sulph, T_cold_1) - CO.a_w_ice(prs, T_cold_1),
            ) > CMI_het.ABIFM_J(
                prs,
                dust,
                CO.a_w_xT(prs, x_sulph, T_warm_1) - CO.a_w_ice(prs, T_warm_1),
            )

            TT.@test CMI_het.ABIFM_J(
                prs,
                dust,
                CO.a_w_eT(prs, e_cold, T_cold_2) - CO.a_w_ice(prs, T_cold_2),
            ) > CMI_het.ABIFM_J(
                prs,
                dust,
                CO.a_w_eT(prs, e_warm, T_warm_2) - CO.a_w_ice(prs, T_warm_2),
            )
        end
    end
end

println("Testing Float64")
test_heterogeneous_ice_nucleation(Float64)

println("Testing Float32")
test_heterogeneous_ice_nucleation(Float32)
