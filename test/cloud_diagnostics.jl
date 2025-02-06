import Test as TT

import CloudMicrophysics as CM
import ClimaParams as CP

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.CloudDiagnostics as CMD

@info "Cloud diagnostics tests"

function test_cloud_diagnostics(FT)

    # Seifert and Beheng 2006 parameters
    override_file = joinpath(
        pkgdir(CM),
        "src",
        "parameters",
        "toml",
        "SB2006_limiters.toml",
    )
    toml_dict = CP.create_toml_dict(FT; override_file)
    SB2006 = CMP.SB2006(toml_dict)
    SB2006_no_limiters = CMP.SB2006(toml_dict, false)

    # Water parameters
    wtr = CMP.WaterProperties(FT)
    rain = CMP.Rain(FT)
    cloud_liquid = CMP.CloudLiquid(FT)
    cloud_ice = CMP.CloudIce(FT)

    TT.@testset "1M microphysics RadarReflectivity" begin

        # some example values
        ρ_air = FT(1)
        q_rai = FT(0.18e-3)

        TT.@test CMD.radar_reflectivity_1M(rain, q_rai, ρ_air) ≈ FT(12.17) atol =
            0.2

        q_rai = FT(0.89e-4)

        TT.@test CMD.radar_reflectivity_1M(rain, q_rai, ρ_air) ≈ FT(6.68) atol =
            0.2

    end

    TT.@testset "2M microphysics - Seifert and Beheng 2006 effective radius and reflectivity" begin
        #setup
        ρₐ = FT(1)

        q_liq = [FT(2.128e-4), FT(2.128e-20), FT(1.6e-12), FT(0), FT(1.037e-25)]
        N_liq = [FT(15053529), FT(3), FT(5512), FT(0), FT(5.225e-12)]
        q_rai = [FT(1.573e-4), FT(1.573e-4), FT(1.9e-15), FT(0), FT(2.448e-27)]
        N_rai = [FT(510859), FT(510859), FT(0), FT(0), FT(5.136e-18)]

        # reference values
        rr = [FT(-12.561951), FT(-12.579899), FT(-150), FT(-150), FT(-150)]
        reff = [FT(2.319383e-5), FT(6.91594e-5), FT(0), FT(0), FT(0)]

        for (qₗ, Nₗ, qᵣ, Nᵣ, rₑ, Z) in zip(q_liq, N_liq, q_rai, N_rai, reff, rr)
            for SB in [SB2006, SB2006_no_limiters]

                #action
                Z_val = CMD.radar_reflectivity_2M(SB, qₗ, qᵣ, Nₗ, Nᵣ, ρₐ)
                rₑ_val = CMD.effective_radius_2M(SB, qₗ, qᵣ, Nₗ, Nᵣ, ρₐ)

                #test
                TT.@test rₑ_val ≈ rₑ atol = FT(1e-6)
                TT.@test Z_val ≈ Z atol = FT(1e-4)
            end
        end
    end

    TT.@testset "Effective radius - '1/3' power law from Liu and Hallett (1997)" begin
        #setup
        ρ_air = FT(1)
        ρ_w = FT(1000)
        q_liq = FT(2.128e-4)
        N_liq = FT(15053529)
        q_rai = FT(1.573e-4)
        N_rai = FT(510859)

        #action
        reff = CMD.effective_radius_Liu_Hallet_97(
            wtr,
            ρ_air,
            q_liq,
            N_liq,
            q_rai,
            N_rai,
        )
        #test
        TT.@test reff ≈ FT(2.66e-05) atol = FT(8e-6)

        TT.@test CMD.effective_radius_Liu_Hallet_97(
            wtr,
            ρ_air,
            q_liq,
            FT(100),
            FT(0),
            FT(0),
        ) == CMD.effective_radius_Liu_Hallet_97(wtr, ρ_air, q_liq)
    end

    TT.@testset "Constant effective radius" begin
        TT.@test CMD.effective_radius_const(cloud_liquid) == FT(14e-6)
        TT.@test CMD.effective_radius_const(cloud_ice) == FT(25e-6)
    end
end

println("Testing Float64")
test_cloud_diagnostics(Float64)

println("Testing Float32")
test_cloud_diagnostics(Float32)
