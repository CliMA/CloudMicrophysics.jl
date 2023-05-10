import Test as TT

import CloudMicrophysics as CM
import CLIMAParameters as CP
import Thermodynamics as TD

const AM = CM.AerosolModel
const AA = CM.AerosolActivation
const CMP = CM.Parameters

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

function test_aerosol_activation(FT)

    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    param_set = cloud_microphysics_parameters(toml_dict)
    thermo_params = CMP.thermodynamics_params(param_set)

    # Atmospheric conditions
    T = FT(294)    # air temperature K
    p = FT(1e5)    # air pressure Pa
    w = FT(0.5)    # vertical velocity m/s

    # We need the phase partition here only so that we can compute the
    # moist air R_m and cp_m in aerosol activation module.
    # We are assuming here saturated conditions and no liquid water or ice.
    # This is consistent with the assumptions of the aerosol activation scheme.
    p_vs = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    q_vs = 1 / (1 - CMP.molmass_ratio(param_set) * (p_vs - p) / p_vs)
    q = TD.PhasePartition(q_vs)

    # Accumulation mode
    r_dry_accum = FT(0.243 * 1e-6) # m
    stdev_accum = FT(1.4)          # -
    N_accum = FT(100 * 1e6)      # 1/m3

    # Coarse Mode
    r_dry_coarse = FT(1.5 * 1e-6)  # m
    stdev_coarse = FT(2.1)         # -
    N_coarse = FT(1e6)             # 1/m3

    # Abdul-Razzak and Ghan 2000 mode
    r_dry_paper = FT(0.05 * 1e-6)  # m
    stdev_paper = FT(2)            # -
    N_1_paper = FT(100 * 1e6)      # 1/m3

    # TODO - what κ values we should use?
    # Sea Salt - universal parameters
    M_seasalt = CMP.molmass_seasalt(param_set)
    ρ_seasalt = CMP.rho_seasalt(param_set)
    ϕ_seasalt = CMP.osm_coeff_seasalt(param_set)
    ν_seasalt = CMP.N_ion_seasalt(param_set)
    ϵ_seasalt = CMP.water_soluble_mass_frac_seasalt(param_set)
    κ_seasalt = CMP.kappa_seasalt(param_set)
    # Sulfate - universal parameters
    M_sulfate = CMP.molmass_sulfate(param_set)
    ρ_sulfate = CMP.rho_sulfate(param_set)
    ϕ_sulfate = CMP.osm_coeff_sulfate(param_set)
    ν_sulfate = CMP.N_ion_sulfate(param_set)
    ϵ_sulfate = CMP.water_soluble_mass_frac_sulfate(param_set)
    κ_sulfate = CMP.kappa_sulfate(param_set)

    # Aerosol size distribution modes
    accum_seasalt_B = AM.Mode_B(
        r_dry_accum,
        stdev_accum,
        N_accum,
        (FT(1.0),),
        (ϵ_seasalt,),
        (ϕ_seasalt,),
        (M_seasalt,),
        (ν_seasalt,),
        (ρ_seasalt,),
        1,
    )
    accum_seasalt_κ = AM.Mode_κ(
        r_dry_accum,
        stdev_accum,
        N_accum,
        (FT(1.0),),
        (FT(1.0),),
        (M_seasalt,),
        (κ_seasalt,),
        1,
    )

    coarse_seasalt_B = AM.Mode_B(
        r_dry_coarse,
        stdev_coarse,
        N_coarse,
        (FT(1.0),),
        (ϵ_seasalt,),
        (ϕ_seasalt,),
        (M_seasalt,),
        (ν_seasalt,),
        (ρ_seasalt,),
        1,
    )
    coarse_seasalt_κ = AM.Mode_κ(
        r_dry_coarse,
        stdev_coarse,
        N_coarse,
        (FT(1.0),),
        (FT(1.0),),
        (M_seasalt,),
        (κ_seasalt,),
        1,
    )

    paper_mode_1_B = AM.Mode_B(
        r_dry_paper,
        stdev_paper,
        N_1_paper,
        (FT(1.0),),
        (ϵ_sulfate,),
        (ϕ_sulfate,),
        (M_sulfate,),
        (ν_sulfate,),
        (ρ_sulfate,),
        1,
    )
    paper_mode_1_κ = AM.Mode_κ(
        r_dry_paper,
        stdev_paper,
        N_1_paper,
        (FT(1.0),),
        (FT(1.0),),
        (M_sulfate,),
        (κ_sulfate,),
        1,
    )

    # Aerosol size distributions
    AM_1_B = AM.AerosolDistribution((accum_seasalt_B,))
    AM_2_B = AM.AerosolDistribution((coarse_seasalt_B, accum_seasalt_B))
    AM_3_B = AM.AerosolDistribution((accum_seasalt_B, coarse_seasalt_B))

    AM_1_κ = AM.AerosolDistribution((accum_seasalt_κ,))
    AM_2_κ = AM.AerosolDistribution((coarse_seasalt_κ, accum_seasalt_κ))
    AM_3_κ = AM.AerosolDistribution((accum_seasalt_κ, coarse_seasalt_κ))

    TT.@testset "callable and positive" begin

        for AM_t in (AM_3_B, AM_3_κ)
            TT.@test all(
                AA.mean_hygroscopicity_parameter(param_set, AM_t) .> 0.0,
            )
            TT.@test AA.max_supersaturation(param_set, AM_t, T, p, w, q) > 0.0
            TT.@test all(
                AA.N_activated_per_mode(param_set, AM_t, T, p, w, q) .> 0.0,
            )
            TT.@test all(
                AA.M_activated_per_mode(param_set, AM_t, T, p, w, q) .> 0.0,
            )
            TT.@test AA.total_N_activated(param_set, AM_t, T, p, w, q) > 0.0
            TT.@test AA.total_M_activated(param_set, AM_t, T, p, w, q) > 0.0
        end
    end

    TT.@testset "same mean hygroscopicity for the same aerosol" begin

        TT.@test AA.mean_hygroscopicity_parameter(param_set, AM_3_B)[1] ==
                 AA.mean_hygroscopicity_parameter(param_set, AM_1_B)[1]

        TT.@test AA.mean_hygroscopicity_parameter(param_set, AM_3_B)[2] ==
                 AA.mean_hygroscopicity_parameter(param_set, AM_1_B)[1]

        TT.@test AA.mean_hygroscopicity_parameter(param_set, AM_3_κ)[1] ==
                 AA.mean_hygroscopicity_parameter(param_set, AM_1_κ)[1]

        TT.@test AA.mean_hygroscopicity_parameter(param_set, AM_3_κ)[2] ==
                 AA.mean_hygroscopicity_parameter(param_set, AM_1_κ)[1]

    end

    TT.@testset "B and kappa hygroscopicities are equivalent" begin

        TT.@test all(
            isapprox(
                AA.mean_hygroscopicity_parameter(param_set, AM_3_κ)[2],
                AA.mean_hygroscopicity_parameter(param_set, AM_3_B)[2],
                rtol = 0.1,
            ),
        )
    end

    TT.@testset "order of modes does not matter" begin

        TT.@test AA.total_N_activated(param_set, AM_3_B, T, p, w, q) ==
                 AA.total_N_activated(param_set, AM_2_B, T, p, w, q)
        TT.@test AA.total_M_activated(param_set, AM_3_B, T, p, w, q) ==
                 AA.total_M_activated(param_set, AM_2_B, T, p, w, q)

        TT.@test AA.total_N_activated(param_set, AM_3_κ, T, p, w, q) ==
                 AA.total_N_activated(param_set, AM_2_κ, T, p, w, q)
        TT.@test AA.total_M_activated(param_set, AM_3_κ, T, p, w, q) ==
                 AA.total_M_activated(param_set, AM_2_κ, T, p, w, q)
    end

    TT.@testset "Abdul-Razzak and Ghan 2000 Fig 1" begin

        # data read from Fig 1 in Abdul-Razzak and Ghan 2000
        # using https://automeris.io/WebPlotDigitizer/
        N_2_obs = [
            FT(18.74716810149539),
            FT(110.41572270049846),
            FT(416.00589034889026),
            FT(918.1014952424102),
            FT(1914.816492976891),
            FT(4919.913910285455),
        ]
        N_act_obs = [
            FT(0.7926937018577255),
            FT(0.7161078386950611),
            FT(0.5953670140462167),
            FT(0.4850589034888989),
            FT(0.34446080652469424),
            FT(0.162630267331219),
        ]

        N_act_frac_B = Vector{FT}(undef, 6)
        N_act_frac_κ = Vector{FT}(undef, 6)

        it = 1
        for N_2_paper in N_2_obs
            paper_mode_2_B = AM.Mode_B(
                r_dry_paper,
                stdev_paper,
                N_2_paper * FT(1e6),
                (FT(1.0),),
                (ϵ_sulfate,),
                (ϕ_sulfate,),
                (M_sulfate,),
                (ν_sulfate,),
                (ρ_sulfate,),
                1,
            )
            paper_mode_2_κ = AM.Mode_κ(
                r_dry_paper,
                stdev_paper,
                N_2_paper * FT(1e6),
                (FT(1.0),),
                (FT(1.0),),
                (M_sulfate,),
                (κ_sulfate,),
                1,
            )

            AD_B = AM.AerosolDistribution((paper_mode_1_B, paper_mode_2_B))
            AD_κ = AM.AerosolDistribution((paper_mode_1_κ, paper_mode_2_κ))

            N_act_frac_B[it] =
                AA.N_activated_per_mode(param_set, AD_B, T, p, w, q)[1] /
                N_1_paper
            N_act_frac_κ[it] =
                AA.N_activated_per_mode(param_set, AD_κ, T, p, w, q)[1] /
                N_1_paper
            it += 1
        end

        TT.@test all(isapprox(N_act_frac_B, N_act_obs, rtol = 0.05))
        TT.@test all(isapprox(N_act_frac_κ, N_act_obs, rtol = 0.1))

    end
end

println("")
println("Testing Float64")
test_aerosol_activation(Float64)

println("")
println("Testing Float32")
test_aerosol_activation(Float32)
