import Test as TT

import ClimaParams as CP

import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.Parameters as CMP

function test_aerosol_activation(FT)

    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    aip = CMP.AirProperties(FT)
    #default ARG2000 parameter values
    ap_default = CMP.AerosolActivationParameters(FT)
    # calibrated values based on PySDM
    override_file =
        joinpath(pkgdir(CM), "src", "parameters", "toml", "ARG2000.toml")
    toml_dict = CP.create_toml_dict(FT; override_file)
    ap_calibrated = CMP.AerosolActivationParameters(toml_dict)

    # Atmospheric conditions
    T = FT(294)    # air temperature K
    p = FT(1e5)    # air pressure Pa
    w = FT(0.5)    # vertical velocity m/s

    # We need the phase partition here only so that we can compute the
    # moist air R_m and cp_m in aerosol activation module.
    # We are assuming here saturated conditions and no liquid water or ice.
    # This is consistent with the assumptions of the aerosol activation scheme.
    p_vs = TDI.saturation_vapor_pressure_over_liquid(tps, T)
    q_vs = 1 / (1 - 1 / TDI.Rd_over_Rv(tps) * (p_vs - p) / p_vs)
    q_tot = q_vs
    q_liq = FT(0)
    q_ice = FT(0)
    N_liq = FT(1000)
    N_ice = FT(1000)

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
    seasalt = CMP.Seasalt(FT)
    # Sulfate - universal parameters
    sulfate = CMP.Sulfate(FT)

    # Aerosol size distribution modes
    accum_seasalt_B = AM.Mode_B(
        r_dry_accum,
        stdev_accum,
        N_accum,
        (FT(1.0),),
        (seasalt.ϵ,),
        (seasalt.ϕ,),
        (seasalt.M,),
        (seasalt.ν,),
        (seasalt.ρ,),
    )
    accum_seasalt_κ = AM.Mode_κ(
        r_dry_accum,
        stdev_accum,
        N_accum,
        (FT(1.0),),
        (FT(1.0),),
        (seasalt.M,),
        (seasalt.κ,),
    )

    coarse_seasalt_B = AM.Mode_B(
        r_dry_coarse,
        stdev_coarse,
        N_coarse,
        (FT(1.0),),
        (seasalt.ϵ,),
        (seasalt.ϕ,),
        (seasalt.M,),
        (seasalt.ν,),
        (seasalt.ρ,),
    )
    coarse_seasalt_κ = AM.Mode_κ(
        r_dry_coarse,
        stdev_coarse,
        N_coarse,
        (FT(1.0),),
        (FT(1.0),),
        (seasalt.M,),
        (seasalt.κ,),
    )

    paper_mode_1_B = AM.Mode_B(
        r_dry_paper,
        stdev_paper,
        N_1_paper,
        (FT(1.0),),
        (sulfate.ϵ,),
        (sulfate.ϕ,),
        (sulfate.M,),
        (sulfate.ν,),
        (sulfate.ρ,),
    )
    paper_mode_1_κ = AM.Mode_κ(
        r_dry_paper,
        stdev_paper,
        N_1_paper,
        (FT(1.0),),
        (FT(1.0),),
        (sulfate.M,),
        (sulfate.κ,),
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
            for ap in (ap_default, ap_calibrated)
                TT.@test all(AA.mean_hygroscopicity_parameter(ap, AM_t) .> 0.0)

                TT.@test AA.max_supersaturation(
                    ap,
                    AM_t,
                    aip,
                    tps,
                    T,
                    p,
                    w,
                    q_tot, q_liq, q_ice,
                ) > 0.0
                TT.@test AA.max_supersaturation(
                    ap,
                    AM_t,
                    aip,
                    tps,
                    T,
                    p,
                    w,
                    q_tot, q_liq, q_ice,
                    N_liq, N_ice,
                ) ≥ 0.0

                TT.@test all(
                    AA.N_activated_per_mode(ap, AM_t, aip, tps, T, p, w, q_tot, q_liq, q_ice) .>
                    0.0,
                )
                TT.@test all(
                    AA.N_activated_per_mode(ap, AM_t, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice) .≥
                    0.0,
                )

                TT.@test all(
                    AA.M_activated_per_mode(ap, AM_t, aip, tps, T, p, w, q_tot, q_liq, q_ice) .>
                    0.0,
                )
                TT.@test all(
                    AA.M_activated_per_mode(ap, AM_t, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice) .≥
                    0.0,
                )

                TT.@test AA.total_N_activated(ap, AM_t, aip, tps, T, p, w, q_tot, q_liq, q_ice) >
                         0.0
                TT.@test AA.total_N_activated(ap, AM_t, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice) ≥
                         0.0

                TT.@test AA.total_M_activated(ap, AM_t, aip, tps, T, p, w, q_tot, q_liq, q_ice) >
                         0.0
                TT.@test AA.total_M_activated(ap, AM_t, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice) ≥
                         0.0
            end
        end
    end

    TT.@testset "same mean hygroscopicity for the same aerosol" begin
        for ap in (ap_default, ap_calibrated)

            TT.@test AA.mean_hygroscopicity_parameter(ap, AM_3_B)[1] ==
                     AA.mean_hygroscopicity_parameter(ap, AM_1_B)[1]

            TT.@test AA.mean_hygroscopicity_parameter(ap, AM_3_B)[2] ==
                     AA.mean_hygroscopicity_parameter(ap, AM_1_B)[1]

            TT.@test AA.mean_hygroscopicity_parameter(ap, AM_3_κ)[1] ==
                     AA.mean_hygroscopicity_parameter(ap, AM_1_κ)[1]

            TT.@test AA.mean_hygroscopicity_parameter(ap, AM_3_κ)[2] ==
                     AA.mean_hygroscopicity_parameter(ap, AM_1_κ)[1]
        end
    end

    TT.@testset "B and kappa hygroscopicities are equivalent" begin
        for ap in (ap_default, ap_calibrated)
            TT.@test all(
                isapprox(
                    AA.mean_hygroscopicity_parameter(ap, AM_3_κ)[2],
                    AA.mean_hygroscopicity_parameter(ap, AM_3_B)[2],
                    rtol = 0.1,
                ),
            )
        end
    end

    TT.@testset "order of modes does not matter" begin
        for ap in (ap_default, ap_calibrated)
            TT.@test AA.total_N_activated(ap, AM_3_B, aip, tps, T, p, w, q_tot, q_liq, q_ice) ==
                     AA.total_N_activated(ap, AM_2_B, aip, tps, T, p, w, q_tot, q_liq, q_ice)
            TT.@test AA.total_M_activated(ap, AM_3_B, aip, tps, T, p, w, q_tot, q_liq, q_ice) ==
                     AA.total_M_activated(ap, AM_2_B, aip, tps, T, p, w, q_tot, q_liq, q_ice)

            TT.@test AA.total_N_activated(ap, AM_3_κ, aip, tps, T, p, w, q_tot, q_liq, q_ice) ==
                     AA.total_N_activated(ap, AM_2_κ, aip, tps, T, p, w, q_tot, q_liq, q_ice)
            TT.@test AA.total_M_activated(ap, AM_3_κ, aip, tps, T, p, w, q_tot, q_liq, q_ice) ==
                     AA.total_M_activated(ap, AM_2_κ, aip, tps, T, p, w, q_tot, q_liq, q_ice)
        end
    end

    TT.@testset "Abdul-Razzak and Ghan 2000 Fig 1" begin
        for ap in (ap_default,)

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
                    (sulfate.ϵ,),
                    (sulfate.ϕ,),
                    (sulfate.M,),
                    (sulfate.ν,),
                    (sulfate.ρ,),
                )
                paper_mode_2_κ = AM.Mode_κ(
                    r_dry_paper,
                    stdev_paper,
                    N_2_paper * FT(1e6),
                    (FT(1.0),),
                    (FT(1.0),),
                    (sulfate.M,),
                    (sulfate.κ,),
                )

                AD_B = AM.AerosolDistribution((paper_mode_1_B, paper_mode_2_B))
                AD_κ = AM.AerosolDistribution((paper_mode_1_κ, paper_mode_2_κ))

                N_act_frac_B[it] =
                    AA.N_activated_per_mode(ap, AD_B, aip, tps, T, p, w, q_tot, q_liq, q_ice)[1] /
                    N_1_paper
                N_act_frac_κ[it] =
                    AA.N_activated_per_mode(ap, AD_κ, aip, tps, T, p, w, q_tot, q_liq, q_ice)[1] /
                    N_1_paper
                it += 1
            end

            TT.@test all(isapprox(N_act_frac_B, N_act_obs, rtol = 0.05))
            TT.@test all(isapprox(N_act_frac_κ, N_act_obs, rtol = 0.1))
        end
    end
end

TT.@testset "Aerosol Activation Tests ($FT)" for FT in (Float64, Float32)
    test_aerosol_activation(FT)
end
nothing
