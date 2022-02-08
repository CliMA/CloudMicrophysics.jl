import Test

import CloudMicrophysics
import Thermodynamics

const TT = Test
const AM = CloudMicrophysics.AerosolModel
const AA = CloudMicrophysics.AerosolActivation
const TD = Thermodynamics

# build the parameter sets
param_set = CloudMicrophysicsParameters(
    src_parameter_dict,
    NoMicrophysicsParameters(),
    Thermodynamics.ThermodynamicsParameters(src_parameter_dict),
)

# we need the molmass ratio for testing - here just use the saved version
test_parameter_dict = ("molmass_ratio" => param_set.molmass_ratio)

# Atmospheric conditions
T = 294.0       # air temperature K
p = 100000.0    # air pressure Pa
w = 0.5         # vertical velocity m/s

# We need the phase partition here only so that we can compute the
# moist air R_m and cp_m in aerosol activation module.
# We are assuming here saturated conditions and no liquid water or ice.
# This is consistent with the assumptions of the aerosol activation scheme.
p_vs = TD.saturation_vapor_pressure(param_set, T, TD.Liquid())
q_vs = 1 / (1 - test_parameter_dict["molmass_ratio"] * (p_vs - p) / p_vs)
q = TD.PhasePartition(q_vs, 0.0, 0.0)

# Accumulation mode
r_dry_accum = 0.243 * 1e-6 # m
stdev_accum = 1.4          # -
N_accum = 100.0 * 1e6      # 1/m3

# Coarse Mode
r_dry_coarse = 1.5 * 1e-6  # m
stdev_coarse = 2.1         # -
N_coarse = 1.0 * 1e6       # 1/m3

# Abdul-Razzak and Ghan 2000 mode
r_dry_paper = 0.05 * 1e-6  # m
stdev_paper = 2.0          # -
N_1_paper = 100.0 * 1e6    # 1/m3

# TODO - move areosol properties to CLIMAParameters
# Sea Salt - universal parameters
M_seasalt = 0.058443
ρ_seasalt = 2170.0
ϕ_seasalt = 0.9
ν_seasalt = 2.0
ϵ_seasalt = 1.0
κ_seasalt = 1.12   # TODO - which value to take?
# Sulfate - universal parameters
M_sulfate = 0.132
ρ_sulfate = 1770.0
ϕ_sulfate = 1.0
ν_sulfate = 3.0
ϵ_sulfate = 1.0
κ_sulfate = 0.53   # TODO - which value to take?

# Aerosol size distribution modes
accum_seasalt_B = AM.Mode_B(
    r_dry_accum,
    stdev_accum,
    N_accum,
    (1.0,),
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
    (1.0,),
    (1.0,),
    (M_seasalt,),
    (κ_seasalt,),
    1,
)

coarse_seasalt_B = AM.Mode_B(
    r_dry_coarse,
    stdev_coarse,
    N_coarse,
    (1.0,),
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
    (1.0,),
    (1.0,),
    (M_seasalt,),
    (κ_seasalt,),
    1,
)

paper_mode_1_B = AM.Mode_B(
    r_dry_paper,
    stdev_paper,
    N_1_paper,
    (1.0,),
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
    (1.0,),
    (1.0,),
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
        TT.@test all(AA.mean_hygroscopicity_parameter(param_set, AM_t) .> 0.0)
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

    TT.@test all(isapprox(
        AA.mean_hygroscopicity_parameter(param_set, AM_3_κ)[2],
        AA.mean_hygroscopicity_parameter(param_set, AM_3_B)[2],
        rtol = 0.1,
    ))
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
        18.74716810149539,
        110.41572270049846,
        416.00589034889026,
        918.1014952424102,
        1914.816492976891,
        4919.913910285455,
    ]
    N_act_obs = [
        0.7926937018577255,
        0.7161078386950611,
        0.5953670140462167,
        0.4850589034888989,
        0.34446080652469424,
        0.162630267331219,
    ]

    N_act_frac_B = Vector{Float64}(undef, 6)
    N_act_frac_κ = Vector{Float64}(undef, 6)

    it = 1
    for N_2_paper in N_2_obs
        paper_mode_2_B = AM.Mode_B(
            r_dry_paper,
            stdev_paper,
            N_2_paper * 1e6,
            (1.0,),
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
            N_2_paper * 1e6,
            (1.0,),
            (1.0,),
            (M_sulfate,),
            (κ_sulfate,),
            1,
        )

        AD_B = AM.AerosolDistribution((paper_mode_1_B, paper_mode_2_B))
        AD_κ = AM.AerosolDistribution((paper_mode_1_κ, paper_mode_2_κ))

        N_act_frac_B[it] =
            AA.N_activated_per_mode(param_set, AD_B, T, p, w, q)[1] / N_1_paper
        N_act_frac_κ[it] =
            AA.N_activated_per_mode(param_set, AD_κ, T, p, w, q)[1] / N_1_paper
        it += 1
    end

    TT.@test all(isapprox(N_act_frac_B, N_act_obs, rtol = 0.05))
    TT.@test all(isapprox(N_act_frac_κ, N_act_obs, rtol = 0.1))

end
