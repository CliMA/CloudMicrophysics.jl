import Test as TT

import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA

import CLIMAParameters as CP

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const param_set = EarthParameterSet()

# Atmospheric conditions
T = 294.0       # air temperature K
p = 100000.0    # air pressure Pa
w = 0.5         # vertical velocity m/s

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
# Sulfate - universal parameters
M_sulfate = 0.132
ρ_sulfate = 1770.0
ϕ_sulfate = 1.0
ν_sulfate = 3.0
ϵ_sulfate = 1.0

# Aerosol size distribution modes
accum_seasalt = AM.Mode(
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

coarse_seasalt = AM.Mode(
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

paper_mode_1 = AM.Mode(
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

# Aerosol size distributions
AM_1 = AM.AerosolDistribution((accum_seasalt,))
AM_2 = AM.AerosolDistribution((coarse_seasalt, accum_seasalt))
AM_3 = AM.AerosolDistribution((accum_seasalt, coarse_seasalt))

@TT.testset "callable and positive" begin

    @TT.test all(AA.mean_hygroscopicity(param_set, AM_3) .> 0.0)
    @TT.test AA.max_supersaturation(param_set, AM_3, T, p, w) > 0.0
    @TT.test all(AA.N_activated_per_mode(param_set, AM_3, T, p, w) .> 0.0)
    @TT.test all(AA.M_activated_per_mode(param_set, AM_3, T, p, w) .> 0.0)
    @TT.test AA.total_N_activated(param_set, AM_3, T, p, w) > 0.0
    @TT.test AA.total_M_activated(param_set, AM_3, T, p, w) > 0.0

end

@TT.testset "same mean hygroscopicity for the same aerosol" begin

    @TT.test AA.mean_hygroscopicity(param_set, AM_3)[1] ==
             AA.mean_hygroscopicity(param_set, AM_1)[1]
    @TT.test AA.mean_hygroscopicity(param_set, AM_3)[2] ==
             AA.mean_hygroscopicity(param_set, AM_1)[1]

end

@TT.testset "order of modes does not matter" begin

    @TT.test AA.total_N_activated(param_set, AM_3, T, p, w) ==
             AA.total_N_activated(param_set, AM_2, T, p, w)
    @TT.test AA.total_M_activated(param_set, AM_3, T, p, w) ==
             AA.total_M_activated(param_set, AM_2, T, p, w)

end

@TT.testset "Abdul-Razzak and Ghan 2000 Fig 1" begin

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

    N_act_frac = Vector{Float64}(undef, 6)

    it = 1
    for N_2_paper in N_2_obs
        paper_mode_2 = AM.Mode(
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

        AD = AM.AerosolDistribution((paper_mode_1, paper_mode_2))
        N_act_frac[it] =
            AA.N_activated_per_mode(param_set, AD, T, p, w)[1] / N_1_paper

        it += 1
    end

    @TT.test all(isapprox(N_act_frac, N_act_obs, rtol = 0.05))
end
