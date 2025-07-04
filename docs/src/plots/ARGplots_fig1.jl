import Plots as PL

import CloudMicrophysics as CM
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ClimaParams

FT = Float64
tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
aip = CMP.AirProperties(FT)
ap = CMP.AerosolActivationParameters(FT)

# Atmospheric conditions
T = 294.0         # air temperature
p = 1000.0 * 1e2   # air pressure
w = 0.5           # vertical velocity

# We need the phase partition here only so that we can compute the
# moist R_m and cp_m in aerosol activation module.
# We are assuming here saturated conditions and no liquid water or ice.
# This is consistent with the assumptions of the aerosol activation scheme.
p_vs = TDI.saturation_vapor_pressure_over_liquid(tps, T)
q_vs = 1 / (1 - 1 / TDI.Rd_over_Rv(tps) * (p_vs - p) / p_vs)

# Abdul-Razzak and Ghan 2000 Figure 1 mode 1
# https://doi.org/10.1029/1999JD901161
r_dry = 0.05 * 1e-6 # um
stdev = 2.0         # -
N_1 = 100.0 * 1e6   # 1/m3

# Sulfate - universal parameters
sulfate = CMP.Sulfate(FT)

mass_fractions_1 = (1.0,)
paper_mode_1_B = AM.Mode_B(
    r_dry,
    stdev,
    N_1,
    mass_fractions_1,
    (sulfate.ϵ,),
    (sulfate.ϕ,),
    (sulfate.M,),
    (sulfate.ν,),
    (sulfate.ρ,),
)

N_2_range = range(0, stop = 5000 * 1e6, length = 100)
N_act_frac_B = Vector{Float64}(undef, 100)

it = 1
for N_2 in N_2_range
    mass_fractions_2 = (1.0,)
    paper_mode_2_B = AM.Mode_B(
        r_dry,
        stdev,
        N_2,
        mass_fractions_2,
        (sulfate.ϵ,),
        (sulfate.ϕ,),
        (sulfate.M,),
        (sulfate.ν,),
        (sulfate.ρ,),
    )
    AD_B = AM.AerosolDistribution((paper_mode_1_B, paper_mode_2_B))
    N_act_frac_B[it] =
        AA.N_activated_per_mode(ap, AD_B, aip, tps, T, p, w, q_vs, FT(0), FT(0))[1] / N_1
    global it += 1
end

# data read from Fig 1 in Abdul-Razzak and Ghan 2000
# using https://automeris.io/WebPlotDigitizer/
include(joinpath(pkgdir(CM), "docs", "src", "plots", "ARGdata.jl"))

PL.plot(
    N_2_range * 1e-6,
    N_act_frac_B,
    label = "CliMA-B",
    xlabel = "Mode 2 aerosol number concentration [1/cm3]",
    ylabel = "Mode 1 number fraction activated",
)
PL.scatter!(
    Fig1_x_obs,
    Fig1_y_obs,
    markercolor = :black,
    label = "paper observations",
)
PL.plot!(
    Fig1_x_param,
    Fig1_y_param,
    linecolor = :black,
    label = "paper parameterization",
)

PL.savefig("Abdul-Razzak_and_Ghan_fig_1.svg")
