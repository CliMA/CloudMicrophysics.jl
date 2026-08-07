import Plots

import CloudMicrophysics
import CLIMAParameters
import Thermodynamics

const PL = Plots
const AM = CloudMicrophysics.AerosolModel
const AA = CloudMicrophysics.AerosolActivation
const CMP = CloudMicrophysics.Parameters
const TD = Thermodynamics

FT = Float64
tps = Thermodynamics.Parameters.ThermodynamicsParameters(FT)
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
p_vs = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
q_vs = 1 / (1 - TD.Parameters.molmass_ratio(tps) * (p_vs - p) / p_vs)
q = TD.PhasePartition(q_vs, 0.0, 0.0)

# Abdul-Razzak and Ghan 2000 Figure 1 mode 1
# https://doi.org/10.1029/1999JD901161
r_dry = 0.05 * 1e-6 # um
stdev = 2.0         # -
N_1 = 100.0 * 1e6   # 1/m3

# Sulfate - universal parameters
sulfate = CMP.Sulfate(FT)

n_components_1 = 1
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
    n_components_1,
)

N_2_range = range(0, stop = 5000 * 1e6, length = 100)
N_act_frac_B = Vector{Float64}(undef, 100)

it = 1
for N_2 in N_2_range
    n_components_2 = 1
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
        n_components_2,
    )
    AD_B = AM.AerosolDistribution((paper_mode_1_B, paper_mode_2_B))
    N_act_frac_B[it] =
        AA.N_activated_per_mode(ap, AD_B, aip, tps, T, p, w, q)[1] / N_1
    global it += 1
end

# data read from Fig 1 in Abdul-Razzak and Ghan 2000
# using https://automeris.io/WebPlotDigitizer/
include("plots/ARGdata.jl")

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
