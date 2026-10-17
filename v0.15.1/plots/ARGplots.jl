import Plots

import CloudMicrophysics
import CLIMAParameters
import Thermodynamics

const PL = Plots
const AM = CloudMicrophysics.AerosolModel
const AA = CloudMicrophysics.AerosolActivation
const CP = CLIMAParameters
const CMP = CloudMicrophysics.Parameters
const TD = Thermodynamics

FT = Float64

tps = CMP.ThermodynamicsParameters(FT)
aip = CMP.AirProperties(FT)
ap = CMP.AerosolActivationParameters(FT)

include("plots/ARGdata.jl")

# Atmospheric conditions
T = 294.0         # air temperature
p = 1000.0 * 1e2   # air pressure

# We need the phase partition here only so that we can compute the
# moist R_m and cp_m in aerosol activation module.
# We are assuming here saturated conditions and no liquid water or ice.
# This is consistent with the assumptions of the aerosol activation scheme.
p_vs = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
q_vs = 1 / (1 - TD.Parameters.molmass_ratio(tps) * (p_vs - p) / p_vs)
q = TD.PhasePartition(q_vs, 0.0, 0.0)

# Sulfate
sulfate = CMP.Sulfate(FT)

# Insoluble material from ARG2000
ϵ_insol = 0.0         # solubility of insol
ϕ_insol = 0.0         # osmotic coeff of insol
M_insol = 0.044       # molar mass of insol
ν_insol = 0.0         # ions of salt in insol
ρ_insol = 1770.0      # density of insol
κ_insol = 0.0         # hygroscopicity of insol

function mass2vol(mass_mixing_ratios)
    if length(mass_mixing_ratios) == 2
        densities = (sulfate.ρ, ρ_insol)
    else
        densities = (sulfate.ρ,)
    end
    volfractions =
        (mass_mixing_ratios ./ densities) ./
        sum(mass_mixing_ratios ./ densities)
    return volfractions
end

# Abdul-Razzak and Ghan 2000
# https://doi.org/10.1029/1999JD901161
function make_ARG_figX(X)
    p1 = PL.plot()
    p2 = PL.plot()

    for v_B in (true, false)
        # mode 1 definitions
        r_dry_1 = 0.05 * 1e-6 # um
        stdev_1 = 2.0         # -
        N_1 = 100.0 * 1e6   # 1/m3

        if X in (1, 4)
            vol_mixing_ratios_1 = (1.0,)
            mass_mixing_ratios_1 = (1.0,)
            n_components_1 = 1
            if v_B
                paper_mode_1 = AM.Mode_B(
                    r_dry_1,
                    stdev_1,
                    N_1,
                    mass_mixing_ratios_1,
                    (sulfate.ϵ,),
                    (sulfate.ϕ,),
                    (sulfate.M,),
                    (sulfate.ν,),
                    (sulfate.ρ,),
                    n_components_1,
                )
            else
                paper_mode_1 = AM.Mode_κ(
                    r_dry_1,
                    stdev_1,
                    N_1,
                    vol_mixing_ratios_1,
                    mass_mixing_ratios_1,
                    (sulfate.M,),
                    (sulfate.κ,),
                    n_components_1,
                )
            end
        end

        if X in (2, 3, 5)
            vol_mixing_ratios_1 = (1.0, 0.0)
            mass_mixing_ratios_1 = (1.0, 0.0)
            n_components_1 = 2
            if v_B
                paper_mode_1 = AM.Mode_B(
                    r_dry_1,
                    stdev_1,
                    N_1,
                    mass_mixing_ratios_1,
                    (sulfate.ϵ, ϵ_insol),
                    (sulfate.ϕ, ϕ_insol),
                    (sulfate.M, M_insol),
                    (sulfate.ν, ν_insol),
                    (sulfate.ρ, ρ_insol),
                    n_components_1,
                )
            else
                paper_mode_1 = AM.Mode_κ(
                    r_dry_1,
                    stdev_1,
                    N_1,
                    vol_mixing_ratios_1,
                    mass_mixing_ratios_1,
                    (sulfate.M, M_insol),
                    (sulfate.κ, κ_insol),
                    n_components_1,
                )
            end
        end

        # parcel and mode 2 definitions
        len = 100
        global it = 1
        act_frac1 = Vector{Float64}(undef, len)
        act_frac2 = Vector{Float64}(undef, len)

        stdev_2 = 2.0   # -
        if X == 1
            w = 0.5                                         # vertical velocity, m/s
            r_dry_2 = 0.05 * 1e-6                           # um
            N_2 = range(100, stop = 5000, length = len) * 1e6   # 1/m3
            n_components_2 = 1                              # 1 mode
            mass_mixing_ratios_2 = (1.0,)                   # all sulfate
            vol_mixing_ratios_2 = (1.0,)                    # all sulfate

            for N2i in N_2
                if v_B
                    paper_mode_2 = AM.Mode_B(
                        r_dry_2,
                        stdev_2,
                        N2i,
                        mass_mixing_ratios_2,
                        (sulfate.ϵ,),
                        (sulfate.ϕ,),
                        (sulfate.M,),
                        (sulfate.ν,),
                        (sulfate.ρ,),
                        n_components_2,
                    )
                else
                    paper_mode_2 = AM.Mode_κ(
                        r_dry_2,
                        stdev_2,
                        N2i,
                        vol_mixing_ratios_2,
                        mass_mixing_ratios_2,
                        (sulfate.M,),
                        (sulfate.κ,),
                        n_components_2,
                    )
                end
                AD = AM.AerosolDistribution((paper_mode_1, paper_mode_2))
                act_frac1[it] =
                    AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q)[1] /
                    N_1
                act_frac2[it] =
                    AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q)[2] /
                    N2i
                global it += 1
            end

            x1_obs = Fig1_x_obs
            y1_obs = Fig1_y_obs
            x1_param = Fig1_x_param
            y1_param = Fig1_y_param

            x2_obs = Fig1_x_obs
            y2_obs = Fig1_y_obs
            x2_param = Fig1_x_param
            y2_param = Fig1_y_param

            xvar = N_2 * 1e-6
            xlabel = "Mode 2 aerosol number concentration [1/cm3]"
        elseif X == 2
            w = 0.5                                         # vertical velocity, m/s
            r_dry_2 = 0.05 * 1e-6                           # um
            N_2 = range(100, stop = 5000, length = len) * 1e6   # 1/m3
            n_components_2 = 2                              # 2 modes
            mass_mixing_ratios_2 = (0.1, 0.9)                # 10% sulfate, 90% insoluble
            vol_mixing_ratios_2 = mass2vol(mass_mixing_ratios_2)

            for N2i in N_2
                if v_B
                    paper_mode_2 = AM.Mode_B(
                        r_dry_2,
                        stdev_2,
                        N2i,
                        mass_mixing_ratios_2,
                        (sulfate.ϵ, ϵ_insol),
                        (sulfate.ϕ, ϕ_insol),
                        (sulfate.M, M_insol),
                        (sulfate.ν, ν_insol),
                        (sulfate.ρ, ρ_insol),
                        n_components_2,
                    )
                else
                    paper_mode_2 = AM.Mode_κ(
                        r_dry_2,
                        stdev_2,
                        N2i,
                        vol_mixing_ratios_2,
                        mass_mixing_ratios_2,
                        (sulfate.M, M_insol),
                        (sulfate.κ, κ_insol),
                        n_components_2,
                    )
                end
                AD = AM.AerosolDistribution((paper_mode_1, paper_mode_2))
                act_frac1[it] =
                    AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q)[1] /
                    N_1
                act_frac2[it] =
                    AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q)[2] /
                    N2i
                global it += 1
            end

            x1_obs = Fig2a_x_obs
            y1_obs = Fig2a_y_obs
            x1_param = Fig2a_x_param
            y1_param = Fig2a_y_param

            x2_obs = Fig2b_x_obs
            y2_obs = Fig2b_y_obs
            x2_param = Fig2b_x_param
            y2_param = Fig2b_y_param

            xvar = N_2 * 1e-6
            xlabel = "Mode 2 aerosol number concentration [1/cm3]"
        elseif X == 3
            w = 0.5                                     # vertical velocity, m/s
            r_dry_2 = 0.05 * 1e-6                       # um
            N_2 = 100 * 1e6                             # 1/m3
            n_components_2 = 2                          # 2 modes
            # ranging from 10% to 100% sulfate, 90% to 0% insoluble
            xvar = range(0.1, stop = 1, length = len)
            mass_mixing_ratios_2 = [(i, 1 - i) for i in xvar]

            for mmr2i in mass_mixing_ratios_2
                vmr2i = mass2vol(mmr2i)
                if v_B
                    paper_mode_2 = AM.Mode_B(
                        r_dry_2,
                        stdev_2,
                        N_2,
                        mmr2i,
                        (sulfate.ϵ, ϵ_insol),
                        (sulfate.ϕ, ϕ_insol),
                        (sulfate.M, M_insol),
                        (sulfate.ν, ν_insol),
                        (sulfate.ρ, ρ_insol),
                        n_components_2,
                    )
                else
                    paper_mode_2 = AM.Mode_κ(
                        r_dry_2,
                        stdev_2,
                        N_2,
                        vmr2i,
                        mmr2i,
                        (sulfate.M, M_insol),
                        (sulfate.κ, κ_insol),
                        n_components_2,
                    )
                end
                AD = AM.AerosolDistribution((paper_mode_1, paper_mode_2))
                act_frac1[it] =
                    AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q)[1] /
                    N_1
                act_frac2[it] =
                    AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q)[2] /
                    N_2
                global it += 1
            end

            x1_obs = Fig3a_x_obs
            y1_obs = Fig3a_y_obs
            x1_param = Fig3a_x_param
            y1_param = Fig3a_y_param

            x2_obs = Fig3b_x_obs
            y2_obs = Fig3b_y_obs
            x2_param = Fig3b_x_param
            y2_param = Fig3b_y_param

            xlabel = "Mode 2 soluble fraction (fraction of sulfate)"
        elseif X == 4
            w = 0.5                                     # vertical velocity, m/s
            r_dry_2 = range(0.01, stop = 0.5, length = len) * 1e-6 # um
            N_2 = 100 * 1e6                             # 1/m3
            n_components_2 = 1                          # 1 mode
            mass_mixing_ratios_2 = (1.0,)               # all sulfate
            vol_mixing_ratios_2 = mass2vol(mass_mixing_ratios_2)

            for rd2i in r_dry_2
                if v_B
                    paper_mode_2 = AM.Mode_B(
                        rd2i,
                        stdev_2,
                        N_2,
                        mass_mixing_ratios_2,
                        (sulfate.ϵ,),
                        (sulfate.ϕ,),
                        (sulfate.M,),
                        (sulfate.ν,),
                        (sulfate.ρ,),
                        n_components_2,
                    )
                else
                    paper_mode_2 = AM.Mode_κ(
                        rd2i,
                        stdev_2,
                        N_2,
                        vol_mixing_ratios_2,
                        mass_mixing_ratios_2,
                        (sulfate.M,),
                        (sulfate.κ,),
                        n_components_2,
                    )
                end
                AD = AM.AerosolDistribution((paper_mode_1, paper_mode_2))
                act_frac1[it] =
                    AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q)[1] /
                    N_1
                act_frac2[it] =
                    AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q)[2] /
                    N_2
                global it += 1
            end

            x1_obs = Fig4a_x_obs
            y1_obs = Fig4a_y_obs
            x1_param = Fig4a_x_param
            y1_param = Fig4a_y_param

            x2_obs = Fig4b_x_obs
            y2_obs = Fig4b_y_obs
            x2_param = Fig4b_x_param
            y2_param = Fig4b_y_param

            xvar = r_dry_2 * 1e6
            xlabel = "Mode 2 mean radius [um]"
        elseif X == 5
            w = range(0.01, stop = 5, length = len)     # vertical velocity, m/s
            r_dry_2 = 0.05 * 1e-6                       # um
            N_2 = 100 * 1e6                             # 1/m3
            n_components_2 = 2                          # 2 modes
            mass_mixing_ratios_2 = (0.1, 0.9)           # 10% sulfate, 90% insoluble
            vol_mixing_ratios_2 = mass2vol(mass_mixing_ratios_2)

            if v_B
                paper_mode_2 = AM.Mode_B(
                    r_dry_2,
                    stdev_2,
                    N_2,
                    mass_mixing_ratios_2,
                    (sulfate.ϵ, ϵ_insol),
                    (sulfate.ϕ, ϕ_insol),
                    (sulfate.M, M_insol),
                    (sulfate.ν, ν_insol),
                    (sulfate.ρ, ρ_insol),
                    n_components_2,
                )
            else
                paper_mode_2 = AM.Mode_κ(
                    r_dry_2,
                    stdev_2,
                    N_2,
                    vol_mixing_ratios_2,
                    mass_mixing_ratios_2,
                    (sulfate.M, M_insol),
                    (sulfate.κ, κ_insol),
                    n_components_2,
                )
            end
            AD = AM.AerosolDistribution((paper_mode_1, paper_mode_2))

            for wi in w
                act_frac1[it] =
                    AA.N_activated_per_mode(ap, AD, aip, tps, T, p, wi, q)[1] /
                    N_1
                act_frac2[it] =
                    AA.N_activated_per_mode(ap, AD, aip, tps, T, p, wi, q)[2] /
                    N_2
                global it += 1
            end

            x1_obs = Fig5a_x_obs
            y1_obs = Fig5a_y_obs
            x1_param = Fig5a_x_param
            y1_param = Fig5a_y_param

            x2_obs = Fig5b_x_obs
            y2_obs = Fig5b_y_obs
            x2_param = Fig5b_x_param
            y2_param = Fig5b_y_param

            xvar = w
            xlabel = "Vertical velocity, w [m/s]"
        end
        v_B ? label = "CliMA-B" : label = "CliMA-κ"
        PL.plot!(
            p1,
            xvar,
            act_frac1,
            label = label,
            ylim = [0, 1],
            ylabel = "Mode 1 act frac",
            title = "ARG2000 Fig " * string(X),
        )
        PL.plot!(
            p2,
            xvar,
            act_frac2,
            legend = false,
            ylim = [0, 1],
            xlabel = xlabel,
            ylabel = "Mode 2 act frac",
        )
        if v_B == false
            PL.scatter!(
                p1,
                x1_obs,
                y1_obs,
                markercolor = :black,
                label = "ARG2000 observations",
            )
            PL.plot!(
                p1,
                x1_param,
                y1_param,
                linecolor = :black,
                label = "ARG2000 parameterization",
            )
            PL.scatter!(p2, x2_obs, y2_obs, markercolor = :black)
            PL.plot!(p2, x2_param, y2_param, linecolor = :black)
        end
    end

    PL.plot(p1, p2, layout = (2, 1))
    PL.savefig("Abdul-Razzak_and_Ghan_fig_" * string(X) * ".svg")
end
