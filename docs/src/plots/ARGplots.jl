import CairoMakie as MK

import CloudMicrophysics as CM
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI

FT = Float64

tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
aip = CMP.AirProperties(FT)
ap = CMP.AerosolActivationParameters(FT)

include(joinpath(pkgdir(CM), "docs", "src", "plots", "ARGdata.jl"))

# Atmospheric conditions
T = 294.0         # air temperature
p = 1000.0 * 1e2   # air pressure

# We need the phase partition here only so that we can compute the
# moist R_m and cp_m in aerosol activation module.
# We are assuming here saturated conditions and no liquid water or ice.
# This is consistent with the assumptions of the aerosol activation scheme.
p_vs = TDI.saturation_vapor_pressure_over_liquid(tps, T)
q_vs = 1 / (1 - 1 / TDI.Rd_over_Rv(tps) * (p_vs - p) / p_vs)

# Sulfate
sulfate = CMP.Sulfate(FT)

# Insoluble material from ARG2000
ϵ_insol = 0.0         # solubility of insol
ϕ_insol = 0.0         # osmotic coeff of insol
M_insol = 0.044       # molar mass of insol
ν_insol = 0.0         # ions of salt in insol
ρ_insol = 1770.0      # density of insol
κ_insol = 0.0         # hygroscopicity of insol

all_components_B = (
    (sulfate.ϵ, ϵ_insol),
    (sulfate.ϕ, ϕ_insol),
    (sulfate.M, M_insol),
    (sulfate.ν, ν_insol),
    (sulfate.ρ, ρ_insol),
)

all_components_κ = (
    (sulfate.M, M_insol), 
    (sulfate.κ, κ_insol),
)

first_n_values(tup, n) = getindex(tup, 1:n)

function make_mode_B(r_dry, stdev, N, mass_mixing_ratios)
    n_components = length(mass_mixing_ratios)
    return AM.Mode_B(
        r_dry, stdev, N, mass_mixing_ratios,
        first_n_values.(all_components_B, n_components)...,
    )
end

function make_mode_κ(r_dry, stdev, N, vol_mixing_ratios, mass_mixing_ratios)
    n_components = length(mass_mixing_ratios)
    return AM.Mode_κ(
        r_dry, stdev, N, vol_mixing_ratios, mass_mixing_ratios,
        first_n_values.(all_components_κ, n_components)...,
    )
end

function mass2vol(mass_mixing_ratios)
    n_components = length(mass_mixing_ratios)
    densities = first_n_values((sulfate.ρ, ρ_insol), n_components)
    volfractions = @. mass_mixing_ratios / densities / $sum(mass_mixing_ratios / densities)
    return volfractions
end

function compute_activation_fractions(mode_1, mode_2, w)
    AD = AM.AerosolDistribution((mode_1, mode_2))
    act_per_mode = AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q_vs, FT(0), FT(0))
    return act_per_mode
end

# function get_fig_params(X, )
#     @assert X in (1, 2, 3, 4, 5)
#     xvar = Dict(
#         1 => :N_2,
#         2 => :N_2,
#         3 => :sulfate_fraction,
#         4 => :r_dry_2,
#         5 => :w,
#     )[X]
#     mass_mixing_ratios_1 = X in (1, 4) ? (1.0,) : (1.0, 0.0)
#     default_params = (;
#         mass_mixing_ratios_1,
#         r_dry_1 = 0.05 * 1e-6,  # um
#         stdev_1 = 2.0,          # -
#         N_1 = 100.0 * 1e6,      # 1/m3
#         stdev_2 = 2.0,
#         r_dry_2 = 0.05 * 1e-6,  # um
#         N_2 = 100.0 * 1e6,      # 1/m3
#         w = 0.5,                # vertical velocity, m/s
#     )

#     if X == 1
#         return (;
#             default_params...,
#             mass_mixing_ratios_2 = (1.0,),
#         )
#     elseif X == 2
#         return (;
#             default_params...,
#             mass_mixing_ratios_2 = (0.1, 0.9),
#         )
#     elseif X == 3
#         return (;
#             default_params...,
#             mass_mixing_ratios_2 = [(i, 1 - i) for i in range(0.1, stop = 1, length = 100)],
#         )
#     elseif X == 4
#     if X in (1, 4)
#         return (;
#             mass_mixing_ratios_1 = (1.0,),
#             w = 0.5,
#             r_dry_2 = 0.05 * 1e-6,
#             N_2 = range(100, stop = 5000, length = 100) * 1e6,
#             mass_mixing_ratios_2 = X == 1 ? (1.0,) : (0.1, 0.9),
#         )
#     elseif X in (2, 3, 5)
#         return (;
#             mass_mixing_ratios_1 = (1.0, 0.0),
#         )
#     end
# end

# Abdul-Razzak and Ghan 2000
# https://doi.org/10.1029/1999JD901161
function make_ARG_figX(X)
    # Create figure with subplots
    fig = MK.Figure(; size = (800, 600))
    ax1 = MK.Axis(fig[1, 1]; ylabel = "Mode 1", title = "ARG2000 Fig $X")
    ax2 = MK.Axis(fig[2, 1]; ylabel = "Mode 2", xlabel = "Mode 2 aerosol number concentration [1/cm3]")
    MK.Label(fig[:, 0], "Activation fraction"; rotation = π / 2)
    MK.linkaxes!(ax1, ax2)
    MK.limits!(ax2, (0, nothing), (0, 1))

    for v_B in (true, false)
        # mode 1 definitions
        r_dry_1 = 0.05 * 1e-6 # um
        stdev_1 = 2.0         # -
        N_1 = 100.0 * 1e6   # 1/m3

        if X in (1, 4)
            vol_mixing_ratios_1 = (1.0,)
            mass_mixing_ratios_1 = (1.0,)
            if v_B
                paper_mode_1 = AM.Mode_B(
                    r_dry_1, stdev_1, N_1, mass_mixing_ratios_1,
                    (sulfate.ϵ,), (sulfate.ϕ,), (sulfate.M,), (sulfate.ν,), (sulfate.ρ,),
                )
            else
                paper_mode_1 = AM.Mode_κ(
                    r_dry_1, stdev_1, N_1, vol_mixing_ratios_1, mass_mixing_ratios_1,
                    (sulfate.M,), (sulfate.κ,),
                )
            end
        end

        if X in (2, 3, 5)
            vol_mixing_ratios_1 = (1.0, 0.0)
            mass_mixing_ratios_1 = (1.0, 0.0)
            if v_B
                paper_mode_1 = make_mode_B(r_dry_1, stdev_1, N_1, mass_mixing_ratios_1)
            else
                paper_mode_1 = make_mode_κ(r_dry_1, stdev_1, N_1, vol_mixing_ratios_1, mass_mixing_ratios_1)
            end
        end

        # parcel and mode 2 definitions
        len = 100
        act_frac1 = Vector{Float64}(undef, len)
        act_frac2 = Vector{Float64}(undef, len)

        stdev_2 = 2.0   # -
        if X == 1
            w = 0.5                                         # vertical velocity, m/s
            r_dry_2 = 0.05 * 1e-6                           # um
            N_2 = range(100, stop = 5000, length = len) * 1e6   # 1/m3
            mass_mixing_ratios_2 = (1.0,)                   # all sulfate
            vol_mixing_ratios_2 = (1.0,)                    # all sulfate

            for (it, N2i) in enumerate(N_2)
                if v_B
                    paper_mode_2 = make_mode_B(r_dry_2, stdev_2, N2i, mass_mixing_ratios_2)
                else
                    paper_mode_2 = make_mode_κ(r_dry_2, stdev_2, N2i, vol_mixing_ratios_2, mass_mixing_ratios_2)
                end
                AD = AM.AerosolDistribution((paper_mode_1, paper_mode_2))
                act_per_mode = AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q_vs, FT(0), FT(0))
                act_frac1[it] = act_per_mode[1] / N_1
                act_frac2[it] = act_per_mode[2] / N2i
            end

            

            xvar = N_2 * 1e-6
            xlabel = "Mode 2 aerosol number concentration [1/cm3]"
        elseif X == 2
            w = 0.5                                         # vertical velocity, m/s
            r_dry_2 = 0.05 * 1e-6                           # um
            N_2 = range(100, stop = 5000, length = len) * 1e6   # 1/m3
            mass_mixing_ratios_2 = (0.1, 0.9)                # 10% sulfate, 90% insoluble
            vol_mixing_ratios_2 = mass2vol(mass_mixing_ratios_2)

            for (it, N2i) in enumerate(N_2)
                if v_B
                    paper_mode_2 = make_mode_B(r_dry_2, stdev_2, N2i, mass_mixing_ratios_2)
                else
                    paper_mode_2 = make_mode_κ(r_dry_2, stdev_2, N2i, vol_mixing_ratios_2, mass_mixing_ratios_2)
                end
                AD = AM.AerosolDistribution((paper_mode_1, paper_mode_2))
                act_per_mode = AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q_vs, FT(0), FT(0))
                act_frac1[it] = act_per_mode[1] / N_1
                act_frac2[it] = act_per_mode[2] / N2i
            end


            xvar = N_2 * 1e-6
            xlabel = "Mode 2 aerosol number concentration [1/cm3]"
        elseif X == 3
            w = 0.5                                     # vertical velocity, m/s
            r_dry_2 = 0.05 * 1e-6                       # um
            N_2 = 100 * 1e6                             # 1/m3
            # ranging from 10% to 100% sulfate, 90% to 0% insoluble
            xvar = range(0.1, stop = 1, length = len)
            mass_mixing_ratios_2 = [(i, 1 - i) for i in xvar]

            for (it, mmr2i) in enumerate(mass_mixing_ratios_2)
                vmr2i = mass2vol(mmr2i)
                if v_B
                    paper_mode_2 = make_mode_B(r_dry_2, stdev_2, N_2, mmr2i)
                else
                    paper_mode_2 = make_mode_κ(r_dry_2, stdev_2, N_2, vmr2i, mmr2i)
                end
                AD = AM.AerosolDistribution((paper_mode_1, paper_mode_2))
                act_per_mode = AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q_vs, FT(0), FT(0))
                act_frac1[it] = act_per_mode[1] / N_1
                act_frac2[it] = act_per_mode[2] / N_2
            end


            xlabel = "Mode 2 soluble fraction (fraction of sulfate)"
        elseif X == 4
            w = 0.5                                     # vertical velocity, m/s
            r_dry_2 = range(0.01, stop = 0.5, length = len) * 1e-6 # um
            N_2 = 100 * 1e6                             # 1/m3
            mass_mixing_ratios_2 = (1.0,)               # all sulfate
            vol_mixing_ratios_2 = mass2vol(mass_mixing_ratios_2)

            for (it, rd2i) in enumerate(r_dry_2)
                if v_B
                    paper_mode_2 = make_mode_B(rd2i, stdev_2, N_2, mass_mixing_ratios_2)
                else
                    paper_mode_2 = make_mode_κ(rd2i, stdev_2, N_2, vol_mixing_ratios_2, mass_mixing_ratios_2)
                end
                AD = AM.AerosolDistribution((paper_mode_1, paper_mode_2))
                act_per_mode = AA.N_activated_per_mode(ap, AD, aip, tps, T, p, w, q_vs, FT(0), FT(0))
                act_frac1[it] = act_per_mode[1] / N_1
                act_frac2[it] = act_per_mode[2] / N_2
            end


            xvar = r_dry_2 * 1e6
            xlabel = "Mode 2 mean radius [um]"
        elseif X == 5
            w = range(0.01, stop = 5, length = len)     # vertical velocity, m/s
            r_dry_2 = 0.05 * 1e-6                       # um
            N_2 = 100 * 1e6                             # 1/m3
            mass_mixing_ratios_2 = (0.1, 0.9)           # 10% sulfate, 90% insoluble
            vol_mixing_ratios_2 = mass2vol(mass_mixing_ratios_2)

            paper_mode_2_B = make_mode_B(r_dry_2, stdev_2, N_2, mass_mixing_ratios_2)
            paper_mode_2_κ = make_mode_κ(r_dry_2, stdev_2, N_2, vol_mixing_ratios_2, mass_mixing_ratios_2)
            AD = AM.AerosolDistribution((paper_mode_1, paper_mode_2))

            for (it, wi) in enumerate(w)
                act_per_mode = AA.N_activated_per_mode(ap, AD, aip, tps, T, p, wi, q_vs, FT(0), FT(0))
                act_frac1[it] = act_per_mode[1] / N_1
                act_frac2[it] = act_per_mode[2] / N_2
            end


            xvar = w
            xlabel = "Vertical velocity, w [m/s]"
        end
        
        # Update xlabel for the second subplot
        ax2.xlabel[] = xlabel
        
        v_B ? label = "CliMA-B" : label = "CliMA-κ"
        
        # Plot lines
        MK.lines!(ax1, xvar, act_frac1, label = label)
        MK.lines!(ax2, xvar, act_frac2)
    end

    x_obs, y_obs, x_param, y_param = get_ARG_data(X)
    MK.scatter!(ax1, x_obs, y_obs, color = :black, label = "ARG2000 observations")
    MK.lines!(ax1, x_param, y_param, color = :black, label = "ARG2000 parameterization")
    MK.scatter!(ax2, x_obs, y_obs, color = :black)
    MK.lines!(ax2, x_param, y_param, color = :black)

    # Add legend to the first subplot
    MK.axislegend(ax1, position = :rt)
    
    # Save the figure
    # MK.save("Abdul-Razzak_and_Ghan_fig_$X.svg", fig)
    
    return fig
end
