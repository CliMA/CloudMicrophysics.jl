import CairoMakie as Plt
import CloudMicrophysics as CM

const P3 = CM.P3Scheme
FT = Float64

# include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))
# toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
# const param_set = cloud_microphysics_parameters(toml_dict)

# const ρ_i::FT = param_set.density_ice_water
const ρ_i::FT = 916.7
const β_va::FT = 1.9
const α_va::FT = (7.38e-11) * 10^((6 * β_va) - 3)
const D_th::FT = ((FT(π) * ρ_i) / (6 * α_va))^(1 / (β_va - 3))

"""
mass(D, ρ_r, F_r)

 - D: maximum particle dimension
 - breakpoints: result of breakpoints() function, i.e.
  a vector containing D_cr, D_gr, ρ_g, ρ_d, in that order
 - F_r: rime mass fraction (q_rim/q_i)

 m(D) regime,
 which computes breakpoints, classifies particles, and returns mass;
 used to create figures for the docs page.
"""
function mass(D::FT, breakpoints::Vector{FT}, F_r::FT) where {FT <: Real}
    # Helper functions:
    """
    m_s(D, ρ)

    - D: maximum particle dimension
    - ρ: bulk ice density (note to self: use ρ_i for small ice and ρ_g for graupel)

    m(D) relation for spherical ice (small ice or completely rimed ice)
    """
    function m_s(D::FT, ρ::FT) where {FT <: Real}
        return (FT(π) / 6) * ρ * D^3
    end

    """
    m_nl(D)
    
     - D: maximum particle dimension
    
    
    m(D) relation for large, nonspherical ice (used for unrimed and dense types)
    """
    function m_nl(D::FT) where {FT <: Real}
        return α_va * D^β_va
    end

    """
    m_r(D, F_r)

    - D: maximum particle dimension
    - F_r: rime mass fraction (q_rim/q_i)

    m(D) relation for partially rimed ice
    """
    function m_r(D::FT, F_r::FT) where {FT <: Real}
        return (α_va / (1 - F_r)) * D^β_va
    end

    # Mass regime:
    if D <= D_th
        return m_s(D, ρ_i) # small spherical ice
    elseif F_r == 0
        return m_nl(D) # large, nonspherical, unrimed ice
    else
        if D >= breakpoints[1]
            return m_r(D, F_r) # partially rimed ice
        elseif D < breakpoints[1]
            if D >= breakpoints[2]
                return m_s(D, breakpoints[3]) # graupel
            elseif D < breakpoints[2]
                return m_nl(D) # dense nonspherical ice
            end
        end
    end
end


D_range = range(1e-5, stop=1e-2, length=10000)
ρ_r = [200.0, 400.0, 800.0]
F_r = [0.0, 0.5, 0.8, 0.95]

function p3_m_plot1(len_D_range::Int64 = 10000, colors::Vector{String})
    D_range = range(1e-5, stop = 1e-2, length = len_D_range)

    fig1_a = Plt.Figure()
        
    ax1_a = Plt.Axis(
        fig1_a[1:7,1:8],
        title = Plt.L"m(D) regime for $ρ_r = 400 kg m^{-3}$",
        xlabel = Plt.L"$D$ (mm)",
        ylabel = Plt.L"$m$ (kg)",
        xscale = Plt.log10, yscale = Plt.log10,
        xminorticksvisible = true,
        xminorticks =Plt.IntervalsBetween(5),
        yminorticksvisible = true,
        yminorticks = Plt.IntervalsBetween(3),
        xticks = [0.01, 0.1, 1, 10],
        aspect = 2
    )
    
    fig1_a_0 = Plt.lines!(
        ax1_a,
        D_range * 1e3,
        [mass(D, P3.breakpoints(400.0, 0.0), 0.0) for D in D_range],
        color = "indigo"
    )
    
    fig1_a_5 = Plt.lines!(
        ax1_a,
        D_range * 1e3,
        [mass(D, P3.breakpoints(400.0, 0.5), 0.5) for D in D_range],
        color = "chartreuse"
    )
    
    fig1_a_8 = Plt.lines!(
        ax1_a,
        D_range * 1e3,
        [mass(D, P3.breakpoints(400.0, 0.8), 0.8) for D in D_range],
        color = "red"
    )
    
    d_cr_5 = Plt.lines!(
        ax1_a,
        [D = P3.breakpoints(400.0, 0.5)[1] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = 10000),
        linestyle = "---",
        color = "chartreuse"
    )
    
    d_cr_8 = Plt.lines!(
        ax1_a,
        [D = P3.breakpoints(400.0, 0.8)[1] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = 10000),
        linestyle = "---",
        color = "red"
    )
    
    d_gr_5 = Plt.lines!(
        ax1_a,
        [D = P3.breakpoints(400.0, 0.5)[2] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = 10000),
        linestyle = "...",
        color = "chartreuse"
    )
    
    d_gr_8 = Plt.lines!(
        ax1_a,
        [D = P3.breakpoints(400.0, 0.8)[2] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = 10000),
        linestyle = "...",
        color = "red"
    )
    
    d_tha = Plt.lines!(
        ax1_a,
        [D = D_th * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = 10000),
        linestyle = "---"
    )
    
    leg1_a = Plt.Legend(
        fig1_a[8:9,1],
        [fig1_a_0, fig1_a_5, fig1_a_8],
        [Plt.L"$F_{r} = 0.0$", Plt.L"$F_{r} = 0.5$", Plt.L"$F_{r} = 0.8$"]
    )
    
    leg1_a_dth = Plt.Legend(
        fig1_a[8:9,3],
        [d_tha],
        [Plt.L"$D_{th}$"]
    )
    
    leg1_a_dcr = Plt.Legend(
        fig1_a[8:9,7], 
        [d_cr_5, d_cr_8],
        [Plt.L"$D_{cr}$ for $F_{r} = 0.5$", Plt.L"$D_{cr}$ for $F_{r} = 0.8$"]
    )
        
    leg1_a_dgr = Plt.Legend(
        fig1_a[8:9,5],
        [d_gr_5, d_gr_8],
        [Plt.L"$D_{gr}$ for $F_{r} = 0.5$", Plt.L"$D_{gr}$ for $F_{r} = 0.8$"]
    )
    
    Plt.save("MorrisonandMilbrandtFig1a.svg", fig1_a)
end

function p3_m_plot2(len_D_range::Int64 = 10000, colors::Vector{String})
    fig1_b = Plt.Figure()
    
    ax1_b = Plt.Axis(
        fig1_b[1:10,1:11],
        title = Plt.L"m(D) regime for $F_r = 0.95$",
        xlabel = Plt.L"$D$ (mm)",
        ylabel = Plt.L"$m$ (kg)",
        xscale = Plt.log10,
        yscale = Plt.log10,
        yminorticksvisible = true,
        yminorticks = Plt.IntervalsBetween(3),
        xminorticksvisible = true,
        xminorticks = Plt.IntervalsBetween(5),
        xticks = [0.01, 0.1, 1, 10], aspect = 1.67
    )
    
    fig1_b200 = Plt.lines!(
        ax1_b,
        D_range * 1e3,
        [P3.m(D, P3.breakpoints(200.0, 0.95), 0.95) for D in D_range],
        color = "indigo"
    )

    fig1_b400 = Plt.lines!(
        ax1_b,
        D_range * 1e3,
        [P3.m(D, P3.breakpoints(400.0, 0.95), 0.95) for D in D_range],
        color = "chartreuse"
    )
    
    fig1_b800 = Plt.lines!(ax1_b, D_range * 1e3,
    [P3.m(D, P3.breakpoints(800.0, 0.95), 0.95) for D in D_range], color = "red")
    d_thb = Plt.lines!(ax1_b, [D = D_th * 1e3 for D in D_range],
    range(1e-14, stop = 1e-4, length = 10000), linestyle = "---")
    d_cr_200 = Plt.lines!(ax1_b, [D = P3.breakpoints(200.0, 0.95)[1] * 1e3 for D in D_range],
    range(1e-14, stop = 1e-4, length = 10000), linestyle = "---", color = "indigo")
    d_cr_400 = Plt.lines!(ax1_b, [D = P3.breakpoints(400.0, 0.95)[1] * 1e3 for D in D_range],
    range(1e-14, stop = 1e-4, length = 10000), linestyle = "---", color = "chartreuse")
    d_cr_800 = Plt.lines!(ax1_b, [D = P3.breakpoints(800.0, 0.95)[1] * 1e3 for D in D_range],
    range(1e-14, stop = 1e-4, length = 10000), linestyle = "---", color = "red")
    d_gr_200 = Plt.lines!(ax1_b, [D = P3.breakpoints(200.0, 0.95)[2] * 1e3 for D in D_range],
    range(1e-14, stop = 1e-4, length = 10000), linestyle = "...", color = "indigo")
    d_gr_400 = Plt.lines!(ax1_b, [D = P3.breakpoints(400.0, 0.8)[2] * 1e3 for D in D_range],
    range(1e-14, stop = 1e-4, length = 10000), linestyle = "...", color = "chartreuse")
    d_gr_800 = Plt.lines!(ax1_b, [D = P3.breakpoints(800.0, 0.95)[2] * 1e3 for D in D_range],
    range(1e-14, stop = 1e-4, length = 10000), linestyle = "...", color = "red")
    leg1_b = Plt.Legend(fig1_b[11:12,2], [fig1_b200, fig1_b400, fig1_b800],
    [Plt.L"$\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$\rho_{r} = 400.0 kg m^{-3}$", Plt.L"$\rho_{r} = 800.0 kg m^{-3}$"])
    leg1_b_dth = Plt.Legend(fig1_b[11:12,3], [d_thb],
    [Plt.L"$D_{th}$"])
    leg1_b_dcr = Plt.Legend(fig1_b[11:12,6], [d_cr_200, d_cr_400, d_cr_800],
    [Plt.L"$D_{cr}$ for $\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$D_{cr}$ for $\rho_{r} = 400.0 kg m^{-3}$", 
    Plt.L"$D_{cr}$ for $\rho_{r} = 800.0 kg m^{-3}$"])
    leg1_b_dgr = Plt.Legend(fig1_b[11:12,4], [d_gr_200, d_gr_400, d_gr_800],
    [Plt.L"$D_{gr}$ for $\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$D_{gr}$ for $\rho_{r} = 400.0 kg m^{-3}$", 
    Plt.L"$D_{gr}$ for $\rho_{r} = 800.0 kg m^{-3}$"])
    Plt.save("P3_2.png", fig1_b, px_per_unit = 2)
end
















#     plt = Plt.Figure()
    
#     if fig == 1        
#         ax = Plt.Axis(
#             plt[1:7,1:8],
#             title = Plt.L"m(D) regime for $ρ_r = 400 kg m^{-3}$",
#             xlabel = Plt.L"$D$ (mm)",
#             ylabel = Plt.L"$m$ (kg)",
#             xscale = Plt.log10, yscale = Plt.log10,
#             xminorticksvisible = true,
#             xminorticks =Plt.IntervalsBetween(5),
#             yminorticksvisible = true,
#             yminorticks = Plt.IntervalsBetween(3),
#             xticks = [0.01, 0.1, 1, 10],
#             aspect = 2
#         )

#         ρ_r = [400.0, 400.0, 400.0]
#         F_r = [0.0, 0.5, 0.8]

#     elseif fig == 2
#         ax = Plt.Axis(
#             plt[1:10,1:11],
#             title = Plt.L"m(D) regime for $F_r = 0.95$",
#             xlabel = Plt.L"$D$ (mm)",
#             ylabel = Plt.L"$m$ (kg)",
#             xscale = Plt.log10,
#             yscale = Plt.log10,
#             xminorticksvisible = true,
#             xminorticks = Plt.IntervalsBetween(5),
#             xticks = [0.01, 0.1, 1, 10],
#             yminorticksvisible = true,
#             yminorticks = Plt.IntervalsBetween(3),
#             aspect = 1.67
#         )

#         ρ_r = [200.0, 400.0, 800.0]
#         F_r = [0.95, 0.95, 0.95]
    
#     end

#     D_range = range(1e-5, stop = 1e-2, length = len_D_range)
#     lines = Any[0, 0, 0]

#     for i in [1:3]
#         ρ_r_i = ρ_r[i]
#         F_r_i = F_r[i]
#         lines[i] = Plt.lines!(
#             ax,
#             D_range * 1e3,
#             [mass(D, P3.thresholds(ρ_r_i, F_r_i), F_r_i) for D in D_range],
#             color = colors[i]
#         )
#     end
# end

# function m_plot(
#         ax_range_1::UnitRange{Int64},
#         ax_range_2::UnitRange{Int64},
#         plt_title::AbstractString,
#         len_D_range::Int64,
#         ρ_r::Vector{Float64},
#         F_r::Vector{Float64},
#         colors::Vector{String},
#         fig::Int64
#     )
    
#     fig = Plt.Figure()
    
#     ax = Plt.Axis(fig[ax_range_1, ax_range_2],
#         title = plt_title,
#         xlabel = Plt.L"$D$ (mm)",
#         ylabel = Plt.L"$m$ (kg)",
#         xscale = Plt.log10, yscale = Plt.log10,
#         xminorticksvisible = true,
#         xminorticks =Plt.IntervalsBetween(5),
#         yminorticksvisible = true,
#         yminorticks = Plt.IntervalsBetween(3),
#         xticks = [0.01, 0.1, 1, 10],
#         aspect = 2
#     )
    
#     D_range = range(1e-5, stop = 1e-2, length = len_D_range)
#     lines = Any[0, 0, 0]
#     thresholds = (d_cr = Any[0, 0, 0], d_gr = Any[0, 0, 0])
    
#     d_th = Plt.lines!(
#         ax,
#         [D = D_th * 1e3 for D in D_range],
#         range(1e-14, stop = 1e-4, length = len_D_range),
#         linestyle = "---"
#     )

#     for i in range(1, length(colors))
#         lines[i] = Plt.lines!(
#             ax,
#             D_range * 1e3,
#             [P3.m(D, P3.breakpoints(ρ_r[i], F_r[i]), F_r[i]) for D in D_range],
#             color = colors[i]
#         )
#     end

#     for j in [1:2]
#         for i in range(1, length(colors))
#             if F_r[i] == 0.0
#                 thresholds[j] = deleteat!(thresholds[j], i)
#                 i -= 1
#             elseif i > length(thresholds[j])
#                 break
#             else
#                 thresholds[j][i] = Plt.lines!(
#                     ax,
#                     [D = P3.breakpoints(ρ_r[i], F_r[i])[j] * 1e3 for D in D_range],
#                     range(1e-14, stop = 1e-4, length = len_D_range),
#                     linestyle = "---",
#                     color = colors[i]
#                 )
#             end
#         end
#     end

#     if fig == 1
#         leg1_a = Plt.Legend(
#             fig[(lastindex(ax_range_1) + 1):(lastindex(ax_range_1) + 2), 1],
#             lines,
#             [Plt.L"$F_{r} = 0.0$", Plt.L"$F_{r} = 0.5$", Plt.L"$F_{r} = 0.8$"]
#         )
        
#         leg1_a_dth = Plt.Legend(
#             fig[(lastindex(ax_range_1) + 1):(lastindex(ax_range_1) + 2), 3],
#             [d_th],
#             [Plt.L"$D_{th}$"]
#         )
        
#         leg1_a_dcr = Plt.Legend(
#             fig[(lastindex(ax_range_1) + 1):(lastindex(ax_range_1) + 2), 7],
#             thresholds.d_cr,
#             [Plt.L"$D_{cr}$ for $F_{r} = 0.5$", Plt.L"$D_{cr}$ for $F_{r} = 0.8$"]
#         )
        
#         leg1_a_dgr = Plt.Legend(
#             fig[(lastindex(ax_range_1) + 1):(lastindex(ax_range_1) + 2), 5],
#             thresholds.d_gr,
#             [Plt.L"$D_{gr}$ for $F_{r} = 0.5$", Plt.L"$D_{gr}$ for $F_{r} = 0.8$"]
#         )
    
#     elseif fig == 2
#         leg1_b = Plt.Legend(
#             fig[(lastindex(ax_range_1) + 1):(lastindex(ax_range_1) + 2), 2],
#             lines,
#             [Plt.L"$\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$\rho_{r} = 400.0 kg m^{-3}$", Plt.L"$\rho_{r} = 800.0 kg m^{-3}$"]
#         )
        
#         leg1_b_dth = Plt.Legend(
#             fig1_b[(lastindex(ax_range_1) + 1):(lastindex(ax_range_1) + 2), 3],
#             d_th,
#             [Plt.L"$D_{th}$"]
#         )
        
#         leg1_b_dcr = Plt.Legend(
#             fig[(lastindex(ax_range_1) + 1):(lastindex(ax_range_1) + 2), 6],
#             thresholds.d_cr,
#             [Plt.L"$D_{cr}$ for $\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$D_{cr}$ for $\rho_{r} = 400.0 kg m^{-3}$", Plt.L"$D_{cr}$ for $\rho_{r} = 800.0 kg m^{-3}$"]
#         )
        
#         leg1_b_dgr = Plt.Legend(
#             fig[(lastindex(ax_range_1) + 1):(lastindex(ax_range_1) + 2), 4],
#             thresholds.d_gr,
#             [Plt.L"$D_{gr}$ for $\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$D_{gr}$ for $\rho_{r} = 400.0 kg m^{-3}$", Plt.L"$D_{gr}$ for $\rho_{r} = 800.0 kg m^{-3}$"]
#         )
    
#     end

#     Plt.save("P3_fig_" * string(fig) * ".svg", fig, px_per_unit = 2)

# end



# # Figure 1a from M&M with CairoMakie
# function p3fig(fig)
#     if fig == 1
#         fig1_a = Plt.Figure()
        
#         ax1_a = Plt.Axis(
#             fig1_a[1:7,1:8],
#             title = Plt.L"m(D) regime for $ρ_r = 400 kg m^{-3}$",
#             xlabel = Plt.L"$D$ (mm)",
#             ylabel = Plt.L"$m$ (kg)",
#             xscale = Plt.log10, yscale = Plt.log10,
#             xminorticksvisible = true,
#             xminorticks =Plt.IntervalsBetween(5),
#             yminorticksvisible = true,
#             yminorticks = Plt.IntervalsBetween(3),
#             xticks = [0.01, 0.1, 1, 10],
#             aspect = 2
#         )
        
#         fig1_a_0 = Plt.lines!(
#             ax1_a, D_range * 1e3,
#             [P3.m(D, P3.breakpoints(400.0, 0.0), 0.0) for D in D_range],
#             color = "indigo"
#             )
        
#         fig1_a_5 = Plt.lines!(
#             ax1_a,
#             D_range * 1e3,
#             [P3.m(D, P3.breakpoints(400.0, 0.5), 0.5) for D in D_range],
#             color = "chartreuse"
#         )
        
#         fig1_a_8 = Plt.lines!(
#             ax1_a,
#             D_range * 1e3,
#             [P3.m(D, P3.breakpoints(400.0, 0.8), 0.8) for D in D_range],
#             color = "red"
#         )
        
#         d_cr_5 = Plt.lines!(
#             ax1_a,
#             [D = P3.breakpoints(400.0, 0.5)[1] * 1e3 for D in D_range],
#             range(1e-14, stop = 1e-4, length = 10000),
#             linestyle = "---",
#             color = "chartreuse"
#         )
        
#         d_cr_8 = Plt.lines!(
#             ax1_a,
#             [D = P3.breakpoints(400.0, 0.8)[1] * 1e3 for D in D_range],
#             range(1e-14, stop = 1e-4, length = 10000),
#             linestyle = "---",
#             color = "red"
#         )
        
#         d_gr_5 = Plt.lines!(
#             ax1_a,
#             [D = P3.breakpoints(400.0, 0.5)[2] * 1e3 for D in D_range],
#             range(1e-14, stop = 1e-4, length = 10000),
#             linestyle = "...",
#             color = "chartreuse"
#         )
        
#         d_gr_8 = Plt.lines!(
#             ax1_a,
#             [D = P3.breakpoints(400.0, 0.8)[2] * 1e3 for D in D_range],
#             range(1e-14, stop = 1e-4, length = 10000),
#             linestyle = "...",
#             color = "red"
#         )
        
#         d_tha = Plt.lines!(
#             ax1_a,
#             [D = D_th * 1e3 for D in D_range],
#             range(1e-14, stop = 1e-4, length = 10000),
#             linestyle = "---"
#         )
        
#         leg1_a = Plt.Legend(
#             fig1_a[8:9,1],
#             [fig1_a_0, fig1_a_5, fig1_a_8],
#  [Plt.L"$F_{r} = 0.0$", Plt.L"$F_{r} = 0.5$", Plt.L"$F_{r} = 0.8$"])
# leg1_a_dth = Plt.Legend(fig1_a[8:9,3], [d_tha],
#  [Plt.L"$D_{th}$"])
# leg1_a_dcr = Plt.Legend(fig1_a[8:9,7], [d_cr_5, d_cr_8],
#  [Plt.L"$D_{cr}$ for $F_{r} = 0.5$", Plt.L"$D_{cr}$ for $F_{r} = 0.8$"])
# leg1_a_dgr = Plt.Legend(fig1_a[8:9,5], [d_gr_5, d_gr_8],
#  [Plt.L"$D_{gr}$ for $F_{r} = 0.5$", Plt.L"$D_{gr}$ for $F_{r} = 0.8$"])
# #Plt.save("P3_1.png", fig1_a, px_per_unit = 2)

# # Figure 1b from M&M with CairoMakie
# fig1_b = Plt.Figure()
# ax1_b = Plt.Axis(fig1_b[1:10,1:11],
#  title = Plt.L"m(D) regime for $F_r = 0.95$", xlabel = Plt.L"$D$ (mm)", ylabel = Plt.L"$m$ (kg)",
#  xscale = Plt.log10, yscale = Plt.log10, 
#  xminorticksvisible = true, xminorticks = Plt.IntervalsBetween(5),
#  xticks = [0.01, 0.1, 1, 10], aspect = 1.67)
# fig1_b200 = Plt.lines!(ax1_b, D_range * 1e3,
#  [P3.m(D, P3.breakpoints(200.0, 0.95), 0.95) for D in D_range], color = "indigo")
# fig1_b400 = Plt.lines!(ax1_b, D_range * 1e3,
#  [P3.m(D, P3.breakpoints(400.0, 0.95), 0.95) for D in D_range], color = "chartreuse")
# fig1_b800 = Plt.lines!(ax1_b, D_range * 1e3,
#  [P3.m(D, P3.breakpoints(800.0, 0.95), 0.95) for D in D_range], color = "red")
# d_thb = Plt.lines!(ax1_b, [D = D_th * 1e3 for D in D_range],
#  range(1e-14, stop = 1e-4, length = 10000), linestyle = "---")
# d_cr_200 = Plt.lines!(ax1_b, [D = P3.breakpoints(200.0, 0.95)[1] * 1e3 for D in D_range],
#  range(1e-14, stop = 1e-4, length = 10000), linestyle = "---", color = "indigo")
# d_cr_400 = Plt.lines!(ax1_b, [D = P3.breakpoints(400.0, 0.95)[1] * 1e3 for D in D_range],
#  range(1e-14, stop = 1e-4, length = 10000), linestyle = "---", color = "chartreuse")
# d_cr_800 = Plt.lines!(ax1_b, [D = P3.breakpoints(800.0, 0.95)[1] * 1e3 for D in D_range],
# range(1e-14, stop = 1e-4, length = 10000), linestyle = "---", color = "red")
# d_gr_200 = Plt.lines!(ax1_b, [D = P3.breakpoints(200.0, 0.95)[2] * 1e3 for D in D_range],
#  range(1e-14, stop = 1e-4, length = 10000), linestyle = "...", color = "indigo")
# d_gr_400 = Plt.lines!(ax1_b, [D = P3.breakpoints(400.0, 0.8)[2] * 1e3 for D in D_range],
# range(1e-14, stop = 1e-4, length = 10000), linestyle = "...", color = "chartreuse")
# d_gr_800 = Plt.lines!(ax1_b, [D = P3.breakpoints(800.0, 0.95)[2] * 1e3 for D in D_range],
#  range(1e-14, stop = 1e-4, length = 10000), linestyle = "...", color = "red")
# leg1_b = Plt.Legend(fig1_b[11:12,2], [fig1_b200, fig1_b400, fig1_b800],
#  [Plt.L"$\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$\rho_{r} = 400.0 kg m^{-3}$", Plt.L"$\rho_{r} = 800.0 kg m^{-3}$"])
# leg1_b_dth = Plt.Legend(fig1_b[11:12,3], [d_thb],
#  [Plt.L"$D_{th}$"])
# leg1_b_dcr = Plt.Legend(fig1_b[11:12,6], [d_cr_200, d_cr_400, d_cr_800],
#  [Plt.L"$D_{cr}$ for $\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$D_{cr}$ for $\rho_{r} = 400.0 kg m^{-3}$", 
#  Plt.L"$D_{cr}$ for $\rho_{r} = 800.0 kg m^{-3}$"])
# leg1_b_dgr = Plt.Legend(fig1_b[11:12,4], [d_gr_200, d_gr_400, d_gr_800],
#  [Plt.L"$D_{gr}$ for $\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$D_{gr}$ for $\rho_{r} = 400.0 kg m^{-3}$", 
#  Plt.L"$D_{gr}$ for $\rho_{r} = 800.0 kg m^{-3}$"])
# fig1_b

# #Plt.save("P3_2.png", fig1_b, px_per_unit = 2)
# #fig1_a