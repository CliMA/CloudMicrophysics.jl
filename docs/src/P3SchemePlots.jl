import CairoMakie as Plt
import CloudMicrophysics as CM
# import ..Parameters as CMP
const CMP = CM.Parameters
const P3 = CM.P3Scheme
const APS = CMP.AbstractCloudMicrophysicsParameters
FT = Float64

include(joinpath("..", "..", "test", "create_parameters.jl"))
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
const param_set = cloud_microphysics_parameters(toml_dict)
thermo_params = CMP.thermodynamics_params(param_set)

# bulk density of ice
const ρ_i::FT = CMP.ρ_cloud_ice(param_set)
const β_va::FT = 1.9
const α_va::FT = (7.38e-11) * 10^((6 * β_va) - 3)
const D_th::FT = ((FT(π) * ρ_i) / (6 * α_va))^(1 / (β_va - 3))

"""
mass(D, thresholds, F_r)

 - D: maximum particle dimension
 - thresholds: result of thresholds() function, i.e.
  a vector containing D_cr, D_gr, ρ_g, ρ_d, in that order;
  pass nothing in the case of no rime, and thresholds defaults to 0.0
 - F_r: rime mass fraction (q_rim/q_i)

 m(D) regime,
 which computes thresholds, classifies particles, and returns mass;
 used to create figures for the docs page.
"""
function mass(
    D::FT,
    F_r::FT,
    thresholds::Vector{FT} = [0.0, 0.0, 0.0, 0.0],
) where {FT <: Real}
    # Helper functions:
    """
    mass_s(D, ρ)

    - D: maximum particle dimension (m)
    - ρ: bulk ice density (note to self: use ρ_i for small ice and ρ_g for graupel) (kg m^-3)

    m(D) relation for spherical ice (small ice or completely rimed ice), returns mass (kg)
    """
    function mass_s(D::FT, ρ::FT) where {FT <: Real}
        return FT(π) / 6 * ρ * D^3
    end

    """
    mass_nl(D)

     - D: maximum particle dimension (m)


    m(D) relation for large, nonspherical ice (used for unrimed and dense types), returns mass (kg)
    """
    function mass_nl(D::FT) where {FT <: Real}
        return α_va * D^β_va
    end

    """
    mass_r(D, F_r)

    - D: maximum particle dimension (m)
    - F_r: rime mass fraction (q_rim/q_i)

    m(D) relation for partially rimed ice, returns mass (kg)
    """
    function mass_r(D::FT, F_r::FT) where {FT <: Real}
        return (α_va / (1 - F_r)) * D^β_va
    end

    # Mass regime:
    if D <= D_th
        return mass_s(D, ρ_i) # small spherical ice
    elseif F_r == 0
        return mass_nl(D) # large, nonspherical, unrimed ice
    else
        if D >= thresholds[1]
            return mass_r(D, F_r) # partially rimed ice
        elseif D < thresholds[1]
            if D >= thresholds[2]
                return mass_s(D, thresholds[3]) # graupel
            elseif D < thresholds[2]
                return mass_nl(D) # dense nonspherical ice
            end
        end
    end
end

"""
p3_m_plot1(colors, threshold_colors, len_D_range = 10000)

 - colors: string vector with three elements corresponding to the line colors
 - threshold_colors: string vector with three elements corresponding to the threshold line colors
 - len_D_range: amount of points to graph, defaults to 10000

Function that allows for the replication of Fig. 1a from Morrison and Milbrandt 2015.
"""
function p3_m_plot1(
    colors::Vector{String},
    threshold_colors::Vector{String},
    len_D_range::Int64 = 10000,
)
    D_range = range(1e-5, stop = 1e-2, length = len_D_range)

    fig1_a = Plt.Figure()

    ax1_a = Plt.Axis(
        fig1_a[1:7, 1:8],
        title = Plt.L"m(D) regime for $ρ_r = 400 kg m^{-3}$",
        xlabel = Plt.L"$D$ (mm)",
        ylabel = Plt.L"$m$ (kg)",
        xscale = Plt.log10,
        yscale = Plt.log10,
        xminorticksvisible = true,
        xminorticks = Plt.IntervalsBetween(5),
        yminorticksvisible = true,
        yminorticks = Plt.IntervalsBetween(3),
        xticks = [0.01, 0.1, 1, 10],
        aspect = 2,
    )

    fig1_a_0 = Plt.lines!(
        ax1_a,
        D_range * 1e3,
        [mass(D, 0.0) for D in D_range],
        color = colors[1],
    )

    fig1_a_5 = Plt.lines!(
        ax1_a,
        D_range * 1e3,
        [mass(D, 0.5, P3.thresholds(400.0, 0.5)) for D in D_range],
        color = colors[2],
    )

    fig1_a_8 = Plt.lines!(
        ax1_a,
        D_range * 1e3,
        [mass(D, 0.8, P3.thresholds(400.0, 0.8)) for D in D_range],
        color = colors[3],
    )

    d_cr_5 = Plt.lines!(
        ax1_a,
        [D = P3.thresholds(400.0, 0.5)[1] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "---",
        color = threshold_colors[2],
    )

    d_cr_8 = Plt.lines!(
        ax1_a,
        [D = P3.thresholds(400.0, 0.8)[1] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "---",
        color = threshold_colors[3],
    )

    d_gr_5 = Plt.lines!(
        ax1_a,
        [D = P3.thresholds(400.0, 0.5)[2] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "...",
        color = threshold_colors[2],
        linewidth = 1.5,
    )

    d_gr_8 = Plt.lines!(
        ax1_a,
        [D = P3.thresholds(400.0, 0.8)[2] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "...",
        color = threshold_colors[3],
        linewidth = 1.5,
    )

    d_tha = Plt.lines!(
        ax1_a,
        [D = D_th * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "---",
        color = "red",
    )

    leg1_a = Plt.Legend(
        fig1_a[8:9, 1],
        [fig1_a_0, fig1_a_5, fig1_a_8],
        [Plt.L"$F_{r} = 0.0$", Plt.L"$F_{r} = 0.5$", Plt.L"$F_{r} = 0.8$"],
    )

    leg1_a_dth = Plt.Legend(fig1_a[8:9, 3], [d_tha], [Plt.L"$D_{th}$"])

    leg1_a_dcr = Plt.Legend(
        fig1_a[8:9, 7],
        [d_cr_5, d_cr_8],
        [Plt.L"$D_{cr}$ for $F_{r} = 0.5$", Plt.L"$D_{cr}$ for $F_{r} = 0.8$"],
    )

    leg1_a_dgr = Plt.Legend(
        fig1_a[8:9, 5],
        [d_gr_5, d_gr_8],
        [Plt.L"$D_{gr}$ for $F_{r} = 0.5$", Plt.L"$D_{gr}$ for $F_{r} = 0.8$"],
    )

    Plt.save("MorrisonandMilbrandtFig1a.svg", fig1_a)
end

"""
p3_m_plot2(colors, threshold_colors, len_D_range = 10000)

 - colors: string vector with three elements corresponding to the line colors
 - threshold_colors: string vector with three elements corresponding to the threshold line colors
 - len_D_range: amount of points to graph, defaults to 10000

Function that allows for the replication of Fig. 1b from Morrison and Milbrandt 2015.
"""
function p3_m_plot2(
    colors::Vector{String},
    threshold_colors::Vector{String},
    len_D_range::Int64 = 10000,
)
    D_range = range(1e-5, stop = 1e-2, length = len_D_range)

    fig1_b = Plt.Figure()

    ax1_b = Plt.Axis(
        fig1_b[1:10, 1:11],
        title = Plt.L"m(D) regime for $F_r = 0.95$",
        xlabel = Plt.L"$D$ (mm)",
        ylabel = Plt.L"$m$ (kg)",
        xscale = Plt.log10,
        yscale = Plt.log10,
        yminorticksvisible = true,
        yminorticks = Plt.IntervalsBetween(3),
        xminorticksvisible = true,
        xminorticks = Plt.IntervalsBetween(5),
        xticks = [0.01, 0.1, 1, 10],
        aspect = 1.67,
    )

    fig1_b200 = Plt.lines!(
        ax1_b,
        D_range * 1e3,
        [mass(D, 0.95, P3.thresholds(200.0, 0.95)) for D in D_range],
        color = colors[1],
    )

    fig1_b400 = Plt.lines!(
        ax1_b,
        D_range * 1e3,
        [mass(D, 0.95, P3.thresholds(400.0, 0.95)) for D in D_range],
        color = colors[2],
    )

    fig1_b800 = Plt.lines!(
        ax1_b,
        D_range * 1e3,
        [mass(D, 0.95, P3.thresholds(800.0, 0.95)) for D in D_range],
        color = colors[3],
    )

    d_thb = Plt.lines!(
        ax1_b,
        [D = D_th * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "---",
        color = "red",
    )

    d_cr_200 = Plt.lines!(
        ax1_b,
        [D = P3.thresholds(200.0, 0.95)[1] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "---",
        color = threshold_colors[1],
    )

    d_cr_400 = Plt.lines!(
        ax1_b,
        [D = P3.thresholds(400.0, 0.95)[1] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "---",
        color = threshold_colors[2],
    )

    d_cr_800 = Plt.lines!(
        ax1_b,
        [D = P3.thresholds(800.0, 0.95)[1] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "---",
        color = threshold_colors[3],
    )

    d_gr_200 = Plt.lines!(
        ax1_b,
        [D = P3.thresholds(200.0, 0.95)[2] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "...",
        color = threshold_colors[1],
        linewidth = 1.5,
    )

    d_gr_400 = Plt.lines!(
        ax1_b,
        [D = P3.thresholds(400.0, 0.8)[2] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "...",
        color = threshold_colors[2],
    )

    d_gr_800 = Plt.lines!(
        ax1_b,
        [D = P3.thresholds(800.0, 0.95)[2] * 1e3 for D in D_range],
        range(1e-14, stop = 1e-4, length = len_D_range),
        linestyle = "...",
        color = threshold_colors[3],
        linewidth = 1.5,
    )

    leg1_b = Plt.Legend(
        fig1_b[11:12, 2],
        [fig1_b200, fig1_b400, fig1_b800],
        [
            Plt.L"$\rho_{r} = 200.0 kg m^{-3}$",
            Plt.L"$\rho_{r} = 400.0 kg m^{-3}$",
            Plt.L"$\rho_{r} = 800.0 kg m^{-3}$",
        ],
    )

    leg1_b_dth = Plt.Legend(fig1_b[11:12, 3], [d_thb], [Plt.L"$D_{th}$"])

    leg1_b_dcr = Plt.Legend(
        fig1_b[11:12, 6],
        [d_cr_200, d_cr_400, d_cr_800],
        [
            Plt.L"$D_{cr}$ for $\rho_{r} = 200.0 kg m^{-3}$",
            Plt.L"$D_{cr}$ for $\rho_{r} = 400.0 kg m^{-3}$",
            Plt.L"$D_{cr}$ for $\rho_{r} = 800.0 kg m^{-3}$",
        ],
    )

    leg1_b_dgr = Plt.Legend(
        fig1_b[11:12, 4],
        [d_gr_200, d_gr_400, d_gr_800],
        [
            Plt.L"$D_{gr}$ for $\rho_{r} = 200.0 kg m^{-3}$",
            Plt.L"$D_{gr}$ for $\rho_{r} = 400.0 kg m^{-3}$",
            Plt.L"$D_{gr}$ for $\rho_{r} = 800.0 kg m^{-3}$",
        ],
    )

    Plt.save("MorrisonandMilbrandtFig1b.svg", fig1_b)
end

"""
p3_heatmap(path, quantity, ρ_r_axis, F_r_axis)

 - path to look-up table file
 - quantity: option that toggles between heat maps of different quantities

Uses P3Scheme.read_threshold_table() to generate matrices from an existing
NetCDF lookup table and plots the matrices as heat map functions of ρ_r and F_r.
"""
function p3_heatmap(
    quantity::String,
    path::String = "/Users/rowan/Desktop/p3_scheme_work/CloudMicrophysics.jl/test.nc",
    ρ_r_axis::AbstractRange{FT} = range(start = 100, stop = 950, length = 10),
    F_r_axis::AbstractRange{FT} = range(start = 0.01, stop = 0.99, length = 10),
) where {FT <: Real}
    matrices = P3.read_threshold_table(path)
    if quantity === "D_cr"
        i = 1
    elseif quantity === "D_gr"
        i = 2
    elseif quantity === "ρ_g"
        i = 3
    elseif quantity === "ρ_d"
        i = 4
    else
        throw(ArgumentError("Specify a valid quantity to graph."))
    end

    fig, ax, hm = Plt.heatmap(ρ_r_axis, F_r_axis, matrices[i])
    Plt.Colorbar(fig[:, end + 1], hm)

    fig
end
