import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CairoMakie: Makie

const PSP3 = CMP.ParametersP3

FT = Float64

params = CMP.ParametersP3(FT; slope_law = :constant)

function get_values(
    params::CMP.ParametersP3,
    Chen2022::CMP.Chen2022VelType,
    L::FT,
    N::FT,
    ρ_a::FT,
    x_resolution::Int,
    y_resolution::Int,
) where {FT}
    F_rims = range(FT(0), stop = FT(1 - eps(FT)), length = x_resolution)
    ρ_rs = range(FT(25), stop = FT(975), length = y_resolution)

    V_m = zeros(x_resolution, y_resolution)
    V_m_ϕ = zeros(x_resolution, y_resolution)
    D_m = zeros(x_resolution, y_resolution)
    D_m_regimes = zeros(x_resolution, y_resolution)
    ϕᵢ = zeros(x_resolution, y_resolution)

    for i in 1:x_resolution
        for j in 1:y_resolution
            F_rim = F_rims[i]
            ρ_r = ρ_rs[j]
            state = P3.get_state(params; F_rim, ρ_r)
            dist = P3.get_distribution_parameters(state; L, N)

            V_m[i, j] = P3.ice_terminal_velocity(
                dist, Chen2022, ρ_a; use_aspect_ratio = false,
            )[2]
            V_m_ϕ[i, j] =
                P3.ice_terminal_velocity(dist, Chen2022, ρ_a; use_aspect_ratio = true)[2]

            D_m[i, j] = P3.D_m(dist)
            D_m_regimes[i, j] = D_m[i, j]
            ϕᵢ[i, j] = P3.ϕᵢ(state, D_m[i, j])

            # plot the regimes
            if state.D_th > D_m[i, j]
                # small spherical ice
                D_m_regimes[i, j] = 0
            elseif state.F_rim == 0
                # large nonspherical unrimed ice
                D_m_regimes[i, j] = 0
            elseif state.D_gr > D_m[i, j] >= state.D_th
                # dense nonspherical ice
                D_m_regimes[i, j] = 1
            elseif state.D_cr > D_m[i, j] >= state.D_gr
                # graupel
                D_m_regimes[i, j] = 2
            else #elseif D >= state.D_cr
                # partially rimed ice
                D_m_regimes[i, j] = 3
            end
        end
    end
    D_m *= 1e3
    return (; F_rims, ρ_rs, D_m_regimes, D_m, ϕᵢ, V_m, V_m_ϕ)
end

theme = Makie.Theme(
    Axis = (;
        width = 350,
        height = 350,
        limits = ((0, 1.0), (25, 975)),
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        yticks = [200, 400, 600, 800],
    ),
    Contour = (;
        color = :black,
        labels = true,
        levels = 3,
        linewidth = 1.5,
        labelsize = 18,
    ),
)

function figure_2()

    ### CALCULATE VALUES ###
    Chen2022 = CMP.Chen2022VelType(FT)
    # density of air in kg/m^3
    ρ_a = FT(1.2) #FT(1.293)
    # small D_m
    L_s = FT(0.0008)
    N_s = FT(1e6)
    # medium D_m
    L_m = FT(0.22)
    N_m = FT(1e6)
    # large D_m
    L_l = FT(0.7)
    N_l = FT(1e6)

    # get V_m and D_m
    xres = 100
    yres = 100

    (F_rims, ρ_rs, D_m_regimes_s, D_m_s, ϕᵢ_s, V_m_s, V_m_ϕ_s) =
        get_values(params, Chen2022, L_s, N_s, ρ_a, xres, yres)
    (F_rimm, ρ_rm, D_m_regimes_m, D_m_m, ϕᵢ_m, V_m_m, V_m_ϕ_m) =
        get_values(params, Chen2022, L_m, N_m, ρ_a, xres, yres)
    (F_riml, ρ_rl, D_m_regimes_l, D_m_l, ϕᵢ_l, V_m_l, V_m_ϕ_l) =
        get_values(params, Chen2022, L_l, N_l, ρ_a, xres, yres)

    ### PLOT ###
    fig = Makie.Figure()

    # Plot velocities as in Fig 2 in Morrison and Milbrandt 2015

    colormap = Makie.cgrad(:PuBuGn_3, 3, categorical = true)
    regime_contour_kwargs = (; levels = 3, colormap)

    row = 1
    ax1 = Makie.Axis(fig[row, 1]; title = "Particle regimes with small Dₘ")
    hm = Makie.contourf!(ax1, F_rims, ρ_rs, D_m_regimes_s; regime_contour_kwargs...)

    ax2 = Makie.Axis(fig[row, 2]; title = "Particle regimes with medium Dₘ")
    hm = Makie.contourf!(ax2, F_rimm, ρ_rm, D_m_regimes_m; regime_contour_kwargs...)

    ax3 = Makie.Axis(fig[row, 3]; title = "Particle regimes with large Dₘ")
    hm = Makie.contourf!(ax3, F_riml, ρ_rl, D_m_regimes_l; regime_contour_kwargs...)

    map(1:3) do col
        ticks = (
            [1, 2, 3],
            ["dense\n nonspherical ice", "graupel", "partially\n rimed ice"],
        )
        Makie.Colorbar(
            fig[row, col];
            colormap,
            ticks,
            vertical = false,
            width = Makie.Relative(0.95),
            height = 10,
            halign = 0.5,
            valign = 0.02,
            tellheight = false,
            colorrange = (0.5, 3.5),
            ticklabelpad = 0,
        )
    end

    function make_plots(
        row,
        col,
        F_rim,
        ρ_r;
        cfvals,
        cvals = nothing,
        title = "",
    )
        gp = fig[row, col]
        ax = Makie.Axis(gp; title)
        row3_opts =
            row == 3 ?
            (;
                ticklabelcolor = :white,
                ticks = 0:0.25:1,
                leftspinecolor = :white,
                rightspinecolor = :white,
                bottomspinecolor = :white,
                topspinecolor = :white,
            ) : (;)
        hm = Makie.contourf!(ax, F_rim, ρ_r, cfvals)
        Makie.Colorbar(
            gp,
            hm;
            halign = 0.05,
            valign = 0.05,
            height = Makie.Relative(0.60),
            width = 10,
            tellwidth = false,
            ticklabelpad = 0,
            row3_opts...,
        )
        !isnothing(cvals) && Makie.contour!(ax, F_rim, ρ_r, cvals)
    end

    row += 1
    title = "Dₘ (mm) (L = 8e-4 kgm⁻³, N = 1e6 m⁻³)"
    make_plots(row, 1, F_rims, ρ_rs; cfvals = D_m_s, title)
    title = "Dₘ (mm) (L = 0.22 kgm⁻³, N = 1e6 m⁻³)"
    make_plots(row, 2, F_rimm, ρ_rm; cfvals = D_m_m, title)
    title = "Dₘ (mm) (L = 0.7 kgm⁻³, N = 1e6 m⁻³)"
    make_plots(row, 3, F_riml, ρ_rl; cfvals = D_m_l, title)

    row += 1
    title = "ϕᵢ with small Dₘ"
    make_plots(row, 1, F_rims, ρ_rs; cfvals = ϕᵢ_s, cvals = D_m_s, title)
    title = "ϕᵢ with medium Dₘ"
    make_plots(row, 2, F_rimm, ρ_rm; cfvals = ϕᵢ_m, cvals = D_m_m, title)
    title = "ϕᵢ with large Dₘ"
    make_plots(row, 3, F_riml, ρ_rl; cfvals = ϕᵢ_l, cvals = D_m_l, title)

    row += 1
    title = "Vₘ (ϕᵢ = 1) with small Dₘ"
    make_plots(row, 1, F_rims, ρ_rs; cfvals = V_m_s, cvals = D_m_s, title)
    title = "Vₘ (ϕᵢ = 1) with medium Dₘ",
    make_plots(row, 2, F_rimm, ρ_rm; cfvals = V_m_m, cvals = D_m_m, title)
    title = "Vₘ (ϕᵢ = 1) with large Dₘ"
    make_plots(row, 3, F_riml, ρ_rl; cfvals = V_m_l, cvals = D_m_l, title)

    row += 1
    title = "Vₘ (using ϕᵢ) with small Dₘ"
    make_plots(row, 1, F_rims, ρ_rs; cfvals = V_m_ϕ_s, cvals = D_m_s, title)
    title = "Vₘ (using ϕᵢ) with medium Dₘ"
    make_plots(row, 2, F_rimm, ρ_rm; cfvals = V_m_ϕ_m, cvals = D_m_m, title)
    title = "Vₘ (using ϕᵢ) with large Dₘ"
    make_plots(row, 3, F_riml, ρ_rl; cfvals = V_m_ϕ_l, cvals = D_m_l, title)

    axs = filter(ax -> ax isa Makie.Axis, fig.content)
    Makie.linkaxes!(axs...)

    Makie.resize_to_layout!(fig)
    Makie.save("MorrisonandMilbrandtFig2.svg", fig)
    fig
end
#! format: on

# Terminal Velocity figure
Makie.with_theme(figure_2, theme)
