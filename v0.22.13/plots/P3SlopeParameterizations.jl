import CairoMakie: Makie
import CloudMicrophysics as CM
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
import Printf

FT = Float64

logλs = @. log(10.0^FT(3:0.01:5))

function make_slope_plot(slope_law, title)
    fig = Makie.Figure(size = (400, 300), figure_padding = 20)

    ax = Makie.Axis(fig[1, 1]; title, xscale = log10, xlabel = "λ", ylabel = "μ")
    Makie.vlines!(λ_bnds, color = (:gray, 0.5))

    Makie.lines!(ax, exp.(logλs), P3.get_μ.(slope_law, logλs))

    return fig
end


# Power law parameterization
params = CMP.ParametersP3(FT)
slope_power_law = params.slope

log_λ_from_μ(spl::CMP.SlopePowerLaw, μ) = log((μ + spl.c) / spl.a) / spl.b
λ_bnds = log_λ_from_μ.(slope_power_law, [0.0, 6.0]) .|> exp

fig = Makie.with_theme(Makie.theme_minimal()) do
    make_slope_plot(slope_power_law, "μ as a function of λ (power law)")
end
Makie.save("P3SlopeParameterizations_power_law.svg", fig)

# Multiple solutions issue

function make_multiple_solutions_plot()
    state = P3.get_state(params; F_rim = 0.0, ρ_r = 400.0)
    L_known = 2.39e-4
    N_known = 1e5
    target_log_L_div_N = log(L_known) - log(N_known)
    shape_problem(log_λ) = P3.log_L_div_N(state, log_λ) - target_log_L_div_N
    shape_sol = shape_problem.(logλs)
    # Brute-force roots
    dists = P3.get_distribution_parameters_all_solutions(state; L = L_known, N = N_known)

    # Make figure
    fig = Makie.Figure(size = (800, 300), figure_padding = 20)
    ax = Makie.Axis(fig[1, 1];
        title = "Valid solutions of shape solver",
        xscale = log10, xlabel = "λ", ylabel = "shape solution", limits = ((10^3.55, 10^4.7), (-0.5, 0.5)),
    )
    Makie.vlines!(λ_bnds, color = (:gray, 0.5))
    Makie.hlines!(ax, 0; color = (:gray, 0.2))

    Makie.lines!(ax, exp.(logλs), shape_sol; color = :black)
    map(dists) do dist
        λ = exp(dist.log_λ)
        Makie.vlines!(ax, λ; label = Printf.@sprintf("λ = %.0f m⁻¹", λ))
    end
    Makie.axislegend(ax; framevisible = true)

    # N_ice axis
    mm_to_m = 1e-3 # m
    m_to_mm = 1e3  # mm
    m_to_cm = 1e2  # cm
    Ds = 10.0 .^ (-4:0.001:1.0) * mm_to_m  # in m
    Ds_mm = Ds * m_to_mm # for plotting
    ax_opts = (;
        limits = (extrema(Ds_mm), nothing),
        xscale = log10,
        xlabel = "D (mm)",
    )
    Makie.Axis(
        fig[1, 2];
        ax_opts...,
        # ylabel = L"$N_{ice}$ (1/cm³)", 
        ylabel = "[1/cm⁴]",
        title = "Number concentration distribution, N'(D)",
    )
    Makie.vlines!(state.D_th * m_to_mm, color = (:gray, 0.5))
    for dist in dists
        N′ = @. exp(P3.log_N′ice(dist, Ds))
        Makie.lines!(Ds_mm, N′ / m_to_cm^4)
        # Makie.lines!(Ds_mm, cumsum(N′.(Ds) .* Ds / m_to_cm^3))  # CUMULATIVE DISTRIBUTION
    end

    return fig
end

fig = Makie.with_theme(Makie.theme_minimal()) do
    make_multiple_solutions_plot()
end
Makie.save("P3SlopeParameterizations_multiple_solutions.svg", fig)
