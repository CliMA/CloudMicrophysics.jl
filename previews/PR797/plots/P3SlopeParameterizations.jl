import CairoMakie: Makie
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

FT = Float64

logλs = @. log(10.0^FT(3:0.01:5))

function make_slope_plot(hard_slope, smooth_slope, λ_bnds)
    fig = Makie.Figure(size = (500, 320), figure_padding = 20)

    ax = Makie.Axis(fig[1, 1]; title = "Slope parameter", xscale = log10,
        xlabel = "λ (m⁻¹)", ylabel = "μ",
    )
    Makie.vlines!(ax, λ_bnds; color = (:gray, 0.5))

    Makie.lines!(ax, exp.(logλs), P3.get_μ.(hard_slope, logλs); color = :black, label = "Hard limits")
    Makie.lines!(ax, exp.(logλs), P3.get_μ.(smooth_slope, logλs); color = :royalblue, label = "Smooth default")
    Makie.axislegend(ax; position = :lt)

    return fig
end


smooth_params = CMP.ParametersP3(FT)
hard_params = CMP.ParametersP3(FT; slope = CMP.SlopePowerLaw(FT))

log_λ_from_μ(spl::CMP.SlopePowerLaw, μ) = log((μ + spl.c) / spl.a) / spl.b
λ_bnds = log_λ_from_μ.(hard_params.slope, [0.0, 6.0]) .|> exp

fig = Makie.with_theme(Makie.theme_minimal()) do
    make_slope_plot(hard_params.slope, smooth_params.slope, λ_bnds)
end
Makie.save("P3SlopeParameterizations_slope_laws.svg", fig)

function make_multiple_solutions_plot(hard_params, smooth_params, λ_bnds)
    L_known = 2.39e-4
    N_known = 1e5
    hard_state = P3.P3State(hard_params, L_known, N_known, 0.0, 400.0)
    smooth_state = P3.P3State(smooth_params, L_known, N_known, 0.0, 400.0)
    target_logLdN = log(L_known) - log(N_known)
    hard_shape_problem(logλ) = P3.logLdivN(hard_state, logλ) - target_logLdN
    smooth_shape_problem(logλ) = P3.logLdivN(smooth_state, logλ) - target_logLdN
    hard_roots = P3.get_distribution_logλ_all_solutions(hard_state)
    smooth_roots = P3.get_distribution_logλ_all_solutions(smooth_state)
    @assert length(hard_roots) == 3
    @assert length(smooth_roots) == 1

    fig = Makie.Figure(size = (950, 350), figure_padding = 20)
    ax = Makie.Axis(fig[1, 1];
        title = "Shape solve for the same ice state",
        xscale = log10, xlabel = "λ (m⁻¹)", ylabel = "shape residual", limits = ((10^3.55, 10^4.7), (-0.5, 0.5)),
    )
    Makie.vlines!(ax, λ_bnds; color = (:gray, 0.5))
    Makie.hlines!(ax, 0; color = (:gray, 0.2))

    Makie.lines!(ax, exp.(logλs), hard_shape_problem.(logλs); color = :black, label = "Hard limits: three roots")
    Makie.lines!(ax, exp.(logλs), smooth_shape_problem.(logλs); color = :royalblue, label = "Smooth default: one root")
    Makie.scatter!(ax, exp.(hard_roots), zeros(length(hard_roots));
        color = :white, strokecolor = :black, strokewidth = 1.5, markersize = 12)
    Makie.scatter!(ax, exp.(smooth_roots), zeros(length(smooth_roots));
        color = :royalblue, marker = :diamond, markersize = 9)
    Makie.axislegend(ax; position = :rt)

    mm_to_m = 1e-3 # m
    m_to_mm = 1e3  # mm
    m_to_cm = 1e2  # cm
    Ds = 10.0 .^ (-4:0.001:1.0) * mm_to_m  # in m
    Ds_mm = Ds * m_to_mm # for plotting
    ax_psd = Makie.Axis(fig[1, 2];
        limits = (extrema(Ds_mm), nothing), xscale = log10,
        xlabel = "D (mm)", ylabel = "N′(D) (cm⁻⁴)",
        title = "Number distributions",
    )
    Makie.vlines!(ax_psd, P3.get_D_th(hard_params) * m_to_mm; color = (:gray, 0.5))
    for (i, logλ) in enumerate(hard_roots)
        N′ = P3.size_distribution(hard_state, logλ)
        Makie.lines!(ax_psd, Ds_mm, N′.(Ds) / m_to_cm^4;
            color = (:black, 0.5), label = i == 1 ? "Hard limits" : nothing)
    end
    N′ = P3.size_distribution(smooth_state, only(smooth_roots))
    Makie.lines!(ax_psd, Ds_mm, N′.(Ds) / m_to_cm^4; color = :royalblue, linewidth = 2, label = "Smooth default")
    Makie.axislegend(ax_psd; position = :rt)

    return fig
end

fig = Makie.with_theme(Makie.theme_minimal()) do
    make_multiple_solutions_plot(hard_params, smooth_params, λ_bnds)
end
Makie.save("P3SlopeParameterizations_multiple_solutions.svg", fig)
