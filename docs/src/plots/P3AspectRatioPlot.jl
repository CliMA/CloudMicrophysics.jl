using CairoMakie
import CloudMicrophysics as CM
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

axis_theme = Theme(
    Axis = (
        xscale = log10,
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(5),
        xticks = [0.01, 0.1, 1, 10],
        limits = ((0.01, 10.0), (0, 1.05)),
        xgridvisible = false,
        ygridvisible = false,
        xlabel = "D (mm)",
    ),
    linewidth = 3,
    VLines = (linewidth = 1.5,),
)

logocolors = Makie.Colors.JULIA_LOGO_COLORS
cl = [logocolors.blue, logocolors.green, logocolors.red]

#! format: off
function p3_relations_plot()

	mm = 1e3 # m to mm, units

	FT = Float64
	params = CMP.ParametersP3(FT)
	
	D_range = 10.0 .^ (-2:0.01:1.0) * 1e-3
	fig = Figure(size=(1200, 400), figure_padding = 20)

	# define plot axis
	ax1L = Axis(fig[1, 1], ylabel = L"$ϕ(D)$ [-]", title = L"Regimes for $ρ_{rim} = 400$ kg/m³")
	ax1R = Axis(fig[1, 2], title = L"Regimes for $F_{rim} = 0.4$")

	colgap!(fig.layout, 40) # add space between columns so xticklabels don't overlap

	function lines_and_vlines!(ax, state, color)
		lines!(ax, D_range * mm, P3.ϕᵢ.(state, D_range); color)

		(; D_th, D_gr, D_cr) = P3.get_thresholds_ρ_g(state)
		vlines!(ax, D_th * mm; linestyle = (:dot, :loose), color = :gray)
		vlines!(ax, D_gr * mm; linestyle = :dash, color)
		vlines!(ax, D_cr * mm; linestyle = :dashdot, color)
	end

	# Make a plot for each rime fraction
	ρ_rim₀ = 400.0
	F_rims = [0.0, 0.5, 0.8]
	for (i, F_rim) in enumerate(F_rims)
		color = cl[i]
		state = P3.get_state(params; F_rim, ρ_rim = ρ_rim₀)
		lines_and_vlines!(ax1L, state, color)
	end

	# Make a plot for each rime density
	F_rim₀ = 0.4
	ρ_rims = [200.0, 400.0, 800.0]
	for (i, ρ_rim) in enumerate(ρ_rims)
		color = cl[i]
		state = P3.get_state(params; F_rim = F_rim₀, ρ_rim)
		lines_and_vlines!(ax1R, state, color)
	end

	linkyaxes!(ax1R, ax1L)

	# Add legend
	leg_elems = [
		map(color -> LineElement(;color), cl), 
		map(linestyle -> LineElement(;linestyle, linewidth = 1.5, color = :gray), [(:dot, :loose), (:dash, :dense), (:dashdot, :dense)])
	]
	thresh_labels = [L"$D_{th}$", L"$D_{gr}$", L"$D_{cr}$"]
	leg_kwargs = (; orientation = :horizontal, nbanks=10, tellheight = false, tellwidth = false, margin = (10, 10, 10, 10), halign = :left, valign = :center)
	Legend(fig[1,1], leg_elems, [Makie.latexstring.(L"F_{rim} = ", F_rims), thresh_labels], ["Rime fraction", "Thresholds"]; leg_kwargs...)
	Legend(fig[1,2], leg_elems, [Makie.latexstring.(L"ρ_{rim} = ", ρ_rims), thresh_labels], ["Rime density", "Thresholds"]; leg_kwargs...)

	resize_to_layout!(fig)
	return fig
end
#! format: on

fig = with_theme(p3_relations_plot, axis_theme)
save("P3Scheme_aspect_ratio.svg", fig)
fig
