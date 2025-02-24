using CairoMakie
import CloudMicrophysics as CM
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

axis_theme = Theme(
	Axis = (
		xscale = log10,
		xminorticksvisible = true, xminorticks = IntervalsBetween(5),
		xticks = [0.01, 0.1, 1, 10],
		limits = ((0.01, 10.0), nothing),
		xgridvisible = false, ygridvisible = false,
	),
	linewidth = 3,
	VLines = (
		linewidth = 1.5,
	)
)

logocolors = Makie.Colors.JULIA_LOGO_COLORS
cl = [logocolors.blue, logocolors.green, logocolors.red]

#! format: off
function p3_relations_plot()

	mm = 1e3 # m to mm, units

	FT = Float64
	F_liq = FT(0)
	params = CMP.ParametersP3(FT)
	
	D_range = 10.0 .^ (-2:0.01:2.0) * 1e-3
	fig = Figure(size=(1200, 900), figure_padding = 20)

	# define plot axis
	ax1L = Axis(fig[1, 1], yscale = log10, ylabel = L"$m(D)$ (kg)", title = L"Regimes for $ρ_r = 400$ kg/m³")
	ax1R = Axis(fig[1, 2], yscale = log10, title = L"Regimes for $F_{rim} = 0.95$")
	ax2L = Axis(fig[2, 1], yscale = log10, ylabel = L"$A(D)$ (m²)")
	ax2R = Axis(fig[2, 2], yscale = log10)
	ax3L = Axis(fig[3, 1], xlabel = L"$D$ (mm)", ylabel = L"$ρ$ (kg/m³)")
	ax3R = Axis(fig[3, 2], xlabel = L"$D$ (mm)")
	map(contents(fig[:, 2])) do ax ax.yticklabelsvisible = false end
	map(contents(fig[1:2, :])) do ax ax.xticklabelsvisible = false end

	colgap!(fig.layout, 40) # add space between columns so xticklabels don't overlap

	function lines_and_vlines!((ax1, ax2, ax3), state, color)
		lines!(ax1, D_range * mm, P3.ice_mass.(state, D_range); color)
		lines!(ax2, D_range * mm, P3.ice_area.(state, D_range); color)
		lines!(ax3, D_range * mm, P3.ice_density.(state, D_range); color)

		for ax in [ax1, ax2, ax3]
			vlines!(ax, state.D_th * mm; linestyle = (:dot, :loose), color = :gray)
			vlines!(ax, state.D_gr * mm; linestyle = :dash, color)
			vlines!(ax, state.D_cr * mm; linestyle = :dashdot, color)
		end
	end

	# Make a plot for each rime fraction
	ρ_r₀ = 400.0
	F_rims = [0.0, 0.5, 0.8]
	for (i, F_rim) in enumerate(F_rims)
		color = cl[i]
		state = P3.get_state(params; F_rim, F_liq, ρ_r = ρ_r₀)
		lines_and_vlines!((ax1L, ax2L, ax3L), state, color)
	end

	# Make a plot for each rime density
	F_rim₀ = 0.95
	ρ_rs = [200.0, 400.0, 800.0]
	for (i, ρ_r) in enumerate(ρ_rs)
		color = cl[i]
		state = P3.get_state(params; F_rim = F_rim₀, F_liq, ρ_r)
		lines_and_vlines!((ax1R, ax2R, ax3R), state, color)
	end

	linkyaxes!(ax1R, ax1L)
	ylims!(ax1R, 1e-14, 1e-4)
	linkyaxes!(ax2R, ax2L)
	linkyaxes!(ax3R, ax3L)

	# Add legend
	leg_elems = [
		map(color -> LineElement(;color), cl), 
		map(linestyle -> LineElement(;linestyle, linewidth = 1.5, color = :gray), [(:dot, :loose), (:dash, :dense), (:dashdot, :dense)])
	]
	thresh_labels = [L"$D_{th}$", L"$D_{gr}$", L"$D_{cr}$"]
	leg_kwargs = (; orientation = :horizontal, nbanks=10, tellheight = false, tellwidth = false, margin = (10, 10, 10, 10), halign = :left, valign = :top)
	Legend(fig[1,1], leg_elems, [Makie.latexstring.(L"F_{rim} = ", F_rims), thresh_labels], ["Rime fraction", "Thresholds"]; leg_kwargs...)
	Legend(fig[1,2], leg_elems, [Makie.latexstring.(L"ρ_r = ", ρ_rs), thresh_labels], ["Rime density", "Thresholds"]; leg_kwargs...)

	resize_to_layout!(fig)
	return fig
end
#! format: on

fig = with_theme(p3_relations_plot, axis_theme)
save("P3Scheme_relations.svg", fig)
fig