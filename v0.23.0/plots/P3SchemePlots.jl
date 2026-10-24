import CairoMakie: Makie
import CloudMicrophysics as CM
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

axis_theme = Makie.Theme(
    Axis = (
        xscale = log10,
        xminorticksvisible = true,
        xminorticks = Makie.IntervalsBetween(5),
        xticks = [0.01, 0.1, 1, 10],
        limits = ((0.01, 10.0), nothing),
        xgridvisible = false,
        ygridvisible = false,
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
	
	D_range = 10.0 .^ (-2:0.01:2.0) * 1e-3
	fig = Makie.Figure(size=(1200, 900), figure_padding = 20)

	# define plot axis
	ax1L = Makie.Axis(fig[1, 1], yscale = log10, ylabel = Makie.L"$m(D)$ (kg)", title = Makie.L"Regimes for $ρ_r = 400$ kg/m³")
	ax1R = Makie.Axis(fig[1, 2], yscale = log10, title = Makie.L"Regimes for $F_{rim} = 0.95$")
	ax2L = Makie.Axis(fig[2, 1], yscale = log10, ylabel = Makie.L"$A(D)$ (m²)")
	ax2R = Makie.Axis(fig[2, 2], yscale = log10)
	ax3L = Makie.Axis(fig[3, 1], xlabel = Makie.L"$D$ (mm)", ylabel = Makie.L"$ρ$ (kg/m³)")
	ax3R = Makie.Axis(fig[3, 2], xlabel = Makie.L"$D$ (mm)")
	map(Makie.contents(fig[:, 2])) do ax ax.yticklabelsvisible = false end
	map(Makie.contents(fig[1:2, :])) do ax ax.xticklabelsvisible = false end

	Makie.colgap!(fig.layout, 40) # add space between columns so xticklabels don't overlap

	function lines_and_vlines!((ax1, ax2, ax3), state, color)
		Makie.lines!(ax1, D_range * mm, P3.ice_mass.(state, D_range); color)
		Makie.lines!(ax2, D_range * mm, P3.ice_area.(state, D_range); color)
		Makie.lines!(ax3, D_range * mm, P3.ice_density.(state, D_range); color)

		for ax in [ax1, ax2, ax3]
			Makie.vlines!(ax, state.D_th * mm; linestyle = (:dot, :loose), color = :gray)
			Makie.vlines!(ax, state.D_gr * mm; linestyle = :dash, color)
			Makie.vlines!(ax, state.D_cr * mm; linestyle = :dashdot, color)
		end
	end

	# Make a plot for each rime fraction
	ρ_r₀ = 400.0
	F_rims = [0.0, 0.5, 0.8]
	for (i, F_rim) in enumerate(F_rims)
		color = cl[i]
		state = P3.get_state(params; F_rim, ρ_r = ρ_r₀)
		lines_and_vlines!((ax1L, ax2L, ax3L), state, color)
	end

	# Make a plot for each rime density
	F_rim₀ = 0.95
	ρ_rs = [200.0, 400.0, 800.0]
	for (i, ρ_r) in enumerate(ρ_rs)
		color = cl[i]
		state = P3.get_state(params; F_rim = F_rim₀, ρ_r)
		lines_and_vlines!((ax1R, ax2R, ax3R), state, color)
	end

	Makie.linkyaxes!(ax1R, ax1L)
	Makie.ylims!(ax1R, 1e-14, 1e-4)
	Makie.linkyaxes!(ax2R, ax2L)
	Makie.linkyaxes!(ax3R, ax3L)

	# Add legend
	leg_elems = [
		map(color -> Makie.LineElement(;color), cl), 
		map(linestyle -> Makie.LineElement(;linestyle, linewidth = 1.5, color = :gray), [(:dot, :loose), (:dash, :dense), (:dashdot, :dense)])
	]
	thresh_labels = [Makie.L"$D_{th}$", Makie.L"$D_{gr}$", Makie.L"$D_{cr}$"]
	leg_kwargs = (; orientation = :horizontal, nbanks=10, tellheight = false, tellwidth = false, margin = (10, 10, 10, 10), halign = :left, valign = :top)
	Makie.Legend(fig[1,1], leg_elems, [Makie.latexstring.(Makie.L"F_{rim} = ", F_rims), thresh_labels], ["Rime fraction", "Thresholds"]; leg_kwargs...)
	Makie.Legend(fig[1,2], leg_elems, [Makie.latexstring.(Makie.L"ρ_r = ", ρ_rs), thresh_labels], ["Rime density", "Thresholds"]; leg_kwargs...)

	Makie.resize_to_layout!(fig)
	return fig
end
#! format: on

fig = Makie.with_theme(p3_relations_plot, axis_theme)
Makie.save("P3Scheme_relations.svg", fig)
fig
