import CairoMakie: Makie
import CloudMicrophysics as CM
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

# Make heatmaps of D_gr, D_cr, ρ_g, ρ_d as a function of F_rim and ρ_r
FT = Float64
params = CMP.ParametersP3(FT)

F_rims = FT.(0.0:0.01:0.9)
ρ_rs = FT.(200.0:10.0:900.0)

get_state2(F_rim, ρ_r) = P3.get_state(params; F_rim, ρ_r)

states = get_state2.(F_rims, ρ_rs')

fig = Makie.Figure(size = (1200, 900), figure_padding = 20)

D_th = round(states[1].D_th * 1e3, digits = 4) # mm
# Make a plot for each threshold
thresh_keys = (:D_gr, :D_cr)
for (i, key) in enumerate(thresh_keys)
    title = string(key) * if i == 1
        " (D_th = $(D_th) mm)"
    else
        ""
    end
    ax = Makie.Axis(fig[i, 1], xlabel = "F_rim", ylabel = "ρ_r", title = title)
    threshold = getproperty.(states, key) * 1e3
    hm = Makie.heatmap!(ax, F_rims, ρ_rs, threshold)
    Makie.Colorbar(fig[i, 2], hm, label = "mm")
end

# Make a plot for each density
ρ_keys = (:ρ_g, :ρ_d)
for (i, key) in enumerate(ρ_keys)
    ax = Makie.Axis(fig[i, 3], xlabel = "F_rim", ylabel = "ρ_r", title = string(key))
    density = if key == :ρ_g
        getproperty.(states, key)
    else # :ρ_d
        P3.get_ρ_d.(params.mass, F_rims, ρ_rs')
    end
    hm = Makie.heatmap!(ax, F_rims, ρ_rs, density)
    Makie.Colorbar(fig[i, 4], hm, label = "kg/m³")
end

Makie.save("P3Thresholds.svg", fig)
fig
