import CairoMakie: CairoMakie, Makie
CairoMakie.activate!(type = "svg")

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI
FT = Float64

# parameters
params = CMP.ParametersP3(FT)
vel = CMP.Chen2022VelType(FT)
quad = CMP.Microphysics2MParams(FT; with_ice = true).ice.quad  # the package default rule
aps = CMP.AirProperties(FT)
tps = TDI.PS(FT)

# ice content [kg/m³] and number concentration [1/m³]
Lᵢ = FT(1e-4)
Nᵢ = FT(2e5)

# temperature above freezing [K]
ΔT_range = range(FT(0), FT(2), length = 200)

# (ρₐ [kg/m³], Fᵣ, ρᵣ [kg/m³]) and line color
cases = (
    (FT(1.2), FT(0.8), FT(800), :skyblue),
    (FT(1.2), FT(0.2), FT(800), :blue3),
    (FT(1.2), FT(0.2), FT(200), :orchid),
    (FT(0.5), FT(0.2), FT(200), :purple),
)

Makie.with_theme(Makie.theme_minimal(), fontsize = 20, linewidth = 3) do
    fig = Makie.Figure(size = (1000, 500))
    ax_L = Makie.Axis(fig[1, 1]; xlabel = "T - T_freeze [K]", ylabel = "ice mass melting rate [g/m³/s]")
    ax_N = Makie.Axis(fig[1, 2]; xlabel = "T - T_freeze [K]", ylabel = "ice number melting rate [1/cm³/s]")
    for (ρₐ, Fᵣ, ρᵣ, color) in cases
        state = P3.P3State(params, Lᵢ, Nᵢ, Fᵣ, ρᵣ)
        logλ = P3.get_distribution_logλ(state)
        melt = P3.ice_melt.(vel, aps, tps, params.T_freeze .+ ΔT_range, ρₐ, state, logλ; quad)
        label = "ρₐ = $ρₐ kg/m³, Fᵣ = $Fᵣ, ρᵣ = $(Int(ρᵣ)) kg/m³"
        Makie.lines!(ax_L, ΔT_range, getfield.(melt, :dLdt) * 1e3; color, label)
        Makie.lines!(ax_N, ΔT_range, getfield.(melt, :dNdt) * 1e-6; color)
    end
    Makie.Legend(fig[2, 1:2], ax_L; orientation = :horizontal, nbanks = 2, framevisible = false)
    Makie.save("P3_ice_melt.svg", fig)
end
