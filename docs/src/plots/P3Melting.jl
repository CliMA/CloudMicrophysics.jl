import CairoMakie: CairoMakie, Makie
CairoMakie.activate!(type = "svg")

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI
FT = Float64

# parameters
params = CMP.ParametersP3(FT; slope_law = :constant)
vel = CMP.Chen2022VelType(FT)
aps = CMP.AirProperties(FT)
tps = TDI.PS(FT)

# model time step (for limiting)
dt = FT(1)

# initial ice content and number concentration
Lᵢ = FT(1e-4)
Nᵢ = FT(2e5)

# tested temperature range
ΔT_range = range(1e-4, stop = 0.025, length = 1000)

# limiters to not melt more mass and number than we have
max_dLdt = Lᵢ / dt
max_dNdt = Nᵢ / dt

#! format: off

ρₐ = FT(1.2)
Fᵣ = FT(0.8)
ρᵣ = FT(800)
label1 = "ρₐ=$ρₐ kg/m³, Fᵣ=$Fᵣ, ρᵣ=$(Int(ρᵣ)) kg/m³"
state = P3.get_state(params; F_rim = Fᵣ, ρ_rim = ρᵣ, L_ice = Lᵢ, N_ice = Nᵢ)
logλ = P3.get_distribution_logλ(state)
melt1 = P3.ice_melt.(vel, aps, tps, params.T_freeze .+ ΔT_range, ρₐ, state, logλ)
dLdt1 = getfield.(melt1, :dLdt)
dNdt1 = getfield.(melt1, :dNdt)

Fᵣ = FT(0.2)
label2 = "ρₐ=$ρₐ kg/m³, Fᵣ=$Fᵣ, ρᵣ=$(Int(ρᵣ)) kg/m³"
state = P3.get_state(params; F_rim = Fᵣ, ρ_rim = ρᵣ, L_ice = Lᵢ, N_ice = Nᵢ)
logλ = P3.get_distribution_logλ(state)
melt2 = P3.ice_melt.(vel, aps, tps, params.T_freeze .+ ΔT_range, ρₐ, state, logλ)
dLdt2 = getfield.(melt2, :dLdt)
dNdt2 = getfield.(melt2, :dNdt)

ρᵣ = FT(200)
label3 = "ρₐ=$ρₐ kg/m³, Fᵣ=$Fᵣ, ρᵣ=$(Int(ρᵣ)) kg/m³"
state = P3.get_state(params; F_rim = Fᵣ, ρ_rim = ρᵣ, L_ice = Lᵢ, N_ice = Nᵢ)
logλ = P3.get_distribution_logλ(state)
melt3 = P3.ice_melt.(vel, aps, tps, params.T_freeze .+ ΔT_range, ρₐ, state, logλ)
dLdt3 = getfield.(melt3, :dLdt)
dNdt3 = getfield.(melt3, :dNdt)

ρₐ = FT(0.5)
label4 = "ρₐ=$ρₐ kg/m³, Fᵣ=$Fᵣ, ρᵣ=$(Int(ρᵣ)) kg/m³"
state = P3.get_state(params; F_rim = Fᵣ, ρ_rim = ρᵣ, L_ice = Lᵢ, N_ice = Nᵢ)
logλ = P3.get_distribution_logλ(state)
melt4 = P3.ice_melt.(vel, aps, tps, params.T_freeze .+ ΔT_range, ρₐ, state, logλ)
dLdt4 = getfield.(melt4, :dLdt)
dNdt4 = getfield.(melt4, :dNdt)

# plotting
Makie.with_theme(Makie.theme_minimal(), fontsize = 22, linewidth = 3) do
    fig = Makie.Figure(size = (800, 600))

    ax1 = Makie.Axis(fig[1, 1]; 
        xlabel = "T [°C]", ylabel = "ice mass melting rate [g/m³/s]",
        title = "P3 ice melting",
    )
    ax2 = Makie.Axis(fig[1, 1]; 
        xlabel = "T [°C]", ylabel = "ice number melting rate [1/cm³/s]",
        yaxisposition = :right, rightspinevisible = true,
    )

    l_max_dLdt = Makie.hlines!(ax1, max_dLdt * 1e3;  color = :gray, linewidth = 1)
    l_max_dNdt = Makie.hlines!(ax2, max_dNdt * 1e-6; color = :gray, linewidth = 1, label = "limit")

    l_dLdt1 = Makie.lines!(ax1, ΔT_range,  dLdt1 * 1e3;  color = :skyblue)
    l_dNdt1 = Makie.lines!(ax2, ΔT_range,  dNdt1 * 1e-6; color = :skyblue, label = label1)

    l_dLdt2 = Makie.lines!(ax1, ΔT_range,  dLdt2 * 1e3;  color = :blue3)
    l_dNdt2 = Makie.lines!(ax2, ΔT_range,  dNdt2 * 1e-6; color = :blue3, label = label2)

    l_dLdt3 = Makie.lines!(ax1, ΔT_range,  dLdt3 * 1e3;  color = :orchid)
    l_dNdt3 = Makie.lines!(ax2, ΔT_range,  dNdt3 * 1e-6; color = :orchid, label = label3)

    l_dLdt4 = Makie.lines!(ax1, ΔT_range,  dLdt4 * 1e3;  color = :purple)
    l_dNdt4 = Makie.lines!(ax2, ΔT_range,  dNdt4 * 1e-6; color = :purple, label = label4)

    Makie.axislegend(ax2; position = :rb, framevisible = false)
    Makie.save("P3_ice_melt.png", fig)
end
