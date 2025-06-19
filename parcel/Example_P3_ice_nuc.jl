import OrdinaryDiffEq as ODE
import CairoMakie as MK

import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ClimaParams as CP

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float32
# get free parameters
tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
wps = CM.Parameters.WaterProperties(FT)

# Initial conditions
Nₐ = FT(2000)
Nₗ = FT(2000)
Nᵢ = FT(0)
rₗ = FT(1.25e-6)
p₀ = FT(20000)
qᵥ = FT(5e-4)
qₗ = FT(Nₗ * 4 / 3 * π * rₗ^3 * wps.ρw / 1.2)
qᵢ = FT(0)
qₜ = qᵥ + qₗ + qᵢ
ln_INPC = FT(0)

# Moisture dependent initial conditions
Rₐ = TDI.Rₘ(tps, qₜ, qₗ, qᵢ)
Rᵥ = TDI.Rᵥ(tps)
e = eᵥ(qᵥ, p₀, Rₐ, Rᵥ)

# Simulation parameters passed into ODE solver
r_nuc = FT(1.25e-6)                     # assumed size of nucleated particles
w = FT(0.5)                             # updraft speed
const_dt = FT(0.1)                      # model timestep
t_max = FT(50)
deposition_growth = "Deposition"
size_distribution = "Monodisperse"

# Plotting
fig = MK.Figure(size = (1000, 350), fontsize = 20)
ax1 = MK.Axis(fig[1, 1], ylabel = "ICNC [cm^-3]", xlabel = "Time [s]")
ax2 = MK.Axis(fig[1, 2], ylabel = "ICNC [cm^-3]", xlabel = "Time [s]")
ax3 = MK.Axis(fig[1, 3], ylabel = "ICNC [cm^-3]", xlabel = "Time [s]")

ice_nucleation_modes_list = ["P3_dep", "P3_het", "P3_hom"]
T₀_list = [FT(235), FT(235), FT(233.2)]
color = [:blue, :red, :orange]
ax_list = [ax1, ax2, ax3]

for it in [1, 2, 3]
    mode = ice_nucleation_modes_list[it]
    local T₀ = T₀_list[it]
    local ρₐ = TDI.air_density(tps, T₀, p₀, qᵥ + qₗ + qᵢ, qₗ, qᵢ)
    local eₛ = TDI.saturation_vapor_pressure_over_liquid(tps, T₀)
    local Sₗ = FT(e / eₛ)
    local IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]

    #! format: off
    if mode == "P3_dep"
        local params = parcel_params{FT}(
            const_dt = const_dt,
            r_nuc = r_nuc,
            w = w,
            deposition = mode,
            deposition_growth = deposition_growth,
            ice_size_distribution = size_distribution,
        )
    elseif mode == "P3_het"
        local params = parcel_params{FT}(
            const_dt = const_dt,
            r_nuc = r_nuc,
            w = w,
            heterogeneous = mode,
            deposition_growth = deposition_growth,
            ice_size_distribution = size_distribution,
        )
    elseif mode == "P3_hom"
        local params = parcel_params{FT}(
            const_dt = const_dt,
            r_nuc = r_nuc,
            w = w,
            homogeneous = mode,
            deposition_growth = deposition_growth,
            ice_size_distribution = size_distribution,
        )
    end
    # solve ODE
    local sol = run_parcel(IC, FT(0), t_max, params)

    MK.lines!(ax_list[it], sol.t, sol[9, :] * 1e-6, label = mode, color = color[it]) # saturation
    MK.axislegend(ax_list[it], framevisible = false, labelsize = 12, orientation = :horizontal, position = :rb)

    #! format: on
end

MK.save("P3_ice_nuc.svg", fig)
nothing
