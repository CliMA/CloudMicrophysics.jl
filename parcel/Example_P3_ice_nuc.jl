import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaParams as CP

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float32
# get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CM.Parameters.WaterProperties(FT)

# Initial conditions
Nₐ = FT(2000)
Nₗ = FT(2000)
Nᵢ = FT(0)
rₗ = FT(1.25e-6)
p₀ = FT(20000)
qᵥ = FT(8.3e-4)
qₗ = FT(Nₗ * 4 / 3 * π * rₗ^3 * wps.ρw / 1.2)
qᵢ = FT(0)
ln_INPC = FT(0)

# Moisture dependent initial conditions
q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
Rₐ = TD.gas_constant_air(tps, q)
Rᵥ = TD.Parameters.R_v(tps)
e = eᵥ(qᵥ, p₀, Rₐ, Rᵥ)

# Simulation parameters passed into ODE solver
r_nuc = FT(1.25e-6)                     # assumed size of nucleated particles
w = FT(0.5)                             # updraft speed
const_dt = FT(0.1)                      # model timestep
t_max = FT(50)
deposition_growth = "Deposition"
size_distribution = "Monodisperse"

# Plotting
fig = MK.Figure(size = (900, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Ice Saturation [-]")
ax2 = MK.Axis(fig[1, 2], ylabel = "ICNC [cm^-3]")
ax7 = MK.Axis(fig[1, 3], ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "Ice Saturation [-]")
ax4 = MK.Axis(fig[2, 2], ylabel = "ICNC [cm^-3]")
ax8 = MK.Axis(fig[2, 3], ylabel = "Temperature [K]")
ax5 = MK.Axis(fig[3, 1], ylabel = "Ice Saturation [-]", xlabel = "Time [s]")
ax6 = MK.Axis(fig[3, 2], ylabel = "ICNC [cm^-3]", xlabel = "Time [s]")
ax9 = MK.Axis(fig[3, 3], ylabel = "Temperature [K]", xlabel = "Time [s]")

ice_nucleation_modes_list = ["P3_dep", "P3_het", "P3_hom"]
T₀_list = [FT(235), FT(235), FT(233.2)]
color = [:blue, :red, :orange]
ax_row = [[ax1, ax2, ax7], [ax3, ax4, ax8], [ax5, ax6, ax9]]

for it in [1, 2, 3]
    mode = ice_nucleation_modes_list[it]
    local T₀ = T₀_list[it]
    local ts = TD.PhaseNonEquil_pTq(tps, p₀, T₀, q)
    local ρₐ = TD.air_density(tps, ts)
    local eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
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
            size_distribution = size_distribution,
        )
    elseif mode == "P3_het"
        local params = parcel_params{FT}(
            const_dt = const_dt,
            r_nuc = r_nuc,
            w = w,
            heterogeneous = mode,
            deposition_growth = deposition_growth,
            size_distribution = size_distribution,
        )
    elseif mode == "P3_hom"
        local params = parcel_params{FT}(
            const_dt = const_dt,
            r_nuc = r_nuc,
            w = w,
            homogeneous = mode,
            deposition_growth = deposition_growth,
            size_distribution = size_distribution,
        )
    end
    # solve ODE
    local sol = run_parcel(IC, FT(0), t_max, params)

    MK.lines!(ax_row[it][1], sol.t, S_i.(tps, sol[3, :], (sol[1, :])), label = mode, color = color[it]) # saturation
    MK.lines!(ax_row[it][2], sol.t, sol[9, :] * 1e-6, color = color[it])                               # ICNC
    MK.lines!(ax_row[it][3], sol.t, sol[3, :], color = color[it])                                      # Temperature
    MK.axislegend(ax_row[it][1], framevisible = false, labelsize = 12, orientation = :horizontal, position = :rt)

    #! format: on
end

MK.save("P3_ice_nuc.svg", fig)
