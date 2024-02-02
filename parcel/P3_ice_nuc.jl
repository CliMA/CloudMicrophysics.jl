import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "parcel.jl"))
FT = Float32
# get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
aps = CMP.AirProperties(FT)
wps = CMP.WaterProperties(FT)
ip = CMP.IceNucleationParameters(FT)

# Constants
ρₗ = wps.ρw
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)

# Initial conditions
Nₐ = FT(2000)
Nₗ = FT(2000)
Nᵢ = FT(0)
rₗ = FT(1.25e-6)
p₀ = FT(20000)
T₀_dep = FT(235)
T₀_het = FT(235)
T₀_hom = FT(233.2)
qᵥ = FT(8.3e-4)
qₗ = FT(Nₗ * 4 / 3 * π * rₗ^3 * ρₗ / 1.2)
qᵢ = FT(0)
x_sulph = FT(0)

# Moisture dependent initial conditions
q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
ts_dep = TD.PhaseNonEquil_pTq(tps, p₀, T₀_dep, q)
ts_het = TD.PhaseNonEquil_pTq(tps, p₀, T₀_het, q)
ts_hom = TD.PhaseNonEquil_pTq(tps, p₀, T₀_hom, q)
ρₐ_dep = TD.air_density(tps, ts_dep)
ρₐ_het = TD.air_density(tps, ts_het)
ρₐ_hom = TD.air_density(tps, ts_hom)
Rₐ = TD.gas_constant_air(tps, q)
eₛ_dep = TD.saturation_vapor_pressure(tps, T₀_dep, TD.Liquid())
eₛ_het = TD.saturation_vapor_pressure(tps, T₀_het, TD.Liquid())
eₛ_hom = TD.saturation_vapor_pressure(tps, T₀_hom, TD.Liquid())
e = qᵥ * p₀ * R_v / Rₐ

Sₗ_dep = FT(e / eₛ_dep)
Sₗ_het = FT(e / eₛ_het)
Sₗ_hom = FT(e / eₛ_hom)
IC_dep = [Sₗ_dep, p₀, T₀_dep, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]
IC_het = [Sₗ_het, p₀, T₀_het, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]
IC_hom = [Sₗ_hom, p₀, T₀_hom, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]

ξ(T) =
    TD.saturation_vapor_pressure(tps, T, TD.Liquid()) /
    TD.saturation_vapor_pressure(tps, T, TD.Ice())
S_i(T, S_liq) = ξ(T) * S_liq

# Simulation parameters passed into ODE solver
r_nuc = FT(1.25e-6)                     # assumed size of nucleated particles
w = FT(0.5)                             # updraft speed
α_m = FT(0.5)                           # accomodation coefficient
const_dt = FT(0.1)                      # model timestep
t_max = FT(50)
aerosol = []
ice_nucleation_modes_list = [["P3_Deposition"], ["P3_het"], ["P3_hom"]]
growth_modes = ["Deposition"]
droplet_size_distribution_list = [["Monodisperse"]]

# Plotting
fig = MK.Figure(resolution = (900, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Ice Saturation [-]")
ax2 = MK.Axis(fig[1, 2], ylabel = "ICNC [cm^-3]")
ax7 = MK.Axis(fig[1, 3], ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "Ice Saturation [-]")
ax4 = MK.Axis(fig[2, 2], ylabel = "ICNC [cm^-3]")
ax8 = MK.Axis(fig[2, 3], ylabel = "Temperature [K]")
ax5 = MK.Axis(fig[3, 1], ylabel = "Ice Saturation [-]", xlabel = "Time [s]")
ax6 = MK.Axis(fig[3, 2], ylabel = "ICNC [cm^-3]", xlabel = "Time [s]")
ax9 = MK.Axis(fig[3, 3], ylabel = "Temperature [K]", xlabel = "Time [s]")

for ice_nucleation_modes in ice_nucleation_modes_list
    nuc_mode = ice_nucleation_modes[1]
    droplet_size_distribution = droplet_size_distribution_list[1]
    p = (;
        wps,
        aps,
        tps,
        ip,
        const_dt,
        r_nuc,
        w,
        α_m,
        aerosol,
        ice_nucleation_modes,
        growth_modes,
        droplet_size_distribution,
    )
    # solve ODE
    #! format: off
    if "P3_Deposition" in ice_nucleation_modes
        sol = run_parcel(IC_dep, FT(0), t_max, p)
        MK.lines!(ax1, sol.t, S_i.(sol[3, :], (sol[1, :])), label = nuc_mode) # saturation
        MK.lines!(ax2, sol.t, sol[9, :] * 1e-6)                               # ICNC
        MK.lines!(ax7, sol.t, sol[3, :])                                      # Temperature
        MK.axislegend(ax1, framevisible = false, labelsize = 12, orientation = :horizontal, position = :rt)
    elseif "P3_het" in ice_nucleation_modes
        sol = run_parcel(IC_het, FT(0), t_max, p)
        MK.lines!(ax3, sol.t, S_i.(sol[3, :], (sol[1, :])), label = nuc_mode, color = :red) # saturation
        MK.lines!(ax4, sol.t, sol[9, :] * 1e-6, color = :red)                               # ICNC
        MK.lines!(ax8, sol.t, sol[3, :], color = :red)                                                    # Temperature
        MK.axislegend(ax3, framevisible = false, labelsize = 12, orientation = :horizontal, position = :lt)
    elseif "P3_hom" in ice_nucleation_modes
        sol = run_parcel(IC_hom, FT(0), t_max, p)
        MK.lines!(ax5, sol.t, S_i.(sol[3, :], (sol[1, :])), label = nuc_mode, color = :orange) # saturation
        MK.lines!(ax6, sol.t, sol[9, :] * 1e-6, color = :orange)                               # ICNC
        MK.lines!(ax9, sol.t, sol[3, :], color = :orange)                                                       # Temperature
        MK.axislegend(ax5, framevisible = false, labelsize = 12, orientation = :horizontal, position = :lt)
    end
    #! format: on
end

MK.save("P3_ice_nuc.svg", fig)
