import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

# boilerplate code to get free parameter values
include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))
# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "parcel.jl"))
# Boiler plate code to have access to model parameters and constants
FT=Float32
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
prs = cloud_microphysics_parameters(toml_dict)
thermo_params = CMP.thermodynamics_params(prs)

# Initial conditions
Nₐ = FT(0)
Nₗ = FT(100*1e3)
Nᵢ = FT(0)
p = FT(1e5)
T = FT(273.15+20)
ρ = FT(1.2)
qₗ = FT(5*1e-3)
qᵢ = FT(0)
x_sulph = FT(0)
Sᵢ = FT(1.2)
eₛ = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
e = Sᵢ * eₛ
qᵥ = eₛ / (ρ * CMP.R_v(prs) * T)
IC = [Sᵢ, p, T, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]

# Simulation parameters passed into ODE solver
r_nuc = FT(0.5 * 1.e-4 * 1e-6)             # assumed size of nucleated particles
w = FT(1.5)                                # updraft speed
α_m = FT(0.5)                              # accomodation coefficient
const_dt = 0.5                             # model timestep
t_max = 10 * 60
ice_nucleation_modes = []                  #
growth_modes = ["Condensation",]           # switch on condensation
p = (; prs, const_dt, r_nuc, w, α_m, ice_nucleation_modes, growth_modes)

# solve ODE
sol = run_parcel(IC, 0, t_max, p)

# Plot results
fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], xlabel = "Supersaturation [-]", ylabel = "Height [m]")
ax2 = MK.Axis(fig[1, 2], xlabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], xlabel = "q_vap [g/kg]", ylabel = "Height [m]")
ax4 = MK.Axis(fig[2, 2], xlabel = "q_liq [g/kg]")
MK.lines!(ax1, sol[1, :], sol.t * w,)
MK.lines!(ax2, sol[3, :], sol.t * w,)
MK.lines!(ax3, sol[4, :] * 1e3, sol.t * w)
MK.lines!(ax4, sol[5, :] * 1e3, sol.t * w)
MK.save("liquid_only_parcel.svg", fig)
