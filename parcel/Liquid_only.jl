import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CloudMicrophysics.CommonTypes as CMT
import CLIMAParameters as CP

# boilerplate code to get free parameter values
include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))
# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "parcel.jl"))
# Boiler plate code to have access to model parameters and constants
FT = Float64
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
prs = cloud_microphysics_parameters(toml_dict)
thermo_params = CMP.thermodynamics_params(prs)
air_props = CMT.AirProperties(FT)

# Constants
ρₗ = FT(CMP.ρ_cloud_liq(prs))
R_v = FT(CMP.R_v(prs))
R_d = FT(CMP.R_d(prs))

# Initial conditions
Nₐ = FT(0)
Nₗ = FT(200 * 1e6)
Nᵢ = FT(0)
r₀ = FT(8e-6)
p₀ = FT(800 * 1e2)
T₀ = FT(273.15 + 7.0)
x_sulph = FT(0)

eₛ = TD.saturation_vapor_pressure(thermo_params, T₀, TD.Liquid())
e = eₛ
ϵ = R_d / R_v

# mass per volume for dry air, vapor and liquid
md_v = (p₀ - e) / R_d / T₀
mv_v = e / R_v / T₀
ml_v = Nₗ * 4 / 3 * π * ρₗ * r₀^3

qᵥ = mv_v / (md_v + mv_v + ml_v)
qₗ = ml_v / (md_v + mv_v + ml_v)
qᵢ = FT(0)

# Moisture dependent initial conditions
q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
ts = TD.PhaseNonEquil_pTq(thermo_params, p₀, T₀, q)
ρₐ = TD.air_density(thermo_params, ts)

Rₐ = TD.gas_constant_air(thermo_params, q)
Sₗ = FT(e / eₛ)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]

# Simulation parameters passed into ODE solver
r_nuc = FT(0.5 * 1.e-4 * 1e-6)             # assumed size of nucleated particles
w = FT(10)                                 # updraft speed
α_m = FT(0.5)                              # accomodation coefficient
const_dt = FT(0.5)                         # model timestep
t_max = FT(20)
ice_nucleation_modes = []                  # no freezing
growth_modes_list = [["Condensation"], ["Condensation_DSD"]] # switch on condensation

# Data from Rogers(1975) Figure 1
# https://www.tandfonline.com/doi/abs/10.1080/00046973.1975.9648397
#! format: off
Rogers_time_supersat = [0.0645, 0.511, 0.883, 1.4, 2.07, 2.72, 3.24, 3.89, 4.53, 5.87, 7.16, 9.79, 16.0, 19.8]
Rogers_supersat = [0.0268, 0.255, 0.393, 0.546, 0.707, 0.805, 0.863, 0.905, 0.938, 0.971, 0.978, 0.963, 0.910, 0.885]
Rogers_time_radius = [0.561, 2, 3.99, 10.7, 14.9, 19.9]
Rogers_radius = [8.0, 8.08, 8.26, 8.91, 9.26, 9.68]
#! format: on

fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Supersaturation [%]")
ax2 = MK.Axis(fig[3, 1], xlabel = "Time [s]", ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "q_vap [g/kg]")
ax4 = MK.Axis(fig[2, 2], xlabel = "Time [s]", ylabel = "q_liq [g/kg]")
ax5 = MK.Axis(fig[1, 2], ylabel = "radius [μm]")

MK.lines!(ax1, Rogers_time_supersat, Rogers_supersat, label = "Rogers_1975")
MK.lines!(ax5, Rogers_time_radius, Rogers_radius)

for growth_modes in growth_modes_list
    p = (;
        prs,
        air_props,
        thermo_params,
        const_dt,
        r_nuc,
        w,
        α_m,
        ice_nucleation_modes,
        growth_modes,
    )
    # solve ODE
    sol = run_parcel(IC, FT(0), t_max, p)

    gm = growth_modes[1]

    # Plot results
    MK.lines!(ax1, sol.t, (sol[1, :] .- 1) * 100.0, label = gm)
    MK.lines!(ax2, sol.t, sol[3, :])
    MK.lines!(ax3, sol.t, sol[4, :] * 1e3)
    MK.lines!(ax4, sol.t, sol[5, :] * 1e3)

    sol_Nₗ = sol[8, :]
    sol_qₗ = sol[5, :]
    if gm == "Condensation"
        rₗ = cbrt.(sol_qₗ ./ sol_Nₗ ./ (4 / 3 * π) / ρₗ * ρₐ)
    end
    if gm == "Condensation_DSD"
        λ = cbrt.(32 .* π .* sol_Nₗ ./ sol_qₗ * ρₗ / ρₐ)
        rₗ = 2 ./ λ
    end
    MK.lines!(ax5, sol.t, rₗ * 1e6)
end

#fig[3, 2] = MK.Legend(fig, ax1, framevisible = false)
MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 12,
    orientation = :horizontal,
    nbanks = 2,
    position = :rb,
)

MK.save("liquid_only_parcel.svg", fig)
