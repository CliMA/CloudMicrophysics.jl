import OrdinaryDiffEq as ODE
import CairoMakie as MK

import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ClimaParams as CP

FT = Float64

include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))

# choosing a value for the ice relaxation timescale
τ = FT(10)
@info("using τ value ", τ)
override_file = Dict(
    "condensation_evaporation_timescale" =>
        Dict("value" => τ, "type" => "float"),
    "sublimation_deposition_timescale" =>
        Dict("value" => τ, "type" => "float"),
)
override_toml_dict = CP.create_toml_dict(FT; override_file)
liquid = CMP.CloudLiquid(override_toml_dict)
ice = CMP.CloudIce(override_toml_dict)
@info("relaxations:", liquid.τ_relax, ice.τ_relax)

# Get free parameters
tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)
# Constants
ρₗ = wps.ρw
ρᵢ = wps.ρi
R_v = TDI.Rᵥ(tps)
R_d = TDI.Rd(tps)

# Initial conditions
Nₐ = FT(0)
Nₗ = FT(200 * 1e6)
Nᵢ = FT(1e6)
r₀ₗ = FT(1e-6)
r₀ᵢ = FT(8e-6)
p₀ = FT(800 * 1e2)
T₀ = FT(243)
ln_INPC = FT(0)
e_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T₀)
Sₗ = FT(1)
e = Sₗ * e_sat
md_v = (p₀ - e) / R_d / T₀
mv_v = e / R_v / T₀
ml_v = Nₗ * 4 / 3 * FT(π) * ρₗ * r₀ₗ^3
mi_v = Nᵢ * 4 / 3 * FT(π) * ρᵢ * r₀ᵢ^3
qᵥ = mv_v / (md_v + mv_v + ml_v + mi_v)
qₗ = ml_v / (md_v + mv_v + ml_v + mi_v)
qᵢ = mi_v / (md_v + mv_v + ml_v + mi_v)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]

# Simulation parameters passed into ODE solver
w = FT(1)                     # updraft speed
const_dt = FT(0.001)          # model timestep
t_max = FT(20)
condensation_growth = "NonEq_Condensation"
deposition_growth = "NonEq_Deposition"

# Setup the plots
fig = MK.Figure(size = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Liquid Supersaturation [%]")
ax2 = MK.Axis(fig[2, 1], ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[3, 1], xlabel = "Time [s]", ylabel = "q_lcl [g/kg]")
ax4 = MK.Axis(fig[1, 2], ylabel = "Ice Supersaturation [%]")
ax5 = MK.Axis(fig[2, 2], ylabel = "q_vap [g/kg]")
ax6 = MK.Axis(fig[3, 2], xlabel = "Time [s]", ylabel = "q_icl [g/kg]")

params = parcel_params{FT}(
    condensation_growth = condensation_growth,
    deposition_growth = deposition_growth,
    const_dt = const_dt,
    w = w,
    liquid = liquid,
    ice = ice,
)

# solve ODE
sol = run_parcel(IC, FT(0), t_max, params)

# Plot results
MK.lines!(ax1, sol.t, (sol[1, :] .- 1) * 100.0)
MK.lines!(ax2, sol.t, sol[3, :])
MK.lines!(ax3, sol.t, sol[5, :] * 1e3)
MK.lines!(ax4, sol.t, (S_i.(tps, sol[3, :], sol[1, :]) .- 1) * 100.0)
MK.lines!(ax5, sol.t, sol[4, :] * 1e3)
MK.lines!(ax6, sol.t, sol[6, :] * 1e3)

MK.save("noneq_parcel.svg", fig)
