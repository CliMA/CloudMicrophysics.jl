import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CloudMicrophysics.CloudDiagnostics as CD
import CloudMicrophysics.Parameters as CMP
import ClimaParams as CP

FT = Float64

include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))


# choosing a value for the ice relaxation timescale
τ = FT(10)
@info("using τ value ", τ)

# Get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)
# Constants
ρₗ = wps.ρw
ρᵢ = wps.ρi
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)

# Initial conditions
Nₐ = FT(0)
Nₗ = FT(200 * 1e6)
Nᵢ = FT(1e6)
r₀ = FT(8e-6)
p₀ = FT(800 * 1e2)
T₀ = FT(243)
ln_INPC = FT(0)
e_sat = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
Sₗ = FT(1)
e = Sₗ * e_sat
md_v = (p₀ - e) / R_d / T₀
mv_v = e / R_v / T₀
ml_v = Nₗ * 4 / 3 * FT(π) * ρₗ * r₀^3
mi_v = Nᵢ * 4 / 3 * FT(π) * ρᵢ * r₀^3
qᵥ = mv_v / (md_v + mv_v + ml_v + mi_v)
qₗ = ml_v / (md_v + mv_v + ml_v + mi_v)
qᵢ = mi_v / (md_v + mv_v + ml_v + mi_v)

override_file = Dict(
    "condensation_evaporation_timescale" =>
        Dict("value" => τ, "type" => "float"),
)

liquid_toml_dict = CP.create_toml_dict(FT; override_file)
liquid = CMP.CloudLiquid(liquid_toml_dict)

override_file = Dict(
    "sublimation_deposition_timescale" =>
        Dict("value" => τ, "type" => "float"),
)

ice_toml_dict = CP.create_toml_dict(FT; override_file)

ice = CMP.CloudIce(ice_toml_dict)
@info("relaxations:", liquid.τ_relax, ice.τ_relax)

IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]
simple = false

# Simulation parameters passed into ODE solver
w = FT(1)                                 # updraft speed
const_dt = FT(0.001)                         # model timestep
t_max = FT(20)#FT(const_dt*1)
size_distribution_list = ["Monodisperse"]

condensation_growth = "NonEq_Condensation"
deposition_growth = "NonEq_Deposition"

# Setup the plots
fig = MK.Figure(size = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Liquid Supersaturation [%]")
ax2 = MK.Axis(fig[2, 1], ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[3, 1], xlabel = "Time [s]", ylabel = "q_liq [g/kg]")
ax4 = MK.Axis(fig[1, 2], ylabel = "Ice Supersaturation [%]")
ax5 = MK.Axis(fig[2, 2], ylabel = "q_vap [g/kg]")
ax6 = MK.Axis(fig[3, 2], xlabel = "Time [s]", ylabel = "q_ice [g/kg]")

for DSD in size_distribution_list
    local params = parcel_params{FT}(
        liq_size_distribution = DSD,
        condensation_growth = condensation_growth,
        deposition_growth = deposition_growth,
        const_dt = const_dt,
        w = w,
        liquid = liquid,
        ice = ice,
    )

    # solve ODE
    local sol = run_parcel(IC, FT(0), t_max, params)

    # Plot results
    MK.lines!(ax1, sol.t, (sol[1, :] .- 1) * 100.0)
    MK.lines!(ax2, sol.t, sol[3, :])
    MK.lines!(ax3, sol.t, sol[5, :] * 1e3)
    MK.lines!(ax4, sol.t, (S_i.(tps, sol[3, :], sol[1, :]) .- 1) * 100.0)
    MK.lines!(ax5, sol.t, sol[4, :] * 1e3)
    MK.lines!(ax6, sol.t, sol[6, :] * 1e3)


    sol_Nₗ = sol[8, :]
    sol_Nᵢ = sol[9, :]
    sol_T = sol[3, :]
    sol_p = sol[2, :]
    sol_qᵥ = sol[4, :]
    sol_qₗ = sol[5, :]
    sol_qᵢ = sol[6, :]

    q = TD.PhasePartition.(sol_qᵥ + sol_qₗ + sol_qᵢ, sol_qₗ, sol_qᵢ)

    local q = TD.PhasePartition.(sol_qᵥ + sol_qₗ + sol_qᵢ, sol_qₗ, sol_qᵢ)
    local ts = TD.PhaseNonEquil_pTq.(tps, sol_p, sol_T, q)
    local ρₐ = TD.air_density.(tps, ts)

end

MK.save("noneq_parcel.svg", fig)
nothing
