import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import ClimaParams as CP

FT = Float64

include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))

# Get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)
liquid = CMP.CloudLiquid(FT)
ice = CMP.CloudIce(FT)
# Constants
ρₗ = wps.ρw
ρᵢ = wps.ρi
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)

# Initial conditions
Nₐ = FT(0)
Nₗ = FT(200 * 1e6)
Nᵢ = FT(200 * 1e3)
r₀ = FT(8e-6)
p₀ = FT(800 * 1e2)
T₀ = FT(273.15 - 30.0)
ln_INPC = FT(0)
e_sat = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
Sₗ = FT(0.8)
e = Sₗ * e_sat
md_v = (p₀ - e) / R_d / T₀
mv_v = e / R_v / T₀
ml_v = Nₗ * 4 / 3 * FT(π) * ρₗ * r₀^3
mi_v = Nᵢ * 4 / 3 * FT(π) * ρᵢ * r₀^3
qᵥ = mv_v / (md_v + mv_v + ml_v + mi_v)
qₗ = ml_v / (md_v + mv_v + ml_v + mi_v)
qᵢ = mi_v / (md_v + mv_v + ml_v + mi_v)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]
simple = false

# Simulation parameters passed into ODE solver
w = FT(1)                                 # updraft speed
const_dt = FT(20)                         # model timestep
t_max = FT(20)#FT(const_dt*1)
size_distribution_list = ["Monodisperse", "Gamma"]

if simple
    condensation_growth = "NonEq_Condensation_Simple_Morrison"
    deposition_growth = "NonEq_Deposition_Simple_Morrison"
else
    condensation_growth = "NonEq_Condensation"
    deposition_growth = "NonEq_Deposition"
end

# Setup the plots
fig = MK.Figure(size = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Liquid Supersaturation [%]")
ax2 = MK.Axis(fig[2, 1], ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[3, 1], xlabel = "Time [s]", ylabel = "q_liq [g/kg]")
ax4 = MK.Axis(fig[1, 2], ylabel = "Ice Supersaturation [%]")
ax5 = MK.Axis(fig[2, 2], ylabel = "q_vap [g/kg]")
ax6 = MK.Axis(fig[3, 2], xlabel = "Time [s]", ylabel = "q_ice [g/kg]")
#ax7 = MK.Axis(fig[1, 3], xlabel = "Time [s]", ylabel = "internal energy")

for DSD in size_distribution_list
    local params = parcel_params{FT}(
        liq_size_distribution = DSD,
        ice_size_distribution = DSD,
        condensation_growth = condensation_growth,
        deposition_growth = deposition_growth,
        const_dt = const_dt,
        w = w,
        liquid = liquid,#,
        #ice = ice
    )
    # solve ODE
    local sol = run_parcel(IC, FT(0), t_max, params)

    # Plot results
    MK.lines!(ax1, sol.t, (sol[1, :] .- 1) * 100.0, label = DSD)
    MK.lines!(ax2, sol.t, sol[3, :])
    MK.lines!(ax3, sol.t, sol[5, :] * 1e3)
    MK.lines!(ax4, sol.t, (S_i.(tps, sol[3, :], sol[1, :]) .- 1) * 100.0, label = DSD)
    MK.lines!(ax5, sol.t, sol[4, :] * 1e3)
    MK.lines!(ax6, sol.t, sol[6, :] * 1e3)


    sol_Nₗ = sol[8, :]
    sol_Nᵢ = sol[9, :]
    # Compute the current air density
    sol_T = sol[3, :]
    sol_p = sol[2, :]
    sol_qᵥ = sol[4, :]
    sol_qₗ = sol[5, :]
    sol_qᵢ = sol[6, :]

    # calculating the saturation vapor pressure from S:
    #e_test = sol_qᵥ.*sol_p .* (R_v/(R_v+R_d))
    q = TD.PhasePartition.(sol_qᵥ + sol_qₗ + sol_qᵢ, sol_qₗ, sol_qᵢ)
    e_test = sol_qᵥ .* sol_p .* (TD.Parameters.R_v(tps) ./ TD.gas_constant_air.(tps, q))
    #e_test = sol_qᵥ.*sol_p .* (TD.Parameters.R_v(tps)/(TD.Parameters.R_v(tps)+TD.Parameters.R_d(tps)))
    e_sat_fromS = e_test ./ sol[1, :]

    # calculating it from T:
    e_sat_fromT = TD.saturation_vapor_pressure.(tps, sol_T, TD.Liquid())

    rel_error = (e_sat_fromT .- e_sat_fromS) ./ (e_sat_fromT)

    # plot them:
    #MK.lines!(ax7, sol.t, rel_error)
    #MK.lines!(ax7, sol.t, e_sat_fromT, label="from T")
    #MK.lines!(ax7, sol.t, e_sat_fromS, label="from S")

    q_sat_fromS = sol_qᵥ .* sol[1, :] #- sol_qᵥ

    #MK.lines!(ax8, sol.t, q_sat_fromS)

    # calculating internal energy as a sanity check

    int_energy = TD.internal_energy.(tps, sol_T, q)

    #MK.lines!(ax7, sol.t, int_energy)

    local q = TD.PhasePartition.(sol_qᵥ + sol_qₗ + sol_qᵢ, sol_qₗ, sol_qᵢ)
    local ts = TD.PhaseNonEquil_pTq.(tps, sol_p, sol_T, q)
    local ρₐ = TD.air_density.(tps, ts)
end

MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 12,
    orientation = :horizontal,
    nbanks = 2,
    position = :rb,
)

#MK.axislegend(
#    ax7,
#    framevisible = false,
#    labelsize = 12,
#    orientation = :horizontal,
#    nbanks = 2,
#    position = :rb,
#)

if simple
    #MK.save("/Users/oliviaalcabes/Documents/research/microphysics/parcel_sims/timestep_experiment_simple_morrison/both.png", fig)
    MK.save("ice_noneq_parcel_simple_morrison.svg", fig)
else
    #MK.save("ice_noneq_parcel_morrison.svg", fig)
    MK.save("/Users/oliviaalcabes/Documents/research/microphysics/parcel_sims/timestep_experiment_morrison/both.png", fig)
end