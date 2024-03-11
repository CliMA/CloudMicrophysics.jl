import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import Cloudy as CL
import Cloudy.ParticleDistributions as CPD
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsFlexible as CMF
import ClimaParams as CP

FT = Float64

"""
    ODE problem definitions
"""
function parcel_model_cloudy(dY, Y, p, t)
    # Numerical precision used in the simulation
    FT = eltype(Y[5:end])
    # Simulation parameters
    (; wps, tps, aps, w, clinfo) = p
    # Y values
    Sₗ = Y[1]
    p_air = Y[2]
    T = Y[3]
    qᵥ = Y[4]
    moments = Y[5:end]

    # Constants
    Rᵥ = TD.Parameters.R_v(tps)
    grav = TD.Parameters.grav(tps)
    ρₗ = wps.ρw

    # Get thermodynamic parameters, phase partition and create thermo state.
    q = TD.PhasePartition(qᵥ, FT(0), FT(0))  # use dry air density to compute ql
    ts = TD.PhaseNonEquil_pTq(tps, p_air, T, q)
    ρ_air = TD.air_density(tps, ts)
    
    qᵢ = FT(0.0)
    qₗ = CPD.get_standard_N_q(clinfo.pdists).M_liq * ρₗ / ρ_air
    q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
    ts = TD.PhaseNonEquil_pTq(tps, p_air, T, q)

    # Constants and variables that depend on the moisture content
    R_air = TD.gas_constant_air(tps, q)
    cp_air = TD.cp_m(tps, q)
    L_subl = TD.latent_heat_sublim(tps, T)
    L_fus = TD.latent_heat_fusion(tps, T)
    L_vap = TD.latent_heat_vapor(tps, T)
    ρ_air = TD.air_density(tps, ts)

    # Adiabatic parcel coefficients
    a1 = L_vap * grav / cp_air / T^2 / Rᵥ - grav / R_air / T
    a2 = 1 / qᵥ
    a3 = L_vap^2 / Rᵥ / T^2 / cp_air
    a4 = L_vap * L_subl / Rᵥ / T^2 / cp_air
    a5 = L_vap * L_fus / Rᵥ / cp_air / (T^2)

    # Cloudy condnsational growth / evaporation
    dmom_ce = CMF.condensation(clinfo, aps, tps, T, Sₗ)
    @show dmom_ce

    # ... water mass budget
    dqₗ_dt_v2l = dmom_ce[2] + dmom_ce[4]

    # Update the tendecies
    dqₗ_dt = dqₗ_dt_v2l
    dqᵥ_dt = -dqₗ_dt

    dSₗ_dt =
        a1 * w * Sₗ - (a2 + a3) * Sₗ * dqₗ_dt_v2l

    dp_air_dt = -p_air * grav / R_air / T * w

    dT_dt =
        -grav / cp_air * w +
        L_vap / cp_air * dqₗ_dt_v2l

    # Set tendencies
    dY[1] = dSₗ_dt      # saturation ratio over liquid water
    dY[2] = dp_air_dt   # pressure
    dY[3] = dT_dt       # temperature
    dY[4] = dqᵥ_dt      # vapor specific humidity
    dY[5:end] = dmom_ce
end

function run_parcel_cloudy(Yinit, clinfo, t_0, t_end, pp)
    FT = typeof(t_0)

    println(" ")
    println("Size distributions: ", clinfo.pdists)
    println("Condensation growth only ")

    # Parameters for the ODE solver
    p = (
        wps = pp.wps, 
        aps = pp.aps,
        tps = pp.tps, 
        w = pp.w, 
        clinfo = clinfo
    )

    problem = ODE.ODEProblem(parcel_model_cloudy, Yinit, (FT(t_0), FT(t_end)), p)
    return ODE.solve(
        problem,
        ODE.Euler(),
        dt = pp.const_dt,
        reltol = 100 * eps(FT),
        abstol = 100 * eps(FT),
    )
end

# Get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)
aps = CMP.AirProperties(FT)
# Constants
ρₗ = wps.ρw
ρᵢ = wps.ρi
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)

# Initial conditions
dist_init = [
    CPD.ExponentialPrimitiveParticleDistribution(FT(100 * 1e6), FT(1e5*1e-18*1e3)), # 100/cm^3; 10^5 µm^3
    CPD.GammaPrimitiveParticleDistribution(FT(1 * 1e6), FT(1e6*1e-18*1e3), FT(1)),   # 1/cm^3; 10^6 µm^3; k=1
]
moments_init = FT.([100.0 * 1e6, 1e-2, 1.0 * 1e6, 1e-3, 2e-12])
p₀ = FT(800 * 1e2)
T₀ = FT(273.15 + 7.0)
e = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
Sₗ = FT(1)
md_v = (p₀ - e) / R_d / T₀
mv_v = e / R_v / T₀
ml_v = moments_init[2] + moments_init[4]
qᵥ = mv_v / (md_v + mv_v + ml_v)
qₗ = ml_v / (md_v + mv_v + ml_v)
qᵢ = FT(0)

# Simulation parameters passed into ODE solver
w = FT(10)                                  # updraft speed
const_dt = FT(0.5)                         # model timestep
t_max = FT(20)
clinfo = CMF.CLSetup{FT}(
    pdists = dist_init,
    mom = moments_init
)

Y0 = [Sₗ, p₀, T₀, qᵥ, moments_init...]
dY = zeros(FT,9)
p = (
    wps = wps,
    aps = aps,
    tps = tps, 
    w = w, 
    clinfo = clinfo, 
    const_dt = const_dt
)

# Setup the plots
fig = MK.Figure(size = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Supersaturation [%]")
ax2 = MK.Axis(fig[3, 1], xlabel = "Time [s]", ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "q_vap [g/kg]")
ax4 = MK.Axis(fig[2, 2], xlabel = "Time [s]", ylabel = "q_liq [g/kg]")
ax5 = MK.Axis(fig[1, 2], ylabel = "radius [μm]")

# solve ODE
sol = run_parcel_cloudy(Y0, clinfo, FT(0), t_max, p)
@show sol.u

# plot results - time series
MK.lines!(ax1, sol.t, (sol[1, :] .- 1) * 100.0)
MK.lines!(ax2, sol.t, sol[3, :])
MK.lines!(ax3, sol.t, sol[4, :] * 1e3)
ρ_air = sol[2, :] / R_v ./ sol[3, :]
M_l = sol[6, :] + sol[8, :]  # kg / m^3 air
N_l = sol[5, :] + sol[7, :]  # number / m^3 air
r_l = (M_l ./ N_l / ρₗ / 4 / π * 3).^(1/3) * 1e6 
MK.lines!(ax4, sol.t, M_l ./ ρ_air * 1e3)
MK.lines!(ax5, sol.t, r_l)

MK.save("cloudy_parcel.svg", fig)