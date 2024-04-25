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

    # ... water mass budget
    dqₗ_dt_v2l = dmom_ce[2] / ρ_air #+ dmom_ce[4]

    # Update the tendecies
    dqₗ_dt = dqₗ_dt_v2l
    dqᵥ_dt = -dqₗ_dt

    dSₗ_dt = a1 * w * Sₗ - (a2 + a3) * Sₗ * dqₗ_dt_v2l

    dp_air_dt = -p_air * grav / R_air / T * w

    dT_dt = -grav / cp_air * w + L_vap / cp_air * dqₗ_dt_v2l

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
    p = (wps = pp.wps, aps = pp.aps, tps = pp.tps, w = pp.w, clinfo = clinfo)

    problem =
        ODE.ODEProblem(parcel_model_cloudy, Yinit, (FT(t_0), FT(t_end)), p)
    return ODE.solve(
        problem,
        ODE.Euler(),
        dt = pp.const_dt,
        reltol = 100 * eps(FT),
        abstol = 100 * eps(FT),
    )
end

function init_conditions(ρₗ, type::String)
    FT = typeof(ρₗ)
    r₀ = 8e-6
    N = 200 * 1e6
    m₀ = FT(4 / 3 * FT(π) * r₀^3 * ρₗ)
    if type == "monodisperse"
        dist_init = [
            CPD.MonodispersePrimitiveParticleDistribution(
                N,
                m₀
            )
        ]
        moments_init = [N, N*m₀]
        ml_v = moments_init[2]
    elseif type == "gamma"
        k = FT(2)
        θ = m₀ / k
        dist_init = [CPD.GammaPrimitiveParticleDistribution(
            N,
            θ,
            k
        )]
        moments_init = CPD.get_moments(dist_init[1])
        ml_v = moments_init[2]
    elseif type == "mixture"
        M0 = [9 * N / 10, N / 10]
        M1 = [N * m₀ / 2, N * m₀ / 2]
        k = FT(2)
        dist_init = [
            CPD.ExponentialPrimitiveParticleDistribution(
                M0[1],
                M1[1] / M0[1]
            )
            CPD.GammaPrimitiveParticleDistribution(
                M0[2],
                M1[2] / M0[2] / k,
                k
            )   
        ]
        moments_init = vcat(CPD.get_moments.(dist_init)...)
        ml_v = moments_init[2] + moments_init[4]
    end
    return (dist_init, moments_init, ml_v)
end

function get_spectrum(clinfo, moments)
    x = 10 .^ (collect(range(-13, -10, 100)))
    r = (x / 1000 * 3 / 4 / π) .^ (1 / 3) * 1e6 # plot in µm
    y = zeros(FT, length(r))

    for (i, dist) in enumerate(clinfo.pdists)
        ind_rng = CL.get_dist_moments_ind_range(clinfo.NProgMoms, i)
        CPD.update_dist_from_moments!(dist, moments[ind_rng])
        y = y .+ (3 * x .^ 2 .* dist.(x))
    end

    return (r, y)
end

##### RUN SIMULATION ###

# Get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)
aps = CMP.AirProperties(FT)
# Constants
ρₗ = wps.ρw
ρᵢ = wps.ρi
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)
# Atmospheric State
p₀ = FT(800 * 1e2)
T₀ = FT(273.15 + 7.0)
e = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
Sₗ = FT(1)
md_v = (p₀ - e) / R_d / T₀
mv_v = e / R_v / T₀
qᵢ = FT(0)
# Simulation parameters passed into ODE solver
w = FT(10)                                  # updraft speed
const_dt = FT(0.5)                         # model timestep
t_max = FT(20)

# Data from Rogers(1975) Figure 1
# https://www.tandfonline.com/doi/abs/10.1080/00046973.1975.9648397
#! format: off
Rogers_time_supersat = [0.0645, 0.511, 0.883, 1.4, 2.07, 2.72, 3.24, 3.89, 4.53, 5.87, 7.16, 9.79, 16.0, 19.8]
Rogers_supersat = [0.0268, 0.255, 0.393, 0.546, 0.707, 0.805, 0.863, 0.905, 0.938, 0.971, 0.978, 0.963, 0.910, 0.885]
Rogers_time_radius = [0.561, 2, 3.99, 10.7, 14.9, 19.9]
Rogers_radius = [8.0, 8.08, 8.26, 8.91, 9.26, 9.68]
#! format: on

# Setup the plots
fig = MK.Figure(size = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Supersaturation [%]")
ax2 = MK.Axis(fig[3, 1], xlabel = "Time [s]", ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "q_vap [g/kg]")
ax4 = MK.Axis(fig[2, 2], xlabel = "Time [s]", ylabel = "q_liq [g/kg]")
ax5 = MK.Axis(fig[1, 2], ylabel = "radius [μm]")
ax6 = MK.Axis(fig[3, 2], xlabel = "radius [μm]", ylabel = "dm / d(ln r) [kg / m^3]", xscale=log10)
MK.lines!(ax1, Rogers_time_supersat, Rogers_supersat, label = "Rogers_1975")
MK.lines!(ax5, Rogers_time_radius, Rogers_radius)

# Initial conditions
size_distribution_list = ["monodisperse", "gamma", "mixture"]
for j in 1:length(size_distribution_list)

    DSD = size_distribution_list[j]
    (dist_init, moments_init, ml_v) = init_conditions(ρₗ, DSD)

    qᵥ = mv_v / (md_v + mv_v + ml_v)
    qₗ = ml_v / (md_v + mv_v + ml_v)

    clinfo = CMF.CLSetup{FT}(pdists = dist_init, mom = moments_init)
    
    (r, y) = get_spectrum(clinfo, moments_init)
    if DSD != "monodisperse"
        MK.lines!(ax6, r, y, color=MK.Cycled(j+1), linestyle=:dash, label = "init")
    end


    Y0 = [Sₗ, p₀, T₀, qᵥ, moments_init...]
    dY = zeros(FT, 9)
    p = (
        wps = wps,
        aps = aps,
        tps = tps,
        w = w,
        clinfo = clinfo,
        const_dt = const_dt,
    )

    # solve ODE
    sol = run_parcel_cloudy(Y0, clinfo, FT(0), t_max, p)

    # plot results - time series
    MK.lines!(ax1, sol.t, (sol[1, :] .- 1) * 100.0, label = DSD)
    MK.lines!(ax2, sol.t, sol[3, :])
    MK.lines!(ax3, sol.t, sol[4, :] * 1e3)

    ρ_air = sol[2, :] / R_v ./ sol[3, :]
    M_l = sol[6, :] #+ sol[8, :]  # kg / m^3 air
    N_l = sol[5, :] #+ sol[7, :]  # number / m^3 air
    r_l = (M_l ./ N_l / ρₗ / 4 / π * 3) .^ (1 / 3) * 1e6
    (r, y) = get_spectrum(clinfo, sol[5:end, end])

    MK.lines!(ax4, sol.t, M_l ./ ρ_air * 1e3)
    MK.lines!(ax5, sol.t, r_l)
    if DSD != "monodisperse"
        MK.lines!(ax6, r, y, color=MK.Cycled(j+1), label = "final")
    end
end

MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 12,
    orientation = :horizontal,
    nbanks = 2,
    position = :rb,
)

MK.axislegend(
    ax6,
    framevisible = false,
    labelsize = 12,
    orientation = :horizontal,
    nbanks = 2,
    position = :rt,
)

MK.save("cloudy_parcel.svg", fig)
