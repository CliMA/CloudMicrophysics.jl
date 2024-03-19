import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions
import Random
import Distributions
import CairoMakie as MK
import LinearAlgebra
import OrdinaryDiffEq as ODE

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import Thermodynamics as TD

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float64

# Random number generator
rng_seed = 24
rng = Random.seed!(Random.GLOBAL_RNG, rng_seed)

# Get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)

# Define other parameters and initial conditions
ρₗ = wps.ρw
Nₐ = FT(8e6) # according to fig 2 panel b deposition nucleation experiment
Nₗ = FT(0)
Nᵢ = FT(0)
r₀ = FT(0.5e-6)
p₀ = FT(987.018 * 1e2)
T₀ = FT(212.978)

R_d = TD.Parameters.R_d(tps)
R_v = TD.Parameters.R_v(tps)
ϵₘ = R_d / R_v
C_v = FT(1.08509 * 1e-6)
qᵥ = ϵₘ / (ϵₘ - 1 + 1 / C_v)
qₗ = FT(0)
qᵢ = FT(0)
x_sulph = FT(0)

q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
Rₐ = TD.gas_constant_air(tps, q)
eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
e = eᵥ(qᵥ, p₀, Rₐ, R_v)
Sₗ = FT(e / eₛ)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]

w = FT(0.4)
const_dt = FT(1)
t_max = FT(500)
aerosol = CMP.Kaolinite(FT)
dep_nucleation = "ABDINM"
deposition_growth = "Deposition"
size_distribution = "Monodisperse"

params = (;
    const_dt,
    w,
    dep_nucleation,
    deposition_growth,
    size_distribution,
)

# Define model which wraps around parcel and overwrites calibrated parameters
function run_model(p, coefficients)
    # grabbing parameters
    ABDINM_m_calibrated, ABDINM_c_calibrated = coefficients
    (; const_dt, w) = p
    (; dep_nucleation, deposition_growth) = p
    (; size_distribution) = p

    # overwrite coefficients
    override_file = Dict(
        "China2017_J_deposition_m_Kaolinite" =>
            Dict("value" => ABDINM_m_calibrated, "type" => "float"),
        "China2017_J_deposition_c_Kaolinite" =>
            Dict("value" => ABDINM_c_calibrated, "type" => "float"),
    )
    kaolinite_calibrated = CP.create_toml_dict(FT; override_file)
    overwrite = CMP.Kaolinite(kaolinite_calibrated)

    # run parcel with new coefficients
    local params = parcel_params{FT}(
        const_dt = const_dt,
        w = w,
        aerosol = overwrite,
        deposition = dep_nucleation ,
        deposition_growth = deposition_growth,
        size_distribution = size_distribution,
    )

    # solve ODE
    local sol = run_parcel(IC, FT(0), t_max, params)
    @info(sol[9,:])         # TODO - figure out what to multiply the return by
    return sol[9, :] * 1e10 # ICNC
end

# Creating noisy pseudo-observations
observation_data_names = ["ABDINM_m", "ABDINM_c"]
coeff_true = [FT(27.551), FT(-2.2209)]
n_samples = 10
G_truth = run_model(params, coeff_true)     # ICNC from running parcel w default values
y_truth = zeros(length(G_truth), n_samples) # where noisy ICNC will be stored

dim_output = length(G_truth)
Γ = 0.01 * LinearAlgebra.I * (maximum(G_truth) - minimum(G_truth))
noise_dist = Distributions.MvNormal(zeros(dim_output), Γ)
y_truth = G_truth .+ rand(noise_dist)

# Define prior distributions for the coefficients
# stats = [mean, std dev, lower bound, upper bound]
ABDINM_m_stats = [FT(20), FT(1), FT(0), Inf]
ABDINM_c_stats = [FT(-1), FT(1), -Inf, Inf]
ABDINM_m_prior = EKP.constrained_gaussian(
    observation_data_names[1],
    ABDINM_m_stats[1],
    ABDINM_m_stats[2],
    ABDINM_m_stats[3],
    ABDINM_m_stats[4],
)
ABDINM_c_prior = EKP.constrained_gaussian(
    observation_data_names[2],
    ABDINM_c_stats[1],
    ABDINM_c_stats[2],
    ABDINM_c_stats[3],
    ABDINM_c_stats[4],
)
prior = EKP.combine_distributions([ABDINM_m_prior, ABDINM_c_prior])

# Generate initial ensember and set up EKI
N_ensemble = 20      # runs N_ensemble trials per iteration
N_iterations = 10    # number of iterations the inverse problem goes through
initial_ensemble = EKP.construct_initial_ensemble(rng, prior, N_ensemble)
EKI_obj = EKP.EnsembleKalmanProcess(
    initial_ensemble,
    y_truth,
    Γ,
    EKP.Inversion();
    rng = rng,
)

# Carry out the EKI calibration
# ϕ_n_values[iteration] stores ensembles of calibrated coeffs in that iteration
ϕ_n_values = []
for n in 1:N_iterations
    ϕ_n = EKP.get_ϕ_final(prior, EKI_obj)
    G_ens = hcat([run_model(params, ϕ_n[:, i]) for i in 1:N_ensemble]...)
    EKP.update_ensemble!(EKI_obj, G_ens)

    global ϕ_n_values = vcat(ϕ_n_values, [ϕ_n])
end

# Plotting ICNC w true params, the initial ensemble, and the final ensemble
fig = MK.Figure(size = (800, 600))
ax1 = MK.Axis(
    fig[1, 1],
    ylabel = "Deposition m [-]",
    xlabel = "iteration number",
)
ax2 = MK.Axis(
    fig[2, 1],
    ylabel = "Deposition c [-]",
    xlabel = "iteration number",
)

iterations = collect(1:N_iterations)
ABDINM_m_mean = zeros(length(iterations))
ABDINM_c_mean = zeros(length(iterations))
for iter in iterations
    ABDINM_m_mean[iter] =
        Distributions.mean(ϕ_n_values[iter][1, i] for i in 1:N_ensemble)
    ABDINM_c_mean[iter] =
        Distributions.mean(ϕ_n_values[iter][2, i] for i in 1:N_ensemble)
end

MK.lines!(
    ax1,
    iterations,
    ABDINM_m_mean,
    label = "ensemble mean",
    color = :orange,
)
MK.lines!(
    ax1,
    iterations,
    zeros(length(ABDINM_m_mean)) .+ coeff_true[1],
    label = "default value",
)
MK.lines!(
    ax2,
    iterations,
    ABDINM_c_mean,
    label = "ensemble mean",
    color = :orange,
)
MK.lines!(
    ax2,
    iterations,
    zeros(length(ABDINM_m_mean)) .+ coeff_true[2],
    label = "default value",
)

MK.axislegend(ax1, framevisible = true, labelsize = 12)
MK.axislegend(ax2, framevisible = true, labelsize = 12)

MK.save("perfect_calibration_dep.svg", fig)

# Check if the calibrated coefficients are similar to the known default values
ABDINM_m_ekp = round(
    Distributions.mean(ϕ_n_values[N_iterations][1, 1:N_ensemble]),
    digits = 6,
)
ABDINM_c_ekp = round(
    Distributions.mean(ϕ_n_values[N_iterations][2, 1:N_ensemble]),
    digits = 6,
)

println("Deposition m [-]: ", ABDINM_m_ekp, " vs ", coeff_true[1])
println("Deposition c [-]: ", ABDINM_c_ekp, " vs ", coeff_true[2])
