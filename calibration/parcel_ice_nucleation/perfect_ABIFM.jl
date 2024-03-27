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

# Define other parameters and initial conditions
ρₗ = wps.ρw
Nₐ = FT(0)
Nₗ = FT(2000)
Nᵢ = FT(0)
r₀ = FT(1e-6)
p₀ = FT(800 * 1e2)
T₀ = FT(251)
qᵥ = FT(8.1e-4)
qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2)) # 1.2 should be ρₐ
qᵢ = FT(0)
x_sulph = FT(0)

q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
R_v = TD.Parameters.R_v(tps)
Rₐ = TD.gas_constant_air(tps, q)
eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
e = eᵥ(qᵥ, p₀, Rₐ, R_v)
Sₗ = FT(e / eₛ)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]

w = FT(0.4)
const_dt = FT(1)
t_max = FT(500)
aerosol = CMP.DesertDust(FT)
heterogeneous = "ABIFM"
condensation_growth = "Condensation"
deposition_growth = "Deposition"
size_distribution = "Monodisperse"

params = (;
    const_dt,
    w,
    heterogeneous,
    condensation_growth,
    deposition_growth,
    size_distribution,
)

# Define model which wraps around parcel and overwrites calibrated parameters
function run_model(p, coefficients)
    # grabbing parameters
    m_coeff_calibrated, c_coeff_calibrated = coefficients
    (; const_dt, w) = p
    (; heterogeneous, condensation_growth, deposition_growth) = p
    (; size_distribution) = p

    # overwrite coefficients
    override_file = Dict(
        "AlpertKnopf2016_J_ABIFM_m_DesertDust" =>
            Dict("value" => m_coeff_calibrated, "type" => "float"),
        "AlpertKnopf2016_J_ABIFM_c_DesertDust" =>
            Dict("value" => c_coeff_calibrated, "type" => "float"),
    )
    desertdust_calibrated = CP.create_toml_dict(FT; override_file)
    overwrite = CMP.DesertDust(desertdust_calibrated)

    # run parcel with new coefficients
    local params = parcel_params{FT}(
        const_dt = const_dt,
        w = w,
        aerosol = overwrite,
        heterogeneous = heterogeneous,
        condensation_growth = condensation_growth,
        deposition_growth = deposition_growth,
        size_distribution = size_distribution,
    )

    # solve ODE
    local sol = run_parcel(IC, FT(0), t_max, params)
    return sol[9, :] # ICNC
end

# Creating noisy pseudo-observations
observation_data_names = ["m_coeff", "c_coeff"]
coeff_true = [FT(22.62), FT(-1.35)]         # [aerosol.ABIFM_m, aerosol.ABIFM_c] ; should this be hard coded?
n_samples = 10                              # what is this?
G_truth = run_model(params, coeff_true)     # ICNC from running parcel w default values
y_truth = zeros(length(G_truth), n_samples) # where noisy ICNC will be stored

dim_output = length(G_truth)
Γ = 0.01 * LinearAlgebra.I * (maximum(G_truth) - minimum(G_truth))
noise_dist = Distributions.MvNormal(zeros(dim_output), Γ)
y_truth = G_truth .+ rand(noise_dist)

# Define prior distributions for the coefficients
# stats = [mean, std dev, lower bound, upper bound]
m_stats = [FT(20), FT(1), FT(0), Inf]
c_stats = [FT(-1), FT(1), -Inf, Inf]
prior_m = EKP.constrained_gaussian(
    observation_data_names[1],
    m_stats[1],
    m_stats[2],
    m_stats[3],
    m_stats[4],
)
prior_c = EKP.constrained_gaussian(
    observation_data_names[2],
    c_stats[1],
    c_stats[2],
    c_stats[3],
    c_stats[4],
)
prior = EKP.combine_distributions([prior_m, prior_c])

# Generate initial ensember and set up EKI
N_ensemble = 20     # runs N_ensemble trials per iteration
N_iterations = 5    # number of iterations the inverse problem goes through
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
    ylabel = "m coefficient [-]",
    xlabel = "iteration number",
)
ax2 = MK.Axis(
    fig[2, 1],
    ylabel = "c coefficient [-]",
    xlabel = "iteration number",
)

iterations = collect(1:N_iterations)
mean_m_coeff = zeros(length(iterations))
mean_c_coeff = zeros(length(iterations))
for iter in iterations
    mean_m_coeff[iter] =
        Distributions.mean(ϕ_n_values[iter][1, i] for i in 1:N_ensemble)
    mean_c_coeff[iter] =
        Distributions.mean(ϕ_n_values[iter][2, j] for j in 1:N_ensemble)
end

MK.lines!(
    ax1,
    iterations,
    mean_m_coeff,
    label = "ensemble mean",
    color = :orange,
)
MK.lines!(
    ax1,
    iterations,
    zeros(length(mean_m_coeff)) .+ coeff_true[1],
    label = "default value",
)
MK.lines!(
    ax2,
    iterations,
    mean_c_coeff,
    label = "ensemble mean",
    color = :orange,
)
MK.lines!(
    ax2,
    iterations,
    zeros(length(mean_c_coeff)) .+ coeff_true[2],
    label = "default value",
)

MK.axislegend(ax1, framevisible = true, labelsize = 12)
MK.axislegend(ax2, framevisible = true, labelsize = 12)

MK.save("calibration_default_IN.svg", fig)

# Check if the calibrated coefficients are similar to the known default values
m_coeff_ekp = round(
    Distributions.mean(ϕ_n_values[N_iterations][1, 1:N_ensemble]),
    digits = 6,
)
c_coeff_ekp = round(
    Distributions.mean(ϕ_n_values[N_iterations][2, 1:N_ensemble]),
    digits = 6,
)

println("m coefficient [-]: ", m_coeff_ekp, " vs ", coeff_true[1])
println("c coefficient [-]: ", c_coeff_ekp, " vs ", coeff_true[2])
