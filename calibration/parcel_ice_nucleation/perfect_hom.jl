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
Nₐ = FT(0)
Nₗ = FT(3.6 * 1e8)
Nᵢ = FT(0)
r₀ = FT(200 * 1e-9)
p₀ = FT(997.345 * 1e2)
T₀ = FT(243.134)

R_d = TD.Parameters.R_d(tps)
R_v = TD.Parameters.R_v(tps)
ϵₘ = R_d / R_v
qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2)) # 1.2 should be ρₐ
C_l = FT(qₗ / ((1 - qₗ) * ϵₘ + qₗ))  # concentration/mol fraction of liquid
C_v = FT(357.096 * 1e-6 - C_l)     # concentration/mol fraction of vapor
qᵥ = ϵₘ / (ϵₘ - 1 + 1 / C_v)
qᵢ = FT(0)
x_sulph = FT(0)

q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
Rₐ = TD.gas_constant_air(tps, q)
eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
e = eᵥ(qᵥ, p₀, Rₐ, R_v)
Sₗ = FT(e / eₛ)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]

w = FT(1.2)
const_dt = FT(1)
t_max = FT(650)
homogeneous = "ABHOM"
condensation_growth = "Condensation"
deposition_growth = "Deposition"
size_distribution = "Monodisperse"

params = (;
    const_dt,
    w,
    homogeneous,
    condensation_growth,
    deposition_growth,
    size_distribution,
)

# Define model which wraps around parcel and overwrites calibrated parameters
function run_model(p, coefficients)
    # grabbing parameters
    coeff1, coeff2, coeff3, coeff4 = coefficients
    (; const_dt, w) = p
    (; homogeneous, condensation_growth, deposition_growth) = p
    (; size_distribution) = p

    # overwrite coefficients
    override_file = Dict(
        "Koop2000_J_hom_coeff1" =>
            Dict("value" => coeff1, "type" => "float"),
        "Koop2000_J_hom_coeff2" =>
            Dict("value" => coeff2, "type" => "float"),
        "Koop2000_J_hom_coeff3" =>
            Dict("value" => coeff3, "type" => "float"),
        "Koop2000_J_hom_coeff4" =>
            Dict("value" => coeff4, "type" => "float"),
    )
    koop_coeffs_calibrated = CP.create_toml_dict(FT; override_file)
    overwrite = CMP.IceNucleationParameters(koop_coeffs_calibrated)

    # run parcel with new coefficients
    local params = parcel_params{FT}(
        const_dt = const_dt,
        w = w,
        homogeneous = homogeneous,
        condensation_growth = condensation_growth,
        deposition_growth = deposition_growth,
        size_distribution = size_distribution,
        ips = overwrite,
    )

    # solve ODE
    local sol = run_parcel(IC, FT(0), t_max, params)
    return sol[9, :] # ICNC
end

# Creating noisy pseudo-observations
observation_data_names = ["coeff1", "coeff2", "coeff3", "coeff4"]
coeff_true = [FT(-906.7), FT(8502), FT(26924), FT(29180)]
n_samples = 10
G_truth = run_model(params, coeff_true)     # ICNC from running parcel w default values
y_truth = zeros(length(G_truth), n_samples) # where noisy ICNC will be stored

dim_output = length(G_truth)
Γ = 0.01 * LinearAlgebra.I * (maximum(G_truth) - minimum(G_truth))
noise_dist = Distributions.MvNormal(zeros(dim_output), Γ)
y_truth = G_truth .+ rand(noise_dist)

# Define prior distributions for the coefficients
# stats = [mean, std dev, lower bound, upper bound]
coeff1_stats = [FT(-900), FT(1), -Inf, Inf]
coeff2_stats = [FT(8000), FT(1), -Inf, Inf]
coeff3_stats = [FT(26000), FT(1), -Inf, Inf]
coeff4_stats = [FT(29000), FT(1), -Inf, Inf]
prior_coeff1 = EKP.constrained_gaussian(
    observation_data_names[1],
    coeff1_stats[1],
    coeff1_stats[2],
    coeff1_stats[3],
    coeff1_stats[4],
)
prior_coeff2 = EKP.constrained_gaussian(
    observation_data_names[2],
    coeff2_stats[1],
    coeff2_stats[2],
    coeff2_stats[3],
    coeff2_stats[4],
)
prior_coeff3 = EKP.constrained_gaussian(
    observation_data_names[3],
    coeff3_stats[1],
    coeff3_stats[2],
    coeff3_stats[3],
    coeff3_stats[4],
)
prior_coeff4 = EKP.constrained_gaussian(
    observation_data_names[4],
    coeff4_stats[1],
    coeff4_stats[2],
    coeff4_stats[3],
    coeff4_stats[4],
)
prior = EKP.combine_distributions([prior_coeff1, prior_coeff2, prior_coeff3, prior_coeff4])

# Generate initial ensember and set up EKI
N_ensemble = 10     # runs N_ensemble trials per iteration
N_iterations = 10   # number of iterations the inverse problem goes through
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
ax1 = MK.(
    fig[1, 1],
    ylabel = "Koop2000 Coeff1 [-]",
    xlabel = "iteration number",
)
ax2 = MK.Axis(
    fig[1, 2],
    ylabel = "Koop2000 Coeff2 [-]",
    xlabel = "iteration number",
)
ax3 = MK.Axis(
    fig[2, 1],
    ylabel = "Koop2000 Coeff3 [-]",
    xlabel = "iteration number",
)
ax4 = MK.Axis(
    fig[2, 2],
    ylabel = "Koop2000 Coeff4 [-]",
    xlabel = "iteration number",
)

iterations = collect(1:N_iterations)
mean_coeff1 = zeros(length(iterations))
mean_coeff2 = zeros(length(iterations))
mean_coeff3 = zeros(length(iterations))
mean_coeff4 = zeros(length(iterations))

for iter in iterations
    mean_coeff1[iter] =
        Distributions.mean(ϕ_n_values[iter][1, i] for i in 1:N_ensemble)
    mean_coeff2[iter] =
        Distributions.mean(ϕ_n_values[iter][2, j] for j in 1:N_ensemble)
    mean_coeff3[iter] =
        Distributions.mean(ϕ_n_values[iter][3, k] for k in 1:N_ensemble)
    mean_coeff4[iter] =
        Distributions.mean(ϕ_n_values[iter][4, l] for l in 1:N_ensemble)
end

MK.lines!(
    ax1,
    iterations,
    mean_coeff1,
    label = "ensemble mean",
    color = :orange,
)
MK.lines!(
    ax1,
    iterations,
    zeros(length(mean_coeff1)) .+ coeff_true[1],
    label = "default value",
)
MK.lines!(
    ax2,
    iterations,
    mean_coeff2,
    label = "ensemble mean",
    color = :orange,
)
MK.lines!(
    ax2,
    iterations,
    zeros(length(mean_coeff2)) .+ coeff_true[2],
    label = "default value",
)
MK.lines!(
    ax3,
    iterations,
    mean_coeff3,
    label = "ensemble mean",
    color = :orange,
)
MK.lines!(
    ax3,
    iterations,
    zeros(length(mean_coeff3)) .+ coeff_true[3],
    label = "default value",
)
MK.lines!(
    ax4,
    iterations,
    mean_coeff4,
    label = "ensemble mean",
    color = :orange,
)
MK.lines!(
    ax4,
    iterations,
    zeros(length(mean_coeff4)) .+ coeff_true[4],
    label = "default value",
)

MK.axislegend(ax1, framevisible = true, labelsize = 12)
MK.axislegend(ax2, framevisible = true, labelsize = 12)
MK.axislegend(ax3, framevisible = true, labelsize = 12)
MK.axislegend(ax4, framevisible = true, labelsize = 12)

MK.save("calibration_default_IN.svg", fig)

# Check if the calibrated coefficients are similar to the known default values
coeff1_ekp = round(
    Distributions.mean(ϕ_n_values[N_iterations][1, 1:N_ensemble]),
    digits = 6,
)
coeff2_ekp = round(
    Distributions.mean(ϕ_n_values[N_iterations][2, 1:N_ensemble]),
    digits = 6,
)
coeff3_ekp = round(
    Distributions.mean(ϕ_n_values[N_iterations][3, 1:N_ensemble]),
    digits = 6,
)
coeff4_ekp = round(
    Distributions.mean(ϕ_n_values[N_iterations][4, 1:N_ensemble]),
    digits = 6,
)

println("coefficient 1 [-]: ", coeff1_ekp, " vs ", coeff_true[1])
println("coefficient 2 [-]: ", coeff2_ekp, " vs ", coeff_true[2])
println("coefficient 3 [-]: ", coeff3_ekp, " vs ", coeff_true[3])
println("coefficient 4 [-]: ", coeff4_ekp, " vs ", coeff_true[4])
