import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import Thermodynamics as TD

#! format: off
function perf_model_params(FT, IN_mode)
    if IN_mode == "ABDINM"
        const_dt = FT(1)
        w = FT(0.4)
        dep_nucleation = "ABDINM"
        heterogeneous = "None"
        homogeneous = "None"
        condensation_growth = "None"
        deposition_growth = "Deposition"
        size_distribution = "Monodisperse"
    elseif IN_mode == "ABIFM"
        const_dt = FT(1)
        w = FT(0.4)
        dep_nucleation = "None"
        heterogeneous = "ABIFM"
        homogeneous = "None"
        condensation_growth = "Condensation"
        deposition_growth = "Deposition"
        size_distribution = "Monodisperse"
    elseif IN_mode == "ABHOM"
        const_dt = FT(1)
        w = FT(1)
        dep_nucleation = "None"
        heterogeneous = "None"
        homogeneous = "ABHOM"
        condensation_growth = "None"
        deposition_growth = "Deposition"
        size_distribution = "Monodisperse"
    end
    params = (; const_dt, w,
        condensation_growth, deposition_growth,     # growth
        size_distribution,                          # size distribution
        dep_nucleation, heterogeneous, homogeneous, # nucleation
    )
    return params
end

function perf_model_IC(FT, IN_mode)
    tps = TD.Parameters.ThermodynamicsParameters(FT)
    wps = CMP.WaterProperties(FT)

    ρₗ = wps.ρw
    R_d = TD.Parameters.R_d(tps)
    R_v = TD.Parameters.R_v(tps)
    ϵₘ = R_d / R_v

    if IN_mode == "ABDINM"
        Nₐ = FT(8e6) # according to fig 2 panel b deposition nucleation experiment
        Nₗ = FT(0)
        Nᵢ = FT(0)
        r₀ = FT(0.5e-6)
        p₀ = FT(987.018 * 1e2)
        T₀ = FT(212.978)
        C_v = FT(10.8509 * 1e-6)
        qᵥ = ϵₘ / (ϵₘ - 1 + 1 / C_v)
        qₗ = FT(0)
        qᵢ = FT(0)
        q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
        Rₐ = TD.gas_constant_air(tps, q)
        eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
        e = eᵥ(qᵥ, p₀, Rₐ, R_v)
        Sₗ = FT(e / eₛ)
    elseif IN_mode == "ABIFM"
        Nₐ = FT(0)
        Nₗ = FT(8e6)
        Nᵢ = FT(0)
        r₀ = FT(0.5e-6)
        p₀ = FT(987.018 * 1e2)
        T₀ = FT(212.978)
        qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2)) # 1.2 should be ρₐ
        C_l = FT(qₗ / ((1 - qₗ) * ϵₘ + qₗ))  # concentration/mol fraction of liquid
        C_v = FT(10.8509 * 1e-6 - C_l)     # concentration/mol fraction of vapor
        qᵥ = ϵₘ / (ϵₘ - 1 + 1 / C_v)
        qᵢ = FT(0)
        q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
        Rₐ = TD.gas_constant_air(tps, q)
        eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
        e = eᵥ(qᵥ, p₀, Rₐ, R_v)
        Sₗ = FT(e / eₛ)
    elseif IN_mode == "ABHOM"
        Nₐ = FT(0)
        Nₗ = FT(300 * 1e6)
        Nᵢ = FT(0)
        r₀ = FT(25 * 1e-9)
        p₀ = FT(9712.183)
        T₀ = FT(190)
        qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2)) # 1.2 should be ρₐ
        C_l = FT(qₗ / ((1 - qₗ) * ϵₘ + qₗ))  # concentration/mol fraction of liquid
        C_v = FT(5 * 1e-6)
        qᵥ = ϵₘ / (ϵₘ - 1 + 1 / C_v)
        qᵢ = FT(0)
        q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
        Rₐ = TD.gas_constant_air(tps, q)
        eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
        e = eᵥ(qᵥ, p₀, Rₐ, R_v)
        Sₗ = FT(e / eₛ)
    end
    return [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, FT(0)]
end

function perf_model_pseudo_data(FT, IN_mode, params, IC)
    n_samples = 10

    if IN_mode == "ABDINM"
        coeff_true = [FT(27.551), FT(-2.2209)]
    elseif IN_mode == "ABIFM"
        coeff_true = [FT(54.58834), FT(-10.54758)]
    elseif IN_mode == "ABHOM"
        coeff_true = [FT(255.927125), FT(-68.553283)]
    end

    G_truth = run_model(params, coeff_true, IN_mode, FT, IC)
    dim_output = length(G_truth)

    Γ = 0.005 * LinearAlgebra.I * (maximum(G_truth) - minimum(G_truth))
    noise_dist = Distributions.MvNormal(zeros(dim_output), Γ)

    y_truth = zeros(length(G_truth), n_samples) # where noisy ICNC will be stored
    y_truth = G_truth .+ rand(noise_dist)
    return [y_truth, Γ, coeff_true]
end
#! format: on
