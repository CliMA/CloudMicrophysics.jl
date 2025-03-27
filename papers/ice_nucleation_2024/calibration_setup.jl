import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import Thermodynamics as TD

#! format: off
function perf_model_params(FT, IN_mode)
    if IN_mode == "ABDINM"
        const_dt = FT(1)
        w = FT(5)
        t_max = FT(100)
        aerosol_act = "None"
        aerosol = CMP.Kaolinite(FT)
        aero_σ_g = FT(0)
        r_nuc = FT(1e-7)
        A_aer = FT(4 * π * r_nuc^2)
        dep_nucleation = "ABDINM"
        heterogeneous = "ABIFM"
        homogeneous = "ABHOM"
        condensation_growth = "None"
        deposition_growth = "Deposition"
        liq_size_distribution = "Monodisperse"
        ice_size_distribution = "Monodisperse"
        prescribed_thermodynamics = false
        t_profile = []
        T_profile = []
        P_profile = []
        ips = CMP.IceNucleationParameters(FT)
    elseif IN_mode == "ABIFM"
        const_dt = FT(1)
        w = FT(1)
        t_max = FT(100)
        aerosol_act = "None"
        aerosol = CMP.Kaolinite(FT)
        aero_σ_g = FT(0)
        r_nuc = FT(1e-7)
        A_aer = FT(4 * π * r_nuc^2)
        dep_nucleation = "ABDINM"
        heterogeneous = "ABIFM"
        homogeneous = "ABHOM"
        condensation_growth = "Condensation"
        deposition_growth = "Deposition"
        liq_size_distribution = "Monodisperse"
        ice_size_distribution = "Monodisperse"
        prescribed_thermodynamics = false
        t_profile = []
        T_profile = []
        P_profile = []
        ips = CMP.IceNucleationParameters(FT)
    elseif IN_mode == "ABHOM"
        const_dt = FT(1)
        w = FT(1)
        t_max = FT(100)
        aerosol_act = "None"
        aerosol = CMP.Sulfate(FT)
        aero_σ_g = FT(0)
        r_nuc = FT(1e-7)
        A_aer = FT(4 * π * r_nuc^2)
        dep_nucleation = "ABDINM"
        heterogeneous = "ABIFM"
        homogeneous = "ABHOM"
        condensation_growth = "None"
        deposition_growth = "Deposition"
        liq_size_distribution = "Monodisperse"
        ice_size_distribution = "Monodisperse"
        prescribed_thermodynamics = false
        t_profile = []
        T_profile = []
        P_profile = []
        ips = CMP.IceNucleationParameters(FT)
    end
    params = (; const_dt, w, t_max, ips,
        prescribed_thermodynamics, t_profile, T_profile, P_profile,
        aerosol_act, aerosol, r_nuc, aero_σ_g, A_aer,   # aerosol activation
        condensation_growth, deposition_growth,         # growth
        liq_size_distribution, ice_size_distribution,   # size distribution
        dep_nucleation, heterogeneous, homogeneous,     # nucleation
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
        Nₐ = FT(2e8)
        Nₗ = FT(0)
        Nᵢ = FT(0)
        r₀ = FT(1e-7)
        p₀ = FT(80000)
        T₀ = FT(230)
        qᵥ = FT(8.8e-5)
        qₗ = FT(0)
        qᵢ = FT(0)
        q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
        Rₐ = TD.gas_constant_air(tps, q)
        eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
        e = eᵥ(qᵥ, p₀, Rₐ, R_v)
        Sₗ = FT(e / eₛ)
    elseif IN_mode == "ABIFM"
        Nₐ = FT(0)
        Nₗ = FT(2000)
        Nᵢ = FT(0)
        r₀ = FT(1e-6)
        p₀ = FT(800 * 1e2)
        T₀ = FT(251)
        qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2)) # 1.2 should be ρₐ
        qᵥ = FT(8.1e-4)
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

function perf_model_pseudo_data(FT, IN_mode, params, IC, end_sim)
    n_samples = 10

    if IN_mode == "ABDINM"
        coeff_true = [FT(27.551), FT(-2.2209)]
    elseif IN_mode == "ABIFM"
        coeff_true = [FT(54.58834), FT(-10.54758)]
    elseif IN_mode == "ABHOM"
        coeff_true = [FT(255.927125), FT(-68.553283)]
    end

    G_truth = run_model(params, IN_mode, coeff_true, FT, IC, end_sim, calibration = true)
    dim_output = length(G_truth)

    Γ = 0.01 * LinearAlgebra.I * (maximum(G_truth) - minimum(G_truth))
    noise_dist = Distributions.MvNormal(zeros(dim_output), Γ)

    y_truth = zeros(length(G_truth), n_samples) # where noisy ICNC will be stored
    y_truth = G_truth .+ rand(noise_dist)
    return [y_truth, Γ, coeff_true]
end

function AIDA_IN05_params(FT, w, t_max, t_profile, T_profile, P_profile)
    IN_mode = "ABHOM"
    const_dt = FT(1)
    prescribed_thermodynamics = true
    aerosol_act = "AeroAct"
    aerosol = CMP.Sulfate(FT)
    dep_nucleation = "ABDINM"
    heterogeneous = "ABIFM"
    homogeneous = "ABHOM"
    condensation_growth = "Condensation"
    deposition_growth = "Deposition"
    liq_size_distribution = "Gamma"
    ice_size_distribution = "Gamma"
    aero_σ_g = FT(2.3)
    r_nuc = FT(1e-7) #FT(3.057e-6)
    A_aer = FT(4 * π * r_nuc^2)
    ips = CMP.IceNucleationParameters(FT)

    params = (; const_dt, w, t_max, ips,
        prescribed_thermodynamics, t_profile, T_profile, P_profile,
        aerosol_act, aerosol, r_nuc, aero_σ_g, A_aer,   # aerosol activation
        condensation_growth, deposition_growth,         # growth
        liq_size_distribution, ice_size_distribution,   # size distribution
        dep_nucleation, heterogeneous, homogeneous,     # ice nucleation
    )
    return params
end

function AIDA_IN05_IC(FT, data_file)
    tps = TD.Parameters.ThermodynamicsParameters(FT)
    wps = CMP.WaterProperties(FT)

    ρₗ = wps.ρw
    ρᵢ = wps.ρi
    R_d = TD.Parameters.R_d(tps)
    R_v = TD.Parameters.R_v(tps)
    ϵₘ = R_d / R_v

    if data_file == "in05_17_aida.edf"
        Nₗ = FT(320.876 * 1e6)
        Nᵢ = FT(0.05 * 1e6)
        Nₐ = FT(360 * 1e6) - Nₗ - Nᵢ
        r₀ = FT(5.18056 / 2 * 1e-6)
        p₀ = FT(894.409 * 1e2)
        T₀ = FT(237.871)
        qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))  # 1.2 should be ρₐ
        qᵢ = FT(Nᵢ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))
        m_l = Nₗ * ρₗ *  4 * π / 3 * r₀^3
        m_i = Nᵢ * ρᵢ *  4 * π / 3 * r₀^3
        e = FT(30.2314)
        qᵥ = (e / R_v / T₀) / ((p₀ - e) / (R_d * T₀) + e / R_v / T₀ + m_l + m_i)
        q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
        Rₐ = TD.gas_constant_air(tps, q)
        eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
        Sₗ = FT(e / eₛ)
    elseif data_file == "in05_18_aida.edf"
        Nₗ = FT(212.885 * 1e6)
        Nᵢ = FT(0.787 * 1e6)
        Nₐ = FT(275 * 1e6) - Nₗ - Nᵢ
        r₀ = FT(6.83925 / 2 * 1e-6)
        p₀ = FT(878.345 * 1e2)
        T₀ = FT(237.234)
        qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2)) # 1.2 should be ρₐ
        qᵢ = FT(Nᵢ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))
        m_l = Nₗ * ρₗ *  4 * π / 3 * r₀^3
        m_i = Nᵢ * ρᵢ *  4 * π / 3 * r₀^3
        e = FT(28.9235)
        qᵥ = (e / R_v / T₀) / ((p₀ - e) / (R_d * T₀) + e / R_v / T₀ + m_l + m_i)
        q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
        Rₐ = TD.gas_constant_air(tps, q)
        eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
        Sₗ = FT(e / eₛ)
    end
    return [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, FT(0)]
end

function AIDA_IN07_params(FT, w, t_max, t_profile, T_profile, P_profile, batch_name)
    IN_mode = "ABDINM"
    const_dt = FT(1)
    prescribed_thermodynamics = true
    aerosol_act = "None"
    if batch_name == "IN0701"
        aerosol = CMP.ArizonaTestDust(FT)
    elseif batch_name == "IN0719"
        aerosol = CMP.Illite(FT)
    elseif batch_name == "DEP"
        aerosol = CMP.Dust(FT)
    end
    dep_nucleation = "ABDINM"
    heterogeneous = "ABIFM"
    homogeneous = "ABHOM"
    condensation_growth = "Condensation"
    deposition_growth = "Deposition"
    liq_size_distribution = "Gamma"
    ice_size_distribution = "Gamma"
    aero_σ_g = FT(2.3)
    r_nuc = FT(0.1742857 * 1e-6)
    A_aer = FT(4 * π * r_nuc^2)
    # r_nuc = r₀ in IC
    ips = CMP.IceNucleationParameters(FT)

    params = (; const_dt, w, t_max, ips,
        prescribed_thermodynamics, t_profile, T_profile, P_profile,
        aerosol_act, aerosol, r_nuc, aero_σ_g, A_aer,   # aerosol activation
        condensation_growth, deposition_growth,         # growth
        liq_size_distribution, ice_size_distribution,   # size distribution
        dep_nucleation, heterogeneous, homogeneous,     # ice nucleation
    )
    return params
end

function AIDA_IN07_IC(FT, data_file)
    tps = TD.Parameters.ThermodynamicsParameters(FT)
    wps = CMP.WaterProperties(FT)

    ρₗ = wps.ρw
    ρᵢ = wps.ρi
    R_d = TD.Parameters.R_d(tps)
    R_v = TD.Parameters.R_v(tps)
    ϵₘ = R_d / R_v
    if data_file == "in07_01_aida.edf"
        Nₗ = FT(14.2922 * 1e6)
        Nᵢ = FT(0.635 * 1e6)
        Nₐ = FT(146 * 1e6) - Nₗ - Nᵢ
        r₀ = FT(0.0146991 / 2 * 1e-6)
        # r₀ is weighted avg of two modes as listed in Table 2 of Mohler et al (2008)
        p₀ = FT(981.610 * 1e2)
        T₀ = FT(208.917)
        qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))  # 1.2 should be ρₐ
        qᵢ = FT(Nᵢ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))
        m_l = Nₗ * ρₗ *  4 * π / 3 * r₀^3
        m_i = Nᵢ * ρᵢ *  4 * π / 3 * r₀^3
        e = FT(0.682376)
        qᵥ = (e / R_v / T₀) / ((p₀ - e) / (R_d * T₀) + e / R_v / T₀ + m_l + m_i)
        q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
        Rₐ = TD.gas_constant_air(tps, q)
        eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
        Sₗ = FT(e / eₛ)
    elseif data_file == "in07_19_aida.edf"
        Nₗ = FT(0)
        Nᵢ = FT(0)
        Nₐ = FT(95 * 1e6) - Nₗ - Nᵢ
        r₀ = FT(0.1742857 * 1e-6)
        # r₀ is weighted avg of two modes as listed in Table 2 of Mohler et al (2008)
        p₀ = FT(982.977 * 1e2)
        T₀ = FT(208.961)
        qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2)) # 1.2 should be ρₐ
        qᵢ = FT(Nᵢ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))
        m_l = Nₗ * ρₗ *  4 * π / 3 * r₀^3
        m_i = Nᵢ * ρᵢ *  4 * π / 3 * r₀^3
        e = FT(0.632723)
        qᵥ = (e / R_v / T₀) / ((p₀ - e) / (R_d * T₀) + e / R_v / T₀ + m_l + m_i)
        q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
        Rₐ = TD.gas_constant_air(tps, q)
        eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
        Sₗ = FT(e / eₛ)
    end
    return [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, FT(0)]
end

function TROPIC04_params(FT, w, t_max, t_profile, T_profile, P_profile)
    IN_mode = "ABHOM"
    const_dt = FT(1)
    prescribed_thermodynamics = true
    aerosol_act = "AeroAct"
    aerosol = CMP.Sulfate(FT)
    dep_nucleation = "ABDINM"
    heterogeneous = "ABIFM"
    homogeneous = "ABHOM"
    condensation_growth = "Condensation"
    deposition_growth = "Deposition"
    liq_size_distribution = "Gamma"
    ice_size_distribution = "Gamma"
    aero_σ_g = FT(0.45605)
    r_nuc = FT(1.15 / 2 * 1e-6)
    A_aer = FT(4 * π * r_nuc^2)
    ips = CMP.IceNucleationParameters(FT)

    params = (; const_dt, w, t_max, ips,
        prescribed_thermodynamics, t_profile, T_profile, P_profile,
        aerosol_act, aerosol, r_nuc, aero_σ_g, A_aer,   # aerosol activation
        condensation_growth, deposition_growth,         # growth
        liq_size_distribution, ice_size_distribution,   # size distribution
        dep_nucleation, heterogeneous, homogeneous,     # ice nucleation
    )
    return params
end

function TROPIC04_IC(FT)
    # refers to exp 4 of campaign TROPIC04
    tps = TD.Parameters.ThermodynamicsParameters(FT)
    wps = CMP.WaterProperties(FT)

    ρₗ = wps.ρw
    ρᵢ = wps.ρi
    R_d = TD.Parameters.R_d(tps)
    R_v = TD.Parameters.R_v(tps)

    Nₗ = FT(0)
    Nᵢ = FT(0)
    Nₐ = FT(66.19233 * 1e6) - Nₗ - Nᵢ
    r₀ = FT(1.15 / 2 * 1e-6)
    p₀ = FT(100477.45358)
    T₀ = FT(209.79503)
    qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))  # 1.2 should be ρₐ
    qᵢ = FT(Nᵢ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))
    m_l = Nₗ * ρₗ *  4 * π / 3 * r₀^3
    m_i = Nᵢ * ρᵢ *  4 * π / 3 * r₀^3
    e = FT(0.6418042)
    qᵥ = (e / R_v / T₀) / ((p₀ - e) / (R_d * T₀) + e / R_v / T₀ + m_l + m_i)
    eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
    Sₗ = FT(e / eₛ)

    return [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, FT(0)]
end

function ACI04_22_params(FT, w, t_max, t_profile, T_profile, P_profile, batch_name)
    # Niemand et al (2012)
    IN_mode = "ABIFM"
    const_dt = FT(1)
    prescribed_thermodynamics = true
    aerosol_act = "None"
    aerosol = batch_name == "ACI04_22" ? CMP.MiddleEasternDust(FT) : CMP.Dust(FT)
    dep_nucleation = "ABDINM"
    heterogeneous = "ABIFM"
    homogeneous = "ABHOM"
    condensation_growth = "Condensation"
    deposition_growth = "Deposition"
    liq_size_distribution = "Gamma"
    ice_size_distribution = "Gamma"
    aero_σ_g = FT(1.47)           # avg of 2 modes
    r_nuc = FT(0.645 / 2 * 1e-6)  # avg of 2 modes
    A_aer = FT(4 * π * r_nuc^2)
    ips = CMP.IceNucleationParameters(FT)

    params = (; const_dt, w, t_max, ips,
        prescribed_thermodynamics, t_profile, T_profile, P_profile,
        aerosol_act, aerosol, r_nuc, aero_σ_g, A_aer,   # aerosol activation
        condensation_growth, deposition_growth,         # growth
        liq_size_distribution, ice_size_distribution,   # size distribution
        dep_nucleation, heterogeneous, homogeneous,     # ice nucleation
    )
    return params
end

function ACI04_22_IC(FT)
    # Niemand et al (2012)
    tps = TD.Parameters.ThermodynamicsParameters(FT)
    wps = CMP.WaterProperties(FT)

    ρₗ = wps.ρw
    ρᵢ = wps.ρi
    R_d = TD.Parameters.R_d(tps)
    R_v = TD.Parameters.R_v(tps)

    Nᵢ = FT(10075.141555)
    Nₗ = FT(252.1 * 1e6) - Nᵢ
    Nₐ = FT(0)
    r₀ = FT(10 * 1e-6)
    p₀ = FT(97785.714286)
    T₀ = FT(253.4983146)
    qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))  # 1.2 should be ρₐ
    qᵢ = FT(Nᵢ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))
    m_l = Nₗ * ρₗ *  4 * π / 3 * r₀^3
    m_i = Nᵢ * ρᵢ *  4 * π / 3 * r₀^3
    e = FT(104.5108)
    qᵥ = (e / R_v / T₀) / ((p₀ - e) / (R_d * T₀) + e / R_v / T₀ + m_l + m_i)
    eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
    Sₗ = FT(e / eₛ)

    return [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, FT(0)]
end

function EXP19_params(FT, w, t_max, t_profile, T_profile, P_profile, batch_name)
    # Cotten et al (2007)
    IN_mode = "ABIFM"
    const_dt = FT(1)
    prescribed_thermodynamics = true
    aerosol_act = "None"
    aerosol = batch_name == "EXP19" ? CMP.AsianDust(FT) : CMP.Dust(FT)
    dep_nucleation = "ABDINM"
    heterogeneous = "ABIFM"
    homogeneous = "ABHOM"
    condensation_growth = "Condensation"
    deposition_growth = "Deposition"
    liq_size_distribution = "Gamma"
    ice_size_distribution = "Gamma"
    aero_σ_g = FT(log(1.75))    # converted std dev to geometric std dev
    r_nuc = FT(0.4 / 2 * 1e-6)   # value is mode radius, not mean
    A_aer = FT(4 * π * r_nuc^2)
    ips = CMP.IceNucleationParameters(FT)

    params = (; const_dt, w, t_max, ips,
        prescribed_thermodynamics, t_profile, T_profile, P_profile,
        aerosol_act, aerosol, r_nuc, aero_σ_g, A_aer,   # aerosol activation
        condensation_growth, deposition_growth,         # growth
        liq_size_distribution, ice_size_distribution,   # size distribution
        dep_nucleation, heterogeneous, homogeneous,     # ice nucleation
    )
    return params
end

function EXP19_IC(FT)
    # Cotten et al (2007)
    tps = TD.Parameters.ThermodynamicsParameters(FT)
    wps = CMP.WaterProperties(FT)

    ρₗ = wps.ρw
    ρᵢ = wps.ρi
    R_d = TD.Parameters.R_d(tps)
    R_v = TD.Parameters.R_v(tps)

    Nᵢ = FT(161515.09922)
    Nₗ = FT(129.2 * 1e6) - Nᵢ
    Nₐ = FT(0)
    r₀ = FT(8 * 1e-6)
    p₀ = FT(99766.35514)
    T₀ = FT(251.721428)
    qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))  # 1.2 should be ρₐ
    qᵢ = FT(Nᵢ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))
    m_l = Nₗ * ρₗ *  4 * π / 3 * r₀^3
    m_i = Nᵢ * ρᵢ *  4 * π / 3 * r₀^3
    e = FT(57.103126)
    qᵥ = (e / R_v / T₀) / ((p₀ - e) / (R_d * T₀) + e / R_v / T₀ + m_l + m_i)
    eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
    Sₗ = FT(e / eₛ)

    return [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, FT(0)]
end

function EXP45_params(FT, w, t_max, t_profile, T_profile, P_profile, batch_name)
    # Cotten et al (2007)
    IN_mode = "ABDINM"
    const_dt = FT(1)
    prescribed_thermodynamics = true
    aerosol_act = "None"
    aerosol = batch_name == "EXP45" ? CMP.SaharanDust(FT) : CMP.Dust(FT)
    dep_nucleation = "ABDINM"
    heterogeneous = "ABIFM"
    homogeneous = "ABHOM"
    condensation_growth = "Condensation"
    deposition_growth = "Deposition"
    liq_size_distribution = "Gamma"
    ice_size_distribution = "Gamma"
    aero_σ_g = FT(log(1.75))     # converted std dev to geometric std dev
    r_nuc = FT(0.4 / 2 * 1e-6)   # value is mode radius, not mean
    A_aer = FT(4 * π * r_nuc^2)
    ips = CMP.IceNucleationParameters(FT)

    params = (; const_dt, w, t_max, ips,
        prescribed_thermodynamics, t_profile, T_profile, P_profile,
        aerosol_act, aerosol, r_nuc, aero_σ_g, A_aer,   # aerosol activation
        condensation_growth, deposition_growth,         # growth
        liq_size_distribution, ice_size_distribution,   # size distribution
        dep_nucleation, heterogeneous, homogeneous,     # ice nucleation
    )
    return params
end

function EXP45_IC(FT)
    # Cotten et al (2007)
    tps = TD.Parameters.ThermodynamicsParameters(FT)
    wps = CMP.WaterProperties(FT)

    ρₗ = wps.ρw
    ρᵢ = wps.ρi
    R_d = TD.Parameters.R_d(tps)
    R_v = TD.Parameters.R_v(tps)

    Nₗ = FT(0)
    Nᵢ = FT(748663.1016)
    Nₐ = FT(9 * 1e6) - Nᵢ- Nₗ
    r₀ = FT(10 * 1e-6)
    p₀ = FT(99583.33333)
    T₀ = FT(213.351)
    qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))  # 1.2 should be ρₐ
    qᵢ = FT(Nᵢ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2))
    m_l = Nₗ * ρₗ *  4 * π / 3 * r₀^3
    m_i = Nᵢ * ρᵢ *  4 * π / 3 * r₀^3
    e = FT(1.08486)
    qᵥ = (e / R_v / T₀) / ((p₀ - e) / (R_d * T₀) + e / R_v / T₀ + m_l + m_i)
    eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
    Sₗ = FT(e / eₛ)

    return [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, FT(0)]
end

#! format: on
