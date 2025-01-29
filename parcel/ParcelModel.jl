import Thermodynamics as TD
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.AerosolModel as AM

include(joinpath(pkgdir(CM), "parcel", "ParcelParameters.jl"))

"""
    Parcel simulation parameters
"""
Base.@kwdef struct parcel_params{FT} <: CMP.ParametersType{FT}
    prescribed_thermodynamics = false
    t_profile = []
    T_profile = []
    P_profile = []
    aerosol_act = "None"
    deposition = "None"
    heterogeneous = "None"
    homogeneous = "None"
    condensation_growth = "None"
    deposition_growth = "None"
    liq_size_distribution = "Monodisperse"
    ice_size_distribution = "Monodisperse"
    aerosol = Empty{FT}()
    aero_σ_g = FT(0)
    wps = CMP.WaterProperties(FT)
    aps = CMP.AirProperties(FT)
    tps = TD.Parameters.ThermodynamicsParameters(FT)
    aap = CMP.AerosolActivationParameters(FT)
    ips = CMP.IceNucleationParameters(FT)
    liquid = CMP.CloudLiquid(FT)
    ice = CMP.CloudIce(FT)
    H₂SO₄ps = CMP.H2SO4SolutionParameters(FT)
    const_dt = 1
    w = FT(1)
    r_nuc = FT(0.5 * 1.e-4 * 1e-6)
    A_aer = FT(1e-9)
    sampling_interval = FT(1)
    γ = FT(1)
    ip = CMP.Frostenberg2023(FT)
end

"""
    If negative, clip the tracers to zero when computing tendencies
"""
clip!(x, lim) = max(x, lim)
clip!(x) = clip!(x, eltype(x)(0))

"""
    ODE problem definitions
"""
function parcel_model(dY, Y, p, t)
    # Numerical precision used in the simulation
    FT = eltype(Y)
    # Simulation parameters
    (; prescribed_thermodynamics, t_profile, T_profile, P_profile) = p
    (; wps, tps, r_nuc, w) = p
    (; liq_distr, ice_distr) = p
    (; aero_act_params, dep_params, imm_params, hom_params) = p
    (; ce_params, ds_params) = p
    # Y values stored in a named tuple for ease of use
    state = (
        Sₗ = Y[1],
        p_air = Y[2],
        T = Y[3],
        qᵥ = clip!(Y[4]),
        qₗ = clip!(Y[5]),
        qᵢ = clip!(Y[6]),
        Nₐ = clip!(Y[7]),
        Nₗ = clip!(Y[8]),
        Nᵢ = clip!(Y[9]),
        ln_INPC = Y[10],   # needed only in stochastic Frostenberg
        t = t,
    )

    # Constants
    Rᵥ = TD.Parameters.R_v(tps)
    grav = TD.Parameters.grav(tps)
    ρᵢ = wps.ρi
    ρₗ = wps.ρw

    # Get the state values
    (; Sₗ, p_air, T, qᵥ, qₗ, qᵢ, Nₗ, Nᵢ, t) = state
    # Get thermodynamic parameters, phase partition and create thermo state.
    q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
    ts = TD.PhaseNonEquil_pTq(tps, p_air, T, q)

    # Constants and variables that depend on the moisture content
    R_air = TD.gas_constant_air(tps, q)
    cp_air = TD.cp_m(tps, q)
    L_subl = TD.latent_heat_sublim(tps, T)
    L_fus = TD.latent_heat_fusion(tps, T)
    L_vap = TD.latent_heat_vapor(tps, T)
    ρ_air = TD.air_density(tps, ts)

    # Vapor pressures
    e = qᵥ * p_air * Rᵥ / R_air
    eₛₗ = TD.saturation_vapor_pressure(tps, T, TD.Liquid())

    # Mean radius, area and volume of liquid droplets and ice crystals
    PSD_liq = distribution_moments(liq_distr, qₗ, Nₗ, ρₗ, ρ_air)
    PSD_ice = distribution_moments(ice_distr, qᵢ, Nᵢ, ρᵢ, ρ_air)

    # Aerosol activation
    dNₗ_dt_act = aerosol_activation(aero_act_params, state)
    dqₗ_dt_act = dNₗ_dt_act * 4 * FT(π) / 3 * r_nuc^3 * ρₗ / ρ_air

    # Deposition ice nucleation
    # (All deposition parameterizations assume monodisperse aerosol size distr)
    dNᵢ_dt_dep = deposition_nucleation(dep_params, state, dY)
    dqᵢ_dt_dep = dNᵢ_dt_dep * 4 / 3 * FT(π) * r_nuc^3 * ρᵢ / ρ_air

    # Heterogeneous ice nucleation
    #@info(imm_params)
    dln_INPC_imm = INPC_model(imm_params, state)
    dNᵢ_dt_imm = immersion_freezing(imm_params, PSD_liq, state)
    dqᵢ_dt_imm = dNᵢ_dt_imm * PSD_liq.V * ρᵢ / ρ_air

    # Homogeneous ice nucleation
    dNᵢ_dt_hom = homogeneous_freezing(hom_params, PSD_liq, state)
    dqᵢ_dt_hom = dNᵢ_dt_hom * PSD_liq.V * ρᵢ / ρ_air

    # Condensation/evaporation
    dqₗ_dt_ce = condensation(ce_params, PSD_liq, state, ρ_air)
    # Deposition/sublimation
    dqᵢ_dt_ds = deposition(ds_params, PSD_ice, state, ρ_air)

    # number concentration and ...
    dNᵢ_dt = dNᵢ_dt_dep + dNᵢ_dt_imm + dNᵢ_dt_hom
    dNₐ_dt = -dNᵢ_dt_dep - dNₗ_dt_act
    dNₗ_dt = dNₗ_dt_act - dNᵢ_dt_imm - dNᵢ_dt_hom
    # ... water mass budget
    dqₗ_dt_v2l = dqₗ_dt_ce + dqₗ_dt_act
    dqᵢ_dt_l2i = dqᵢ_dt_imm + dqᵢ_dt_hom
    dqᵢ_dt_v2i = dqᵢ_dt_dep + dqᵢ_dt_ds

    # Update the tendecies
    dqᵢ_dt = dqᵢ_dt_v2i + dqᵢ_dt_l2i
    dqₗ_dt = dqₗ_dt_v2l - dqᵢ_dt_l2i
    dqᵥ_dt = -dqₗ_dt_v2l - dqᵢ_dt_v2i

    dp_air_dt =
        prescribed_thermodynamics ? AIDA_rate(t, t_profile, P_profile) :
        -p_air * grav / R_air / T * w

    dT_dt =
        prescribed_thermodynamics ? AIDA_rate(t, t_profile, T_profile) :
        -grav / cp_air * w +
        L_vap / cp_air * dqₗ_dt_v2l +
        L_fus / cp_air * dqᵢ_dt_l2i +
        L_subl / cp_air * dqᵢ_dt_v2i

    de_dt = dqᵥ_dt * p_air * Rᵥ / R_air + qᵥ * dp_air_dt * Rᵥ / R_air
    deₛₗ_dt = L_vap * eₛₗ / Rᵥ / T^2 * dT_dt
    dSₗ_dt = 1 / eₛₗ * de_dt - e / (eₛₗ)^2 * deₛₗ_dt

    # Set tendencies
    dY[1] = dSₗ_dt      # saturation ratio over liquid water
    dY[2] = dp_air_dt   # pressure
    dY[3] = dT_dt       # temperature
    dY[4] = dqᵥ_dt      # vapor specific humidity
    dY[5] = dqₗ_dt      # liquid water specific humidity
    dY[6] = dqᵢ_dt      # ice specific humidity
    dY[7] = dNₐ_dt      # number concentration of interstitial aerosol
    dY[8] = dNₗ_dt      # mumber concentration of droplets
    dY[9] = dNᵢ_dt      # number concentration of activated particles
    dY[10] = dln_INPC_imm
end

"""
    run_parcel(IC, t_0, t_end, pp)

Returns the solution of an ODE probelm defined by the parcel model.

Inputs:
 - IC - Vector with initial conditions: [Sₗ, p_air, T, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln(INPC)]
 - t_0 - simulation start time
 - t_end - simulation end time
 - pp - struct with parcel simulation parameters

Parcel state vector (all variables are in base SI units):
 - Sₗ    - saturation ratio over liquid water
 - p_air - air pressure
 - T     - temperature
 - qᵥ    - vapor specific humidity
 - qₗ    - liquid water specific humidity
 - qᵢ    - ice specific humidity
 - Nₐ    - number concentration of interstitial aerosol
 - Nₗ    - number concentration of existing water droplets
 - Nᵢ    - concentration of activated ice crystals
 - ln(INPC) - logarithm of ice nucleation particles for the stochastic het. freezing option
 - t     - time in seconds since the beginning of the simulation

The parcel parameters struct comes with default values that can be overwritten:
 - deposition - string with deposition ice nucleation parameterization choice ["None" (default), "MohlerAF", "MohlerRate", "ActivityBased", "P3_dep"]
 - heterogeneous - string with heterogeneous ice nucleation parameterization choice ["None" (default), "ImmersionFreezing", "P3_het"]
 - homogeneous - string with homogeneous ice nucleation parameterization choice ["None" (default), "ActivityBased", "P3_hom"]
 - condensation_growth - string with condensation/evaporation options for cloud droplets ["None" (default), "Condensation"]
 - deposition_growth - string with deposition/sublimation options for ice crystals ["None" (default), "Deposition"]
 - size_distribution - string with cloud droplet and ice crystal size disribution choice ["Monodisperse" (default), "Gamma"]
 - aerosol - a struct with aerosol parameters required by the nucleation parameterizations, see CloudMicrophysics documentation for all the options. The default is an Empty struct.
 - wps, aps, tps, ips, H₂SO₄ps - structs with additional free parameters needed by the parameterizations. By default we use the values stored in ClimaParams.jl. See CloudMicrophysics docs for more details.
 - const_dt - model timestep [s]. Parcel model is using a simple Euler timestepper. Default value is 1 s
 - w - parcel vertical velocity [m/s]. Default value is 1 m/s
 - r_nuc - assumed size of nucleating ice crystals. Default value is 5e-11 [m]
 - A_aer - assumed surface area of ice nucleating aerosol. Default value assumes radius of 5e-11 [m]
 - ip - parameters of INPC frequency distribution for Frostenberg parametrization of immerison freezing
 - sampling_interval - number of time steps between random draws in Frostenberg_random parametrization of immerison freezing
 - γ - inverse timescale for the stochastic process in Frostenberg_stochastic parametrization of immerison freezing
"""
function run_parcel(IC, t_0, t_end, pp)

    FT = eltype(IC)

    info = "\nSize distribution (liq): $(pp.liq_size_distribution)\n"
    if pp.liq_size_distribution == "Monodisperse"
        liq_distr = Monodisperse{FT}()
    elseif pp.liq_size_distribution == "Gamma"
        liq_distr = Gamma{FT}()
    else
        throw("Unrecognized size distribution")
    end

    info *= "Size distribution (ice): $(pp.ice_size_distribution)\n"
    if pp.ice_size_distribution == "Monodisperse"
        ice_distr = Monodisperse{FT}()
    elseif pp.ice_size_distribution == "Gamma"
        ice_distr = Gamma{FT}()
    else
        throw("Unrecognized size distribution")
    end

    info *= "Aerosol: $(chop( string(typeof(pp.aerosol)), head = 29, tail = 9))\n"

    info *= "Aerosol Activation: $(pp.aerosol_act)\n"
    if pp.aerosol_act == "None"
        aero_act_params = Empty{FT}()
    elseif pp.aerosol_act == "AeroAct"
        aero_act_params =
            AeroAct{FT}(pp.aap, pp.aerosol, pp.aero_σ_g, pp.r_nuc, pp.const_dt)
    else
        throw("Unrecognized aerosol activation mode")
    end

    info *= "Deposition: $(pp.deposition)\n"
    if pp.deposition == "None"
        dep_params = Empty{FT}()
    elseif pp.deposition == "MohlerAF"
        dep_params = MohlerAF{FT}(pp.ips, pp.aerosol, pp.tps, pp.const_dt)
    elseif pp.deposition == "MohlerRate"
        dep_params = MohlerRate{FT}(pp.ips, pp.aerosol, pp.tps, pp.const_dt)
    elseif pp.deposition == "ABDINM"
        dep_params = ABDINM{FT}(pp.tps, pp.aerosol, pp.r_nuc, pp.const_dt)
    elseif pp.deposition == "P3_dep"
        dep_params = P3_dep{FT}(pp.ips, pp.const_dt)
    else
        throw("Unrecognized deposition mode")
    end

    info *= "Heterogeneous: $(pp.heterogeneous)\n"
    if pp.heterogeneous == "None"
        imm_params = Empty{FT}()
    elseif pp.heterogeneous == "ABIFM"
        imm_params = ABIFM{FT}(pp.tps, pp.aerosol, pp.A_aer, pp.const_dt)
    elseif pp.heterogeneous == "P3_het"
        imm_params = P3_het{FT}(pp.ips, pp.const_dt)
    elseif pp.heterogeneous == "Frostenberg_random"
        imm_params =
            Frostenberg_random{FT}(pp.ip, pp.sampling_interval, pp.const_dt)
    elseif pp.heterogeneous == "Frostenberg_mean"
        imm_params = Frostenberg_mean{FT}(pp.ip, pp.const_dt)
    elseif pp.heterogeneous == "Frostenberg_stochastic"
        imm_params = Frostenberg_stochastic{FT}(pp.ip, pp.γ, pp.const_dt)
    else
        throw("Unrecognized heterogeneous mode")
    end

    info *= "Homogeneous: $(pp.homogeneous)\n"
    if pp.homogeneous == "None"
        hom_params = Empty{FT}()
    elseif pp.homogeneous == "ABHOM"
        hom_params = ABHOM{FT}(pp.tps, pp.ips, pp.const_dt)
    elseif pp.homogeneous == "P3_hom"
        hom_params = P3_hom{FT}(pp.const_dt)
    else
        throw("Unrecognized homogeneous mode")
    end

    info *= "Condensation growth: $(pp.condensation_growth)\n"
    if pp.condensation_growth == "None"
        ce_params = Empty{FT}()
    elseif pp.condensation_growth == "Condensation"
        ce_params = CondParams{FT}(pp.aps, pp.tps, pp.const_dt)
    elseif pp.condensation_growth == "NonEq_Condensation_Anna"
        ce_params = NonEqCondParams_Anna{FT}(pp.tps, pp.liquid)
    elseif pp.condensation_growth == "NonEq_Condensation"
        ce_params = NonEqCondParams{FT}(pp.tps, pp.liquid, pp.ice)
    else
        throw("Unrecognized condensation growth mode")
    end

    info *= "Deposition growth: $(pp.deposition_growth)\n"
    if pp.deposition_growth == "None"
        ds_params = Empty{FT}()
    elseif pp.deposition_growth == "Deposition"
        ds_params = DepParams{FT}(pp.aps, pp.tps)
    elseif pp.deposition_growth == "NonEq_Deposition_Anna"
        ds_params = NonEqDepParams_Anna{FT}(pp.tps, pp.ice)
    elseif pp.deposition_growth == "NonEq_Deposition"
        ds_params = NonEqDepParams{FT}(pp.tps, pp.liquid, pp.ice)
    else
        throw("Unrecognized deposition growth mode")
    end
    # To see parcel parameters in output, uncomment
    # @info info

    # Parameters for the ODE solver
    p = (
        prescribed_thermodynamics = pp.prescribed_thermodynamics,
        t_profile = pp.t_profile,
        T_profile = pp.T_profile,
        P_profile = pp.P_profile,
        liq_distr = liq_distr,
        ice_distr = ice_distr,
        aero_act_params = aero_act_params,
        dep_params = dep_params,
        imm_params = imm_params,
        hom_params = hom_params,
        ce_params = ce_params,
        ds_params = ds_params,
        wps = pp.wps,
        tps = pp.tps,
        aps = pp.aps,
        aap = pp.aap,
        r_nuc = pp.r_nuc,
        w = pp.w,
    )

    problem = ODE.ODEProblem(parcel_model, IC, (FT(t_0), FT(t_end)), p)
    return ODE.solve(
        problem,
        ODE.Euler(),
        dt = pp.const_dt,
        reltol = 100 * eps(FT),
        abstol = 100 * eps(FT),
    )
end
