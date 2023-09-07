import OrdinaryDiffEq as ODE
import CairoMakie as MK

import Thermodynamics as TD
import CloudMicrophysics as CM

import CloudMicrophysics.CommonTypes as CMT
import CloudMicrophysics.Common as CMO
import CloudMicrophysics.HetIceNucleation as CMI_het
import CloudMicrophysics.Parameters as CMP

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

"""
    ODE problem definition
"""
function parcel_model(dY, Y, p, t)

    # Get simulation parameters
    (; prs, const_dt, r_nuc, w, α_m, options) = p
    (; freeze_mode, deposition_growth, condensation_growth) = options
    # Numerical precision used in the simulation
    FT = eltype(Y)

    # Our state vector
    S_i = Y[1]        # supersaturation over ice
    N_act = Y[2]      # number concentration of activated particles
    p_a = Y[3]        # pressure
    T = Y[4]          # temperature
    q_vap = Y[5]      # vapor specific humidity
    q_liq = Y[6]      # liquid water specific humidity
    q_ice = Y[7]      # ice specific humidity
    N_aerosol = Y[8]  # number concentration of interstitial aerosol
    N_droplets = Y[9] # number concentration of existing droplets
    x_sulph = Y[10]   # percent mass sulphuric acid

    # Constants
    R_v = CMP.R_v(prs)
    grav = CMP.grav(prs)
    ρ_ice = CMP.ρ_cloud_ice(prs)
    ρ_liq = CMP.ρ_cloud_liq(prs)


    # Get thermodynamic parameters, phase partition and create thermo state.
    thermo_params = CMP.thermodynamics_params(prs)
    q = TD.PhasePartition(q_vap + q_ice, FT(0), q_ice) # TODO add liquid
    ts = TD.PhaseNonEquil_pTq(thermo_params, p_a, T, q)

    # Constants and variables that depend on the moisture content
    R_a = TD.gas_constant_air(thermo_params, q)
    cp_a = TD.cp_m(thermo_params, q)
    L_subl = TD.latent_heat_sublim(thermo_params, T)
    L_fus = TD.latent_heat_fusion(thermo_params, T)
    ρ = TD.air_density(thermo_params, ts)

    # We are solving for both q_ and supersaturation. This is unneccessary.
    # I'm checking here if we are consistent. We should decide
    # if we want to solve for S (as in the cirrucs box model) or for q_
    # which is closer to what ClimaAtmos is doing.
    e_si = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())
    e = S_i * e_si
    e_check = q_vap * p_a * R_v / R_a
    @assert isapprox(e, e_check; rtol = 0.01)

    # Adiabatic parcel coefficients
    a1 = L_subl * grav / cp_a / T^2 / R_v - grav / R_a / T
    a2 = p_a / e_si * R_v / R_a
    a3 = L_subl^2 / R_v / T^2 / cp_a
    a4 = L_subl * L_fus / R_v / T^2 / cp_a

    dN_act_dt_deposition = FT(0)
    dN_act_dt_ABIFM = FT(0)
    dN_act_dt_homogeneous = FT(0)

    # Activating new crystals
    τ_relax = const_dt
    if freeze_mode == "no freezing"
        @warn "No freezing will occur."

    elseif freeze_mode == "deposition"
        AF = CMI_het.dust_activated_number_fraction(
            prs,
            S_i,
            T,
            CMT.DesertDustType(),
        )

        dN_act_dt_deposition = max(FT(0), AF * N_aerosol - N_act) / τ_relax
    else
        @warn "Invalid freezing mode in run_parcel argument. Running without freezing.\nPlease choose between \"no freezing\" or \"deposition\" (include quotation marks)."

    end
    dN_act_dt = dN_act_dt_deposition
    dN_aerosol_dt = -dN_act_dt
    dqi_dt_new_particles = dN_act_dt * 4 / 3 * π * r_nuc^3 * ρ_ice / ρ

    # TODO - Activating new droplets
    dN_droplets_dt = FT(0)

    # Growth
    # Deposition on existing crystals (assuming all are the same...)
    G_ice = CMO.G_func(prs, T, TD.Ice())
    r_ice = N_act > 0 ? cbrt(q_ice / N_act / (4 / 3 * π) / ρ_ice * ρ) : FT(0)
    C_ice = r_ice
    dqi_dt_deposition =
        deposition_growth == true ?
        1 / ρ * N_act * α_m * 4 * π * C_ice * (S_i - 1) * G_ice : FT(0)
    # Condensation on existing droplets (assuming all are the same)
    G_liq = CMO.G_func(prs, T, TD.Liquid())
    r_liq = N_droplets > 0 ? cbrt(q_liq / N_droplets / (4 / 3 * π) / ρ_liq * ρ) : FT(0)
    C_liq = r_liq
    dql_dt_condensation =
        condensation_growth == true ?
        1 / ρ * N_droplets * α_m * 4 * π * C_liq * (S_i - 1) * G_liq : FT(0)

    # Sum of all phase changes
    dqi_dt = dqi_dt_new_particles + dqi_dt_deposition
    dql_dt = dN_droplets_dt + dql_dt_condensation

    # Update the tendecies
    dS_i_dt = a1 * w * S_i - (a2 + a3 * S_i) * dqi_dt - (a2 + a4 * S_i) * dql_dt
    dp_a_dt = -p_a * grav / R_a / T * w
    dT_dt = -grav / cp_a * w + L_subl / cp_a * dqi_dt
    dq_vap_dt = -dqi_dt    # TODO - balance with liquid
    dq_ice_dt = dqi_dt
    dq_liq_dt = dql_dt
    x_sulph_dt = FT(0)

    # Set tendencies
    dY[1] = dS_i_dt        # supersaturation over ice
    dY[2] = dN_act_dt      # number concentration of activated particles
    dY[3] = dp_a_dt        # pressure
    dY[4] = dT_dt          # temperature
    dY[5] = dq_vap_dt      # vapor specific humidity
    dY[6] = dq_liq_dt      # liquid water specific humidity
    dY[7] = dq_ice_dt      # ice specific humidity
    dY[8] = dN_aerosol_dt  # number concentration of interstitial aerosol
    dY[9] = dN_droplets_dt # mumber concentration of droplets
    dY[10] = x_sulph_dt    # sulphuric acid concentration

    # TODO - add diagnostics output (radius, S, etc)
end
