import Thermodynamics as TD
import CloudMicrophysics as CM

const CMT = CM.CommonTypes
const CMO = CM.Common
const CMI = CM.HetIceNucleation
const CMP = CM.Parameters

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

"""
    ODE problem definition
"""
function parcel_model(dY, Y, p, t)

    # Get simulation parameters
    (; prs, const_dt, r_nuc, w, α_m) = p
    # Numerical precision used in the simulation
    FT = eltype(Y)

    # Our state vector
    S_i = Y[1]        # supersaturation over ice
    N_act = Y[2]      # number concentration of activated particles
    p_a = Y[3]        # pressure
    T = Y[4]          # temperature
    q_vap = Y[5]      # vapor specific humidity
    q_ice = Y[6]      # ice specific humidity
    N_aerosol = Y[7]  # number concentration of interstitial aerosol

    # Constants
    R_v = CMP.R_v(prs)
    grav = CMP.grav(prs)
    ρ_ice = CMP.ρ_cloud_liq(prs)

    # Get thermodynamic parameters, phase partition and create thermo state.
    thermo_params = CMP.thermodynamics_params(prs)
    q = TD.PhasePartition(q_vap + q_ice, FT(0), q_ice)
    ts = TD.PhaseNonEquil_pTq(thermo_params, p_a, T, q)

    # Constants and variables that depend on the moisture content
    R_a = TD.gas_constant_air(thermo_params, q)
    cp_a = TD.cp_m(thermo_params, q)
    L_subl = TD.latent_heat_sublim(thermo_params, T)
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

    # Activating new crystals
    AF = CMI.dust_activated_number_fraction(prs, S_i, T, CMT.DesertDustType())
    τ_relax = const_dt
    dN_act_dt = max(FT(0), AF * N_aerosol - N_act) / τ_relax
    dN_aerosol_dt = -dN_act_dt
    dqi_dt_new_particles = dN_act_dt * 4 / 3 * π * r_nuc^3 * ρ_ice / ρ

    # Growing existing crystals (assuming all are the same...)
    G = CMO.G_func(prs, T, TD.Ice())
    r = N_act > 0 ? cbrt(q_ice / N_act / (4 / 3 * π) / ρ_ice * ρ) : 0
    C = r
    dqi_dt_deposition = 1 / ρ * N_act * α_m * 4 * π * C * (S_i - 1) * G
    #dqi_dt_deposition = FT(0.0) # Use this if you don't want to grow

    # Sum of all phase changes
    dqi_dt = dqi_dt_new_particles + dqi_dt_deposition
    # TODO - update dqw_dt when implementing homo. and immersion freezing

    # Update the tendecies
    dS_i_dt = a1 * w * S_i - (a2 + a3 * S_i) * dqi_dt
    dp_a_dt = -p_a * grav / R_a / T * w
    dT_dt = -grav / cp_a * w + L_subl / cp_a * dqi_dt
    dq_vap_dt = -dqi_dt
    dq_ice_dt = dqi_dt
    # dq_liq_dt = dqw_dt # Use this when introducing liquid water

    # Set tendencies
    dY[1] = dS_i_dt        # supersaturation over ice
    dY[2] = dN_act_dt      # number concentration of activated particles
    dY[3] = dp_a_dt        # pressure
    dY[4] = dT_dt          # temperature
    dY[5] = dq_vap_dt      # vapor specific humidity
    dY[6] = dq_ice_dt      # ice specific humidity
    dY[7] = dN_aerosol_dt  # number concentration of interstitial aerosol
    # add dY state for dq_liq_dt when introducing liquid

    # TODO - add diagnostics output (radius, S, etc)
end
