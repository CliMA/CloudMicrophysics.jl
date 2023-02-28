import Thermodynamics
import CloudMicrophysics
import CLIMAParameters

const CP = CLIMAParameters
const TD = Thermodynamics
const CMT = CloudMicrophysics.CommonTypes
const CMI = CloudMicrophysics.HetIceNucleation
const CMP = CloudMicrophysics.Parameters
const APS = CMP.AbstractCloudMicrophysicsParameters

# Boiler plate code needed to have access to model parameters and constants
include(joinpath(pkgdir(CloudMicrophysics), "test", "create_parameters.jl"))
FT = Float64
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
const param_set = cloud_microphysics_parameters(toml_dict)

#q_tot = FT(10 * 1e-3)
#q_liq = FT(0)
#q_ice = FT(1 * 1e-3)
#q = TD.PhasePartition(q_tot, q_liq, q_ice)
#R_air = TD.gas_constant_air(thermo_params, q)
#print(R_air, "\n")

"""
    cirrus_box_model(param_set, T, p, S_i, w, n_ik, r_ik, n_step, dt)

 - `param_set` - set containing model parameters
 - `T` - temperature
 - `p` - pressure
 - `S_i` - supersaturation
 - `w` - vertical velocity
 - `n_ik` - ice nuclei number density
 - `r_ik` - ice particle radius
 - `n_step` - number of time steps
 - `dt` - simulation time step

Computes the cirrus box model evolution.
"""
function cirrus_box_model(
    param_set::APS,
    T::FT,
    p::FT,
    S_i::FT,
    w::FT,
    n_ik::FT,
    r_ik::FT,
    n_step::Int,
    dt::FT,
)
    # TODO - make dt a function of Courant number and w

    # all constants, variables, etc in SI units
    # careful with pressure units (i.e. edit if pressure input will be in hPa)

    thermo_params = CMP.thermodynamics_params(param_set)

    # constants
    R = CMP.gas_constant(param_set)
    k_B = CMP.k_Boltzmann(param_set)
    N_A = CMP.avogad(param_set)
    grav = CMP.grav(param_set)
    cp_water = CMP.cp_l(param_set)
    cp_dry_air = CMP.cp_d(param_set)
    M_dry_air = CMP.molmass_dryair(param_set) # TODO - should this be for dry air or total air (including)?
    M_water = CMP.molmass_water(param_set)
    D_vapor = CMP.D_vapor(param_set) # TODO -diffusion coefficient of water in air, check if this is the correct value

    for it in range(start=1, stop=n_step, step=1)

        # thermodynamic functions
        L_subl = TD.latent_heat_sublim(thermo_params, T)
        # TODO -which one do we need here?
        #p_sat = exp(16.3872)*exp(-3885.7/(T-237.15+230.17)); # saturated/vapor pressure of water, temp dependent. Antoine's Eqn, valid for 0-200C.
        p_sat_over_liquid = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
        p_sat_over_ice = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())

        ## TODO - which of those should be moved to ClimaParameters?
        #depCoeff = 0.5             # deposition coefficient of water molecules impinging on ice surface
        #capacitance = 1            # capacitance factor (accounts for geometry of ice crystal)
        #m_water = 3*10^(-23)/1000  # mass of one water molecule, kg/molecule
        #specVol_water = 0.001 / 1000 * M_water / N_A # specific volume of water, m^3/molecule
        #v_th_water = sqrt(8/pi*k_B/m_water*T)        # thermal speed of water, temp dependent, m/s
        #n_sat = p_sat_over_ice / T / k_B   # water vapor number density @ ice saturation, P_sat and temp dependent, from IG law, molecules/m^3
        #ventFactor = FT(0) # TODO - what value should we use for the ventilation factor?

        ## Coefficients from Karcher et al (2006)
        #a1 = ((L_subl * M_water * grav) / (cp_water * R* T^2)) - ((M_dry_air * grav) / (R * T)) # temperature dependence, m^-1
        #a2 = (n_sat)^-1;      # m^3/molecules
        #a3 = ((L_subl^2 * M_water * m_water) / (cp_water * M_dry_air)) / (T * p);         # temperature & pressure dependence, m^3/molecule
        #b1 = specVol_water * depCoeff * v_th_water * n_sat*((Si - 1) / 4);                # S_i dependence
        #b2 = depCoeff * v_th_water / 4 / capacitance / D_vapor

        # temperature change from dry adiabatic parcel
        T -= grav / cp_dry_air * w * dt

        ## number density of ice nuclei change
        #n_ik_new  = n_ik * (T_new / T)^(cp_water/R - 1) # TODO - is this cp_water, dry air, ice?

        ## radius of ice particle change
        #r_ik = ((1 + b2*prev_r_ik)*sqrt(1 + (2*ventFactor*b1*b2*dt)/((1+b2*prev_r_ik)^2)) - 1) / b2

        ## latent heat effects from phase change are incorporated as a downdraft that modifies the parcel velocity
        #w_p = (a2 + a3*prev_Si)/(a1*prev_Si) * (4*pi/specVol_water) * ((ventFactor*b1*r_ik^2)/(1+b2*prev_r_ik)) * n_ik; # for one type (if mulitple, must sum everything after 1st term)

        #effective_w = w - w_p          # effective velocity = updraft - downdraft

        #ΔSi = a1 * Si * effective_w * dt

        print("iter = ", it, " temperature = ", T, "\n")
    end
    return nothing
end

T = FT(200)
p = FT(80000.0)
S_i = FT(1.2)
w = FT(0.5)
n_ik = FT(0)
r_ik = FT(0)
n_step = 5
dt = FT(0.1)

cirrus_box_model(param_set, T, p, S_i, w, n_ik, r_ik, n_step, dt)
