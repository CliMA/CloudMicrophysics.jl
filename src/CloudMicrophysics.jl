module CloudMicrophysics

import Thermodynamics
import Thermodynamics.ThermodynamicsParameters

import CLIMAParameters

#types of Microphysics
abstract type AbstractMicrophysicsParameters end

struct NoMicrophysicsParameters <: AbstractMicrophysicsParameters end #just for testing

struct Microphysics_0M_Parameters{FT} <: AbstractMicrophysicsParameters
    τ_precip::FT
    qc_0::FT
    S_0::FT
end
function Microphysics_0M_Parameters(param_set)
    aliases = ["τ_precip","qc_0","S_0"]
    (τ_precip, qc_0, S_0) = CLIMAParameters.get_parameter_values!(param_set,aliases, Microphysics_0M_Parameters)
    return Microphysics_0M_Parameters{CLIMAParameters.get_parametric_type(param_set)}(τ_precip, qc_0, S_0)
end


struct Microphysics_1M_Parameters{FT} <: AbstractMicrophysicsParameters
    C_drag::FT
    K_therm::FT
    D_vapor::FT
    ν_air::FT
    τ_cond_evap::FT
    τ_sub_dep::FT
    r_ice_snow::FT
    n0_ice::FT
    r0_ice::FT
    me_ice::FT
    χm_ice::FT
    Δm_ice::FT
    a_vent_rai::FT
    b_vent_rai::FT
    n0_rai::FT
    r0_rai::FT
    me_rai::FT
    ae_rai::FT
    ve_rai::FT
    χm_rai::FT
    Δm_rai::FT
    χa_rai::FT
    Δa_rai::FT
    χv_rai::FT
    Δv_rai::FT
    τ_acnv_rai::FT
    q_liq_threshold::FT
    q_ice_threshold::FT
    τ_acnv_sno::FT
    a_vent_sno::FT
    b_vent_sno::FT
    ν_sno::FT
    μ_sno::FT
    r0_sno::FT
    me_sno::FT
    ae_sno::FT
    ve_sno::FT
    χm_sno::FT
    Δm_sno::FT
    χa_sno::FT
    Δa_sno::FT
    χv_sno::FT
    Δv_sno::FT
    E_liq_rai::FT
    E_liq_sno::FT
    E_ice_rai::FT
    E_ice_sno::FT
    E_rai_sno::FT
    ρ_cloud_ice::FT
    ρ_cloud_liq::FT
    grav::FT
    T_freeze::FT
    gas_constant::FT
    molmass_water::FT
    N_Sc::FT
    m0_ice::FT
    m0_rai::FT
    a0_rai::FT
    m0_sno::FT
    a0_sno::FT
    v0_sno::FT
    R_v::FT
    TPS::ThermodynamicsParameters{FT}
end
function Microphysics_1M_Parameters(
    param_set,
    TPS::ThermodynamicsParameters{FT},
) where {FT}
    
    aliases = ["C_drag", "K_therm", "D_vapor", "ν_air", "τ_cond_evap",
               "τ_sub_dep", "r_ice_snow", "n0_ice", "r0_ice", "me_ice",
               "ρ_cloud_ice", "χm_ice", "Δm_ice", "q_liq_threshold",
               "τ_acnv_rai", "a_vent_rai", "b_vent_rai", "n0_rai",
               "r0_rai", "me_rai", "ae_rai", "ve_rai", "χm_rai", "Δm_rai",
               "χa_rai", "Δa_rai", "χv_rai", "Δv_rai", "q_ice_threshold",
               "τ_acnv_sno", "a_vent_sno", "b_vent_sno", "ν_sno", "μ_sno",
               "r0_sno", "me_sno", "ae_sno", "ve_sno", "χm_sno", "Δm_sno",
               "χa_sno", "Δa_sno", "χv_sno", "Δv_sno", "E_liq_rai", "E_liq_sno",
               "E_ice_rai", "E_ice_sno", "E_rai_sno", "ρ_cloud_liq",
               "ρ_cloud_ice", "grav", "T_freeze", "gas_constant", "molmass_water"] 

    (C_drag, K_therm, D_vapor, ν_air, τ_cond_evap, τ_sub_dep,
     r_ice_snow, n0_ice, r0_ice, me_ice, ρ_cloud_ice, χm_ice,
     Δm_ice, q_liq_threshold, τ_acnv_rai, a_vent_rai, b_vent_rai,
     n0_rai, r0_rai, me_rai, ae_rai, ve_rai, χm_rai, Δm_rai, χa_rai,
     Δa_rai, χv_rai, Δv_rai, q_ice_threshold, τ_acnv_sno, a_vent_sno,
     b_vent_sno, ν_sno, μ_sno, r0_sno, me_sno, ae_sno, ve_sno,
     χm_sno, Δm_sno, χa_sno, Δa_sno, χv_sno, Δv_sno, E_liq_rai,
     E_liq_sno, E_ice_rai, E_ice_sno, E_rai_sno, ρ_cloud_liq,
     ρ_cloud_ice, grav, T_freeze, gas_constant, molmass_water) = CLIMAParameters.get_parameter_values!(param_set,aliases,Microphysics_1M_Parameters)

    #derived parameters
    N_Sc = ν_air / D_vapor
    m0_ice = 4 / 3 * π * ρ_cloud_ice * r0_ice^me_ice
    m0_rai = 4 / 3 * π * ρ_cloud_liq * r0_rai^me_rai
    a0_rai = π * r0_rai^ae_rai
    m0_sno = 1e-1 * r0_sno^me_sno
    a0_sno = 0.3 * π * r0_sno^ae_sno
    v0_sno = 2^(9 / 4) * r0_sno^ve_sno
    R_v = gas_constant / molmass_water


    return Microphysics_1M_Parameters{CLIMAParameters.get_parametric_type(param_set)}(
        C_drag,
        K_therm,
        D_vapor,
        ν_air,
        τ_cond_evap,
        τ_sub_dep,
        r_ice_snow,
        n0_ice,
        r0_ice,
        me_ice,
        χm_ice,
        Δm_ice,
        a_vent_rai,
        b_vent_rai,
        n0_rai,
        r0_rai,
        me_rai,
        ae_rai,
        ve_rai,
        χm_rai,
        Δm_rai,
        χa_rai,
        Δa_rai,
        χv_rai,
        Δv_rai,
        τ_acnv_rai,
        q_liq_threshold,
        q_ice_threshold,
        τ_acnv_sno,
        a_vent_sno,
        b_vent_sno,
        ν_sno,
        μ_sno,
        r0_sno,
        me_sno,
        ae_sno,
        ve_sno,
        χm_sno,
        Δm_sno,
        χa_sno,
        Δa_sno,
        χv_sno,
        Δv_sno,
        E_liq_rai,
        E_liq_sno,
        E_ice_rai,
        E_ice_sno,
        E_rai_sno,
        ρ_cloud_ice,
        ρ_cloud_liq,
        grav,
        T_freeze,
        gas_constant,
        molmass_water,
        N_Sc,
        m0_ice,
        m0_rai,
        a0_rai,
        m0_sno,
        a0_sno,
        v0_sno,
        R_v,
        TPS,
    )


end


#General parameters outside of modular
struct CloudMicrophysicsParameters{FT, AMPS <: AbstractMicrophysicsParameters}
    K_therm::FT
    D_vapor::FT
    molmass_dryair::FT
    molmass_water::FT
    gas_constant::FT
    ρ_cloud_liq::FT
    surface_tension_coeff::FT
    grav::FT
    molmass_ratio::FT
    R_v::FT
    MPS::AMPS
    TPS::ThermodynamicsParameters{FT}
end

# For example:
# CloudMicrophysicsParameters(
#     param_set,
#     Microphysics_1M_Parameters(param_set),
#     ThermodynamicsParameters(param_set)
# )

function CloudMicrophysicsParameters(
    param_set,
    MPS::AMPS,
    TPS::ThermodynamicsParameters{FT},
) where {FT, AMPS <: AbstractMicrophysicsParameters}

    aliases = ["K_therm", "D_vapor", "molmass_dryair", "molmass_water",
               "gas_constant", "ρ_cloud_liq", "surface_tension_coeff", "grav"]

    (K_therm, D_vapor, molmass_dryair, molmass_water,
     gas_constant,ρ_cloud_liq, surface_tension_coeff, grav) = CLIMAParameters.get_parameter_values!(param_set,aliases,CloudMicrophysicsParameters)
    
    #derived parameters 
    molmass_ratio = molmass_dryair / molmass_water
    R_v = gas_constant / molmass_water

    return CloudMicrophysicsParameters{CLIMAParameters.get_parametric_type(param_set), AMPS}(
        K_therm,
        D_vapor,
        molmass_dryair,
        molmass_water,
        gas_constant,
        ρ_cloud_liq,
        surface_tension_coeff,
        grav,
        molmass_ratio,
        R_v,
        MPS,
        TPS,
    )

end


include("Common.jl")
include("Microphysics_0M.jl")
include("Microphysics_1M.jl")
include("AerosolModel.jl")
include("AerosolActivation.jl")

end # module
