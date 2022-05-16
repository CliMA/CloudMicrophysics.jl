module CloudMicrophysics

import Thermodynamics
import Thermodynamics.ThermodynamicsParameters

import CLIMAParameters
const CP = CLIMAParameters

#parameters for moisture Microphysics
abstract type AbstractMoistureParameters end
struct EqMoistureParameters <: AbstractMoistureParameters end # unused
struct NonEqMoistureParameters{FT} <: AbstractMoistureParameters
    τ_cond_evap::FT
    τ_sub_dep::FT
end
function NonEqMoistureParameters(param_struct)
    aliases = ["τ_cond_evap","τ_sub_dep"]
    (τ_cond_evap,τ_sub_dep) = CP.get_parameter_values!(
        param_struct,
        aliases,
        "NonEqMoisture",
    )
    return NonEqMoistureParameters{
        CP.get_parametric_type(param_struct),
    }(
        τ_cond_evap,
        τ_sub_dep,
    )
end
    

#types of Microphysics
abstract type AbstractPrecipitationParameters end

struct NoMicrophysicsParameters <: AbstractPrecipitationParameters end #just for testing

struct Microphysics_0M_Parameters{FT} <: AbstractPrecipitationParameters
    τ_precip::FT
    qc_0::FT
    S_0::FT
end
function Microphysics_0M_Parameters(param_struct)
    aliases = ["τ_precip", "qc_0", "S_0"]
    (τ_precip, qc_0, S_0) = CP.get_parameter_values!(
        param_struct,
        aliases,
        "Microphysics_0M",
    )
    return Microphysics_0M_Parameters{
        CP.get_parametric_type(param_struct),
    }(
        τ_precip,
        qc_0,
        S_0,
    )
end


struct Microphysics_1M_Parameters{FT} <: AbstractPrecipitationParameters
    C_drag::FT
    K_therm::FT
    D_vapor::FT
    ν_air::FT
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
    param_struct,
    TPS::ThermodynamicsParameters{FT},
) where {FT}

    aliases = [
        "C_drag",
        "K_therm",
        "D_vapor",
        "ν_air",
        "r_ice_snow",
        "n0_ice",
        "r0_ice",
        "me_ice",
        "ρ_cloud_ice",
        "χm_ice",
        "Δm_ice",
        "q_liq_threshold",
        "τ_acnv_rai",
        "a_vent_rai",
        "b_vent_rai",
        "n0_rai",
        "r0_rai",
        "me_rai",
        "ae_rai",
        "ve_rai",
        "χm_rai",
        "Δm_rai",
        "χa_rai",
        "Δa_rai",
        "χv_rai",
        "Δv_rai",
        "q_ice_threshold",
        "τ_acnv_sno",
        "a_vent_sno",
        "b_vent_sno",
        "ν_sno",
        "μ_sno",
        "r0_sno",
        "me_sno",
        "ae_sno",
        "ve_sno",
        "χm_sno",
        "Δm_sno",
        "χa_sno",
        "Δa_sno",
        "χv_sno",
        "Δv_sno",
        "E_liq_rai",
        "E_liq_sno",
        "E_ice_rai",
        "E_ice_sno",
        "E_rai_sno",
        "ρ_cloud_liq",
        "ρ_cloud_ice",
        "grav",
        "T_freeze",
        "gas_constant",
        "molmass_water",
    ]

    (
        C_drag,
        K_therm,
        D_vapor,
        ν_air,
        r_ice_snow,
        n0_ice,
        r0_ice,
        me_ice,
        ρ_cloud_ice,
        χm_ice,
        Δm_ice,
        q_liq_threshold,
        τ_acnv_rai,
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
        ρ_cloud_liq,
        ρ_cloud_ice,
        grav,
        T_freeze,
        gas_constant,
        molmass_water,
    ) = CP.get_parameter_values!(
        param_struct,
        aliases,
        "Microphysics_1M",
    )

    #derived parameters
    N_Sc = ν_air / D_vapor
    m0_ice = 4 / 3 * π * ρ_cloud_ice * r0_ice^me_ice
    m0_rai = 4 / 3 * π * ρ_cloud_liq * r0_rai^me_rai
    a0_rai = π * r0_rai^ae_rai
    m0_sno = 1e-1 * r0_sno^me_sno
    a0_sno = 0.3 * π * r0_sno^ae_sno
    v0_sno = 2^(9 / 4) * r0_sno^ve_sno
    R_v = gas_constant / molmass_water


    return Microphysics_1M_Parameters{
        CP.get_parametric_type(param_struct),
    }(
        C_drag,
        K_therm,
        D_vapor,
        ν_air,
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
struct CloudMicrophysicsParameters{FT, APPS, AMPS}
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
    PPS::APPS
    MPS::AMPS
    TPS::ThermodynamicsParameters{FT}
end

# For example:
# CloudMicrophysicsParameters(
#     param_struct,
#     Microphysics_1M_Parameters(param_struct),
#     ThermodynamicsParameters(param_struct)
# )

function CloudMicrophysicsParameters(
    param_struct,
    PPS::APPS,
    MPS::AMPS,
    TPS::ThermodynamicsParameters{FT},
) where {FT,
         APPS <: AbstractPrecipitationParameters,
         AMPS <: AbstractMoistureParameters}

    aliases = [
        "K_therm",
        "D_vapor",
        "molmass_dryair",
        "molmass_water",
        "gas_constant",
        "ρ_cloud_liq",
        "surface_tension_coeff",
        "grav",
    ]

    (
        K_therm,
        D_vapor,
        molmass_dryair,
        molmass_water,
        gas_constant,
        ρ_cloud_liq,
        surface_tension_coeff,
        grav,
    ) = CP.get_parameter_values!(
        param_struct,
        aliases,
        "CloudMicrophysics",
    )

    #derived parameters 
    molmass_ratio = molmass_dryair / molmass_water
    R_v = gas_constant / molmass_water

    return CloudMicrophysicsParameters{
        CP.get_parametric_type(param_struct),
        APPS,
        AMPS,
    }(
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
        PPS,
        MPS,
        TPS,
    )

end

include("CommonTypes.jl")

include("Common.jl")
include("Microphysics0M.jl")
include("Microphysics1M.jl")
include("MicrophysicsNonEq.jl")
include("AerosolModel.jl")
include("AerosolActivation.jl")

end # module
