module CloudMicrophysics

import Thermodynamics
import Thermodynamics.ThermodynamicsParameters

#types of Microphysics
abstract type AbstractMicrophysicsParameters end

struct NoMicrophysicsParameters <: AbstractMicrophysicsParameters end #just for testing

struct Microphysics_0M_Parameters{FT} <: AbstractMicrophysicsParameters
    τ_precip::FT
    qc_0::FT
    S_0::FT
end
function Microphysics_0M_Parameters(param_set::Dict)
    τ_precip = param_set["τ_precip"]
    qc_0 = param_set["qc_0"]
    S_0 = param_set["S_0"]
    return Microphysics_0M_Parameters{typeof(param_set["τ_precip"])}(
        τ_precip,
        qc_0,
        S_0,
    )
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
    N_Sc::FT
    m0_ice::FT
    m0_rai::FT
    a0_rai::FT       
    m0_sno::FT
    a0_sno::FT
    v0_sno::FT
end
function Microphysics_1M_Parameters(param_set::Dict)
    C_drag = param_set["C_drag"]
    K_therm = param_set["K_therm"]
    D_vapor = param_set["D_vapor"]
    ν_air = param_set["ν_air"]
    
    τ_cond_evap = param_set["τ_cond_evap"]
    τ_sub_dep = param_set["τ_sub_dep"]
    r_ice_snow = param_set["r_ice_snow"]
    n0_ice = param_set["n0_ice"]
    r0_ice = param_set["r0_ice"]
    me_ice = param_set["me_ice"]
    ρ_cloud_ice = param_set["ρ_cloud_ice"]
    χm_ice = param_set["χm_ice"]
    Δm_ice = param_set["Δm_ice"]

    q_liq_threshold = param_set["q_liq_threshold"]
    τ_acnv_rai = param_set["τ_acnv_rai"]
    a_vent_rai = param_set["a_vent_rai"]
    b_vent_rai = param_set["b_vent_rai"]
    n0_rai = param_set["n0_rai"]
    r0_rai = param_set["r0_rai"]
    me_rai = param_set["me_rai"]
    ae_rai = param_set["ae_rai"]
    ve_rai = param_set["ve_rai"]
    χm_rai = param_set["χm_rai"]
    Δm_rai = param_set["Δm_rai"]
    χa_rai = param_set["χa_rai"]
    Δa_rai = param_set["Δa_rai"]
    χv_rai = param_set["χv_rai"]
    Δv_rai = param_set["Δv_rai"]

    q_ice_threshold = param_set["q_ice_threshold"]
    τ_acnv_sno = param_set["τ_acnv_sno"]
    a_vent_sno = param_set["a_vent_sno"]
    b_vent_sno = param_set["b_vent_sno"]
    ν_sno = param_set["ν_sno"]
    μ_sno = param_set["μ_sno"]
    r0_sno = param_set["r0_sno"]
    me_sno = param_set["me_sno"]
    ae_sno = param_set["ae_sno"]
    ve_sno = param_set["ve_sno"]
    χm_sno = param_set["χm_sno"]
    Δm_sno = param_set["Δm_sno"]
    χa_sno = param_set["χa_sno"]
    Δa_sno = param_set["Δa_sno"]
    χv_sno = param_set["χv_sno"]
    Δv_sno = param_set["Δv_sno"]
    
    E_liq_rai = param_set["E_liq_rai"]
    E_liq_sno = param_set["E_liq_sno"]
    E_ice_rai = param_set["E_ice_rai"]
    E_ice_sno = param_set["E_ice_sno"]
    E_rai_sno = param_set["E_rai_sno"]

    ρ_cloud_liq = param_set["ρ_cloud_liq"]
    ρ_cloud_ice = param_set["ρ_cloud_ice"]
    grav = param_set["grav"]
    T_freeze = param_set["T_freeze"]
    
    #derived parameters
    N_Sc = ν_air / D_vapor
    m0_ice = 4/3 * π * ρ_cloud_ice * r0_ice ^ me_ice
    m0_rai = 4/3 * π * ρ_cloud_liq * r0_rai ^ me_rai
    a0_rai = π * r0_rai ^ ae_rai
    m0_sno = 1e-1 * r0_sno ^ me_sno
    a0_sno = 0.3 * π * r0_sno ^ ae_sno
    v0_sno = 2^(9/4) * r0_sno ^ ve_sno

    
    return Microphysics_1M_Parameters{valtype(param_set)}(
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
        N_Sc,
        m0_ice,
        m0_rai,
        a0_rai,       
        m0_sno,
        a0_sno,
        v0_sno,
    )
    
    
end


#General parameters outside of modular
struct CloudMicrophysicsParameters{FT, AMPS <: AbstractMicrophysicsParameters}
    K_therm::FT
    D_vapor::FT
    molmass_water::FT
    gas_constant::FT
    ρ_cloud_liq::FT
    surface_tension_coeff::FT    
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
    param_set::Dict,
    MPS::AMPS,
    TPS::ThermodynamicsParameters{FT}) where {FT, AMPS <: AbstractMicrophysicsParameters}
    
    K_therm = param_set["K_therm"]
    D_vapor = param_set["D_vapor"]
    molmass_water = param_set["molmass_water"]
    gas_constant = param_set["gas_constant"]
    ρ_cloud_liq = param_set["ρ_cloud_liq"]
    surface_tension_coeff = param_set["surface_tension_coeff"]
    
    #derived parameters (one could also get this from thermodynamics)
    R_v = gas_constant / molmass_ratio
    molmass_ratio = molmass_dryair / molmass_water

    #derived_parameters
    
    return CloudMicrophysicsParameters{valtype(param_set), AMPS}(
        K_therm,
        D_vapor,
        molmass_water,
        gas_constant,
        ρ_cloud_liq,
        surface_tension_coeff,
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
