import CloudMicrophysics.Parameters as CMP
import Thermodynamics.Parameters as TDP

struct Empty{FT} <: CMP.ParametersType{FT} end

struct AeroAct{FT} <: CMP.ParametersType{FT}
    aap::CMP.ParametersType{FT}
    aerosol::CMP.AerosolType{FT}
    aero_σ_g::FT
    r_nuc::FT
    const_dt::FT
end

struct MohlerAF{FT} <: CMP.ParametersType{FT}
    ips::CMP.ParametersType{FT}
    aerosol::CMP.AerosolType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    const_dt::FT
end

struct MohlerRate{FT} <: CMP.ParametersType{FT}
    ips::CMP.ParametersType{FT}
    aerosol::CMP.AerosolType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    const_dt::FT
end

struct ABDINM{FT} <: CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    aerosol::CMP.AerosolType{FT}
    r_nuc::FT
    const_dt::FT
end

struct P3_dep{FT} <: CMP.ParametersType{FT}
    ips::CMP.ParametersType{FT}
    const_dt::FT
end

struct ABIFM{FT} <: CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    aerosol::CMP.AerosolType{FT}
    A_aer::FT
    const_dt::FT
end

struct P3_het{FT} <: CMP.ParametersType{FT}
    ips::CMP.ParametersType{FT}
    const_dt::FT
end

struct Frostenberg_random{FT} <: CMP.ParametersType{FT}
    ip::CMP.ParametersType{FT}
    sampling_interval::FT
    const_dt::FT
end

struct Frostenberg_stochastic{FT} <: CMP.ParametersType{FT}
    ip::CMP.ParametersType{FT}
    γ::FT
    const_dt::FT
end

struct Frostenberg_mean{FT} <: CMP.ParametersType{FT}
    ip::CMP.ParametersType{FT}
    const_dt::FT
end

struct ABHOM{FT} <: CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    ips::CMP.ParametersType{FT}
    const_dt::FT
end

struct P3_hom{FT} <: CMP.ParametersType{FT}
    const_dt::FT
end

struct CondParams{FT} <: CMP.ParametersType{FT}
    aps::CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    const_dt::FT
end

struct NonEqCondParams_Anna{FT} <: CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    liquid::CMP.CloudLiquid{FT}
end

struct NonEqCondParams{FT} <: CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    liquid::CMP.CloudLiquid{FT}
    ice::CMP.CloudIce{FT}
end

struct DepParams{FT} <: CMP.ParametersType{FT}
    aps::CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
end

struct NonEqDepParams_Anna{FT} <: CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    ice::CMP.CloudIce{FT}
end

struct NonEqDepParams{FT} <: CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    liquid::CMP.CloudLiquid{FT}
    ice::CMP.CloudIce{FT}
end