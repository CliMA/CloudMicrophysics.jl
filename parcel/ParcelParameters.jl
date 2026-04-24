import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface.TD.Parameters as TDP

struct Empty <: CMP.ParametersType end

struct AeroAct{FT} <: CMP.ParametersType
    aap::CMP.ParametersType
    aerosol::CMP.AerosolType
    aero_σ_g::FT
    r_nuc::FT
    const_dt::FT
    Nₐ::FT
end

struct MohlerAF{FT} <: CMP.ParametersType
    ips::CMP.ParametersType
    aerosol::CMP.AerosolType
    tps::TDP.ThermodynamicsParameters{FT}
    const_dt::FT
end

struct MohlerRate{FT} <: CMP.ParametersType
    ips::CMP.ParametersType
    aerosol::CMP.AerosolType
    tps::TDP.ThermodynamicsParameters{FT}
    const_dt::FT
end

struct ABDINM{FT} <: CMP.ParametersType
    tps::TDP.ThermodynamicsParameters{FT}
    aerosol::CMP.AerosolType
    r_nuc::FT
    const_dt::FT
end

struct P3_dep{FT} <: CMP.ParametersType
    ips::CMP.ParametersType
    const_dt::FT
end

struct ABIFM{FT} <: CMP.ParametersType
    tps::TDP.ThermodynamicsParameters{FT}
    aerosol::CMP.AerosolType
    A_aer::FT
    const_dt::FT
end

struct P3_het{FT} <: CMP.ParametersType
    ips::CMP.ParametersType
    const_dt::FT
end

struct Frostenberg_random{FT} <: CMP.ParametersType
    ip::CMP.ParametersType
    sampling_interval::FT
    const_dt::FT
end

struct Frostenberg_stochastic{FT} <: CMP.ParametersType
    ip::CMP.ParametersType
    γ::FT
    const_dt::FT
end

struct Frostenberg_mean{FT} <: CMP.ParametersType
    ip::CMP.ParametersType
    const_dt::FT
end

struct ABHOM{FT} <: CMP.ParametersType
    tps::TDP.ThermodynamicsParameters{FT}
    ips::CMP.ParametersType
    const_dt::FT
end

struct P3_hom{FT} <: CMP.ParametersType
    const_dt::FT
end

struct CondParams{FT} <: CMP.ParametersType
    aps::CMP.ParametersType
    tps::TDP.ThermodynamicsParameters{FT}
    const_dt::FT
end

struct NonEqCondParams_simple{FT} <: CMP.ParametersType
    tps::TDP.ThermodynamicsParameters{FT}
    liquid::CMP.CloudLiquid{FT}
end

struct NonEqCondParams{FT} <: CMP.ParametersType
    tps::TDP.ThermodynamicsParameters{FT}
    liquid::CMP.CloudLiquid{FT}
    dt::FT
end

struct DepParams{FT} <: CMP.ParametersType
    aps::CMP.ParametersType
    tps::TDP.ThermodynamicsParameters{FT}
    const_dt::FT
end

struct NonEqDepParams_simple{FT} <: CMP.ParametersType
    tps::TDP.ThermodynamicsParameters{FT}
    ice::CMP.CloudIce{FT}
end

struct NonEqDepParams{FT} <: CMP.ParametersType
    tps::TDP.ThermodynamicsParameters{FT}
    ice::CMP.CloudIce{FT}
    dt::FT
end

struct NonEqDepFrostenbergParams{FT} <: CMP.ParametersType{FT}
    aps::CMP.AirProperties{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    ip::CMP.Frostenberg2023{FT}
    ice::CMP.CloudIce{FT}
    dt::FT
end
