import CloudMicrophysics.Parameters as CMP
import Thermodynamics.Parameters as TDP

struct Empty{FT} <: CMP.ParametersType{FT} end

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
end

struct ABDINM{FT} <: CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    aerosol::CMP.AerosolType{FT}
    r_nuc::FT
end

struct P3_dep{FT} <: CMP.ParametersType{FT}
    ips::CMP.ParametersType{FT}
    const_dt::FT
end

struct ABIFM{FT} <: CMP.ParametersType{FT}
    H₂SO₄ps::CMP.H2SO4SolutionParameters{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    aerosol::CMP.AerosolType{FT}
    A_aer::FT
end

struct P3_het{FT} <: CMP.ParametersType{FT}
    ips::CMP.ParametersType{FT}
    const_dt::FT
end

struct Frostenberg_random{FT} <: CMP.ParametersType{FT}
    ip::CMP.ParametersType{FT}
    drawing_interval::FT
end

struct Frostenberg_stochastic{FT} <: CMP.ParametersType{FT}
    ip::CMP.ParametersType{FT}
    γ::FT
end

struct Frostenberg_mean{FT} <: CMP.ParametersType{FT}
    ip::CMP.ParametersType{FT}
end

struct ABHOM{FT} <: CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    ips::CMP.ParametersType{FT}
end

struct P3_hom{FT} <: CMP.ParametersType{FT}
    const_dt::FT
end

struct CondParams{FT} <: CMP.ParametersType{FT}
    aps::CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
end

struct DepParams{FT} <: CMP.ParametersType{FT}
    aps::CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
end
