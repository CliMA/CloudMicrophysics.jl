export KK2000, B1994, TC1980, LD2004, SB2006

"""
    AcnvKK2000

Khairoutdinov and Kogan (2000) autoconversion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AcnvKK2000{FT} <: ParametersType
    "Autoconversion coefficient A"
    A::FT
    "Autoconversion coefficient a"
    a::FT
    "Autoconversion coefficient b"
    b::FT
    "Autoconversion coefficient c"
    c::FT
end

"""
    AccrKK2000

Khairoutdinov and Kogan (2000) accretion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AccrKK2000{FT} <: ParametersType
    "Accretion coefficient A"
    A::FT
    "Accretion coefficient a"
    a::FT
    "Accretion coefficient b"
    b::FT
end

"""
    KK2000

The type and parameters for 2-moment precipitation formation by
Khairoutdinov and Kogan (2000)

DOI:10.1175/1520-0493(2000)128<0229:ANCPPI>2.0.CO;2

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct KK2000{AV, AR} <: Precipitation2MType
    "Autoconversion parameters"
    acnv::AV
    "Accretion parameters"
    accr::AR
end

KK2000(toml_dict::CP.ParamDict) =
    KK2000(; acnv = AcnvKK2000(toml_dict), accr = AccrKK2000(toml_dict))

"""
    AcnvB1994

Beheng (1994) autoconversion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AcnvB1994{FT} <: ParametersType
    "Autoconversion coeff C"
    C::FT
    "Autoconversion coeff a"
    a::FT
    "Autoconversion coeff b"
    b::FT
    "Autoconversion coeff c"
    c::FT
    "Autoconversion coeff N_0"
    N_0::FT
    "Autoconversion coeff d_low"
    d_low::FT
    "Autoconversion coeff d_high"
    d_high::FT
    "Threshold for smooth transition steepness"
    k::FT
end

"""
    AccrB1994

Beheng (1994) accretion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AccrB1994{FT} <: ParametersType
    "Accretion coefficient A"
    A::FT
end

"""
    B1994

The type and parameter for 2-moment precipitation formation by Beheng (1994)
DOI: 10.1016/0169-8095(94)90020-5

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct B1994{AV, AR} <: Precipitation2MType
    "Autoconversion parameters"
    acnv::AV
    "Accretion parameters"
    accr::AR
end

B1994(toml_dict::CP.ParamDict) =
    B1994(; acnv = AcnvB1994(toml_dict), accr = AccrB1994(toml_dict))

"""
    AcnvTC1980

Tripoli and Cotton (1980) autoconversion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AcnvTC1980{FT} <: ParametersType
    "Autoconversion coefficient a"
    a::FT
    "Autoconversion coefficient b"
    b::FT
    "Autoconversion coefficient D"
    D::FT
    "Autoconversion coefficient r_0"
    r_0::FT
    "Autoconversion coefficient me_liq"
    me_liq::FT
    "Autoconversion coefficient m0_liq_coeff"
    m0_liq_coeff::FT
    "Threshold for smooth transition steepness"
    k::FT
end

function AcnvTC1980(td::CP.ParamDict)
    parameters = make_params(td, name_map(AcnvTC1980))
    m0_liq_coeff = parameters.m0_liq_coeff * 4 / 3 * π
    return AcnvTC1980(; parameters..., m0_liq_coeff)
end

"""
    AccrTC1980

Tripoli and Cotton (1980) accretion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AccrTC1980{FT} <: ParametersType
    "Accretion coefficient A"
    A::FT
end

"""
    TC1980

The type and parameters for 2-moment precipitation formation by
Tripoli and Cotton (1980)

DOI: 10.1175/1520-0450(1980)019<1037:ANIOSF>2.0.CO;2

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct TC1980{AV, AR} <: Precipitation2MType
    "Autoconversion parameters"
    acnv::AV
    "Accretion parameters"
    accr::AR
end

TC1980(toml_dict::CP.ParamDict) =
    TC1980(; acnv = AcnvTC1980(toml_dict), accr = AccrTC1980(toml_dict))

"""
    LD2004

The type and parameters for 2-moment precipitation formation by
Liu and Daum (2004)

DOI: 10.1175/1520-0469(2004)061<1539:POTAPI>2.0.CO;2

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct LD2004{FT} <: Precipitation2MType
    "Autoconversion coefficient R_6C_0"
    R_6C_0::FT
    "Autoconversion coefficient E_0"
    E_0::FT
    "liquid water density [kg/m3]"
    ρ_w::FT
    "Threshold for smooth transition steepness"
    k::FT
end

"""
    RainParticlePDF_SB2006

Abstract type for the size distribution parameters of rain particles

See [`RainParticlePDF_SB2006_limited`](@ref) and [`RainParticlePDF_SB2006_notlimited`](@ref)
for the concrete types. These can be constructed with:
```julia
RainParticlePDF_SB2006(toml_dict; is_limited = true) # -> RainParticlePDF_SB2006_limited

RainParticlePDF_SB2006(toml_dict; is_limited = false) # -> RainParticlePDF_SB2006_notlimited
```
where `toml_dict` is a `CP.ParamDict` containing the parameters for the size
distribution, and `is_limited` is a boolean indicating whether to use the
limited or not limited version of the size distribution.
"""
abstract type RainParticlePDF_SB2006 <: ParametersType end
RainParticlePDF_SB2006(toml_dict::CP.ParamDict; is_limited = true) =
    is_limited ?
    RainParticlePDF_SB2006_limited(toml_dict) :
    RainParticlePDF_SB2006_notlimited(toml_dict)

Base.show(io::IO, mime::MIME"text/plain", x::RainParticlePDF_SB2006) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)

"""
    RainParticlePDF_SB2006_limited

Rain size distribution parameters from SB2006 including the limiters
on drop maximum mass and the size distribution coefficients N0 and lambda

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct RainParticlePDF_SB2006_limited{FT} <: RainParticlePDF_SB2006
    "Raindrop size distribution coefficient νr"
    νr::FT
    "Raindrop size distribution coefficient μr"
    μr::FT
    "Raindrop minimal mass"
    xr_min::FT
    "Raindrop maximum mass"
    xr_max::FT
    "Raindrop size distribution coefficient N0 min"
    N0_min::FT
    "Raindrop size distribution coefficient N0 max"
    N0_max::FT
    "Raindrop size distribution coefficient lambda min"
    λ_min::FT
    "Raindrop size distribution coefficient lambda max"
    λ_max::FT
    "Cloud liquid water density [kg/m3]"
    ρw::FT
    "Reference air density [kg/m3]"
    ρ0::FT
end

"""
    RainParticlePDF_SB2006_notlimited

Rain size distribution parameters from SB2006 but without the limiters

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct RainParticlePDF_SB2006_notlimited{FT} <: RainParticlePDF_SB2006
    "Raindrop size distribution coefficient νr"
    νr::FT
    "Raindrop size distribution coefficient μr"
    μr::FT
    "Raindrop minimum mass"
    xr_min::FT
    "Raindrop maximum mass"
    xr_max::FT
    "Cloud liquid water density [kg/m3]"
    ρw::FT
    "Reference air density [kg/m3]"
    ρ0::FT
end

islimited(::RainParticlePDF_SB2006_limited) = true
islimited(::RainParticlePDF_SB2006_notlimited) = false

"""
    CloudParticlePDF_SB2006

Cloud droplets size distribution parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct CloudParticlePDF_SB2006{FT} <: ParametersType
    "Cloud droplet size distribution coefficient νc"
    νc::FT
    "Cloud droplet size distribution coefficient μc"
    μc::FT
    "Cloud droplets minimum mass"
    xc_min::FT
    "Cloud droplets maximum mass (equals the SB2006 minimum raindrop mass)"
    xc_max::FT
    "Cloud liquid water density [kg/m3]"
    ρw::FT
    "pre-computed loggamma((νc + 1) / μc)"
    loggamma_z1::FT
    "pre-computed loggamma((νc + 2) / μc)"
    loggamma_z2::FT
end

function CloudParticlePDF_SB2006(td::CP.ParamDict)
    parameters = make_params(td, name_map(CloudParticlePDF_SB2006))
    (; νc, μc) = parameters
    z1 = (νc + 1) / μc
    z2 = (νc + 2) / μc
    loggamma_z1 = SF.loggamma(z1)
    loggamma_z2 = SF.loggamma(z2)
    return CloudParticlePDF_SB2006(; parameters..., loggamma_z1, loggamma_z2)
end

"""
    AcnvSB2006

Autoconversion parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AcnvSB2006{FT} <: ParametersType
    "Collection kernel coefficient"
    kcc::FT
    "Minimum mass of rain droplets"
    x_star::FT
    "Reference air density [kg/m3]"
    ρ0::FT
    "Autoconversion correcting function coeff A"
    A::FT
    "Autoconversion correcting function coeff a"
    a::FT
    "Autoconversion correcting function coeff b"
    b::FT
end

"""
    AccrSB2006

Accretion parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AccrSB2006{FT} <: ParametersType
    "Collection kernel coefficient Kcr"
    kcr::FT
    "Accretion correcting function coefficient τ_0"
    τ0::FT
    "Reference air density [kg/m3]"
    ρ0::FT
    "Accretion correcting function coefficient c"
    c::FT
end

"""
    SelfColSB2006

Rain selfcollection parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct SelfColSB2006{FT} <: ParametersType
    "Collection kernel coefficient krr"
    krr::FT
    "Collection kernel coefficient kappa rr"
    κrr::FT
    "Raindrop self collection coefficient d"
    d::FT
end

"""
    BreakupSB2006

Rain breakup parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct BreakupSB2006{FT} <: ParametersType
    "Raindrop equilibrium mean diameter"
    Deq::FT
    "Raindrop breakup mean diameter threshold"
    Dr_th::FT
    "Raindrops breakup coefficient kbr"
    kbr::FT
    "Raindrops breakup coefficient kappa br"
    κbr::FT
end

"""
    EvaporationSB2006

Rain evaporation parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct EvaporationSB2006{FT} <: ParametersType
    "Ventilation coefficient a"
    av::FT
    "Ventilation coefficient b"
    bv::FT
    "Rain evaporation coefficient α"
    α::FT
    "Rain evaporation coefficient β"
    β::FT
    "Reference air density [kg/m3]"
    ρ0::FT
    "pre-computed av * Γ(2) / 6^(1/3) = av / cbrt(6)"
    a_vent_1::FT
    "pre-computed bv * Γ(5/2 + 3β/2) / 6^(β/2 + 1/2)"
    b_vent_1::FT
    "pre-computed av * 6^(2/3) = av * cbrt(36)"
    a_vent_0_coeff::FT
    "pre-computed bv / 6^(β/2 - 1/2)"
    b_vent_0_coeff::FT
    "pre-computed incomplete-gamma exponent for the 0th ventilation moment: -1/2 + 3β/2"
    β_vent_0::FT
end

function EvaporationSB2006(td::CP.ParamDict)
    parameters = make_params(td, name_map(EvaporationSB2006))
    (; av, bv, β) = parameters
    FT = typeof(av)
    a_vent_1 = av / cbrt(FT(6))
    b_vent_1 = bv * SF.gamma(FT(5 // 2) + FT(3 // 2) * β) / FT(6)^(β / 2 + FT(1 // 2))
    a_vent_0_coeff = av * cbrt(FT(36))
    b_vent_0_coeff = bv / FT(6)^(β / 2 - FT(0.5))
    β_vent_0 = FT(-0.5) + FT(1.5) * β
    return EvaporationSB2006(; parameters..., a_vent_1, b_vent_1, a_vent_0_coeff, b_vent_0_coeff, β_vent_0)
end

"""
    NumberAdjustmentHorn2012

Number concentration adjustment parameter from Horn (2012, DOI: 10.5194/gmd-5-345-2012)

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct NumberAdjustmentHorn2012{FT} <: ParametersType
    "Number concentration adjustment timescale [s]"
    τ::FT
end

"""
    SB2006

The type and parameters for 2-moment precipitation formation by
Seifert and Beheng (2006). The pdf_r type chooses between running with or without
limiters on raindrop size distribution parameters

DOI: 10.1007/s00703-005-0112-4

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct SB2006{PDc, PDr, AV, AR, SC, BR, EV, NA} <: Precipitation2MType
    "Cloud particle size distribution parameters"
    pdf_c::PDc
    "Rain particle size distribution parameters"
    pdf_r::PDr
    "Autoconversion parameters"
    acnv::AV
    "Accretion parameters"
    accr::AR
    "Rain selfcollection parameters"
    self::SC
    "Rain breakup parameters"
    brek::BR
    "Rain evaporation parameters"
    evap::EV
    "Number concentration adjustment parameter"
    numadj::NA
end

# Construct SB2006 from a ClimaParams TOML dict
SB2006(toml_dict::CP.ParamDict; is_limited = true) =
    SB2006(;
        pdf_c = CloudParticlePDF_SB2006(toml_dict),
        pdf_r = RainParticlePDF_SB2006(toml_dict; is_limited),
        acnv = AcnvSB2006(toml_dict),
        accr = AccrSB2006(toml_dict),
        self = SelfColSB2006(toml_dict),
        brek = BreakupSB2006(toml_dict),
        evap = EvaporationSB2006(toml_dict),
        numadj = NumberAdjustmentHorn2012(toml_dict),
    )

"""
    CondEvap2M{FT}

Relaxation parameters for 2-moment vapor condensation / liquid evaporation

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct CondEvap2M{FT} <: ParametersType
    "condensation/evaporation relaxation timescale [s]"
    τ_relax::FT
end

"""
    SubDep2M{FT}

Relaxation parameters for 2-moment vapor deposition / ice sublimation

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct SubDep2M{FT} <: ParametersType
    "deposition/sublimation relaxation timescale [s]"
    τ_relax::FT
end
