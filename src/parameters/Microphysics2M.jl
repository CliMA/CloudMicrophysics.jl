export KK2000, B1994, TC1980, LD2004, VarTimescaleAcnv, SB2006

"""
    AcnvKK2000

Khairoutdinov and Kogan (2000) autoconversion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AcnvKK2000{FT} <: ParametersType{FT}
    "Autoconversion coefficient A"
    A::FT
    "Autoconversion coefficient a"
    a::FT
    "Autoconversion coefficient b"
    b::FT
    "Autoconversion coefficient c"
    c::FT
end

function AcnvKK2000(td::CP.AbstractTOMLDict)
    name_map = (;
        :KK2000_autoconversion_coeff_A => :A,
        :KK2000_autoconversion_coeff_a => :a,
        :KK2000_autoconversion_coeff_b => :b,
        :KK2000_autoconversion_coeff_c => :c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return AcnvKK2000{FT}(; parameters...)
end

"""
    AccrKK2000

Khairoutdinov and Kogan (2000) accretion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AccrKK2000{FT} <: ParametersType{FT}
    "Accretion coefficient A"
    A::FT
    "Accretion coefficient a"
    a::FT
    "Accretion coefficient b"
    b::FT
end

function AccrKK2000(td::CP.AbstractTOMLDict)
    name_map = (;
        :KK2000_accretion_coeff_A => :A,
        :KK2000_accretion_coeff_a => :a,
        :KK2000_accretion_coeff_b => :b,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return AccrKK2000{FT}(; parameters...)
end

"""
    KK2000

The type and parameters for 2-moment precipitation formation by
Khairoutdinov and Kogan (2000)

DOI:10.1175/1520-0493(2000)128<0229:ANCPPI>2.0.CO;2

# Fields
$(DocStringExtensions.FIELDS)
"""
struct KK2000{FT, AV, AR} <: Precipitation2MType{FT}
    "Autoconversion parameters"
    acnv::AV
    "Accretion parameters"
    accr::AR
end

KK2000(::Type{FT}) where {FT <: AbstractFloat} = KK2000(CP.create_toml_dict(FT))


function KK2000(toml_dict::CP.AbstractTOMLDict)
    acnv = AcnvKK2000(toml_dict)
    accr = AccrKK2000(toml_dict)
    FT = CP.float_type(toml_dict)
    return KK2000{FT, typeof(acnv), typeof(accr)}(acnv, accr)
end

"""
    AcnvB1994

Beheng (1994) autoconversion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AcnvB1994{FT} <: ParametersType{FT}
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
    "Threshold for smooth tranistion steepness"
    k::FT
end

function AcnvB1994(td::CP.AbstractTOMLDict)
    name_map = (;
        :B1994_autoconversion_coeff_C => :C,
        :B1994_autoconversion_coeff_a => :a,
        :B1994_autoconversion_coeff_b => :b,
        :B1994_autoconversion_coeff_c => :c,
        :B1994_autoconversion_coeff_N_0 => :N_0,
        :B1994_autoconversion_coeff_d_low => :d_low,
        :B1994_autoconversion_coeff_d_high => :d_high,
        :threshold_smooth_transition_steepness => :k,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return AcnvB1994{FT}(; parameters...)
end

"""
    AccrB1994

Beheng (1994) accretion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AccrB1994{FT} <: ParametersType{FT}
    "Accretion coefficient A"
    A::FT
end

function AccrB1994(toml_dict::CP.AbstractTOMLDict)
    (; B1994_accretion_coeff_A) = CP.get_parameter_values(
        toml_dict,
        "B1994_accretion_coeff_A",
        "CloudMicrophysics",
    )
    FT = CP.float_type(toml_dict)
    return AccrB1994{FT}(B1994_accretion_coeff_A)
end

"""
    B1994

The type and parameter for 2-moment precipitation formation by Beheng (1994)
DOI: 10.1016/0169-8095(94)90020-5

# Fields
$(DocStringExtensions.FIELDS)
"""
struct B1994{FT, AV, AR} <: Precipitation2MType{FT}
    "Autoconversion coeff C"
    acnv::AV
    "Autoconversion coeff a"
    accr::AR
end

B1994(::Type{FT}) where {FT <: AbstractFloat} = B1994(CP.create_toml_dict(FT))

function B1994(toml_dict::CP.AbstractTOMLDict)
    acnv = AcnvB1994(toml_dict)
    accr = AccrB1994(toml_dict)
    FT = CP.float_type(toml_dict)
    return B1994{FT, typeof(acnv), typeof(accr)}(acnv, accr)
end

"""
    AcnvTC1980

Tripoli and Cotton (1980) autoconversion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AcnvTC1980{FT} <: ParametersType{FT}
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
    "Threshold for smooth tranistion steepness"
    k::FT
end

function AcnvTC1980(td::CP.AbstractTOMLDict)
    name_map = (;
        :TC1980_autoconversion_coeff_a => :a,
        :TC1980_autoconversion_coeff_b => :b,
        :TC1980_autoconversion_coeff_D => :D,
        :TC1980_autoconversion_coeff_r_0 => :r_0,
        :TC1980_autoconversion_coeff_me_liq => :me_liq,
        :threshold_smooth_transition_steepness => :k,
        :density_liquid_water => :m0_liq_coeff,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    m0_liq_coeff = 4 / 3 * π * parameters.m0_liq_coeff

    FT = CP.float_type(td)
    return AcnvTC1980{FT}(; parameters..., m0_liq_coeff)
end

"""
    AccrTC1980

Tripoli and Cotton (1980) accretion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AccrTC1980{FT} <: ParametersType{FT}
    "Accretion coefficient A"
    A::FT
end

function AccrTC1980(toml_dict::CP.AbstractTOMLDict)
    (; TC1980_accretion_coeff_A) = CP.get_parameter_values(
        toml_dict,
        "TC1980_accretion_coeff_A",
        "CloudMicrophysics",
    )
    FT = CP.float_type(toml_dict)
    return AccrTC1980{FT}(TC1980_accretion_coeff_A)
end

"""
    TC1980

The type and parameters for 2-moment precipitation formation by
Tripoli and Cotton (1980)

DOI: 10.1175/1520-0450(1980)019<1037:ANIOSF>2.0.CO;2

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TC1980{FT, AV, AR} <: Precipitation2MType{FT}
    "Autoconversion parameters"
    acnv::AV
    "Accretion parameters"
    accr::AR
end

TC1980(::Type{FT}) where {FT <: AbstractFloat} = TC1980(CP.create_toml_dict(FT))

function TC1980(toml_dict::CP.AbstractTOMLDict)
    acnv = AcnvTC1980(toml_dict)
    accr = AccrTC1980(toml_dict)
    FT = CP.float_type(toml_dict)
    return TC1980{FT, typeof(acnv), typeof(accr)}(acnv, accr)
end

"""
    LD2004

The type and parameters for 2-moment precipitation formation by
Liu and Daum (2004)

DOI: 10.1175/1520-0469(2004)061<1539:POTAPI>2.0.CO;2

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct LD2004{FT} <: Precipitation2MType{FT}
    "Autoconversion coefficient R_6C_0"
    R_6C_0::FT
    "Autoconversion coefficient E_0"
    E_0::FT
    "liquid water density [kg/m3]"
    ρ_w::FT
    "Threshold for smooth tranistion steepness"
    k::FT
end

LD2004(::Type{FT}) where {FT <: AbstractFloat} = LD2004(CP.create_toml_dict(FT))

function LD2004(td::CP.AbstractTOMLDict)
    name_map = (;
        :LD2004_R_6C_coeff => :R_6C_0,
        :LD2004_E_0_coeff => :E_0,
        :density_liquid_water => :ρ_w,
        :threshold_smooth_transition_steepness => :k,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return LD2004{FT}(; parameters...)
end

"""
    VarTimescaleAcnv

The type for 2-moment precipitation formation based on the
1-moment parameterization with variable time scale Azimi et al (2023)

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VarTimescaleAcnv{FT} <: Precipitation2MType{FT}
    "Timescale [s]"
    τ::FT
    "Powerlaw coefficient [-]"
    α::FT
end

VarTimescaleAcnv(::Type{FT}) where {FT <: AbstractFloat} =
    VarTimescaleAcnv(CP.create_toml_dict(FT))

function VarTimescaleAcnv(td::CP.AbstractTOMLDict)
    name_map = (;
        :rain_autoconversion_timescale => :τ,
        :Variable_time_scale_autoconversion_coeff_alpha => :α,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return VarTimescaleAcnv{FT}(; parameters...)
end

"""
    RainParticlePDF_SB2006

Abstract type for the size distribution parameters of rain particles

See [`RainParticlePDF_SB2006_limited`](@ref) and [`RainParticlePDF_SB2006_notlimited`](@ref)
for the concrete types.
"""
abstract type RainParticlePDF_SB2006{FT} <: ParametersType{FT} end

"""
    RainParticlePDF_SB2006_limited

Rain size distribution parameters from SB2006 including the limiters
on drop maximum mass and the size distribution coefficinets N0 and lambda

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct RainParticlePDF_SB2006_limited{FT} <: RainParticlePDF_SB2006{FT}
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

function RainParticlePDF_SB2006_limited(td::CP.AbstractTOMLDict)
    name_map = (;
        :SB2006_rain_distribution_coeff_nu => :νr,
        :SB2006_rain_distribution_coeff_mu => :μr,
        :SB2006_raindrops_min_mass => :xr_min,
        :SB2006_raindrops_max_mass => :xr_max,
        :SB2006_raindrops_size_distribution_coeff_N0_min => :N0_min,
        :SB2006_raindrops_size_distribution_coeff_N0_max => :N0_max,
        :SB2006_raindrops_size_distribution_coeff_lambda_min => :λ_min,
        :SB2006_raindrops_size_distribution_coeff_lambda_max => :λ_max,
        :density_liquid_water => :ρw,
        :SB2006_reference_air_density => :ρ0,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return RainParticlePDF_SB2006_limited{FT}(; parameters...)
end

"""
    RainParticlePDF_SB2006

Rain size distribution parameters from SB2006 but without the limiters

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct RainParticlePDF_SB2006_notlimited{FT} <: RainParticlePDF_SB2006{FT}
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

function RainParticlePDF_SB2006_notlimited(td::CP.AbstractTOMLDict)
    name_map = (;
        :SB2006_rain_distribution_coeff_nu => :νr,
        :SB2006_rain_distribution_coeff_mu => :μr,
        :SB2006_raindrops_min_mass => :xr_min,
        :SB2006_raindrops_max_mass => :xr_max,
        :density_liquid_water => :ρw,
        :SB2006_reference_air_density => :ρ0,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return RainParticlePDF_SB2006_notlimited{FT}(; parameters...)
end

islimited(pdf::RainParticlePDF_SB2006_limited) = true
islimited(pdf::RainParticlePDF_SB2006_notlimited) = false

"""
    CloudParticlePDF_SB2006

Cloud droplets size distribution parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct CloudParticlePDF_SB2006{FT} <: ParametersType{FT}
    "Cloud droplet size distribution coefficient νc"
    νc::FT
    "Cloud droplet size distribution coefficient μc"
    μc::FT
    "Cloud droplets minimum mass"
    xc_min::FT
    "Cloud droplets maximum mass"
    xc_max::FT
    "Cloud liquid water density [kg/m3]"
    ρw::FT
end

function CloudParticlePDF_SB2006(td::CP.AbstractTOMLDict)
    name_map = (;
        :SB2006_cloud_gamma_distribution_parameter => :νc,
        :SB2006_cloud_gamma_distribution_coeff_mu => :μc,
        :SB2006_cloud_droplets_min_mass => :xc_min,
        :SB2006_raindrops_min_mass => :xc_max,
        :density_liquid_water => :ρw,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return CloudParticlePDF_SB2006{FT}(; parameters...)
end

"""
    AcnvSB2006

Autoconversion parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AcnvSB2006{FT} <: ParametersType{FT}
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

function AcnvSB2006(td::CP.AbstractTOMLDict)
    name_map = (;
        :SB2006_collection_kernel_coeff_kcc => :kcc,
        :SB2006_raindrops_min_mass => :x_star,
        :SB2006_reference_air_density => :ρ0,
        :SB2006_autoconversion_correcting_function_coeff_A => :A,
        :SB2006_autoconversion_correcting_function_coeff_a => :a,
        :SB2006_autoconversion_correcting_function_coeff_b => :b,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return AcnvSB2006{FT}(; parameters...)
end


"""
    AccrSB2006

Accretion parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AccrSB2006{FT} <: ParametersType{FT}
    "Collection kernel coefficient Kcr"
    kcr::FT
    "Accretion correcting function coefficient τ_0"
    τ0::FT
    "Reference air density [kg/m3]"
    ρ0::FT
    "Accretion correcting function coefficient c"
    c::FT
end

function AccrSB2006(td::CP.AbstractTOMLDict)
    name_map = (;
        :SB2006_collection_kernel_coeff_kcr => :kcr,
        :SB2006_accretion_correcting_function_coeff_tau0 => :τ0,
        :SB2006_reference_air_density => :ρ0,
        :SB2006_accretion_correcting_function_coeff_c => :c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return AccrSB2006{FT}(; parameters...)
end

"""
    SelfColSB2006

Rain selfcollection parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SelfColSB2006{FT} <: ParametersType{FT}
    "Collection kernel coefficient krr"
    krr::FT
    "Collection kernel coefficient kappa rr"
    κrr::FT
    "Raindrop self collection coefficient d"
    d::FT
end

function SelfColSB2006(td::CP.AbstractTOMLDict)
    name_map = (;
        :SB2006_collection_kernel_coeff_krr => :krr,
        :SB2006_collection_kernel_coeff_kapparr => :κrr,
        Symbol("SB2006_raindrops_self-collection_coeff_d") => :d,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return SelfColSB2006{FT}(; parameters...)
end

"""
    BreakupSB2006

Rain breakup parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct BreakupSB2006{FT} <: ParametersType{FT}
    "Raindrop equilibrium mean diamater"
    Deq::FT
    "Raindrop breakup mean diamater threshold"
    Dr_th::FT
    "Raindrops breakup coefficient kbr"
    kbr::FT
    "Raindrops breakup coefficient kappa br"
    κbr::FT
end

function BreakupSB2006(td::CP.AbstractTOMLDict)
    name_map = (;
        :SB2006_raindrops_equlibrium_mean_diameter => :Deq,
        :SB2006_raindrops_breakup_mean_diameter_threshold => :Dr_th,
        :SB2006_raindrops_breakup_coeff_kbr => :kbr,
        :SB2006_raindrops_breakup_coeff_kappabr => :κbr,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return BreakupSB2006{FT}(; parameters...)
end

"""
    EvaporationSB2006

Rain evaporation parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct EvaporationSB2006{FT} <: ParametersType{FT}
    "Ventillation coefficient a"
    av::FT
    "Ventillation coefficient b"
    bv::FT
    "Rain evapoartion coefficient α"
    α::FT
    "Rain evapoartion coefficient β"
    β::FT
    "Reference air density [kg/m3]"
    ρ0::FT
end

function EvaporationSB2006(td::CP.AbstractTOMLDict)
    name_map = (;
        :SB2006_ventilation_factor_coeff_av => :av,
        :SB2006_ventilation_factor_coeff_bv => :bv,
        :SB2006_rain_evaportation_coeff_alpha => :α,
        :SB2006_rain_evaportation_coeff_beta => :β,
        :SB2006_reference_air_density => :ρ0,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return EvaporationSB2006{FT}(; parameters...)
end

"""
    NumberAdjustmentHorn2012

Number concentration adjustment parameter from Horn (2012, DOI: 10.5194/gmd-5-345-2012)

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct NumberAdjustmentHorn2012{FT} <: ParametersType{FT}
    "Number concentration adjustment timescale [s]"
    τ::FT
end

function NumberAdjustmentHorn2012(td::CP.AbstractTOMLDict)
    name_map = (;
        :Horn2012_number_concentration_adjustment_timescale => :τ,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return NumberAdjustmentHorn2012{FT}(; parameters...)
end

"""
    SB2006

The type and parameters for 2-moment precipitation formation by
Seifert and Beheng (2006). The pdf_r type choses between running with or without
limiters on raindrop size distribution parameters

DOI: 10.1007/s00703-005-0112-4

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SB2006{FT, PDc, PDr, AV, AR, SC, BR, EV, NA} <: Precipitation2MType{FT}
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

SB2006(::Type{FT}, is_limited = true) where {FT <: AbstractFloat} =
    SB2006(CP.create_toml_dict(FT), is_limited)

function SB2006(toml_dict::CP.AbstractTOMLDict, is_limited = true)
    pdf_c = CloudParticlePDF_SB2006(toml_dict)
    pdf_r =
        is_limited ? RainParticlePDF_SB2006_limited(toml_dict) :
        RainParticlePDF_SB2006_notlimited(toml_dict)
    acnv = AcnvSB2006(toml_dict)
    accr = AccrSB2006(toml_dict)
    self = SelfColSB2006(toml_dict)
    brek = BreakupSB2006(toml_dict)
    evap = EvaporationSB2006(toml_dict)
    numadj = NumberAdjustmentHorn2012(toml_dict)
    FT = CP.float_type(toml_dict)
    PDc = typeof(pdf_c)
    PDr = typeof(pdf_r)
    AN = typeof(acnv)
    AR = typeof(accr)
    SE = typeof(self)
    BR = typeof(brek)
    EV = typeof(evap)
    NA = typeof(numadj)
    return SB2006{FT, PDc, PDr, AN, AR, SE, BR, EV, NA}(
        pdf_c,
        pdf_r,
        acnv,
        accr,
        self,
        brek,
        evap,
        numadj,
    )
end
