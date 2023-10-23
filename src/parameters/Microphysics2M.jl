export KK2000, B1994, TC1980, LD2004, VarTimescaleAcnv, SB2006

"""
    AcnvKK2000

Khairoutdinov and Kogan (2000) autoconversion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AcnvKK2000{FT} <: ParametersType{FT}
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
struct AccrKK2000{FT} <: ParametersType{FT}
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
struct KK2000{FT, AV, AR} <: Precipitation2MType{FT}
    "Autoconversion parameters"
    acnv::AV
    "Accretion parameters"
    accr::AR
end

function KK2000(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    acnv = AcnvKK2000(
        FT(data["KK2000_auctoconversion_coeff_A"]["value"]),
        FT(data["KK2000_auctoconversion_coeff_a"]["value"]),
        FT(data["KK2000_auctoconversion_coeff_b"]["value"]),
        FT(data["KK2000_auctoconversion_coeff_c"]["value"]),
    )
    accr = AccrKK2000(
        FT(data["KK2000_accretion_coeff_A"]["value"]),
        FT(data["KK2000_accretion_coeff_a"]["value"]),
        FT(data["KK2000_accretion_coeff_b"]["value"]),
    )
    return KK2000{FT, typeof(acnv), typeof(accr)}(acnv, accr)
end

"""
    AcnvB1994

Beheng (1994) autoconversion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AcnvB1994{FT} <: ParametersType{FT}
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

function B1994(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict

    acnv = AcnvB1994(
        FT(data["B1994_auctoconversion_coeff_C"]["value"]),
        FT(data["B1994_auctoconversion_coeff_a"]["value"]),
        FT(data["B1994_auctoconversion_coeff_b"]["value"]),
        FT(data["B1994_auctoconversion_coeff_c"]["value"]),
        FT(data["B1994_auctoconversion_coeff_N_0"]["value"]),
        FT(data["B1994_auctoconversion_coeff_d_low"]["value"]), # typo aucto => auto
        FT(data["B1994_auctoconversion_coeff_d_high"]["value"]), # type aucto => auto
        FT(data["threshold_smooth_transition_steepness"]["value"]),
    )
    accr = AccrB1994(FT(data["B1994_accretion_coeff_A"]["value"]))

    return B1994{FT, typeof(acnv), typeof(accr)}(acnv, accr)
end

"""
    AcnvTC1980

Tripoli and Cotton (1980) autoconversion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AcnvTC1980{FT} <: ParametersType{FT}
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

function TC1980(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    ρ_cloud_liq = FT(data["density_liquid_water"]["value"])
    acnv = AcnvTC1980(
        FT(data["TC1980_autoconversion_coeff_a"]["value"]),
        FT(data["TC1980_autoconversion_coeff_b"]["value"]),
        FT(data["TC1980_autoconversion_coeff_D"]["value"]),
        FT(data["TC1980_autoconversion_coeff_r_0"]["value"]),
        FT(data["TC1980_autoconversion_coeff_me_liq"]["value"]),
        FT(4 / 3) * π * ρ_cloud_liq,
        FT(data["threshold_smooth_transition_steepness"]["value"]),
    )
    accr = AccrTC1980(FT(data["TC1980_accretion_coeff_A"]["value"]))

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
struct LD2004{FT} <: Precipitation2MType{FT}
    "Autoconversion coefficient R_6C_0"
    R_6C_0::FT
    "Autoconversion coefficient E_0"
    E_0::FT
    "liquid water density [kg/m3]"
    ρ_w::FT
    "Threshold for smooth tranistion steepness"
    k::FT
end

function LD2004(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return LD2004(
        FT(data["LD2004_R_6C_coeff"]["value"]),
        FT(data["LD2004_E_0_coeff"]["value"]),
        FT(data["density_liquid_water"]["value"]),
        FT(data["threshold_smooth_transition_steepness"]["value"]),
    )
end

"""
    VarTimescaleAcnv

The type for 2-moment precipitation formation based on the
1-moment parameterization with variable time scale Azimi et al (2023)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct VarTimescaleAcnv{FT} <: Precipitation2MType{FT}
    "Timescale [s]"
    τ::FT
    "Powerlaw coefficient [-]"
    α::FT
end

function VarTimescaleAcnv(
    ::Type{FT};
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return VarTimescaleAcnv(
        FT(data["rain_autoconversion_timescale"]["value"]),
        FT(data["Variable_time_scale_autoconversion_coeff_alpha"]["value"]),
    )
end

"""
    ParticlePDF_SB2006

Rain size distribution parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ParticlePDF_SB2006{FT} <: ParametersType{FT}
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
    AcnvSB2006

Autoconversion parameters from SB2006

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AcnvSB2006{FT} <: ParametersType{FT}
    "Collection kernel coefficient"
    kcc::FT
    "Cloud gamma distribution coefficient"
    νc::FT
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
struct AccrSB2006{FT} <: ParametersType{FT}
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
struct SelfColSB2006{FT} <: ParametersType{FT}
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
struct BreakupSB2006{FT} <: ParametersType{FT}
    "Raindrop equilibrium mean diamater"
    Deq::FT
    "Raindrop breakup mean diamater threshold"
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
struct EvaporationSB2006{FT} <: ParametersType{FT}
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

"""
    SB2006

The type and parameters for 2-moment precipitation formation by
Seifert and Beheng (2006)

DOI: 10.1007/s00703-005-0112-4

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SB2006{FT, PD, AV, AR, SC, BR, EV} <: Precipitation2MType{FT}
    "Rain particle size distribution parameters"
    pdf::PD
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
end

function SB2006(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict

    pdf = ParticlePDF_SB2006(
        FT(data["SB2006_raindrops_min_mass"]["value"]),
        FT(data["SB2006_raindrops_max_mass"]["value"]),
        FT(data["SB2006_raindrops_size_distribution_coeff_N0_min"]["value"]),
        FT(data["SB2006_raindrops_size_distribution_coeff_N0_max"]["value"]),
        FT(
            data["SB2006_raindrops_size_distribution_coeff_lambda_min"]["value"],
        ),
        FT(
            data["SB2006_raindrops_size_distribution_coeff_lambda_max"]["value"],
        ),
        FT(data["density_liquid_water"]["value"]),
        FT(data["SB2006_reference_air_density"]["value"]),
    )
    acnv = AcnvSB2006(
        FT(data["SB2006_collection_kernel_coeff_kcc"]["value"]),
        FT(data["SB2006_cloud_gamma_distribution_parameter"]["value"]),
        FT(data["SB2006_raindrops_min_mass"]["value"]),
        FT(data["SB2006_reference_air_density"]["value"]),
        FT(data["SB2006_autoconversion_correcting_function_coeff_A"]["value"]),
        FT(data["SB2006_autoconversion_correcting_function_coeff_a"]["value"]),
        FT(data["SB2006_autoconversion_correcting_function_coeff_b"]["value"]),
    )
    accr = AccrSB2006(
        FT(data["SB2006_collection_kernel_coeff_kcr"]["value"]),
        FT(data["SB2006_accretion_correcting_function_coeff_tau0"]["value"]),
        FT(data["SB2006_reference_air_density"]["value"]),
        FT(data["SB2006_accretion_correcting_function_coeff_c"]["value"]),
    )
    self = SelfColSB2006(
        FT(data["SB2006_collection_kernel_coeff_krr"]["value"]),
        FT(data["SB2006_collection_kernel_coeff_kapparr"]["value"]),
        FT(data["SB2006_raindrops_self-collection_coeff_d"]["value"]),
    )
    brek = BreakupSB2006(
        FT(data["SB2006_raindrops_equlibrium_mean_diameter"]["value"]),
        FT(data["SB2006_raindrops_breakup_mean_diameter_threshold"]["value"]),
        FT(data["SB2006_raindrops_breakup_coeff_kbr"]["value"]),
        FT(data["SB2006_raindrops_breakup_coeff_kappabr"]["value"]),
    )
    evap = EvaporationSB2006(
        FT(data["SB2006_ventilation_factor_coeff_av"]["value"]),
        FT(data["SB2006_ventilation_factor_coeff_bv"]["value"]),
        FT(data["SB2006_rain_evaportation_coeff_alpha"]["value"]),
        FT(data["SB2006_rain_evaportation_coeff_beta"]["value"]),
        FT(data["SB2006_reference_air_density"]["value"]),
    )
    return SB2006{
        FT,
        typeof(pdf),
        typeof(acnv),
        typeof(accr),
        typeof(self),
        typeof(brek),
        typeof(evap),
    }(
        pdf,
        acnv,
        accr,
        self,
        brek,
        evap,
    )
end
