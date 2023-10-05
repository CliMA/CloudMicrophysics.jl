"""

Additional type hierarchy to dispatch over for some microphysics parameters

"""
module CommonTypes

using DocStringExtensions
import CLIMAParameters as CP

export AbstractCloudType
export AbstractPrecipType
export LiquidType
export IceType
export RainType
export SnowType
export AbstractAerosolDistribution
export Abstract2MPrecipType
export KK2000Type
export B1994Type
export TC1980Type
export LD2004Type
export SB2006Type
export AbstractTerminalVelocityType
export Blk1MVelType
export SB2006VelType
export Chen2022Type
export AbstractAerosolType
export ArizonaTestDustType
export DesertDustType
export KaoliniteType
export IlliteType
export CollisionEfficiency

"""
    AbstractAerosolDistribution

The top-level super-type for all aerosol distribution types.
"""
abstract type AbstractAerosolDistribution{T} end
Base.broadcastable(x::AbstractAerosolDistribution) = tuple(x)

"""
    AbstractCloudType

The top-level super-type for cloud liquid water and cloud ice types
"""
abstract type AbstractCloudType end
Base.broadcastable(x::AbstractCloudType) = tuple(x)

"""
    AbstractPrecipType

The top-level super-type for precipitation types (rain and snow)
"""
abstract type AbstractPrecipType end
Base.broadcastable(x::AbstractPrecipType) = tuple(x)

"""
    LiquidType

The type for cloud liquid water condensate

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LiquidType{FT} <: AbstractCloudType
    "condensation evaporation non_equil microphysics relaxation timescale [s]"
    τ_relax::FT
end

function LiquidType(::Type{FT}, toml_dict = CP.create_toml_dict(FT)) where {FT}
    (; data) = toml_dict
    τ_relax = FT(data["condensation_evaporation_timescale"]["value"])
    return LiquidType(τ_relax)
end
Base.broadcastable(x::LiquidType) = tuple(x)

"""
    IceType

The type for cloud ice condensate

# Fields
$(DocStringExtensions.FIELDS)
"""
struct IceType{FT} <: AbstractCloudType
    "cloud ice crystals length scale"
    r0::FT
    "cloud ice mass size relation coefficient me"
    me::FT
    "cloud ice mass size relation coefficient delm"
    Δm::FT
    "cloud ice mass_size_relation_coefficient chim"
    χm::FT
    "density ice water"
    ρ::FT
    "cloud ice size distribution coefficient n0"
    n0::FT
    "m0"
    m0::FT
    "ice snow threshold radius"
    r_ice_snow::FT
    "deposition sublimation non_equil microphysics relaxation timescale [s]"
    τ_relax::FT
end

function IceType(::Type{FT}, toml_dict = CP.create_toml_dict(FT)) where {FT}
    (; data) = toml_dict
    r0 = FT(data["cloud_ice_crystals_length_scale"]["value"])
    me = FT(data["cloud_ice_mass_size_relation_coefficient_me"]["value"])
    Δm = FT(data["cloud_ice_mass_size_relation_coefficient_delm"]["value"])
    χm = FT(data["cloud_ice_mass_size_relation_coefficient_chim"]["value"])
    ρ = FT(data["density_ice_water"]["value"])
    n0 = FT(data["cloud_ice_size_distribution_coefficient_n0"]["value"])
    m0 = FT(4 / 3) * π * ρ * r0^me
    r_ice_snow = FT(data["ice_snow_threshold_radius"]["value"])
    τ_relax = FT(data["sublimation_deposition_timescale"]["value"])
    return IceType(r0, me, Δm, χm, ρ, n0, m0, r_ice_snow, τ_relax)
end
Base.broadcastable(x::IceType) = tuple(x)

"""
    RainType

The type for rain

# Fields
$(DocStringExtensions.FIELDS)
"""
struct RainType{FT} <: AbstractPrecipType
    "rain drop length scale"
    r0::FT
    "rain mass size relation coefficient me"
    me::FT
    "rain mass size relation coefficient delm"
    Δm::FT
    "rain mass_size_relation_coefficient chim"
    χm::FT
    "rain cross section size relation coefficient ae"
    ae::FT
    "rain cross section size relation coefficient delta"
    Δa::FT
    "rain cross section size relation coefficient chia"
    χa::FT
    "rain drop size distribution coefficient n0"
    n0::FT
    "density_cloud_liquid"
    ρ_cloud_liq::FT
    "rain drop drag coefficient"
    C_drag::FT
    "grav"
    grav::FT
    "m0"
    m0::FT
    "a0"
    a0::FT
    "rain terminal velocity size relation coefficient ve"
    ve::FT
    "rain terminal velocity size relation coefficient delv"
    Δv::FT
    "rain terminal velocity size relation coefficient chiv"
    χv::FT
    "rain ventilation coefficient `a`"
    a_vent::FT
    "rain ventilation coefficient `b`"
    b_vent::FT
end

function RainType(::Type{FT}, toml_dict = CP.create_toml_dict(FT)) where {FT}
    (; data) = toml_dict
    r0 = FT(data["rain_drop_length_scale"]["value"])
    me = FT(data["rain_mass_size_relation_coefficient_me"]["value"])
    Δm = FT(data["rain_mass_size_relation_coefficient_delm"]["value"])
    χm = FT(data["rain_mass_size_relation_coefficient_chim"]["value"])
    ae = FT(data["rain_cross_section_size_relation_coefficient_ae"]["value"])
    Δa = FT(data["rain_cross_section_size_relation_coefficient_dela"]["value"])# typo in .toml file
    χa = FT(data["rain_cross_section_size_relation_coefficient_chia"]["value"])
    n0 = FT(data["rain_drop_size_distribution_coefficient_n0"]["value"])
    ρ_cloud_liq = FT(data["density_liquid_water"]["value"])
    grav = FT(data["gravitational_acceleration"]["value"])
    C_drag = FT(data["rain_drop_drag_coefficient"]["value"])
    m0 = FT(4 / 3) * π * ρ_cloud_liq * r0^me
    a0 = FT(π) * r0^ae
    ve =
        FT(data["rain_terminal_velocity_size_relation_coefficient_ve"]["value"])
    Δv = FT(
        data["rain_terminal_velocity_size_relation_coefficient_delv"]["value"],
    )
    χv = FT(
        data["rain_terminal_velocity_size_relation_coefficient_chiv"]["value"],
    )
    a_vent = FT(data["rain_ventillation_coefficient_a"]["value"])
    b_vent = FT(data["rain_ventillation_coefficient_b"]["value"])
    return RainType(
        r0,
        me,
        Δm,
        χm,
        ae,
        Δa,
        χa,
        n0,
        ρ_cloud_liq,
        C_drag,
        grav,
        m0,
        a0,
        ve,
        Δv,
        χv,
        a_vent,
        b_vent,
    )
end
Base.broadcastable(x::RainType) = tuple(x)

"""
    SnowType

The type for snow

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SnowType{FT} <: AbstractPrecipType
    "snow flake length scale"
    r0::FT
    "snow mass size relation coefficient me"
    me::FT
    "snow mass size relation coefficient delm"
    Δm::FT
    "snow mass size relation coefficient chim"
    χm::FT
    "snow cross section size relation coefficient ae"
    ae::FT
    "snow cross section size relation coefficient delta"
    Δa::FT
    "snow cross section size relation coefficient chia"
    χa::FT
    "density ice water"
    ρ::FT
    "snow terminal velocity size relation coefficient"
    ve::FT
    "snow terminal velocity size relation coefficient delv"
    Δv::FT
    "snow terminal velocity size relation coefficient chiv"
    χv::FT
    "m0"
    m0::FT
    "a0"
    a0::FT
    "v0"
    v0::FT
    "snow flake size distribution coefficient μ"
    μ::FT
    "snow flake size distribution coefficient ν"
    ν::FT
    "freezing temperature of water"
    T_freeze::FT
    "snow ventilation coefficient `a`"
    a_vent::FT
    "snow ventilation coefficient `b`"
    b_vent::FT
end

function SnowType(::Type{FT}, toml_dict = CP.create_toml_dict(FT)) where {FT}
    (; data) = toml_dict
    r0 = FT(data["snow_flake_length_scale"]["value"])
    me = FT(data["snow_mass_size_relation_coefficient_me"]["value"])
    Δm = FT(data["snow_mass_size_relation_coefficient_delm"]["value"])
    χm = FT(data["snow_mass_size_relation_coefficient_chim"]["value"])
    ae = FT(data["snow_cross_section_size_relation_coefficient"]["value"])
    Δa = FT(data["snow_cross_section_size_relation_coefficient_dela"]["value"])# typo in .toml file
    χa = FT(data["snow_cross_section_size_relation_coefficient_chia"]["value"])
    ρ = FT(data["density_ice_water"]["value"])
    ve = FT(data["snow_terminal_velocity_size_relation_coefficient"]["value"])
    Δv = FT(
        data["snow_terminal_velocity_size_relation_coefficient_delv"]["value"],
    )
    χv = FT(
        data["snow_terminal_velocity_size_relation_coefficient_chiv"]["value"],
    )
    μ = FT(data["snow_flake_size_distribution_coefficient_mu"]["value"])
    ν = FT(data["snow_flake_size_distribution_coefficient_nu"]["value"])
    m0 = FT(1e-1) * r0^me
    a0 = FT(0.3) * π * r0^ae
    v0 = FT(2^(9 / 4)) * r0^ve
    T_freeze = FT(data["temperature_water_freeze"]["value"])
    a_vent = FT(data["snow_ventillation_coefficient_a"]["value"])
    b_vent = FT(data["snow_ventillation_coefficient_b"]["value"])
    return SnowType(
        r0,
        me,
        Δm,
        χm,
        ae,
        Δa,
        χa,
        ρ,
        ve,
        Δv,
        χv,
        m0,
        a0,
        v0,
        μ,
        ν,
        T_freeze,
        a_vent,
        b_vent,
    )
end
Base.broadcastable(x::SnowType) = tuple(x)

"""
    Abstract2MPrecipType

The top-level super-type for 2-moment precipitation formation types
"""
abstract type Abstract2MPrecipType end
Base.broadcastable(x::Abstract2MPrecipType) = tuple(x)

"""
    KK2000Type

The type for 2-moment precipitation formation by Khairoutdinov and Kogan (2000)
"""
struct KK2000Type <: Abstract2MPrecipType end
Base.broadcastable(x::KK2000Type) = tuple(x)

"""
    B1994Type

The type for 2-moment precipitation formation by Beheng (1994)
"""
struct B1994Type <: Abstract2MPrecipType end
Base.broadcastable(x::B1994Type) = tuple(x)

"""
    TC1980Type

The type for 2-moment precipitation formation by Tripoli and Cotton (1980)
"""
struct TC1980Type <: Abstract2MPrecipType end
Base.broadcastable(x::TC1980Type) = tuple(x)

"""
    LD2004Type

The type for 2-moment precipitation formation by Liu and Daum (2004)
"""
struct LD2004Type <: Abstract2MPrecipType end
Base.broadcastable(x::LD2004Type) = tuple(x)

"""
    VarTimeScaleAcnvType

The type for 2-moment precipitation formation based on the 1-moment parameterization
"""
struct VarTimeScaleAcnvType <: Abstract2MPrecipType end
Base.broadcastable(x::VarTimeScaleAcnvType) = tuple(x)

"""
    SB2006Type

The type for 2-moment precipitation formation by Seifert and Beheng (2006)
"""
struct SB2006Type <: Abstract2MPrecipType end
Base.broadcastable(x::SB2006Type) = tuple(x)


struct AccretionKK2000{FT}
    A::FT
    a::FT
    b::FT
end

function AccretionKK2000(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    A = FT(data["KK2000_accretion_coeff_A"]["value"])
    a = FT(data["KK2000_accretion_coeff_a"]["value"])
    b = FT(data["KK2000_accretion_coeff_b"]["value"])
    return AccretionKK2000(A, a, b)
end

struct AutoconversionKK2000{FT}
    A::FT
    a::FT
    b::FT
    c::FT
end

function AutoconversionKK2000(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    A = FT(data["KK2000_auctoconversion_coeff_A"]["value"])
    a = FT(data["KK2000_auctoconversion_coeff_a"]["value"])
    b = FT(data["KK2000_auctoconversion_coeff_b"]["value"])
    c = FT(data["KK2000_auctoconversion_coeff_c"]["value"])
    return AutoconversionKK2000(A, a, b, c)
end

struct AccretionB1994{FT}
    A::FT
end

function AccretionB1994(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    A = FT(data["B1994_accretion_coeff_A"]["value"])
    return AccretionB1994(A)
end

struct AutoconversionB1994{FT}
    C::FT
    a::FT
    b::FT
    c::FT
    N_0::FT
    d_low::FT
    d_high::FT
    k::FT
end

function AutoconversionB1994(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    C = FT(data["B1994_auctoconversion_coeff_C"]["value"])
    a = FT(data["B1994_auctoconversion_coeff_a"]["value"])
    b = FT(data["B1994_auctoconversion_coeff_b"]["value"])
    c = FT(data["B1994_auctoconversion_coeff_c"]["value"])
    N_0 = FT(data["B1994_auctoconversion_coeff_N_0"]["value"])
    d_low = FT(data["B1994_auctoconversion_coeff_d_low"]["value"]) # typo aucto => auto
    d_high = FT(data["B1994_auctoconversion_coeff_d_high"]["value"]) # type aucto => auto
    k_threshold_steepness =
        FT(data["threshold_smooth_transition_steepness"]["value"])
    return AutoconversionB1994(
        C,
        a,
        b,
        c,
        N_0,
        d_low,
        d_high,
        k_threshold_steepness,
    )
end

struct AutoconversionVarTimescale{FT}
    τ::FT
    α::FT
end

function AutoconversionVarTimescale(
    ::Type{FT};
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    τ = FT(data["rain_autoconversion_timescale"]["value"])
    α = FT(data["Variable_time_scale_autoconversion_coeff_alpha"]["value"])
    return AutoconversionVarTimescale(τ, α)
end

struct AccretionTC1980{FT}
    A::FT
end

function AccretionTC1980(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    A = FT(data["TC1980_accretion_coeff_A"]["value"])
    return AccretionTC1980(A)
end

struct AutoconversionTC1980{FT}
    a::FT
    b::FT
    D::FT
    r_0::FT
    me_liq::FT
    m0_liq_coeff::FT
    k::FT
end

function AutoconversionTC1980(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    a = FT(data["TC1980_autoconversion_coeff_a"]["value"])
    b = FT(data["TC1980_autoconversion_coeff_b"]["value"])
    D = FT(data["TC1980_autoconversion_coeff_D"]["value"])
    r_0 = FT(data["TC1980_autoconversion_coeff_r_0"]["value"])
    me_liq = FT(data["TC1980_autoconversion_coeff_me_liq"]["value"])
    ρ_cloud_liq = FT(data["density_liquid_water"]["value"])
    m0_liq_coeff = FT(4 / 3) * π * ρ_cloud_liq
    k_threshold_steepness =
        FT(data["threshold_smooth_transition_steepness"]["value"])
    return AutoconversionTC1980(
        a,
        b,
        D,
        r_0,
        me_liq,
        m0_liq_coeff,
        k_threshold_steepness,
    )
end

struct AutoconversionLD2004{FT}
    R_6C_0::FT
    E_0::FT
    ρ_w::FT
    k::FT
end

function AutoconversionLD2004(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    R_6C_0 = FT(data["LD2004_R_6C_coeff"]["value"])
    E_0 = FT(data["LD2004_E_0_coeff"]["value"])
    ρ_cloud_liq = FT(data["density_liquid_water"]["value"])
    k_threshold_steepness =
        FT(data["threshold_smooth_transition_steepness"]["value"])
    return AutoconversionLD2004(R_6C_0, E_0, ρ_cloud_liq, k_threshold_steepness)
end

struct AccretionSB2006{FT}
    kcr::FT
    τ0::FT
    ρ0::FT
    c::FT
end

function AccretionSB2006(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    kcr = FT(data["SB2006_collection_kernel_coeff_kcr"]["value"])
    τ0 = FT(data["SB2006_accretion_correcting_function_coeff_tau0"]["value"])
    ρ0 = FT(data["SB2006_reference_air_density"]["value"])
    c = FT(data["SB2006_accretion_correcting_function_coeff_c"]["value"])
    return AccretionSB2006(kcr, τ0, ρ0, c)
end

abstract type Autoconversion end

struct AutoconversionSB2006{FT} <: Autoconversion
    kcc::FT
    νc::FT
    x_star::FT
    ρ0::FT
    A::FT
    a::FT
    b::FT
end

function AutoconversionSB2006(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    kcc = FT(data["SB2006_collection_kernel_coeff_kcc"]["value"])
    νc = FT(data["SB2006_cloud_gamma_distribution_parameter"]["value"]) # shoud this be γ?
    x_star = FT(data["SB2006_raindrops_min_mass"]["value"])
    ρ0 = FT(data["SB2006_reference_air_density"]["value"])
    A = FT(data["SB2006_autoconversion_correcting_function_coeff_A"]["value"])
    a = FT(data["SB2006_autoconversion_correcting_function_coeff_a"]["value"])
    b = FT(data["SB2006_autoconversion_correcting_function_coeff_b"]["value"])
    return AutoconversionSB2006(kcc, νc, x_star, ρ0, A, a, b)
end

struct TerminalVelocitySB2006{FT}
    ρ0::FT
    aR::FT
    bR::FT
    cR::FT
    xr_min::FT
    xr_max::FT
    N0_min::FT
    N0_max::FT
    λ_min::FT
    λ_max::FT
    ρ_cloud_liq::FT
end

function TerminalVelocitySB2006(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    ρ0 = FT(data["SB2006_reference_air_density"]["value"])
    aR = FT(data["SB2006_raindrops_terminal_velocity_coeff_aR"]["value"])
    bR = FT(data["SB2006_raindrops_terminal_velocity_coeff_bR"]["value"])
    cR = FT(data["SB2006_raindrops_terminal_velocity_coeff_cR"]["value"])
    xr_min = FT(data["SB2006_raindrops_min_mass"]["value"])
    xr_max = FT(data["SB2006_raindrops_max_mass"]["value"])
    N0_min =
        FT(data["SB2006_raindrops_size_distribution_coeff_N0_min"]["value"])
    N0_max =
        FT(data["SB2006_raindrops_size_distribution_coeff_N0_max"]["value"])
    λ_min =
        FT(data["SB2006_raindrops_size_distribution_coeff_lambda_min"]["value"])
    λ_max =
        FT(data["SB2006_raindrops_size_distribution_coeff_lambda_max"]["value"])
    ρ_cloud_liq = FT(data["density_liquid_water"]["value"])
    return TerminalVelocitySB2006(
        ρ0,
        aR,
        bR,
        cR,
        xr_min,
        xr_max,
        N0_min,
        N0_max,
        λ_min,
        λ_max,
        ρ_cloud_liq,
    )
end

struct EvaporationSB2006{FT, TV}
    av::FT
    bv::FT
    α::FT
    β::FT
    tv::TV
end

function EvaporationSB2006(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    av = FT(data["SB2006_ventilation_factor_coeff_av"]["value"])
    bv = FT(data["SB2006_ventilation_factor_coeff_bv"]["value"])
    α = FT(data["SB2006_rain_evaportation_coeff_alpha"]["value"])
    β = FT(data["SB2006_rain_evaportation_coeff_beta"]["value"])
    tv = TerminalVelocitySB2006(FT)
    return EvaporationSB2006{FT, typeof(tv)}(av, bv, α, β, tv)
end

struct BreakupSB2006{FT, TV}
    Deq::FT
    Dr_th::FT
    kbr::FT
    κbr::FT
    tv::TV
end

function BreakupSB2006(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    Deq = FT(data["SB2006_raindrops_equlibrium_mean_diameter"]["value"])
    Dr_th =
        FT(data["SB2006_raindrops_breakup_mean_diameter_threshold"]["value"])
    kbr = FT(data["SB2006_raindrops_breakup_coeff_kbr"]["value"])
    κbr = FT(data["SB2006_raindrops_breakup_coeff_kappabr"]["value"])
    tv = TerminalVelocitySB2006(FT)
    return BreakupSB2006{FT, typeof(tv)}(Deq, Dr_th, kbr, κbr, tv)
end

struct SelfCollectionSB2006{FT, TV}
    krr::FT
    κrr::FT
    d::FT
    tv::TV
end

function SelfCollectionSB2006(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    krr = FT(data["SB2006_collection_kernel_coeff_krr"]["value"])
    κrr = FT(data["SB2006_collection_kernel_coeff_kapparr"]["value"])
    d = FT(data["SB2006_raindrops_self-collection_coeff_d"]["value"])
    tv = TerminalVelocitySB2006(FT)
    return SelfCollectionSB2006{FT, typeof(tv)}(krr, κrr, d, tv)
end

"""
    AbstractTerminalVelocityType

The top-level super-type for terminal velocity parameterizations
"""
abstract type AbstractTerminalVelocityType end
Base.broadcastable(x::AbstractTerminalVelocityType) = tuple(x)

"""
    Blk1MVelType

The type for precipitation terminal velocity from the simple 1-moment scheme

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Blk1MVelType{FT} <: AbstractTerminalVelocityType
    "rain/snow terminal velocity size relation coefficient chiv"
    χv::FT
    "rain/snow terminal velocity size relation coefficient ve"
    ve::FT
    "rain/snow terminal velocity size relation coefficient delv"
    Δv::FT
end

function Blk1MVelType(
    ::Type{FT},
    ::RainType,
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict

    χv = FT(
        data["rain_terminal_velocity_size_relation_coefficient_chiv"]["value"],
    )
    ve =
        FT(data["rain_terminal_velocity_size_relation_coefficient_ve"]["value"])
    Δv = FT(
        data["rain_terminal_velocity_size_relation_coefficient_delv"]["value"],
    )
    return Blk1MVelType{FT}(χv, ve, Δv)
end

function Blk1MVelType(
    ::Type{FT},
    ::SnowType,
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    χv = FT(
        data["snow_terminal_velocity_size_relation_coefficient_chiv"]["value"],
    )
    ve = FT(data["snow_terminal_velocity_size_relation_coefficient"]["value"])
    Δv = FT(
        data["snow_terminal_velocity_size_relation_coefficient_delv"]["value"],
    )
    return Blk1MVelType{FT}(χv, ve, Δv)
end

Base.broadcastable(x::Blk1MVelType) = tuple(x)

"""
    SB2006VelType

The type for precipitation terminal velocity from Seifert and Beheng (2006)
"""
struct SB2006VelType <: AbstractTerminalVelocityType end
Base.broadcastable(x::SB2006VelType) = tuple(x)

struct Chen2022TypeSnowIce{FT}
    As::FT
    Bs::FT
    Cs::FT
    Es::FT
    Fs::FT
    Gs::FT
end

function Chen2022TypeSnowIce(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    A = (
        FT(data["Chen2022_table_B3_As_coeff_1"]["value"]),
        FT(data["Chen2022_table_B3_As_coeff_2"]["value"]),
        FT(data["Chen2022_table_B3_As_coeff_3"]["value"]),
    )
    B = (
        FT(data["Chen2022_table_B3_Bs_coeff_1"]["value"]),
        FT(data["Chen2022_table_B3_Bs_coeff_2"]["value"]),
        FT(data["Chen2022_table_B3_Bs_coeff_3"]["value"]),
    )
    C = (
        FT(data["Chen2022_table_B3_Cs_coeff_1"]["value"]),
        FT(data["Chen2022_table_B3_Cs_coeff_2"]["value"]),
        FT(data["Chen2022_table_B3_Cs_coeff_3"]["value"]),
        FT(data["Chen2022_table_B3_Cs_coeff_4"]["value"]),
    )
    E = (
        FT(data["Chen2022_table_B3_Es_coeff_1"]["value"]),
        FT(data["Chen2022_table_B3_Es_coeff_2"]["value"]),
        FT(data["Chen2022_table_B3_Es_coeff_3"]["value"]),
    )
    F = (
        FT(data["Chen2022_table_B3_Fs_coeff_1"]["value"]),
        FT(data["Chen2022_table_B3_Fs_coeff_2"]["value"]),
        FT(data["Chen2022_table_B3_Fs_coeff_3"]["value"]),
    )
    G = (
        FT(data["Chen2022_table_B3_Gs_coeff_1"]["value"]),
        FT(data["Chen2022_table_B3_Gs_coeff_2"]["value"]),
        FT(data["Chen2022_table_B3_Gs_coeff_3"]["value"]),
    )
    ρᵢ = FT(data["density_ice_water"]["value"])
    return Chen2022TypeSnowIce(
        A[1] * (log(ρᵢ))^2 − A[2] * log(ρᵢ) - A[3],
        FT(1) / (B[1] + B[2] * log(ρᵢ) + B[3] / sqrt(ρᵢ)),
        C[1] + C[2] * exp(C[3] * ρᵢ) + C[4] * sqrt(ρᵢ),
        E[1] - E[2] * (log(ρᵢ))^2 + E[3] * sqrt(ρᵢ),
        -exp(F[1] - F[2] * (log(ρᵢ))^2 + F[3] * log(ρᵢ)),
        FT(1) / (G[1] + G[2] / (log(ρᵢ)) - G[3] * log(ρᵢ) / ρᵢ),
    )
end

struct Chen2022TypeRain{FT}
    ρ0::FT
    a1::FT
    a2::FT
    a3::FT
    a3_pow::FT
    b1::FT
    b2::FT
    b3::FT
    b_ρ::FT
    c1::FT
    c2::FT
    c3::FT
end

function Chen2022TypeRain(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    ρ0 = FT(data["Chen2022_table_B1_q_coeff"]["value"])
    a1 = FT(data["Chen2022_table_B1_a1_coeff"]["value"])
    a2 = FT(data["Chen2022_table_B1_a2_coeff"]["value"])
    a3 = FT(data["Chen2022_table_B1_a3_coeff"]["value"])
    a3_pow = FT(data["Chen2022_table_B1_a3_pow_coeff"]["value"])

    b1 = FT(data["Chen2022_table_B1_b1_coeff"]["value"])
    b2 = FT(data["Chen2022_table_B1_b2_coeff"]["value"])
    b3 = FT(data["Chen2022_table_B1_b3_coeff"]["value"])
    b_ρ = FT(data["Chen2022_table_B1_b_rho_coeff"]["value"])

    c1 = FT(data["Chen2022_table_B1_c1_coeff"]["value"])
    c2 = FT(data["Chen2022_table_B1_c2_coeff"]["value"])
    c3 = FT(data["Chen2022_table_B1_c3_coeff"]["value"])
    return Chen2022TypeRain(ρ0, a1, a2, a3, a3_pow, b1, b2, b3, b_ρ, c1, c2, c3)
end

"""
    Chen2022Type

The type for precipitation terminal velocity from Chen et. al. 2022
"""
struct Chen2022Type{FT, R, SI} <: AbstractTerminalVelocityType
    rain::R
    snowice::SI
end

function Chen2022Type(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    rain = Chen2022TypeRain(FT, toml_dict)
    snowice = Chen2022TypeSnowIce(FT, toml_dict)
    return Chen2022Type{FT, typeof(rain), typeof(snowice)}(rain, snowice)
end
Base.broadcastable(x::Chen2022Type) = tuple(x)

"""
    AbstractAerosolType

The top-level super-type for all aerosol types.
"""
abstract type AbstractAerosolType end
Base.broadcastable(x::AbstractAerosolType) = tuple(x)

"""
    ArizonaTestDustType

The type for Arizona test dust for deposition activated fraction
"""
struct ArizonaTestDustType <: AbstractAerosolType end
Base.broadcastable(x::ArizonaTestDustType) = tuple(x)

"""
    DesertDustType

The type for desert dust for deposition activated fraction
"""
struct DesertDustType <: AbstractAerosolType end
Base.broadcastable(x::DesertDustType) = tuple(x)

"""
    KaoliniteType

The type for Kaolinite for ABIFM nucleation rate
"""
struct KaoliniteType <: AbstractAerosolType end
Base.broadcastable(x::KaoliniteType) = tuple(x)

"""
    IlliteType

The type for illite for ABIFM nucleation rate
"""
struct IlliteType <: AbstractAerosolType end
Base.broadcastable(x::IlliteType) = tuple(x)

"""
    CollisionEfficiency{FT}

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CollisionEfficiency{FT}
    "cloud liquid-rain collision efficiency"
    e_liq_rai::FT
    "cloud liquid-snow collision efficiency"
    e_liq_sno::FT
    "cloud ice-rain collision efficiency"
    e_ice_rai::FT
    "cloud ice-snow collision efficiency"
    e_ice_sno::FT
    "rain-snow collision efficiency"
    e_rai_sno::FT
end

function CollisionEfficiency(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    e_liq_rai = FT(data["cloud_liquid_rain_collision_efficiency"]["value"])
    e_liq_sno = FT(data["cloud_liquid_snow_collision_efficiency"]["value"])
    e_ice_rai = FT(data["cloud_ice_rain_collision_efficiency"]["value"])
    e_ice_sno = FT(data["cloud_ice_snow_collision_efficiency"]["value"])
    e_rai_sno = FT(data["rain_snow_collision_efficiency"]["value"])
    return CollisionEfficiency(
        e_liq_rai,
        e_liq_sno,
        e_ice_rai,
        e_ice_sno,
        e_rai_sno,
    )
end

"""
    RainAutoconversion1M{FT}

Rain autoconversion parameters for 1M microphysics.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Autoconversion1M{FT} <: Autoconversion
    "rain/snow autoconversion timescale"
    τ::FT
    "cloud liquid water/ice specific humidity autoconversion threshold"
    q_threshold::FT
    "threshold smooth transition steepness"
    k::FT
end

function Autoconversion1M(
    ::Type{FT},
    ::RainType{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    τ = FT(data["rain_autoconversion_timescale"]["value"])
    q_liq_threshold = FT(
        data["cloud_liquid_water_specific_humidity_autoconversion_threshold"]["value"],
    )
    k_threshold_steepness =
        FT(data["threshold_smooth_transition_steepness"]["value"])
    return Autoconversion1M(τ, q_liq_threshold, k_threshold_steepness)
end

function Autoconversion1M(
    ::Type{FT},
    ::SnowType{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    τ = FT(data["snow_autoconversion_timescale"]["value"])
    q_ice_threshold = FT(
        data["cloud_ice_specific_humidity_autoconversion_threshold"]["value"],
    )
    k_threshold_steepness =
        FT(data["threshold_smooth_transition_steepness"]["value"])
    return Autoconversion1M(τ, q_ice_threshold, k_threshold_steepness)
end

"""
    AirProperties{FT}

Air properties.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AirProperties{FT}
    "thermal conductivity of air"
    K_therm::FT
    "diffusivity of water vapor"
    D_vapor::FT
    "kinematic viscosity of air"
    ν_air::FT
end

function AirProperties(
    ::Type{FT},
    toml_dict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    K_therm = FT(data["thermal_conductivity_of_air"]["value"])
    D_vapor = FT(data["diffusivity_of_water_vapor"]["value"])
    ν_air = FT(data["kinematic_viscosity_of_air"]["value"])
    return AirProperties(K_therm, D_vapor, ν_air)
end
end #module CommoniTypes.jl
