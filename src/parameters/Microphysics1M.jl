export CloudLiquid, CloudIce, Rain, Snow, VarTimescaleAcnv

"""
    ParticlePDFSnow{FT}

A struct with snow size distribution parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct ParticlePDFSnow{FT} <: ParametersType
    "snow size distribution coefficient [1/m4]"
    μ::FT
    "snow size distribution coefficient [-]"
    ν::FT
end

"""
    ParticlePDFIceRain{FT}

A struct with snow size distribution parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ParticlePDFIceRain{FT} <: ParametersType
    "Size distribution coefficient [1/m4]"
    n0::FT
end

"""
    ParticleMass{FT}

A struct with coefficients of the assumed mass(size) relationship for particles

m(r) = m0 χm (r/r0)^(me + Δm)

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct ParticleMass{FT} <: ParametersType
    "particle length scale [m]"
    r0::FT
    "mass size relation coefficient [kg]"
    m0::FT
    "mass size relation coefficient [-]"
    me::FT
    "mass size relation coefficient [-]"
    Δm::FT
    "mass size relation coefficient [-]"
    χm::FT
    "pre-computed gamma(me + Δm + 1) for performance [-]"
    gamma_coeff::FT
end
ShowMethods.field_units(::ParticleMass) = (; r0 = "m", m0 = "kg")

"""
    ParticleArea{FT}

A struct with coefficients of the assumed cross_section_area(size)
relationship for particles

a(r) = a0 χa (r/r0)^(ae + Δa)

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct ParticleArea{FT} <: ParametersType
    "cross section size relation coefficient [m²]"
    a0::FT
    "cross section size relation coefficient [-]"
    ae::FT
    "cross section size relation coefficient [-]"
    Δa::FT
    "cross section size relation coefficient [-]"
    χa::FT
end
ShowMethods.field_units(::ParticleArea) = (; a0 = "m²")

"""
    Ventilation{FT}

A struct with ventilation coefficients

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Ventilation{FT} <: ParametersType
    "ventilation coefficient `a` [-]"
    a::FT
    "ventilation coefficient `b` [-]"
    b::FT
end

"""
    SnowAspectRatio{FT}

A struct with aspect ratio coefficients

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SnowAspectRatio{FT} <: ParametersType
    "aspect ratio [-]"
    ϕ::FT
    "power law coefficient in terminal velocity parameterization from Chen et al 2022 [-]"
    κ::FT
end

"""
    Acnv1M{FT}

A struct with autoconversion parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Acnv1M{FT} <: ParametersType
    "autoconversion timescale [s]"
    τ::FT
    "condensate specific content autoconversion threshold [-]"
    q_threshold::FT
    "threshold smooth transition steepness [-]"
    k::FT
end

"""
    CloudLiquid{FT}

The parameters and type for cloud liquid water condensate

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct CloudLiquid{FT} <: CloudCondensateType
    "water density [kg/m³]"
    ρw::FT
    "effective radius [m]"
    r_eff::FT
    "assumed number concentration for cloud sedimentation [1/m³]"
    N_0::FT
end
ShowMethods.field_units(::CloudLiquid) =
    (; ρw = "kg/m³", r_eff = "m", N_0 = "1/m³")

function CloudLiquid(toml_dict::CP.ParamDict)
    name_map = (;
        :density_liquid_water => :ρw,
        :liquid_cloud_effective_radius => :r_eff,
        :cloud_liquid_sedimentation_number_concentration => :N_0,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    return CloudLiquid(; parameters...)
end

"""
    CloudIce{FT, PD, MS}

The parameters and type for cloud ice condensate

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct CloudIce{FT, PD, MS} <: CloudCondensateType
    "a struct with size distribution parameters"
    pdf::PD
    "a struct with mass size relation parameters"
    mass::MS
    "cloud ice apparent density [kg/m³]"
    ρᵢ::FT
    "effective radius [m]"
    r_eff::FT
    "assumed number concentration for cloud sedimentation [1/m³]"
    N_0::FT
end
ShowMethods.field_units(::CloudIce) =
    (; ρᵢ = "kg/m³", r_eff = "m", N_0 = "1/m³")

function CloudIce(toml_dict::CP.ParamDict)
    name_map = (;
        :cloud_ice_apparent_density => :ρᵢ,
        :cloud_ice_size_distribution_coefficient_n0 => :n0,
        :ice_cloud_effective_radius => :r_eff,
        :cloud_ice_sedimentation_number_concentration => :N_0,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    mass = ParticleMass(CloudIce, toml_dict)
    pdf = ParticlePDFIceRain(p.n0)
    return CloudIce(; pdf, mass, p.ρᵢ, p.r_eff, p.N_0)
end

function ParticleMass(::Type{CloudIce}, td::CP.ParamDict)
    name_map = (;
        :cloud_ice_apparent_density => :ρᵢ,
        :cloud_ice_crystals_length_scale => :r0,
        :cloud_ice_mass_size_relation_coefficient_me => :me,
        :cloud_ice_mass_size_relation_coefficient_delm => :Δm,
        :cloud_ice_mass_size_relation_coefficient_chim => :χm,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    m0 = p.ρᵢ * p.r0^p.me * π * 4 / 3
    gamma_coeff = SF.gamma(p.me + p.Δm + 1)
    return ParticleMass(; p.r0, m0, p.me, p.Δm, p.χm, gamma_coeff)
end

"""
    Rain{FT, PD, MS, AR, VT}

The parameters and type for rain

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Rain{PD, MS, AR, VT} <: PrecipitationType
    "a struct with size distribution parameters"
    pdf::PD
    "a struct with mass size relation parameters"
    mass::MS
    "a struct with cross section size relation parameters"
    area::AR
    "a struct with ventilation coefficients"
    vent::VT
end

function Rain(toml_dict::CP.ParamDict)
    name_map = (;
        :rain_drop_size_distribution_coefficient_n0 => :n0,
        :rain_ventilation_coefficient_a => :a,
        :rain_ventilation_coefficient_b => :b,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    return Rain(;
        mass = ParticleMass(Rain, toml_dict),
        area = ParticleArea(Rain, toml_dict),
        pdf = ParticlePDFIceRain(p.n0),
        vent = Ventilation(p.a, p.b),
    )
end

function ParticleMass(::Type{Rain}, td::CP.ParamDict)
    name_map = (;
        :density_liquid_water => :ρ,
        :rain_drop_length_scale => :r0,
        :rain_mass_size_relation_coefficient_me => :me,
        :rain_mass_size_relation_coefficient_delm => :Δm,
        :rain_mass_size_relation_coefficient_chim => :χm,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    m0 = p.ρ * p.r0^p.me * π * 4 / 3
    gamma_coeff = SF.gamma(p.me + p.Δm + 1)
    return ParticleMass(; p.r0, m0, p.me, p.Δm, p.χm, gamma_coeff)
end

function ParticleArea(::Type{Rain}, td::CP.ParamDict)
    name_map = (;
        :rain_drop_length_scale => :r0,
        :rain_cross_section_size_relation_coefficient_ae => :ae,
        :rain_cross_section_size_relation_coefficient_dela => :Δa,
        :rain_cross_section_size_relation_coefficient_chia => :χa,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    a0 = π * p.r0^p.ae
    return ParticleArea(; a0, p.ae, p.Δa, p.χa)
end

"""
    Snow{FT, PD, MS, AR, VT, AP}

The parameters and type for snow

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Snow{FT, PD, MS, AR, VT, AP} <: PrecipitationType
    "a struct with size distribution parameters"
    pdf::PD
    "a struct with mass size relation parameters"
    mass::MS
    "a struct with cross section size relation parameters"
    area::AR
    "a struct with ventilation coefficients"
    vent::VT
    "a struct with aspect ratio parameters"
    aspr::AP
    "snow apparent density [kg/m3]"
    ρᵢ::FT
    "pre-computed gamma(α+4)/6 for oblate aspect ratio [-]"
    gamma_aspect_oblate::FT
    "pre-computed gamma(α+4)/6 for prolate aspect ratio [-]"
    gamma_aspect_prolate::FT
end

function Snow(toml_dict::CP.ParamDict)
    name_map = (;
        :snow_apparent_density => :ρᵢ,
        :snow_flake_size_distribution_coefficient_mu => :μ,
        :snow_flake_size_distribution_coefficient_nu => :ν,
        :snow_ventilation_coefficient_a => :a,
        :snow_ventilation_coefficient_b => :b,
        :snow_aspect_ratio => :ϕ,
        :snow_aspect_ratio_coefficient => :κ,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")

    mass = ParticleMass(Snow, toml_dict)
    area = ParticleArea(Snow, toml_dict)
    FT = CP.float_type(toml_dict)

    # Pre-compute gamma aspect ratio for oblate and prolate shapes
    # Oblate: α = me + Δm - 3/2 * (ae + Δa)
    # Prolate: α = 3 * (ae + Δa) - 2 * (me + Δm)
    α_oblate = mass.me + mass.Δm - (3 // 2) * (area.ae + area.Δa)
    α_prolate = 3 * (area.ae + area.Δa) - 2 * (mass.me + mass.Δm)

    return Snow(;
        pdf = ParticlePDFSnow(p.μ, p.ν),
        mass,
        area,
        vent = Ventilation(p.a, p.b),
        aspr = SnowAspectRatio(p.ϕ, p.κ),
        p.ρᵢ,
        gamma_aspect_oblate = SF.gamma(α_oblate + 4) / SF.gamma(FT(4)),
        gamma_aspect_prolate = SF.gamma(α_prolate + 4) / SF.gamma(FT(4)),
    )
end

function ParticleMass(::Type{Snow}, td::CP.ParamDict)
    name_map = (;
        :snow_flake_length_scale => :r0,
        :snow_mass_size_relation_coefficient_me => :me,
        :snow_mass_size_relation_coefficient_delm => :Δm,
        :snow_mass_size_relation_coefficient_chim => :χm,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    m0 = p.r0^p.me / 10
    gamma_coeff = SF.gamma(p.me + p.Δm + 1)
    return ParticleMass(; p.r0, m0, p.me, p.Δm, p.χm, gamma_coeff)
end

function ParticleArea(::Type{Snow}, toml_dict::CP.ParamDict)
    name_map = (;
        :snow_flake_length_scale => :r0,
        :snow_cross_section_size_relation_coefficient => :ae,
        :snow_cross_section_size_relation_coefficient_dela => :Δa,
        :snow_cross_section_size_relation_coefficient_chia => :χa,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    FT = CP.float_type(toml_dict)
    a0 = FT(0.3 * π * p.r0^p.ae)
    return ParticleArea(; a0, p.ae, p.Δa, p.χa)
end



"""
    VarTimescaleAcnv{FT}

Parameters for the variable-timescale rain autoconversion scheme used in the 1-moment
microphysics scheme, following Azimi et al. (2023).
Active when `RainAutoconversionPrescribedNd` is selected in `Microphysics1MOptions`.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct VarTimescaleAcnv{FT} <: Precipitation2MType
    "Autoconversion timescale [s]"
    τ::FT
    "Powerlaw coefficient relating timescale to droplet number [-]"
    α::FT
    "Prescribed cloud droplet number concentration [1/m³]"
    Nc::FT
end

function VarTimescaleAcnv(td::CP.ParamDict)
    name_map = (;
        :rain_autoconversion_timescale => :τ,
        :Variable_time_scale_autoconversion_coeff_alpha => :α,
        :prescribed_cloud_droplet_number_concentration => :Nc,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return VarTimescaleAcnv(; parameters...)
end
