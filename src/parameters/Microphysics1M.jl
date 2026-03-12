export CloudLiquid, CloudIce, Rain, Snow, CollisionEff

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

"""
    ParticleArea{FT}

A struct with coefficients of the assumed cross_section_area(size)
relationship for particles

a(r) = a0 χa (r/r0)^(ae + Δa)

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct ParticleArea{FT} <: ParametersType
    "cross section size relation coefficient [m2]"
    a0::FT
    "cross section size relation coefficient [-]"
    ae::FT
    "cross section size relation coefficient [-]"
    Δa::FT
    "cross section size relation coefficient [-]"
    χa::FT
end

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
    "condensation evaporation non_equil microphysics relaxation timescale [s]"
    τ_relax::FT
    "water density [kg/m3]"
    ρw::FT
    "effective radius [m]"
    r_eff::FT
    "assumed number concentration for cloud sedimentation [1/m3]"
    N_0::FT
end

function CloudLiquid(toml_dict::CP.ParamDict)
    name_map = (;
        :condensation_evaporation_timescale => :τ_relax,
        :density_liquid_water => :ρw,
        :liquid_cloud_effective_radius => :r_eff,
        :cloud_liquid_sedimentation_number_concentration => :N_0,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    return CloudLiquid(; parameters...)
end

"""
    CloudIce{FT, MS}

The parameters and type for cloud ice condensate

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct CloudIce{FT, PD, MS} <: CloudCondensateType
    "a struct with size distribution parameters"
    pdf::PD
    "a struct with mass size relation parameters"
    mass::MS
    "particle length scale [m]"
    r0::FT
    "ice snow threshold radius [m]"
    r_ice_snow::FT
    "deposition sublimation non_equil microphysics relaxation timescale [s]"
    τ_relax::FT
    "cloud ice apparent density [kg/m3]"
    ρᵢ::FT
    "effective radius [m]"
    r_eff::FT
    "assumed number concentration for cloud sedimentation [1/m3]"
    N_0::FT
end

function CloudIce(toml_dict::CP.ParamDict)
    name_map = (;
        :cloud_ice_apparent_density => :ρᵢ,
        :cloud_ice_crystals_length_scale => :r0,
        :cloud_ice_size_distribution_coefficient_n0 => :n0,
        :cloud_ice_mass_size_relation_coefficient_me => :me,
        :ice_snow_threshold_radius => :r_ice_snow,
        :sublimation_deposition_timescale => :τ_relax,
        :ice_cloud_effective_radius => :r_eff,
        :cloud_ice_sedimentation_number_concentration => :N_0,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    mass = ParticleMass(CloudIce, toml_dict)
    pdf = ParticlePDFIceRain(p.n0)
    return CloudIce(; pdf, mass, p.r0, p.r_ice_snow, p.τ_relax, p.ρᵢ, p.r_eff, p.N_0)
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
    Rain{FT, MS, AR, VT, AC}

The parameters and type for rain

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Rain{FT, PD, MS, AR, VT, AC} <: PrecipitationType
    "a struct with size distribution parameters"
    pdf::PD
    "a struct with mass size relation parameters"
    mass::MS
    "a struct with cross section size relation parameers"
    area::AR
    "a struct with ventilation coefficients"
    vent::VT
    "a struct with cloud water to rain autoconversion parameters"
    acnv1M::AC
    "particle length scale [m]"
    r0::FT
end

function Rain(toml_dict::CP.ParamDict)
    name_map = (;
        :density_liquid_water => :ρ,
        :rain_drop_length_scale => :r0,
        :rain_drop_size_distribution_coefficient_n0 => :n0,
        :rain_mass_size_relation_coefficient_me => :me,
        :rain_cross_section_size_relation_coefficient_ae => :ae,
        :rain_autoconversion_timescale => :τ,
        :cloud_liquid_water_specific_humidity_autoconversion_threshold => :q_threshold,
        :threshold_smooth_transition_steepness => :k,
        :rain_ventilation_coefficient_a => :a,
        :rain_ventilation_coefficient_b => :b,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    return Rain(;
        mass = ParticleMass(Rain, toml_dict),
        area = ParticleArea(Rain, toml_dict),
        pdf = ParticlePDFIceRain(p.n0),
        vent = Ventilation(p.a, p.b),
        acnv1M = Acnv1M(p.τ, p.q_threshold, p.k),
        p.r0,
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
    Snow{FT, PD, MS, AR, VT, AP, AC}

The parameters and type for snow

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Snow{FT, PD, MS, AR, VT, AP, AC} <: PrecipitationType
    "a struct with size distribution parameters"
    pdf::PD
    "a struct with mass size relation parameters"
    mass::MS
    "a struct with cross section size relation parameers"
    area::AR
    "a struct with ventilation coefficients"
    vent::VT
    "a struct with aspect ratio parameters"
    aspr::AP
    "a struct with ice to snow autoconversion parameters"
    acnv1M::AC
    "particle length scale [m]"
    r0::FT
    "freezing temperature of water [K]"
    T_freeze::FT
    "snow apparent density [kg/m3]"
    ρᵢ::FT
    "pre-computed gamma(α+4)/6 for oblate aspect ratio [-]"
    gamma_aspect_oblate::FT
    "pre-computed gamma(α+4)/6 for prolate aspect ratio [-]"
    gamma_aspect_prolate::FT
end

function Snow(toml_dict::CP.ParamDict)
    name_map = (;
        :cloud_ice_crystals_length_scale => :r0,
        :snow_apparent_density => :ρᵢ,
        :snow_autoconversion_timescale => :τ,
        :cloud_ice_specific_humidity_autoconversion_threshold => :q_threshold,
        :threshold_smooth_transition_steepness => :k,
        :temperature_water_freeze => :T_freeze,
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
        acnv1M = Acnv1M(p.τ, p.q_threshold, p.k),
        p.r0, p.T_freeze, p.ρᵢ,
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
    CollisionEff{FT}

Collision efficiency parameters for the 1-moment scheme

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct CollisionEff{FT} <: ParametersType
    "cloud liquid-rain collision efficiency [-]"
    e_lcl_rai::FT
    "cloud liquid-snow collision efficiency [-]"
    e_lcl_sno::FT
    "cloud ice-rain collision efficiency [-]"
    e_icl_rai::FT
    "cloud ice-snow collision efficiency [-]"
    e_icl_sno::FT
    "rain-snow collision efficiency [-]"
    e_rai_sno::FT
    "rain-snow velocity dispersion coefficient [-]"
    coeff_disp::FT
end

function CollisionEff(td::CP.ParamDict)
    name_map = (;
        :cloud_liquid_rain_collision_efficiency => :e_lcl_rai,
        :cloud_liquid_snow_collision_efficiency => :e_lcl_sno,
        :cloud_ice_rain_collision_efficiency => :e_icl_rai,
        :cloud_ice_snow_collision_efficiency => :e_icl_sno,
        :rain_snow_collision_efficiency => :e_rai_sno,
        :rain_snow_velocity_dispersion_coefficient => :coeff_disp,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return CollisionEff(; parameters...)
end
