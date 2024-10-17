export CloudLiquid, CloudIce, Rain, Snow, CollisionEff

"""
    ParticlePDFSnow{FT}

A struct with snow size distribution parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ParticlePDFSnow{FT} <: ParametersType{FT}
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
struct ParticlePDFIceRain{FT} <: ParametersType{FT}
    "snow size distribution coefficient [1/m4]"
    n0::FT
end

"""
    ParticleMass{FT}

A struct with coefficients of the assumed mass(size) relationship for particles

m(r) = m0 χm (r/r0)^(me + Δm)

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ParticleMass{FT} <: ParametersType{FT}
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
end

"""
    ParticleArea{FT}

A struct with coefficients of the assumed cross_section_area(size)
relationship for particles

a(r) = a0 χa (r/r0)^(ae + Δa)

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ParticleArea{FT} <: ParametersType{FT}
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
struct Ventilation{FT} <: ParametersType{FT}
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
struct SnowAspectRatio{FT} <: ParametersType{FT}
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
struct Acnv1M{FT} <: ParametersType{FT}
    "autoconversion timescale [s]"
    τ::FT
    "specific humidity autoconversion threshold [-]"
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
Base.@kwdef struct CloudLiquid{FT} <: CloudCondensateType{FT}
    "condensation evaporation non_equil microphysics relaxation timescale [s]"
    τ_relax::FT
    "water density [kg/m3]"
    ρw::FT
end

CloudLiquid(::Type{FT}) where {FT <: AbstractFloat} =
    CloudLiquid(CP.create_toml_dict(FT))

function CloudLiquid(toml_dict::CP.AbstractTOMLDict)
    name_map = (;
        :condensation_evaporation_timescale => :τ_relax,
        :density_liquid_water => :ρw,
    )
    parameters =
        CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    FT = CP.float_type(toml_dict)
    return CloudLiquid{FT}(; parameters...)
end

"""
    CloudIce{FT, MS}

The parameters and type for cloud ice condensate

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CloudIce{FT, PD, MS} <: CloudCondensateType{FT}
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
    "ice density [kg/m3]"
    ρi::FT
end

CloudIce(::Type{FT}) where {FT <: AbstractFloat} =
    CloudIce(CP.create_toml_dict(FT))

function CloudIce(toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT))
    name_map = (;
        :density_ice_water => :ρi,
        :cloud_ice_crystals_length_scale => :r0,
        :cloud_ice_size_distribution_coefficient_n0 => :n0,
        :cloud_ice_mass_size_relation_coefficient_me => :me,
        :ice_snow_threshold_radius => :r_ice_snow,
        :sublimation_deposition_timescale => :τ_relax,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    mass = ParticleMass(CloudIce, toml_dict)
    pdf = ParticlePDFIceRain(p.n0)
    FT = CP.float_type(toml_dict)
    P = typeof(pdf)
    M = typeof(mass)
    return CloudIce{FT, P, M}(pdf, mass, p.r0, p.r_ice_snow, p.τ_relax, p.ρi)
end

function ParticleMass(::Type{CloudIce}, td::CP.AbstractTOMLDict)
    name_map = (;
        :density_ice_water => :ρi,
        :cloud_ice_crystals_length_scale => :r0,
        :cloud_ice_mass_size_relation_coefficient_me => :me,
        :cloud_ice_mass_size_relation_coefficient_delm => :Δm,
        :cloud_ice_mass_size_relation_coefficient_chim => :χm,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    m0 = 4 / 3 * π * p.ρi * p.r0^p.me
    FT = CP.float_type(td)
    return ParticleMass{FT}(; p.r0, m0, p.me, p.Δm, p.χm)
end

"""
    Rain{FT, MS, AR, VT, AC}

The parameters and type for rain

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Rain{FT, PD, MS, AR, VT, AC} <: PrecipitationType{FT}
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

Rain(::Type{FT}) where {FT <: AbstractFloat} = Rain(CP.create_toml_dict(FT))

function Rain(toml_dict::CP.AbstractTOMLDict)
    name_map = (;
        :density_liquid_water => :ρ,
        :rain_drop_length_scale => :r0,
        :rain_drop_size_distribution_coefficient_n0 => :n0,
        :rain_mass_size_relation_coefficient_me => :me,
        :rain_cross_section_size_relation_coefficient_ae => :ae,
        :rain_autoconversion_timescale => :τ,
        :cloud_liquid_water_specific_humidity_autoconversion_threshold =>
            :q_threshold,
        :threshold_smooth_transition_steepness => :k,
        :rain_ventillation_coefficient_a => :a,
        :rain_ventillation_coefficient_b => :b,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")

    mass = ParticleMass(Rain, toml_dict)
    area = ParticleArea(Rain, toml_dict)
    pdf = ParticlePDFIceRain(p.n0)
    vent = Ventilation(p.a, p.b)
    acnv1M = Acnv1M(p.τ, p.q_threshold, p.k)

    FT = CP.float_type(toml_dict)
    P = typeof(pdf)
    M = typeof(mass)
    A = typeof(area)
    V = typeof(vent)
    AC = typeof(acnv1M)
    return Rain{FT, P, M, A, V, AC}(pdf, mass, area, vent, acnv1M, p.r0)
end

function ParticleMass(::Type{Rain}, td::CP.AbstractTOMLDict)
    name_map = (;
        :density_liquid_water => :ρ,
        :rain_drop_length_scale => :r0,
        :rain_mass_size_relation_coefficient_me => :me,
        :rain_mass_size_relation_coefficient_delm => :Δm,
        :rain_mass_size_relation_coefficient_chim => :χm,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    m0 = 4 / 3 * π * p.ρ * p.r0^p.me
    FT = CP.float_type(td)
    return ParticleMass{FT}(; p.r0, m0, p.me, p.Δm, p.χm)
end

function ParticleArea(::Type{Rain}, td::CP.AbstractTOMLDict)
    name_map = (;
        :rain_drop_length_scale => :r0,
        :rain_cross_section_size_relation_coefficient_ae => :ae,
        :rain_cross_section_size_relation_coefficient_dela => :Δa,
        :rain_cross_section_size_relation_coefficient_chia => :χa,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    a0 = π * p.r0^p.ae
    FT = CP.float_type(td)
    return ParticleArea{FT}(; a0, p.ae, p.Δa, p.χa)
end

"""
    Snow{FT, PD, MS, AR, VT, AP, AC}

The parameters and type for snow

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Snow{FT, PD, MS, AR, VT, AP, AC} <: PrecipitationType{FT}
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
    "ice density [kg/m3]"
    ρᵢ::FT
end

Snow(::Type{FT}) where {FT <: AbstractFloat} = Snow(CP.create_toml_dict(FT))

function Snow(toml_dict::CP.AbstractTOMLDict)
    name_map = (;
        :cloud_ice_crystals_length_scale => :r0,
        :density_ice_water => :ρi,
        :snow_autoconversion_timescale => :τ,
        :cloud_ice_specific_humidity_autoconversion_threshold =>
            :q_threshold,
        :threshold_smooth_transition_steepness => :k,
        :temperature_water_freeze => :T_freeze,
        :snow_flake_size_distribution_coefficient_mu => :μ,
        :snow_flake_size_distribution_coefficient_nu => :ν,
        :snow_ventillation_coefficient_a => :a,
        :snow_ventillation_coefficient_b => :b,
        :snow_aspect_ratio => :ϕ,
        :snow_aspect_ratio_coefficient => :κ,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")

    mass = ParticleMass(Snow, toml_dict)
    area = ParticleArea(Snow, toml_dict)
    pdf = ParticlePDFSnow(p.μ, p.ν)
    vent = Ventilation(p.a, p.b)
    aspr = SnowAspectRatio(p.ϕ, p.κ)
    acnv1M = Acnv1M(p.τ, p.q_threshold, p.k)
    FT = CP.float_type(toml_dict)
    return Snow{
        FT,
        typeof(pdf),
        typeof(mass),
        typeof(area),
        typeof(vent),
        typeof(aspr),
        typeof(acnv1M),
    }(
        pdf,
        mass,
        area,
        vent,
        aspr,
        acnv1M,
        p.r0,
        p.T_freeze,
        p.ρi,
    )
end

function ParticleMass(::Type{Snow}, td::CP.AbstractTOMLDict)
    name_map = (;
        :density_liquid_water => :ρ,
        :snow_flake_length_scale => :r0,
        :snow_mass_size_relation_coefficient_me => :me,
        :snow_mass_size_relation_coefficient_delm => :Δm,
        :snow_mass_size_relation_coefficient_chim => :χm,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    m0 = 1e-1 * p.r0^p.me
    FT = CP.float_type(td)
    return ParticleMass{FT}(; p.r0, m0, p.me, p.Δm, p.χm)
end

function ParticleArea(::Type{Snow}, toml_dict::CP.AbstractTOMLDict)
    name_map = (;
        :snow_flake_length_scale => :r0,
        :snow_cross_section_size_relation_coefficient => :ae,
        :snow_cross_section_size_relation_coefficient_dela => :Δa,
        :snow_cross_section_size_relation_coefficient_chia => :χa,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    a0 = 0.3 * π * p.r0^p.ae
    FT = CP.float_type(toml_dict)
    return ParticleArea{FT}(; a0, p.ae, p.Δa, p.χa)
end

"""
    CollisionEff{FT}

Collision efficiency parameters for the 1-moment scheme

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct CollisionEff{FT} <: ParametersType{FT}
    "cloud liquid-rain collision efficiency [-]"
    e_liq_rai::FT
    "cloud liquid-snow collision efficiency [-]"
    e_liq_sno::FT
    "cloud ice-rain collision efficiency [-]"
    e_ice_rai::FT
    "cloud ice-snow collision efficiency [-]"
    e_ice_sno::FT
    "rain-snow collision efficiency [-]"
    e_rai_sno::FT
end

CollisionEff(::Type{FT}) where {FT <: AbstractFloat} =
    CollisionEff(CP.create_toml_dict(FT))

function CollisionEff(td::CP.AbstractTOMLDict)
    name_map = (;
        :cloud_liquid_rain_collision_efficiency => :e_liq_rai,
        :cloud_liquid_snow_collision_efficiency => :e_liq_sno,
        :cloud_ice_rain_collision_efficiency => :e_ice_rai,
        :cloud_ice_snow_collision_efficiency => :e_ice_sno,
        :rain_snow_collision_efficiency => :e_rai_sno,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return CollisionEff{FT}(; parameters...)
end
