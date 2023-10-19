export CloudLiquid, CloudIce, Rain, Snow, CollisionEff

"""
    ParticlePDFSnow{FT}

A struct with snow size distribution parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ParticlePDFSnow{FT} <: ParametersType{FT}
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
struct ParticleMass{FT} <: ParametersType{FT}
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
struct ParticleArea{FT} <: ParametersType{FT}
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
struct CloudLiquid{FT} <: CloudCondensateType{FT}
    "condensation evaporation non_equil microphysics relaxation timescale [s]"
    τ_relax::FT
end

function CloudLiquid(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return CloudLiquid(FT(data["condensation_evaporation_timescale"]["value"]))
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
end

function CloudIce(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict

    n0 = FT(data["cloud_ice_size_distribution_coefficient_n0"]["value"])
    pdf = ParticlePDFIceRain(n0)

    ρi = FT(data["density_ice_water"]["value"])
    r0 = FT(data["cloud_ice_crystals_length_scale"]["value"])

    me = FT(data["cloud_ice_mass_size_relation_coefficient_me"]["value"])
    m0 = FT(4 / 3) * π * ρi * r0^me
    Δm = FT(data["cloud_ice_mass_size_relation_coefficient_delm"]["value"])
    χm = FT(data["cloud_ice_mass_size_relation_coefficient_chim"]["value"])
    mass = ParticleMass(r0, m0, me, Δm, χm)

    return CloudIce{FT, typeof(pdf), typeof(mass)}(
        pdf,
        mass,
        r0,
        FT(data["ice_snow_threshold_radius"]["value"]),
        FT(data["sublimation_deposition_timescale"]["value"]),
    )
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

function Rain(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict

    ρ = FT(data["density_liquid_water"]["value"])
    r0 = FT(data["rain_drop_length_scale"]["value"])

    n0 = FT(data["rain_drop_size_distribution_coefficient_n0"]["value"])
    pdf = ParticlePDFIceRain(n0)

    me = FT(data["rain_mass_size_relation_coefficient_me"]["value"])
    m0 = FT(4 / 3) * π * ρ * r0^me
    Δm = FT(data["rain_mass_size_relation_coefficient_delm"]["value"])
    χm = FT(data["rain_mass_size_relation_coefficient_chim"]["value"])
    mass = ParticleMass(r0, m0, me, Δm, χm)

    ae = FT(data["rain_cross_section_size_relation_coefficient_ae"]["value"])
    Δa = FT(data["rain_cross_section_size_relation_coefficient_dela"]["value"]) # typo in .toml file
    χa = FT(data["rain_cross_section_size_relation_coefficient_chia"]["value"])
    a0 = FT(π) * r0^ae
    area = ParticleArea(a0, ae, Δa, χa)

    a = FT(data["rain_ventillation_coefficient_a"]["value"])
    b = FT(data["rain_ventillation_coefficient_b"]["value"])
    vent = Ventilation(a, b)

    τ = FT(data["rain_autoconversion_timescale"]["value"])
    q_threshold = FT(
        data["cloud_liquid_water_specific_humidity_autoconversion_threshold"]["value"],
    )
    k = FT(data["threshold_smooth_transition_steepness"]["value"])
    acnv1M = Acnv1M(τ, q_threshold, k)

    return Rain{
        FT,
        typeof(pdf),
        typeof(mass),
        typeof(area),
        typeof(vent),
        typeof(acnv1M),
    }(
        pdf,
        mass,
        area,
        vent,
        acnv1M,
        r0,
    )
end

"""
    Snow{FT, PD, MS, AR, VT, AC}

The parameters and type for snow

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Snow{FT, PD, MS, AR, VT, AC} <: PrecipitationType{FT}
    "a struct with size distribution parameters"
    pdf::PD
    "a struct with mass size relation parameters"
    mass::MS
    "a struct with cross section size relation parameers"
    area::AR
    "a struct with ventilation coefficients"
    vent::VT
    "a struct with ice to snow autoconversion parameters"
    acnv1M::AC
    "particle length scale [m]"
    r0::FT
    "freezing temperature of water [K]"
    T_freeze::FT
    "ice density [kg/m3]"
    ρᵢ::FT
end

function Snow(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict

    r0 = FT(data["snow_flake_length_scale"]["value"])
    ρi = FT(data["density_ice_water"]["value"])

    μ = FT(data["snow_flake_size_distribution_coefficient_mu"]["value"])
    ν = FT(data["snow_flake_size_distribution_coefficient_nu"]["value"])
    pdf = ParticlePDFSnow(μ, ν)

    me = FT(data["snow_mass_size_relation_coefficient_me"]["value"])
    m0 = FT(1e-1) * r0^me
    Δm = FT(data["snow_mass_size_relation_coefficient_delm"]["value"])
    χm = FT(data["snow_mass_size_relation_coefficient_chim"]["value"])
    mass = ParticleMass(r0, m0, me, Δm, χm)

    ae = FT(data["snow_cross_section_size_relation_coefficient"]["value"])
    Δa = FT(data["snow_cross_section_size_relation_coefficient_dela"]["value"]) # typo in .toml file
    χa = FT(data["snow_cross_section_size_relation_coefficient_chia"]["value"])
    a0 = FT(0.3) * π * r0^ae
    area = ParticleArea(a0, ae, Δa, χa)

    a = FT(data["snow_ventillation_coefficient_a"]["value"])
    b = FT(data["snow_ventillation_coefficient_b"]["value"])
    vent = Ventilation(a, b)

    τ = FT(data["snow_autoconversion_timescale"]["value"])
    q_threshold = FT(
        data["cloud_ice_specific_humidity_autoconversion_threshold"]["value"],
    )
    k = FT(data["threshold_smooth_transition_steepness"]["value"])
    acnv1M = Acnv1M(τ, q_threshold, k)

    return Snow{
        FT,
        typeof(pdf),
        typeof(mass),
        typeof(area),
        typeof(vent),
        typeof(acnv1M),
    }(
        pdf,
        mass,
        area,
        vent,
        acnv1M,
        r0,
        FT(data["temperature_water_freeze"]["value"]),
        ρi,
    )
end

"""
    CollisionEff{FT}

Collision efficiency parameters for the 1-moment scheme

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CollisionEff{FT} <: ParametersType{FT}
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

function CollisionEff(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return CollisionEff(
        FT(data["cloud_liquid_rain_collision_efficiency"]["value"]),
        FT(data["cloud_liquid_snow_collision_efficiency"]["value"]),
        FT(data["cloud_ice_rain_collision_efficiency"]["value"]),
        FT(data["cloud_ice_snow_collision_efficiency"]["value"]),
        FT(data["rain_snow_collision_efficiency"]["value"]),
    )
end
