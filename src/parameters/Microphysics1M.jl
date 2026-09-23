export CloudLiquid, CloudIce, Rain, Snow, VarTimescaleAcnv, KesslerAcnv

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

A struct with cloud ice or rain size distribution parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ParticlePDFIceRain{FT} <: ParametersType
    "size distribution coefficient [1/m⁴]"
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

A struct with threshold autoconversion parameters (used for the snow autoconversion of
cloud ice, `NoSupersaturation`). Rain autoconversion (`Kessler1M`) uses [`KesslerAcnv`](@ref).

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
    p = make_params(toml_dict, name_map(CloudIce))
    mass = ParticleMass(CloudIce, toml_dict)
    pdf = ParticlePDFIceRain(p.n0)
    return CloudIce(; pdf, mass, p.ρᵢ, p.r_eff, p.N_0)
end

"""
    IceNumberTemperatureFit{FT}

Parameters of the prescribed exponential temperature dependence of the cloud ice
number concentration used by the `TemperatureDependentIceNumber` cloud ice
formation option (see `MicrophysicsNonEq.ice_number_concentration`):

    N_ice(T) = min(N_ref exp(a + b max(T_freeze - T, 0)), N_max)

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct IceNumberTemperatureFit{FT} <: ParametersType
    "prefactor of the fit [1/m³]"
    N_ref::FT
    "intercept of the exponent [-]"
    a::FT
    "slope of the exponent [1/K]"
    b::FT
    "upper bound of the ice number concentration [1/m³]"
    N_max::FT
    "freezing temperature [K]"
    T_freeze::FT
end

function ParticleMass(::Type{CloudIce}, td::CP.ParamDict)
    p = make_params(td, name_map(ParticleMass, CloudIce))
    m0 = p.ρᵢ * p.r0^p.me * π * 4 / 3
    gamma_coeff = SF.gamma(p.me + p.Δm + 1)
    return ParticleMass(; p.r0, m0, p.me, p.Δm, p.χm, gamma_coeff)
end

"""
    Rain{PD, MS, AR, VT}

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
    p = make_params(toml_dict, name_map(Rain))
    return Rain(;
        mass = ParticleMass(Rain, toml_dict),
        area = ParticleArea(Rain, toml_dict),
        pdf = ParticlePDFIceRain(p.n0),
        vent = Ventilation(p.a, p.b),
    )
end

function ParticleMass(::Type{Rain}, td::CP.ParamDict)
    p = make_params(td, name_map(ParticleMass, Rain))
    m0 = p.ρ * p.r0^p.me * π * 4 / 3
    gamma_coeff = SF.gamma(p.me + p.Δm + 1)
    return ParticleMass(; p.r0, m0, p.me, p.Δm, p.χm, gamma_coeff)
end

function ParticleArea(::Type{Rain}, td::CP.ParamDict)
    p = make_params(td, name_map(ParticleArea, Rain))
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
    "snow apparent density [kg/m³], used in the Chen et al. (2022) coefficients"
    ρᵢ::FT
    "bulk ice density [kg/m³], used in the spheroid aspect-ratio relation"
    ρᵢ_bulk::FT
    "pre-computed gamma(α+4)/6 for oblate aspect ratio [-]"
    gamma_aspect_oblate::FT
    "pre-computed gamma(α+4)/6 for prolate aspect ratio [-]"
    gamma_aspect_prolate::FT
end

function Snow(toml_dict::CP.ParamDict)
    p = make_params(toml_dict, name_map(Snow))

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
        p.ρᵢ_bulk,
        gamma_aspect_oblate = SF.gamma(α_oblate + 4) / SF.gamma(FT(4)),
        gamma_aspect_prolate = SF.gamma(α_prolate + 4) / SF.gamma(FT(4)),
    )
end

function ParticleMass(::Type{Snow}, td::CP.ParamDict)
    p = make_params(td, name_map(ParticleMass, Snow))
    m0 = p.r0^p.me / 10
    gamma_coeff = SF.gamma(p.me + p.Δm + 1)
    return ParticleMass(; p.r0, m0, p.me, p.Δm, p.χm, gamma_coeff)
end

function ParticleArea(::Type{Snow}, toml_dict::CP.ParamDict)
    p = make_params(toml_dict, name_map(ParticleArea, Snow))
    FT = CP.float_type(toml_dict)
    a0 = FT(0.3 * π * p.r0^p.ae)
    return ParticleArea(; a0, p.ae, p.Δa, p.χa)
end

"""
    VarTimescaleAcnv{FT}

Parameters for the variable-timescale rain autoconversion scheme used in the 1-moment
microphysics scheme, following Azimi et al. (2023).
Active when `PrescribedNd` is selected in `Microphysics1MOptions`.

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

"""
    KesslerAcnv{FT}

Parameters of the 1-moment Kessler rain autoconversion.
Active when `Kessler1M` is selected in in `Microphysics1MOptions`.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct KesslerAcnv{FT} <: ParametersType
    "Autoconversion timescale for quiescent (stratiform) conditions [s]"
    τ_slow::FT
    "Autoconversion timescale for strong vertical motions (convective) [s]"
    τ_fast::FT
    "Autoconversion threshold for quiescent (stratiform) conditions [kg/kg]"
    q_threshold_slow::FT
    "Autoconversion threshold for strong vertical motions (convective) [kg/kg]"
    q_threshold_fast::FT
    "Velocity scale of the blending [m/s]"
    w_0::FT
    "Threshold smooth transition steepness [-]"
    k::FT
end

function KesslerAcnv(td::CP.ParamDict)
    parameters = make_params(td, name_map(KesslerAcnv))
    parameters.w_0 > 0 ||
        throw(ArgumentError("rain_autoconversion_velocity_scale (w_0) must be positive, got $(parameters.w_0)"))
    (parameters.τ_slow > 0 && parameters.τ_fast > 0) ||
        throw(
            ArgumentError(
                "rain autoconversion timescales must be positive, got τ_slow = $(parameters.τ_slow), τ_fast = $(parameters.τ_fast)",
            ),
        )
    (parameters.q_threshold_slow >= 0 && parameters.q_threshold_fast >= 0) ||
        throw(ArgumentError("rain autoconversion thresholds must be non-negative"))
    return KesslerAcnv(; parameters...)
end
