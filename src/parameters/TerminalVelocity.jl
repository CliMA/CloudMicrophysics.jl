export Blk1MVelType, StokesRegimeVelType, SB2006VelType, Chen2022VelType, TerminalVelocityParams

"""
    Blk1MVelTypeRain

The type for precipitation terminal velocity from the simple 1-moment scheme
for rain

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Blk1MVelTypeRain{FT} <: ParametersType{FT}
    "particle length scale [m]"
    r0::FT
    "rain terminal velocity size relation coefficient [-]"
    ve::FT
    "rain terminal velocity size relation coefficient [-]"
    Δv::FT
    "rain terminal velocity size relation coefficient [-]"
    χv::FT
    "cloud water density [kg/m3]"
    ρw::FT
    "rain drop drag coefficient [-]"
    C_drag::FT
    "gravitational acceleration [m/s2]"
    grav::FT
    "pre-computed gamma((ve + Δv + 5) / 2) for ventilation [-]"
    gamma_vent::FT
    "pre-computed gamma(me + ve + Δm + Δv + 1) for terminal velocity [-]"
    gamma_term::FT
    "pre-computed gamma(ae + ve + Δa + Δv + 1) for accretion [-]"
    gamma_accr::FT
end

function Blk1MVelTypeRain(td::CP.ParamDict)
    vel_map = (;
        :snow_flake_length_scale => :r0,
        :rain_terminal_velocity_size_relation_coefficient_ve => :ve,
        :rain_terminal_velocity_size_relation_coefficient_delv => :Δv,
        :rain_terminal_velocity_size_relation_coefficient_chiv => :χv,
        :density_liquid_water => :ρw,
        :rain_drop_drag_coefficient => :C_drag,
        :gravitational_acceleration => :grav,
    )
    mass_map = (;
        :rain_mass_size_relation_coefficient_me => :me,
        :rain_mass_size_relation_coefficient_delm => :Δm,
    )
    area_map = (;
        :rain_cross_section_size_relation_coefficient_ae => :ae,
        :rain_cross_section_size_relation_coefficient_dela => :Δa,
    )
    parameters = CP.get_parameter_values(td, vel_map, "CloudMicrophysics")
    mass_p = CP.get_parameter_values(td, mass_map, "CloudMicrophysics")
    area_p = CP.get_parameter_values(td, area_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    gamma_vent = FT(SF.gamma((parameters.ve + parameters.Δv + 5) / 2))
    gamma_term = FT(SF.gamma(mass_p.me + parameters.ve + mass_p.Δm + parameters.Δv + 1))
    gamma_accr = FT(SF.gamma(area_p.ae + parameters.ve + area_p.Δa + parameters.Δv + 1))
    return Blk1MVelTypeRain{FT}(; parameters..., gamma_vent, gamma_term, gamma_accr)
end

"""
    Blk1MVelTypeSnow

The type for precipitation terminal velocity from the simple 1-moment scheme
for snow

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Blk1MVelTypeSnow{FT} <: ParametersType{FT}
    "particle length scale [m]"
    r0::FT
    "snow terminal velocity size relation coefficient [-]"
    ve::FT
    "snow terminal velocity size relation coefficient [-]"
    Δv::FT
    "snow terminal velocity size relation coefficient [-]"
    χv::FT
    "snow terminal velocity size relation coefficient [m/s]"
    v0::FT
    "pre-computed gamma((ve + Δv + 5) / 2) for ventilation [-]"
    gamma_vent::FT
    "pre-computed gamma(me + ve + Δm + Δv + 1) for terminal velocity [-]"
    gamma_term::FT
    "pre-computed gamma(ae + ve + Δa + Δv + 1) for accretion [-]"
    gamma_accr::FT
end

function Blk1MVelTypeSnow(td::CP.ParamDict)
    vel_map = (;
        :snow_flake_length_scale => :r0,
        :snow_terminal_velocity_size_relation_coefficient => :ve,
        :snow_terminal_velocity_size_relation_coefficient_delv => :Δv,
        :snow_terminal_velocity_size_relation_coefficient_chiv => :χv,
    )
    mass_map = (;
        :snow_mass_size_relation_coefficient_me => :me,
        :snow_mass_size_relation_coefficient_delm => :Δm,
    )
    area_map = (;
        :snow_cross_section_size_relation_coefficient => :ae,
        :snow_cross_section_size_relation_coefficient_dela => :Δa,
    )
    parameters = CP.get_parameter_values(td, vel_map, "CloudMicrophysics")
    mass_p = CP.get_parameter_values(td, mass_map, "CloudMicrophysics")
    area_p = CP.get_parameter_values(td, area_map, "CloudMicrophysics")

    v0 = 2^(9 / 4) * parameters.r0^parameters.ve
    FT = CP.float_type(td)
    gamma_vent = FT(SF.gamma((parameters.ve + parameters.Δv + 5) / 2))
    gamma_term = FT(SF.gamma(mass_p.me + parameters.ve + mass_p.Δm + parameters.Δv + 1))
    gamma_accr = FT(SF.gamma(area_p.ae + parameters.ve + area_p.Δa + parameters.Δv + 1))
    return Blk1MVelTypeSnow{FT}(; parameters..., v0, gamma_vent, gamma_term, gamma_accr)
end

"""
    Blk1MVelType

The type for precipitation terminal velocity from the simple 1-moment scheme
(defined for rain and snow)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Blk1MVelType{FT, R, S} <: TerminalVelocityType{FT}
    rain::R
    snow::S
end

Blk1MVelType(::Type{FT}) where {FT <: AbstractFloat} =
    Blk1MVelType(CP.create_toml_dict(FT))

function Blk1MVelType(toml_dict::CP.ParamDict)
    rain = Blk1MVelTypeRain(toml_dict)
    snow = Blk1MVelTypeSnow(toml_dict)
    FT = CP.float_type(toml_dict)
    return Blk1MVelType{FT, typeof(rain), typeof(snow)}(rain, snow)
end

"""
    StokesRegimeVelType

The type for precipitation terminal velocity in the Stokes regime (Re < 1)

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct StokesRegimeVelType{FT} <: TerminalVelocityType{FT}
    ρw::FT
    ν_air::FT
    grav::FT
end

StokesRegimeVelType(::Type{FT}) where {FT <: AbstractFloat} =
    StokesRegimeVelType(CP.create_toml_dict(FT))

function StokesRegimeVelType(td::CP.ParamDict)
    name_map = (;
        :density_liquid_water => :ρw,
        :kinematic_viscosity_of_air => :ν_air,
        :gravitational_acceleration => :grav,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return StokesRegimeVelType{FT}(; parameters...)
end

"""
    SB2006VelType

The type for precipitation terminal velocity from Seifert and Beheng 2006

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SB2006VelType{FT} <: TerminalVelocityType{FT}
    ρ0::FT
    aR::FT
    bR::FT
    cR::FT
    ρw::FT
    ν_air::FT
    grav::FT
end

SB2006VelType(::Type{FT}) where {FT <: AbstractFloat} =
    SB2006VelType(CP.create_toml_dict(FT))

function SB2006VelType(td::CP.ParamDict)
    name_map = (;
        :SB2006_reference_air_density => :ρ0,
        :SB2006_raindrops_terminal_velocity_coeff_aR => :aR,
        :SB2006_raindrops_terminal_velocity_coeff_bR => :bR,
        :SB2006_raindrops_terminal_velocity_coeff_cR => :cR,
        :density_liquid_water => :ρw,
        :kinematic_viscosity_of_air => :ν_air,
        :gravitational_acceleration => :grav,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return SB2006VelType{FT}(; parameters...)
end

"""
    Chen2022VelTypeSmallIce

The type for precipitation terminal velocity from Chen et al 2022 for small ice.
See Table B3 for parameter definitions. DOI: 10.1016/j.atmosres.2022.106171

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Chen2022VelTypeSmallIce{FT, N, M} <: TerminalVelocityType{FT}
    A::NTuple{N, FT}
    B::NTuple{N, FT}
    C::NTuple{M, FT}
    E::NTuple{N, FT}
    F::NTuple{N, FT}
    G::NTuple{N, FT}
    "cutoff for small vs large ice particle dimension [m]"
    cutoff::FT
end

function Chen2022VelTypeSmallIce(td::CP.ParamDict)
    # TODO: These should be array parameters.
    name_map = (;
        :Chen2022_table_B3_As => :A,
        :Chen2022_table_B3_Bs => :B,
        :Chen2022_table_B3_Cs => :C,
        :Chen2022_table_B3_Es => :E,
        :Chen2022_table_B3_Fs => :F,
        :Chen2022_table_B3_Gs => :G,
        :Chen2022_ice_cutoff => :cutoff,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    # hack!
    parameters = map(p -> p isa Vector ? Tuple(p) : p, parameters)
    FT = CP.float_type(td)
    return Chen2022VelTypeSmallIce{FT, 3, 4}(; parameters...)
end

"""
    Chen2022VelTypeLargeIce

The type for precipitation terminal velocity from Chen et al 2022 for large ice.
See Table B4 for parameter definitions. DOI: 10.1016/j.atmosres.2022.106171

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Chen2022VelTypeLargeIce{FT, N} <: TerminalVelocityType{FT}
    A::NTuple{N, FT}
    B::NTuple{N, FT}
    C::NTuple{N, FT}
    E::NTuple{N, FT}
    F::NTuple{N, FT}
    G::NTuple{N, FT}
    H::NTuple{N, FT}
    "cutoff for small vs large ice particle dimension [m]"
    cutoff::FT
end

function Chen2022VelTypeLargeIce(td::CP.ParamDict)
    # TODO: These should be array parameters.
    name_map = (;
        :Chen2022_table_B5_Al => :A,
        :Chen2022_table_B5_Bl => :B,
        :Chen2022_table_B5_Cl => :C,
        :Chen2022_table_B5_El => :E,
        :Chen2022_table_B5_Fl => :F,
        :Chen2022_table_B5_Gl => :G,
        :Chen2022_table_B5_Hl => :H,
        :Chen2022_ice_cutoff => :cutoff,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    # hack!
    parameters = map(p -> p isa Vector ? Tuple(p) : p, parameters)
    FT = CP.float_type(td)
    return Chen2022VelTypeLargeIce{FT, 3}(; parameters...)
end

"""
    Chen2022VelTypeRain

The type for precipitation terminal velocity from Chen et al 2022 for
rain. See Table B1 for parameter definitions.
DOI: 10.1016/j.atmosres.2022.106171

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Chen2022VelTypeRain{FT, N} <: TerminalVelocityType{FT}
    ρ0::FT
    a::NTuple{N, FT}
    a3_pow::FT
    b::NTuple{N, FT}
    b_ρ::FT
    c::NTuple{N, FT}
end

Chen2022VelTypeRain(::Type{FT}) where {FT <: AbstractFloat} =
    Chen2022VelTypeRain(CP.create_toml_dict(FT))

function Chen2022VelTypeRain(td::CP.ParamDict)
    name_map = (;
        :Chen2022_table_B1_q_coeff => :ρ0,
        :Chen2022_table_B1_ai => :a,
        :Chen2022_table_B1_a3_pow_coeff => :a3_pow,
        :Chen2022_table_B1_bi => :b,
        :Chen2022_table_B1_b_rho_coeff => :b_ρ,
        :Chen2022_table_B1_ci => :c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    # hack!
    parameters = map(p -> p isa Vector ? Tuple(p) : p, parameters)
    FT = CP.float_type(td)
    return Chen2022VelTypeRain{FT, 3}(; parameters...)
end

"""
    Chen2022VelType

The type for precipitation terminal velocity from Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
(defined for rain, snow and cloud ice)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Chen2022VelType{FT, R, SI, LI} <: TerminalVelocityType{FT}
    rain::R
    small_ice::SI
    large_ice::LI
end

Chen2022VelType(::Type{FT}) where {FT <: AbstractFloat} =
    Chen2022VelType(CP.create_toml_dict(FT))

function Chen2022VelType(toml_dict::CP.ParamDict)
    rain = Chen2022VelTypeRain(toml_dict)
    small_ice = Chen2022VelTypeSmallIce(toml_dict)
    large_ice = Chen2022VelTypeLargeIce(toml_dict)
    FT = CP.float_type(toml_dict)
    return Chen2022VelType{
        FT,
        typeof(rain),
        typeof(small_ice),
        typeof(large_ice),
    }(
        rain,
        small_ice,
        large_ice,
    )
end

"""
    TerminalVelocityParams{FT, STOKES, CHEN, BLK1M}

Unified container for all terminal velocity parameterizations used in CloudMicrophysics.

This struct consolidates terminal velocity parameters for sedimentation calculations
across all microphysics schemes (0M, 1M, 2M, NonEq). It provides a single source of
truth for velocity parameterizations, eliminating duplication and enabling flexible
mixing of velocity types across schemes.

# Fields
- `stokes::STOKES`: [`StokesRegimeVelType`](@ref) — terminal velocity for cloud liquid droplets (Stokes regime)
- `chen2022::CHEN`: [`Chen2022VelType`](@ref) — terminal velocity for rain and ice (Chen et al. 2022)
- `blk1m::BLK1M`: [`Blk1MVelType`](@ref) — terminal velocity for rain and snow (1-moment scheme)

# Usage

Terminal velocities are used for sedimentation (advection of hydrometeors), not for
microphysical process rates (which are handled by scheme-specific parameter structs).

Different schemes use different velocity parameterizations:
- **Cloud liquid**: `stokes` (all schemes)
- **Cloud ice**: `chen2022.small_ice` (all schemes)
- **Rain (1M)**: `blk1m.rain`
- **Rain (2M)**: `chen2022.rain`
- **Snow**: `blk1m.snow` (used by both 1M and 2M)

# Example
```julia
using CloudMicrophysics.Parameters as CMP
using CloudMicrophysics.MicrophysicsNonEq as CMNe
using CloudMicrophysics.Microphysics1M as CM1

# Create unified terminal velocity parameters
tv = CMP.TerminalVelocityParams(Float64)

# Cloud liquid sedimentation (Stokes regime)
v_liq = CMNe.terminal_velocity(cloud_liquid, tv.stokes, ρ, q_liq)

# Cloud ice sedimentation (Chen 2022)
v_ice = CMNe.terminal_velocity(cloud_ice, tv.chen2022.small_ice, ρ, q_ice)

# Rain sedimentation (1M scheme uses Blk1M)
v_rain_1m = CM1.terminal_velocity(rain, tv.blk1m.rain, ρ, q_rai)

# Snow sedimentation (shared by 1M and 2M)
v_snow = CM1.terminal_velocity(snow, tv.blk1m.snow, ρ, q_sno)
```

# See Also
- [`StokesRegimeVelType`](@ref): Stokes regime terminal velocity
- [`Chen2022VelType`](@ref): Chen et al. (2022) terminal velocity
- [`Blk1MVelType`](@ref): 1-moment bulk terminal velocity
"""
struct TerminalVelocityParams{FT, STOKES, CHEN, BLK1M} <: ParametersType{FT}
    stokes::STOKES
    chen2022::CHEN
    blk1m::BLK1M
end

"""
    TerminalVelocityParams(::Type{FT}) where {FT <: AbstractFloat}

Create a `TerminalVelocityParams` object from a floating point type.

# Arguments
- `FT`: Floating point type (e.g., Float64, Float32)
"""
TerminalVelocityParams(::Type{FT}) where {FT <: AbstractFloat} =
    TerminalVelocityParams(CP.create_toml_dict(FT))

"""
    TerminalVelocityParams(toml_dict::CP.ParamDict)

Create a `TerminalVelocityParams` object from a ClimaParams TOML dictionary.

# Arguments
- `toml_dict`: ClimaParams parameter dictionary
"""
function TerminalVelocityParams(toml_dict::CP.ParamDict)
    FT = CP.float_type(toml_dict)

    stokes = StokesRegimeVelType(toml_dict)
    chen2022 = Chen2022VelType(toml_dict)
    blk1m = Blk1MVelType(toml_dict)

    return TerminalVelocityParams{FT, typeof(stokes), typeof(chen2022), typeof(blk1m)}(
        stokes,
        chen2022,
        blk1m,
    )
end
