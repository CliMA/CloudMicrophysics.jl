export Blk1MVelType, StokesRegimeVelType, SB2006VelType, Chen2022VelType

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
end

function Blk1MVelTypeRain(td::CP.ParamDict)
    name_map = (;
        :snow_flake_length_scale => :r0,
        :rain_terminal_velocity_size_relation_coefficient_ve => :ve,
        :rain_terminal_velocity_size_relation_coefficient_delv => :Δv,
        :rain_terminal_velocity_size_relation_coefficient_chiv => :χv,
        :density_liquid_water => :ρw,
        :rain_drop_drag_coefficient => :C_drag,
        :gravitational_acceleration => :grav,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return Blk1MVelTypeRain{FT}(; parameters...)
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
end

function Blk1MVelTypeSnow(td::CP.ParamDict)
    name_map = (;
        :snow_flake_length_scale => :r0,
        :snow_terminal_velocity_size_relation_coefficient => :ve,
        :snow_terminal_velocity_size_relation_coefficient_delv => :Δv,
        :snow_terminal_velocity_size_relation_coefficient_chiv => :χv,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")

    v0 = 2^(9 / 4) * parameters.r0^parameters.ve
    FT = CP.float_type(td)
    return Blk1MVelTypeSnow{FT}(; parameters..., v0)
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
struct Chen2022VelType{FT, R, SI, LI}
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
