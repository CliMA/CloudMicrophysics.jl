export Blk1MVelType, SB2006VelType, Chen2022VelType

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

function Blk1MVelTypeRain(td::CP.AbstractTOMLDict)
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

function Blk1MVelTypeSnow(td::CP.AbstractTOMLDict)
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

function Blk1MVelType(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    rain = Blk1MVelTypeRain(toml_dict)
    snow = Blk1MVelTypeSnow(toml_dict)
    return Blk1MVelType{FT, typeof(rain), typeof(snow)}(rain, snow)
end


"""
    SB2006VelType

The type for precipitation terminal velocity from Seifert and Beheng 2006
(Defined only for rain)

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SB2006VelType{FT} <: TerminalVelocityType{FT}
    ρ0::FT
    aR::FT
    bR::FT
    cR::FT
end

SB2006VelType(::Type{FT}) where {FT <: AbstractFloat} =
    SB2006VelType(CP.create_toml_dict(FT))

function SB2006VelType(td::CP.AbstractTOMLDict)
    name_map = (;
        :SB2006_reference_air_density => :ρ0,
        :SB2006_raindrops_terminal_velocity_coeff_aR => :aR,
        :SB2006_raindrops_terminal_velocity_coeff_bR => :bR,
        :SB2006_raindrops_terminal_velocity_coeff_cR => :cR,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return SB2006VelType{FT}(; parameters...)
end

"""
    Chen2022VelTypeSnowIce

The type for precipitation terminal velocity from Chen et al 2022 for
snow and ice. See Table B3 for parameter definitions.
DOI: 10.1016/j.atmosres.2022.106171

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Chen2022VelTypeSnowIce{FT} <: ParametersType{FT}
    As::FT
    Bs::FT
    Cs::FT
    Es::FT
    Fs::FT
    Gs::FT
    "density of cloud ice [kg/m3]"
    ρᵢ::FT
end

function Chen2022VelTypeSnowIce(toml_dict::CP.AbstractTOMLDict)
    # TODO: These should be array parameters.
    name_map = (;
        :Chen2022_table_B3_As_coeff_1 => :A1,
        :Chen2022_table_B3_As_coeff_2 => :A2,
        :Chen2022_table_B3_As_coeff_3 => :A3,
        :Chen2022_table_B3_Bs_coeff_1 => :B1,
        :Chen2022_table_B3_Bs_coeff_2 => :B2,
        :Chen2022_table_B3_Bs_coeff_3 => :B3,
        :Chen2022_table_B3_Cs_coeff_1 => :C1,
        :Chen2022_table_B3_Cs_coeff_2 => :C2,
        :Chen2022_table_B3_Cs_coeff_3 => :C3,
        :Chen2022_table_B3_Cs_coeff_4 => :C4,
        :Chen2022_table_B3_Es_coeff_1 => :E1,
        :Chen2022_table_B3_Es_coeff_2 => :E2,
        :Chen2022_table_B3_Es_coeff_3 => :E3,
        :Chen2022_table_B3_Fs_coeff_1 => :F1,
        :Chen2022_table_B3_Fs_coeff_2 => :F2,
        :Chen2022_table_B3_Fs_coeff_3 => :F3,
        :Chen2022_table_B3_Gs_coeff_1 => :G1,
        :Chen2022_table_B3_Gs_coeff_2 => :G2,
        :Chen2022_table_B3_Gs_coeff_3 => :G3,
        :density_ice_water => :ρᵢ,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    FT = CP.float_type(toml_dict)
    return Chen2022VelTypeSnowIce{FT}(
        p.A1 * (log(p.ρᵢ))^2 − p.A2 * log(p.ρᵢ) - p.A3,
        FT(1) / (p.B1 + p.B2 * log(p.ρᵢ) + p.B3 / sqrt(p.ρᵢ)),
        p.C1 + p.C2 * exp(p.C3 * p.ρᵢ) + p.C4 * sqrt(p.ρᵢ),
        p.E1 - p.E2 * (log(p.ρᵢ))^2 + p.E3 * sqrt(p.ρᵢ),
        -exp(p.F1 - p.F2 * (log(p.ρᵢ))^2 + p.F3 * log(p.ρᵢ)),
        FT(1) / (p.G1 + p.G2 / (log(p.ρᵢ)) - p.G3 * log(p.ρᵢ) / p.ρᵢ),
        p.ρᵢ,
    )
end

"""
    Chen2022VelTypeRain

The type for precipitation terminal velocity from Chen et al 2022 for
rain. See Table B1 for parameter definitions.
DOI: 10.1016/j.atmosres.2022.106171

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Chen2022VelTypeRain{FT} <: ParametersType{FT}
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

Chen2022VelTypeRain(::Type{FT}) where {FT <: AbstractFloat} =
    Chen2022VelTypeRain(CP.create_toml_dict(FT))

function Chen2022VelTypeRain(td::CP.AbstractTOMLDict)
    name_map = (;
        :Chen2022_table_B1_q_coeff => :ρ0,
        :Chen2022_table_B1_a1_coeff => :a1,
        :Chen2022_table_B1_a2_coeff => :a2,
        :Chen2022_table_B1_a3_coeff => :a3,
        :Chen2022_table_B1_a3_pow_coeff => :a3_pow,
        :Chen2022_table_B1_b1_coeff => :b1,
        :Chen2022_table_B1_b2_coeff => :b2,
        :Chen2022_table_B1_b3_coeff => :b3,
        :Chen2022_table_B1_b_rho_coeff => :b_ρ,
        :Chen2022_table_B1_c1_coeff => :c1,
        :Chen2022_table_B1_c2_coeff => :c2,
        :Chen2022_table_B1_c3_coeff => :c3,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return Chen2022VelTypeRain{FT}(; parameters...)
end

"""
    Chen2022VelType

The type for precipitation terminal velocity from Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
(defied for rain, snow and cloud ice)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Chen2022VelType{FT, R, SI} <: TerminalVelocityType{FT}
    rain::R
    snow_ice::SI
end

function Chen2022VelType(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    rain = Chen2022VelTypeRain(toml_dict)
    snow_ice = Chen2022VelTypeSnowIce(toml_dict)
    return Chen2022VelType{FT, typeof(rain), typeof(snow_ice)}(rain, snow_ice)
end
