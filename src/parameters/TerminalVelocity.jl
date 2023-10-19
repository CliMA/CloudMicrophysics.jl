export Blk1MVelType, SB2006VelType, Chen2022VelType

"""
    Blk1MVelTypeRain

The type for precipitation terminal velocity from the simple 1-moment scheme
for rain

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Blk1MVelTypeRain{FT} <: ParametersType{FT}
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

function Blk1MVelTypeRain(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return Blk1MVelTypeRain(
        FT(data["snow_flake_length_scale"]["value"]),
        FT(
            data["rain_terminal_velocity_size_relation_coefficient_ve"]["value"],
        ),
        FT(
            data["rain_terminal_velocity_size_relation_coefficient_delv"]["value"],
        ),
        FT(
            data["rain_terminal_velocity_size_relation_coefficient_chiv"]["value"],
        ),
        FT(data["density_liquid_water"]["value"]),
        FT(data["rain_drop_drag_coefficient"]["value"]),
        FT(data["gravitational_acceleration"]["value"]),
    )
end

"""
    Blk1MVelTypeSnow

The type for precipitation terminal velocity from the simple 1-moment scheme
for snow

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Blk1MVelTypeSnow{FT} <: ParametersType{FT}
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

function Blk1MVelTypeSnow(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    r0 = FT(data["snow_flake_length_scale"]["value"])
    ve = FT(data["snow_terminal_velocity_size_relation_coefficient"]["value"])
    return Blk1MVelTypeSnow(
        r0,
        ve,
        FT(
            data["snow_terminal_velocity_size_relation_coefficient_delv"]["value"],
        ),
        FT(
            data["snow_terminal_velocity_size_relation_coefficient_chiv"]["value"],
        ),
        FT(2^(9 / 4)) * r0^ve,
    )
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
    rain = Blk1MVelTypeRain(FT, toml_dict)
    snow = Blk1MVelTypeSnow(FT, toml_dict)
    return Blk1MVelType{FT, typeof(rain), typeof(snow)}(rain, snow)
end


"""
    SB2006VelType

The type for precipitation terminal velocity from Seifert and Beheng 2006
(Defined only for rain)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SB2006VelType{FT} <: TerminalVelocityType{FT}
    ρ0::FT
    aR::FT
    bR::FT
    cR::FT
end

function SB2006VelType(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return SB2006VelType(
        FT(data["SB2006_reference_air_density"]["value"]),
        FT(data["SB2006_raindrops_terminal_velocity_coeff_aR"]["value"]),
        FT(data["SB2006_raindrops_terminal_velocity_coeff_bR"]["value"]),
        FT(data["SB2006_raindrops_terminal_velocity_coeff_cR"]["value"]),
    )
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

function Chen2022VelTypeSnowIce(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
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
    return Chen2022VelTypeSnowIce(
        A[1] * (log(ρᵢ))^2 − A[2] * log(ρᵢ) - A[3],
        FT(1) / (B[1] + B[2] * log(ρᵢ) + B[3] / sqrt(ρᵢ)),
        C[1] + C[2] * exp(C[3] * ρᵢ) + C[4] * sqrt(ρᵢ),
        E[1] - E[2] * (log(ρᵢ))^2 + E[3] * sqrt(ρᵢ),
        -exp(F[1] - F[2] * (log(ρᵢ))^2 + F[3] * log(ρᵢ)),
        FT(1) / (G[1] + G[2] / (log(ρᵢ)) - G[3] * log(ρᵢ) / ρᵢ),
        ρᵢ,
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
struct Chen2022VelTypeRain{FT} <: ParametersType{FT}
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

function Chen2022VelTypeRain(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
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
    return Chen2022VelTypeRain(
        ρ0,
        a1,
        a2,
        a3,
        a3_pow,
        b1,
        b2,
        b3,
        b_ρ,
        c1,
        c2,
        c3,
    )
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
    rain = Chen2022VelTypeRain(FT, toml_dict)
    snow_ice = Chen2022VelTypeSnowIce(FT, toml_dict)
    return Chen2022VelType{FT, typeof(rain), typeof(snow_ice)}(rain, snow_ice)
end
