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

Blk1MVelType(::Type{FT}) where {FT <: AbstractFloat} =
    Blk1MVelType(CP.create_toml_dict(FT))

function Blk1MVelType(toml_dict::CP.AbstractTOMLDict)
    rain = Blk1MVelTypeRain(toml_dict)
    snow = Blk1MVelTypeSnow(toml_dict)
    FT = CP.float_type(toml_dict)
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
    Al::FT
    Bl::FT
    Cl::FT
    El::FT
    Fl::FT
    Gl::FT
    Hl::FT
    "density of cloud ice [kg/m3]"
    ρᵢ::FT
end

function Chen2022VelTypeSnowIce(toml_dict::CP.AbstractTOMLDict)
    # TODO: These should be array parameters.
    name_map = (;
        :Chen2022_table_B3_As => :As,
        :Chen2022_table_B3_Bs => :Bs,
        :Chen2022_table_B3_Cs => :Cs,
        :Chen2022_table_B3_Es => :Es,
        :Chen2022_table_B3_Fs => :Fs,
        :Chen2022_table_B3_Gs => :Gs,
        :Chen2022_table_B5_Al => :Al,
        :Chen2022_table_B5_Bl => :Bl,
        :Chen2022_table_B5_Cl => :Cl,
        :Chen2022_table_B5_El => :El,
        :Chen2022_table_B5_Fl => :Fl,
        :Chen2022_table_B5_Gl => :Gl,
        :Chen2022_table_B5_Hl => :Hl,
        :density_ice_water => :ρᵢ,
    )
    p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    FT = CP.float_type(toml_dict)
    return Chen2022VelTypeSnowIce{FT}(
        p.As[2] * (log(p.ρᵢ))^2 − p.As[3] * log(p.ρᵢ) + p.As[1],
        FT(1) / (p.Bs[1] + p.Bs[2] * log(p.ρᵢ) + p.Bs[3] / sqrt(p.ρᵢ)),
        p.Cs[1] + p.Cs[2] * exp(p.Cs[3] * p.ρᵢ) + p.Cs[4] * sqrt(p.ρᵢ),
        p.Es[1] - p.Es[2] * (log(p.ρᵢ))^2 + p.Es[3] * sqrt(p.ρᵢ),
        -exp(p.Fs[1] - p.Fs[2] * (log(p.ρᵢ))^2 + p.Fs[3] * log(p.ρᵢ)),
        FT(1) / (p.Gs[1] + p.Gs[2] / (log(p.ρᵢ)) - p.Gs[3] * log(p.ρᵢ) / p.ρᵢ),
        p.Al[1] + p.Al[2] * log(p.ρᵢ) + p.Al[3] * (p.ρᵢ)^(-3 / 2),
        exp(p.Bl[1] + p.Bl[2] * log(p.ρᵢ)^2 + p.Bl[3] * log(p.ρᵢ)),
        exp(p.Cl[1] + p.Cl[2] / log(p.ρᵢ) + p.Cl[3] / p.ρᵢ),
        p.El[1] + p.El[2] * log(p.ρᵢ) * sqrt(p.ρᵢ) + p.El[3] * sqrt(p.ρᵢ),
        p.Fl[1] + p.Fl[2] * log(p.ρᵢ) + p.Fl[3] * exp(-p.ρᵢ),
        (
            p.Gl[1] + p.Gl[2] * log(p.ρᵢ) * sqrt(p.ρᵢ) + p.Gl[3] / sqrt(p.ρᵢ)
        )^(-1),
        p.Hl[1] + p.Hl[2] * (p.ρᵢ)^(5 / 2) - p.Hl[3] * exp(-p.ρᵢ),
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
Base.@kwdef struct Chen2022VelTypeRain{FT, N} <: ParametersType{FT}
    ρ0::FT
    a::NTuple{N, FT}
    a3_pow::FT
    b::NTuple{N, FT}
    b_ρ::FT
    c::NTuple{N, FT}
end

Chen2022VelTypeRain(::Type{FT}) where {FT <: AbstractFloat} =
    Chen2022VelTypeRain(CP.create_toml_dict(FT))

function Chen2022VelTypeRain(td::CP.AbstractTOMLDict)
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
(defied for rain, snow and cloud ice)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Chen2022VelType{FT, R, SI} <: TerminalVelocityType{FT}
    rain::R
    snow_ice::SI
end

Chen2022VelType(::Type{FT}) where {FT <: AbstractFloat} =
    Chen2022VelType(CP.create_toml_dict(FT))

function Chen2022VelType(toml_dict::CP.AbstractTOMLDict)
    rain = Chen2022VelTypeRain(toml_dict)
    snow_ice = Chen2022VelTypeSnowIce(toml_dict)
    FT = CP.float_type(toml_dict)
    return Chen2022VelType{FT, typeof(rain), typeof(snow_ice)}(rain, snow_ice)
end
