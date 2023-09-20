"""
    H2S04NucleationParameters{FT}

Parameters for pure sulfuric acid nucleation

# Fields
$(DocStringExtensions.FIELDS)
"""

Base.@kwdef struct H2S04NucleationParameters{FT} <:
                   AbstractCloudMicrophysicsParameters
    p_b_n::FT
    p_b_i::FT
    u_b_n::FT
    u_b_i::FT
    v_b_n::FT
    v_b_i::FT
    w_b_n::FT
    w_b_i::FT
    p_t_n::FT
    p_t_i::FT
    u_t_n::FT
    u_t_i::FT
    v_t_n::FT
    v_t_i::FT
    w_t_n::FT
    w_t_i::FT
    p_A_n::FT
    p_A_i::FT
    a_n::FT
    a_i::FT
end

"""
OrganicNucleationParameters{FT}

Parameters for pure organic nucleation

# Fields
$(DocStringExtensions.FIELDS)
"""

Base.@kwdef struct OrganicNucleationParameters{FT} <:
                   AbstractCloudMicrophysicsParameters
    a_1::FT
    a_2::FT
    a_3::FT
    a_4::FT
    a_5::FT
    Y_MTO3::FT
    Y_MTOH::FT
    k_MTO3::FT
    k_MTOH::FT
    exp_MTO3::FT
    exp_MTOH::FT
end

"""
MixedNucleationParameters{FT}

Parameters for mixed organic and sulfuric acid nucleation

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MixedNucleationParameters{FT} <:
                   AbstractCloudMicrophysicsParameters
    k_H2SO4org::FT
    k_MTOH::FT
    exp_MTOH::FT
end

for var in [
    :H2S04NucleationParameters,
    :OrganicNucleationParameters,
    :MixedNucleationParameters,
]
    @eval function $var(::Type{FT}) where {FT}
        toml_dict = CP.create_toml_dict(FT)
        aliases = string.(fieldnames($var))
        pairs = CP.get_parameter_values!(toml_dict, aliases)
        return $var{FT}(; pairs...)
    end
end
