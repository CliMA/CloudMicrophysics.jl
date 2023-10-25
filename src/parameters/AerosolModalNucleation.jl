export H2S04NucleationParameters,
    OrganicNucleationParameters, MixedNucleationParameters

"""
    H2S04NucleationParameters{FT}

Parameters for pure sulfuric acid nucleation from Dunne et al 1016
DOI:10.1126/science.aaf2649

# Fields
$(DocStringExtensions.FIELDS)
"""
struct H2S04NucleationParameters{FT} <: ParametersType{FT}
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

function H2S04NucleationParameters(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return H2S04NucleationParameters(
        FT(data["mam3_nucleation_p_b_n_neutral"]["value"]),
        FT(data["mam3_nucleation_p_b_i_ion_induced"]["value"]),
        FT(data["mam3_nucleation_u_b_n_neutral"]["value"]),
        FT(data["mam3_nucleation_u_b_i_ion_induced"]["value"]),
        FT(data["mam3_nucleation_v_b_n_neutral"]["value"]),
        FT(data["mam3_nucleation_v_b_i_ion_induced"]["value"]),
        FT(data["mam3_nucleation_w_b_n_neutral"]["value"]),
        FT(data["mam3_nucleation_w_b_i_ion_induced"]["value"]),
        FT(data["mam3_nucleation_p_t_n_neutral"]["value"]),
        FT(data["mam3_nucleation_p_t_i_ion_induced"]["value"]),
        FT(data["mam3_nucleation_u_t_n_neutral"]["value"]),
        FT(data["mam3_nucleation_u_t_i_ion_induced"]["value"]),
        FT(data["mam3_nucleation_v_t_n_neutral"]["value"]),
        FT(data["mam3_nucleation_v_t_i_ion_induced"]["value"]),
        FT(data["mam3_nucleation_w_t_n_neutral"]["value"]),
        FT(data["mam3_nucleation_w_t_i_ion_induced"]["value"]),
        FT(data["mam3_nucleation_p_A_n_neutral"]["value"]),
        FT(data["mam3_nucleation_p_A_i_ion_induced"]["value"]),
        FT(data["mam3_nucleation_a_n_neutral"]["value"]),
        FT(data["mam3_nucleation_a_i_ion_induced"]["value"]),
    )
end

"""
OrganicNucleationParameters{FT}

Parameters for pure organic nucleation from Kirkby 2016
DOI: 10.1038/nature17953

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OrganicNucleationParameters{FT} <: ParametersType{FT}
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

function OrganicNucleationParameters(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return OrganicNucleationParameters(
        FT(data["mam3_nucleation_a_1_neutral"]["value"]),
        FT(data["mam3_nucleation_a_2_neutral"]["value"]),
        FT(data["mam3_nucleation_a_3_ion_induced"]["value"]),
        FT(data["mam3_nucleation_a_4_ion_induced"]["value"]),
        FT(data["mam3_nucleation_a_5"]["value"]),
        FT(data["mam3_nucleation_Y_MTO3_percent"]["value"]),
        FT(data["mam3_nucleation_Y_MTOH_percent"]["value"]),
        FT(data["mam3_nucleation_k_MTO3_organic_factor"]["value"]),
        FT(data["mam3_nucleation_k_MTOH_organic_factor"]["value"]),
        FT(data["mam3_nucleation_exp_MTO3_organic_factor"]["value"]),
        FT(data["mam3_nucleation_exp_MTOH_organic_factor"]["value"]),
    )
end

"""
MixedNucleationParameters{FT}

Parameters for mixed organic and sulfuric acid nucleation from Riccobono et al 2014
DOI:10.1126/science.1243527

# Fields
$(DocStringExtensions.FIELDS)
"""
struct MixedNucleationParameters{FT} <: ParametersType{FT}
    k_H2SO4org::FT
    k_MTOH::FT
    exp_MTOH::FT
end

function MixedNucleationParameters(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return MixedNucleationParameters(
        FT(
            data["mam3_nucleation_k_H2SO4_mixed_organic_sulfuric_acid_factor"]["value"],
        ),
        FT(data["mam3_nucleation_k_MTOH_organic_factor"]["value"]),
        FT(data["mam3_nucleation_exp_MTOH_organic_factor"]["value"]),
    )
end
