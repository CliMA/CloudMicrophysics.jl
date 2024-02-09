export H2S04NucleationParameters,
    OrganicNucleationParameters, MixedNucleationParameters

"""
    H2S04NucleationParameters{FT}

Parameters for pure sulfuric acid nucleation from Dunne et al 1016
DOI:10.1126/science.aaf2649

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct H2S04NucleationParameters{FT} <: ParametersType{FT}
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

H2S04NucleationParameters(::Type{FT}) where {FT <: AbstractFloat} =
    H2S04NucleationParameters(CP.create_toml_dict(FT))

function H2S04NucleationParameters(td::CP.AbstractTOMLDict)
    name_map = (;
        :mam3_nucleation_p_b_n_neutral => :p_b_n,
        :mam3_nucleation_p_b_i_ion_induced => :p_b_i,
        :mam3_nucleation_u_b_n_neutral => :u_b_n,
        :mam3_nucleation_u_b_i_ion_induced => :u_b_i,
        :mam3_nucleation_v_b_n_neutral => :v_b_n,
        :mam3_nucleation_v_b_i_ion_induced => :v_b_i,
        :mam3_nucleation_w_b_n_neutral => :w_b_n,
        :mam3_nucleation_w_b_i_ion_induced => :w_b_i,
        :mam3_nucleation_p_t_n_neutral => :p_t_n,
        :mam3_nucleation_p_t_i_ion_induced => :p_t_i,
        :mam3_nucleation_u_t_n_neutral => :u_t_n,
        :mam3_nucleation_u_t_i_ion_induced => :u_t_i,
        :mam3_nucleation_v_t_n_neutral => :v_t_n,
        :mam3_nucleation_v_t_i_ion_induced => :v_t_i,
        :mam3_nucleation_w_t_n_neutral => :w_t_n,
        :mam3_nucleation_w_t_i_ion_induced => :w_t_i,
        :mam3_nucleation_p_A_n_neutral => :p_A_n,
        :mam3_nucleation_p_A_i_ion_induced => :p_A_i,
        :mam3_nucleation_a_n_neutral => :a_n,
        :mam3_nucleation_a_i_ion_induced => :a_i,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return H2S04NucleationParameters{FT}(; parameters...)
end


"""
OrganicNucleationParameters{FT}

Parameters for pure organic nucleation from Kirkby 2016
DOI: 10.1038/nature17953

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct OrganicNucleationParameters{FT} <: ParametersType{FT}
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

OrganicNucleationParameters(::Type{FT}) where {FT <: AbstractFloat} =
    OrganicNucleationParameters(CP.create_toml_dict(FT))

function OrganicNucleationParameters(td::CP.AbstractTOMLDict)
    name_map = (;
        :mam3_nucleation_a_1_neutral => :a_1,
        :mam3_nucleation_a_2_neutral => :a_2,
        :mam3_nucleation_a_3_ion_induced => :a_3,
        :mam3_nucleation_a_4_ion_induced => :a_4,
        :mam3_nucleation_a_5 => :a_5,
        :mam3_nucleation_Y_MTO3_percent => :Y_MTO3,
        :mam3_nucleation_Y_MTOH_percent => :Y_MTOH,
        :mam3_nucleation_k_MTO3_organic_factor => :k_MTO3,
        :mam3_nucleation_k_MTOH_organic_factor => :k_MTOH,
        :mam3_nucleation_exp_MTO3_organic_factor => :exp_MTO3,
        :mam3_nucleation_exp_MTOH_organic_factor => :exp_MTOH,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return OrganicNucleationParameters{FT}(; parameters...)
end

"""
MixedNucleationParameters{FT}

Parameters for mixed organic and sulfuric acid nucleation from Riccobono et al 2014
DOI:10.1126/science.1243527

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MixedNucleationParameters{FT} <: ParametersType{FT}
    k_H2SO4org::FT
    k_MTOH::FT
    exp_MTOH::FT
end

MixedNucleationParameters(::Type{FT}) where {FT <: AbstractFloat} =
    MixedNucleationParameters(CP.create_toml_dict(FT))

function MixedNucleationParameters(td::CP.AbstractTOMLDict)
    name_map = (;
        :mam3_nucleation_k_H2SO4_mixed_organic_sulfuric_acid_factor =>
            :k_H2SO4org,
        :mam3_nucleation_k_MTOH_organic_factor => :k_MTOH,
        :mam3_nucleation_exp_MTOH_organic_factor => :exp_MTOH,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return MixedNucleationParameters{FT}(; parameters...)
end
