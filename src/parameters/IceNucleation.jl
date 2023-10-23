export IceNucleationParameters

"""
    Mohler2006{FT}

Parameters for ice nucleation from Mohler et al 2006
DOI: 10.5194/acp-6-3007-2006

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Mohler2006{FT} <: ParametersType{FT}
    "max allowed supersaturation [-]"
    Sᵢ_max::FT
    "threshold temperature [K]"
    T_thr::FT
end

function Mohler2006(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return Mohler2006(
        FT(data["Mohler2006_maximum_allowed_Si"]["value"]),
        FT(data["Mohler2006_threshold_T"]["value"]),
    )
end

"""
    Koop2000{FT}

Parameters for ice nucleation from Koop et al 2000
DOI: 10.1038/35020537

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Koop2000{FT} <: ParametersType{FT}
    "min Δaw [-]"
    Δa_w_min::FT
    "max Δaw [-]"
    Δa_w_max::FT
    "coefficient [-]"
    c₁::FT
    "coefficient [-]"
    c₂::FT
    "coefficient [-]"
    c₃::FT
    "coefficient [-]"
    c₄::FT
end

function Koop2000(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return Koop2000(
        FT(data["Koop2000_min_delta_aw"]["value"]),
        FT(data["Koop2000_max_delta_aw"]["value"]),
        FT(data["Koop2000_J_hom_coeff1"]["value"]),
        FT(data["Koop2000_J_hom_coeff2"]["value"]),
        FT(data["Koop2000_J_hom_coeff3"]["value"]),
        FT(data["Koop2000_J_hom_coeff4"]["value"]),
    )
end

"""
    IceNucleationParameters{FT, DEP, HOM}

Parameters for ice nucleation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct IceNucleationParameters{FT, DEP, HOM} <: ParametersType{FT}
    deposition::DEP
    homogeneous::HOM
end

function IceNucleationParameters(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    deposition = Mohler2006(FT, toml_dict)
    homogeneous = Koop2000(FT, toml_dict)
    return IceNucleationParameters{FT, typeof(deposition), typeof(homogeneous)}(
        deposition,
        homogeneous,
    )
end
