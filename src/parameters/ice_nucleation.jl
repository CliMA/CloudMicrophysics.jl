"""
    IceNucleationParameters{FT}

Parameters for ice nucleation

# Fields
$(DocStringExtensions.FIELDS)
"""

struct IceNucleationParameters{FT} <: AbstractCloudMicrophysicsParameters
    "max allowed supersaturation in Mohler 2006 [-]"
    Sᵢ_max::FT
    "threshold temperature in Mohler 2006 [K]"
    T_thr::FT
    "min Δaw from Koop 2000 [-]"
    Δa_w_min::FT
    "max Δaw from Koop 2000 [-]"
    Δa_w_max::FT
    "coefficient from Koop 2000 [-]"
    c₁::FT
    "coefficient from Koop 2000 [-]"
    c₂::FT
    "coefficient from Koop 2000 [-]"
    c₃::FT
    "coefficient from Koop 2000 [-]"
    c₄::FT
end

function IceNucleationParameters(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT)
    (; data) = toml_dict
    return IceNucleationParameters(
        FT(data["Mohler2006_maximum_allowed_Si"]["value"]),
        FT(data["Mohler2006_threshold_T"]["value"]),
        FT(data["Koop2000_min_delta_aw"]["value"]),
        FT(data["Koop2000_max_delta_aw"]["value"]),
        FT(data["Koop2000_J_hom_coeff1"]["value"]),
        FT(data["Koop2000_J_hom_coeff2"]["value"]),
        FT(data["Koop2000_J_hom_coeff3"]["value"]),
        FT(data["Koop2000_J_hom_coeff4"]["value"]),
    )
end
