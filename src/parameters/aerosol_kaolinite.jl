"""
    Kaolinite{FT}

Parameters for kaolinite

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Kaolinite{FT} <: AbstractAerosolProperties
    "m coefficient for immersion freezing J [-]"
    m::FT
    "c coefficient for immersion freezing J [-]"
    c::FT
end
Base.broadcastable(x::Kaolinite) = tuple(x)

function Kaolinite(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT)
    (; data) = toml_dict
    return Kaolinite(
        FT(data["KnopfAlpert2013_J_ABIFM_m_Kaolinite"]["value"]),
        FT(data["KnopfAlpert2013_J_ABIFM_c_Kaolinite"]["value"]),
    )
end
