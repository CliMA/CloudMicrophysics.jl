"""
    Illite{FT}

Parameters for illite

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Illite{FT} <: AbstractAerosolProperties
    "m coefficient for immersion freezing J [-]"
    m::FT
    "c coefficient for immersion freezing J [-]"
    c::FT
end
Base.broadcastable(x::Illite) = tuple(x)

function Illite(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT)
    (; data) = toml_dict
    return Illite(
        FT(data["KnopfAlpert2013_J_ABIFM_m_Illite"]["value"]),
        FT(data["KnopfAlpert2013_J_ABIFM_c_Illite"]["value"]),
    )
end
