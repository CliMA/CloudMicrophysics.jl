export Feldspar

"""
    Feldspar{FT}

Parameters for Feldspar from Alpert et al 2022
DOI: 10.1039/D1EA00077B

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Feldspar{FT} <: AerosolType{FT}
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
end

Feldspar(::Type{FT}) where {FT <: AbstractFloat} =
    Feldspar(CP.create_toml_dict(FT))

function Feldspar(td::CP.AbstractTOMLDict)
    name_map = (;
        :Alpert2022_J_deposition_m_Feldspar => :deposition_m,
        :Alpert2022_J_deposition_c_Feldspar => :deposition_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return Feldspar{FT}(; parameters...)
end
