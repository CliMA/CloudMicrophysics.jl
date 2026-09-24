export Seasalt

"""
    Seasalt{FT}

Parameters for seasalt

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Seasalt{FT} <: AerosolType
    "molar mass [kg/mol]"
    M::FT
    "density [kg/m3]"
    ρ::FT
    "osmotic coefficient [-]"
    ϕ::FT
    "ion number [-]"
    ν::FT
    "water soluble mass fraction [-]"
    ϵ::FT
    "hygroscopicity parameter [-]"
    κ::FT
end
