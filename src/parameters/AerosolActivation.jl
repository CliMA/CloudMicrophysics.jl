export AerosolActivationParameters

"""
    AerosolActivationParameters{FT}

Parameters for Abdul-Razzak and Ghan 2000 aerosol activation scheme
DOI: 10.1029/1999JD901161

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AerosolActivationParameters{FT} <: ParametersType{FT}
    "molar mass of water [kg/mol]"
    M_w::FT
    "gas constant [J/mol/K]"
    R::FT
    "cloud water density [kg/m3]"
    ρ_w::FT
    "surface tension of water [N/m]"
    σ::FT
    "gravitational_acceleration [m/s2]"
    g::FT
end

AerosolActivationParameters(::Type{FT}) where {FT <: AbstractFloat} =
    AerosolActivationParameters(CP.create_toml_dict(FT))

function AerosolActivationParameters(td::CP.AbstractTOMLDict)
    name_map = (;
        :molar_mass_water => :M_w,
        :gas_constant => :R,
        :density_liquid_water => :ρ_w,
        :surface_tension_water => :σ,
        :gravitational_acceleration => :g,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return AerosolActivationParameters{FT}(; parameters...)
end
