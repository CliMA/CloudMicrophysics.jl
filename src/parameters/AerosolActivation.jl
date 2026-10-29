export AerosolActivationParameters

"""
    AerosolActivationParameters{FT}

Parameters for Abdul-Razzak and Ghan 2000 aerosol activation scheme
DOI: 10.1029/1999JD901161

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AerosolActivationParameters{FT} <: ParametersType
    "molar mass of water [kg/mol]"
    M_w::FT
    "gas constant [J/mol/K]"
    R::FT
    "cloud water density [kg/m3]"
    ρ_w::FT
    "cloud ice density [kg/m3]"
    ρ_i::FT
    "surface tension of water [N/m]"
    σ::FT
    "gravitational acceleration [m/s2]"
    g::FT
    "scaling coefficient in Abdul-Razzak and Ghan 2000 [-]"
    f1::FT
    "scaling coefficient in Abdul-Razzak and Ghan 2000 [-]"
    f2::FT
    "scaling coefficient in Abdul-Razzak and Ghan 2000 [-]"
    g1::FT
    "scaling coefficient in Abdul-Razzak and Ghan 2000 [-]"
    g2::FT
    "power of (zeta / eta) in Abdul-Razzak and Ghan 2000 [-]"
    p1::FT
    "power of (S_m^2 / (zeta + 3 * eta)) in Abdul-Razzak and Ghan 2000 [-]"
    p2::FT
end
