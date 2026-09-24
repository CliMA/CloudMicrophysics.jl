export H2SO4SolutionParameters

"""
    H2SO4SolutionParameters{FT}

Parameters for water activity of H2SO4 solutions
from Luo et al 1995. DOI: 10.1029/94GL02988

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct H2SO4SolutionParameters{FT} <: ParametersType
    "max temperature for which the parameterization is valid [K]"
    T_max::FT
    "min temperature for which the parameterization is valid [K]"
    T_min::FT
    "coefficient [-]"
    w_2::FT
    "coefficient [-]"
    c1::FT
    "coefficient [-]"
    c2::FT
    "coefficient [-]"
    c3::FT
    "coefficient [-]"
    c4::FT
    "coefficient [-]"
    c5::FT
    "coefficient [-]"
    c6::FT
    "coefficient [-]"
    c7::FT
end
