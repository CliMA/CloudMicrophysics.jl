export H2SO4SolutionParameters

"""
    H2SO4SolutionParameters{FT}

Parameters for water activity of H2SO4 solutions
from Luo et al 1995. DOI: 10.1029/94GL02988

# Fields
$(DocStringExtensions.FIELDS)
"""
struct H2SO4SolutionParameters{FT} <: ParametersType{FT}
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

function H2SO4SolutionParameters(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return H2SO4SolutionParameters(
        FT(data["p_over_sulphuric_acid_solution_T_max"]["value"]),
        FT(data["p_over_sulphuric_acid_solution_T_min"]["value"]),
        FT(data["p_over_sulphuric_acid_solution_w_2"]["value"]),
        FT(data["p_over_sulphuric_acid_solution_c1"]["value"]),
        FT(data["p_over_sulphuric_acid_solution_c2"]["value"]),
        FT(data["p_over_sulphuric_acid_solution_c3"]["value"]),
        FT(data["p_over_sulphuric_acid_solution_c4"]["value"]),
        FT(data["p_over_sulphuric_acid_solution_c5"]["value"]),
        FT(data["p_over_sulphuric_acid_solution_c6"]["value"]),
        FT(data["p_over_sulphuric_acid_solution_c7"]["value"]),
    )
end
