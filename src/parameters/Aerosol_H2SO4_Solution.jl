export H2SO4SolutionParameters

"""
    H2SO4SolutionParameters{FT}

Parameters for water activity of H2SO4 solutions
from Luo et al 1995. DOI: 10.1029/94GL02988

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct H2SO4SolutionParameters{FT}
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

H2SO4SolutionParameters(::Type{FT}) where {FT <: AbstractFloat} =
    H2SO4SolutionParameters(CP.create_toml_dict(FT))

function H2SO4SolutionParameters(td::CP.AbstractTOMLDict)
    name_map = (;
        :p_over_sulphuric_acid_solution_T_max => :T_max,
        :p_over_sulphuric_acid_solution_T_min => :T_min,
        :p_over_sulphuric_acid_solution_w_2 => :w_2,
        :p_over_sulphuric_acid_solution_c1 => :c1,
        :p_over_sulphuric_acid_solution_c2 => :c2,
        :p_over_sulphuric_acid_solution_c3 => :c3,
        :p_over_sulphuric_acid_solution_c4 => :c4,
        :p_over_sulphuric_acid_solution_c5 => :c5,
        :p_over_sulphuric_acid_solution_c6 => :c6,
        :p_over_sulphuric_acid_solution_c7 => :c7,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return H2SO4SolutionParameters{FT}(; parameters...)
end
