export IceNucleationParameters

"""
    Mohler2006{FT}

Parameters for ice nucleation from Mohler et al 2006
DOI: 10.5194/acp-6-3007-2006

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Mohler2006{FT} <: ParametersType{FT}
    "max allowed supersaturation [-]"
    Sᵢ_max::FT
    "threshold temperature [K]"
    T_thr::FT
end

function Mohler2006(td::CP.AbstractTOMLDict)
    name_map = (;
        :Mohler2006_maximum_allowed_Si => :Sᵢ_max,
        :Mohler2006_threshold_T => :T_thr,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return Mohler2006{FT}(; parameters...)
end

"""
    Koop2000{FT}

Parameters for ice nucleation from Koop et al 2000
DOI: 10.1038/35020537

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Koop2000{FT} <: ParametersType{FT}
    "min Δaw [-]"
    Δa_w_min::FT
    "max Δaw [-]"
    Δa_w_max::FT
    "coefficient [-]"
    c₁::FT
    "coefficient [-]"
    c₂::FT
    "coefficient [-]"
    c₃::FT
    "coefficient [-]"
    c₄::FT
end

function Koop2000(td::CP.AbstractTOMLDict)
    name_map = (;
        :Koop2000_min_delta_aw => :Δa_w_min,
        :Koop2000_max_delta_aw => :Δa_w_max,
        :Koop2000_J_hom_coeff1 => :c₁,
        :Koop2000_J_hom_coeff2 => :c₂,
        :Koop2000_J_hom_coeff3 => :c₃,
        :Koop2000_J_hom_coeff4 => :c₄,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return Koop2000{FT}(; parameters...)
end

"""
    IceNucleationParameters{FT, DEP, HOM}

Parameters for ice nucleation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct IceNucleationParameters{FT, DEP, HOM} <: ParametersType{FT}
    deposition::DEP
    homogeneous::HOM
end

function IceNucleationParameters(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    deposition = Mohler2006(toml_dict)
    homogeneous = Koop2000(toml_dict)
    DEP = typeof(deposition)
    HOM = typeof(homogeneous)
    return IceNucleationParameters{FT, DEP, HOM}(deposition, homogeneous)
end
