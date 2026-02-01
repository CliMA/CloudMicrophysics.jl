export Microphysics0MParams

"""
    Microphysics0MParams{FT, P}

Unified parameter container for 0-moment microphysics.

# Fields
- `precip::P`: Parameters0M — precipitation removal parameters
"""
struct Microphysics0MParams{FT, P} <: ParametersType{FT}
    precip::P
end

"""
    Microphysics0MParams(::Type{FT}) where {FT <: AbstractFloat}

Create a `Microphysics0MParams` object from a floating point type.
"""
Microphysics0MParams(::Type{FT}) where {FT <: AbstractFloat} =
    Microphysics0MParams(CP.create_toml_dict(FT))

"""
    Microphysics0MParams(toml_dict::CP.ParamDict)

Create a `Microphysics0MParams` object from a ClimaParams TOML dictionary.
"""
function Microphysics0MParams(toml_dict::CP.ParamDict)
    precip = Parameters0M(toml_dict)
    FT = CP.float_type(toml_dict)
    return Microphysics0MParams{FT, typeof(precip)}(precip)
end
