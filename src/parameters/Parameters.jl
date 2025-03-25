"""
    Parameters

A module for CloudMicrophysics.jl free parameters.
"""
module Parameters

using DocStringExtensions

import ClimaParams as CP

# Super-types (dispatch, broadcasting, etc...)
include("AbstractTypes.jl")

# Parameters for moist air and water
include("AirProperties.jl")
include("WaterProperties.jl")

# Parameters for different aerosol species
include("Aerosol_H2SO4_Solution.jl")
include("AerosolATD.jl")
include("AerosolSeasalt.jl")
include("AerosolSulfate.jl")
include("AerosolIllite.jl")
include("AerosolKaolinite.jl")
include("AerosolDesertDust.jl")
include("AerosolFeldspar.jl")
include("AerosolFerrihydrite.jl")
include("AerosolAsianDust.jl")
include("AerosolMiddleEasternDust.jl")
include("AerosolSaharanDust.jl")
include("AerosolDust.jl")

# Parameters for aerosol specific parameterizations
include("AerosolActivation.jl")
include("IceNucleation.jl")
include("AerosolModalNucleation.jl")

# Cloud microphysics parameters
include("Microphysics0M.jl")
include("Microphysics1M.jl")
include("Microphysics2M.jl")
include("MicrophysicsP3.jl")
# Terminal velocity parameters  (can be used with different microph. schemes)
include("TerminalVelocity.jl")

for T in (
    Chen2022VelTypeRain,
    Chen2022VelTypeSmallIce,
    Chen2022VelTypeLargeIce,
    Chen2022VelType,
    CloudLiquid,
)
    @eval Base.Broadcast.broadcastable(x::$T) = x
    @eval Base.ndims(::Type{<:$T}) = 0
    @eval Base.size(::$T) = ()
    @eval Base.@propagate_inbounds Base.getindex(x::$T, i) = x
end

end # module
