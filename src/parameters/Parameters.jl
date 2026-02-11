"""
    Parameters

A module for CloudMicrophysics.jl free parameters.
"""
module Parameters

using DocStringExtensions

import ClimaParams as CP
import SpecialFunctions as SF

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

# Unified parameter containers
include("Microphysics0MParams.jl")
include("Microphysics1MParams.jl")
include("Microphysics2MParams.jl")


end # module
