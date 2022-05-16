module CloudMicrophysics

import CLIMAParameters
const CP = CLIMAParameters
const APS = CP.AbstractParameterSet

include("InternalClimaParams.jl")
import .InternalClimaParams
const ICP = InternalClimaParams

Base.broadcastable(param_set::APS) = Ref(param_set)

include("CommonTypes.jl")
include("Common.jl")
include("Microphysics0M.jl")
include("Microphysics1M.jl")
include("MicrophysicsNonEq.jl")
include("AerosolModel.jl")
include("AerosolActivation.jl")

end # module
