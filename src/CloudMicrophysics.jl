module CloudMicrophysics

include("Parameters.jl")
import .Parameters
const CMP = Parameters

# function stubs to be re-defined inside Microphysics1M and Microphysics2M
export conv_q_liq_to_q_rai
export accretion
function conv_q_liq_to_q_rai end
function accretion end

include("CommonTypes.jl")
include("Common.jl")
include("Microphysics0M.jl")
include("Microphysics1M.jl")
include("Microphysics2M.jl")
include("MicrophysicsNonEq.jl")
include("AerosolModel.jl")
include("AerosolActivation.jl")
include("Nucleation.jl")
include("IceNucleation.jl")
include("PrecipitationSusceptibility.jl")
include("PreprocessAerosolData.jl")
include("ActivationEmulatorModels.jl")

end # module
