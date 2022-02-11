#get parameters from file
# read parameters needed for tests
import CLIMAParameters
src_parameter_dict = CLIMAParameters.create_parameter_struct(dict_type = "alias")

#include("parameter_tests.jl")
include("aerosol_activation_tests.jl")
include("microphysics_tests.jl")
include("gpu_tests.jl")
