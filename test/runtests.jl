#get parameters from file
# read parameters needed for tests
import CLIMAParameters
src_parameter_dict = CLIMAParameters.create_parameter_struct(dict_type = "alias")

aerosol_parameter_dict = src_parameter_dict
include("aerosol_activation_tests.jl")

micro_parameter_dict = src_parameter_dict
include("microphysics_tests.jl")

gpu_parameter_dict = src_parameter_dict
include("gpu_tests.jl")
