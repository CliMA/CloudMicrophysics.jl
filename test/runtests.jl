#get parameters from file
using TOML
parameter_parse = TOML.parsefile(joinpath(@__DIR__, "parameters.toml"))
src_parameter_dict = Dict{Symbol, Float64}()
for (key, val) in parameter_parse
    # In the future - we will use the full names,
    # src_parameter_dict[Symbol(key)] = val["value"]

    # for now we use the aliases
    src_parameter_dict[Symbol(val["alias"])] = val["value"]
end

#parameters for testing (not in source):
test_parameter_dict = parameter_parse["molmass_ratio"]


include("aerosol_activation_tests.jl")
include("microphysics_tests.jl")
include("gpu_tests.jl")
