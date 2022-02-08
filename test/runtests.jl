#get parameters from file
using TOML
parameter_parse = TOML.parsefile(joinpath(@__DIR__, "parameters.toml"))
src_parameter_dict = Dict{String, Float64}()
for (key, val) in parameter_parse
    # In the future - we will use the full names,
    # src_parameter_dict[key] = val["value"]

    # for now we use the aliases
    src_parameter_dict[val["alias"]] = val["value"]
end

include("aerosol_activation_tests.jl")
include("microphysics_tests.jl")
include("gpu_tests.jl")
