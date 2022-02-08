using Test


### read parameters needed for tests
using TOML
parameter_parse = TOML.parsefile(joinpath(@__DIR__, "parameters.toml"))
param_dict = Dict{String, Any}()
for (key, val) in parameter_parse
    # In the future - we will use the full names,
    # param_dict[key] = val["value"]

    # for now we use the aliases
    param_dict[val["alias"]] = val["value"]
end
full_parameter_set = param_dict
universal_constant_aliases = [
    "gas_constant",
    "light_speed",
    "h_Planck",
    "k_Boltzmann",
    "Stefan",
    "astro_unit",
    "avogad",
]



#get CLIMAParameters to check consistency
import CLIMAParameters
const CP = CLIMAParameters
const CP_planet = CLIMAParameters.Planet
const CP_micro = CLIMAParameters.Atmos.Microphysics
const CP_micro_0M = CLIMAParameters.Atmos.Microphysics_0M

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const param_set_cpp = EarthParameterSet()

using Thermodynamics


@testset "parameter read tests" begin

    #check parameter_file agrees with current CP defaults
    for (k, v) in full_parameter_set
        if ~(k in universal_constant_aliases)
            try 
                cp_k = getfield(CP_planet, Symbol(k))
            catch
                try
                    cp_k = getfield(CP_micro, Symbol(k))
                catch
                    cp_k = getfield(CP_micro_0M, Symbol(k))
                end
            end
            @test (v ≈ cp_k(param_set_cpp))
        else
            cp_k = getfield(CP, Symbol(k))
            @test (v ≈ cp_k())
        end
    end

    #check new local implementation agrees with CP
    param_set_by_alias =
        Thermodynamics.ThermodynamicsParameters(full_parameter_set)
    for fn in fieldnames(typeof(param_set_by_alias))
        v = getfield(param_set_by_alias, fn)
        if ~(String(fn) in universal_constant_aliases)
            try
                cp_fn = getfield(CP_planet, fn)
            catch
                try
                    cp_fn = getfield(CP_micro, fn)
                catch
                    cp_fn = getfield(CP_micro_0M, fn)
                end
            end
            @test (v ≈ cp_fn(param_set_cpp))
        else
            cp_fn = getfield(CP, fn)
            @test (v ≈ cp_fn())
        end
    end
end
