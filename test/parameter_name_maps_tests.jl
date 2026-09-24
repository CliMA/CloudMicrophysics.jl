import Test as TT

import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP

"""
Every parameter struct and 1-moment process option registered in `Parameters.parameter_groups()`
is constructed on a fresh ClimaParams dictionary, and the parameters it reads (the `used_in` tags
ClimaParams records) are compared with its `name_map`:
  - every name in the map is actually read by the constructor (no stale map entries),
  - every name read by the constructor is listed in some registered map (no undocumented
    parameters; nested constructors such as `ParticleMass(Rain, td)` account for the rest).
This is what keeps the generated parameter documentation in sync with the code.
"""
function test_name_maps(FT)
    TT.@testset "ClimaParams name maps ($FT)" begin
        summary = CMP.name_map_summary()
        all_names = union((Set(keys(nm)) for nm in values(summary))...)
        TT.@test length(summary) == sum(length(ks) for (_, ks) in CMP.parameter_groups())
        for (title, group_keys) in CMP.parameter_groups(), key in group_keys
            td = CP.create_toml_dict(FT)
            obj = build_from_key(key, td)
            TT.@test obj !== nothing
            used = Set(Symbol(n) for (n, v) in td.data if haskey(v, "used_in"))
            nm = CMP.name_map(key)
            # local names are unique within a map
            TT.@test allunique(values(nm))
            # mapped-but-unread names would make the documentation lie
            TT.@test isempty(setdiff(Set(keys(nm)), used))
            # read-but-unmapped names would be missing from the documentation
            TT.@test isempty(setdiff(used, all_names))
            if key isa Type && key <: CMP.ParametersType && key in CMP.PLAIN_PARAMETER_TYPES
                # generated constructors: every local name is a struct field (fields with a
                # `@kwdef` default, e.g. `Frostenberg2023.log_a`, need not be in the map)
                TT.@test Set(values(nm)) ⊆ Set(fieldnames(key))
            end
        end
    end
end

build_from_key(T::Type, td) = T <: CMP.MicrophysicsOption ? CMP.process_params_for(T(), td) : T(td)
build_from_key(key::Tuple, td) = key[1](key[2], td)

TT.@testset "Parameter name maps" begin
    test_name_maps(Float64)
    test_name_maps(Float32)
end
