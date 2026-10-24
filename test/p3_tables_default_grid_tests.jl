# The tables on the default grid, against the targets of the accuracy section of the P3 documentation.
# The suite in runtests.jl does not include this file.
include("p3_tables_helpers.jl")

@testset "P3 tables on the default grid ($FT)" for FT in (Float64, Float32)
    ice = CMP.P3IceParams(CP.create_toml_dict(FT))
    tables = P3.with_tables(ice).integration
    test_ice_table_accuracy(ice, tables; n_states = 2000, median = 0.01, p99 = 0.03)
    q = collision_error_quantiles(ice, tables, 2000)
    for name in (:cloud_mass, :rain_mass, :cloud_number, :rain_number, :frozen_mass, :rime_volume)
        @test q[name].median < 0.01
        @test q[name].p99 < 0.03
        @test q[name].mismatch == 0
    end
    for name in (:cloud_frozen, :cloud_shed, :rain_frozen, :rain_shed)
        @test q[name].median < 0.01
        @test q[name].p99 < 0.03
    end
    @test q.shed_fraction.p99 < 0.01
    for name in (:ice_mass_tendency, :rime_mass_tendency, :rime_volume_tendency)
        @test q[name].median < 0.01
        @test q[name].p99 < 0.03
    end
end
nothing
