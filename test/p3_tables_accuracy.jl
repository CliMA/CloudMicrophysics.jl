# The errors of the P3 lookup tables against the Gauss-Legendre quadrature of order 16 of the
# parameters. The test suite does not include this file; run it with
# `julia --project=test test/p3_tables_accuracy.jl`. To measure a proposed grid, build tables on it
# and print their errors:
#
#     ice = CMP.P3IceParams(CP.create_toml_dict(Float64))
#     grid = P3.P3TableGrid(ice; slope = (; logλ = 256, w = 16, ρ_rim = 16))
#     print_table_accuracy(table_accuracy(ice, P3.with_tables(ice; grid).integration))
include("p3_tables_helpers.jl")

# The quantiles of the errors of `tables` against the quadrature of `ice` on `n_states` seeded
# states, relative errors except the absolute error of `logλ`
table_accuracy(ice, tables; n_states = 2000) = ice_error_quantiles(ice, tables, n_states)

function print_table_accuracy(q)
    println(rpad("Quantity", 24), lpad("Median", 12), lpad("99th pct.", 12), lpad("Maximum", 12))
    for (name, v) in pairs(q)
        values = (lpad(string(round(v[k]; sigdigits = 2)), 12) for k in (:median, :p99, :max))
        println(rpad(string(name), 24), values...)
    end
end

@testset "P3 tables on the default grid ($FT)" for FT in (Float64, Float32)
    ice = CMP.P3IceParams(CP.create_toml_dict(FT))
    q = table_accuracy(ice, P3.with_tables(ice).integration)
    println("Errors of the P3 tables on the default grid in $FT")
    print_table_accuracy(q)
    @test q.logλ.p99 < 1e-3
    for name in filter(!=(:logλ), keys(q))
        @test q[name].median < 0.01
        @test q[name].p99 < 0.03
        @test q[name].mismatch == 0
    end
end
nothing
