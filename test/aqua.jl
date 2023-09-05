using Test
using CloudMicrophysics
using Aqua

@testset "Quality assurance - Aqua" begin

    @testset "Aqua tests (performance)" begin
        # This tests that we don't accidentally run into
        # https://github.com/JuliaLang/julia/issues/29393
        ua = Aqua.detect_unbound_args_recursively(CloudMicrophysics)
        # Uncomment for debugging:
        # for unbound_arg in ua
        #     @show unbound_arg
        # end

        # There are unbound ambiguities in CM that match this pattern:
        # https://github.com/JuliaTesting/Aqua.jl/issues/86
        @test length(ua) <= 4

        # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
        # Test that we're not introducing method ambiguities across deps
        ambs = Aqua.detect_ambiguities(CloudMicrophysics; recursive = true)
        pkg_match(pkgname, pkdir::Nothing) = false
        pkg_match(pkgname, pkdir::AbstractString) = occursin(pkgname, pkdir)
        filter!(
            x -> pkg_match("CloudMicrophysics", pkgdir(last(x).module)),
            ambs,
        )

        # If the number of ambiguities is less than the limit below,
        # then please lower the limit based on the new number of ambiguities.
        # We're trying to drive this number down to zero to reduce latency.
        # Uncomment for debugging:
        # for method_ambiguity in ambs
        #     @show method_ambiguity
        # end
        @test length(ambs) ≤ 0
    end
    @testset "Aqua tests (additional)" begin
        Aqua.test_undefined_exports(CloudMicrophysics)
        Aqua.test_stale_deps(CloudMicrophysics)
        Aqua.test_deps_compat(CloudMicrophysics)
        Aqua.test_project_extras(CloudMicrophysics)
        Aqua.test_project_toml_formatting(CloudMicrophysics)
        Aqua.test_piracy(CloudMicrophysics)
    end

end
nothing
