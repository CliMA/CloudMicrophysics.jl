include("p3_tables_helpers.jl")

function test_uniform_table(FT)
    @testset "UniformTable lookup" begin
        axes = (P3.TableAxis(-1, 2, 4), P3.TableAxis(0, 1, 3), P3.TableAxis(10, 20, 5))
        f(x, y, z) = 1 + x - 2y + z / 10 + x * y * z / 100
        g(x, y, z) = exp(x - y + z / 20)
        point((x, y, z)) = P3.SA.SVector(f(x, y, z), g(x, y, z))
        t = P3.UniformTable(P3.tabulate(point, axes), axes, FT; logs = (false, true))
        @test t isa P3.UniformTable{FT, 3, 2} && size(t.data) == (2, 4, 3, 5)
        # Multilinear interpolation reproduces a multilinear function, and of the logarithm an
        # exponential of one
        rtol = 10eps(FT)
        for (x, y, z) in ((-1, 0, 10), (2, 1, 20), (0.3, 0.7, 12.5), (1.9, 0.01, 19.9))
            (v_f, v_g) = P3.lookup(t, FT(x), FT(y), FT(z))
            @test v_f ≈ f(x, y, z) rtol = rtol
            @test exp(v_g) ≈ g(x, y, z) rtol = 10rtol
        end
        # Coordinates outside the grid move to the nearest edge, NaN to the first point
        @test P3.lookup(t, FT(-5), FT(2), FT(Inf))[1] ≈ f(-1, 1, 20) rtol = rtol
        @test P3.lookup(t, FT(-Inf), FT(NaN), FT(NaN))[1] ≈ f(-1, 0, 10) rtol = rtol
        @test P3.uniform_index(FT(NaN), FT(0), FT(1), 5) == (1, 0)
        @test P3.uniform_index(FT(Inf), FT(0), FT(1), 5) == (4, 1)
        # A coordinate with partials gives the derivative of the interpolation along its axis, and
        # none outside the grid
        @test FD.derivative(x -> P3.lookup(t, x, FT(0.5), FT(15))[1], FT(0.3)) ≈ 1 + FT(0.5 * 15 / 100) rtol = rtol
        @test FD.derivative(x -> P3.lookup(t, x, FT(0.5), FT(15))[1], FT(3)) == 0
        @test (@inferred P3.lookup(t, FT(0.3), FT(0.7), FT(12.5))) isa P3.SA.SVector{2, FT}
        @test (BT.@ballocated P3.lookup($t, $(FT(0.3)), $(FT(0.7)), $(FT(12.5))) samples = 1 evals = 1) == 0

        @test_throws ArgumentError P3.TableAxis(0, 1, 1)
        @test_throws ArgumentError P3.TableAxis(1, 0, 4)
        @test P3.axis_points(P3.TableAxis(0, 1, 3)) == [0, 0.5, 1]
        @test P3.grid_point(axes, CartesianIndex(2, 3, 5)) == (0.0, 1.0, 20.0)
        values = P3.tabulate(x -> f(x...), axes)
        @test size(values) == (4, 3, 5) && values[2, 3, 5] == f(0, 1, 20)
        @test values == map(I -> f(P3.grid_point(axes, I)...), CartesianIndices((4, 3, 5)))
        @test_throws DimensionMismatch P3.UniformTable(zeros(4, 3, 4), axes, FT)
        # A log component must be positive, and every value finite
        @test_throws DomainError P3.UniformTable(zeros(4, 3, 5), axes, FT; logs = (true,))
        @test_throws DomainError P3.UniformTable(fill(NaN, 4, 3, 5), axes, FT)
        @test all(n -> P3.bisection_steps(n) == ceil(Int, log2(n - 1)), 2:300)
    end

    @testset "Table of a point function on an inner axis" begin
        inner = P3.TableAxis(0, 1, 5)
        axes = (P3.TableAxis(-1, 1, 3), P3.TableAxis(0, 2, 4))
        point((x, y)) = P3.SA.SVector(ntuple(k -> x * y + P3.axis_point(inner, k), 5))
        t = P3.inner_axis_table(P3.tabulate(point, axes), inner, axes, FT)
        @test t isa P3.UniformTable{FT, 3, 1} && size(t.data) == (1, 5, 3, 4)
        @test P3.lookup(t, FT(0.3), FT(0.5), FT(1.5))[1] ≈ FT(0.3 + 0.75) rtol = 10eps(FT)
        @test_throws DimensionMismatch P3.inner_axis_table(P3.tabulate(point, axes), P3.TableAxis(0, 1, 4), axes, FT)
    end
end

function test_table_grid(FT)
    @testset "P3TableGrid bounds from the parameters" begin
        ice = CMP.P3IceParams(CP.create_toml_dict(FT))
        params = ice.scheme
        grid = P3.P3TableGrid(ice)
        @test (grid.slope.logλ.lo, grid.slope.logλ.hi) == Float64.(P3.logλ_bounds(params))
        @test (grid.slope.w.lo, grid.slope.w.hi) == (0, 1)
        @test (grid.ice.ρ_rim.lo, grid.ice.ρ_rim.hi) == Float64.(P3.rime_density_bounds(params))
        # The integral tables start at the slope of the largest unrimed mean mass, the smallest slope
        @test (grid.ice.logλ.lo, grid.ice.logλ.hi) == P3.integral_logλ_bounds(params)
        unrimed = P3.P3Particles(params, FT(0), params.ρ_i)
        @test exp(P3.logLdivN(unrimed, grid.ice.logλ.lo)) ≈ P3.ice_mean_particle_mass_max(params) rtol = 1e-3
        rimed = P3.P3Particles(params, FT(0.99), FT(500))
        @test P3.get_distribution_logλ(rimed, P3.ice_mean_particle_mass_max(params), FT(1)) > grid.ice.logλ.lo
        @test grid.liquid_ice.logλ.lo == grid.ice.logλ.lo == grid.wet_set.logλ.lo
        @test all(exp.((grid.ice.logρₐ.lo, grid.ice.logρₐ.hi)) .≈ (0.005, 1.7))
        @test exp(grid.liquid_ice.cloud_logx̄.lo) ≈ ice.cloud_pdf.xc_min
        # The limited rain distribution has its slope within its limits
        (; ρw, λ_min, λ_max) = ice.rain_pdf
        @test exp(grid.liquid_ice.rain_logx̄.hi) ≈ π * ρw / λ_min^3
        @test exp(grid.liquid_ice.rain_logx̄.lo) ≈ π * ρw / λ_max^3
        @test grid.rime_density.logR.lo == -6 && grid.rime_density.logR.hi ≈ log(12) + 6
        @test (grid.particle_rime_density.logR.lo, grid.particle_rime_density.logR.hi) ==
              (grid.rime_density.logR.lo, grid.rime_density.logR.hi)
        @test grid.particle_rime_density.logD.lo == grid.wet_set.logD.lo
        # The number logit spans the bounds of the collision integrals
        p = P3.liquid_ice_collision_quantile(Float64)
        @test grid.wet_set.number_logit.hi ≈ -grid.wet_set.number_logit.lo ≈ log((1 - p) / p)
        for logλ in (grid.wet_set.logλ.lo, grid.wet_set.logλ.hi)
            μ = P3.get_μ(params.slope, logλ)
            @test grid.wet_set.μ.lo ≤ μ ≤ grid.wet_set.μ.hi
            @test grid.wet_set.logD.lo ≤ P3.collision_logξ_lo(μ) - logλ + 1e-5
            @test P3.collision_logξ_hi(μ) - logλ ≤ grid.wet_set.logD.hi + 1e-5
        end
        # The axes of each table group
        axes = P3.table_axes(grid)
        @test axes.ice == values(grid.ice)
        @test axes.cloud_collision == (values(grid.liquid_ice)[1:4]..., grid.liquid_ice.cloud_logx̄)
        @test length(axes.rain_particle_rime_density) == 5
        @test P3.P3TableGrid(ice; slope = (; logλ = 32, w = 4, ρ_rim = 3)).slope.ρ_rim.n == 3
        mp = CMP.Microphysics2MParams(CP.create_toml_dict(FT); with_ice = true)
        @test P3.P3TableGrid(mp) == grid

        @test_throws ArgumentError P3.P3TableGrid(CMP.P3IceParams(CP.create_toml_dict(FT); is_limited = false))
        @test_throws ArgumentError P3.P3TableGrid(CMP.Microphysics2MParams(CP.create_toml_dict(FT)))

        # The rime coordinate and its inverse
        @test P3.rime_coordinate(FT(0)) == 0 && P3.rime_coordinate(FT(1)) == 1
        @test P3.rime_mass_fraction_at(P3.rime_coordinate(0.3)) ≈ 0.3
        @test P3.rime_mass_fraction_at(1.0) < 1
    end
end

function test_inverse_slope(FT, ice, tables)
    params = ice.scheme
    params64 = CMP.ParametersP3(CP.create_toml_dict(Float64))
    (; groups) = tables

    @testset "P3Tables type" begin
        @test eltype(tables) == FT
        @test keys(groups) == (:slope,)
        @test groups.slope isa P3.UniformTable{FT, 3, 1}
        @test size(groups.slope.data) == (1, 128, 16, 16)
        @test tables.fingerprint == P3.table_fingerprint(ice)
        @test !isempty(sprint(show, MIME"text/plain"(), tables))

        # Adaptation reaches the tables inside the parameters and inside a state
        mp = CMP.Microphysics2MParams(CP.create_toml_dict(FT); with_ice = true, integration = tables)
        mp_view = Adapt.adapt(ViewAdaptor(), mp)
        @test mp_view.ice.integration.groups.slope.data isa SubArray
        @test mp_view.ice.integration.fingerprint == tables.fingerprint
        @test mp_view.ice.scheme == mp.ice.scheme
        st = P3.P3State(params, tables, FT(1e-4), FT(1e5), FT(3e-5), FT(1e-7))
        st_view = Adapt.adapt(ViewAdaptor(), st)
        @test st_view.integration.groups.slope.data isa SubArray
        @test st_view.dist == st.dist && st_view.particles == st.particles
        @test P3.P3State(params, mp_view.ice.integration, FT(1e-4), FT(1e5), FT(3e-5), FT(1e-7)).dist == st.dist
    end

    @testset "Inverse slope table against the solve" begin
        states = sample_ice_states(params64, 2000)
        tabulated(logx̄, F_rim, ρ_rim) =
            P3.tabulated_distribution_logλ(
                groups.slope,
                P3.P3Particles(params, FT(F_rim), FT(ρ_rim)),
                FT(exp(logx̄)),
                FT(1),
            )
        errors = [abs(tabulated(s...) - reference_logλ(params64, s...)) for s in states]
        @test maximum(errors) < 3e-3
        @test sort(errors)[end - 19] < 1e-3

        # At a table point the inversion before the Newton step is exact
        grid = P3.P3TableGrid(ice)
        logλs, ws, ρs = map(P3.axis_points, values(grid.slope))
        for (k, j, l) in ((5, 1, 1), (40, 7, 3), (100, 16, 16), (127, 12, 9))
            particles = P3.P3Particles(params, FT(P3.rime_mass_fraction_at(ws[j])), FT(ρs[l]))
            target = groups.slope.data[1, k, j, l]
            logλ = P3.tabulated_distribution_logλ(groups.slope, particles, exp(target), one(FT))
            @test logλ ≈ logλs[k] atol = FT === Float64 ? 1e-6 : 1e-3
        end
    end

    @testset "P3State under tables" begin
        ρq_ice, ρn_ice, ρq_rim, ρb_rim = FT(1e-4), FT(1e5), FT(3e-5), FT(1e-7)
        st = @inferred P3.P3State(params, tables, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
        @test (BT.@ballocated P3.P3State($params, $tables, $ρq_ice, $ρn_ice, $ρq_rim, $ρb_rim) samples = 1 evals = 1) ==
              0
        sq = P3.P3State(params, ice.integration, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
        @test st.integration === tables
        @test st.particles == sq.particles
        @test st.dist.ρn == sq.dist.ρn
        @test st.dist.logλ ≈ sq.dist.logλ atol = 3e-3
        @test st.dist.μ ≈ sq.dist.μ atol = 1e-2
        @test st.dist.logN₀ ≈ sq.dist.logN₀ rtol = 1e-3
        # Broadcasting over cells
        sts = P3.P3State.(params, tables, [ρq_ice, 2ρq_ice], ρn_ice, FT(0), FT(0))
        @test sts[2].dist.logλ < sts[1].dist.logλ

        # An absent population: the bounds of the solve and no number
        for (q, n) in ((FT(0), ρn_ice), (ρq_ice, FT(0)), (FT(0), FT(0)))
            st₀ = P3.P3State(params, tables, q, n, FT(0), FT(0))
            sq₀ = P3.P3State(params, ice.integration, q, n, FT(0), FT(0))
            @test st₀.dist.ρn == 0 && st₀.dist.logλ == sq₀.dist.logλ
            @test st₀.dist.logN₀ == -Inf
        end
        # Nonfinite input
        @test isfinite(P3.P3State(params, tables, FT(NaN), ρn_ice, FT(0), FT(0)).dist.logλ)

        # The shape is held under differentiation, and the number scales the intercept
        for rule in (tables, ice.integration)
            ∂logλ_∂q = FD.derivative(q -> P3.P3State(params, rule, q, ρn_ice, ρq_rim, ρb_rim).dist.logλ, ρq_ice)
            @test ∂logλ_∂q == 0
        end
        logN₀(rule, n) = P3.P3State(params, rule, ρq_ice, n, ρq_rim, ρb_rim).dist.logN₀
        @test FD.derivative(n -> logN₀(tables, n), ρn_ice) ≈ FD.derivative(n -> logN₀(ice.integration, n), ρn_ice)
    end
end

@testset "P3 tables ($FT)" for FT in (Float64, Float32)
    test_uniform_table(FT)
    test_table_grid(FT)
    ice = CMP.P3IceParams(CP.create_toml_dict(FT))
    grid = P3.P3TableGrid(ice; TEST_GRID...)
    tables = tables_of(ice, grid, (; slope = P3.SlopePoint(ice.scheme)))
    test_inverse_slope(FT, ice, tables)
end
nothing
