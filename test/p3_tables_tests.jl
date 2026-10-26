include("p3_tables_helpers.jl")

function test_uniform_table(FT)
    @testset "UniformTable lookup" begin
        axes = (P3.TableAxis(-1, 2, 4), P3.TableAxis(0, 1, 3), P3.TableAxis(10, 20, 5))
        f(x, y, z) = 1 + x - 2y + z / 10 + x * y * z / 100
        g(x, y, z) = exp(x - y + z / 20)
        point((x, y, z)) = P3.SA.SVector(f(x, y, z), g(x, y, z))
        t = P3.UniformTable(P3.tabulate(point, axes), axes, FT; logs = (false, true))
        @test t isa P3.UniformTable{FT, 3, 2} && size(t.data) == (2, 4, 3, 5)
        rtol = 10eps(FT)
        for (x, y, z) in ((-1, 0, 10), (2, 1, 20), (0.3, 0.7, 12.5), (1.9, 0.01, 19.9))
            (v_f, v_g) = P3.lookup(t, FT(x), FT(y), FT(z))
            @test v_f ≈ f(x, y, z) rtol = rtol
            @test exp(v_g) ≈ g(x, y, z) rtol = 10rtol
        end
        @test P3.lookup(t, FT(-5), FT(2), FT(Inf))[1] ≈ f(-1, 1, 20) rtol = rtol
        @test P3.lookup(t, FT(-Inf), FT(NaN), FT(NaN))[1] ≈ f(-1, 0, 10) rtol = rtol
        @test P3.uniform_index(FT(NaN), FT(0), FT(1), 5) == (1, 0)
        @test P3.uniform_index(FT(Inf), FT(0), FT(1), 5) == (4, 1)
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
        @test_throws DomainError P3.UniformTable(zeros(4, 3, 5), axes, FT; logs = (true,))
        @test_throws DomainError P3.UniformTable(fill(NaN, 4, 3, 5), axes, FT)
        @test all(n -> P3.bisection_steps(n) == ceil(Int, log2(n - 1)), 2:300)
        axes₁ = (P3.TableAxis(-1, 2, 4),)
        t₁ = P3.UniformTable(P3.tabulate(((x,),) -> 1 + 2x, axes₁), axes₁, FT)
        @test P3.lookup(t₁, FT(0.3))[1] ≈ FT(1.6) rtol = rtol
    end

    @testset "Table of a point function on an inner axis" begin
        inner = P3.TableAxis(0, 1, 5)
        axes = (P3.TableAxis(-1, 1, 3), P3.TableAxis(0, 2, 4))
        point((x, y)) = P3.SA.SVector(ntuple(k -> x * y + P3.axis_point(inner, k), 5))
        t = P3.UniformTable(P3.tabulate(point, axes), (inner, axes...), FT; inner = true)
        @test t isa P3.UniformTable{FT, 3, 1} && size(t.data) == (1, 5, 3, 4)
        @test P3.lookup(t, FT(0.3), FT(0.5), FT(1.5))[1] ≈ FT(0.3 + 0.75) rtol = 10eps(FT)
        @test_throws DimensionMismatch P3.UniformTable(
            P3.tabulate(point, axes),
            (P3.TableAxis(0, 1, 4), axes...),
            FT;
            inner = true,
        )
    end
end

function test_table_grid(FT)
    @testset "P3TableGrid bounds from the parameters" begin
        ice = CMP.P3IceParams(CP.create_toml_dict(FT))
        params = ice.scheme
        grid = P3.P3TableGrid(ice)
        (; slope) = grid.groups
        ice_axes = grid.groups.ice
        @test all(a -> a isa P3.TableAxis{FT}, Iterators.flatten(values(grid.groups)))
        @test (slope[1].lo, slope[1].hi) == P3.logλ_bounds(params)
        @test (slope[2].lo, slope[2].hi) == (0, 1)
        @test (ice_axes[3].lo, ice_axes[3].hi) == P3.rime_density_bounds(params)
        @test (ice_axes[1].lo, ice_axes[1].hi) == P3.integral_logλ_bounds(params)
        unrimed = P3.P3Particles(params, FT(0), params.ρ_i)
        @test exp(P3.logLdivN(unrimed, ice_axes[1].lo)) ≈ P3.ice_mean_particle_mass_max(params) rtol = 1e-3
        rimed = P3.P3Particles(params, FT(0.99), FT(500))
        @test P3.get_distribution_logλ(rimed, P3.ice_mean_particle_mass_max(params), FT(1)) > ice_axes[1].lo
        @test all(exp.((ice_axes[4].lo, ice_axes[4].hi)) .≈ (0.005, 1.7))
        @test P3.P3TableGrid(ice; slope = (; logλ = 32, w = 4, ρ_rim = 3)).groups.slope[3].n == 3
        @test_throws ArgumentError P3.P3TableGrid(ice; liquid_ice = TEST_ICE_AXES)
        mp = CMP.Microphysics2MParams(CP.create_toml_dict(FT); with_ice = true)
        @test P3.P3TableGrid(mp) == grid
        @test_throws ArgumentError P3.P3TableGrid(CMP.Microphysics2MParams(CP.create_toml_dict(FT)))

        @test P3.rime_coordinate(FT(0)) == 0 && P3.rime_coordinate(FT(1)) == 1
        @test P3.rime_mass_fraction_at(P3.rime_coordinate(0.3)) ≈ 0.3
        @test P3.rime_mass_fraction_at(1.0) < 1
    end
end

function test_table_construction(FT, ice, tables)
    grid = P3.P3TableGrid(ice; TEST_GRID...)
    groups = P3.table_groups(ice, grid, P3.GaussLegendre(FT, 16))

    @testset "Point functions of the table groups" begin
        @test keys(groups) == keys(grid.groups)
        for ((; point, inner), axes) in zip(groups, grid.groups)
            a = inner ? Base.tail(axes) : axes
            x_first = P3.grid_point(a, CartesianIndex(map(_ -> 1, a)))
            x_last = P3.grid_point(a, CartesianIndex(map(ax -> ax.n, a)))
            @test isbits(point)
            @test (BT.@ballocated $point($x_first) samples = 1 evals = 1) == 0
            @test (@inferred point(x_last)) isa Union{FT, P3.SA.SVector{<:Any, FT}}
            inner && @test length(point(x_first)) == first(axes).n
        end
    end

    @testset "Tables for the whole scheme" begin
        @test P3.tables_are_current(ice) && P3.tables_are_current(P3.with_integration(ice, tables))
        mp = CMP.Microphysics2MParams(CP.create_toml_dict(FT); with_ice = true)
        mp_tab = P3.with_tables(mp; grid = P3.P3TableGrid(mp; TINY_GRID...))
        @test mp_tab.ice.integration isa P3.P3Tables
        @test mp_tab.warm_rain == mp.warm_rain && mp_tab.ice.scheme == mp.ice.scheme
        @test P3.tables_are_current(mp_tab) && P3.tables_are_current(mp)
        @test P3.tables_are_current(CMP.Microphysics2MParams(CP.create_toml_dict(FT)))
        @test_throws ArgumentError P3.with_tables(CMP.Microphysics2MParams(CP.create_toml_dict(FT)))
    end

    @testset "Tables of other parameters" begin
        toml = CP.create_toml_dict(FT)
        @test CMP.P3IceParams(toml; integration = tables).integration === tables
        @test CMP.Microphysics2MParams(toml; with_ice = true, integration = tables).ice.integration === tables
        @test_throws ArgumentError CMP.P3IceParams(toml; is_limited = false, integration = tables)
        @test_throws ArgumentError CMP.P3IceParams(
            CP.create_toml_dict(Float64 === FT ? Float32 : Float64);
            integration = tables,
        )
        not_limited = CMP.P3IceParams(toml; is_limited = false)
        @test !P3.tables_are_current(P3.with_integration(not_limited, tables))
    end
end

function test_inverse_slope(FT, ice, tables)
    params = ice.scheme
    params64 = CMP.ParametersP3(CP.create_toml_dict(Float64))
    (; groups) = tables

    @testset "P3Tables type" begin
        @test eltype(tables) == FT
        @test keys(groups) == keys(P3.P3TableGrid(ice).groups)
        @test all(t -> t isa P3.UniformTable && eltype(t) == FT, groups)
        @test groups.slope isa P3.UniformTable{FT, 3, 1}
        @test size(groups.slope.data) == (1, 128, 16, 16)
        @test groups.ice isa P3.UniformTable{FT, 4, 4}
        @test size(groups.ice.data) == (4, values(TEST_ICE_AXES)...)
        @test !isempty(sprint(show, MIME"text/plain"(), tables))

        mp = CMP.Microphysics2MParams(CP.create_toml_dict(FT); with_ice = true, integration = tables)
        mp_view = Adapt.adapt(ViewAdaptor(), mp)
        @test mp_view.ice.integration.groups.slope.data isa SubArray
        @test mp_view.ice.integration.groups.ice.data isa SubArray
        @test mp_view.ice.integration.fingerprint == tables.fingerprint
        @test mp_view.ice.scheme == mp.ice.scheme
        st = P3.P3State(params, tables, FT(1e-4), FT(1e5), FT(3e-5), FT(1e-7))
        st_view = Adapt.adapt(ViewAdaptor(), st)
        @test st_view.integration.groups.ice.data isa SubArray
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

        grid = P3.P3TableGrid(ice)
        logλs, ws, ρs = map(P3.axis_points, grid.groups.slope)
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
        sts = P3.P3State.(params, tables, [ρq_ice, 2ρq_ice], ρn_ice, FT(0), FT(0))
        @test sts[2].dist.logλ < sts[1].dist.logλ

        for (q, n) in ((FT(0), ρn_ice), (ρq_ice, FT(0)), (FT(0), FT(0)))
            st₀ = P3.P3State(params, tables, q, n, FT(0), FT(0))
            sq₀ = P3.P3State(params, ice.integration, q, n, FT(0), FT(0))
            @test st₀.dist.ρn == 0 && st₀.dist.logλ == sq₀.dist.logλ
            @test st₀.dist.logN₀ == -Inf
        end
        @test isfinite(P3.P3State(params, tables, FT(NaN), ρn_ice, FT(0), FT(0)).dist.logλ)

        for rule in (tables, ice.integration)
            ∂logλ_∂q = FD.derivative(q -> P3.P3State(params, rule, q, ρn_ice, ρq_rim, ρb_rim).dist.logλ, ρq_ice)
            @test ∂logλ_∂q == 0
        end
        logN₀(rule, n) = P3.P3State(params, rule, ρq_ice, n, ρq_rim, ρb_rim).dist.logN₀
        @test FD.derivative(n -> logN₀(tables, n), ρn_ice) ≈ FD.derivative(n -> logN₀(ice.integration, n), ρn_ice)
    end
end

# The elasticities `∂ log y / ∂ log x` of the outputs `f(x)` at `x`, with zero for a zero output
elasticities(f, x) = FD.jacobian(f, x) .* x' ./ map(y -> iszero(y) ? one(y) : y, f(x))

function test_ice_tables(FT, ice, tables)
    params = ice.scheme
    vel = ice.terminal_velocity
    aps = CMP.AirProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    ρq_ice, ρn_ice, ρq_rim, ρb_rim, ρₐ, T = FT(1e-4), FT(1e5), FT(3e-5), FT(1e-7), FT(0.8), FT(275.15)

    @testset "Ice tables against the quadrature" begin
        test_ice_table_accuracy(ice, tables; n_states = 400, median = 0.01, p99 = 0.1)
    end

    @testset "Ice tables at the grid points" begin
        points = map(P3.axis_points, P3.P3TableGrid(ice; ice = TEST_ICE_AXES).groups.ice)
        rtol = FT === Float64 ? 1e-10 : 1e-4
        for (i, j, k, l) in ((1, 1, 1, 1), (20, 4, 2, 3), (40, 9, 5, 4), (63, 7, 3, 2))
            logλ, w, ρ_rim, logρₐ = points[1][i], points[2][j], points[3][k], points[4][l]
            st = P3.table_state(params, P3.GaussLegendre(FT, 16), logλ, w, ρ_rim)
            (M_v, self_collection, v_n, v_m) = exp.(P3.lookup(tables.groups.ice, FT.((logλ, w, ρ_rim, logρₐ))...))
            @test M_v ≈ P3.ventilation_moment(vel, st, exp(logρₐ)) rtol = rtol
            @test self_collection ≈ P3.ice_self_collection(vel, st, exp(logρₐ)).dNdt rtol = rtol
            @test P3.node_self_collection(vel, st.integration, st, exp(logρₐ)) ≈
                  P3.ice_self_collection(vel, st, exp(logρₐ)).dNdt rtol = FT === Float64 ? 1e-12 : 100eps(FT)
            velocities = P3.ice_terminal_velocities(vel, st, exp(logρₐ))
            @test v_n ≈ velocities.v_n rtol = rtol
            @test v_m ≈ velocities.v_m rtol = rtol
        end
    end

    @testset "Ventilation from its two moments" begin
        st = P3.P3State(params, ice.integration, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
        (a, b) = CO.ventilation_coefficients(params.vent, aps)
        (; ρn, μ, logλ) = st.dist
        @test a * ρn * (μ + 1) * exp(-logλ) + b * P3.ventilation_moment(vel, st, ρₐ) ≈
              P3.ice_ventilation_integral(vel, aps, st, ρₐ) rtol = FT === Float64 ? 1e-5 : 1e-4
    end

    @testset "Absent ice gives zero under tables" begin
        for (q, n) in ((FT(0), ρn_ice), (ρq_ice, FT(0)), (FT(0), FT(0)))
            st = P3.P3State(params, tables, q, n, FT(0), FT(0))
            @test P3.ice_terminal_velocities(vel, st, ρₐ) == (; v_n = 0, v_m = 0)
            @test P3.ice_self_collection(vel, st, ρₐ).dNdt == 0
            @test P3.ice_ventilation_integral(vel, aps, st, ρₐ) == 0
            @test P3.ice_melt(vel, aps, tps, st, ρₐ, T) == (; dNdt = 0, dLdt = 0)
        end
    end

    @testset "The table path evaluates no quadrature node" begin
        throwing =
            (; small_ice = ThrowingVelocity(vel.small_ice.cutoff), large_ice = ThrowingVelocity(vel.large_ice.cutoff))
        sq = P3.P3State(params, ice.integration, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
        @test_throws Exception P3.ice_terminal_velocities(throwing, sq, ρₐ)
        st = P3.P3State(params, tables, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
        @test P3.ice_terminal_velocities(throwing, st, ρₐ) == P3.ice_terminal_velocities(vel, st, ρₐ)
        @test P3.ice_self_collection(throwing, st, ρₐ) == P3.ice_self_collection(vel, st, ρₐ)
        @test P3.ice_ventilation_integral(throwing, aps, st, ρₐ) == P3.ice_ventilation_integral(vel, aps, st, ρₐ)
        @test P3.ice_melt(throwing, aps, tps, st, ρₐ, T) == P3.ice_melt(vel, aps, tps, st, ρₐ, T)
    end

    @testset "Liquid-ice collisions of a tabulated state by quadrature" begin
        (; cloud_pdf, rain_pdf) = ice
        st = P3.P3State(params, tables, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
        sq = P3.P3State(st.particles, ice.integration, st.ρq_ice, st.ρn_ice, st.ρb_rim, st.dist)
        ρw = FT(1000)
        mass = D -> ρw * CO.volume_sphere_D(D)
        args = (FT(2e-4), FT(1e8), FT(1e-4), FT(1e4), ρₐ, FT(265))
        @test P3.∫liquid_ice_collisions(cloud_pdf, rain_pdf, aps, tps, vel, mass, st, args...) ==
              P3.∫liquid_ice_collisions(cloud_pdf, rain_pdf, aps, tps, vel, mass, sq, args...)
    end

    @testset "Inference and allocations under tables" begin
        st = P3.P3State(params, tables, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
        @test (@inferred P3.ice_terminal_velocities(vel, st, ρₐ)) isa @NamedTuple{v_n::FT, v_m::FT}
        @test (@inferred P3.ice_self_collection(vel, st, ρₐ)).dNdt isa FT
        @test (@inferred P3.ice_ventilation_integral(vel, aps, st, ρₐ)) isa FT
        @test (@inferred P3.ice_melt(vel, aps, tps, st, ρₐ, T)).dLdt isa FT
        @test (BT.@ballocated P3.ice_terminal_velocities($vel, $st, $ρₐ) samples = 1 evals = 1) == 0
        @test (BT.@ballocated P3.ice_self_collection($vel, $st, $ρₐ) samples = 1 evals = 1) == 0
        @test (BT.@ballocated P3.ice_melt($vel, $aps, $tps, $st, $ρₐ, $T) samples = 1 evals = 1) == 0
    end

    @testset "Derivatives of the ice integrals under tables" begin
        x = [ρq_ice, ρn_ice, ρq_rim, ρb_rim, ρₐ]
        integrals(rule) =
            x -> begin
                st = P3.P3State(params, rule, x[1:4]...)
                [
                    P3.ice_terminal_velocities(vel, st, x[5]).v_m,
                    P3.ice_self_collection(vel, st, x[5]).dNdt,
                    P3.ice_melt(vel, aps, tps, st, x[5], T).dLdt,
                ]
            end
        E_tab = elasticities(integrals(tables), x)
        E_quad = elasticities(integrals(ice.integration), x)
        @test all(isfinite, E_tab)
        @test E_tab[:, 2] ≈ [0, 2, 0] atol = 10eps(FT)
        @test maximum(abs, E_tab - E_quad) < 0.1
    end
end

@testset "P3 tables ($FT)" for FT in (Float64, Float32)
    test_uniform_table(FT)
    test_table_grid(FT)
    ice = CMP.P3IceParams(CP.create_toml_dict(FT))
    # The comparisons with the quadrature run in Float32 only, see #834
    grid = FT === Float32 ? TEST_GRID : TINY_GRID
    tables = P3.with_tables(ice; grid = P3.P3TableGrid(ice; grid...)).integration
    test_table_construction(FT, ice, tables)
    if FT === Float32
        test_inverse_slope(FT, ice, tables)
        test_ice_tables(FT, ice, tables)
    end
end
nothing
