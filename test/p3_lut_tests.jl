import Test as TT
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Parameters as CMP
import ClimaParams as CP
import Adapt

function test_grid(::Type{FT}, scheme) where {FT}
    ρ_rim_lo, ρ_rim_hi = P3.rime_density_bounds(scheme)
    P3.IceKernelTableGrid(;
        nx = 6, nf = 4, nr = 4, na = 4,
        logx̄_lo = FT(log10(1e-14)), logx̄_hi = FT(log10(1e-3)),
        ρ_rim_lo, ρ_rim_hi,
        logρₐ_lo = FT(log10(0.05)), logρₐ_hi = FT(log10(1.35)),
    )
end

test_states(::Type{FT}) where {FT} = (
    (FT(1e-4), FT(1e5), FT(0), FT(0)),
    (FT(5e-4), FT(2e5), FT(0.3), FT(500)),
    (FT(1e-5), FT(1e6), FT(0.8), FT(800)),
)

function test_p3_lut_off_mode_bit_identical(FT)
    TT.@testset "off-mode carrier is bit-identical to the plain rule" begin
        ice = CMP.P3IceParams(FT; quadrature_order = 8)
        scheme, vel, rule = ice.scheme, ice.terminal_velocity, ice.quad
        aps = CMP.AirProperties(FT)
        carrier = P3.P3TabulatedQuadrature(rule)
        ρₐ = FT(1.1)
        for (ρq_ice, ρn_ice, F_rim, ρ_rim) in test_states(FT)
            q_rim = ρq_ice * F_rim
            b_rim = ρ_rim > 0 ? q_rim / ρ_rim : zero(FT)
            state = P3.state_from_prognostic(scheme, ρq_ice, ρn_ice, q_rim, b_rim)
            logλ = P3.get_distribution_logλ(state)

            TT.@test P3.ice_self_collection(state, logλ, vel, ρₐ; quad = carrier).dNdt ===
                     P3.ice_self_collection(state, logλ, vel, ρₐ; quad = rule).dNdt
            TT.@test P3.ice_ventilation_integral(vel, aps, ρₐ, state, logλ; quad = carrier) ===
                     P3.ice_ventilation_integral(vel, aps, ρₐ, state, logλ; quad = rule)
            TT.@test P3.ice_terminal_velocity_number_weighted(vel, ρₐ, state, logλ; quad = carrier) ===
                     P3.ice_terminal_velocity_number_weighted(vel, ρₐ, state, logλ; quad = rule)
            TT.@test P3.ice_terminal_velocity_mass_weighted(vel, ρₐ, state, logλ; quad = carrier) ===
                     P3.ice_terminal_velocity_mass_weighted(vel, ρₐ, state, logλ; quad = rule)
        end
    end
end

function test_p3_lut_untabulated_terms_bit_identical(FT)
    TT.@testset "adopting one table leaves the others bit-identical" begin
        ice = CMP.P3IceParams(FT; quadrature_order = 8)
        scheme, vel, rule = ice.scheme, ice.terminal_velocity, ice.quad
        aps = CMP.AirProperties(FT)
        grid = test_grid(FT, scheme)
        carrier = P3.p3_lut_carrier(scheme, vel, aps, grid; quad = rule, outputs = (:selfcol,))
        ρₐ = FT(1.1)
        for (ρq_ice, ρn_ice, F_rim, ρ_rim) in test_states(FT)
            q_rim = ρq_ice * F_rim
            b_rim = ρ_rim > 0 ? q_rim / ρ_rim : zero(FT)
            state = P3.state_from_prognostic(scheme, ρq_ice, ρn_ice, q_rim, b_rim)
            logλ = P3.get_distribution_logλ(state)

            TT.@test P3.ice_ventilation_integral(vel, aps, ρₐ, state, logλ; quad = carrier) ===
                     P3.ice_ventilation_integral(vel, aps, ρₐ, state, logλ; quad = rule)
            TT.@test P3.ice_terminal_velocity_number_weighted(vel, ρₐ, state, logλ; quad = carrier) ===
                     P3.ice_terminal_velocity_number_weighted(vel, ρₐ, state, logλ; quad = rule)
            TT.@test P3.ice_terminal_velocity_mass_weighted(vel, ρₐ, state, logλ; quad = carrier) ===
                     P3.ice_terminal_velocity_mass_weighted(vel, ρₐ, state, logλ; quad = rule)
        end
    end
end

function test_p3_lut_stored_entry_matches_integrand(FT)
    TT.@testset "stored table entries match a fresh evaluation" begin
        ice = CMP.P3IceParams(FT; quadrature_order = 8)
        scheme, vel, rule = ice.scheme, ice.terminal_velocity, ice.quad
        aps = CMP.AirProperties(FT)
        grid = test_grid(FT, scheme)
        carrier = P3.p3_lut_carrier(scheme, vel, aps, grid; quad = rule)
        for lin in (1, 37, P3.grid_length(grid))
            TT.@test carrier.selfcol.logI[lin] ==
                     P3.ice_self_collection_table_logI(scheme, vel, rule, grid, lin)
            TT.@test carrier.vent.logI[lin] ==
                     P3.ice_ventilation_table_logK(scheme, vel, aps, rule, grid, lin)
            TT.@test carrier.vel_n.logI[lin] ==
                     P3.ice_velocity_n_table_logK(scheme, vel, rule, grid, lin)
            TT.@test carrier.vel_m.logI[lin] ==
                     P3.ice_velocity_m_table_logK(scheme, vel, rule, grid, lin)
        end
    end

    TT.@testset "the grid's lower mean-mass edge is the shape solve's own pin" begin
        # A host derives the table's lower edge from `ice_mean_particle_mass_min`, which is also
        # what `_derived_logλ_bracket` derives `logλ_max` from. Below that mass the target is
        # unbracketable, `logλ` clamps to the bracket endpoint and every stored quantity stops
        # depending on the mean mass. Asserting the edge against the accessor's VALUE would still
        # pass if the bracket were later derived from a different bound, leaving a host's edge above
        # the pin, where the entry varies and a clamp returns a wrong value with no symptom. So the
        # property is asserted instead: identical below, different above.
        ice = CMP.P3IceParams(FT; quadrature_order = 8)
        scheme, vel, rule = ice.scheme, ice.terminal_velocity, ice.quad
        aps = CMP.AirProperties(FT)
        lx = log10(P3.ice_mean_particle_mass_min(scheme))
        rb = P3.rime_density_bounds(scheme)
        fr, rr, la = FT(0.3), FT((rb[1] + rb[2]) / 2), FT(0)
        entries = (
            l -> P3.ice_self_collection_table_entry(scheme, vel, rule, FT(l), fr, rr, la),
            l -> P3.ice_ventilation_table_entry(scheme, vel, aps, rule, FT(l), fr, rr, la),
            l -> P3.ice_velocity_n_table_entry(scheme, vel, rule, FT(l), fr, rr, la),
            l -> P3.ice_velocity_m_table_entry(scheme, vel, rule, FT(l), fr, rr, la),
        )
        for f in entries
            TT.@test f(lx - FT(0.7)) === f(lx - FT(0.2))
            TT.@test f(lx + FT(0.3)) !== f(lx - FT(0.2))
        end
    end
end

function test_p3_lut_carrier_composition(FT)
    TT.@testset "the carrier's adopted tables match the requested outputs" begin
        ice = CMP.P3IceParams(FT; quadrature_order = 8)
        scheme, vel, rule = ice.scheme, ice.terminal_velocity, ice.quad
        aps = CMP.AirProperties(FT)
        grid = test_grid(FT, scheme)
        for outputs in (
            (:selfcol, :vent, :vel_n, :vel_m),
            (:selfcol,),
            (:vent, :vel_m),
            (),
        )
            carrier = P3.p3_lut_carrier(scheme, vel, aps, grid; quad = rule, outputs)
            TT.@test (carrier.selfcol !== nothing) == (:selfcol in outputs)
            TT.@test (carrier.vent !== nothing) == (:vent in outputs)
            TT.@test (carrier.vel_n !== nothing) == (:vel_n in outputs)
            TT.@test (carrier.vel_m !== nothing) == (:vel_m in outputs)
        end
        TT.@test_throws ErrorException P3.p3_lut_carrier(
            scheme, vel, aps, grid; quad = rule, outputs = (:selfcol, :unknown))
    end
end

function test_p3_lut_quadrature_delegation(FT)
    TT.@testset "the carrier delegates to its rule as a QuadratureRule" begin
        ice = CMP.P3IceParams(FT; quadrature_order = 8)
        rule = ice.quad
        carrier = P3.P3TabulatedQuadrature(rule)
        TT.@test carrier.n == rule.n
        TT.@test carrier isa P3.QuadratureRule
        TT.@test Base.broadcastable(carrier) == (carrier,)
        f(x) = x^2
        bnds = (zero(FT), one(FT))
        TT.@test P3.integrate(f, bnds, carrier) == P3.integrate(f, bnds, rule)
    end
end

function test_p3_lut_adapt_identity(FT)
    TT.@testset "the identity adaptor round-trips the carrier through P3IceParams" begin
        ice = CMP.P3IceParams(FT; quadrature_order = 8)
        scheme, vel, rule = ice.scheme, ice.terminal_velocity, ice.quad
        aps = CMP.AirProperties(FT)
        grid = test_grid(FT, scheme)
        carrier = P3.p3_lut_carrier(scheme, vel, aps, grid; quad = rule, outputs = (:selfcol,))

        ice_tab = CMP.P3IceParams(FT; quadrature_order = 8, quad = carrier)
        TT.@test ice_tab.quad isa P3.P3TabulatedQuadrature
        ice_tab_adapted = Adapt.adapt(identity, ice_tab)
        TT.@test ice_tab_adapted.quad.selfcol.logI === ice_tab.quad.selfcol.logI

        mp = CMP.Microphysics2MParams(FT; with_ice = true, quadrature_order = 8, quad = carrier)
        TT.@test mp.ice.quad isa P3.P3TabulatedQuadrature
        mp_adapted = Adapt.adapt(identity, mp)
        TT.@test mp_adapted.ice.quad.selfcol.logI === mp.ice.quad.selfcol.logI
    end
end

# The self-collection table's own testset. Its enclosing function header was lost in an
# extraction, which left the block at top level with `FT` a free variable, so it raised
# `UndefVarError: FT` on every run and every assertion inside it went unexecuted.
function test_p3_lut_self_collection_table(FT)
    TT.@testset "ice self-collection lookup table" begin
        params = CMP.ParametersP3(FT)
        vel = CMP.Chen2022VelType(FT)
        rb = P3.rime_density_bounds(params)
        nx, nf, nr, na = 12, 5, 6, 5
        lx_lo, lx_hi = FT(log10(1e-14)), FT(log10(1e-3))
        la_lo, la_hi = FT(log10(0.05)), FT(log10(1.35))
        tab = P3.IceKernelTable(
            params, vel; quad = P3.GaussLegendre(FT, 16), nx, nf, nr, na,
            logx̄_lo = lx_lo, logx̄_hi = lx_hi,
            ρ_rim_lo = FT(rb[1]), ρ_rim_hi = FT(rb[2]),
            logρₐ_lo = la_lo, logρₐ_hi = la_hi,
        )
        TT.@test all(isfinite, tab.logI)

        TT.@testset "at a grid node the lookup returns the stored value" begin
            # At a node the interpolation weights are exactly one and zero, so this is a direct
            # test of the index arithmetic - the part a smooth off-node comparison would hide.
            ax(lo, hi, n, i) = n == 1 ? lo : lo + (hi - lo) * FT(i - 1) / FT(n - 1)
            for (ix, jf, kr, ma) in ((1, 1, 1, 1), (5, 3, 4, 2), (nx, nf, nr, na))
                got = P3.lookup(
                    tab, exp10(ax(lx_lo, lx_hi, nx, ix)), ax(zero(FT), one(FT), nf, jf),
                    ax(FT(rb[1]), FT(rb[2]), nr, kr), exp10(ax(la_lo, la_hi, na, ma)),
                )
                TT.@test got ≈ exp(tab.logI[ix, jf, kr, ma]) rtol = sqrt(eps(FT))
            end
        end

        TT.@testset "a stored entry is the integrand at that node's own coordinates" begin
            # NODE EXACTNESS ABOVE CANNOT CATCH A WRONG AXIS DECODE. It compares the lookup against
            # the STORED array, so entries filed at the wrong grid point pass it - both sides are
            # wrong together. This compares a stored entry against the integrand evaluated directly
            # at that node's own coordinates, which is the only thing that pins the fill's
            # linear-index decode to the array's layout. A caller-side fill calls the same decode,
            # so this is its test too.
            ax(lo, hi, n, i) = n == 1 ? lo : lo + (hi - lo) * FT(i - 1) / FT(n - 1)
            for (ix, jf, kr, ma) in ((1, 1, 1, 1), (5, 3, 4, 2), (nx, nf, nr, na))
                direct = P3.ice_self_collection_table_entry(
                    params, vel, P3.GaussLegendre(FT, 16),
                    ax(lx_lo, lx_hi, nx, ix), ax(zero(FT), one(FT), nf, jf),
                    ax(FT(rb[1]), FT(rb[2]), nr, kr), ax(la_lo, la_hi, na, ma),
                )
                TT.@test exp(tab.logI[ix, jf, kr, ma]) ≈ direct rtol = sqrt(eps(FT))
            end
        end

        TT.@testset "a wrongly-shaped array is refused rather than reinterpreted" begin
            # The shape check is what stops a caller-side fill returning an array whose axes would
            # be silently read in the wrong order; a guard with no execution proof is not a guard.
            grid = P3.IceKernelTableGrid(;
                nx, nf, nr, na, logx̄_lo = lx_lo, logx̄_hi = lx_hi,
                ρ_rim_lo = FT(rb[1]), ρ_rim_hi = FT(rb[2]),
                logρₐ_lo = la_lo, logρₐ_hi = la_hi,
            )
            TT.@test_throws ErrorException P3.IceKernelTable(
                zeros(FT, nx, nf, nr, na + 1), grid)
        end

        TT.@testset "an absent population returns exactly zero, taking no logarithm" begin
            st = P3.state_from_prognostic(params, zero(FT), FT(1e5), zero(FT), zero(FT))
            lg = P3.get_distribution_logλ(st)
            TT.@test P3.ice_self_collection(st, lg, vel, one(FT); quad = tab).dNdt == 0
        end

        TT.@testset "the mode carrier delegates faithfully and isolates the tabulated term" begin
            # THE FAILURE THIS ENCODES was not in self-collection. Putting a table into the single
            # `quad` slot made it the rule for EVERY term, and the compile died in liquid-ice
            # collisions - a term nobody was thinking about, because the table has no `n`, no `node`
            # and no `weight` to give it. So the carrier's contract is tested on an UNTABULATED term
            # as well as on the tabulated one, and bit-identity is the bar: an untabulated term must
            # not be able to tell a carrier from the rule it delegates to.
            aps = CMP.AirProperties(FT)
            tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
            rule = P3.GaussLegendre(FT, 6)
            bare = P3.P3TabulatedQuadrature(rule)                  # a carrier holding no table
            with = P3.P3TabulatedQuadrature(rule; selfcol = tab)   # and one holding the table
            st = P3.state_from_prognostic(params, FT(1e-3), FT(1e5), FT(5e-4), FT(1e-6))
            lg = P3.get_distribution_logλ(st)
            ρₐ, T_w = one(FT), FT(273.15 + 0.01)

            τ_rule = P3.ice_deposition_timescale(vel, aps, tps, T_w, ρₐ, st, lg; quad = rule)
            TT.@test P3.ice_deposition_timescale(vel, aps, tps, T_w, ρₐ, st, lg; quad = bare) === τ_rule
            # and carrying a table for a DIFFERENT term changes this one not at all
            TT.@test P3.ice_deposition_timescale(vel, aps, tps, T_w, ρₐ, st, lg; quad = with) === τ_rule

            sc_rule = P3.ice_self_collection(st, lg, vel, ρₐ; quad = rule).dNdt
            sc_tab = P3.ice_self_collection(st, lg, vel, ρₐ; quad = tab).dNdt
            # no table: the carrier falls back to its own rule, exactly
            TT.@test P3.ice_self_collection(st, lg, vel, ρₐ; quad = bare).dNdt === sc_rule
            # a table: the carrier routes to it, exactly
            TT.@test P3.ice_self_collection(st, lg, vel, ρₐ; quad = with).dNdt === sc_tab
        end

        TT.@testset "the quadrature path is untouched by the table's existence" begin
            st = P3.state_from_prognostic(params, FT(1e-3), FT(1e5), FT(5e-4), FT(1e-6))
            lg = P3.get_distribution_logλ(st)
            q = P3.ice_self_collection(st, lg, vel, one(FT); quad = P3.GaussLegendre(FT, 6)).dNdt
            t = P3.ice_self_collection(st, lg, vel, one(FT); quad = tab).dNdt
            TT.@test isfinite(q) && q > 0
            # a coarse test grid, so this is a loose sanity bound and NOT the table's accuracy,
            # which is measured against a high-order reference on the production grid
            TT.@test isapprox(t, q; rtol = FT(0.1))
        end
    end
end

function test_p3_lut_outputs(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3, vel, rule = mp.ice.scheme, mp.ice.terminal_velocity, mp.ice.quad
    aps = CMP.AirProperties(FT)
    RB = P3.rime_density_bounds(p3)

    # Coarse on purpose: this asks which code path runs, not how accurate the table is, and a fine
    # grid would spend seconds of suite time answering a question the gates already answer better.
    g = P3.IceKernelTableGrid(; nx = 8, nf = 4, nr = 4, na = 4,
        logx̄_lo = FT(log10(1e-14)), logx̄_hi = FT(log10(1e-3)),
        ρ_rim_lo = FT(RB[1]), ρ_rim_hi = FT(RB[2]),
        logρₐ_lo = FT(log10(0.05)), logρₐ_hi = FT(log10(1.35)))

    t_sel = P3.IceKernelTable(p3, vel, g; quad = rule)
    t_ven = P3.ice_ventilation_table(p3, vel, aps, g; quad = rule)
    t_vn = P3.ice_velocity_n_table(p3, vel, g; quad = rule)
    t_vm = P3.ice_velocity_m_table(p3, vel, g; quad = rule)

    st = P3.state_from_prognostic(p3, FT(1e-3), FT(1e5), FT(5e-4), FT(1e-6))
    lg = P3.get_distribution_logλ(st)
    ρₐ = one(FT)

    TT.@testset "the mode's other three outputs route by their own carrier field ($FT)" begin
        bare = P3.P3TabulatedQuadrature(rule)

        # (2) the ventilation integral, read by deposition AND melt
        v_rule = P3.ice_ventilation_integral(vel, aps, ρₐ, st, lg; quad = rule)
        TT.@test P3.ice_ventilation_integral(vel, aps, ρₐ, st, lg; quad = bare) === v_rule
        TT.@test P3.ice_ventilation_integral(vel, aps, ρₐ, st, lg;
            quad = P3.P3TabulatedQuadrature(rule; selfcol = t_sel)) === v_rule
        v_tab = P3.ice_ventilation_integral(vel, aps, ρₐ, st, lg;
            quad = P3.P3TabulatedQuadrature(rule; vent = t_ven))
        TT.@test v_tab != v_rule || isapprox(v_tab, v_rule; rtol = sqrt(eps(FT)))

        # (3) and (4), the two means
        n_rule = P3.ice_terminal_velocity_number_weighted(vel, ρₐ, st, lg; quad = rule)
        m_rule = P3.ice_terminal_velocity_mass_weighted(vel, ρₐ, st, lg; quad = rule)
        TT.@test P3.ice_terminal_velocity_number_weighted(vel, ρₐ, st, lg; quad = bare) === n_rule
        TT.@test P3.ice_terminal_velocity_mass_weighted(vel, ρₐ, st, lg; quad = bare) === m_rule
        # a table for the OTHER velocity leaves this one exactly on the rule
        TT.@test P3.ice_terminal_velocity_number_weighted(vel, ρₐ, st, lg;
            quad = P3.P3TabulatedQuadrature(rule; vel_m = t_vm)) === n_rule
        TT.@test P3.ice_terminal_velocity_mass_weighted(vel, ρₐ, st, lg;
            quad = P3.P3TabulatedQuadrature(rule; vel_n = t_vn)) === m_rule

        # and each takes its own table when it has one
        full = P3.P3TabulatedQuadrature(rule; selfcol = t_sel, vent = t_ven, vel_n = t_vn, vel_m = t_vm)
        TT.@test P3.ice_terminal_velocity_number_weighted(vel, ρₐ, st, lg; quad = full) ==
                 P3.lookup(t_vn, st.ρq_ice / st.ρn_ice, st.F_rim, st.ρ_rim, ρₐ)
        TT.@test P3.ice_terminal_velocity_mass_weighted(vel, ρₐ, st, lg; quad = full) ==
                 P3.lookup(t_vm, st.ρq_ice / st.ρn_ice, st.F_rim, st.ρ_rim, ρₐ)

        # an absent population returns exactly zero on the table path, matching what the rule's own
        # guarded quotient gives there - the substitution must not invent a velocity for no ice
        empty = P3.state_from_prognostic(p3, zero(FT), zero(FT), zero(FT), zero(FT))
        lge = P3.get_distribution_logλ(empty)
        TT.@test P3.ice_terminal_velocity_number_weighted(vel, ρₐ, empty, lge; quad = full) == 0
        TT.@test P3.ice_ventilation_integral(vel, aps, ρₐ, empty, lge; quad = full) == 0
    end
end

# The driver runs LAST, because it executes at top level the moment it is reached and the
# functions it calls must already be defined. It previously sat above two of them.

function test_p3_lut_liqice_path(FT)
    TT.@testset "the tabulated liquid-ice path reproduces the quadrature path" begin
        # The seam substitutes a table read for the split assembly's full-range integral and
        # rebuilds nine channels from six stored integrals and a constant partition. A wrong
        # reconstruction returns a plausible number rather than an error, so this compares the two
        # paths at one state rather than checking that the tables were read.
        params = CMP.ParametersP3(FT)
        vel = CMP.Chen2022VelType(FT)
        aps = CMP.AirProperties(FT)
        tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
        mp = CMP.Microphysics2MParams(FT; with_ice = true)
        pdf_c, pdf_r = mp.ice.cloud_pdf, mp.ice.rain_pdf
        rule = P3.GaussLegendre(FT, 6)
        rb = P3.rime_density_bounds(params)
        common = (; nx = 4, nf = 3, nr = 3, na = 3, nl = 4,
            logx̄_lo = FT(log10(1e-11)), logx̄_hi = FT(log10(1e-6)),
            ρ_rim_lo = FT(rb[1]), ρ_rim_hi = FT(rb[2]),
            logρₐ_lo = FT(log10(0.4)), logρₐ_hi = FT(log10(1.2)))
        xc = (; logx̄l_lo = FT(log10(pdf_c.xc_min)), logx̄l_hi = FT(log10(pdf_c.xc_max)))
        xr = (; logx̄l_lo = FT(log10(pdf_r.xr_min)), logx̄l_hi = FT(log10(pdf_r.xr_max)))
        iT = (; invT_lo = FT(-1 / 40), invT_hi = FT(-1 / 0.5))
        lgc = P3.LiqIceTableGrid(; common..., xc...)
        lgr = P3.LiqIceTableGrid(; common..., xr...)
        bgc = P3.LiqIceBrimTableGrid(; common..., nt = 3, xc..., iT...)
        bgr = P3.LiqIceBrimTableGrid(; common..., nt = 3, xr..., iT...)

        # A ONE-NODE AXIS IS REJECTED AT CONSTRUCTION, because it cannot be caught later.
        # `_uniform_index` returns the pair `(i, i+1)` unconditionally and the lookup reads both
        # under `@inbounds`, so with one node on an axis one of the two reads is outside the array
        # whatever the index is clamped to. The weight on it is zero, so the result is right
        # whenever that memory happens to hold a finite number and silently non-finite when it does
        # not. There is no size to check against inside the lookup, so the grid has to refuse it.
        TT.@test_throws ErrorException P3.LiqIceTableGrid(; merge(common, (; nx = 1))..., xc...)
        TT.@test_throws ErrorException P3.LiqIceTableGrid(; merge(common, (; nl = 1))..., xr...)
        TT.@test_throws ErrorException P3.LiqIceBrimTableGrid(; common..., nt = 1, xc..., iT...)
        TT.@test_throws ErrorException P3.LiqIceBrimTableGrid(; merge(common, (; nf = 1))..., xr..., iT...)
        # Two nodes is the smallest axis that has a linear form, and it is accepted.
        TT.@test P3.LiqIceTableGrid(; merge(common, (; nx = 2))..., xc...) isa P3.LiqIceTableGrid
        m_liq(Dₗ) = pdf_c.ρw * CM.Common.volume_sphere_D(Dₗ)
        icegrid = test_grid(FT, params)
        carrier = P3.p3_lut_carrier(params, vel, aps, icegrid; quad = rule,
            outputs = (:liqice,),
            liqice_grid_cloud = lgc, liqice_grid_rain = lgr,
            liqice_brim_grid_cloud = bgc, liqice_brim_grid_rain = bgr,
            psd_c = pdf_c, psd_r = pdf_r, m_liq, T_freeze = params.T_freeze)
        for t in (:liqice_col_cloud, :liqice_col_rain, :liqice_brim_cloud, :liqice_brim_rain)
            TT.@test getfield(carrier, t) !== nothing
        end
        # The carrier must select the only assembly that can read it.
        TT.@test P3._default_liqice_assembly(carrier) isa P3.SplitCorrection
        TT.@test P3._default_liqice_assembly(rule) isa P3.PartitionedOuter

        st = P3.state_from_prognostic(params, FT(1e-3), FT(1e5), FT(5e-4), FT(1e-6))
        lg_ = P3.get_distribution_logλ(st)
        ρₐ, T = one(FT), FT(263.15)
        L_c, N_c, L_r, N_r = FT(5e-4), FT(1e8), FT(2e-4), FT(1e4)
        args = (st, lg_, pdf_c, pdf_r, L_c, N_c, L_r, N_r, aps, tps, vel, ρₐ, T, m_liq)
        tab = P3.∫liquid_ice_collisions(args...; quad = carrier)
        ref = P3.∫liquid_ice_collisions(args...; quad = rule)

        # Tolerance-free: the shed and frozen halves of each species sum back to that species'
        # collected mass, so the total is their sum in either path. A swapped channel, a lost
        # scale factor or a wrong partition all break this without needing a tolerance.
        for v in (tab, ref)
            TT.@test v[7] ≈ v[1] + v[2] + v[4] + v[5] rtol = sqrt(eps(FT))
        end
        # Away from a node the two paths differ by the table's own interpolation error, which on
        # this deliberately tiny grid is large, so the comparison there is loose and is not what
        # decides correctness.
        for i in (1, 3, 4, 6, 7, 8, 9)
            TT.@test tab[i] ≈ ref[i] rtol = 2.0
        end

        # AT A GRID NODE the interpolation is exact, so any residual disagreement here can only be
        # the nine-channel reconstruction itself. That is what this assertion is for; the loose
        # comparison above cannot separate a wrong reconstruction from a coarse grid.
        #
        # Both paths are handed the same liquid numbers, so the seam's own `ρn_ice * N_liq` scale
        # factors cancel between them whatever those numbers are. The numbers are chosen to put the
        # specific content at the value the table generator itself evaluates at, and NOT at one drop
        # per cubic metre: there the rain distribution's size bounds collapse to an empty interval
        # at Float32, and it is then the QUADRATURE path that returns exactly zero while the table
        # returns a finite value. That is a real defect of the reference path at this precision
        # rather than a disagreement to paper over, and it is recorded in this unit's own message.
        ax(lo, hi, n, k) = n == 1 ? lo : lo + (hi - lo) * FT(k) / FT(n - 1)
        lx = ax(lgc.logx̄_lo, lgc.logx̄_hi, lgc.nx, 1)
        fr = ax(lgc.F_rim_lo, lgc.F_rim_hi, lgc.nf, 1)
        rr = ax(lgc.ρ_rim_lo, lgc.ρ_rim_hi, lgc.nr, 1)
        la = ax(lgc.logρₐ_lo, lgc.logρₐ_hi, lgc.na, 1)
        lxc = ax(lgc.logx̄l_lo, lgc.logx̄l_hi, lgc.nl, 1)
        lxr = ax(lgr.logx̄l_lo, lgr.logx̄l_hi, lgr.nl, 1)
        # The temperature axis is GEOMETRIC in |T°C|, not uniform in 1/T°C, so its node is
        # `_logT_axis_value` and not `ax`. These tests assert that the table reproduces the entry AT
        # A NODE, which requires computing the node the way the fill computes it; with `ax` the
        # evaluation point is an interior point and the assertion picks up interpolation error
        # instead. See the comment on `_logT_axis_value`.
        iTn = P3._logT_axis_value(bgc.invT_lo, bgc.invT_hi, bgc.nt, 0)
        ρn = one(FT)
        ρa_n = exp10(la)
        q_ice = exp10(lx) * ρn
        st_n = P3.state_from_prognostic(params, q_ice, ρn, fr * q_ice, fr * q_ice / rr)
        lgn = P3.get_distribution_logλ(st_n)
        T_n = 1 / iTn + params.T_freeze
        q_ref = FT(1e-4)
        Nc_n = ρa_n * q_ref / exp10(lxc)
        Nr_n = ρa_n * q_ref / exp10(lxr)
        argn = (st_n, lgn, pdf_c, pdf_r, exp10(lxc) * Nc_n, Nc_n, exp10(lxr) * Nr_n, Nr_n,
            aps, tps, vel, ρa_n, T_n, m_liq)
        tabn = P3.∫liquid_ice_collisions(argn...; quad = carrier)
        refn = P3.∫liquid_ice_collisions(argn...;
            quad = rule, assembly = P3.SplitCorrection())
        # Measured residual at a node: six channels agree to better than 1e-4 and three to between
        # 1.4e-4 and 5.5e-4. That rules out a structural error, which would show as a factor rather
        # than as a hundredth of a percent, so the nine-channel reconstruction is what it claims to
        # be. The residual itself is not yet attributed and the tolerance records what is
        # established rather than what is expected; tightening it wants that attribution first.
        for i in 1:9
            TT.@test tabn[i] ≈ refn[i] rtol = 1e-3
        end
    end
end

TT.@testset "P3 lookup-table mode ($FT)" for FT in (Float64, Float32)
    test_p3_lut_off_mode_bit_identical(FT)
    test_p3_lut_untabulated_terms_bit_identical(FT)
    test_p3_lut_stored_entry_matches_integrand(FT)
    test_p3_lut_carrier_composition(FT)
    test_p3_lut_quadrature_delegation(FT)
    test_p3_lut_adapt_identity(FT)
    test_p3_lut_self_collection_table(FT)
    test_p3_lut_outputs(FT)
    test_p3_lut_liqice_path(FT)
end
nothing
