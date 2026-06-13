#####
##### 2M+P3 Rosenbrock-Euler substepping (`RosenbrockAverage`)
#####

"""
    MicroState2MP3{FT}

The eight prognostic 2M+P3 species as a `StaticArrays.FieldVector`: vector
algebra for substep updates and linear solves, named fields for readability.
Internal to the [`RosenbrockAverage`](@ref) implementation.
"""
struct MicroState2MP3{FT} <: SA.FieldVector{8, FT}
    q_lcl::FT
    n_lcl::FT
    q_rai::FT
    n_rai::FT
    q_ice::FT
    n_ice::FT
    q_rim::FT
    b_rim::FT
end
SA.similar_type(::Type{<:MicroState2MP3}, ::Type{FT}, ::SA.Size{(8,)}) where {FT} =
    MicroState2MP3{FT}

"""
    _instantaneous_2mp3_tendency(mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ)

The raw instantaneous 2M+P3 tendency projected onto the eight prognostic
species: the unlimited process rates of the `Microphysics2Moment` entry,
without timestep-dependent clipping.

This is the function the Rosenbrock step differentiates. It is the model
physics — an analytic time-averaged relaxation (Morrison and Milbrandt
2015, J. Atmos. Sci., 72, 287-311, Appendix C) carries no tendency clip;
the L-stable one-stage Rosenbrock update damps the stiff vapor-exchange
subsystem monotonically, so the overshoot an explicit-Euler supersaturation
cap suppresses cannot arise, while a 1/dt clip would add dt-independent
error and break convergence under refinement (Wan et al. 2020, J. Adv.
Model. Earth Syst., 12, e2019MS001982).
"""
@inline function _instantaneous_2mp3_tendency(mp, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
)
    full = bulk_microphysics_tendencies(Microphysics2Moment(), mp, tps,
        ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    )
    return (; full.dq_lcl_dt, full.dn_lcl_dt, full.dq_rai_dt, full.dn_rai_dt,
        full.dq_ice_dt, full.dn_ice_dt, full.dq_rim_dt, full.db_rim_dt)
end

"""
    Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)

Callable bundling the frozen per-substep context; applying it to the species
vector evaluates [`_instantaneous_2mp3_tendency`](@ref). A top-level struct
rather than a closure, so `ForwardDiff` differentiates a concretely-typed
callable.

The frozen `q_tot` is promoted to the state's element type at the call: a
zero-partial Dual is exact for a constant, and it keeps
`eltype(q_tot)`-keyed fallback returns in the tendency functions concretely
typed (a no-op in the primal pass). `logλ` stays plain so distribution-shape
slots never carry derivatives (the forward `gamma_inc` has no shape rule);
`T` and `ρ` stay plain because they feed working-type computations that must
remain floats.
"""
struct Instantaneous2MP3Tendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    T::F
    q_tot::F
    logλ::F
end
@inline function (g::Instantaneous2MP3Tendency)(x::SA.StaticVector{8})
    (q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim) = x
    tend = _instantaneous_2mp3_tendency(g.mp, g.tps,
        g.ρ, g.T, eltype(x)(g.q_tot),
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, g.logλ,
    )
    return MicroState2MP3(values(tend)...)
end

"""
    _per_process_2mp3(mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ)

Decompose the raw instantaneous 2M+P3 tendency into per-process contributions,
each a [`MicroState2MP3`](@ref) over the eight prognostic species, returned as
a `NamedTuple`. Replays the exact process calls of the warm-rain + P3 ice
[`bulk_microphysics_tendencies`](@ref) entry in the same order with the same
intermediate quantities, but routes each process's contribution into its own
species vector instead of accumulating all of them into one running sum. The
sum over the returned processes therefore equals the full raw tendency
[`_instantaneous_2mp3_tendency`](@ref) evaluates (asserted in the verbose
tests), so no physics is recomputed differently — the parts are merely exposed.

The processes are the warm-rain set (activation, cloud condensation/evaporation,
rain evaporation, autoconversion, cloud self-collection, accretion, rain
self-collection, rain breakup, cloud and rain number adjustment) and the P3 ice
set (liquid-ice collision/riming, ice aggregation, ice melting, F23 deposition
nucleation, Bigg immersion freezing of cloud drops, ice deposition/sublimation,
ice number adjustment, rain heterogeneous freezing). The ice processes are
gated identically to the entry, so each is exactly zero when its gate is shut.

Evaluated at the primal state only: it supplies the per-process right-hand
sides `f_p` for the linear post-solve attribution and is not differentiated.
"""
@inline function _per_process_2mp3(mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
) where {WR, ICE <: CMP.P3IceParams}
    FT = eltype(ρ)
    ϵₘ = UT.ϵ_numerics_2M_M(FT)
    ϵₙ = UT.ϵ_numerics_2M_N(FT)
    # Clamp negative inputs to zero, matching the entry.
    ρ = UT.clamp_to_nonneg(ρ)
    q_tot = UT.clamp_to_nonneg(q_tot)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_rai = UT.clamp_to_nonneg(q_rai)
    n_lcl = UT.clamp_to_nonneg(n_lcl)
    n_rai = UT.clamp_to_nonneg(n_rai)
    q_ice = UT.clamp_to_nonneg(q_ice)
    n_ice = UT.clamp_to_nonneg(n_ice)
    q_rim = UT.clamp_to_nonneg(q_rim)
    b_rim = UT.clamp_to_nonneg(b_rim)

    o = zero(FT)
    # builder for a per-process vector in the (q_lcl, n_lcl, q_rai, n_rai,
    # q_ice, n_ice, q_rim, b_rim) order shared with the entry's accumulators
    Z() = MicroState2MP3(o, o, o, o, o, o, o, o)

    # Volumetric quantities for P3 functions (entry convention).
    L_lcl = q_lcl * ρ
    L_rai = q_rai * ρ
    N_lcl = n_lcl * ρ
    N_rai = n_rai * ρ
    L_ice = q_ice * ρ
    N_ice = n_ice * ρ
    L_rim = q_rim * ρ
    B_rim = b_rim * ρ
    state = CMP3.state_from_prognostic(mp.ice.scheme, L_ice, N_ice, L_rim, B_rim)

    aps = mp.warm_rain.air_properties
    subdep = mp.warm_rain.subdep

    #####
    ##### Warm-rain processes (mirrors `warm_rain_tendencies_2m`)
    #####
    warm_rain = mp.warm_rain
    sb = warm_rain.seifert_beheng
    condevap = warm_rain.condevap
    N_lcl_wr = ρ * n_lcl
    N_rai_wr = ρ * n_rai

    # activation (cloud number only)
    dn_lcl_activation_dt = activation_source(
        warm_rain.activation_scheme, tps, ρ, T, q_tot, q_lcl, q_ice, n_lcl, o, o,
    )
    activation = MicroState2MP3(o, dn_lcl_activation_dt, o, o, o, o, o, o)

    # cloud condensation / evaporation (cloud mass only; number neglected)
    micro_mock = (; q_tot, q_lcl, q_icl = q_ice, q_rai, q_sno = zero(q_ice))
    thermo_mock = (; ρ, T)
    ∂ₜq_lcl_cond = CMNonEq.conv_q_vap_to_q_lcl(
        CMP.CloudLiquidFormation(condevap.τ_relax), nothing, tps, micro_mock, thermo_mock,
    )
    cloud_condevap = MicroState2MP3(∂ₜq_lcl_cond, o, o, o, o, o, o, o)

    # rain evaporation (rain mass + number)
    evap = CM2.rain_evaporation(sb, aps, tps, q_tot, q_lcl, q_ice, q_rai, zero(q_ice), ρ, N_rai_wr, T)
    rain_evap = MicroState2MP3(o, o, evap.∂ₜq_rai, evap.∂ₜρn_rai / ρ, o, o, o, o)

    # autoconversion (cloud → rain, mass + number)
    acnv = CM2.autoconversion(sb.acnv, sb.pdf_c, q_lcl, q_rai, ρ, N_lcl_wr)
    autoconv = MicroState2MP3(
        acnv.dq_lcl_dt, acnv.dN_lcl_dt / ρ, acnv.dq_rai_dt, acnv.dN_rai_dt / ρ, o, o, o, o,
    )

    # cloud self-collection (cloud number only)
    ∂ₜN_lcl_sc = CM2.cloud_liquid_self_collection(sb.acnv, sb.pdf_c, q_lcl, ρ, acnv.dN_lcl_dt)
    cloud_selfcol = MicroState2MP3(o, ∂ₜN_lcl_sc / ρ, o, o, o, o, o, o)

    # accretion (cloud → rain, mass; cloud number)
    accr = CM2.accretion(sb, q_lcl, q_rai, ρ, N_lcl_wr)
    accretion_wr = MicroState2MP3(accr.dq_lcl_dt, accr.dN_lcl_dt / ρ, accr.dq_rai_dt, o, o, o, o, o)

    # rain self-collection (rain number only)
    ∂ₜN_rai_sc = CM2.rain_self_collection(sb.pdf_r, sb.self, q_rai, ρ, N_rai_wr)
    rain_selfcol = MicroState2MP3(o, o, o, ∂ₜN_rai_sc / ρ, o, o, o, o)

    # rain breakup (rain number only)
    ∂ₜN_rai_br = CM2.rain_breakup(sb.pdf_r, sb.brek, q_rai, ρ, N_rai_wr, ∂ₜN_rai_sc)
    rain_breakup = MicroState2MP3(o, o, o, ∂ₜN_rai_br / ρ, o, o, o, o)

    # number adjustment for mass limits (cloud, then rain)
    numadj_lcl = (; sb.numadj.τ, x_min = sb.pdf_c.xc_min, x_max = sb.pdf_c.xc_max)
    ∂ₜn_lcl_numadj = CM2.number_tendency_from_mass_limits(numadj_lcl, q_lcl, n_lcl)
    cloud_numadj = MicroState2MP3(o, ∂ₜn_lcl_numadj, o, o, o, o, o, o)
    numadj_rai = (; sb.numadj.τ, x_min = sb.pdf_r.xr_min, x_max = sb.pdf_r.xr_max)
    ∂ₜn_rai_numadj = CM2.number_tendency_from_mass_limits(numadj_rai, q_rai, n_rai)
    rain_numadj = MicroState2MP3(o, o, o, ∂ₜn_rai_numadj, o, o, o, o)

    #####
    ##### P3 ice processes (mirrors the warm-rain + P3 ice entry)
    #####
    p3 = mp.ice.scheme
    vel = mp.ice.terminal_velocity
    pdf_c = mp.ice.cloud_pdf
    pdf_r = mp.ice.rain_pdf
    ice_nucleation = mp.ice.ice_nucleation
    inp_depletion_model = mp.ice.inp_depletion_model
    quad = mp.ice.quad

    # gated cluster: liquid-ice collision, ice aggregation, ice melting
    if q_ice > ϵₘ && n_ice > ϵₙ
        coll = CMP3.bulk_liquid_ice_collision_sources(
            state, logλ, pdf_c, pdf_r, L_lcl, N_lcl, L_rai, N_rai, aps, tps, vel, ρ, T;
            quad,
        )
        liquid_ice_collision = MicroState2MP3(
            coll.∂ₜq_c, coll.∂ₜN_c / ρ, coll.∂ₜq_r, coll.∂ₜN_r / ρ,
            coll.∂ₜL_ice / ρ, o, coll.∂ₜL_rim / ρ, coll.∂ₜB_rim / ρ,
        )

        S_ice_agg = CMP3.ice_self_collection(state, logλ, vel, ρ; quad)
        ice_aggregation = MicroState2MP3(o, o, o, o, o, -S_ice_agg.dNdt / ρ, o, o)

        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        melt = ifelse(T > T_freeze,
            CMP3.ice_melt(vel, aps, tps, T, ρ, state, logλ; quad),
            (; dNdt = zero(ρ), dLdt = zero(ρ)),
        )
        ∂ₜq_ice_melt = melt.dLdt / ρ
        ∂ₜn_ice_melt = melt.dNdt / ρ
        ∂ₜq_rim_melt = -∂ₜq_ice_melt * state.F_rim
        ∂ₜb_rim_melt = ifelse(state.ρ_rim > 0, -∂ₜq_ice_melt * state.F_rim / state.ρ_rim, zero(FT))
        ice_melting = MicroState2MP3(
            o, o, ∂ₜq_ice_melt, ∂ₜn_ice_melt, -∂ₜq_ice_melt, -∂ₜn_ice_melt,
            ∂ₜq_rim_melt, ∂ₜb_rim_melt,
        )
    else
        liquid_ice_collision = Z()
        ice_aggregation = Z()
        ice_melting = Z()
    end

    # F23 deposition nucleation (pristine ice, F_rim = 0)
    τ_act = inp_depletion_model.τ_act
    D_nuc = FT(10e-6)
    m_nuc = p3.ρ_i * CO.volume_sphere_D(D_nuc)
    n_active = CM_HetIce.n_active(inp_depletion_model, n_ice)
    f23_dep = CM_HetIce.f23_deposition_rate(
        ice_nucleation, tps, T, ρ, q_tot, q_lcl + q_rai, q_ice, n_active;
        m_nuc, τ_act, inpc_log_shift = zero(ρ),
    )
    f23_deposition = MicroState2MP3(o, o, o, o, f23_dep.∂ₜq_frz, f23_dep.∂ₜn_frz, o, o)

    # Bigg immersion freezing of cloud drops (fully-rimed embryo graupel)
    cld_bigg = CM_HetIce.liquid_freezing_rate(mp.ice.rain_freezing, pdf_c, tps, q_lcl, ρ, N_lcl, T)
    cld_cap = CM_HetIce.f23_immersion_limit_rate(
        ice_nucleation, T, ρ; τ = τ_act, inpc_log_shift = zero(ρ), n_active,
    )
    ∂ₜn_imm = min(cld_bigg.∂ₜn_frz, cld_cap.∂ₜn_frz)
    ∂ₜq_imm = ifelse(cld_bigg.∂ₜn_frz > 0, cld_bigg.∂ₜq_frz * ∂ₜn_imm / cld_bigg.∂ₜn_frz, zero(FT))
    bigg_immersion = MicroState2MP3(
        -∂ₜq_imm, -∂ₜn_imm, o, o, ∂ₜq_imm, ∂ₜn_imm, ∂ₜq_imm, ∂ₜq_imm / p3.ρ_i,
    )

    # ice deposition / sublimation (rim drains on the sublimation branch only)
    n_per_q_ice = ifelse(q_ice > ϵₘ, n_ice / q_ice, zero(n_ice))
    micro_mock_ice = (; q_tot, q_lcl, q_icl = q_ice, q_rai, q_sno = zero(q_ice))
    ∂ₜq_ice_dep = CMNonEq.conv_q_vap_to_q_icl(
        CMP.ConstantTimescale(subdep.τ_relax), nothing, tps, micro_mock_ice, thermo_mock,
    )
    ∂ₜq_ice_dep = ifelse(T > tps.T_freeze, min(∂ₜq_ice_dep, zero(T)), ∂ₜq_ice_dep)
    ∂ₜn_ice_dep = ifelse(∂ₜq_ice_dep < 0, n_per_q_ice * ∂ₜq_ice_dep, zero(∂ₜq_ice_dep))
    ∂ₜq_ice_sub = min(∂ₜq_ice_dep, 0)
    ∂ₜq_rim_sub = ∂ₜq_ice_sub * state.F_rim
    ∂ₜb_rim_sub = ifelse(state.ρ_rim > 0, ∂ₜq_ice_sub * state.F_rim / state.ρ_rim, zero(FT))
    ice_depsub = MicroState2MP3(o, o, o, o, ∂ₜq_ice_dep, ∂ₜn_ice_dep, ∂ₜq_rim_sub, ∂ₜb_rim_sub)

    # ice number adjustment for mass limits
    numadj = (; τ = FT(100), x_min = FT(1e-12), x_max = FT(1e-5))
    ∂ₜn_ice_numadj = CM2.number_tendency_from_mass_limits(numadj, q_ice, n_ice)
    ice_numadj = MicroState2MP3(o, o, o, o, o, ∂ₜn_ice_numadj, o, o)

    # rain heterogeneous freezing (Bigg; frozen rain fully rimed)
    rain_frz = CM_HetIce.liquid_freezing_rate(mp.ice.rain_freezing, pdf_r, tps, q_rai, ρ, N_rai, T)
    rain_freezing = MicroState2MP3(
        o, o, -rain_frz.∂ₜq_frz, -rain_frz.∂ₜn_frz,
        rain_frz.∂ₜq_frz, rain_frz.∂ₜn_frz, rain_frz.∂ₜq_frz, rain_frz.∂ₜq_frz / p3.ρ_i,
    )

    return (;
        activation, cloud_condevap, rain_evap, autoconv, cloud_selfcol,
        accretion = accretion_wr, rain_selfcol, rain_breakup, cloud_numadj, rain_numadj,
        liquid_ice_collision, ice_aggregation, ice_melting, f23_deposition,
        bigg_immersion, ice_depsub, ice_numadj, rain_freezing,
    )
end

"""
    Verbose2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)

Per-process companion to [`Instantaneous2MP3Tendency`](@ref): applying it to
the species vector returns a `NamedTuple` of per-process tendency contributions
(each a [`MicroState2MP3`](@ref)) via [`_per_process_2mp3`](@ref), instead of
only their sum. The frozen context and the primal physics are identical to
[`Instantaneous2MP3Tendency`](@ref); this functor is evaluated at the primal
state only (it supplies the right-hand sides `f_p` for the linear post-solve
attribution and is not itself differentiated).
"""
struct Verbose2MP3Tendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    T::F
    q_tot::F
    logλ::F
end
@inline function (g::Verbose2MP3Tendency)(x::SA.StaticVector{8, FT}) where {FT}
    (q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim) = x
    return _per_process_2mp3(g.mp, g.tps,
        g.ρ, g.T, FT(g.q_tot),
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, g.logλ,
    )
end

"""
    _rosenbrock_channel_mask(x)

Diagonal of the channel projection matrix `P` used by
[`_rosenbrock_update`](@ref): 1 for healthy channels, 0 for near-empty ones
(condensed mass below `1e-10`, per channel: liquid, rain, ice+rime).

This is a scheme partition, not a preconditioner: with `P J P` in place of
`J`, a masked channel's row of `(I/h - P J P)` is the identity row, so the
solve returns exactly the forward-Euler update for that channel while
healthy channels stay implicit (an IMEX-style splitting at channel
granularity). The gate is load-bearing: near-empty channels produce finite
but enormous Jacobian entries whose linearized steady state fabricates
phantom number concentrations that substep refinement cannot remove.
"""
@inline function _rosenbrock_channel_mask(x::MicroState2MP3{FT}) where {FT}
    ϵ_empty = FT(1e-10)
    liq = ifelse(x.q_lcl < ϵ_empty, zero(FT), one(FT))
    rai = ifelse(x.q_rai < ϵ_empty, zero(FT), one(FT))
    ice = ifelse(x.q_ice < ϵ_empty, zero(FT), one(FT))
    return MicroState2MP3(liq, liq, rai, rai, ice, ice, ice, ice)
end

"""
    _euler_update(x, f, h)

Forward-Euler substep, floored at zero.
"""
@inline _euler_update(x, f, h) = max.(x .+ h .* f, 0)

"""
    _rosenbrock_system(x, f, J, z, h)

Build the equilibrated linear system of one linearized-implicit
(Rosenbrock-Euler) substep at state `x` with raw tendency `f`, Jacobian `J`,
channel mask `z`, and substep `h`. Returns `(S, S⁻¹, A)`, the equilibration
matrix `S = Diagonal(|x| + h |f| + ϵ)`, its inverse, and the equilibrated
system matrix `A = I/h - S⁻¹ (P J P) S`, where `P = Diagonal(z)` is the
channel projection built from the per-scheme channel mask `z` (e.g.
[`_rosenbrock_channel_mask`](@ref) for 2M+P3).

Equilibration makes the similarity transform `S⁻¹ A S` O(1)-conditioned: the
raw rows span ~9 orders of magnitude (number vs mass species), and an
unscaled Float32 factorization bleeds roundoff from the large rows into empty
species as phantom mass. It is exact in exact arithmetic and keeps roundoff
relative to each species' own scale.

Splitting the system build from the solve lets multiple right-hand sides
reuse one factorization context: the full-step update and the per-process
attribution (see [`_rosenbrock_solve`](@ref)) solve against the SAME `S`,
`S⁻¹`, `A`, so per-process increments sum exactly to the full increment.
"""
@inline function _rosenbrock_system(x::SA.StaticVector{N, FT}, f, J, z, h) where {N, FT}
    Iₙ = one(SA.SMatrix{N, N, FT})
    s = abs.(x) .+ h .* abs.(f) .+ eps(FT)
    # dense diagonal matrices: an `SDiagonal` wrapper here defeats the
    # optimizer's static-array stack allocation (heap spills per substep)
    P = Iₙ .* z'
    S = Iₙ .* s'
    S⁻¹ = Iₙ .* inv.(s)'
    A = Iₙ / h - S⁻¹ * (P * J * P) * S
    return S, S⁻¹, A
end

"""
    _rosenbrock_solve(S, S⁻¹, A, v)

Solve the equilibrated Rosenbrock system from [`_rosenbrock_system`](@ref) for
the unclamped increment of right-hand side `v`: `Δ = S (A \\ (S⁻¹ v))`, the
equilibrated form of `(I/h - P J P)⁻¹ v`.

This is LINEAR in `v`: with `S`, `S⁻¹`, `A` fixed across calls, `Σ_p
solve(v_p) = solve(Σ_p v_p)` exactly. Per-process attribution exploits this —
each process tendency is pushed through the same linear correction and the
increments sum to the full-step increment with no residual beyond the
positivity clamp (handled separately by the caller).
"""
@inline _rosenbrock_solve(S, S⁻¹, A, v) = S * (A \ (S⁻¹ * v))

"""
    _rosenbrock_update(x, f, J, z, h)

One linearized-implicit (Rosenbrock-Euler) substep: build the equilibrated
system ([`_rosenbrock_system`](@ref)) for `(I/h - P J P) Δx = f`, solve it
([`_rosenbrock_solve`](@ref)), and return `max.(x + Δx, 0)`.
"""
@inline function _rosenbrock_update(x::SA.StaticVector{N, FT}, f, J, z, h) where {N, FT}
    S, S⁻¹, A = _rosenbrock_system(x, f, J, z, h)
    Δx = _rosenbrock_solve(S, S⁻¹, A, f)
    return max.(x .+ Δx, 0)
end

"""
    bulk_microphysics_tendencies(::RosenbrockAverage, ::Microphysics2Moment,
        mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub = 1)

Compute average 2M+P3 microphysics tendencies over `Δt` using `nsub`
linearized-implicit (Rosenbrock-Euler) substeps of the raw instantaneous
tendency.

# Algorithm

For each substep of `h = Δt / nsub`:

1. Evaluate the raw tendency `f` ([`_instantaneous_2mp3_tendency`](@ref)) and
   its exact 8×8 Jacobian `J` via `ForwardDiff` at the current state.
2. Advance with [`_rosenbrock_update`](@ref): solve `(I/h - P J P) Δx = f`
   in equilibrated variables, where the projection `P` routes near-empty
   channels to forward Euler.
3. Update the local temperature from the latent heating of the realized
   increments.

A non-finite state cannot be differentiated and a non-finite Jacobian (an
exotic state escaping the channel mask) cannot be linearized; both fall back
to a forward-Euler substep of the raw tendency. `logλ` and `q_tot` are held
fixed across substeps, matching the explicit-substepping semantics. The
discrete safeguards are h-free conditioning and projection (channel mask,
equilibration, positivity clamp), not model terms: the differentiated
tendency is the unmodified relaxation physics.

Returns the net change in the species over `Δt` divided by `Δt`, in the same
fields as the `Instantaneous` entry (without the activation diagnostic).
"""
@inline function bulk_microphysics_tendencies(::RosenbrockAverage, cm::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    Δt, nsub = 1,
) where {WR, ICE <: CMP.P3IceParams}
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    cp_d = TDI.TD.Parameters.cp_d(tps)

    x = MicroState2MP3{FT}(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    x₀ = x
    Tsub = T
    for _ in 1:nsub_eff
        g = Instantaneous2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ)
        f = g(x)
        x_prev = x
        x = if all(isfinite, x)
            J = FD.jacobian(g, x)
            z = _rosenbrock_channel_mask(x)
            all(isfinite, J) ? _rosenbrock_update(x, f, J, z, h) : _euler_update(x, f, h)
        else
            _euler_update(x, f, h)
        end
        Δ = x - x_prev
        T_safe = max(150, Tsub)
        Tsub += (TDI.Lᵥ(tps, T_safe) * (Δ.q_lcl + Δ.q_rai) + TDI.Lₛ(tps, T_safe) * Δ.q_ice) / cp_d
    end

    rates = (x - x₀) / Δt
    return NamedTuple{(
        :dq_lcl_dt, :dn_lcl_dt, :dq_rai_dt, :dn_rai_dt,
        :dq_ice_dt, :dn_ice_dt, :dq_rim_dt, :db_rim_dt,
    )}(
        Tuple(rates),
    )
end

#####
##### 1M Rosenbrock-Euler substepping (`RosenbrockAverage`)
#####

"""
    MicroState1M{FT}

The four prognostic 1M species as a `StaticArrays.FieldVector`: vector algebra
for substep updates and linear solves, named fields for readability. Internal
to the [`RosenbrockAverage`](@ref) implementation, mirroring
[`MicroState2MP3`](@ref) at dimension 4 (1M carries no number species).
"""
struct MicroState1M{FT} <: SA.FieldVector{4, FT}
    q_lcl::FT
    q_icl::FT
    q_rai::FT
    q_sno::FT
end
SA.similar_type(::Type{<:MicroState1M}, ::Type{FT}, ::SA.Size{(4,)}) where {FT} =
    MicroState1M{FT}

"""
    _instantaneous_1m_tendency(mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno)

The raw instantaneous 1M tendency projected onto the four prognostic species:
the unlimited process rates of the `Microphysics1Moment` `Instantaneous` entry
(`_microphysics_source_terms` aggregated by `_aggregate_tendencies`), without
timestep-dependent clipping.

This is the function the 1M Rosenbrock step differentiates — the same raw RHS
the hand-built `LinearizedAverage` donor-linearizes (the hand operator is the
system matrix of a donor-modified problem, not `∂f/∂q`). The shipped 1M
`LinearizedAverage` applies no supersaturation cap or coupled-sink limiter, so
differentiating the raw tendency is the faithful counterpart of that scheme.
"""
@inline function _instantaneous_1m_tendency(mp, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    src = _microphysics_source_terms(Microphysics1Moment(), mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )
    return _aggregate_tendencies(src)
end

"""
    Raw1MTendency(mp, tps, ρ, T, q_tot)

Callable bundling the frozen per-substep context; applying it to the species
vector evaluates [`_instantaneous_1m_tendency`](@ref). A top-level struct
rather than a closure, so `ForwardDiff` differentiates a concretely-typed
callable, mirroring [`Instantaneous2MP3Tendency`](@ref).

The frozen `q_tot` is promoted to the state's element type at the call so a
zero-partial Dual is exact for the constant and the promotion-keyed fallback
returns in the 1M kernels stay concretely typed (a no-op in the primal pass);
`T` and `ρ` stay plain because they feed working-type computations that must
remain floats. The differentiability of the 1M kernels w.r.t. the mass
channels (with `ρ`/`T`/params held float) is provided by their AD-readiness:
the collected/collecting masses are unconstrained relative to the parameter
type, so a Dual working type flows through without a parameter-type converter.
"""
struct Raw1MTendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    T::F
    q_tot::F
end
@inline function (g::Raw1MTendency)(x::SA.StaticVector{4})
    (q_lcl, q_icl, q_rai, q_sno) = x
    tend = _instantaneous_1m_tendency(g.mp, g.tps,
        g.ρ, g.T, eltype(x)(g.q_tot),
        q_lcl, q_icl, q_rai, q_sno,
    )
    return MicroState1M(tend.dq_lcl_dt, tend.dq_icl_dt, tend.dq_rai_dt, tend.dq_sno_dt)
end

"""
    _per_process_1m(src)

Project each individual 1M source term onto the four prognostic species as a
[`MicroState1M`](@ref) `(q_lcl, q_icl, q_rai, q_sno)`, returning a `NamedTuple`
of these per-process contribution vectors. The signs mirror exactly the single
aggregation point [`_aggregate_tendencies`](@ref), so by construction the sum
over the returned processes equals the aggregated raw tendency the
[`Raw1MTendency`](@ref) functor evaluates — no physics is recomputed, the same
`src` terms are merely re-grouped by process instead of by species.

The processes are the eighteen terms of `_microphysics_source_terms`: vapor ↔
cloud condensate phase change, autoconversion, the accretion family (with its
cold/warm collision arms and thermal-melt by-products pre-routed by
temperature in `src`, so each appears here with a fixed sign), precipitation ↔
vapor phase change, and the two melting transfers. Used by the verbose
post-solve attribution to supply the per-process right-hand sides `f_p`.
"""
@inline function _per_process_1m(src)
    FT = typeof(src.S_phase_change_vap_lcl)
    o = zero(FT)
    # each `MicroState1M(q_lcl, q_icl, q_rai, q_sno)` holds one process's signed
    # contribution; the signs are copied verbatim from `_aggregate_tendencies`
    return (;
        phase_change_vap_lcl = MicroState1M(src.S_phase_change_vap_lcl, o, o, o),
        phase_change_vap_icl = MicroState1M(o, src.S_phase_change_vap_icl, o, o),
        acnv_lcl_rai = MicroState1M(-src.S_acnv_lcl_rai, o, src.S_acnv_lcl_rai, o),
        acnv_icl_sno = MicroState1M(o, -src.S_acnv_icl_sno, o, src.S_acnv_icl_sno),
        accr_lcl_rai = MicroState1M(-src.S_accr_lcl_rai, o, src.S_accr_lcl_rai, o),
        accr_lcl_sno_cold = MicroState1M(-src.S_accr_lcl_sno_cold, o, o, src.S_accr_lcl_sno_cold),
        accr_lcl_sno_warm = MicroState1M(-src.S_accr_lcl_sno_warm, o, src.S_accr_lcl_sno_warm, o),
        accr_melt_lcl_sno = MicroState1M(o, o, src.S_accr_melt_lcl_sno, -src.S_accr_melt_lcl_sno),
        accr_icl_rai = MicroState1M(o, -src.S_accr_icl_rai, o, src.S_accr_icl_rai),
        accr_freeze_icl_rai = MicroState1M(o, o, -src.S_accr_freeze_icl_rai, src.S_accr_freeze_icl_rai),
        accr_icl_sno = MicroState1M(o, -src.S_accr_icl_sno, o, src.S_accr_icl_sno),
        accr_rai_sno_cold = MicroState1M(o, o, -src.S_accr_rai_sno_cold, src.S_accr_rai_sno_cold),
        accr_rai_sno_warm = MicroState1M(o, o, src.S_accr_rai_sno_warm, -src.S_accr_rai_sno_warm),
        accr_melt_rai_sno = MicroState1M(o, o, src.S_accr_melt_rai_sno, -src.S_accr_melt_rai_sno),
        phase_change_vap_rai = MicroState1M(o, o, src.S_phase_change_vap_rai, o),
        phase_change_vap_sno = MicroState1M(o, o, o, src.S_phase_change_vap_sno),
        melt_icl_lcl = MicroState1M(src.S_melt_icl_lcl, -src.S_melt_icl_lcl, o, o),
        melt_sno_rai = MicroState1M(o, o, src.S_melt_sno_rai, -src.S_melt_sno_rai),
    )
end

"""
    Verbose1MTendency(mp, tps, ρ, T, q_tot)

Per-process companion to [`Raw1MTendency`](@ref): applying it to the species
vector returns a `NamedTuple` of per-process tendency contributions (each a
[`MicroState1M`](@ref)) via [`_per_process_1m`](@ref), instead of only their
sum. The frozen context and the primal physics are identical to
[`Raw1MTendency`](@ref); this functor is evaluated at the primal state only (it
supplies the right-hand sides `f_p` for the linear post-solve attribution and
is not itself differentiated).
"""
struct Verbose1MTendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    T::F
    q_tot::F
end
@inline function (g::Verbose1MTendency)(x::SA.StaticVector{4, FT}) where {FT}
    (q_lcl, q_icl, q_rai, q_sno) = x
    src = _microphysics_source_terms(Microphysics1Moment(), g.mp, g.tps,
        g.ρ, g.T, FT(g.q_tot),
        q_lcl, q_icl, q_rai, q_sno,
    )
    return _per_process_1m(src)
end

"""
    _rosenbrock_channel_mask(x::MicroState1M)

Diagonal of the channel projection `P` for the 1M state: 1 for healthy
channels, 0 for near-empty ones (mass below `1e-10`, per channel: liquid, ice,
rain, snow). With no number species the mask is the four mass channels
directly. See the [`MicroState2MP3`](@ref) method for the role of `P` in
[`_rosenbrock_update`](@ref).
"""
@inline function _rosenbrock_channel_mask(x::MicroState1M{FT}) where {FT}
    ϵ_empty = FT(1e-10)
    lcl = ifelse(x.q_lcl < ϵ_empty, zero(FT), one(FT))
    icl = ifelse(x.q_icl < ϵ_empty, zero(FT), one(FT))
    rai = ifelse(x.q_rai < ϵ_empty, zero(FT), one(FT))
    sno = ifelse(x.q_sno < ϵ_empty, zero(FT), one(FT))
    return MicroState1M(lcl, icl, rai, sno)
end

"""
    bulk_microphysics_tendencies(::RosenbrockAverage, ::Microphysics1Moment,
        mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1)

Compute average 1M microphysics tendencies over `Δt` using `nsub`
linearized-implicit (Rosenbrock-Euler) substeps of the raw instantaneous
tendency. A separate option from the hand-built `LinearizedAverage`: this
linearizes with the exact `ForwardDiff` Jacobian of the raw tendency at each
substep instead of a donor-cell-modified system matrix.

# Algorithm

For each substep of `h = Δt / nsub`:

1. Evaluate the raw tendency `f` ([`_instantaneous_1m_tendency`](@ref)) and its
   exact 4×4 Jacobian `J` via `ForwardDiff` at the current state.
2. Advance with [`_rosenbrock_update`](@ref): solve `(I/h - P J P) Δx = f` in
   equilibrated variables, where the projection `P` from
   [`_rosenbrock_channel_mask`](@ref) routes near-empty channels to forward
   Euler.
3. Update the local temperature from the latent heating of the realized
   increments, matching the `LinearizedAverage` 1M convention (constant latent
   heats, liquid+rain on `L_v`, ice+snow on `L_s`).

A non-finite state or Jacobian falls back to a forward-Euler substep of the raw
tendency; `q_tot` is held fixed across substeps. The discrete safeguards are
h-free conditioning and projection (channel mask, equilibration, positivity
clamp), not model terms.

Returns the net change in the species over `Δt` divided by `Δt`, in the same
fields as the `Instantaneous` entry.
"""
@inline function bulk_microphysics_tendencies(::RosenbrockAverage, cm::Microphysics1Moment,
    mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1,
)
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)

    x = MicroState1M{FT}(q_lcl, q_icl, q_rai, q_sno)
    x₀ = x
    Tsub = T
    for _ in 1:nsub_eff
        g = Raw1MTendency(mp, tps, ρ, Tsub, q_tot)
        f = g(x)
        x_prev = x
        x = if all(isfinite, x)
            J = FD.jacobian(g, x)
            z = _rosenbrock_channel_mask(x)
            all(isfinite, J) ? _rosenbrock_update(x, f, J, z, h) : _euler_update(x, f, h)
        else
            _euler_update(x, f, h)
        end
        Δ = x - x_prev
        Tsub += Lv_over_cp * (Δ.q_lcl + Δ.q_rai) + Ls_over_cp * (Δ.q_icl + Δ.q_sno)
    end

    rates = (x - x₀) / Δt
    return (;
        dq_lcl_dt = rates.q_lcl, dq_icl_dt = rates.q_icl,
        dq_rai_dt = rates.q_rai, dq_sno_dt = rates.q_sno,
    )
end

#####
##### Verbose post-solve per-process attribution (`RosenbrockAverageVerbose`)
#####

"""
    _rosenbrock_substep_verbose(g, g_verbose, x, h)

One Rosenbrock-Euler substep with post-solve per-process attribution. Returns
`(x_new, Δx_processes, Δx_clamp)`:

- `x_new` — the realized next state, identical to the non-verbose
  [`_rosenbrock_update`](@ref) / [`_euler_update`](@ref) at the same inputs.
- `Δx_processes` — a `NamedTuple` of per-process realized increments `Δx_p`
  (each a state vector), keyed by the verbose functor `g_verbose`'s processes.
- `Δx_clamp` — the positivity clamp correction `(x_new − x) − Σ_p Δx_p` as a
  state vector.

The substep update is LINEAR in the raw tendency `f` (the equilibrated solve
[`_rosenbrock_solve`](@ref), or forward Euler on the fallback path), so pushing
each process's instantaneous contribution `f_p` through the SAME operator gives
`Σ_p Δx_p` equal to the unclamped increment `Δx` exactly. The realized state
update `max.(x + Δx, 0)` is NOT linear, so the clamp's effect is not
attributable to any process: it is returned separately in `Δx_clamp`, never
folded into `Δx_processes`. The fallbacks (non-finite state or Jacobian) take
forward Euler, where the same linearity holds with `Δx_p = h f_p`.

`g` is the summed tendency functor (differentiated for the Jacobian); for the
solve to be the same operator the per-process functor `g_verbose` must compute
the same physics at the same state — its per-process parts summing to `g(x)`.

Not `@inline`d on purpose: kept as its own specialization so the verbose
functor call is analyzed with the concrete `MicroState{N,FT}` state type, which
sidesteps a JET tuple-broadcast false report in the shared P3 state
constructor that surfaces only when the functor is inlined into a
construct-and-call context. The verbose path is diagnostic, not the hot loop,
so the call boundary costs nothing that matters.
"""
function _rosenbrock_substep_verbose(g, g_verbose, x::SA.StaticVector{N, FT}, h) where {N, FT}
    f = g(x)
    fp = g_verbose(x)
    if all(isfinite, x)
        J = FD.jacobian(g, x)
        z = _rosenbrock_channel_mask(x)
        if all(isfinite, J)
            S, S⁻¹, A = _rosenbrock_system(x, f, J, z, h)
            Δx = _rosenbrock_solve(S, S⁻¹, A, f)
            Δxp = map(fp_i -> _rosenbrock_solve(S, S⁻¹, A, fp_i), fp)
            x_new = max.(x .+ Δx, 0)
            return x_new, Δxp, (x_new - x) - Δx
        end
    end
    Δx = h .* f
    Δxp = map(fp_i -> h .* fp_i, fp)
    x_new = max.(x .+ Δx, 0)
    return x_new, Δxp, (x_new - x) - Δx
end

"""
    _per_process_zero_seed(g_verbose, x)

Zero-valued per-process accumulator matching the `NamedTuple` shape the verbose
functor `g_verbose` returns at state `x`: each process slot set to `zero(x)`.
Used to seed the per-substep accumulation in the verbose averaged entries.

Not `@inline`d, for the same reason as [`_rosenbrock_substep_verbose`](@ref):
the single verbose-functor call stays in its own specialization with the
concrete `MicroState{N,FT}` state type, avoiding a JET tuple-broadcast false
report in the shared P3 state constructor under inlining.
"""
function _per_process_zero_seed(g_verbose, x::SA.StaticVector)
    return map(_ -> zero(x), g_verbose(x))
end

"""
    bulk_microphysics_tendencies(::RosenbrockAverageVerbose, ::Microphysics2Moment,
        mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub = 1)

Diagnostic 2M+P3 Rosenbrock averaged tendency with post-solve per-process
attribution. Runs the identical substep loop as the
[`RosenbrockAverage`](@ref) entry (same state trajectory, same T update, same
net result) and, in addition, accumulates the per-process realized increments
([`_per_process_2mp3`](@ref)) and the positivity clamp correction through
[`_rosenbrock_substep_verbose`](@ref).

Returns a `NamedTuple` with:

- `dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt, dq_ice_dt, dn_ice_dt, dq_rim_dt,
  db_rim_dt` — the net averaged tendencies, identical to
  [`RosenbrockAverage`](@ref).
- `processes` — a `NamedTuple` of per-process realized averaged tendencies
  (accumulated `Σ_substeps Δx_p / Δt`), each a [`MicroState2MP3`](@ref).
- `clamp_correction` — the non-attributable positivity-clamp tendency
  (accumulated `Σ_substeps Δx_clamp / Δt`), a [`MicroState2MP3`](@ref).

By construction `Σ_p processes_p + clamp_correction` equals the net averaged
state change `(x − x₀) / Δt` to the roundoff of the per-substep linear solve.
A diagnostic path: it may allocate (the per-process `map` builds intermediate
state vectors); the non-verbose hot loop is untouched.
"""
@inline function bulk_microphysics_tendencies(::RosenbrockAverageVerbose, cm::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    Δt, nsub = 1,
) where {WR, ICE <: CMP.P3IceParams}
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    cp_d = TDI.TD.Parameters.cp_d(tps)

    x = MicroState2MP3{FT}(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    x₀ = x
    Tsub = T
    Δxp_sum = _per_process_zero_seed(Verbose2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ), x)
    Δx_clamp_sum = zero(x)
    for _ in 1:nsub_eff
        g = Instantaneous2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ)
        gv = Verbose2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ)
        x_prev = x
        x, Δxp, Δx_clamp = _rosenbrock_substep_verbose(g, gv, x, h)
        Δxp_sum = map(+, Δxp_sum, Δxp)
        Δx_clamp_sum += Δx_clamp
        Δ = x - x_prev
        T_safe = max(150, Tsub)
        Tsub += (TDI.Lᵥ(tps, T_safe) * (Δ.q_lcl + Δ.q_rai) + TDI.Lₛ(tps, T_safe) * Δ.q_ice) / cp_d
    end

    rates = (x - x₀) / Δt
    net = NamedTuple{(
        :dq_lcl_dt, :dn_lcl_dt, :dq_rai_dt, :dn_rai_dt,
        :dq_ice_dt, :dn_ice_dt, :dq_rim_dt, :db_rim_dt,
    )}(
        Tuple(rates),
    )
    processes = map(Δxp_i -> Δxp_i / Δt, Δxp_sum)
    clamp_correction = Δx_clamp_sum / Δt
    return merge(net, (; processes, clamp_correction))
end

"""
    bulk_microphysics_tendencies(::RosenbrockAverageVerbose, ::Microphysics1Moment,
        mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1)

Diagnostic 1M Rosenbrock averaged tendency with post-solve per-process
attribution. Runs the identical substep loop as the
[`RosenbrockAverage`](@ref) 1M entry and, in addition, accumulates the
per-process realized increments ([`_per_process_1m`](@ref)) and the positivity
clamp correction through [`_rosenbrock_substep_verbose`](@ref).

Returns a `NamedTuple` with:

- `dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt` — the net averaged tendencies,
  identical to [`RosenbrockAverage`](@ref).
- `processes` — a `NamedTuple` of per-process realized averaged tendencies
  (accumulated `Σ_substeps Δx_p / Δt`), each a [`MicroState1M`](@ref).
- `clamp_correction` — the non-attributable positivity-clamp tendency
  (accumulated `Σ_substeps Δx_clamp / Δt`), a [`MicroState1M`](@ref).

By construction `Σ_p processes_p + clamp_correction` equals the net averaged
state change `(x − x₀) / Δt` to the roundoff of the per-substep linear solve.
A diagnostic path: it may allocate; the non-verbose hot loop is untouched.
"""
@inline function bulk_microphysics_tendencies(::RosenbrockAverageVerbose, cm::Microphysics1Moment,
    mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1,
)
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)

    x = MicroState1M{FT}(q_lcl, q_icl, q_rai, q_sno)
    x₀ = x
    Tsub = T
    Δxp_sum = _per_process_zero_seed(Verbose1MTendency(mp, tps, ρ, Tsub, q_tot), x)
    Δx_clamp_sum = zero(x)
    for _ in 1:nsub_eff
        g = Raw1MTendency(mp, tps, ρ, Tsub, q_tot)
        gv = Verbose1MTendency(mp, tps, ρ, Tsub, q_tot)
        x_prev = x
        x, Δxp, Δx_clamp = _rosenbrock_substep_verbose(g, gv, x, h)
        Δxp_sum = map(+, Δxp_sum, Δxp)
        Δx_clamp_sum += Δx_clamp
        Δ = x - x_prev
        Tsub += Lv_over_cp * (Δ.q_lcl + Δ.q_rai) + Ls_over_cp * (Δ.q_icl + Δ.q_sno)
    end

    rates = (x - x₀) / Δt
    net = (;
        dq_lcl_dt = rates.q_lcl, dq_icl_dt = rates.q_icl,
        dq_rai_dt = rates.q_rai, dq_sno_dt = rates.q_sno,
    )
    processes = map(Δxp_i -> Δxp_i / Δt, Δxp_sum)
    clamp_correction = Δx_clamp_sum / Δt
    return merge(net, (; processes, clamp_correction))
end

#####
##### Implicit-temperature Rosenbrock-Euler substepping (`RosenbrockAverageImplicitT`)
#####

"""
    MicroState2MP3T{FT}

The eight prognostic 2M+P3 species plus the local temperature `T` as a
`StaticArrays.FieldVector`, the augmented state of the
[`RosenbrockAverageImplicitT`](@ref) 2M+P3 entry. Mirrors
[`MicroState2MP3`](@ref) with `T` appended as the final component, so the
species ordering of the leading eight-component block is unchanged and the
species-only kernels see the same layout.

Carrying `T` in the state makes it a differentiation variable: the
species+temperature Jacobian then holds the latent-heating row and the
Clausius-Clapeyron column `∂f/∂T` that the species-only solve splits off into
an explicit between-substep update.
"""
struct MicroState2MP3T{FT} <: SA.FieldVector{9, FT}
    q_lcl::FT
    n_lcl::FT
    q_rai::FT
    n_rai::FT
    q_ice::FT
    n_ice::FT
    q_rim::FT
    b_rim::FT
    T::FT
end
SA.similar_type(::Type{<:MicroState2MP3T}, ::Type{FT}, ::SA.Size{(9,)}) where {FT} =
    MicroState2MP3T{FT}

"""
    _latent_heating_rate_2mp3(tps, T, dq_lcl_dt, dq_rai_dt, dq_ice_dt)

The temperature tendency `dT/dt` from the latent heat of the 2M+P3 phase-change
rates: `(Lᵥ·(dq_lcl_dt + dq_rai_dt) + Lₛ·dq_ice_dt) / cp_d`, the rate form of
the explicit between-substep update of [`RosenbrockAverage`](@ref) (which adds
`(Lᵥ·(Δq_lcl + Δq_rai) + Lₛ·Δq_ice) / cp_d` from the realized increment `Δ`).
Liquid and rain mass condensing/evaporating release/absorb the latent heat of
vaporization `Lᵥ`; ice mass depositing/sublimating (net of melting/freezing,
which the P3 ice mass balance already nets) releases/absorbs the latent heat of
sublimation `Lₛ`. The latent heats are evaluated at `max(150, T)`, matching the
floor the explicit update applies, so an explicit forward-Euler integration of
the augmented tendency reproduces the [`RosenbrockAverage`](@ref) trajectory.

Written as a function of `T` (through `Lᵥ`, `Lₛ`) and of the rates (which depend
on `T` through the saturation specific humidities) so that, with `T` carried in
the augmented state, `ForwardDiff` differentiates the full
temperature → saturation → exchange-flux → latent-heating loop into the
Jacobian.
"""
@inline function _latent_heating_rate_2mp3(tps, T, dq_lcl_dt, dq_rai_dt, dq_ice_dt)
    cp_d = TDI.TD.Parameters.cp_d(tps)
    T_safe = max(oftype(T, 150), T)
    return (TDI.Lᵥ(tps, T_safe) * (dq_lcl_dt + dq_rai_dt) + TDI.Lₛ(tps, T_safe) * dq_ice_dt) / cp_d
end

"""
    Instantaneous2MP3TTendency(mp, tps, ρ, q_tot, logλ)

Callable bundling the frozen per-substep context for the augmented 2M+P3
state; applying it to a [`MicroState2MP3T`](@ref) evaluates the eight species
tendencies ([`_instantaneous_2mp3_tendency`](@ref)) at the state's own `T` and
appends the latent-heating tendency `dT/dt`
([`_latent_heating_rate_2mp3`](@ref)) as the ninth component. A top-level struct
rather than a closure, so `ForwardDiff` differentiates a concretely-typed
callable, mirroring [`Instantaneous2MP3Tendency`](@ref).

Unlike [`Instantaneous2MP3Tendency`](@ref), `T` is NOT frozen context: it is
taken from the state, so the partial w.r.t. the ninth component flows through
the saturation specific humidities (the Clausius-Clapeyron feedback) and the
temperature-dependent latent heats. `ρ`, `q_tot`, and `logλ` are frozen across
substeps exactly as in the species-only entry; `q_tot` is promoted to the
state's element type at the call so the constant carries a zero partial and the
promotion-keyed fallback returns stay concretely typed.
"""
struct Instantaneous2MP3TTendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    q_tot::F
    logλ::F
end
@inline function (g::Instantaneous2MP3TTendency)(x::SA.StaticVector{9})
    (q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, T) = x
    tend = _instantaneous_2mp3_tendency(g.mp, g.tps,
        g.ρ, T, eltype(x)(g.q_tot),
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, g.logλ,
    )
    dT_dt = _latent_heating_rate_2mp3(g.tps, T, tend.dq_lcl_dt, tend.dq_rai_dt, tend.dq_ice_dt)
    return MicroState2MP3T(
        tend.dq_lcl_dt, tend.dn_lcl_dt, tend.dq_rai_dt, tend.dn_rai_dt,
        tend.dq_ice_dt, tend.dn_ice_dt, tend.dq_rim_dt, tend.db_rim_dt, dT_dt,
    )
end

"""
    _rosenbrock_channel_mask(x::MicroState2MP3T)

Diagonal of the channel projection `P` for the augmented 2M+P3 state: the eight
species channels gated exactly as the species-only
[`MicroState2MP3`](@ref) method (liquid, rain, ice+rime below `1e-10` → forward
Euler), and the temperature channel ALWAYS present (`1`). `T` is ~250 K, far
above the `1e-10` near-empty threshold, so it is never masked: its row of the
Rosenbrock system always participates implicitly, which is the entire point of
this mode.
"""
@inline function _rosenbrock_channel_mask(x::MicroState2MP3T{FT}) where {FT}
    ϵ_empty = FT(1e-10)
    liq = ifelse(x.q_lcl < ϵ_empty, zero(FT), one(FT))
    rai = ifelse(x.q_rai < ϵ_empty, zero(FT), one(FT))
    ice = ifelse(x.q_ice < ϵ_empty, zero(FT), one(FT))
    return MicroState2MP3T(liq, liq, rai, rai, ice, ice, ice, ice, one(FT))
end

"""
    bulk_microphysics_tendencies(::RosenbrockAverageImplicitT, ::Microphysics2Moment,
        mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub = 1)

Compute average 2M+P3 microphysics tendencies over `Δt` using `nsub`
linearized-implicit (Rosenbrock-Euler) substeps with the temperature promoted
into the implicitly-solved state.

# Algorithm

For each substep of `h = Δt / nsub`:

1. Evaluate the augmented tendency `f` ([`Instantaneous2MP3TTendency`](@ref)) —
   the eight species rates plus the latent-heating `dT/dt` — and its exact 9×9
   Jacobian `J` via `ForwardDiff` at the current `(species, T)` state. The
   Jacobian now carries `∂f/∂T` (the saturation feedback) and the
   latent-heating row.
2. Advance with [`_rosenbrock_update`](@ref) over the 9-vector: the kernel is
   dimension-generic, so `N = 9` plugs in unchanged. The channel projection
   ([`_rosenbrock_channel_mask`](@ref)) keeps `T` always implicit; equilibration
   `S = |x| + h|f| + ϵ` gives `T` its own ~250-scale row, so the species/T scale
   split (now spanning ~9-12 orders) stays conditioned even in Float32.

No separate between-substep `T` update is done: `T` evolves inside the solve, so
the vapor → latent-heat → temperature → saturation → exchange-flux feedback is
integrated implicitly rather than operator-split. A non-finite state or Jacobian
falls back to a forward-Euler substep of the augmented tendency. `logλ` and
`q_tot` are held fixed across substeps. The discrete safeguards are h-free
conditioning and projection (channel mask, equilibration, positivity clamp on
the species), not model terms.

Returns the net change in the species over `Δt` divided by `Δt`, in the same
fields as the [`RosenbrockAverage`](@ref) entry; the implicit `T` evolution is
internal and not returned (the macro scheme owns the energy update).
"""
@inline function bulk_microphysics_tendencies(::RosenbrockAverageImplicitT, cm::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    Δt, nsub = 1,
) where {WR, ICE <: CMP.P3IceParams}
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)

    g = Instantaneous2MP3TTendency(mp, tps, ρ, q_tot, logλ)
    x = MicroState2MP3T{FT}(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, T)
    x₀ = x
    for _ in 1:nsub_eff
        f = g(x)
        x = if all(isfinite, x)
            J = FD.jacobian(g, x)
            z = _rosenbrock_channel_mask(x)
            all(isfinite, J) ? _rosenbrock_update(x, f, J, z, h) : _euler_update(x, f, h)
        else
            _euler_update(x, f, h)
        end
    end

    rates = (x - x₀) / Δt
    return NamedTuple{(
        :dq_lcl_dt, :dn_lcl_dt, :dq_rai_dt, :dn_rai_dt,
        :dq_ice_dt, :dn_ice_dt, :dq_rim_dt, :db_rim_dt,
    )}(
        (rates.q_lcl, rates.n_lcl, rates.q_rai, rates.n_rai,
        rates.q_ice, rates.n_ice, rates.q_rim, rates.b_rim),
    )
end

#####
##### Implicit-temperature 1M Rosenbrock-Euler substepping (`RosenbrockAverageImplicitT`)
#####

"""
    MicroState1MT{FT}

The four prognostic 1M species plus the local temperature `T` as a
`StaticArrays.FieldVector`, the augmented state of the
[`RosenbrockAverageImplicitT`](@ref) 1M entry. Mirrors [`MicroState1M`](@ref)
with `T` appended as the final component, so the four-species block keeps its
ordering and the species-only kernels see the same layout.
"""
struct MicroState1MT{FT} <: SA.FieldVector{5, FT}
    q_lcl::FT
    q_icl::FT
    q_rai::FT
    q_sno::FT
    T::FT
end
SA.similar_type(::Type{<:MicroState1MT}, ::Type{FT}, ::SA.Size{(5,)}) where {FT} =
    MicroState1MT{FT}

"""
    Raw1MTTendency(mp, tps, ρ, q_tot, Lv_over_cp, Ls_over_cp)

Callable bundling the frozen per-substep context for the augmented 1M state;
applying it to a [`MicroState1MT`](@ref) evaluates the four species tendencies
([`_instantaneous_1m_tendency`](@ref)) at the state's own `T` and appends the
latent-heating tendency `dT/dt = Lv_over_cp·(dq_lcl_dt + dq_rai_dt) +
Ls_over_cp·(dq_icl_dt + dq_sno_dt)` as the fifth component. This is the rate form
of the explicit between-substep update of the [`RosenbrockAverage`](@ref) 1M
entry (which adds `Lv_over_cp·(Δq_lcl + Δq_rai) + Ls_over_cp·(Δq_icl + Δq_sno)`
from the realized increment `Δ`), so an explicit forward-Euler integration of the
augmented tendency reproduces the [`RosenbrockAverage`](@ref) 1M trajectory. The
constant `Lv_over_cp`, `Ls_over_cp` mirror the 1M convention (constant latent
heats, liquid+rain on `L_v`, ice+snow on `L_s`); the temperature dependence of
the feedback enters only through the saturation-dependent rates, exactly as in the
explicit update.

A top-level struct rather than a closure, so `ForwardDiff` differentiates a
concretely-typed callable, mirroring [`Raw1MTendency`](@ref). Unlike that functor
`T` is NOT frozen context: it is taken from the state, so the fifth-component
partial flows through the saturation specific humidities. `ρ` and `q_tot` are
frozen across substeps; `q_tot` is promoted to the state's element type at the
call so the constant carries a zero partial.
"""
struct Raw1MTTendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    q_tot::F
    Lv_over_cp::F
    Ls_over_cp::F
end
@inline function (g::Raw1MTTendency)(x::SA.StaticVector{5})
    (q_lcl, q_icl, q_rai, q_sno, T) = x
    tend = _instantaneous_1m_tendency(g.mp, g.tps,
        g.ρ, T, eltype(x)(g.q_tot),
        q_lcl, q_icl, q_rai, q_sno,
    )
    dT_dt =
        eltype(x)(g.Lv_over_cp) * (tend.dq_lcl_dt + tend.dq_rai_dt) +
        eltype(x)(g.Ls_over_cp) * (tend.dq_icl_dt + tend.dq_sno_dt)
    return MicroState1MT(tend.dq_lcl_dt, tend.dq_icl_dt, tend.dq_rai_dt, tend.dq_sno_dt, dT_dt)
end

"""
    _rosenbrock_channel_mask(x::MicroState1MT)

Diagonal of the channel projection `P` for the augmented 1M state: the four mass
channels gated exactly as the species-only [`MicroState1M`](@ref) method (below
`1e-10` → forward Euler), and the temperature channel ALWAYS present (`1`). `T`
is ~250 K, far above the near-empty threshold, so it is never masked and its row
always participates implicitly.
"""
@inline function _rosenbrock_channel_mask(x::MicroState1MT{FT}) where {FT}
    ϵ_empty = FT(1e-10)
    lcl = ifelse(x.q_lcl < ϵ_empty, zero(FT), one(FT))
    icl = ifelse(x.q_icl < ϵ_empty, zero(FT), one(FT))
    rai = ifelse(x.q_rai < ϵ_empty, zero(FT), one(FT))
    sno = ifelse(x.q_sno < ϵ_empty, zero(FT), one(FT))
    return MicroState1MT(lcl, icl, rai, sno, one(FT))
end

"""
    bulk_microphysics_tendencies(::RosenbrockAverageImplicitT, ::Microphysics1Moment,
        mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1)

Compute average 1M microphysics tendencies over `Δt` using `nsub`
linearized-implicit (Rosenbrock-Euler) substeps with the temperature promoted
into the implicitly-solved state. The 1M counterpart of the 2M+P3
[`RosenbrockAverageImplicitT`](@ref) entry, and the implicit-`T` companion of the
species-only [`RosenbrockAverage`](@ref) 1M entry.

# Algorithm

For each substep of `h = Δt / nsub`:

1. Evaluate the augmented tendency `f` ([`Raw1MTTendency`](@ref)) — the four
   species rates plus the latent-heating `dT/dt` — and its exact 5×5 Jacobian `J`
   via `ForwardDiff` at the current `(species, T)` state, so the Jacobian carries
   `∂f/∂T` and the latent-heating row.
2. Advance with [`_rosenbrock_update`](@ref) over the 5-vector (the kernel is
   dimension-generic, `N = 5`). The channel projection keeps `T` always implicit;
   equilibration gives `T` its own ~250-scale row.

No separate between-substep `T` update is done. A non-finite state or Jacobian
falls back to a forward-Euler substep; `q_tot` is held fixed across substeps.

Returns the net change in the species over `Δt` divided by `Δt`, in the same
fields as the [`RosenbrockAverage`](@ref) 1M entry; the implicit `T` evolution is
internal and not returned.
"""
@inline function bulk_microphysics_tendencies(::RosenbrockAverageImplicitT, cm::Microphysics1Moment,
    mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1,
)
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)

    g = Raw1MTTendency(mp, tps, ρ, q_tot, Lv_over_cp, Ls_over_cp)
    x = MicroState1MT{FT}(q_lcl, q_icl, q_rai, q_sno, T)
    x₀ = x
    for _ in 1:nsub_eff
        f = g(x)
        x = if all(isfinite, x)
            J = FD.jacobian(g, x)
            z = _rosenbrock_channel_mask(x)
            all(isfinite, J) ? _rosenbrock_update(x, f, J, z, h) : _euler_update(x, f, h)
        else
            _euler_update(x, f, h)
        end
    end

    rates = (x - x₀) / Δt
    return (;
        dq_lcl_dt = rates.q_lcl, dq_icl_dt = rates.q_icl,
        dq_rai_dt = rates.q_rai, dq_sno_dt = rates.q_sno,
    )
end
