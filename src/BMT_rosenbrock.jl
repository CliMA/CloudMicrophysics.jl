#####
##### 2M+P3 Rosenbrock-Euler substepping (`RosenbrockAverage`)
#####

"""
    MicroState2MP3{FT}

The eight prognostic 2M+P3 species as a `StaticArrays.FieldVector`. Internal to
the [`RosenbrockAverage`](@ref) implementation.
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
vector evaluates [`_instantaneous_2mp3_tendency`](@ref). `q_tot` is promoted to
the state's element type at the call; `logλ`, `T`, and `ρ` stay plain.
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
each a [`MicroState2MP3`](@ref) over the eight prognostic species, returned as a
`NamedTuple`. Replays the process calls of the warm-rain + P3 ice
[`bulk_microphysics_tendencies`](@ref) entry in the same order with the same
intermediate quantities, routing each process into its own species vector. The
sum over the returned processes equals the full raw tendency
[`_instantaneous_2mp3_tendency`](@ref) evaluates.

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

    # activation (cloud number only): no activation source
    dn_lcl_activation_dt = o
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

    # liquid-ice collision, ice aggregation, ice melting
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
    f23_dep = CM_HetIce.deposition_rate(
        ice_nucleation, tps, T, ρ, q_tot, q_lcl + q_rai, q_ice, n_active;
        m_nuc, τ_act, inpc_log_shift = zero(ρ),
    )
    f23_deposition = MicroState2MP3(o, o, o, o, f23_dep.∂ₜq_frz, f23_dep.∂ₜn_frz, o, o)

    # Bigg immersion freezing of cloud drops (fully-rimed embryo graupel)
    cld_bigg = CM_HetIce.liquid_freezing_rate(mp.ice.rain_freezing, pdf_c, tps, q_lcl, ρ, N_lcl, T)
    cld_cap = CM_HetIce.immersion_limit_rate(
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
only their sum. Evaluated at the primal state only, it supplies the right-hand
sides `f_p` for the linear post-solve attribution.
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
    _rosenbrock_species_mask(x)

Diagonal of the species projection matrix `P` used by
[`_rosenbrock_update`](@ref): 1 for active species, 0 for near-empty ones
(condensed mass below `1e-10`, per species: liquid, rain, ice+rime). A masked
species takes the forward-Euler update while active species stay implicit.
"""
@inline function _rosenbrock_species_mask(x::MicroState2MP3{FT}) where {FT}
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
species mask `z`, and substep `h`. Returns `(S, S⁻¹, A)`, the equilibration
matrix `S = Diagonal(|x| + h |f| + ϵ)`, its inverse, and the equilibrated
system matrix `A = I/h - S⁻¹ B S`, where `B = P J P` is the masked Jacobian and
`P = Diagonal(z)` is the species projection built from the per-scheme species
mask `z` (e.g. [`_rosenbrock_species_mask`](@ref) for 2M+P3).

The system build is separated from the solve so the full-step update and the
per-process attribution can reuse one factorization, solving against the same
`S`, `S⁻¹`, `A`.
"""
@inline function _rosenbrock_system(
    x::SA.StaticVector{N, FT}, f, J, z, h,
) where {N, FT}
    Iₙ = one(SA.SMatrix{N, N, FT})
    s = abs.(x) .+ h .* abs.(f) .+ eps(FT)
    P = Iₙ .* z'
    S = Iₙ .* s'
    S⁻¹ = Iₙ .* inv.(s)'
    B = P * J * P
    A = Iₙ / h - S⁻¹ * B * S
    return S, S⁻¹, A
end

"""
    _rosenbrock_solve(S, S⁻¹, A, v)

Solve the equilibrated Rosenbrock system from [`_rosenbrock_system`](@ref) for
the unclamped increment of right-hand side `v`: `Δ = S (A \\ (S⁻¹ v))`, the
equilibrated form of `(I/h - P J P)⁻¹ v`. Linear in `v`, so per-process
increments sum to the full-step increment.
"""
@inline _rosenbrock_solve(S, S⁻¹, A, v) = S * (A \ (S⁻¹ * v))

"""
    _rosenbrock_update(x, f, J, z, h)

One linearized-implicit (Rosenbrock-Euler) substep: build the equilibrated
system ([`_rosenbrock_system`](@ref)) for `(I/h - B) Δx = f` with the masked
Jacobian `B = P J P`, solve it ([`_rosenbrock_solve`](@ref)), and return
`max.(x + Δx, 0)`.
"""
@inline function _rosenbrock_update(
    x::SA.StaticVector{N, FT}, f, J, z, h,
) where {N, FT}
    S, S⁻¹, A = _rosenbrock_system(x, f, J, z, h)
    Δx = _rosenbrock_solve(S, S⁻¹, A, f)
    return max.(x .+ Δx, 0)
end

bulk_microphysics_tendencies(::RosenbrockAverage, ::Microphysics2Moment, args...) = throw(
    ArgumentError(
        "RosenbrockAverage on the 2M+P3 model supports only ExactJacobian; use rosenbrock_exact()",
    ),
)

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
   species to forward Euler.
3. Update the local temperature from the latent heating of the realized
   increments.

A non-finite state or Jacobian falls back to a forward-Euler substep of the raw
tendency. `logλ` and `q_tot` are held fixed across substeps.

The 2M+P3 model supports only [`ExactJacobian`](@ref); the donor-based matrix is
1M-only.

Returns the net change in the species over `Δt` divided by `Δt`, in the same
fields as the `Instantaneous` entry (without the activation diagnostic).
"""
@inline function bulk_microphysics_tendencies(mode::RosenbrockAverage{<:ExactJacobian}, cm::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    Δt, nsub = 1,
) where {WR, ICE <: CMP.P3IceParams}
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    cp_d = TDI.TD.Parameters.cp_d(tps)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / cp_d
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / cp_d

    x = MicroState2MP3{FT}(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    x₀ = x
    Tsub = T
    for _ in 1:nsub_eff
        g = Instantaneous2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ)
        x_prev = x
        if all(isfinite, x)
            f, J_raw = _tendency_and_jacobian(mode.jacobian, g, x)
            J = _apply_growth(mode.growth, J_raw)
            z = _species_mask(mode.jacobian, mode.growth)(x)
            d = if all(isfinite, J)
                _rosenbrock_update(x, f, J, z, h) - x
            else
                _euler_update(x, f, h) - x
            end
            d = _apply_limiter(mode.limiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps)
            x = max.(x .+ d, 0)
        else
            f = g(x)
            x = _euler_update(x, f, h)
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

The four prognostic 1M species as a `StaticArrays.FieldVector`. Internal to the
[`RosenbrockAverage`](@ref) implementation, mirroring [`MicroState2MP3`](@ref)
at dimension 4 (1M carries no number species).
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
vector evaluates [`_instantaneous_1m_tendency`](@ref), mirroring
[`Instantaneous2MP3Tendency`](@ref). `q_tot` is promoted to the state's element
type at the call; `T` and `ρ` stay plain.
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
of these per-process contribution vectors. The signs match
[`_aggregate_tendencies`](@ref), so the sum over the returned processes equals
the aggregated raw tendency the [`Raw1MTendency`](@ref) functor evaluates. Used
by the verbose post-solve attribution to supply the per-process right-hand sides
`f_p`.
"""
@inline function _per_process_1m(src)
    FT = typeof(src.S_phase_change_vap_lcl)
    o = zero(FT)
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
sum. Evaluated at the primal state only, it supplies the right-hand sides `f_p`
for the linear post-solve attribution.
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
    _rosenbrock_species_mask(x::MicroState1M)

Diagonal of the species projection `P` for the 1M state: 1 for active species,
0 for near-empty ones (mass below `1e-10`, per species: liquid, ice, rain,
snow). See the [`MicroState2MP3`](@ref) method for the role of `P` in
[`_rosenbrock_update`](@ref).
"""
@inline function _rosenbrock_species_mask(x::MicroState1M{FT}) where {FT}
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
tendency. The [`Jacobian`](@ref), [`GrowthTreatment`](@ref), and
[`TendencyLimiter`](@ref) options of `mode` select the substep matrix, the
growth-diagonal treatment, and the increment limiter.

# Algorithm

For each substep of `h = Δt / nsub`:

1. Build the substep Jacobian `J` from the [`Jacobian`](@ref) provider and apply
   the [`GrowthTreatment`](@ref) to it.
2. Advance with [`_rosenbrock_update`](@ref): solve `(I/h - P J P) Δx = f` in
   equilibrated variables, where the projection `P` is selected by
   [`_species_mask`](@ref).
3. Apply the [`TendencyLimiter`](@ref) to the increment.
4. Update the local temperature from the latent heating of the realized
   increments (constant latent heats, liquid+rain on `L_v`, ice+snow on `L_s`).

A non-finite state or Jacobian falls back to a forward-Euler substep of the raw
tendency; `q_tot` is held fixed across substeps.

Returns the net change in the species over `Δt` divided by `Δt`, in the same
fields as the `Instantaneous` entry.
"""
@inline bulk_microphysics_tendencies(mode::RosenbrockAverage, ::Microphysics1Moment,
    mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1) =
    _rosenbrock_average_1m(mode, mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub)

# ---- Shared 1M driver + Jacobian/mask/growth/limiter resolvers ----

"""
    _jacobian_provider(jacobian)

The substep Jacobian provider `(g, x) -> SMatrix` for a [`Jacobian`](@ref)
option: [`_jacobian_1m_linearized`](@ref) for [`DonorJacobian`](@ref),
[`_jacobian_1m_relinearized`](@ref) for [`CoupledDonorJacobian`](@ref), and
[`_ad_jacobian_1m`](@ref) for [`ExactJacobian`](@ref).
"""
@inline _jacobian_provider(::DonorJacobian) = _jacobian_1m_linearized
@inline _jacobian_provider(::CoupledDonorJacobian) = _jacobian_1m_relinearized
@inline _jacobian_provider(::ExactJacobian) = _ad_jacobian_1m

"""
    _species_mask(jacobian, growth)

The species projection `x -> z` for a [`Jacobian`](@ref) and
[`GrowthTreatment`](@ref) pair. The donor-based matrices are bounded by their rate
flooring, so every species stays implicit ([`_full_species_mask`](@ref)). An
explicit growth diagonal removes the unbounded growth of the exact Jacobian, so
its species stay implicit too; the exact Jacobian with implicit growth uses the
near-empty mask ([`_rosenbrock_species_mask`](@ref)).
"""
@inline _species_mask(::DonorJacobian, ::GrowthTreatment) = _full_species_mask
@inline _species_mask(::CoupledDonorJacobian, ::GrowthTreatment) = _full_species_mask
@inline _species_mask(::ExactJacobian, ::ExplicitGrowthDiagonal) = _full_species_mask
@inline _species_mask(::ExactJacobian, ::ImplicitGrowth) = _rosenbrock_species_mask

"""
    _apply_growth(growth, J)

Apply a [`GrowthTreatment`](@ref) to the Jacobian `J`. [`ImplicitGrowth`](@ref)
returns `J` unchanged; [`ExplicitGrowthDiagonal`](@ref) removes the positive
diagonal entries, leaving the off-diagonals and the negative diagonals.
"""
@inline _apply_growth(::ImplicitGrowth, J) = J
@inline function _apply_growth(::ExplicitGrowthDiagonal, J::SA.SMatrix{N, N, FT}) where {N, FT}
    Iₙ = one(SA.SMatrix{N, N, FT})
    return J - max.(Iₙ .* J, zero(FT))
end

"""
    _apply_limiter(limiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps)

Limit the substep increment `d` at state `x`. [`NoLimiter`](@ref) returns `d`.
[`EndStateSaturationAdjustment`](@ref), for a cell at or above saturation over its
more-supersaturated phase whose full-increment end state would drop below it, scales
`d` by `s ∈ [0, 1]` to keep that latent-heated end state at or above saturation over
that phase; otherwise it returns `d`.
"""
@inline _apply_limiter(::NoLimiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps) = d

"""
    _saturation_bisection_count(FT)

Number of bisection iterations to resolve a fraction in `[0, 1]` to the
precision of `FT`.
"""
@inline _saturation_bisection_count(::Type{FT}) where {FT} = ceil(Int, -log2(eps(FT)))

@inline function _apply_limiter(::EndStateSaturationAdjustment,
    x::MicroState1M{FT}, d::MicroState1M{FT},
    ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps,
) where {FT}
    Ssat(xx, TT) = max(
        TDI.supersaturation_over_ice(tps, q_tot, xx.q_lcl + xx.q_rai, xx.q_icl + xx.q_sno, ρ, TT),
        TDI.supersaturation_over_liquid(tps, q_tot, xx.q_lcl + xx.q_rai, xx.q_icl + xx.q_sno, ρ, TT),
    )
    latent(dd) = Lv_over_cp * (dd.q_lcl + dd.q_rai) + Ls_over_cp * (dd.q_icl + dd.q_sno)
    xf = max.(x .+ d, 0)
    if Ssat(x, Tsub) >= 0 && Ssat(xf, Tsub + latent(xf .- x)) < 0
        lo = zero(FT)
        hi = one(FT)
        for _ in 1:_saturation_bisection_count(FT)
            s = (lo + hi) / 2
            xs = max.(x .+ s .* d, 0)
            if Ssat(xs, Tsub + latent(xs .- x)) >= 0
                lo = s
            else
                hi = s
            end
        end
        return lo .* d
    end
    return d
end

@inline function _apply_limiter(::EndStateSaturationAdjustment,
    x::MicroState2MP3{FT}, d::MicroState2MP3{FT},
    ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps,
) where {FT}
    Ssat(xx, TT) = max(
        TDI.supersaturation_over_ice(tps, q_tot, xx.q_lcl + xx.q_rai, xx.q_ice, ρ, TT),
        TDI.supersaturation_over_liquid(tps, q_tot, xx.q_lcl + xx.q_rai, xx.q_ice, ρ, TT),
    )
    latent(dd) = Lv_over_cp * (dd.q_lcl + dd.q_rai) + Ls_over_cp * dd.q_ice
    xf = max.(x .+ d, 0)
    if Ssat(x, Tsub) >= 0 && Ssat(xf, Tsub + latent(xf .- x)) < 0
        lo = zero(FT)
        hi = one(FT)
        for _ in 1:_saturation_bisection_count(FT)
            s = (lo + hi) / 2
            xs = max.(x .+ s .* d, 0)
            if Ssat(xs, Tsub + latent(xs .- x)) >= 0
                lo = s
            else
                hi = s
            end
        end
        return lo .* d
    end
    return d
end

"Exact ForwardDiff Jacobian provider for [`_rosenbrock_average_1m`](@ref)."
@inline _ad_jacobian_1m(g, x) = FD.jacobian(g, x)

"""
    _tendency_and_jacobian(jacobian, g, x)

The raw substep tendency `f = g(x)` and the substep Jacobian (before the growth
treatment) for a [`Jacobian`](@ref) option, returned as `(f, J)`.

For [`ExactJacobian`](@ref) the primal `f` and the `N×N` Jacobian are both
obtained from `ForwardDiff`. The donor-based matrices ([`DonorJacobian`](@ref),
[`CoupledDonorJacobian`](@ref)) produce no tendency by-product, so `f = g(x)` is
evaluated separately.
"""
@inline function _tendency_and_jacobian(::ExactJacobian{C}, g, x::SA.FieldVector{N, FT}) where {C, N, FT}
    Tag = typeof(FD.Tag(g, FT))
    J = SA.MMatrix{N, N, FT}(undef)
    fv = SA.MVector{N, FT}(undef)
    for c in 1:(N ÷ C)
        c₀ = (c - 1) * C
        dx = SA.SVector(
            ntuple(
                i -> FD.Dual{Tag}(x[i], ntuple(s -> ifelse(c₀ + s == i, one(FT), zero(FT)), Val(C))...),
                Val(N),
            ),
        )
        y = g(dx)
        for row in 1:N
            c == 1 && (@inbounds fv[row] = FD.value(y[row]))
            for s in 1:C
                @inbounds J[row, c₀ + s] = FD.partials(y[row], s)
            end
        end
    end
    return typeof(x)(ntuple(i -> @inbounds(fv[i]), Val(N))...), SA.SMatrix(J)
end
@inline _tendency_and_jacobian(::DonorJacobian, g, x) = (g(x), _jacobian_1m_linearized(g, x))
@inline _tendency_and_jacobian(::CoupledDonorJacobian, g, x) =
    (g(x), _jacobian_1m_relinearized(g, x))

"""
    _full_species_mask(x)

The all-ones species projection: every species stays in the implicit
solve.
"""
@inline _full_species_mask(x::SA.StaticVector{N, FT}) where {N, FT} =
    ones(SA.SVector{N, FT})

"""
    _rosenbrock_average_1m(mode, mp, tps, ρ, T, q_tot,
                           q_lcl, q_icl, q_rai, q_sno, Δt, nsub)

Shared 1M Rosenbrock-average substep driver. The substep loop, equilibration,
positivity clamp, and explicit between-substep `T` update are fixed; the
[`RosenbrockAverage`](@ref) `mode` selects the Jacobian, the growth treatment,
and the increment limiter through [`_jacobian_provider`](@ref),
[`_species_mask`](@ref), [`_apply_growth`](@ref), and [`_apply_limiter`](@ref).
"""
@inline function _rosenbrock_average_1m(
    mode::RosenbrockAverage, mp::CMP.Microphysics1MParams, tps, ρ, T, q_tot,
    q_lcl, q_icl, q_rai, q_sno, Δt, nsub,
)
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)
    mask = _species_mask(mode.jacobian, mode.growth)

    x = MicroState1M{FT}(q_lcl, q_icl, q_rai, q_sno)
    x₀ = x
    Tsub = T
    for _ in 1:nsub_eff
        g = Raw1MTendency(mp, tps, ρ, Tsub, q_tot)
        x_prev = x
        if all(isfinite, x)
            f, J_raw = _tendency_and_jacobian(mode.jacobian, g, x)
            J = _apply_growth(mode.growth, J_raw)
            z = mask(x)
            d = if all(isfinite, J)
                _rosenbrock_update(x, f, J, z, h) - x
            else
                _euler_update(x, f, h) - x
            end
            d = _apply_limiter(mode.limiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps)
            x = max.(x .+ d, 0)
        else
            f = g(x)
            x = _euler_update(x, f, h)
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

"""
    _jacobian_1m(FT; lcl_lcl, lcl_icl, ..., sno_sno)

Assemble the 4×4 Jacobian of the 1M tendency over the state
`(q_lcl, q_icl, q_rai, q_sno)` from named species-pair entries: keyword
`<receiver>_<donor>` is `∂(dq_receiver/dt)/∂q_donor`. Entries not supplied are
zero. Centralizes the index layout so a hand-built Jacobian is filled by physical
coupling name rather than numeric `[i, j]` position.
"""
@inline _jacobian_1m(::Type{FT};
    lcl_lcl = zero(FT), lcl_icl = zero(FT), lcl_rai = zero(FT), lcl_sno = zero(FT),
    icl_lcl = zero(FT), icl_icl = zero(FT), icl_rai = zero(FT), icl_sno = zero(FT),
    rai_lcl = zero(FT), rai_icl = zero(FT), rai_rai = zero(FT), rai_sno = zero(FT),
    sno_lcl = zero(FT), sno_icl = zero(FT), sno_rai = zero(FT), sno_sno = zero(FT),
) where {FT} = SA.SMatrix{4, 4, FT}(
    # column-major: each column is a donor q_*, each row a receiver dq_*/dt
    lcl_lcl, icl_lcl, rai_lcl, sno_lcl,  # ∂/∂q_lcl
    lcl_icl, icl_icl, rai_icl, sno_icl,  # ∂/∂q_icl
    lcl_rai, icl_rai, rai_rai, sno_rai,  # ∂/∂q_rai
    lcl_sno, icl_sno, rai_sno, sno_sno,  # ∂/∂q_sno
)

"""
    _jacobian_1m_linearized(g::Raw1MTendency, x::MicroState1M)

Donor-based Jacobian provider for [`_rosenbrock_average_1m`](@ref): the donor-based
linearized system matrix `M` that [`LinearizedAverage`](@ref) uses (built by
`_linearize`), assembled by name via [`_jacobian_1m`](@ref).

`M` is the donor-based linearization, not the exact derivative of the raw
tendency: each donor to receiver transfer is linearized only in its donor species,
so collector-species couplings are absent. It is evaluated at the same
`(ρ, Tsub, q_tot, x)` as the raw tendency `f = g(x)`.
"""
@inline function _jacobian_1m_linearized(g::Raw1MTendency, x::MicroState1M{FT}) where {FT}
    (; q_lcl, q_icl, q_rai, q_sno) = x
    src = _microphysics_source_terms(
        Microphysics1Moment(), g.mp, g.tps, g.ρ, g.T, FT(g.q_tot),
        q_lcl, q_icl, q_rai, q_sno,
    )
    q_min = TDI.TD.Parameters.q_min(g.tps)
    M = _linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min)
    return _jacobian_1m(FT;
        lcl_lcl = M.M11, lcl_icl = M.M12,
        icl_icl = M.M22,
        rai_lcl = M.M31, rai_rai = M.M33, rai_sno = M.M34,
        sno_lcl = M.M41, sno_icl = M.M42, sno_rai = M.M43, sno_sno = M.M44,
    )
end

"""
    _vapor_exchange_rates(g::Raw1MTendency, x::MicroState1M, q_tot)

The four vapor-to-species phase-change rates `(S_phase_change_vap_lcl,
S_phase_change_vap_icl, S_phase_change_vap_rai, S_phase_change_vap_sno)` of the
1M source terms at state `x` and total water `q_tot`, as a `MicroState1M` in the
`(lcl, icl, rai, sno)` receiver order. Used to differentiate the vapor-exchange
coupling with respect to `q_tot`.
"""
@inline function _vapor_exchange_rates(g::Raw1MTendency, x::MicroState1M, q_tot)
    (; q_lcl, q_icl, q_rai, q_sno) = x
    src = _microphysics_source_terms(
        Microphysics1Moment(), g.mp, g.tps, g.ρ, g.T, q_tot,
        q_lcl, q_icl, q_rai, q_sno,
    )
    return MicroState1M(
        src.S_phase_change_vap_lcl, src.S_phase_change_vap_icl,
        src.S_phase_change_vap_rai, src.S_phase_change_vap_sno,
    )
end

"""
    _jacobian_1m_relinearized(g::Raw1MTendency, x::MicroState1M)

Coupled donor-based Jacobian provider for [`_rosenbrock_average_1m`](@ref): the
donor-based matrix [`_jacobian_1m_linearized`](@ref) with the vapor-competition
(Wegener-Bergeron-Findeisen) coupling restored.

The coupling is approximated by its vapor part: each rate depends on the
vapor specific content `q_vap = q_tot - q_lcl - q_icl - q_rai - q_sno`, giving a
condensate derivative `-∂(rate)/∂q_vap`, obtained from the derivative of the
vapor-to-species rates ([`_vapor_exchange_rates`](@ref)) with respect to `q_tot`
and added to each receiver row. The direct condensate dependence of the rates
(for example rain ventilation and the condensate availability terms) is not
recovered; use [`ExactJacobian`](@ref) for the full derivative.
"""
@inline function _jacobian_1m_relinearized(g::Raw1MTendency, x::MicroState1M{FT}) where {FT}
    Jdonor = _jacobian_1m_linearized(g, x)
    dS_dq_tot = FD.derivative(qt -> _vapor_exchange_rates(g, x, qt), FT(g.q_tot))
    wbf = -dS_dq_tot
    coupling = SA.SMatrix{4, 4, FT}(
        wbf[1], wbf[2], wbf[3], wbf[4],
        wbf[1], wbf[2], wbf[3], wbf[4],
        wbf[1], wbf[2], wbf[3], wbf[4],
        wbf[1], wbf[2], wbf[3], wbf[4],
    )
    return Jdonor + coupling
end

#####
##### Verbose post-solve per-process attribution (`Verbose`)
#####

"""
    _rosenbrock_substep_verbose(g, g_verbose, J, z, x, h)

One Rosenbrock-Euler substep with post-solve per-process attribution. Returns
`(x_new, Δx_processes, Δx_clamp)`:

- `x_new` — the realized next state, identical to the non-verbose
  [`_rosenbrock_update`](@ref) / [`_euler_update`](@ref) at the same inputs.
- `Δx_processes` — a `NamedTuple` of per-process realized increments `Δx_p`
  (each a state vector), keyed by the verbose functor `g_verbose`'s processes.
- `Δx_clamp` — the positivity clamp correction `(x_new − x) − Σ_p Δx_p` as a
  state vector.

The substep update is linear in the raw tendency, so each process's
contribution `f_p` (from `g_verbose`, summing to `g(x)`) pushed through the same
solve gives `Σ_p Δx_p` equal to the unclamped increment `Δx`. The positivity
clamp is not linear, so its correction `Δx_clamp` is returned separately. `J`
and `z` are the substep Jacobian (with the growth treatment applied) and the
species mask.

Not `@inline`d: its own specialization keeps the verbose functor call analyzed
with the concrete `MicroState{N, FT}` state type, avoiding a JET tuple-broadcast
false report in the shared P3 state constructor under inlining.
"""
function _rosenbrock_substep_verbose(g, g_verbose, J, z, x::SA.StaticVector{N, FT}, h) where {N, FT}
    f = g(x)
    fp = g_verbose(x)
    if all(isfinite, x) && all(isfinite, J)
        S, S⁻¹, A = _rosenbrock_system(x, f, J, z, h)
        Δx = _rosenbrock_solve(S, S⁻¹, A, f)
        Δxp = map(fp_i -> _rosenbrock_solve(S, S⁻¹, A, fp_i), fp)
        x_new = max.(x .+ Δx, 0)
        return x_new, Δxp, (x_new - x) - Δx
    end
    Δx = h .* f
    Δxp = map(fp_i -> h .* fp_i, fp)
    x_new = max.(x .+ Δx, 0)
    return x_new, Δxp, (x_new - x) - Δx
end

"""
    _per_process_zero_accumulator(g_verbose, x)

Zero-valued per-process accumulator matching the `NamedTuple` shape the verbose
functor `g_verbose` returns at state `x`: each process slot set to `zero(x)`.
Initializes the per-substep accumulation in the verbose averaged entries. Not
`@inline`d, for the same reason as [`_rosenbrock_substep_verbose`](@ref).
"""
function _per_process_zero_accumulator(g_verbose, x::SA.StaticVector)
    return map(_ -> zero(x), g_verbose(x))
end

"""
    bulk_microphysics_tendencies(v::Verbose, ::Microphysics2Moment,
        mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub = 1)

Diagnostic 2M+P3 Rosenbrock averaged tendency with post-solve per-process
attribution. Runs the substep loop of the wrapped
[`RosenbrockAverage`](@ref) mode and accumulates the per-process realized
increments ([`_per_process_2mp3`](@ref)) and the positivity clamp correction
through [`_rosenbrock_substep_verbose`](@ref). The per-process attribution uses
the unlimited solve; the increment limiter is not applied here.

Returns a `NamedTuple` with:

- `dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt, dq_ice_dt, dn_ice_dt, dq_rim_dt,
  db_rim_dt` — the net averaged tendencies.
- `processes` — a `NamedTuple` of per-process realized averaged tendencies
  (accumulated `Σ_substeps Δx_p / Δt`), each a [`MicroState2MP3`](@ref).
- `clamp_correction` — the non-attributable positivity-clamp tendency
  (accumulated `Σ_substeps Δx_clamp / Δt`), a [`MicroState2MP3`](@ref).

By construction `Σ_p processes_p + clamp_correction` equals the net averaged
state change `(x − x₀) / Δt` to the roundoff of the per-substep linear solve.
This is a diagnostic path, separate from the non-verbose entry.
"""
@inline function bulk_microphysics_tendencies(
    v::Verbose{<:RosenbrockAverage{<:ExactJacobian}}, cm::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    Δt, nsub = 1,
) where {WR, ICE <: CMP.P3IceParams}
    FT = typeof(q_tot)
    mode = v.mode
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    cp_d = TDI.TD.Parameters.cp_d(tps)

    x = MicroState2MP3{FT}(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    x₀ = x
    Tsub = T
    Δxp_sum = _per_process_zero_accumulator(Verbose2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ), x)
    Δx_clamp_sum = zero(x)
    for _ in 1:nsub_eff
        g = Instantaneous2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ)
        gv = Verbose2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ)
        J = _apply_growth(mode.growth, FD.jacobian(g, x))
        z = _species_mask(mode.jacobian, mode.growth)(x)
        x_prev = x
        x, Δxp, Δx_clamp = _rosenbrock_substep_verbose(g, gv, J, z, x, h)
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
    bulk_microphysics_tendencies(v::Verbose, ::Microphysics1Moment,
        mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1)

Diagnostic 1M Rosenbrock averaged tendency with post-solve per-process
attribution. Runs the substep loop of the wrapped [`RosenbrockAverage`](@ref)
mode and accumulates the per-process realized increments
([`_per_process_1m`](@ref)) and the positivity clamp correction through
[`_rosenbrock_substep_verbose`](@ref). The per-process attribution uses the
unlimited solve; the increment limiter is not applied here.

Returns a `NamedTuple` with:

- `dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt` — the net averaged tendencies.
- `processes` — a `NamedTuple` of per-process realized averaged tendencies
  (accumulated `Σ_substeps Δx_p / Δt`), each a [`MicroState1M`](@ref).
- `clamp_correction` — the non-attributable positivity-clamp tendency
  (accumulated `Σ_substeps Δx_clamp / Δt`), a [`MicroState1M`](@ref).

By construction `Σ_p processes_p + clamp_correction` equals the net averaged
state change `(x − x₀) / Δt` to the roundoff of the per-substep linear solve.
This is a diagnostic path, separate from the non-verbose entry.
"""
@inline function bulk_microphysics_tendencies(
    v::Verbose{<:RosenbrockAverage}, cm::Microphysics1Moment,
    mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1,
)
    FT = typeof(q_tot)
    mode = v.mode
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)
    jacobian = _jacobian_provider(mode.jacobian)
    mask = _species_mask(mode.jacobian, mode.growth)

    x = MicroState1M{FT}(q_lcl, q_icl, q_rai, q_sno)
    x₀ = x
    Tsub = T
    Δxp_sum = _per_process_zero_accumulator(Verbose1MTendency(mp, tps, ρ, Tsub, q_tot), x)
    Δx_clamp_sum = zero(x)
    for _ in 1:nsub_eff
        g = Raw1MTendency(mp, tps, ρ, Tsub, q_tot)
        gv = Verbose1MTendency(mp, tps, ρ, Tsub, q_tot)
        J = _apply_growth(mode.growth, jacobian(g, x))
        z = mask(x)
        x_prev = x
        x, Δxp, Δx_clamp = _rosenbrock_substep_verbose(g, gv, J, z, x, h)
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
