#####
##### The 2M+P3 manual Jacobian
#####

"""
    _jacobian_2mp3(FT; lcl_lcl, lcl_rai, ..., brim_brim)

Assemble the 8×8 Jacobian of the 2M+P3 tendency over the state
`(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)` from named
species-pair entries: keyword `<receiver>_<donor>` is
`∂(d<receiver>/dt)/∂<donor>`. The species short names are
`lcl`/`nlcl`/`rai`/`nrai`/`ice`/`nice`/`rim`/`brim`. Entries not supplied are
zero. Centralizes the index layout so the hand-built Jacobian is filled by
physical coupling name rather than numeric `[i, j]` position, mirroring the
1M assembler `_jacobian_1m`.
"""
@inline _jacobian_2mp3(::Type{FT};
    lcl_lcl = zero(FT), lcl_nlcl = zero(FT), lcl_rai = zero(FT), lcl_nrai = zero(FT),
    lcl_ice = zero(FT), lcl_nice = zero(FT), lcl_rim = zero(FT), lcl_brim = zero(FT),
    nlcl_lcl = zero(FT), nlcl_nlcl = zero(FT), nlcl_rai = zero(FT), nlcl_nrai = zero(FT),
    nlcl_ice = zero(FT), nlcl_nice = zero(FT), nlcl_rim = zero(FT), nlcl_brim = zero(FT),
    rai_lcl = zero(FT), rai_nlcl = zero(FT), rai_rai = zero(FT), rai_nrai = zero(FT),
    rai_ice = zero(FT), rai_nice = zero(FT), rai_rim = zero(FT), rai_brim = zero(FT),
    nrai_lcl = zero(FT), nrai_nlcl = zero(FT), nrai_rai = zero(FT), nrai_nrai = zero(FT),
    nrai_ice = zero(FT), nrai_nice = zero(FT), nrai_rim = zero(FT), nrai_brim = zero(FT),
    ice_lcl = zero(FT), ice_nlcl = zero(FT), ice_rai = zero(FT), ice_nrai = zero(FT),
    ice_ice = zero(FT), ice_nice = zero(FT), ice_rim = zero(FT), ice_brim = zero(FT),
    nice_lcl = zero(FT), nice_nlcl = zero(FT), nice_rai = zero(FT), nice_nrai = zero(FT),
    nice_ice = zero(FT), nice_nice = zero(FT), nice_rim = zero(FT), nice_brim = zero(FT),
    rim_lcl = zero(FT), rim_nlcl = zero(FT), rim_rai = zero(FT), rim_nrai = zero(FT),
    rim_ice = zero(FT), rim_nice = zero(FT), rim_rim = zero(FT), rim_brim = zero(FT),
    brim_lcl = zero(FT), brim_nlcl = zero(FT), brim_rai = zero(FT), brim_nrai = zero(FT),
    brim_ice = zero(FT), brim_nice = zero(FT), brim_rim = zero(FT), brim_brim = zero(FT),
) where {FT} = SA.SMatrix{8, 8, FT}(
    # column-major: each column is a donor, each row a receiver d<recv>/dt
    lcl_lcl, nlcl_lcl, rai_lcl, nrai_lcl, ice_lcl, nice_lcl, rim_lcl, brim_lcl,       # ∂/∂q_lcl
    lcl_nlcl, nlcl_nlcl, rai_nlcl, nrai_nlcl, ice_nlcl, nice_nlcl, rim_nlcl, brim_nlcl, # ∂/∂n_lcl
    lcl_rai, nlcl_rai, rai_rai, nrai_rai, ice_rai, nice_rai, rim_rai, brim_rai,       # ∂/∂q_rai
    lcl_nrai, nlcl_nrai, rai_nrai, nrai_nrai, ice_nrai, nice_nrai, rim_nrai, brim_nrai, # ∂/∂n_rai
    lcl_ice, nlcl_ice, rai_ice, nrai_ice, ice_ice, nice_ice, rim_ice, brim_ice,       # ∂/∂q_ice
    lcl_nice, nlcl_nice, rai_nice, nrai_nice, ice_nice, nice_nice, rim_nice, brim_nice, # ∂/∂n_ice
    lcl_rim, nlcl_rim, rai_rim, nrai_rim, ice_rim, nice_rim, rim_rim, brim_rim,       # ∂/∂q_rim
    lcl_brim, nlcl_brim, rai_brim, nrai_brim, ice_brim, nice_brim, rim_brim, brim_brim, # ∂/∂b_rim
)

"""
    _condevap_derivs(τ, sat_excess, Γ, cp_air, L, dqs_dT, q_limit, limit_is_ice,
        dcp_dliq, dcp_dice, dlog_τ_dq_limit = zero(τ))

The three closed-form derivatives `(∂s_liq, ∂s_rai, ∂s_ice)` of a relaxation
condensate tendency with respect to the liquid donor, the rain donor, and the
ice donor. The tendency is the substep march's clamped relaxation

```
∂ₜq = max(sat_excess, -max(0, q_limit)) / (τ·Γ)
```

so the limited branch is the one on which the donor bound is the larger of the
two arguments, and that is the branch test taken below. The inner
`max(0, q_limit)` is kept as its own expression, as it is in the march, because
it carries the semantic meaning: it is the donor mass actually available to
evaporate or sublimate.

Passing `Γ = 1` (with `dcp_dliq = dcp_dice = 0`, so the now-inert Γ-coupling
terms below vanish identically) linearizes the BARE rate `sat_excess/τ` that the
2M/P3 primal computes at its two phase-change call sites; passing the real `Γ`
linearizes the FOLDED rate the shared `CMNonEq` helper computes for the 1M twin.
The same function serves both, since Γ is only ever a multiplicative factor on
`τ`, never a control-flow branch.

`q_limit` is the donor's own mass (`q_lcl` for cloud, `q_ice` for ice, selected
by `limit_is_ice`); `limit_binds` selects between the vapor and limited branches
of `∂ₜq`.

`dlog_τ_dq_limit`, `(1/τ)·(∂τ/∂q_limit)`, is the donor's own fractional rate of
change of `τ` with `q_limit`, and is zero where `τ` has no such dependence. It
contributes a SELF-entry-only correction (`-rate·dlog_τ_dq_limit`, `rate` the
same value the primal returns on the active branch), the second half of the
`∂(1/(τΓ))/∂q` chain rule that the Γ-coupling terms below carry the first half
of.
"""
@inline function _condevap_derivs(
    τ, sat_excess, Γ, cp_air, L, dqs_dT, q_limit, limit_is_ice, dcp_dliq, dcp_dice,
    dlog_τ_dq_limit = zero(τ),
)
    FT = typeof(sat_excess)
    τΓ = τ * Γ
    # ∂(1/(τΓ))/∂q = -(τ/(τΓ)²)·∂Γ/∂q, ∂Γ/∂q = -(L·dqs_dT/cp²)·∂cp/∂q
    dΓ_dliq = -(L * dqs_dT / cp_air^2) * dcp_dliq
    dΓ_dice = -(L * dqs_dT / cp_air^2) * dcp_dice
    dinv_dliq = -dΓ_dliq * τ / τΓ^2
    dinv_dice = -dΓ_dice * τ / τΓ^2
    # The march's own clamp, in the march's own form: the limited branch is the one where
    # `-max(0, q_limit)` is the larger argument of the `max`.
    q_available = max(zero(FT), q_limit)
    limit_binds = sat_excess < -q_available
    # Vapor branch: ∂ₜq = sat_excess/(τΓ), ∂sat_excess/∂q = -1 on every donor.
    v_liq = -1 / τΓ + sat_excess * dinv_dliq
    v_ice = -1 / τΓ + sat_excess * dinv_dice
    # Limited branch: ∂ₜq = -q_limit/(τΓ). Only the limit species (q_lcl for
    # cloud, q_ice for ice) carries the -1/(τΓ) self term; q_rai is never the
    # limit. Every donor keeps the -q_limit·∂(1/(τΓ)) Γ term.
    limit_self = 1 / τΓ
    c_liq = -q_limit * dinv_dliq - ifelse(limit_is_ice, zero(FT), limit_self)
    c_rai = -q_limit * dinv_dliq
    c_ice = -q_limit * dinv_dice - ifelse(limit_is_ice, limit_self, zero(FT))
    # τ(q) correction: d(N·κ)/dq_limit with κ = 1/(τΓ) picks up an extra -(τ'/τ)·κ term
    # beyond the Γ-only chain rule above whenever τ' ≠ 0, and `N` times that is
    # -dlog_τ_dq_limit·rate, `rate` the value the primal returns on the active branch. SELF
    # entry only: `τ` depends on the limit species' own mass, not on the other two donors.
    rate = ifelse(limit_binds, -q_limit, sat_excess) / τΓ
    self_correction = -rate * dlog_τ_dq_limit
    ∂s_liq = ifelse(limit_binds, c_liq, v_liq) + ifelse(limit_is_ice, zero(FT), self_correction)
    ∂s_rai = ifelse(limit_binds, c_rai, v_liq)
    ∂s_ice = ifelse(limit_binds, c_ice, v_ice) + ifelse(limit_is_ice, self_correction, zero(FT))
    return (; ∂s_liq, ∂s_rai, ∂s_ice)
end

"""
    _riming_jacobian_block(rs, dlcl, drai)

The liquid-ice collision block of [`_jacobian_2mp3_manual`](@ref), donor-linearized
from the per-donor split `rs` of [`p3_2m_process_rates`](@ref) with the donor
factors `dlcl`, `drai`.

Returned as a `NamedTuple` of `<receiver>_<donor>` contributions to be accumulated
into the matrix entries of the same names.

Riming conserves condensate exactly in the primal
([`CMP3.bulk_liquid_ice_collision_sources`](@ref): the cloud, rain and ice mass
rates sum to zero), and the block reproduces that in the linearization by giving
each donor's loss to the receiver that takes it: the frozen cloud mass to ice, the
shed cloud mass to rain, the frozen rain mass to ice. So the condensate rows of each
donor column sum to zero, as they already do for every other condensate-to-condensate
transfer in the matrix, and the receivers are built as the exact negatives of the
sink terms rather than recomputed, so the cancellation is not merely algebraic.

Routing the shed mass to the *cloud* donor rather than to the net rain rate also
keeps the rain diagonal a sink: the aggregated `∂ₜq_r = -QRFRZ + QCSHD` is
sign-indeterminate and would be dropped wherever shedding exceeds rain freezing,
while the rain donor's own loss `∂ₜq_rai_frz` is a sink unconditionally.
"""
@inline function _riming_jacobian_block(rs, dlcl, drai)
    D_lcl_frz = rs.∂ₜq_lcl_frz * dlcl   # ≤ 0: cloud mass onto ice
    D_lcl_shd = rs.∂ₜq_lcl_shd * dlcl   # ≤ 0: cloud mass collected and shed as rain
    D_rai_frz = rs.∂ₜq_rai_frz * drai   # ≤ 0: rain mass onto ice
    return (;
        lcl_lcl = D_lcl_frz + D_lcl_shd,
        ice_lcl = -D_lcl_frz,
        rai_lcl = -D_lcl_shd,
        rai_rai = D_rai_frz,
        ice_rai = -D_rai_frz,
        # Rime mass tracks the frozen liquid one for one; rime volume has its own
        # per-donor collection rates. Neither enters the condensate total: `q_rim` is a
        # fraction of `q_ice`, not an addition to it.
        rim_lcl = -D_lcl_frz,
        rim_rai = -D_rai_frz,
        brim_lcl = rs.∂ₜb_rim_lcl * dlcl,
        brim_rai = rs.∂ₜb_rim_rai * drai,
    )
end

"""
    _jacobian_2mp3_manual(g, x::MicroState2MP3, pp, rs)

The hand-built 2M+P3 substep Jacobian for [`ManualJacobian`](@ref), evaluated at
the same state as the tendency it linearizes.

`g` is the frozen per-substep context, any object carrying the fields `mp`, `tps`,
`ρ`, `T`, `q_tot`, `logλ`, `w` and `p`; `x` is the eight-species state; `pp` is the
per-process breakdown [`p3_2m_process_rates`](@ref) returns at that state and `rs`
the liquid-ice collision transfer split by donor that it returns alongside.

ONE QUANTITY IS SHARED BY THE TENDENCY AND BOTH JACOBIAN MODES: `pp`, the 21-slot
per-process breakdown. The tendency the substep integrates is `sum(values(pp))`;
[`ExactJacobian`](@ref) differentiates a callable that computes that sum, and
`ManualJacobian` assembles this matrix from the same `pp` at the same state. This
function evaluates no process rate of its own, so the two modes provably linearize
the tendency the substep integrates and cannot drift apart. What it does recompute
are the relaxation TIMESCALES and the closed-form derivative bundles that the rates
alone do not carry (`τ_l`, `τ_i`, the orphan inverse timescales, the mean-mass-bounded
populations and the droplet-activation derivative bundle), each at the same state and
under the same branch tests the primal takes.

The entries are tiered:

  - Tier 1 (closed-form): cloud condensation/evaporation, ice
    deposition/sublimation (with the ice-number sublimation pathway and the rim
    drain), the four number adjustments, droplet activation, and the
    deposition-nucleation number pathway. The supersaturation block carries the
    autocatalytic vapor brake `∂q_vap/∂q_condensate = −1` and the relaxation
    `−1/(τ·Γ)` couplings. Both relaxation timescales are population-dependent and
    are differentiated as such: `1/τ_l ∝ N_lcl^{2/3} q_lcl^{1/3}` from the droplet
    diameter moment and `1/τ_i ∝ n_ice` from the P3 capacitance integral at frozen
    `logλ`, which supply the two number couplings of the phase-change rows and make
    the sublimation number pathway degree-2 homogeneous in `n_ice`.
  - Tier 2 (donor-diagonal): autoconversion, accretion, cloud/rain
    self-collection, breakup, rain evaporation, immersion freezing, and rain
    freezing, linearized in their donor species through
    `D = rate / max(floor, q_donor)` reusing the per-process rates of `pp`.
  - Tier 3 (coupled donor): the mixed-phase quadrature transfers, donor-linearized
    through `D = rate / max(floor, x_donor)` with `−D` on the donor diagonal and
    `+D` on the receiver off-diagonal (the 1M `_linearize` recipe), reusing the
    per-process rates of `pp`. The ice to rain melt transfer (mass donor `q_ice`,
    number donor `n_ice`, rim mass and volume on their own donors) is the dominant
    one: without it the rain melt source integrates approximately explicitly and
    diverges with the step. Liquid-ice collision carries both sides through
    [`_riming_jacobian_block`](@ref), and ice aggregation its own quadratic diagonal.
    No `gamma_inc` shape derivative is taken: the quadrature rate is frozen and only
    the donor dependence is linearized.
"""
@inline function _jacobian_2mp3_manual(g, x::MicroState2MP3{FT}, pp, rs) where {FT}
    mp = g.mp
    tps = g.tps
    # The condensation timescale is state-dependent and reaches a `log(ρ q / N)`, so a
    # host-delivered ρ < 0 would throw a `DomainError` here, which inside a GPU kernel aborts
    # the whole kernel and masks the state that caused it. See [`AIR_DENSITY_FLOOR`](@ref) for
    # why the floor is positive rather than a non-negative clamp.
    ρ = _floored_air_density(g.ρ)
    T = g.T
    q_tot = FT(g.q_tot)
    logλ = g.logλ

    (; q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim) = x
    o = zero(FT)
    # ϵₘ matches no guard in the primal: the `n_ice/q_ice` guard is at zero, and the number
    # adjustments test `q > 0` in both the primal and `_numadj_derivs`. Its one use here is the
    # donor floor of the sublimation cross-derivatives.
    qmin = UT.ϵ_numerics_2M_M(FT)
    # donor floor for the Tier-2 linearizations (the 1M donor recipe's `q_min`)
    q_floor = FT(TDI.TD.Parameters.q_min(tps))
    n_floor = q_floor

    # --- shared thermodynamic constants (T, ρ, q_tot frozen in the substep) ---
    Rᵥ = TDI.Rᵥ(tps)
    Lᵥ = TDI.Lᵥ(tps, T)
    Lₛ = TDI.Lₛ(tps, T)
    cp_v = TDI.TD.Parameters.cp_v(tps)
    cp_l = TDI.TD.Parameters.cp_l(tps)
    cp_i = TDI.TD.Parameters.cp_i(tps)
    T_freeze = TDI.T_freeze(tps)
    cp_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_ice)
    # ∂cp_m/∂q: q_liq = q_lcl + q_rai, q_ice = q_ice (primal convention)
    dcp_dliq = cp_l - cp_v  # ∂cp_m/∂q_lcl = ∂cp_m/∂q_rai
    dcp_dice = cp_i - cp_v  # ∂cp_m/∂q_ice
    qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_ice)

    #####
    ##### Tier 1: closed-form stiff couplings
    #####

    # cloud condensation / evaporation (row q_lcl), branch matched to the primal;
    # τ matches the primal's capacitance-integral timescale
    pdf_c_j = mp.warm_rain.seifert_beheng.pdf_c
    # The zero-mass droplet arm's deactivation test, taken on the CLAMPED state the primal uses so
    # that this matrix takes the same branch its own primal does. `sat_excess_l` below is the
    # UNCLAMPED excess the condensation derivative is built on; the two coincide everywhere except
    # at negative condensate, where the rate being linearized is the clamped one and a branch test
    # has to follow the rate rather than the derivative.
    sat_excess_arm = _liquid_sat_excess(
        tps, ρ, T, q_tot,
        UT.clamp_to_nonneg(q_lcl), UT.clamp_to_nonneg(q_rai), UT.clamp_to_nonneg(q_ice))
    n_lcl_j = CM2.number_bounded_by_mass_limits(
        (; x_min = pdf_c_j.xc_min, x_max = pdf_c_j.xc_max),
        UT.clamp_to_nonneg(q_lcl), UT.clamp_to_nonneg(n_lcl), sat_excess_arm;
        invent_from_zero = false)
    τ_l = CM2.cloud_condensation_timescale(
        pdf_c_j, mp.warm_rain.air_properties, tps, T, ρ,
        UT.clamp_to_nonneg(q_lcl), n_lcl_j * ρ)
    qᵥ_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
    dqsl_dT = CMNonEq.dqcld_dT(qᵥ_sat_liq, Lᵥ, Rᵥ, T)
    sat_excess_l = qᵥ - qᵥ_sat_liq
    # The orphan drain's own branch and rate, on the CLAMPED state the primal uses, so that this
    # matrix takes the same branch its primal does at the same state.
    q_lcl_c = UT.clamp_to_nonneg(q_lcl)
    q_rai_c = UT.clamp_to_nonneg(q_rai)
    sat_excess_orphan_j =
        _liquid_sat_excess(tps, ρ, T, q_tot, q_lcl_c, q_rai_c, UT.clamp_to_nonneg(q_ice))
    lcl_is_orphan_j = !(FD.value(n_lcl) > o) & (q_lcl_c > o)
    rai_is_orphan_j = !(FD.value(n_rai) > o) & (q_rai_c > o)
    pdf_r_j = mp.warm_rain.seifert_beheng.pdf_r
    inv_τ_orphan_lcl = CM2.orphan_mass_inv_timescale(
        mp.warm_rain.air_properties, tps, T, ρ, pdf_c_j.xc_min, pdf_c_j.ρw, sat_excess_orphan_j)
    inv_τ_orphan_rai = CM2.orphan_mass_inv_timescale(
        mp.warm_rain.air_properties, tps, T, ρ, pdf_r_j.xr_min, pdf_r_j.ρw, sat_excess_orphan_j)
    # The primal computes the BARE rate `sat_excess/τ_l` at this call site, so the linearization
    # passes `Γ = 1` with the Γ-coupling inputs zeroed and `_condevap_derivs` reduces to
    # differentiating `sat_excess/τ_l` directly, matching the primal it linearizes.
    # `dlog_τ_dq_liq` carries `τ_l`'s own `q_lcl` dependence into the self entry (see
    # [`_condevap_derivs`](@ref)); it is guarded to zero where the mean-mass floor or
    # `clamp_to_nonneg` already makes `τ_l` inert to `q_lcl`, the same precondition
    # `n_lcl_unbounded` guards for `lcl_nlcl` below.
    dlog_τ_dq_liq = UT.guarded_quotient(-one(FT), 3 * q_lcl)
    cl = _condevap_derivs(τ_l, sat_excess_l, one(FT), cp_air, Lᵥ, dqsl_dT, q_lcl, false,
        zero(FT), zero(FT), dlog_τ_dq_liq)
    # A degenerate τ_l means the capacitance integral underflowed, and the primal zeroes the whole
    # `cloud_condevap` slot there. The derivative of an identically zero rate is zero, so `cl` has
    # to be dropped with it, for the reason `ci` is dropped below.
    cond_active = !CM2.cloud_condensation_is_degenerate(τ_l)
    lcl_lcl = ifelse(cond_active, cl.∂s_liq, o)
    lcl_rai = ifelse(cond_active, cl.∂s_rai, o)
    lcl_ice = ifelse(cond_active, cl.∂s_ice, o)

    # Orphan-mass drain (rows q_lcl and q_rai, own diagonals). `∂ₜq = -q/τ_orphan` with a
    # timescale that does not depend on the mass, so the entry is exactly `-1/τ_orphan`: the
    # linearization is exact rather than a donor recipe. Left explicit, a cell whose entire
    # condensate is orphaned would take its own removal as a forward-Euler step.
    #
    # The coupling through the subsaturation inside `τ_orphan` is DROPPED, deliberately: it is the
    # derivative of a rate whose whole purpose is to remove an unphysical state within a few
    # steps, its sign opposes the diagonal, and carrying it would couple three columns for a
    # correction of order `q/q_sat`.
    lcl_lcl += ifelse(lcl_is_orphan_j, -inv_τ_orphan_lcl, o)
    rai_rai_orphan = ifelse(rai_is_orphan_j, -inv_τ_orphan_rai, o)
    # The ice orphan drain's branch and rate, on the clamped state as above. The predicate is
    # folded into the value, so a populated cell contributes exactly zero; the q_ice diagonal
    # and the rime-ray self-derivatives are applied where those accumulators exist.
    q_ice_c = UT.clamp_to_nonneg(q_ice)
    ice_is_orphan_j = !(FD.value(n_ice) > o) & (q_ice_c > o)
    sat_excess_i_orphan_j = _ice_sat_excess(tps, ρ, T, q_tot, q_lcl_c, q_rai_c, q_ice_c)
    inv_τ_orphan_ice_j = ifelse(
        ice_is_orphan_j,
        CM2.orphan_mass_inv_timescale_ice(
            mp.warm_rain.air_properties, tps, T, ρ,
            CMP3.ice_mean_particle_mass_min(mp.ice.scheme), mp.ice.scheme.ρ_i,
            sat_excess_i_orphan_j),
        o,
    )
    # Droplet-number coupling of the condensation rate. `τ_l` is not a constant in `n_lcl`: it is
    # `ρ q_{v,sl}` over the diameter moment of the droplet PSD, and that moment carries
    # `∂log∫D n dD/∂log N = (ν_c + 2) − μ_c (ν_cD + 2)/μ_cD = 2/3` for every `(ν_c, μ_c)`, i.e.
    # `1/τ_l ∝ N_lcl^{2/3} q_lcl^{1/3}` (the surface a fixed mass presents grows as the drops are
    # divided). Both branches of the bare rate carry no OTHER `n_lcl` dependence, so the
    # donor-linearization recipe gives the entry exactly: homogeneity degree 2/3 times rate over
    # donor, on the vapor and the evaporation-limited branch alike. `pp.cloud_condevap.q_lcl`
    # below is read directly from the primal, so this entry is automatically consistent with it.
    #
    # The rate reads the mean-mass-BOUNDED number, so where that bound binds the primal has no
    # `n_lcl` dependence at all and the entry must vanish with it (the f/J consistency doctrine
    # applied to a clamp rather than to a gate). `n_lcl_j == n_lcl` is the bound's own inertness
    # test, and it also excludes the non-positive numbers, where `clamp_to_nonneg` is the clamp.
    n_lcl_unbounded = n_lcl_j == n_lcl
    lcl_nlcl = ifelse(
        cond_active & n_lcl_unbounded,
        2 * pp.cloud_condevap.q_lcl / (3 * max(n_lcl, floatmin(FT))),
        o,
    )

    # ice deposition / sublimation (rows q_ice, n_ice, q_rim, b_rim);
    # τ matches the primal's capacitance-integral timescale (inputs clamped
    # to nonnegative as in the primal)
    state_i = CMP3.state_from_prognostic(
        mp.ice.scheme,
        UT.clamp_to_nonneg(q_ice) * ρ, UT.clamp_to_nonneg(n_ice) * ρ,
        UT.clamp_to_nonneg(q_rim) * ρ, UT.clamp_to_nonneg(b_rim) * ρ)
    τ_i = CMP3.ice_deposition_timescale(
        mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps, T, ρ,
        state_i, logλ; quad = mp.ice.quad)
    qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
    dqsi_dT = CMNonEq.dqcld_dT(qᵥ_sat_ice, Lₛ, Rᵥ, T)
    sat_excess_i = qᵥ - qᵥ_sat_ice
    # The primal computes the BARE rate `sat_excess/τ_i` at this call site, so `Γ = 1` with the
    # Γ-coupling inputs zeroed, matching the liquid call above. The τ(q) self term is scoped to
    # the LIQUID donor only: `τ_i`'s own `n_ice` dependence is carried separately by `ice_nice`
    # below, and its `q_ice` dependence through the shape-solve output has no closed form. The
    # `zero(FT)` is passed explicitly rather than left to the default, so this call reads as a
    # stated decision rather than an oversight.
    ci = _condevap_derivs(τ_i, sat_excess_i, one(FT), cp_air, Lₛ, dqsi_dT, q_ice, true,
        zero(FT), zero(FT), zero(FT))
    # `wet_surface_deposition_limiter` zeros a positive (deposition) tendency inside the primal's
    # ice relaxation, so the derivative survives only on the sublimation branch
    # (sat_excess_i ≤ 0); that branch of `ci` is unchanged.
    dep_suppressed = (T > T_freeze) & (sat_excess_i > 0)
    # A degenerate τ_i means the capacitance integral underflowed, and the primal zeroes the
    # whole `ice_depsub` slot there. The derivative of an identically zero rate is zero, so `ci`
    # has to be dropped with it: `_condevap_derivs` at the capped timescale returns a
    # self-derivative of order −1/τ_max ≈ −5e-11, which against `f = 0` exactly is an f/J
    # inconsistency rather than a linearization.
    dep_degenerate = CMP3.ice_deposition_is_degenerate(τ_i)
    dep_active = !(dep_suppressed | dep_degenerate)
    ice_lcl = ifelse(dep_active, ci.∂s_liq, o)
    ice_rai = ifelse(dep_active, ci.∂s_rai, o)
    ice_ice = ifelse(dep_active, ci.∂s_ice, o)
    ∂ₜq_ice_dep = pp.ice_depsub.q_ice
    # Ice-number coupling of the deposition rate. `τ_i` is not a constant in `n_ice` either: below
    # the cap it is `ρ q_{v,si}` over the capacitance integral `∫D F_v N′ dD`, and at frozen `logλ`
    # the P3 size distribution carries the number as a single additive `log(n_ice)` in `log N₀`
    # while the quadrature bounds are quantiles of the shape alone, so the integral is exactly
    # linear in `n_ice` and `1/τ_i ∝ n_ice`. The bare rate carries no OTHER `n_ice` dependence, so
    # this entry is the degree-1 donor linearization, exact rather than approximate: rate over
    # donor. `∂ₜq_ice_dep` is read directly from the primal, so this entry is automatically
    # consistent with whatever it computes.
    #
    # It is POSITIVE wherever the air is supersaturated over ice, and it is off-diagonal. That is
    # the physics and not a sign error: more crystals deposit more vapor, so ice mass grows with
    # ice number, and the growth coupling belongs in the linearization the implicit update
    # navigates by. It is largest exactly where the loading is small and the supersaturation is
    # not, the fresh-crystal regime, since `rate/n_ice` is bounded by the per-crystal capacitance
    # and does not diverge as the population thins. `ExplicitGrowthDiagonal` removes positive
    # DIAGONAL entries only, so production keeps this one.
    ice_nice = ifelse(dep_active, ∂ₜq_ice_dep / max(n_ice, floatmin(FT)), o)

    # ice orphan drain (row q_ice, own diagonal): exact, as the warm-phase entries above
    ice_ice += -inv_τ_orphan_ice_j

    # ice number sublimation pathway: ∂ₜn_ice_dep = n_ice·(∂ₜq_ice_dep/q_ice), active on the
    # sublimation branch (∂ₜq_ice_dep < 0) at every positive q_ice - the primal's own guard is at
    # zero, not at `qmin`, so carrying `q_ice > qmin` here would switch the pathway off in J while
    # f still runs it.
    n_sub_active = ∂ₜq_ice_dep < 0 && dep_active
    # The self-derivative is the fractional mass loss the primal forms, bounded by 1/τ_dep
    # whatever q_ice is, so it is carried exactly. The donor cross-derivatives carry n_ice/q_ice,
    # which overflows Float32 on the number-without-mass states this domain carries, so they keep
    # the Tier-1 donor floor: the linearization under-damps there rather than returning Inf.
    sub_frac = UT.guarded_quotient(∂ₜq_ice_dep, q_ice)
    n_per_q = n_ice / max(qmin, q_ice)
    # ∂(n·∂ₜq/q)/∂q = (n/q)·∂(∂ₜq)/∂q − (n/q)·(∂ₜq/q)
    nice_lcl = ifelse(n_sub_active, n_per_q * ice_lcl, o)
    nice_rai = ifelse(n_sub_active, n_per_q * ice_rai, o)
    nice_ice = ifelse(n_sub_active, n_per_q * (ice_ice - sub_frac), o)
    # ∂(n·∂ₜq/q)/∂n = ∂ₜq/q + (n/q)·∂(∂ₜq)/∂n, and the second term is not zero: the sublimation
    # rate is itself linear in n_ice through τ_i, so `(n/q)·(∂ₜq/n) = ∂ₜq/q` again and the pathway
    # is degree-2 homogeneous in the ice number. The self-derivative is therefore 2·∂ₜq/q, twice
    # what a rate-over-donor reading of the pathway alone gives. It is carried exactly rather than
    # as `sub_frac + n_per_q·ice_nice`, because the sublimation branch bounds it by 2/τ_i at
    # any positive q_ice while `n_per_q` keeps the Tier-1 floor against Float32 overflow.
    nice_nice = ifelse(n_sub_active, 2 * sub_frac, o)

    # Rim-drain CROSS couplings on the sublimation branch are dropped from J (Tier 3), the
    # ice-number one among them: `∂ₜq_rim_sub` and `∂ₜb_rim_sub` are the ray-form multiples
    # `q_rim` and `b_rim` of the fractional ice-mass loss, so they carry the same `1/τ_i` number
    # dependence as `∂ₜq_ice_sub`, but adding `rim_nice`/`brim_nice` alone would linearize the rim
    # rows in one donor while their other three stay explicit. The two SELF-derivatives ARE
    # carried below, where that same ray form makes them exact.
    rim_lcl = o
    rim_rai = o
    rim_ice = o
    brim_lcl = o
    brim_rai = o
    brim_ice = o

    # Deposition-nucleation number pathway: the crystal number relaxes toward the closure's
    # target, `∂ₜn_frz = max(0, n_target − n_ice)/τ_act`, so the active arm contributes
    # `∂/∂n_ice = −1/τ_act` on the ice-number diagonal.
    #
    # The mass donors are treated explicitly for the reason `_numadj_derivs` drops its own
    # `1/(x_min·τ)` mass coupling: the target's mass coupling is a gain of order `1/m_nuc`, which
    # is nine orders above the implicit diagonal `I/h = 0.5`. The relaxation diagonal itself is
    # cheap by comparison and is what keeps a relaxation the substep would otherwise integrate
    # explicitly inside the implicit solve.
    #
    # The delivery rate is the inverse of the seed delivery time and depends on the state, so it
    # is evaluated here rather than read from a parameter. Its supersaturation factor is the
    # explicit part of the coupling: the rate falls linearly to zero as the vapor approaches ice
    # saturation and is exactly zero below it, but that dependence is not carried in `J`.
    S_i = TDI.supersaturation_over_ice(tps, q_tot, q_lcl + q_rai, q_ice, ρ, T)
    inv_τ_dep = CM_HetIce.delivery_rate(mp.ice.ice_nucleation, mp, tps, T, S_i)
    dep_nuc_n = pp.ice_deposition.n_ice
    dep_nuc_active = dep_nuc_n > 0
    nice_nice += ifelse(dep_nuc_active, -inv_τ_dep, o)

    # number adjustments (cloud, rain, ice): closed-form, see
    # [`CM2.number_tendency_from_mass_limits`](@ref). Interior ⇒ 0; clamped low/high ⇒
    # ∂/∂n = −1/τ, ∂/∂q = 1/(x_bound·τ).
    sb = mp.warm_rain.seifert_beheng
    (nlcl_lcl, nlcl_nlcl) = _numadj_derivs(
        FT, q_lcl, n_lcl, sb.pdf_c.xc_min, sb.pdf_c.xc_max, sb.numadj.τ, sat_excess_arm)
    (nrai_rai, nrai_nrai) =
        _numadj_derivs(FT, q_rai, n_rai, sb.pdf_r.xr_min, sb.pdf_r.xr_max, sb.numadj.τ)
    numadj_ice = _ice_numadj_params(mp.ice.scheme)
    (nice_ice_adj, nice_nice_adj) =
        _numadj_derivs(FT, q_ice, n_ice, numadj_ice.x_min, numadj_ice.x_max, numadj_ice.τ)
    nice_ice += nice_ice_adj
    nice_nice += nice_nice_adj

    # the rain number the primal evaluates its warm-rain rates at, needed by the relaxation entry
    # below (matching `n_lcl_j` above and the primal's own bounded rain number)
    n_rai_j = CM2.number_bounded_by_mass_limits(
        (; x_min = sb.pdf_r.xr_min, x_max = sb.pdf_r.xr_max),
        UT.clamp_to_nonneg(q_rai), UT.clamp_to_nonneg(n_rai))

    # droplet activation (rows n_lcl and q_lcl) contributes no entries until the warm-rain
    # parameters carry an aerosol and an activation closure.
    nlcl_rai = o
    nlcl_ice = o

    #####
    ##### Tier 2: donor-diagonal linearizations of the warm-rain and freezing transfers
    #####
    # accumulators that receive only Tier-2 contributions (Tier-1 left them zero)
    rai_lcl = o
    rai_rai = rai_rai_orphan
    nrai_nlcl = o
    # each process column = process tendency vector / max(floor, donor); mass
    # pathways routed by the mass donor, number pathways by the number donor.
    dlcl = 1 / max(q_floor, q_lcl)
    drai = 1 / max(q_floor, q_rai)
    dnlcl = 1 / max(n_floor, n_lcl)
    dnrai = 1 / max(n_floor, n_rai)

    # rain evaporation: q_rai, n_rai sinks (donors q_rai, n_rai)
    rai_rai += pp.rain_evap.q_rai * drai
    nrai_nrai += pp.rain_evap.n_rai * dnrai

    # autoconversion: mass donor q_lcl, number donor n_lcl
    lcl_lcl += pp.autoconv.q_lcl * dlcl
    rai_lcl += pp.autoconv.q_rai * dlcl
    nlcl_nlcl += pp.autoconv.n_lcl * dnlcl
    nrai_nlcl += pp.autoconv.n_rai * dnlcl

    # cloud self-collection: n_lcl sink, donor n_lcl
    nlcl_nlcl += pp.cloud_selfcol.n_lcl * dnlcl

    # accretion: mass donor q_lcl, number donor n_lcl
    lcl_lcl += pp.accretion.q_lcl * dlcl
    rai_lcl += pp.accretion.q_rai * dlcl
    nlcl_nlcl += pp.accretion.n_lcl * dnlcl

    # rain self-collection + breakup: the one process whose rate legitimately changes SIGN
    # through the population's shape rather than a thermodynamic attractor, so the donor recipe
    # `(sc + br)/n_rai` is positive - anti-damping - wherever breakup dominates. The pair is a
    # relaxation of the number toward the collisional equilibrium `n_eq = q_rai/x_eq`, and
    # written that way its frozen-shape diagonal is `-1/τ_eff ≤ 0` by construction. See
    # [`CM2.rain_number_relaxation`](@ref). Evaluated at the mean-mass-bounded number the rate
    # was, so f and the deviation are the same function's.
    #
    # Where that bound BINDS the entry is exactly zero, not merely small: the rate then depends
    # on n_rai only through `clamp`, whose derivative is zero, so the chain rule gives zero and
    # the row's damping there is the number adjustment's own `-1/τ`, which is the designated
    # restorer for that regime.
    n_rai_c = UT.clamp_to_nonneg(n_rai)
    nrai_nrai += ifelse(
        n_rai_j == n_rai_c,
        -CM2.rain_number_relaxation(
            sb.pdf_r, sb.self, sb.brek, UT.clamp_to_nonneg(q_rai), ρ, n_rai_j * ρ,
        ).inv_τ_eff,
        o,
    )

    # immersion freezing (cloud → ice): mass donor q_lcl, number donor n_lcl
    lcl_lcl += pp.immersion_freezing.q_lcl * dlcl
    ice_lcl += pp.immersion_freezing.q_ice * dlcl
    rim_lcl += pp.immersion_freezing.q_rim * dlcl
    brim_lcl += pp.immersion_freezing.b_rim * dlcl
    nlcl_nlcl += pp.immersion_freezing.n_lcl * dnlcl
    nice_nlcl = pp.immersion_freezing.n_ice * dnlcl

    # rain freezing (rain → ice): mass donor q_rai, number donor n_rai
    rai_rai += pp.rain_freezing.q_rai * drai
    ice_rai += pp.rain_freezing.q_ice * drai
    rim_rai += pp.rain_freezing.q_rim * drai
    brim_rai += pp.rain_freezing.b_rim * drai
    nrai_nrai += pp.rain_freezing.n_rai * dnrai
    nice_nrai = pp.rain_freezing.n_ice * dnrai

    #####
    ##### Tier 3: coupled donor linearization of the mixed-phase quadrature transfers
    #####
    # The dominant mixed-phase coupling is the ice→rain melt source; without it the
    # rain source integrates approximately explicitly and runs away. Donor-linearize each
    # transfer of primal rate `S` from donor `d` to receiver `r` as `D = S/max(floor, x_d)`
    # with `−D` on the donor diagonal and `+D` on the (r, d) off-diagonal (the 1M
    # `_linearize` recipe), reusing the per-process rates of `pp`. No `gamma_inc` shape
    # derivative is taken: the quadrature rate is held frozen and only the donor dependence
    # is linearized, so the receiver source self-limits as the donor empties within the
    # implicit step.
    rai_ice = o
    nrai_nice = o
    rim_rim = o
    brim_brim = o
    # The rime volume acquires a dependence on the rime mass through the melting densification
    # below; every other rime term is a rate along the ray, on which the pair moves together and
    # neither is a function of the other.
    brim_rim = o

    # ice melt (ice → rain): mass donor q_ice, number donor n_ice; the rim pair drains
    # along its own ray. The melt vector is signed +source into rain / −sink out of ice,
    # so each `pp.ice_melting.<species>` already carries the rate.
    dice = 1 / max(q_floor, q_ice)
    dnice = 1 / max(n_floor, n_ice)
    D_melt_q = pp.ice_melting.q_rai * dice   # ≥ 0
    ice_ice -= D_melt_q
    rai_ice += D_melt_q
    D_melt_n = pp.ice_melting.n_rai * dnice  # ≥ 0
    nice_nice -= D_melt_n
    nrai_nice += D_melt_n
    # The ray-form drain is `−q_rim·frac` / `−b_rim·frac` with the fractional ice-mass loss
    # `frac` frozen at the Tier-3 quadrature rate, so BOTH rim self-derivatives are exactly
    # `−frac` and neither needs a donor floor. Their being EQUAL is what makes the implicit
    # solve rescale the pair proportionally instead of reshaping it: unequal floored diagonals
    # (`1/max(q_floor, q_rim)` against `1/max(n_floor, b_rim)`, a NUMBER floor on a VOLUME
    # moment) would damp the two moments differently and so move the quotient even where `f`
    # preserves it.
    # Bounded exactly as the primal's fraction is, so this entry cannot exceed the rate it
    # linearizes where the conduction limit decides both.
    melt_frac_ice = min(
        -UT.guarded_quotient(pp.ice_melting.q_ice, UT.clamp_to_nonneg(q_ice)),
        CMP3.ice_melt_fraction_limit(
            mp.warm_rain.air_properties, tps, mp.ice.scheme, T).inv_τ,
    )  # ≥ 0
    rim_rim -= melt_frac_ice
    # The rime mass drains along the ray, so its self-derivative is exactly `-melt_frac_ice`. The
    # rime volume does not: the melting densification relaxes the rime density toward solid ice, so
    # `∂ₜb_rim = -melt_frac_ice ρ_i b_rim² / q_rim` and its two partials are
    # `∂/∂b_rim = -2 melt_frac_ice ρ_i b_rim / q_rim` and
    # `∂/∂q_rim = +melt_frac_ice ρ_i b_rim² / q_rim²`, both exact at frozen `melt_frac_ice`.
    # Their Euler combination `b_rim ∂/∂b_rim + q_rim ∂/∂q_rim` returns the tendency itself, so the
    # pair is degree-one homogeneous as every other rime term on this row is.
    #
    # `brim_rim` is an off-diagonal rather than a diagonal, so `ExplicitGrowthDiagonal` does not
    # remove it; it is carried in every mode.
    #
    # `q_rim` is guarded exactly as the primal guards it, so the Jacobian cannot carry a slope at a
    # state where the tendency took the ray branch.
    q_rim_pos = FD.value(q_rim) > o
    ρ_i_melt = mp.ice.scheme.ρ_i
    q_rim_safe = max(q_rim, UT.ϵ_numerics(FD.value(q_rim)))
    brim_brim -= ifelse(
        q_rim_pos, 2 * melt_frac_ice * ρ_i_melt * b_rim / q_rim_safe, melt_frac_ice)
    brim_rim += ifelse(
        q_rim_pos, melt_frac_ice * ρ_i_melt * b_rim * b_rim / (q_rim_safe * q_rim_safe), o)

    # ice aggregation (n_ice sink, donor n_ice). The rate is a double integral of
    # `n(D₁)·n(D₂)·σ·|Δv|` and `n(D) ∝ n_ice` exactly at frozen `logλ` - the shape
    # `(μ, λ)`, the kernel and the integration bounds all come from the particle geometry
    # and carry no `n_ice` - so the rate is exactly second-order homogeneous in the donor
    # and `2·rate/n_ice` is its derivative, not an approximation of one. The generic
    # `rate/donor` recipe is the first-order case of the same linearization and would
    # carry exactly half of it here; the recipe is the donor-linearized derivative, so
    # consistency with the file means taking the exponent the rate actually has. Negative,
    # so `ExplicitGrowthDiagonal` keeps it.
    nice_nice += 2 * pp.ice_aggregation.n_ice * dnice

    # ice sublimation: the same ray form, so the same exact self-derivative on both rim rows,
    # with the rate frozen as Tier 3 freezes the melt quadrature. It is a pure sink, so it only
    # adds diagonal damping, and it is the fractional loss `f` itself carries - the pathway test
    # mirrors the ice-number row's so J cannot run a branch f does not.
    sub_frac_rim = ifelse(n_sub_active, sub_frac, o)  # ≤ 0
    rim_rim += sub_frac_rim
    brim_brim += sub_frac_rim

    # ice orphan drain: the same ray form, so the same exact self-derivative on both rim rows
    rim_rim += -inv_τ_orphan_ice_j
    brim_brim += -inv_τ_orphan_ice_j

    # liquid-ice collision (cloud/rain → ice): both sides, so each donor column's
    # condensate rows sum to zero as the primal's do. See [`_riming_jacobian_block`](@ref).
    rb = _riming_jacobian_block(rs, dlcl, drai)
    lcl_lcl += rb.lcl_lcl
    ice_lcl += rb.ice_lcl
    rai_lcl += rb.rai_lcl
    rai_rai += rb.rai_rai
    ice_rai += rb.ice_rai
    rim_lcl += rb.rim_lcl
    rim_rai += rb.rim_rai
    brim_lcl += rb.brim_lcl
    brim_rai += rb.brim_rai

    # Wet-growth densification (BIWET). `_riming_jacobian_block`'s donor-linearization does not
    # cover it: unlike the donor terms above, this is not a transfer FROM cloud/rain TO ice, it
    # is a self-relaxation of the `(L_rim, B_rim)` pair toward its fully-soaked endpoint, so it
    # belongs on the `b_rim` diagonal directly rather than routed through a donor.
    # `BIWET = f_shd·(ρq_ice/ρ_i − B_rim)/τ_wet`
    # ([`CMP3.bulk_liquid_ice_collision_sources`](@ref)) is EXACTLY linear in `B_rim` once
    # `state`, `logλ` and the quadrature rates are held fixed - `f_shd`, `ρq_ice` and `ρ_i` all
    # come from `state`/`rates`, never from `B_rim` itself - so `∂(BIWET)/∂B_rim = −f_shd/τ_wet`
    # is exact, not a donor-recipe approximation, mirroring the ray-form exactness of the
    # sublimation and orphan self-terms above. `rs.f_shd` is exactly zero where the riming block
    # is inactive (the same gate as `_riming_jacobian_block`'s own donor terms), so this needs no
    # separate presence test. Measured, it restores the fraction of the manual-versus-AD
    # `b_rim`/`b_rim` gap attributable to BIWET's own self-term, which ranges from negligible at
    # light riming to about 93 percent at heavy riming; the remainder is
    # `_riming_jacobian_block`'s own geometry-mediated dependence (the collection rates varying
    # with `F_rim`/`ρ_rim` through `state`), a documented Tier-3 limitation, NOT restored here.
    brim_brim -= rs.f_shd / mp.ice.scheme.τ_wet
    # Wet-growth densification of the rime MASS (QIWET), the twin of the BIWET self-term above.
    # `QIWET = f_shd·ρq_ice·(1 - F_rim)/τ_wet` is `f_shd·(L_ice - L_rim)/τ_wet`, so in the
    # per-kg-air units of this Jacobian it is `f_shd·(q_ice - q_rim)/τ_wet` and
    # `∂(QIWET)/∂q_rim = -f_shd/τ_wet`, exact by the same argument the BIWET comment makes.
    #
    # It differs from BIWET in one respect, which is the reason for the gate. BIWET reads the
    # prognostic `B_rim`; QIWET reads the state's `F_rim`, which the `P3State` constructor clamps
    # to `[0, 1 - eps]`. On either saturated branch `F_rim` no longer moves with `q_rim`, so the
    # exact derivative there is zero rather than `-f_shd/τ_wet`. The term is inert on both of those
    # branches anyway, its factor `1 - F_rim` being `eps` at the upper clamp, so the gate costs
    # nothing and keeps `J` from carrying a slope the primal does not have.
    rim_is_interior = (FD.value(q_rim) > o) & (FD.value(q_rim) < FD.value(q_ice))
    rim_rim -= ifelse(rim_is_interior, rs.f_shd / mp.ice.scheme.τ_wet, o)
    # The number sinks keep their own donors. Ice number takes no source from riming -
    # collected drops merge into crystals that already exist - so the number rows need no
    # receiver, and `min(·, o)` guards the rain slot, whose shed source is
    # sign-indeterminate and is not linearized on the rain donor that does not supply it.
    nlcl_nlcl += min(pp.liquid_ice_collision.n_lcl, o) * dnlcl
    nrai_nrai += min(pp.liquid_ice_collision.n_rai, o) * dnrai

    return _jacobian_2mp3(FT;
        lcl_lcl, lcl_nlcl, lcl_rai, lcl_ice,
        nlcl_lcl, nlcl_nlcl, nlcl_rai, nlcl_ice,
        rai_lcl, rai_rai, rai_ice,
        nrai_nlcl, nrai_rai, nrai_nrai, nrai_nice,
        ice_lcl, ice_rai, ice_ice, ice_nice,
        nice_lcl, nice_nlcl, nice_rai, nice_nrai, nice_ice, nice_nice,
        rim_lcl, rim_rai, rim_ice, rim_rim,
        brim_lcl, brim_rai, brim_ice, brim_rim, brim_brim,
    )
end

"""
    _numadj_derivs(FT, q, n, x_min, x_max, τ, sat_excess = 0)

The implicit derivatives `(∂q, ∂n)` of
[`CM2.number_tendency_from_mass_limits`](@ref) `∂ₜn = (n_target − n)/τ` with
`n_target = clamp(n, q/x_max, q/x_min)`: the relaxation diagonal `∂n = −1/τ`
on the clamped and drained arms, zero in the interior and on the retained arm.
The target's mass coupling (`1/(x_min·τ)` on the high clamp, a gain up to
`1/x_min` per unit `τ`) is treated explicitly: `∂q = 0`.

The zero-mass arm is `q ≤ 0` and splits on `sat_excess` exactly as the primal
does: retained where the vapor is in excess, so `n_target = n` and the whole
entry vanishes with the rate, and drained where it is not. Passing the excess
is what keeps that split on one side of the f/J pair from differing from the
other. The function takes no smallness threshold on the mass, deliberately, so
that a caller cannot introduce one on the derivative side of the pair only.
"""
@inline function _numadj_derivs(
    ::Type{FT}, q, n, x_min, x_max, τ, sat_excess = zero(FT),
) where {FT}
    # The zero-mass arm has two branches and they linearize differently: a population retained
    # under supersaturation has `n_target = n`, so the tendency is identically zero there and so
    # is its derivative, while a drained one relaxes at `-1/τ`. Matching the primal's branch here
    # is the f/J consistency requirement, not a refinement - an entry of `-1/τ` against a rate of
    # zero is the strongest wrong entry in this row, since `1/τ = 0.01` sits against `I/h = 0.5`.
    zero_mass = !(q > zero(FT))
    retained = zero_mass & (FD.value(sat_excess) > 0)
    drained = zero_mass & !retained
    lo = q / x_max
    hi = q / x_min
    clamp_low = !zero_mass && n < lo
    clamp_high = !zero_mass && n > hi
    ∂n = ifelse(drained || clamp_low || clamp_high, -1 / τ, zero(FT))
    return (zero(FT), ∂n)
end

#####
##### The temperature-coupled 9×9 Jacobian
#####

"""
    _jacobian_2mp3t_manual(g, y, pp, rs, ctx)

The temperature-coupled 9×9 substep Jacobian, over the state
`(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, T)`.

The species block is [`_jacobian_2mp3_manual`](@ref)'s own 8×8 directly: that function
linearizes the bare phase-change relaxation (`−1/τ`, no `Γ`) at both its condensation and
deposition call sites, which is exactly what the bare primal the temperature-coupled tendency
sums needs, so the 9-vector embeds the 8-vector block unchanged and no folded-to-bare
correction remains. What the 9-vector adds is a temperature column and a temperature row:

  - the column carries the saturation shift of condensation and deposition
    (`−∂q_sat/∂T / τ` on the active, unlimited branches, with the ice-number sublimation
    pathway riding it) and the melting contribution's temperature derivative, taken from
    `ctx` rather than from a linearization, being a genuinely separate quantity the 8×8 has
    no column for;
  - the row is the latent-heating combination of the species rows.

The freezing rates' exponential temperature dependence is not carried; it stays explicit.

`g` is the frozen per-substep context (fields `mp`, `tps`, `ρ`, `q_tot`, `logλ`, `w`, `p`;
`T` is a state slot here, not a context field), `pp` and `rs` are the primal's per-process
breakdown and its liquid-ice collision split by donor, and `ctx` the shared phase-relaxation
context. All three are computed once at this state by the caller and shared with the tendency
rather than recomputed here, which is what makes the 9-vector obey the same one-tendency
property [`_jacobian_2mp3_manual`](@ref) states.
"""
function _jacobian_2mp3t_manual(g, y::SA.StaticVector{9, FT}, pp, rs, ctx) where {FT}
    (q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, T) = y
    o = zero(FT)
    q_tot = FT(g.q_tot)
    x8 = MicroState2MP3(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    # The species context is this context with the state's own temperature: `w` and `p` ride
    # the shared context at both state sizes, so the species block linearizes droplet
    # activation on the branch the primal fires it on.
    g8 = (; g.mp, g.tps, g.ρ, T, q_tot, g.logλ, g.w, g.p)
    J8c = _jacobian_2mp3_manual(g8, x8, pp, rs)
    (; Lᵥ, Lₛ, cp_air, τ_l, dqsl_dT, sat_excess_l,
        τ_i, dqsi_dT, sat_excess_i, ∂melting_∂T, T_freeze) = ctx

    qmin = UT.ϵ_numerics_2M_M(FT)

    # Own gating, matching `_jacobian_2mp3_manual`'s. It is needed for the temperature COLUMN
    # below, which that function has no equivalent of; it is not a correction on `J8c`, which is
    # bare-consistent and gated by construction.
    cond_active = !CM2.cloud_condensation_is_degenerate(τ_l)
    dep_active =
        !(((T > T_freeze) & (sat_excess_i > 0)) | CMP3.ice_deposition_is_degenerate(τ_i))
    ∂ₜq_ice_dep = pp.ice_depsub.q_ice
    n_sub_active = ∂ₜq_ice_dep < 0 && dep_active
    n_per_q = n_ice / max(qmin, q_ice)

    # temperature column: saturation shift on the unlimited phase-change
    # branches, the ice-number sublimation pathway, and the melting rate's
    # linear dependence on the temperature excess
    limit_l = (sat_excess_l < 0) & (-sat_excess_l > max(o, q_lcl))
    limit_i = (sat_excess_i < 0) & (-sat_excess_i > max(o, q_ice))
    t1 = (cond_active && !limit_l) ? -dqsl_dT / τ_l : o
    t5 = (dep_active && !limit_i) ? -dqsi_dT / τ_i : o
    t6 = n_sub_active ? n_per_q * t5 : o
    melt_col = Tuple(∂melting_∂T)
    col9 = SA.SVector{8, FT}(
        t1 + melt_col[1], melt_col[2], melt_col[3], melt_col[4],
        t5 + melt_col[5], t6 + melt_col[6], melt_col[7], melt_col[8])

    lh =
        (
            SA.SVector{8, FT}(ntuple(c -> J8c[1, c] + J8c[3, c], 8)) .* Lᵥ .+
            SA.SVector{8, FT}(ntuple(c -> J8c[5, c], 8)) .* Lₛ
        ) ./ cp_air
    tTT = (Lᵥ * (col9[1] + col9[3]) + Lₛ * col9[5]) / cp_air
    row9 = SA.SVector{9, FT}(Tuple(lh)..., tTT)
    return vcat(hcat(J8c, col9), row9')
end
