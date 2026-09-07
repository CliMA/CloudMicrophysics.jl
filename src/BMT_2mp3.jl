#####
##### 2M + P3 prognostic state
#####

"""
    MicroState2MP3{FT}

The eight prognostic 2M+P3 species as a `StaticArrays.FieldVector`, in the order
`(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)`.

The order is shared by the per-process breakdown
([`p3_2m_process_rates`](@ref)), the substep state vector and the manual
Jacobian, so an index means the same species everywhere.
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

#####
##### Shared saturation and phase-change helpers
#####

"""
    _liquid_sat_excess(tps, ρ, T, q_tot, q_lcl, q_rai, q_ice)

The vapor specific content in excess of saturation over liquid [kg/kg]. Its sign
and magnitude set the orphan-mass drain of [`CM2.orphan_mass_drain`](@ref), and its
sign alone is the deactivation test of the zero-mass droplet arm of
[`CM2.number_tendency_from_mass_limits`](@ref): droplet number carrying no mass
yet is retained while the air is supersaturated and drained when it is not.

One definition, because the primal, the manual Jacobian and the
temperature-coupled context have to take the same branch at the same state, and
they each reach it from different local variables.
"""
@inline _liquid_sat_excess(tps, ρ, T, q_tot, q_lcl, q_rai, q_ice) =
    TDI.q_vap(q_tot, q_lcl + q_rai, q_ice) -
    TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)

"""
    _ice_sat_excess(tps, ρ, T, q_tot, q_lcl, q_rai, q_ice)

The vapor specific content in excess of saturation over ice [kg/kg], mirroring
[`_liquid_sat_excess`](@ref) for the ice orphan-mass drain
([`CM2.orphan_mass_drain_ice`](@ref)). One definition, so the primal and the
manual Jacobian take the same branch at the same state.
"""
@inline _ice_sat_excess(tps, ρ, T, q_tot, q_lcl, q_rai, q_ice) =
    TDI.q_vap(q_tot, q_lcl + q_rai, q_ice) -
    TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)

"""
    _bare_rate_conv_q_vap_to_q_icl(τ, tps, micro, thermo)

The phase-change rate `sat_excess / τ` toward saturation over ice, with the sublimation branch
bounded by the donor content. It also applies
[`CMNonEq.wet_surface_deposition_limiter`](@ref), so the rate is non-positive above the
freezing point. The liquid twin, [`_bare_rate_conv_q_vap_to_q_lcl`](@ref), is the same form and
arrives with the warm-rain surface that first needs it.

`τ` here is the DERIVED diffusional growth timescale,
[`CM2.cloud_condensation_timescale`](@ref) and [`CMP3.ice_deposition_timescale`](@ref), rather
than a phenomenological closure parameter. The cell-scale relaxation `ds/dt = -Γ s / τ` then
comes out of the coupled mass and energy bookkeeping on its own, so folding `Γ` into `τ` at
this call site would apply that correction a second time. The shared
[`CMNonEq.conv_q_vap_to_q_lcl`](@ref) and [`CMNonEq.conv_q_vap_to_q_icl`](@ref) close over a
phenomenological relaxation timescale instead and keep the fold.
"""
@inline function _bare_rate_conv_q_vap_to_q_icl(τ, tps, micro, thermo)
    (; q_icl) = micro
    (; ρ, T) = thermo
    qᵥ = TDI.q_vap(micro.q_tot, micro.q_lcl + micro.q_rai, micro.q_icl + micro.q_sno)
    qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
    sat_excess = qᵥ - qᵥ_sat_ice
    tendency = ifelse(
        sat_excess < 0,
        -min(-sat_excess, max(0, q_icl)) / τ,
        sat_excess / τ,
    )
    limiter = CMNonEq.wet_surface_deposition_limiter(tendency, tps, T)
    return ifelse(limiter, zero(tendency), tendency)
end

"""
    _ice_numadj_params(p3)

The `(τ, x_min, x_max)` bounds for the ice number adjustment, shared by the per-process primal
([`p3_2m_process_rates`](@ref)) and the manual Jacobian, so the rate and its linearization
cannot disagree.

`x_min` is the scheme's own nucleation starter mass
([`CMP3.ice_mean_particle_mass_min`](@ref)) rather than an independent literal, so the bound
cannot sit above the size at which the scheme injects crystals.

`τ` is `p3.τ_numadj`, read from the `P3_ice_number_adjustment_timescale` TOML key. It equals the
warm phase's `Horn2012_number_concentration_adjustment_timescale` value today; the two keys are
independent and may diverge.
"""
@inline function _ice_numadj_params(p3)
    FT = typeof(p3.ρ_i)
    return (;
        τ = p3.τ_numadj,
        x_min = CMP3.ice_mean_particle_mass_min(p3),
        x_max = CMP3.ice_mean_particle_mass_max(FT),
    )
end

"""
    _ice_melting_species(ρ, ρ_i, q_rim, b_rim, dNdt, dLdt, melt_frac)

The eight-species [`MicroState2MP3`](@ref) contribution of ice melting, from the volumetric
number and mass melting rates of [`CMP3.ice_melt`](@ref) and the bounded fractional ice-mass
loss they share ([`CMP3.ice_melt_fraction`](@ref)): ice becomes rain, rime mass drains along
the ray through the origin with that same fraction, and rime volume leaves the ray.

Rime volume leaving the ray is the melting densification. The porous fraction of a melting
particle drains first, so the rime left behind is denser, and the rime density relaxes toward
solid ice at the rate the ice mass melts. In these variables that is
`∂ₜb_rim = -melt_frac ρ_i b_rim² / q_rim`, which collapses to the ray drain exactly at
`ρ_rim = ρ_i`, where already-solid rime only drains. `q_rim` is guarded rather than divided
into: with no rime mass there is no density to carry, and the ray drain is the only value the
quotient could take by continuity. The guard is a division backstop rather than a physical
scale: the rime pair is projected onto the density cone, so `b_rim² / q_rim` vanishes with
`q_rim` and the quotient needs no floor of its own.

Linear in `(dNdt, dLdt, melt_frac)`, so applying it to the temperature derivatives of the
three rates returns the temperature derivative of the contribution.
"""
@inline function _ice_melting_species(ρ, ρ_i, q_rim, b_rim, dNdt, dLdt, melt_frac)
    ∂ₜq_ice_melt = dLdt / ρ
    ∂ₜn_ice_melt = dNdt / ρ
    o = zero(∂ₜq_ice_melt)
    ∂ₜq_rim_melt = -q_rim * melt_frac
    ∂ₜb_rim_melt =
        -ifelse(
            FD.value(q_rim) > zero(FD.value(q_rim)),
            melt_frac * ρ_i * b_rim * b_rim / max(q_rim, UT.ϵ_numerics(FD.value(q_rim))),
            b_rim * melt_frac,
        )
    return MicroState2MP3(
        o, o, ∂ₜq_ice_melt, ∂ₜn_ice_melt, -∂ₜq_ice_melt, -∂ₜn_ice_melt,
        ∂ₜq_rim_melt, ∂ₜb_rim_melt,
    )
end

#####
##### The per-process primal
#####

"""
    p3_2m_process_rates(mp::CMP.Microphysics2MParams{WR, ICE}, tps, micro, thermo)

The 2M warm rain + P3 ice tendency decomposed into its per-process contributions, returned as
`(pp, rs)`.

`micro` is a `NamedTuple` of specific contents,
`(; q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)`, and `thermo` is a
`NamedTuple` of thermodynamic state, `(; ρ, T, w, p, logλ)`. The vertical velocity `w` and the
air pressure `p` reach droplet activation through `thermo` and nothing else in the tendency
reads them; at zero only the adiabatic-parcel branch of the activation supersaturation is
switched off, and the ambient branch still fires wherever the air is supersaturated over
liquid.

`pp` is a `NamedTuple` of 21 slots, each a [`MicroState2MP3`](@ref) over the eight prognostic
species:

    activation, cloud_condevap, rain_evap, autoconv, cloud_selfcol,
    accretion, rain_selfcol, rain_breakup, cloud_numadj, rain_numadj,
    cloud_orphan, rain_orphan, ice_orphan,
    liquid_ice_collision, ice_aggregation, ice_melting, ice_deposition,
    immersion_freezing, ice_depsub, ice_numadj, rain_freezing

This is the ONLY evaluator of the 2M+P3 rates. The net tendency is the component-wise sum
`sum(values(pp))` at every call site, so the parts-sum-to-total identity

    sum(values(first(p3_2m_process_rates(mp, tps, micro, thermo)))) ==
        Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ, w, p)(x)

holds by construction rather than by two implementations being kept in step, and a test asserts
the slot names and their order.

The primal owns the ancillary treatment once: the canonicalization prelude (non-negative clamps
on the masses and numbers, a positive floor on the air density), the ice presence gate, the
mean-mass number bounds, the orphan drains, and the number-adjustment parameters.

`rs` is the liquid-ice collision transfer split by the donor that supplies it,
`(; ∂ₜq_lcl_frz, ∂ₜq_lcl_shd, ∂ₜq_rai_frz, ∂ₜb_rim_lcl, ∂ₜb_rim_rai, f_shd)`, in the per-kg-air
units of the breakdown (`f_shd` dimensionless), all zero where the transfer is inactive. `pp`
cannot carry the split: it aggregates the whole transfer into one species vector per process,
whose cloud, rain and ice mass slots sum to zero, so two of them determine the third. A
linearization that routes the transfer back to its donor therefore needs
[`CMP3.bulk_liquid_ice_collision_sources`](@ref)'s per-donor rates, and this is where the
quadrature is already paid for. `f_shd` carries the bulk shed fraction, so the manual
Jacobian can differentiate the wet-growth self-term `-B_rim/τ_wet` without recomputing the
quadrature.
"""
@inline function p3_2m_process_rates(
    mp::CMP.Microphysics2MParams{WR, ICE}, tps, micro, thermo,
) where {WR, ICE <: CMP.P3IceParams}
    FT = eltype(thermo.ρ)
    (; T, w, p, logλ) = thermo
    # The AIR DENSITY takes a POSITIVE floor where the masses and numbers take a non-negative
    # clamp: it appears in denominators, inside logs and under negative fractional powers, so
    # zero merely trades a throw for a NaN. See [`AIR_DENSITY_FLOOR`](@ref).
    ρ = _floored_air_density(thermo.ρ)
    q_tot = UT.clamp_to_nonneg(micro.q_tot)
    q_lcl = UT.clamp_to_nonneg(micro.q_lcl)
    q_rai = UT.clamp_to_nonneg(micro.q_rai)
    n_lcl = UT.clamp_to_nonneg(micro.n_lcl)
    n_rai = UT.clamp_to_nonneg(micro.n_rai)
    q_ice = UT.clamp_to_nonneg(micro.q_ice)
    n_ice = UT.clamp_to_nonneg(micro.n_ice)
    q_rim = UT.clamp_to_nonneg(micro.q_rim)
    b_rim = UT.clamp_to_nonneg(micro.b_rim)

    o = zero(FT)
    # builder for a per-process vector in the species order of [`MicroState2MP3`](@ref)
    Z() = MicroState2MP3(o, o, o, o, o, o, o, o)

    # The species and thermodynamic views the vapor-exchange helpers read, built once from the
    # canonicalized state. P3 carries a single ice category, so the snow slot is empty.
    micro_vap = (; q_tot, q_lcl, q_icl = q_ice, q_rai, q_sno = zero(q_ice))
    thermo_vap = (; ρ, T)

    # Process rates are evaluated at the mean-mass-bounded populations; the
    # number adjustments relax the prognostic numbers toward the same bounds.
    sb = mp.warm_rain.seifert_beheng
    # the zero-mass droplet arm's deactivation test, taken once and used by the bounded
    # population the rates see and by the number adjustment that drains it, so the two cannot
    # disagree
    sat_excess_l = _liquid_sat_excess(tps, ρ, T, q_tot, q_lcl, q_rai, q_ice)
    # `false` is the warm phase refusing to invent a population where there is none: at the
    # phantom corner the rates below read the cell as empty, and its orphan mass is removed as
    # mass by the drain rather than given a manufactured carrier. The ice number adjustment
    # below follows the same doctrine.
    n_lcl_b = CM2.number_bounded_by_mass_limits(
        (; x_min = sb.pdf_c.xc_min, x_max = sb.pdf_c.xc_max), q_lcl, n_lcl, sat_excess_l;
        invent_from_zero = false)
    n_rai_b = CM2.number_bounded_by_mass_limits(
        (; x_min = sb.pdf_r.xr_min, x_max = sb.pdf_r.xr_max), q_rai, n_rai; invent_from_zero = false)

    # Volumetric quantities for the P3 functions.
    L_lcl = q_lcl * ρ
    L_rai = q_rai * ρ
    N_lcl = n_lcl_b * ρ
    N_rai = n_rai_b * ρ
    L_ice = q_ice * ρ
    N_ice = n_ice * ρ
    L_rim = q_rim * ρ
    B_rim = b_rim * ρ
    state = CMP3.state_from_prognostic(mp.ice.scheme, L_ice, N_ice, L_rim, B_rim)

    aps = mp.warm_rain.air_properties

    #####
    ##### Warm-rain processes
    #####
    warm_rain = mp.warm_rain

    # Droplet activation: cloud number AND the mass those droplets carry.
    #
    # A nucleation-class source, so it obeys the same rule the ice-side deposition nucleation
    # source obeys: both moments together, at the per-particle mass of the new particles. Each
    # droplet arrives at `CM2.activation_droplet_mass`, which IS the size distribution's minimum
    # droplet mass by derivation rather than by coincidence - the smallest resolvable droplet is
    # a freshly activated one - and which is exactly a 1 μm droplet and exactly the activation
    # radius the rate's diffusional timescale is built on. Supplying number alone leaves the
    # category in a state the scheme has no size for, and the number adjustment below then
    # correctly drains it, so the population never establishes itself.
    #
    # Inside the substep this co-evolves with condensation in one implicit solve, so the
    # supersaturation it activates at is the one the newly grown droplets have already drawn
    # down.
    #
    # `w` and `p` reach the parcel branch of the activation supersaturation only. The ambient
    # branch does not depend on them, so this fires at their defaults wherever the air is
    # supersaturated over liquid, and it is inert only where it is not or where the prescribed
    # aerosol has no particles.
    act = CMAA.cloud_droplet_activation_rate(
        warm_rain.activation, warm_rain.aerosol, aps, tps,
        T, p, w, ρ, q_tot, q_lcl + q_rai, q_ice, n_lcl,
        CM2.activation_droplet_mass(sb.pdf_c),
    )
    activation = MicroState2MP3(act.∂ₜq_lcl, act.∂ₜn_lcl, o, o, o, o, o, o)

    # cloud condensation / evaporation (cloud mass only; number neglected).
    #
    # The relaxation timescale follows the droplet population's capacitance integral, so the
    # unbounded timescale diverges as the population vanishes. But
    # `cloud_condensation_timescale` caps it at `CLOUD_COND_TIMESCALE_MAX` to stay finite, and
    # `sat_excess / τ_max` is not zero, so an existence threshold IS required. The cap binds
    # exactly, and it binds on a switch of its own: `log_pdf_cloud_parameters_mass` returns
    # `logA = -Inf` when the droplet NUMBER is below presence, which sends the diameter moment
    # to zero and the quotient to the cap. A population with number but no mass is not
    # degenerate: its mean droplet mass is floored at the activation droplet mass, it has
    # surface area, and it condenses. Gating on the degeneracy of τ itself needs no new tuned
    # parameter and fires exactly when the capacitance integral underflowed. The number slot is
    # zero either way, so the ungated rate would be a pure mass source with no droplets to carry
    # it.
    τ_cond = CM2.cloud_condensation_timescale(sb.pdf_c, aps, tps, T, ρ, q_lcl, N_lcl)
    ∂ₜq_lcl_cond = _bare_rate_conv_q_vap_to_q_lcl(τ_cond, tps, micro_vap, thermo_vap)
    # No droplets, no surfaces: condensation AND evaporation are both exactly zero.
    ∂ₜq_lcl_cond = ifelse(
        CM2.cloud_condensation_is_degenerate(τ_cond), zero(∂ₜq_lcl_cond), ∂ₜq_lcl_cond)
    cloud_condevap = MicroState2MP3(∂ₜq_lcl_cond, o, o, o, o, o, o, o)

    # rain evaporation (rain mass + number)
    evap = CM2.rain_evaporation(sb, aps, tps, q_tot, q_lcl, q_ice, q_rai, zero(q_ice), ρ, N_rai, T)
    rain_evap = MicroState2MP3(o, o, evap.∂ₜq_rai, evap.∂ₜρn_rai / ρ, o, o, o, o)

    # autoconversion (cloud → rain, mass + number)
    acnv = CM2.autoconversion(sb.acnv, sb.pdf_c, q_lcl, q_rai, ρ, N_lcl)
    autoconv = MicroState2MP3(
        acnv.dq_lcl_dt, acnv.dN_lcl_dt / ρ, acnv.dq_rai_dt, acnv.dN_rai_dt / ρ, o, o, o, o,
    )

    # cloud self-collection (cloud number only), evaluated at the true N_lcl
    ∂ₜN_lcl_sc = CM2.cloud_liquid_self_collection(sb.acnv, sb.pdf_c, q_lcl, ρ, ρ * n_lcl, acnv.dN_lcl_dt)
    cloud_selfcol = MicroState2MP3(o, ∂ₜN_lcl_sc / ρ, o, o, o, o, o, o)

    # accretion (cloud → rain, mass; cloud number)
    accr = CM2.accretion(sb, q_lcl, q_rai, ρ, N_lcl)
    accretion_wr = MicroState2MP3(accr.dq_lcl_dt, accr.dN_lcl_dt / ρ, accr.dq_rai_dt, o, o, o, o, o)

    # rain self-collection (rain number only)
    ∂ₜN_rai_sc = CM2.rain_self_collection(sb.pdf_r, sb.self, q_rai, ρ, N_rai)
    rain_selfcol = MicroState2MP3(o, o, o, ∂ₜN_rai_sc / ρ, o, o, o, o)

    # rain breakup (rain number only)
    ∂ₜN_rai_br = CM2.rain_breakup(sb.pdf_r, sb.brek, q_rai, ρ, N_rai, ∂ₜN_rai_sc)
    rain_breakup = MicroState2MP3(o, o, o, ∂ₜN_rai_br / ρ, o, o, o, o)

    # number adjustment for mass limits (cloud, then rain)
    numadj_lcl = (; sb.numadj.τ, x_min = sb.pdf_c.xc_min, x_max = sb.pdf_c.xc_max)
    ∂ₜn_lcl_numadj = CM2.number_tendency_from_mass_limits(
        numadj_lcl, q_lcl, n_lcl, sat_excess_l; invent_from_zero = false)
    cloud_numadj = MicroState2MP3(o, ∂ₜn_lcl_numadj, o, o, o, o, o, o)
    numadj_rai = (; sb.numadj.τ, x_min = sb.pdf_r.xr_min, x_max = sb.pdf_r.xr_max)
    ∂ₜn_rai_numadj = CM2.number_tendency_from_mass_limits(numadj_rai, q_rai, n_rai; invent_from_zero = false)
    rain_numadj = MicroState2MP3(o, o, o, ∂ₜn_rai_numadj, o, o, o, o)

    # Orphan mass: condensate whose number is absent, so no particle carries it. The number
    # adjustment above does not invent one, which leaves the mass to be dealt with as mass.
    # Below saturation it evaporates at the rate a population of minimum-mass particles would;
    # at or above saturation it is left alone, because activation supplies a real number and the
    # population that arrives adopts it. See `CM2.orphan_mass_drain`.
    #
    # BOTH halves of each predicate read `FD.value`, and the mass half is not optional even
    # though it reads like one. A `Dual` with value zero and a nonzero partial compares GREATER
    # than zero, so a bare `q > o` reads an empty species as PRESENT under AD and selects the
    # drain at a state that has no orphan. The VALUE is still right, because the rate is
    # proportional to the mass and the mass is zero; only the DERIVATIVE leaks, which is why
    # this is invisible to any check on the tendency itself.
    lcl_is_orphan = !(FD.value(n_lcl) > o) & (FD.value(q_lcl) > o)
    rai_is_orphan = !(FD.value(n_rai) > o) & (FD.value(q_rai) > o)
    ∂ₜq_lcl_orphan = ifelse(
        lcl_is_orphan,
        CM2.orphan_mass_drain(aps, tps, T, ρ, q_lcl, sb.pdf_c.xc_min, sb.pdf_c.ρw,
            sat_excess_l),
        o,
    )
    ∂ₜq_rai_orphan = ifelse(
        rai_is_orphan,
        CM2.orphan_mass_drain(aps, tps, T, ρ, q_rai, sb.pdf_r.xr_min, sb.pdf_r.ρw,
            sat_excess_l),
        o,
    )
    cloud_orphan = MicroState2MP3(∂ₜq_lcl_orphan, o, o, o, o, o, o, o)
    rain_orphan = MicroState2MP3(o, o, ∂ₜq_rai_orphan, o, o, o, o, o)

    #####
    ##### P3 ice processes
    #####
    p3 = mp.ice.scheme
    vel = mp.ice.terminal_velocity
    pdf_c = mp.ice.cloud_pdf
    pdf_r = mp.ice.rain_pdf
    ice_nucleation = mp.ice.ice_nucleation
    quad = mp.ice.quad

    # liquid-ice collision, ice aggregation, ice melting.
    #
    # Gated on the PRESENCE of both ice moments rather than on a mass epsilon. An `eps(Float32)`
    # gate is 1.1920929e-7 kg/kg, so below it melting is identically zero at every temperature
    # and switches off discontinuously mid-melt, leaving trace ice unremovable in warm air;
    # riming and aggregation go with it, and none of it applies at Float64. See
    # [`CMP3.ice_population_is_present`](@ref) for why the test is on both moments.
    if CMP3.ice_population_is_present(state)
        coll = CMP3.bulk_liquid_ice_collision_sources(
            state, logλ, pdf_c, pdf_r, L_lcl, N_lcl, L_rai, N_rai, aps, tps, vel, ρ, T;
            B_rim, quad, assembly = CMP3._liqice_partition(mp.ice.liqice_partition, quad),
        )
        liquid_ice_collision = MicroState2MP3(
            coll.∂ₜq_c, coll.∂ₜN_c / ρ, coll.∂ₜq_r, coll.∂ₜN_r / ρ,
            coll.∂ₜL_ice / ρ, o, coll.∂ₜL_rim / ρ, coll.∂ₜB_rim / ρ,
        )
        riming_split = (;
            ∂ₜq_lcl_frz = coll.∂ₜq_c_frz, ∂ₜq_lcl_shd = coll.∂ₜq_c_shd,
            ∂ₜq_rai_frz = coll.∂ₜq_r_frz,
            ∂ₜb_rim_lcl = coll.∂ₜB_rim_c / ρ, ∂ₜb_rim_rai = coll.∂ₜB_rim_r / ρ,
            f_shd = coll.f_shd,
        )

        # The sticking efficiency is an exact per-cell scalar on the rate, applied HERE rather
        # than inside the kernel: it comes out of the collision integral exactly, its kinks are
        # evaluated exactly rather than interpolated, and applying it above the quadrature/table
        # dispatch is what keeps those two paths from computing different physics.
        S_ice_agg = CMP3.ice_self_collection(state, logλ, vel, ρ; quad)
        E_stick = CMP3.ice_sticking_efficiency(state.params, T, state.F_rim)
        ice_aggregation = MicroState2MP3(o, o, o, o, o, -E_stick * S_ice_agg.dNdt / ρ, o, o)

        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        melt = ifelse(T > T_freeze,
            CMP3.ice_melt(vel, aps, tps, T, ρ, state, logλ; quad),
            CMP3.zero_ice_melt(ρ),
        )
        ice_melting =
            _ice_melting_species(ρ, mp.ice.scheme.ρ_i, q_rim, b_rim, melt.dNdt, melt.dLdt, melt.melt_frac)
    else
        liquid_ice_collision = Z()
        riming_split = (;
            ∂ₜq_lcl_frz = o, ∂ₜq_lcl_shd = o, ∂ₜq_rai_frz = o,
            ∂ₜb_rim_lcl = o, ∂ₜb_rim_rai = o, f_shd = o,
        )
        ice_aggregation = Z()
        ice_melting = Z()
    end

    # Deposition nucleation of pristine ice (F_rim = 0), a nucleation-class source that supplies
    # both moments at the nascent-crystal mass. The depletion proxy, the seed delivery time and
    # the liquid grouping for the saturation ratio are all owned by the rate function itself; see
    # `HetIceNucleation.deposition_rate`.
    micro_dep = (; q_tot, q_lcl, q_rai, q_ice, n_ice)
    dep = CM_HetIce.deposition_rate(ice_nucleation, mp, tps, micro_dep, thermo_vap)
    ice_deposition = MicroState2MP3(o, o, o, o, dep.∂ₜq_frz, dep.∂ₜn_frz, o, o)

    # Freezing of cloud drops (fully-rimed embryo graupel), on the SAME composition rain gets.
    #
    # Liquid is liquid: cloud drops and rain drops are not different substances, they are
    # different sizes, so `cloud_freezing_rate` is `rain_freezing_rate`'s composition with the
    # cloud PSD and the Stokes fall speed substituted. Two parallel nucleation pathways whose
    # coefficients add - Bigg immersion and Koop homogeneous - then the same three sequential
    # kinetic stages, with no ice-nucleating-particle budget above either: the reference P3
    # implementation runs cloud-droplet immersion as Bigg alone, with no bound and no aerosol
    # dependence; see `HetIceNucleation.cloud_freezing_rate`.
    cld_frz = CM_HetIce.cloud_freezing_rate(
        mp.ice.rain_freezing, mp.ice.homogeneous, p3.vent, aps, tps,
        pdf_c, q_lcl, ρ, N_lcl, T, TDI.q_vap(q_tot, q_lcl + q_rai, q_ice),
    )
    ∂ₜn_imm = cld_frz.∂ₜn_frz
    ∂ₜq_imm = cld_frz.∂ₜq_frz
    immersion_freezing = MicroState2MP3(
        -∂ₜq_imm, -∂ₜn_imm, o, o, ∂ₜq_imm, ∂ₜn_imm, ∂ₜq_imm, ∂ₜq_imm / p3.ρ_i,
    )

    # ice deposition / sublimation (rime drains on the sublimation branch only).
    #
    # The relaxation timescale follows the population's capacitance integral, so the unbounded
    # timescale diverges as the population vanishes. But `ice_deposition_timescale` caps it at
    # `ICE_DEP_TIMESCALE_MAX` to stay finite, and `deficit / τ_max` is not zero. At the empty
    # ice state the rate is exactly `(qᵥ - qᵥ_sat_ice) / (ICE_DEP_TIMESCALE_MAX ⋅ Γᵢ)`, of order
    # 1e-13 kg/kg/s, while a populated state gives a timescale four orders below the cap. So the
    # bound binds only where the population is absent, and there it is a mass source with no
    # accompanying number source, which manufactures exactly the mass-without-number states the
    # scheme has no size for. It is slow, but it applies to every ice-free supersaturated cell
    # continuously.
    #
    # An existence threshold is therefore required, not as a workaround for a broken timescale
    # but because the safety bound makes the "rate vanishes with the population" argument false
    # at precisely the state it needed to hold for. Gating on the degeneracy of τ itself needs
    # no new tuned parameter and self-scales: it fires exactly when the capacitance integral
    # underflowed.
    τ_dep = CMP3.ice_deposition_timescale(vel, aps, tps, T, ρ, state, logλ; quad)
    ∂ₜq_ice_dep = _bare_rate_conv_q_vap_to_q_icl(τ_dep, tps, micro_vap, thermo_vap)
    # No particles, no surfaces: deposition AND sublimation are both exactly zero. Zeroing the
    # whole rate rather than the growth branch alone is deliberate - the sublimation branch is
    # driven by the same capacitance integral, so it is equally unsupported here, and the rime
    # drains below are formed from it and would otherwise drain rime from a state with none.
    ∂ₜq_ice_dep =
        ifelse(CMP3.ice_deposition_is_degenerate(τ_dep), zero(∂ₜq_ice_dep), ∂ₜq_ice_dep)
    # Sublimation removes crystals in proportion to the mass it removes, at every positive ice
    # loading. The quotient's only hazard is 0/0, so it is guarded at zero and not at a mass
    # epsilon: the fractional loss is formed FIRST, and the sublimation branch bounds
    # `|∂ₜq_ice_dep|` by `q_ice / (τ_dep ⋅ Γᵢ)`, so the quotient is finite for any positive
    # `q_ice` and the number rate is bounded by `n_ice / (τ_dep ⋅ Γᵢ)`. Guarding at `eps(FT)`
    # instead returns exactly zero, so ice below that threshold sublimates its mass away and
    # leaves every crystal behind, at a physical loading and only at Float32.
    sub_frac_ice = UT.guarded_quotient(∂ₜq_ice_dep, q_ice)
    ∂ₜn_ice_dep = ifelse(∂ₜq_ice_dep < 0, n_ice * sub_frac_ice, zero(∂ₜq_ice_dep))
    # Sublimation removes rime along the same ray through the origin that melting does, reusing
    # the fractional loss formed just above. The `min` picks out the sublimation branch, so
    # depositional growth adds exactly zero rime and dilutes `F_rim` instead.
    sub_frac_rim = min(sub_frac_ice, zero(FT))
    ∂ₜq_rim_sub = q_rim * sub_frac_rim
    ∂ₜb_rim_sub = b_rim * sub_frac_rim
    ice_depsub = MicroState2MP3(o, o, o, o, ∂ₜq_ice_dep, ∂ₜn_ice_dep, ∂ₜq_rim_sub, ∂ₜb_rim_sub)

    # ice number adjustment for mass limits; `false` as for the warm phase, so mass with no
    # number is not given a manufactured carrier
    numadj = _ice_numadj_params(p3)
    ∂ₜn_ice_numadj =
        CM2.number_tendency_from_mass_limits(numadj, q_ice, n_ice; invent_from_zero = false)
    ice_numadj = MicroState2MP3(o, o, o, o, o, ∂ₜn_ice_numadj, o, o)

    # Orphan ice: mass whose number is absent, removed as mass since the adjustment above does
    # not invent a carrier for it - the warm-phase doctrine, one phase over. Below ice
    # saturation it sublimates at the rate a population of nucleation-mass crystals would. At
    # or above ice saturation the drain is zero and the mass is left in place; deposition
    # nucleation supplies a number only below its own temperature threshold, so warmer
    # supersaturated orphan mass waits on transport or on the air drying. The rime pair drains
    # along the ray through the origin with the same fraction, so the drain cannot leave rime
    # behind ice. See `CM2.orphan_mass_drain_ice`.
    sat_excess_i_orphan = _ice_sat_excess(tps, ρ, T, q_tot, q_lcl, q_rai, q_ice)
    ice_is_orphan = !(FD.value(n_ice) > o) & (FD.value(q_ice) > o)
    inv_τ_orphan_ice = ifelse(
        ice_is_orphan,
        CM2.orphan_mass_inv_timescale_ice(
            aps, tps, T, ρ, CMP3.ice_mean_particle_mass_min(p3), p3.ρ_i,
            sat_excess_i_orphan),
        o,
    )
    ice_orphan = MicroState2MP3(
        o, o, o, o, -q_ice * inv_τ_orphan_ice, o,
        -q_rim * inv_τ_orphan_ice, -b_rim * inv_τ_orphan_ice,
    )

    # rain heterogeneous freezing (frozen rain fully rimed), limited by how fast the latent heat
    # of fusion can leave the drop.
    #
    # `liquid_freezing_rate` returns the NUCLEATION rate: the Bigg volumetric rate
    # `J = B exp(a (T_freeze - T))`, fitted over roughly 0 to 40 K of supercooling, integrated
    # over the raindrop PSD and gated only on presence. Nothing in it says a drop that has
    # nucleated still has to shed `Lf` before it is ice, and at 82 K of supercooling the
    # nucleation rate alone reaches 1e18 kg⁻¹ s⁻¹ of number.
    #
    # Two PARALLEL nucleation pathways whose coefficients add - Bigg immersion and Koop
    # homogeneous - followed by three SEQUENTIAL conversion stages whose timescales add:
    # nucleate, shed `Lf`, then let the dendrites cross the drop. Adding Koop is what stops the
    # scheme relying on Bigg's out-of-range extrapolation being accidentally right; the
    # physically supported pathway takes over at `ΔT ≈ 35` K, before the extrapolated one runs
    # away. `τ_heat` binds from roughly 20 K of supercooling for millimeter drops and vanishes
    # past `ΔT* ≈ 53` K, where a drop carries enough of its own cold to absorb all of `Lf`;
    # `τ_dend` is what remains there, and it cuts the 82 K rate by about 13 orders. The result
    # is still far above the donor per step, and correctly so: the boundedness is supplied by
    # the substep's implicit update, not by the rate.
    #
    # `τ_heat` is the full Musil balance, conduction plus evaporation at the melting-point drop
    # surface, so it needs the ambient humidity. `TDI.q_vap` clamps a negative remainder to
    # zero, and `drop_freezing_heat_timescale` clamps the resulting vapor deficit at zero as
    # well, so a corrupt `q_tot` can only remove the evaporative term, never turn it into a heat
    # source.
    rain_frz = CM_HetIce.rain_freezing_rate(
        mp.ice.rain_freezing, mp.ice.homogeneous, p3.vent, aps, tps, sb.evap,
        pdf_r, q_rai, ρ, N_rai, T, TDI.q_vap(q_tot, q_lcl + q_rai, q_ice),
    )
    rain_freezing = MicroState2MP3(
        o, o, -rain_frz.∂ₜq_frz, -rain_frz.∂ₜn_frz,
        rain_frz.∂ₜq_frz, rain_frz.∂ₜn_frz, rain_frz.∂ₜq_frz, rain_frz.∂ₜq_frz / p3.ρ_i,
    )

    pp = (;
        activation, cloud_condevap, rain_evap, autoconv, cloud_selfcol,
        accretion = accretion_wr, rain_selfcol, rain_breakup, cloud_numadj, rain_numadj,
        cloud_orphan, rain_orphan, ice_orphan,
        liquid_ice_collision, ice_aggregation, ice_melting, ice_deposition,
        immersion_freezing, ice_depsub, ice_numadj, rain_freezing,
    )
    return (pp, riming_split)
end

#####
##### Tendency functors: the frozen substep context applied to a species vector
#####

"""
    _2mp3_context(g, x)

The `(micro, thermo)` pair [`p3_2m_process_rates`](@ref) takes, built from a frozen substep
context `g` and the species vector `x`. `q_tot` is promoted to the state's element type; the
frozen entries `ρ`, `T`, `w`, `p` and `logλ` stay in their own type, so a differentiated `x`
does not seed them.
"""
@inline function _2mp3_context(g, x::SA.StaticVector{8, FT}) where {FT}
    (q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim) = x
    micro = (; q_tot = FT(g.q_tot), q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    thermo = (; g.ρ, g.T, g.w, g.p, g.logλ)
    return micro, thermo
end

"""
    Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ, w = zero(ρ), p = zero(ρ))

Callable bundling the frozen per-substep context; applying it to the species vector returns the
net 2M+P3 tendency as a [`MicroState2MP3`](@ref), the component-wise sum
`sum(values(pp))` of [`p3_2m_process_rates`](@ref)'s 21 slots.

`w` and `p` are frozen across the substep exactly as `T` and `ρ` are: droplet activation reads
them and the state does not evolve them. Their defaults of zero switch off only the
adiabatic-parcel branch of the activation supersaturation. The ambient branch does not depend on
the updraft, so a caller that supplies no host vertical velocity still activates in
supersaturated air.
"""
struct Instantaneous2MP3Tendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    T::F
    q_tot::F
    logλ::F
    w::F
    p::F
end
Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ) =
    Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ, zero(ρ), zero(ρ))
@inline function (g::Instantaneous2MP3Tendency)(x::SA.StaticVector{8})
    micro, thermo = _2mp3_context(g, x)
    pp, _ = p3_2m_process_rates(g.mp, g.tps, micro, thermo)
    return sum(values(pp))
end

"""
    Verbose2MP3Tendency(mp, tps, ρ, T, q_tot, logλ, w = zero(ρ), p = zero(ρ))

Per-process companion to [`Instantaneous2MP3Tendency`](@ref): applying it to the species vector
returns [`p3_2m_process_rates`](@ref)'s `NamedTuple` of per-process contributions, each a
[`MicroState2MP3`](@ref), instead of only their sum. It supplies the per-process right-hand
sides `f_p` of the linear post-solve attribution.
"""
struct Verbose2MP3Tendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    T::F
    q_tot::F
    logλ::F
    w::F
    p::F
end
Verbose2MP3Tendency(mp, tps, ρ, T, q_tot, logλ) =
    Verbose2MP3Tendency(mp, tps, ρ, T, q_tot, logλ, zero(ρ), zero(ρ))
@inline function (g::Verbose2MP3Tendency)(x::SA.StaticVector{8})
    micro, thermo = _2mp3_context(g, x)
    return first(p3_2m_process_rates(g.mp, g.tps, micro, thermo))
end

#####
##### The instantaneous 2M + P3 entry
#####

"""
    bulk_microphysics_tendencies(
        ::Microphysics2Moment,
        mp::Microphysics2MParams{WR, <:P3IceParams}, tps,
        ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai,
        q_ice, n_ice, q_rim, b_rim, logλ, w = zero(ρ), p = zero(ρ),
    )

Compute 2-moment **warm rain + P3 ice** microphysics tendencies.

The tendency is the component-wise sum of [`p3_2m_process_rates`](@ref)'s 21 per-process slots,
which is the only evaluator of these rates. This method is type-stable and GPU-optimized; the
P3 ice parameters are guaranteed to be non-`Nothing`, eliminating runtime type checks and
dynamic dispatch.

# Arguments
- `mp`: `Microphysics2MParams` with P3 ice parameters present
- `tps`: Thermodynamics parameters
- `ρ`: Air density (kg/m³)
- `T`: Temperature (K)
- `q_tot`: Total water specific content (kg/kg)
- `q_lcl`: Cloud liquid specific content (kg/kg)
- `n_lcl`: Cloud droplet number per kg air (1/kg)
- `q_rai`: Rain specific content (kg/kg)
- `n_rai`: Rain number per kg air (1/kg)
- `q_ice`: Ice specific content (kg/kg)
- `n_ice`: Ice number per kg air (1/kg)
- `q_rim`: Rime mass (kg/kg)
- `b_rim`: Rime volume (m³/kg)
- `logλ`: Log of the P3 distribution slope parameter, log(1/m)
- `w`: Host vertical velocity (m/s), read by droplet activation only
- `p`: Air pressure (Pa), read by droplet activation only

# Returns
`NamedTuple` with all tendency fields:
- `dq_lcl_dt`: Cloud liquid tendency (kg/kg/s)
- `dn_lcl_dt`: Cloud number tendency (1/kg/s)
- `dq_rai_dt`: Rain tendency (kg/kg/s)
- `dn_rai_dt`: Rain number tendency (1/kg/s)
- `dq_ice_dt`: Ice tendency (kg/kg/s)
- `dn_ice_dt`: Ice number tendency (1/kg/s)
- `dq_rim_dt`: Rime mass tendency (kg/kg/s)
- `db_rim_dt`: Rime volume tendency (m³/kg/s)
- `dn_lcl_activation_dt`: the droplet-activation part of `dn_lcl_dt`, which the host reports
  separately and which is already included in `dn_lcl_dt`
"""
@inline function bulk_microphysics_tendencies(
    ::Microphysics2Moment, mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai,
    q_ice, n_ice, q_rim, b_rim, logλ,
    w = zero(ρ), p = zero(ρ),
) where {WR, ICE <: CMP.P3IceParams}
    micro = (; q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    thermo = (; ρ, T, w, p, logλ)
    pp, _ = p3_2m_process_rates(mp, tps, micro, thermo)
    f = sum(values(pp))
    return (;
        dq_lcl_dt = f.q_lcl, dn_lcl_dt = f.n_lcl,
        dq_rai_dt = f.q_rai, dn_rai_dt = f.n_rai,
        dq_ice_dt = f.q_ice, dn_ice_dt = f.n_ice,
        dq_rim_dt = f.q_rim, db_rim_dt = f.b_rim,
        dn_lcl_activation_dt = pp.activation.n_lcl,
    )
end

#####
##### Substep hooks carried by the 2M + P3 state
#####

"""
    _rosenbrock_species_mask(x::MicroState2MP3)

Diagonal of the species projection matrix `P` used by the linearized-implicit substep: 1 for
active species, 0 for near-empty ones (condensed mass below `1e-10`, for rain and for
ice+rime). A masked species takes the forward-Euler update while active species stay implicit.

The LIQUID species is never masked. The mask reads a small condensed mass as "no stiff dynamics
here", which holds for growth and decay of an existing population and fails for a NUCLEATION
source relaxing INTO an empty species: droplet activation carries `1/τ_act` of order 100 per
second at a strongly supersaturated state, so a masked liquid species takes an explicit step
with `h/τ_act` in the hundreds and overshoots the relaxation target by that factor. On the
near-empty band `q_lcl = 1e-13`, `n_lcl = 100` at 92 percent supersaturation over liquid the
explicit path returns 4.2e9 droplets per kg against a relaxation target of 1e7; the implicit
path with the `-1/τ_act` diagonal returns the target.

This method is reached only by the [`ExactJacobian`](@ref) + [`ImplicitGrowth`](@ref)
combination; both named presets use [`_full_species_mask`](@ref).
"""
@inline function _rosenbrock_species_mask(x::MicroState2MP3{FT}) where {FT}
    ϵ_empty = FT(1e-10)
    liq = one(FT)
    rai = ifelse(x.q_rai < ϵ_empty, zero(FT), one(FT))
    ice = ifelse(x.q_ice < ϵ_empty, zero(FT), one(FT))
    return MicroState2MP3(liq, liq, rai, rai, ice, ice, ice, ice)
end

"""
    _rime_pair_b_new(q_new, b_trial, ρ_min, ρ_max)

The rime volume to pair with the already-floored rime mass `q_new`. `b_trial` is `b_rim`'s own
unfloored target, returned unchanged, bit for bit, whenever the implied density `q_new / b_trial`
already lies inside `[ρ_min, ρ_max]`, so the projection is inert on an admissible pair and
moves the pair to the nearest admissible one otherwise.
"""
@inline _rime_pair_b_new(q_new, b_trial, ρ_min, ρ_max) =
    UT.nearest_admissible_b(q_new, b_trial, ρ_min, ρ_max)

"""
    _condensate_phases(x::MicroState2MP3)

The condensed water of a [`MicroState2MP3`](@ref) split by phase: cloud liquid and rain
together, then ice. `q_rim` is a fraction of the ice content rather than an addition to it,
so it does not enter the split.
"""
@inline _condensate_phases(x::MicroState2MP3) = (x.q_lcl + x.q_rai, x.q_ice)

"""
    _apply_positivity_floor(x::MicroState2MP3, Δx, ρ_min, ρ_max)

Positivity floor for the substep increment `Δx` against state `x`, applied species-by-species
EXCEPT the rime mass/volume pair `(q_rim, b_rim)` (indices 7 and 8 of
[`MicroState2MP3`](@ref)), which is projected onto its admissible density cone as a pair rather
than floored independently ([`_rime_pair_b_new`](@ref)).

The pair projection is a property of this state, not a mode: the two rime moments are two
coordinates of one physical object, and flooring them independently admits pairs whose implied
density lies outside `[ρ_min, ρ_max]`. States that carry no such pair take the generic
positivity floor by dispatch.
"""
@inline function _apply_positivity_floor(
    x::MicroState2MP3{FT}, Δx, ρ_min, ρ_max,
) where {FT}
    idx_q, idx_b = 7, 8
    x_new = max.(x .+ Δx, 0)
    b_trial = x[idx_b] + Δx[idx_b]
    b_new = _rime_pair_b_new(x_new[idx_q], b_trial, ρ_min, ρ_max)
    return SA.setindex(x_new, b_new, idx_b)
end

#####
##### 2M+P3 Rosenbrock-Euler substepping (`RosenbrockAverage`)
#####
