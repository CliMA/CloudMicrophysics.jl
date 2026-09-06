#####
##### The 2M+P3 Rosenbrock-Euler march
#####

"""
    _substep_context(g, x)

The `(micro, thermo)` argument pair [`p3_2m_process_rates`](@ref) takes, built once
from a substep tendency callable `g` and the substep state `x`.

`micro` holds the state's own specific contents with `q_tot` promoted to the state's
element type; `thermo` holds the frozen substep context `(ρ, w, p, logλ)` together
with the temperature. Which temperature that is depends on the state: the
eight-species state has none of its own and takes the callable's frozen `T`, while
the temperature-coupled state supplies its ninth component. Every consumer of the
temperature downstream of this function, the per-process primal, the phase
relaxation context and the Jacobian alike, therefore reads it from one place.
"""
@inline function _substep_context(
    g::Instantaneous2MP3Tendency, x::SA.StaticVector{8, FT},
) where {FT}
    (q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim) = x
    micro = (; q_tot = FT(g.q_tot), q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    thermo = (; g.ρ, g.T, g.w, g.p, g.logλ)
    return micro, thermo
end

"""
    _tendency_and_jacobian(::ManualJacobian, g::Instantaneous2MP3Tendency, x)

The raw tendency and the [`ManualJacobian`](@ref) substep matrix from a single
[`p3_2m_process_rates`](@ref) evaluation: the tendency is the component-wise sum of
the per-process breakdown, and [`_jacobian_2mp3_manual`](@ref) consumes the same
breakdown for its donor coefficients. This shares the mixed-phase quadrature kernels
between the tendency and the Jacobian instead of replaying them once for each.
"""
@inline function _tendency_and_jacobian(
    ::ManualJacobian, g::Instantaneous2MP3Tendency, x::SA.StaticVector{8},
)
    micro, thermo = _substep_context(g, x)
    pp, rs = p3_2m_process_rates(g.mp, g.tps, micro, thermo)
    return (sum(values(pp)), _jacobian_2mp3_manual(g, x, pp, rs))
end

"""
    _per_process_rates(g::Instantaneous2MP3Tendency, x)

The bare [`p3_2m_process_rates`](@ref) breakdown at `x`, in [`MicroState2MP3`](@ref)'s own
eight-species space. See [`_per_process_rates`](@ref) for when this is reached.
"""
@inline function _per_process_rates(g::Instantaneous2MP3Tendency, x::SA.StaticVector{8})
    micro, thermo = _substep_context(g, x)
    pp, _ = p3_2m_process_rates(g.mp, g.tps, micro, thermo)
    return pp
end

#####
##### Temperature as a state variable
#####

"""
    MicroState2MP3T{FT}

The eight prognostic 2M+P3 species and the air temperature as a
`StaticArrays.FieldVector`. Internal to the [`RosenbrockAverage`](@ref)
implementation.

Carrying the temperature as a state component makes the latent-heat feedback part of
the same linear solve as the phase change that drives it, rather than an explicit
update applied after the solve has already committed to an increment. The entry
returns species tendencies only: the host owns the energy budget and derives its own
temperature from the mass tendencies.
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
    _positivity_mask(x)

Which slots of a substep state are floored at zero, as a per-slot boolean belonging
to the state type.

Masses and numbers vanish physically, so every species slot is floored. The
temperature is not a content that can run out: it is a thermodynamic coordinate whose
admissible range has nothing to do with zero, so the temperature slot of
[`MicroState2MP3T`](@ref) is exempt. Making the exemption a property of the state
type keeps the substep driver free of state-specific special cases, and keeps a new
state type from silently inheriting the wrong treatment.
"""
@inline _positivity_mask(::SA.StaticVector{N}) where {N} =
    SA.SVector(ntuple(_ -> true, Val(N)))
@inline _positivity_mask(::MicroState2MP3) = SA.SVector(ntuple(_ -> true, Val(8)))
@inline _positivity_mask(::MicroState2MP3T) =
    SA.SVector(ntuple(i -> i != 9, Val(9)))

"""
    _floored_state(x, Δx)

The state `x + Δx` with the floored slots of [`_positivity_mask`](@ref) clipped at
zero and the exempt slots left as they are.

Every place that commits an increment to a substep state goes through here, so the
positivity treatment is one implementation rather than a `max.(x .+ Δx, 0)` repeated
at each site with its own idea of which slots it applies to.
"""
@inline function _floored_state(x::SA.StaticVector{N}, Δx) where {N}
    y = x .+ Δx
    m = _positivity_mask(x)
    FT = eltype(y)
    return typeof(y)(
        ntuple(i -> @inbounds(ifelse(m[i], max(y[i], zero(FT)), y[i])), Val(N))...,
    )
end

"""
    _apply_positivity_floor(x::MicroState2MP3T, Δx, ρ_min, ρ_max)

The accepted increment `Δx` committed to the 2M+P3 state `x`: the per-slot
positivity floor of [`_floored_state`](@ref), except that the rime mass/volume pair
`(q_rim, b_rim)` at indices 7 and 8 is projected onto its admissible density cone as
a pair rather than floored independently.

`ρ_min`/`ρ_max` are [`CMP3.rime_density_bounds`](@ref). The projection returns `b_rim`
unchanged, bit for bit, whenever `q_rim / b_rim` already lies inside the cone, so it
is inert on admissible pairs and acts only where independent flooring would have left
the pair describing a rime density the scheme has no particle for.
"""
@inline function _apply_positivity_floor(x::MicroState2MP3T, Δx, ρ_min, ρ_max)
    idx_q, idx_b = 7, 8
    x_new = _floored_state(x, Δx)
    b_trial = x[idx_b] + Δx[idx_b]
    b_new = UT.nearest_admissible_b(x_new[idx_q], b_trial, ρ_min, ρ_max)
    return SA.setindex(x_new, b_new, idx_b)
end

"""
    _condensate_phases(x::MicroState2MP3T)

The condensed water of a [`MicroState2MP3T`](@ref) split by phase, as for
[`MicroState2MP3`](@ref). `q_rim` is a fraction of the ice content rather than an
addition to it, and the temperature is not a water content, so neither enters the split.
"""
@inline _condensate_phases(x::MicroState2MP3T) = (x.q_lcl + x.q_rai, x.q_ice)

"""
    _rosenbrock_species_mask(x::MicroState2MP3T)

The near-empty species projection of [`MicroState2MP3`](@ref) extended to the
temperature-coupled state. The temperature is never masked out of the implicit solve:
its coupling to the phase-change rows is the reason the state includes it.
"""
@inline function _rosenbrock_species_mask(x::MicroState2MP3T{FT}) where {FT}
    ϵ_empty = FT(1e-10)
    liq = one(FT)
    rai = ifelse(x.q_rai < ϵ_empty, zero(FT), one(FT))
    ice = ifelse(x.q_ice < ϵ_empty, zero(FT), one(FT))
    return MicroState2MP3T(liq, liq, rai, rai, ice, ice, ice, ice, one(FT))
end

"""
    _latent_heating(Lᵥ, Lₛ, cp_air, s)

The temperature response to a 2M+P3 species change `s`: the latent heat released by
the liquid part of `s` (cloud plus rain) and by its ice part, divided by the moist
heat capacity.

This is the whole temperature policy of the 2M+P3 substep, in one expression. It is
linear in `s`, which is what lets the two ways of using it stay one policy rather than
two: applied to the tendency it is the rate of the temperature slot
([`_temperature_2mp3_tendency`](@ref)), and applied to an accepted increment it is the
between-substeps temperature march of the eight-species state
([`_marched_temperature`](@ref)).
"""
@inline _latent_heating(Lᵥ, Lₛ, cp_air, s::MicroState2MP3) =
    (Lᵥ * (s.q_lcl + s.q_rai) + Lₛ * s.q_ice) / cp_air

"""
    _phase_relaxation_context(mp, tps, micro, thermo)

The shared phase-change quantities of the temperature-coupled tendency and its
Jacobian: capacitance timescales, saturation excesses, saturation-slope derivatives,
psychrometric factors, latent heats, the moist heat capacity, and the melting
contribution's temperature derivative. Evaluated at the mean-mass-bounded cloud
population and the nonnegative-clamped ice state, so it takes the same branches at the
same state as the per-process primal does.
"""
@inline function _phase_relaxation_context(mp, tps, micro, thermo)
    (; q_tot, q_lcl, n_lcl, q_rai, q_ice, n_ice, q_rim, b_rim) = micro
    (; T, logλ) = thermo
    FT = typeof(q_tot)
    # The air density appears in denominators, inside logs and under negative
    # fractional powers, so it takes the positive floor rather than a non-negative
    # clamp. See [`AIR_DENSITY_FLOOR`](@ref).
    ρ = _floored_air_density(thermo.ρ)

    Rᵥ = TDI.Rᵥ(tps)
    Lᵥ = TDI.Lᵥ(tps, T)
    Lₛ = TDI.Lₛ(tps, T)
    cp_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_ice)
    qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_ice)

    pdf_c = mp.warm_rain.seifert_beheng.pdf_c
    # the zero-mass droplet arm's deactivation test, on the clamped state, as in the primal
    sat_excess_arm = _liquid_sat_excess(
        tps, ρ, T, q_tot,
        UT.clamp_to_nonneg(q_lcl), UT.clamp_to_nonneg(q_rai), UT.clamp_to_nonneg(q_ice))
    n_lcl_b = CM2.number_bounded_by_mass_limits(
        (; x_min = pdf_c.xc_min, x_max = pdf_c.xc_max),
        UT.clamp_to_nonneg(q_lcl), UT.clamp_to_nonneg(n_lcl), sat_excess_arm;
        invent_from_zero = false)
    τ_l = CM2.cloud_condensation_timescale(
        pdf_c, mp.warm_rain.air_properties, tps, T, ρ,
        UT.clamp_to_nonneg(q_lcl), n_lcl_b * ρ)
    qᵥ_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
    dqsl_dT = CMNonEq.dqcld_dT(qᵥ_sat_liq, Lᵥ, Rᵥ, T)
    Γₗ = CMNonEq.gamma_helper(Lᵥ, cp_air, dqsl_dT)
    sat_excess_l = qᵥ - qᵥ_sat_liq

    state_i = CMP3.state_from_prognostic(
        mp.ice.scheme,
        UT.clamp_to_nonneg(q_ice) * ρ, UT.clamp_to_nonneg(n_ice) * ρ,
        UT.clamp_to_nonneg(q_rim) * ρ, UT.clamp_to_nonneg(b_rim) * ρ)
    τ_i = CMP3.ice_deposition_timescale(
        mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps, T, ρ,
        state_i, logλ; quad = mp.ice.quad)
    qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
    dqsi_dT = CMNonEq.dqcld_dT(qᵥ_sat_ice, Lₛ, Rᵥ, T)
    Γᵢ = CMNonEq.gamma_helper(Lₛ, cp_air, dqsi_dT)
    sat_excess_i = qᵥ - qᵥ_sat_ice

    # The melting contribution's temperature column, evaluated at the clamped ice state
    # the primal builds and under the same two conditions it applies,
    # `ice_population_is_present` and `T > T_freeze`, so the column differentiates the
    # contribution the tendency sums and cannot run a branch the rate does not.
    # `ice_melt` returns the derivatives of its rates alongside them, and
    # `_ice_melting_species` is linear in the rates.
    T_freeze = TDI.T_freeze(tps)
    state_melt = CMP3.state_from_prognostic(
        mp.ice.scheme,
        UT.clamp_to_nonneg(q_ice) * ρ, UT.clamp_to_nonneg(n_ice) * ρ,
        UT.clamp_to_nonneg(q_rim) * ρ, UT.clamp_to_nonneg(b_rim) * ρ)
    ∂melting_∂T = if CMP3.ice_population_is_present(state_melt)
        melt = ifelse(T > T_freeze,
            CMP3.ice_melt(mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps,
                T, ρ, state_melt, logλ; quad = mp.ice.quad),
            CMP3.zero_ice_melt(ρ),
        )
        _ice_melting_species(ρ, mp.ice.scheme.ρ_i, UT.clamp_to_nonneg(q_rim),
            UT.clamp_to_nonneg(b_rim), melt.∂dNdt_∂T, melt.∂dLdt_∂T, melt.∂melt_frac_∂T)
    else
        MicroState2MP3{FT}(0, 0, 0, 0, 0, 0, 0, 0)
    end

    return (; Rᵥ, Lᵥ, Lₛ, cp_air, τ_l, dqsl_dT, Γₗ, sat_excess_l,
        τ_i, dqsi_dT, Γᵢ, sat_excess_i, ∂melting_∂T, T_freeze)
end

"""
    Temperature2MP3Tendency(mp, tps, ρ, q_tot, logλ, w = zero(ρ), p = zero(ρ))

Callable bundling the frozen per-substep context of the temperature-coupled state;
applying it to a [`MicroState2MP3T`](@ref) returns the species rates together with the
latent-heating temperature rate.

The temperature is absent from the frozen context because it is a component of the
state the callable is applied to. `ρ`, `q_tot`, `logλ`, `w` and `p` are frozen across
the substep, as they are for [`Instantaneous2MP3Tendency`](@ref); `w` and `p` reach
the adiabatic-parcel branch of the activation supersaturation and nothing else.
"""
struct Temperature2MP3Tendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    q_tot::F
    logλ::F
    w::F
    p::F
end
Temperature2MP3Tendency(mp, tps, ρ, q_tot, logλ) =
    Temperature2MP3Tendency(mp, tps, ρ, q_tot, logλ, zero(ρ), zero(ρ))

@inline function _substep_context(
    g::Temperature2MP3Tendency, y::SA.StaticVector{9, FT},
) where {FT}
    (q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, T) = y
    micro = (; q_tot = FT(g.q_tot), q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    thermo = (; g.ρ, T, g.w, g.p, g.logλ)
    return micro, thermo
end

"""
    _temperature_2mp3_tendency(pp, ctx)

The temperature-coupled tendency [`MicroState2MP3T`](@ref) from a
[`p3_2m_process_rates`](@ref) breakdown `pp` and a
[`_phase_relaxation_context`](@ref) `ctx`: the component-wise sum of `pp`, with the
latent-heating temperature rate of [`_latent_heating`](@ref) appended.

`pp`'s condensation and deposition slots carry the bare relaxation rate, without the
psychrometric factor `Γ`, because the coupled temperature carries that feedback
instead. See [`_bare_rate_conv_q_vap_to_q_lcl`](@ref).
"""
@inline function _temperature_2mp3_tendency(pp, ctx)
    f8 = sum(values(pp))
    fT = _latent_heating(ctx.Lᵥ, ctx.Lₛ, ctx.cp_air, f8)
    return MicroState2MP3T(Tuple(f8)..., fT)
end

@inline function (g::Temperature2MP3Tendency)(y::SA.StaticVector{9})
    micro, thermo = _substep_context(g, y)
    pp, _ = p3_2m_process_rates(g.mp, g.tps, micro, thermo)
    ctx = _phase_relaxation_context(g.mp, g.tps, micro, thermo)
    return _temperature_2mp3_tendency(pp, ctx)
end

"""
    _tendency_and_jacobian(::TemperatureCoupledJacobian, g::Temperature2MP3Tendency, y)

The temperature-coupled tendency and 9x9 substep Jacobian from a single
[`p3_2m_process_rates`](@ref) and [`_phase_relaxation_context`](@ref) evaluation,
shared between [`_temperature_2mp3_tendency`](@ref) and
[`_jacobian_2mp3t_manual`](@ref) instead of each replaying the mixed-phase quadrature
kernels independently.
"""
@inline function _tendency_and_jacobian(
    ::TemperatureCoupledJacobian, g::Temperature2MP3Tendency, y::SA.StaticVector{9},
)
    micro, thermo = _substep_context(g, y)
    pp, rs = p3_2m_process_rates(g.mp, g.tps, micro, thermo)
    ctx = _phase_relaxation_context(g.mp, g.tps, micro, thermo)
    return (_temperature_2mp3_tendency(pp, ctx), _jacobian_2mp3t_manual(g, y, pp, rs, ctx))
end

"""
    _per_process_rates(g::Temperature2MP3Tendency, y)

The [`p3_2m_process_rates`](@ref) breakdown at `y`, extended one component per process with
that process's own latent-heating row ([`_latent_heating`](@ref)), so the temperature
contribution of a process reaches the sink alongside its mass contribution. Linear in the
per-process rate, exactly as [`_temperature_2mp3_tendency`](@ref)'s summed one is, so the
per-process temperature rows sum to the same net latent heating. See
[`_per_process_rates`](@ref) for when this is reached.
"""
@inline function _per_process_rates(g::Temperature2MP3Tendency, y::SA.StaticVector{9})
    micro, thermo = _substep_context(g, y)
    pp, _ = p3_2m_process_rates(g.mp, g.tps, micro, thermo)
    ctx = _phase_relaxation_context(g.mp, g.tps, micro, thermo)
    return map(
        s -> MicroState2MP3T(Tuple(s)..., _latent_heating(ctx.Lᵥ, ctx.Lₛ, ctx.cp_air, s)),
        pp,
    )
end

"""
    _apply_limiter(::EndStateSaturationAdjustment, x::MicroState2MP3T, d, ρ, Tsub,
                   q_tot, Lv_over_cp, Ls_over_cp, tps)

The end-state saturation adjustment at the temperature-coupled state.

The eight-species method has to reconstruct the end-state temperature from the
increment's latent heating, at constant latent heats over dry air. Here the increment
already includes its own temperature component, marched by the same linear solve that
produced the mass increment, so the adjustment reads the end-state temperature off the
state instead of approximating it. `Lv_over_cp` and `Ls_over_cp` are consequently
unused and are kept only so the driver calls every limiter the same way.

See [`_apply_limiter_diag`](@ref) for the scale-factor-exposing form.
"""
@inline _apply_limiter(::EndStateSaturationAdjustment, x::MicroState2MP3T, d::SA.StaticVector{9},
    ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps,
) = first(_apply_limiter_diag(
    EndStateSaturationAdjustment(), x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps))

@inline function _apply_limiter_diag(::EndStateSaturationAdjustment,
    x::MicroState2MP3T{FT}, d::SA.StaticVector{9},
    ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps,
) where {FT}
    Ssat(xx, TT) = max(
        TDI.supersaturation_over_ice(tps, q_tot, xx.q_lcl + xx.q_rai, xx.q_ice, ρ, TT),
        TDI.supersaturation_over_liquid(tps, q_tot, xx.q_lcl + xx.q_rai, xx.q_ice, ρ, TT),
    )
    latent(dd) = dd.T
    return _saturation_bisection_diag(Ssat, latent, x, d, Tsub)
end

#####
##### The substep march
#####

"""
    _refreshed_logλ(mp, ρ, x)

The ice size-distribution shape parameter re-solved from the marched state `x`.

The substep march moves the ice slots, and the shape parameter is a function of them;
a value carried from the entry state therefore describes a distribution that is not
the one the later substeps are acting on. This is safe at any state the march can
reach, including degenerate ones: [`CMP3.state_from_prognostic`](@ref) projects the
rime pair onto its admissible cone, and [`CMP3.get_distribution_logλ`](@ref) returns
the correct bracket end where the mean mass is degenerate rather than an arbitrary one.
"""
@inline _refreshed_logλ(mp, ρ, x::Union{MicroState2MP3, MicroState2MP3T}) =
    last(_ice_state_and_logλ(mp, ρ, x))

"""
    _ice_state(mp, ρ, x)

The P3 ice state of the substep state `x`, from its four ice slots.

Built once per substep and shared by the shape refresh and the fall speeds, so a rate and
a velocity cannot be evaluated on two different ice states.
"""
@inline _ice_state(mp, ρ, x::Union{MicroState2MP3, MicroState2MP3T}) =
    CMP3.state_from_prognostic(
        mp.ice.scheme, ρ * x.q_ice, ρ * x.n_ice, ρ * x.q_rim, ρ * x.b_rim)

"""
    _ice_state_and_logλ(mp, ρ, x)

The ice state of `x` and its re-solved shape parameter, as a pair, so the state the solve
already built is returned rather than rebuilt by the caller.
"""
@inline function _ice_state_and_logλ(mp, ρ, x::Union{MicroState2MP3, MicroState2MP3T})
    state = _ice_state(mp, ρ, x)
    return (state, CMP3.get_distribution_logλ(state))
end

"""
    _substep_fall_speeds(mp, ρ, x, state, logλ)

The six sedimentation velocities [m/s] at one substep state, as
`(v_lcl_n, v_lcl_m, v_rai_n, v_rai_m, v_ice_n, v_ice_m)`: cloud liquid, rain and P3 ice,
each number-weighted then mass-weighted, all positive downward.

They are evaluated on the substep's OWN state and its own ice shape parameter, which is
the whole reason they are computed here rather than by the host. A host fill reads the
step-entry state, so at any substep count above one the velocity a species sediments with
belongs to a distribution the microphysics has already moved away from; and it reads the
state as the host holds it rather than as the kernel canonicalizes it, so a velocity and a
rate could be evaluated on two different states in the same cell. Both are structural
rather than incidental: no host-side fill can see a substep, and none can see the clamps.

The canonicalization is the primal's own
([`p3_2m_process_rates`](@ref)): the air density takes the positive floor
[`AIR_DENSITY_FLOOR`](@ref), because it sits in a denominator of the Stokes prefactor and
of every specific-to-volumetric conversion here, and the four warm moments take the
non-negative clamp. The ice quartet needs neither, [`CMP3.state_from_prognostic`](@ref)
clamping the moments and projecting the rime pair onto its density cone already.
"""
@inline function _substep_fall_speeds(mp, ρ, x, state, logλ)
    sb = mp.warm_rain.seifert_beheng
    ρₐ = _floored_air_density(ρ)
    q_lcl = UT.clamp_to_nonneg(x.q_lcl)
    q_rai = UT.clamp_to_nonneg(x.q_rai)
    # The warm entries take a NUMBER CONCENTRATION where the substep state carries a
    # specific number, which is the same conversion the primal makes at every warm rate.
    N_lcl = ρₐ * UT.clamp_to_nonneg(x.n_lcl)
    N_rai = ρₐ * UT.clamp_to_nonneg(x.n_rai)
    v_lcl_n, v_lcl_m =
        CM2.cloud_terminal_velocity(sb.pdf_c, mp.warm_rain.cloud_velocity, q_lcl, ρₐ, N_lcl)
    v_rai_n, v_rai_m =
        CM2.rain_terminal_velocity(sb, mp.warm_rain.rain_velocity, q_rai, ρₐ, N_rai)
    v_ice_n = CMP3.ice_terminal_velocity_number_weighted(
        mp.ice.terminal_velocity, ρₐ, state, logλ; quad = mp.ice.quad)
    v_ice_m = CMP3.ice_terminal_velocity_mass_weighted(
        mp.ice.terminal_velocity, ρₐ, state, logλ; quad = mp.ice.quad)
    return (v_lcl_n, v_lcl_m, v_rai_n, v_rai_m, v_ice_n, v_ice_m)
end

"""
    _substep_tendency(x, mp, tps, ρ, Tsub, q_tot, logλ, w, p)

The tendency callable for one substep at state `x`: [`Instantaneous2MP3Tendency`](@ref)
for the eight-species state, which needs the substep temperature frozen into it, and
[`Temperature2MP3Tendency`](@ref) for the temperature-coupled state, which reads its
temperature from the state instead.
"""
@inline _substep_tendency(::MicroState2MP3, mp, tps, ρ, Tsub, q_tot, logλ, w, p) =
    Instantaneous2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ, w, p)
@inline _substep_tendency(::MicroState2MP3T, mp, tps, ρ, Tsub, q_tot, logλ, w, p) =
    Temperature2MP3Tendency(mp, tps, ρ, q_tot, logλ, w, p)

"""
    _marched_temperature(x, x_prev, tps, q_tot, Tsub)

The substep temperature after the increment `x - x_prev` has been committed.

The temperature-coupled state has already marched it: the ninth component came out of
the same linear solve as the mass increment, was scaled by the same water bound and
limiter, and is returned as it stands. The eight-species state has no temperature
component, so it applies the same latent-heating map ([`_latent_heating`](@ref)) to the
realized increment afterwards, at the moist heat capacity of the pre-substep state.

The two are one policy, evaluated at different arguments: the coupled state integrates
the map's output as a rate inside the solve, the eight-species state evaluates it on
the accepted increment after the solve. They do not produce the same number, and they
are not meant to: an increment that the limiter or the water bound has reduced gives the
reduced latent heating in both, but only the coupled state feeds the temperature back
into the phase-change rows it came from.
"""
@inline function _marched_temperature(
    x::MicroState2MP3, x_prev::MicroState2MP3, tps, q_tot, Tsub,
)
    # The latent heats are evaluated at a temperature floored into the domain over which
    # their linear fits are meaningful; the substep temperature can leave it transiently
    # on a fallback branch.
    T_safe = max(150, Tsub)
    cp_air = TDI.cpₘ(tps, q_tot, x_prev.q_lcl + x_prev.q_rai, x_prev.q_ice)
    return Tsub +
           _latent_heating(TDI.Lᵥ(tps, T_safe), TDI.Lₛ(tps, T_safe), cp_air, x - x_prev)
end
@inline _marched_temperature(
    x::MicroState2MP3T, ::MicroState2MP3T, tps, q_tot, Tsub,
) = x.T

"""
    record!(sink, ctx)

Hand one substep's context to a diagnostic sink. [`_rosenbrock_substep_diag`](@ref) calls
this once per substep, with `sink` unchanged from whatever [`_march_2mp3`](@ref) was given.
The production entries pass [`NullSink`](@ref); this method covers `_march_2mp3`'s own
`sink = nothing` default, for a caller that marches without importing `NullSink` at all.
Either way `record!`'s body is empty, so it compiles away and the physics path and the
diagnostic path are the same code and cannot diverge.
"""
@inline record!(::Nothing, ctx) = nothing

"""
    _march_2mp3(mode, mp, tps, ρ, T, q_tot, x₀, logλ, Δt, nsub, w, p, sink = nothing)

The 2M+P3 Rosenbrock-Euler march from `x₀` over `Δt` in `nsub` substeps, returning the
final state and the substep-averaged sedimentation velocities
([`_substep_fall_speeds`](@ref)) as an `SVector{6}`.

One implementation serves both substep states. The state type selects the tendency
callable ([`_substep_tendency`](@ref)), the temperature policy
([`_marched_temperature`](@ref)), the positivity treatment
([`_positivity_mask`](@ref)) and the accepted-increment projection
([`_apply_positivity_floor`](@ref)); everything else, the substep solve
([`_rosenbrock_substep_diag`](@ref)) with its water bound, increment limiter and
fallback branches, and the per-substep shape refresh, is shared by construction. A
state size therefore cannot acquire or lose an ancillary treatment by being marched in
a loop of its own.

`q_tot`, `ρ`, `w` and `p` are held fixed across substeps.
"""
@inline function _march_2mp3(
    mode::RosenbrockAverage, mp, tps, ρ, T, q_tot, x₀::SA.StaticVector{N, FT}, logλ,
    Δt, nsub, w, p, sink = nothing,
) where {N, FT}
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    cp_d = TDI.TD.Parameters.cp_d(tps)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / cp_d
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / cp_d
    ρ_min, ρ_max = CMP3.rime_density_bounds(mp.ice.scheme)

    x = x₀
    Tsub = T
    v_acc = zero(SA.SVector{6, FT})
    for i in 1:nsub_eff
        # THE SHAPE PARAMETER IS REFRESHED FROM THE MARCHED STATE, every substep.
        #
        # Held fixed, it imposes an error the march cannot remove: more substeps only
        # re-use the same stale shape more often, so the error plateaus instead of
        # converging. Measured at the box's 2 s step that floor is 1.6 to 1.8 percent;
        # measured at the coupled 180 s step it is a FACTOR of 2.8 to 4 on converged ice
        # number, which is what makes it worth a solve per substep.
        #
        # It is an accuracy floor and NOT a convergence fix, and the two must not be
        # conflated: the ice-number march fails to converge at a large minority of
        # carrying cells whether the shape parameter is refreshed or not, and that
        # pathology lives in the rates.
        #
        # The first substep re-uses the caller's value rather than re-solving: it already
        # belongs to the entry state, so this spends no solve where the answer is in hand
        # and is exactly inert at `nsub == 1`, where there is no later substep for a stale
        # shape to reach.
        state_sub, logλ_sub =
            i == 1 ? (_ice_state(mp, ρ, x), logλ) : _ice_state_and_logλ(mp, ρ, x)
        # THE FALL SPEEDS ARE SAMPLED HERE, on the state the substep is about to act on and
        # with that substep's own shape parameter, and averaged over the march. The substeps
        # are uniform in length, `h = Δt / nsub_eff`, so the arithmetic mean IS the time
        # average over `Δt` and no weighting is owed. Sampling at loop entry rather than at
        # loop exit is what pairs each velocity with the state its own tendency saw.
        v_acc = v_acc + SA.SVector{6, FT}(_substep_fall_speeds(mp, ρ, x, state_sub, logλ_sub))
        g = _substep_tendency(x, mp, tps, ρ, Tsub, q_tot, logλ_sub, w, p)
        x_prev = x
        x, diag = _rosenbrock_substep_diag(
            mode, g, x, h, q_tot, ρ, Tsub, Lv_over_cp, Ls_over_cp, tps, ρ_min, ρ_max, sink)
        Tsub = _marched_temperature(x, x_prev, tps, q_tot, Tsub)
    end
    return (x, v_acc / FT(nsub_eff))
end

"""
    _species_increment(d)

The eight species components of a substep increment, dropping the temperature
component of a [`MicroState2MP3T`](@ref) increment. The entry returns species
tendencies only: the host owns the energy budget and derives its own temperature from
the mass tendencies, so returning a temperature tendency here would give it a second,
uncoordinated source.
"""
@inline _species_increment(d::MicroState2MP3) = d
@inline _species_increment(d::MicroState2MP3T) =
    MicroState2MP3(d.q_lcl, d.n_lcl, d.q_rai, d.n_rai, d.q_ice, d.n_ice, d.q_rim, d.b_rim)

"""
    _rosenbrock_average_carrier(rates, extras)

The fixed-shape return value of the 2M+P3 [`RosenbrockAverage`](@ref) entries: the
eight species tendencies, the `dn_lcl_activation_dt` slot of the host's microphysics
tendency cache, and an `extras` `NamedTuple`.

`dn_lcl_activation_dt` is zero because droplet activation is a per-process slot of the
substep: its number source is already inside `dn_lcl_dt` and its paired mass inside
`dq_lcl_dt`, so a host that added anything of its own here would double the source. The
field is kept rather than removed so the returned shape still matches the host's cache
type.

`extras` holds substep by-products that are not tendencies: the six substep-averaged
sedimentation velocities `v_lcl_n`, `v_lcl_m`, `v_rai_n`, `v_rai_m`, `v_ice_n` and
`v_ice_m` [m/s], all positive downward, from [`_substep_fall_speeds`](@ref).

They ride here rather than beside the tendencies because they are not tendencies: a host
applies them to a state, where it integrates the eight rates. Keeping them in a
`NamedTuple` of its own is also what let them arrive without changing this function's
signature or any host unpacking, which is what the slot was cut for.
"""
@inline _rosenbrock_average_carrier(rates::MicroState2MP3{FT}, extras::NamedTuple) where {FT} =
    (;
        dq_lcl_dt = rates.q_lcl, dn_lcl_dt = rates.n_lcl,
        dq_rai_dt = rates.q_rai, dn_rai_dt = rates.n_rai,
        dq_ice_dt = rates.q_ice, dn_ice_dt = rates.n_ice,
        dq_rim_dt = rates.q_rim, db_rim_dt = rates.b_rim,
        dn_lcl_activation_dt = zero(FT),
        extras,
    )

"""
    _rosenbrock_average_2mp3(mode, mp, tps, ρ, T, q_tot, x₀, logλ, Δt, nsub, w, p,
        sink = nothing)

March `x₀` with [`_march_2mp3`](@ref) and return the average tendency over `Δt` in the
carrier of [`_rosenbrock_average_carrier`](@ref). Shared by the entries of both
substep states. `sink` passes straight through to [`_march_2mp3`](@ref); the production
entries pass [`NullSink`](@ref).
"""
@inline function _rosenbrock_average_2mp3(
    mode, mp, tps, ρ, T, q_tot, x₀, logλ, Δt, nsub, w, p, sink = nothing,
)
    x, v̄ = _march_2mp3(mode, mp, tps, ρ, T, q_tot, x₀, logλ, Δt, nsub, w, p, sink)
    return _rosenbrock_average_carrier(
        _species_increment(x - x₀) / Δt,
        (;
            v_lcl_n = v̄[1], v_lcl_m = v̄[2],
            v_rai_n = v̄[3], v_rai_m = v̄[4],
            v_ice_n = v̄[5], v_ice_m = v̄[6],
        ),
    )
end

#####
##### Entries
#####

bulk_microphysics_tendencies(::RosenbrockAverage, ::Microphysics2Moment, args...) = throw(
    ArgumentError(
        "RosenbrockAverage on the 2M+P3 model supports only TemperatureCoupledJacobian, ManualJacobian or ExactJacobian; use rosenbrock_manual_temperature(), rosenbrock_manual() or rosenbrock_exact()",
    ),
)

bulk_microphysics_tendencies(
    ::RosenbrockAverage, ::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, Nothing}, args...,
) where {WR} = throw(
    ArgumentError(
        "RosenbrockAverage on Microphysics2Moment requires P3 ice parameters (with_ice = true)",
    ),
)

"""
    bulk_microphysics_tendencies(mode::RosenbrockAverage{TemperatureCoupledJacobian},
        ::Microphysics2Moment, mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub = 1, w = zero(ρ), p = zero(ρ))

Average 2M+P3 microphysics tendencies over `Δt` from `nsub` linearized-implicit
(Rosenbrock-Euler) substeps of the temperature-coupled state
`(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, T)`. This is the production
configuration: the temperature is a state variable of the substepper, so the latent
heating of a phase change and the saturation it changes are solved together instead of
in sequence.

`q_tot` is held fixed across substeps, and so are the host vertical velocity `w` and
air pressure `p`, which only droplet activation reads and only through the
adiabatic-parcel branch of the supersaturation it activates at. `logλ` is the entry
state's ice shape parameter and is refreshed from the marched state at every later
substep. See the [Rosenbrock-average microphysics substepping](@ref) documentation page
for the substep algorithm.

Returns the fixed-shape carrier of [`_rosenbrock_average_carrier`](@ref), whose `extras`
carry the six substep-averaged sedimentation velocities. The temperature the substeps
marched is not returned: the host derives its own from the species tendencies.

`sink` is production's [`NullSink`](@ref) by default; [`Verbose`](@ref) and
[`Trace`](@ref) pass their own through this same keyword.
"""
@inline function bulk_microphysics_tendencies(
    mode::RosenbrockAverage{TemperatureCoupledJacobian}, cm::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    Δt, nsub = 1, w = zero(ρ), p = zero(ρ); sink = NullSink(),
) where {WR, ICE <: CMP.P3IceParams}
    FT = typeof(q_tot)
    y₀ = MicroState2MP3T{FT}(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, T)
    return _rosenbrock_average_2mp3(
        mode, mp, tps, ρ, FT(T), q_tot, y₀, FT(logλ), Δt, nsub, FT(w), FT(p), sink)
end

"""
    bulk_microphysics_tendencies(mode::RosenbrockAverage{<:Union{ExactJacobian, ManualJacobian}},
        ::Microphysics2Moment, mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub = 1, w = zero(ρ), p = zero(ρ))

Average 2M+P3 microphysics tendencies over `Δt` from `nsub` linearized-implicit
(Rosenbrock-Euler) substeps of the eight-species state, with the temperature frozen
within a substep and marched from the realized increment between substeps
([`_marched_temperature`](@ref)).

The substep march is the one
[`bulk_microphysics_tendencies(::RosenbrockAverage{TemperatureCoupledJacobian}, ...)`](@ref)
runs, at eight components rather than nine, so the two differ in the temperature
treatment alone. `q_tot`, `w` and `p` are held fixed across substeps and `logλ` is
refreshed from the marched state, as they are there. The donor-based matrices
([`DonorJacobian`](@ref), [`CoupledDonorJacobian`](@ref)) are 1M-only.

Returns the fixed-shape carrier of [`_rosenbrock_average_carrier`](@ref).

`sink` is production's [`NullSink`](@ref) by default; [`Verbose`](@ref) and
[`Trace`](@ref) pass their own through this same keyword.
"""
@inline function bulk_microphysics_tendencies(
    mode::RosenbrockAverage{<:Union{ExactJacobian, ManualJacobian}}, cm::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    Δt, nsub = 1, w = zero(ρ), p = zero(ρ); sink = NullSink(),
) where {WR, ICE <: CMP.P3IceParams}
    FT = typeof(q_tot)
    x₀ = MicroState2MP3{FT}(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    return _rosenbrock_average_2mp3(
        mode, mp, tps, ρ, T, q_tot, x₀, logλ, Δt, nsub, FT(w), FT(p), sink)
end
