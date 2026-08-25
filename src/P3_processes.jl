"""
    het_ice_nucleation(aerosol, tps, q_lcl, N_lcl, RH, T, ρₐ)

Calculate the ice nucleation rate from heterogeneous freezing due to some `aerosol`

# Arguments
  - `aerosol`: aerosol parameters (supported types: desert dust, illite, kaolinite)
  - `tps`: thermodynamics parameters
  - `q_lcl`: cloud liquid water specific content
  - `N_lcl`: cloud droplet number concentration
  - `RH`: relative humidity
  - `T`: temperature
  - `ρₐ`: air density

# Returns
- A `NamedTuple` with the fields:
  - `dNdt`: ice number concentration change rate [m⁻³ s⁻¹]
  - `dLdt`: ice content change rate [kg m⁻³ s⁻¹]
"""
function het_ice_nucleation(
    aerosol::Union{CMP.DesertDust, CMP.Illite, CMP.Kaolinite},
    tps::TDI.PS,
    q_lcl, N_lcl, RH, T, ρₐ,
)
    FT = eltype(tps)
    #TODO - Also consider rain freezing

    # Immersion freezing nucleation rate coefficient
    J = CM_HetIce.ABIFM_J(aerosol, RH - CO.a_w_ice(tps, T))

    # Assumed erosol surface area
    # TODO - Make it a parameter of ABIFM scheme
    # We could consider making it a function of the droplet size distribution
    A_aer = FT(1e-10)

    # NaN guard: if ABIFM returns a non-finite J for a degenerate input,
    # treat non-finite J as "no nucleation" (0 rate), consistent with
    # the physical limit of an unresolvable state.
    JA_aer = ifelse(isfinite(J), J * A_aer, zero(J))

    dNdt = max(0, JA_aer * N_lcl)
    dLdt = max(0, JA_aer * q_lcl * ρₐ)

    return (; dNdt, dLdt)
end

"""
    _ventilation_from(velocity_params, aps, ρₐ, state, logλ, quad)

Dispatch point for the ventilation integral: evaluate it by quadrature when `quad` is a
quadrature rule, or read it from a lookup table when the mode's carrier holds one. Both
consumers of the ventilation integral call this rather than integrating directly, because a
quadrature rule integrates an integrand, while a table replaces the whole integral and is
keyed on the state, which by the point of integration has already been captured into closures.
"""
@inline _ventilation_from(velocity_params, aps, ρₐ, state, logλ, quad) =
    ice_ventilation_integral(velocity_params, aps, ρₐ, state, logλ; quad)

"""
    ice_ventilation_integral(velocity_params, aps, ρₐ, state, logλ; quad)

The ventilation integral `∫ D · F_v(D) · N′(D) dD` over the ice size distribution, shared by
[`ice_melt`](@ref) and [`ice_deposition_timescale`](@ref), which differ only in the scalar
prefactor each applies to it: `2πK/L_f · ΔT` for melt, `2πG` and `qᵥ_sat_ice` for deposition.

Temperature enters neither consumer through this integral: it enters only through those
prefactors, which is why a four-axis table in `(x̄, F_rim, ρ_rim, ρₐ)` can serve both without a
temperature axis.
"""
@inline function ice_ventilation_integral(velocity_params, aps::CMP.AirProperties, ρₐ, state::P3State, logλ; quad)
    (; vent) = state.params
    v_term = ice_particle_terminal_velocity(velocity_params, ρₐ, state)
    F_v = CO.ventilation_factor(vent, aps, v_term)
    N′ = size_distribution(state, logλ)
    bnds = velocity_integral_bounds(state, logλ, v_term; p = 1e-6)
    integrate(D -> D * F_v(D) * N′(D), bnds, quad)
end

"""
    ice_melt(velocity_params, aps, tps, Tₐ, ρₐ, state, logλ; ∫kwargs...)

# Arguments
 - `velocity_params`: [`CMP.Chen2022VelType`](@ref)
 - `aps`: [`CMP.AirProperties`](@ref)
 - `tps`: thermodynamics parameters
 - `Tₐ`: temperature (K)
 - `ρₐ`: air density
 - `state`: a [`P3State`](@ref) object
 - `logλ`: the log of the slope parameter [log(1/m)]

# Keyword arguments
 - `quad`: quadrature rule (a `Quadrature.QuadratureRule`)

Returns the melting rate of ice (QIMLT in [MorrisonMilbrandt2015](@cite)) as
`(; dNdt, dLdt, melt_frac, ∂dNdt_∂T, ∂dLdt_∂T, ∂melt_frac_∂T)`, the number and mass melting
rates, the fractional ice-mass melting rate they share, and the derivatives of all three with
respect to `Tₐ`.

`melt_frac` is [`ice_melt_fraction`](@ref): `dLdt / ρq_ice`, bounded above by
[`ice_melt_fraction_limit`](@ref). The number melting rate is `ρn_ice * melt_frac`, and the
same fraction drains the rime pair at the caller, so number and rime leave at one bounded
fraction.

The mass rate is the conduction approximation with spherical capacitance ``C = D/2``,

```math
\\frac{dL}{dt} = \\frac{2 π K_{therm} (Tₐ - T_{freeze})}{L_f} ∫ D \\, F_v(D) \\, N'(D) \\, dD,
```

the same capacitance integral as `ice_deposition_timescale`, following
Morrison and Milbrandt (2015). The vapor diffusion contribution in subsaturated air of
Morrison and Milbrandt (2015) is not included.

The conduction integral carries no temperature dependence at all, and the prefactor
`2 π K_therm (Tₐ - T_freeze) / L_f(Tₐ)` depends on temperature only through the excess and
through `L_f`, which is linear with slope `cp_l - cp_i`. Both derivatives are therefore exact
and need no division by the temperature excess. Dropping the `L_f` term instead would cost
`(cp_l - cp_i)/L_f ≈ 6e-3` per K of excess, which is 6% at 10 K.
"""
@inline function ice_melt(
    velocity_params, aps::CMP.AirProperties, tps::TDI.PS,
    Tₐ, ρₐ, state::P3State, logλ;
    quad,
)
    # Note: process not dependent on `F_liq`
    # (we want ice core shape params)
    # Get constants
    FT = eltype(state)
    (; K_therm) = aps
    L_f = TDI.Lf(tps, Tₐ)

    (; ρq_ice, ρn_ice) = state
    (; T_freeze) = state.params

    # The ventilation integral is shared with `ice_deposition_timescale`; the temperature
    # dependence is entirely in the prefactors below.
    ΔT = Tₐ - T_freeze
    fac₀ = 2 * FT(π) * K_therm / L_f
    fac = fac₀ * ΔT
    ∂L_f_∂T = TDI.TD.Parameters.cp_l(tps) - TDI.TD.Parameters.cp_i(tps)
    ∂fac_∂T = fac₀ * (1 - ΔT * ∂L_f_∂T / L_f)
    ∫melt = _ventilation_from(velocity_params, aps, ρₐ, state, logλ, quad)
    dLdt_unclamped = fac * ∫melt

    # only consider melting (not fusion)
    dLdt = max(0, dLdt_unclamped)
    ∂dLdt_∂T = ifelse(dLdt > 0, ∂fac_∂T * ∫melt, zero(dLdt))

    # One fractional loss for all four ice slots: number here, the rime pair at the caller.
    # The number rate follows the identity `dNdt = ρn_ice * (dLdt / ρq_ice)`, with the
    # fraction bounded by the conduction-limited melt rate of a nucleation-size particle.
    (; frac, ∂frac_∂T) =
        ice_melt_fraction(aps, tps, state.params, Tₐ, ρq_ice, dLdt, ∂dLdt_∂T)
    dNdt = ρn_ice * frac
    ∂dNdt_∂T = ρn_ice * ∂frac_∂T

    return (;
        dNdt, dLdt, melt_frac = frac,
        ∂dNdt_∂T, ∂dLdt_∂T, ∂melt_frac_∂T = ∂frac_∂T,
    )
end

"""
    zero_ice_melt(ρₐ)

The [`ice_melt`](@ref) return with every rate and derivative zero, for temperatures at or
below `T_freeze`.
"""
@inline function zero_ice_melt(ρₐ)
    o = zero(ρₐ)
    return (;
        dNdt = o, dLdt = o, melt_frac = o,
        ∂dNdt_∂T = o, ∂dLdt_∂T = o, ∂melt_frac_∂T = o,
    )
end

"""
    ICE_NUCLEATION_DIAMETER(FT)

Diameter of a nascent deposition-nucleation ice crystal [m]: the small-`D` tail of the P3
distribution, and the size at which `HetIceNucleation.deposition_rate` injects new crystals.

Single source of truth for the nucleation size, so that every quantity derived from it moves
together. TODO: put into ClimaParams.
"""
@inline ICE_NUCLEATION_DIAMETER(::Type{FT}) where {FT} = FT(10e-6)

"""
    ice_nucleation_mass(p3)

Mass of one nascent deposition-nucleation ice crystal [kg], `ρ_i (π/6) D_nuc³` with `D_nuc` from
[`ICE_NUCLEATION_DIAMETER`](@ref). This is the smallest particle mass the scheme can create.
"""
@inline ice_nucleation_mass(p3) =
    p3.ρ_i * CO.volume_sphere_D(ICE_NUCLEATION_DIAMETER(typeof(p3.ρ_i)))

"""
    ice_mean_particle_mass_min(p3)
    ice_mean_particle_mass_max(FT)

Bounds of the physical mean ice particle mass range [kg], shared by the ice number adjustment and
by the melt number rate.

The lower bound is [`ice_nucleation_mass`](@ref) rather than a literal of its own, so a numerical
guard cannot end up above the size the scheme nucleates at. It previously read `1e-12` kg,
annotated "~10 μm crystal" but in fact a 12.77 μm solid-ice sphere and 2.08x the nucleation mass
`4.7998e-13` kg, which made the guard active on freshly nucleated populations: reading
`n > q / x_min`, the adjustment relaxed the number toward 48% of what nucleation had just supplied
while conserving the mass. Derived, a population of fresh crystals sits exactly ON the bound, where
`clamp` is inert.

The upper bound is a REGULARIZATION TARGET rather than a physical ceiling: it is the mean mass the
number adjustment relaxes an over-massive population toward, and real populations exceed it. It is
raised here from `1e-5` kg, which is a solid-ice sphere of 2.75 mm, to `1e-4` kg, which is 5.93 mm,
because graupel and small hail routinely carry a mean mass above the smaller value and the
adjustment was therefore acting on ordinary ice rather than on degenerate states. Two independent
records agree on the size of the effect: a gen-2 production pin census over `219,086` populated
cells found mean masses to `0.0011` kg, over a hundred times the old target, and on a 48 km
restart `6.4` percent of populated cells sat above `1e-5` kg against `0.1` percent above `1e-4`,
so the raise removes the adjustment from about ninety-eight percent of the cells it was reaching.
Nothing else reads this bound: the shape solver's own bracket is derived from the LOWER bound
(see [`_derived_logλ_bracket`](@ref), whose docstring records why the two edges are different
kinds of quantity), so this value moves the number adjustment and nothing else.
"""
@inline ice_mean_particle_mass_min(p3) = ice_nucleation_mass(p3)
@inline ice_mean_particle_mass_max(::Type{FT}) where {FT} = FT(1e-4)

"""
    ice_moments_are_admissible(p3, ρq_ice, ρn_ice)

Whether the ice mass and ice number describe a population the scheme can represent: both
positive, and the mean particle mass they imply inside
`[ice_mean_particle_mass_min, ice_mean_particle_mass_max]`.

The single source of the admissibility question, so that a host masking its prognostic state and
[`admissible_ice_moments`](@ref) repairing it cannot come to different answers. Written as two
products rather than as a quotient, so no division is taken and no guard on a vanishing number is
needed; the one state it therefore misreads is a mean mass below the floor in a cell whose mass is
itself below the smallest normal float, more than twenty orders under a single micron-sized
crystal.
"""
@inline function ice_moments_are_admissible(p3, ρq_ice, ρn_ice)
    FT = UT.promote_typeof(ρq_ice, ρn_ice)
    o = zero(FT)
    x_min = FT(ice_mean_particle_mass_min(p3))
    x_max = ice_mean_particle_mass_max(FT)
    return (ρq_ice > o) & (ρn_ice > o) &
           (ρn_ice * x_min <= ρq_ice) & (ρq_ice <= ρn_ice * x_max)
end

"""
    admissible_ice_moments(p3, ρq_ice, ρn_ice, ρq_rim, ρb_rim)

The four P3 ice moments projected onto the admissible set, returned in the order they were given.

A quartet is admissible when it is entirely zero, or when the mass and the number are both
positive and the mean particle mass they imply lies inside the physical range
`[ice_mean_particle_mass_min, ice_mean_particle_mass_max]`. Where it is not, all four moments are
set to zero: the mass has a reservoir and the number does not, so evacuating the category is the
only repair that invents nothing. In a host carrying total water and total energy the ice becomes
vapour through the `q_tot` identity and the sublimation enthalpy leaves the sensible heat by
itself, so a cell repaired this way cools by the right amount with no bookkeeping of its own.

The rime pair is a SUB-PARTITION and is never the trigger: a negative rime mass says the rime
partition is wrong, not that the ice is not there, so it is projected the way
[`state_from_prognostic`](@ref) projects it rather than deleting the ice that carries it. On a
48 km record the two are not close: the pair test fires on cells holding `0.008` percent of the
domain ice, while evacuating on a negative rime mass instead would have taken `0.125` percent.

**This is a HOST-side repair and is deliberately not what the kernel does with the same state.**
`p3_2m_process_rates` clamps every moment it reads and then treats positive mass with no number as
ORPHAN ICE, draining it as mass at the rate a population of nucleation-mass crystals would and only
where the air is subsaturated with respect to ice. That gate is right, and the substep march can
still produce an orphan mid-step, so the drain keeps its work. What the kernel cannot do is
guarantee the invariant to everything else that reads the prognostic state: the diagnostics, the
presence masks and transport all read `Y.c` directly, which is why the ice-number diagnostic
carries a clamp of its own "for display only". This function is what lets that clamp be the
scheme's contract instead of each reader's local patch.

The admissibility question itself is [`ice_moments_are_admissible`](@ref), so a host that masks
its own prognostic state with the predicate and a caller that repairs a quartet with this function
cannot disagree about WHICH STATES ARE ADMISSIBLE. That is a statement about classification and
not about arithmetic: the two agree on every finite state, and on a non-finite one this function
leaves the moment non-finite rather than repairing it, because a non-finite moment reports a
numerical failure rather than an inadmissible physical state and the host's own whole-state check
is what should see it.

Branchless, so it costs the same on every lane of a warp.
"""
@inline function admissible_ice_moments(p3, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    FT = UT.promote_typeof(ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    o = zero(FT)
    # A MULTIPLIER rather than a select, so a non-finite moment stays non-finite: `NaN * 0` is
    # `NaN` while `ifelse` would write a real zero over it. A non-finite moment reports a numerical
    # failure rather than an inadmissible physical state, and this function projects a physical
    # state; laundering the first into a plausible zero destroys the only evidence that something
    # upstream broke, and defeats the host's own whole-state non-finite check. Every other
    # sanitizing step on this surface already propagates it: `clamp_to_nonneg` is `max(zero, x)`,
    # which is what `state_from_prognostic` and the kernel entry apply.
    mask = ifelse(ice_moments_are_admissible(p3, ρq_ice, ρn_ice), one(FT), o)
    q = FT(ρq_ice) * mask
    n = FT(ρn_ice) * mask
    # The rime mass cannot exceed the ice mass that carries it, and cannot be negative; the rime
    # volume then follows its pair onto the density cone, which returns zero of its own accord
    # once the rime mass is zero.
    qr = clamp(FT(ρq_rim) * mask, o, q)
    (ρ_rim_min, ρ_rim_max) = rime_density_bounds(p3)
    br = UT.nearest_admissible_b(qr, FT(ρb_rim), ρ_rim_min, ρ_rim_max)
    return (q, n, qr, br)
end

"""
    ice_melt_fraction_limit(aps, tps, p3, Tₐ)

Upper bound on the fractional ice melt rate [1/s] and its exact temperature derivative, as
`(; inv_τ, ∂inv_τ_∂T)`: the conduction-limited melt rate of a solid-ice sphere at the
nucleation size,

```math
1/τ = \\frac{3 K_{therm} (Tₐ - T_{freeze})}{ρ_i r_{min}^2 L_f(Tₐ)},
\\qquad r_{min} = D_{nuc} / 2,
```

zero at and below `T_freeze`. A fractional melt rate above this bound removes latent heat
faster than conduction supplies it to the smallest particle the scheme creates by nucleation,
so the bound is a property of the air and of [`ICE_NUCLEATION_DIAMETER`](@ref).
"""
@inline function ice_melt_fraction_limit(aps::CMP.AirProperties, tps::TDI.PS, p3, Tₐ)
    (; K_therm) = aps
    L_f = TDI.Lf(tps, Tₐ)
    FT = typeof(p3.ρ_i)
    r_min = ICE_NUCLEATION_DIAMETER(FT) / 2
    ΔT = max(Tₐ - p3.T_freeze, 0)
    fac = 3 * K_therm / (p3.ρ_i * r_min^2 * L_f)
    inv_τ = fac * ΔT
    ∂L_f_∂T = TDI.TD.Parameters.cp_l(tps) - TDI.TD.Parameters.cp_i(tps)
    ∂inv_τ_∂T = ifelse(ΔT > 0, fac * (1 - ΔT * ∂L_f_∂T / L_f), zero(inv_τ))
    return (; inv_τ, ∂inv_τ_∂T)
end

"""
    ice_melt_fraction(aps, tps, p3, Tₐ, ρq_ice, dLdt, ∂dLdt_∂T = zero(dLdt))

Fractional ice melt rate [1/s] and its temperature derivative, `(; frac, ∂frac_∂T)`:
`dLdt / ρq_ice` bounded above by [`ice_melt_fraction_limit`](@ref), and zero where the ice
mass is absent. The quotient alone is unbounded as `ρq_ice` vanishes while the melt integral
stays positive; the bound keeps the fraction at what conduction can melt of the smallest
particle. Shared by the melt number rate and the rime drain, so the number and both rime
moments leave at one fraction.

The value is bounded on both branches. Its derivative with respect to `ρq_ice` on the
quotient branch is `-dLdt / ρq_ice²`, which grows without bound as the mass vanishes and can
overflow Float32 at trace mass even where the value sits below the bound; a consumer
differentiating through this function there relies on its own non-finite handling.
"""
@inline function ice_melt_fraction(
    aps::CMP.AirProperties, tps::TDI.PS, p3, Tₐ, ρq_ice, dLdt, ∂dLdt_∂T = zero(dLdt),
)
    lim = ice_melt_fraction_limit(aps, tps, p3, Tₐ)
    quot = UT.guarded_quotient(dLdt, ρq_ice)
    ∂quot_∂T = UT.guarded_quotient(∂dLdt_∂T, ρq_ice)
    limited = quot > lim.inv_τ
    # `oftype` keeps the selection type stable under automatic differentiation: the limit is
    # independent of the species state, so its state partials are zero, which `convert`
    # supplies.
    frac = ifelse(limited, oftype(quot, lim.inv_τ), quot)
    ∂frac_∂T = ifelse(limited, oftype(∂quot_∂T, lim.∂inv_τ_∂T), ∂quot_∂T)
    return (; frac, ∂frac_∂T)
end

"""
    ice_deposition_timescale(velocity_params, aps, tps, Tₐ, ρₐ, state, logλ; quad)

Compute the vapor deposition relaxation timescale of the ice population from
its capacitance integral,

```math
τ_{dep} = \\frac{ρₐ q_{v,si}}{2π G_i ∫ D F_v(D) N'(D) dD},
```

with spherical capacitance `C = D/2`, following Morrison and Milbrandt (2015).
The timescale diverges as the population vanishes and shrinks as the
integrated particle surface grows.

# Arguments
 - `velocity_params`: [`CMP.Chen2022VelType`](@ref)
 - `aps`: [`CMP.AirProperties`](@ref)
 - `tps`: thermodynamics parameters
 - `Tₐ`: temperature (K)
 - `ρₐ`: air density
 - `state`: a [`P3State`](@ref) object
 - `logλ`: the log of the slope parameter [log(1/m)]

# Keyword arguments
 - `quad`: quadrature rule (a `Quadrature.QuadratureRule`)

# Returns
- Deposition timescale [s], bounded above at [`ICE_DEP_TIMESCALE_MAX`](@ref) to stay finite.
  **A returned value equal to that bound means the capacitance integral underflowed, i.e. the
  population cannot support the relaxation at all.** It does NOT mean "a very slow but real
  relaxation": the unbounded quotient diverges there, and callers that form `deficit / τ` will
  otherwise mint condensate at `deficit / ICE_DEP_TIMESCALE_MAX` on a state with no particles.
  Test with [`ice_deposition_is_degenerate`](@ref) rather than comparing to a literal.
"""
@inline function ice_deposition_timescale(
    velocity_params, aps::CMP.AirProperties, tps::TDI.PS,
    Tₐ, ρₐ, state::P3State, logλ;
    quad,
)
    FT = eltype(state)

    G = CO.G_func_ice(aps, tps, Tₐ)
    qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, Tₐ, ρₐ)

    ∫DFvN = _ventilation_from(velocity_params, aps, ρₐ, state, logλ, quad)

    denom = 2 * FT(π) * G * ∫DFvN
    return min(ρₐ * qᵥ_sat_ice / max(denom, floatmin(FT)), ICE_DEP_TIMESCALE_MAX(FT))
end

"""
    ICE_DEP_TIMESCALE_MAX(FT)

Upper bound on the ice deposition/sublimation relaxation timescale [s], applied by
[`ice_deposition_timescale`](@ref) to keep it finite as the ice population vanishes.

Reaching this bound is a **degenerate** outcome, not a slow one. The unbounded timescale is
`ρₐ q_{v,sat,ice} / (2π G ∫ D F_v N′ dD)`, and the integral goes to zero with the population, so the
quotient diverges. Capping it makes the value finite but leaves `deficit / τ` finite too, which is a
mass source on a state that has no particles to deposit onto. Use
[`ice_deposition_is_degenerate`](@ref) to detect it and zero the rate.
"""
@inline ICE_DEP_TIMESCALE_MAX(::Type{FT}) where {FT} = FT(1e10)

"""
    ice_deposition_is_degenerate(τ_dep)

`true` when `τ_dep` from [`ice_deposition_timescale`](@ref) sits at
[`ICE_DEP_TIMESCALE_MAX`](@ref), meaning the capacitance integral underflowed and there is no
population to deposit onto or sublimate from, so the physical rate is exactly zero.

Named rather than written inline so the bound is compared against its own definition at every call
site instead of a repeated literal.
"""
@inline ice_deposition_is_degenerate(τ_dep::FT) where {FT} =
    τ_dep >= ICE_DEP_TIMESCALE_MAX(FT)

"""
    ice_population_is_present(state)

`true` when the ice population can support the mixed-phase process rates that integrate over it -
riming and the other liquid-ice collisions, aggregation, and melting - i.e. when the number
moment is strictly positive and the mass moment carries at least one nucleated crystal's worth of
mass for that number:

    (ρq_ice / ice_nucleation_mass > ρn_ice) & (ρn_ice > 0)

This is a presence test on the two moments, not a smallness threshold on the mass. A fixed
mixing-ratio threshold such as `ϵ_numerics_2M_M(FT) = eps(FT)` sits at 1.1920929e-7 kg/kg at
Float32 and at 2.2e-16 at Float64, so it switches melting, riming and aggregation off
discontinuously at a physical loading at one precision and at an unreachable one at the other -
precision dependence that is not rounding. The scaled mass test keeps precision independence: it
moves with `ρn_ice` and with [`ice_nucleation_mass`](@ref) rather than fixing a mixing-ratio
threshold, so a trace population passes as soon as its crystals carry their birth mass, at either
precision. The mass side is formed as a quotient because the product
`ρn_ice * ice_nucleation_mass` underflows Float32 at trace number, where the test would degrade
to `ρq_ice > 0`; the quotient stays normal over the full Float32 range of `ρq_ice`.

Both moments enter because all three rates are integrals of the ice size distribution
`N′(D) ∝ ρn_ice`. With `ρn_ice = 0` the distribution is identically zero, so the rates are zero by
construction - and the `ρn_ice > 0` conjunct is what keeps mass-without-number states absent,
where the mass test alone reads positive at any `ρq_ice`; that mass is the ice orphan drain's to
remove, not the integrals' to melt or rime. With less than one crystal's mass per particle the
mean particle mass sits below the smallest particle the scheme can create, the shape solve has no
physical target, and the number adjustment already relaxes the number toward
`ρq_ice / ice_mean_particle_mass_min`; a population of fresh crystals sits on the boundary, where
every process this predicate controls is well posed on either side.
"""
@inline ice_population_is_present(state::P3State) =
    (state.ρq_ice / ice_nucleation_mass(state.params) > state.ρn_ice) &
    (state.ρn_ice > 0)

"""
    collision_cross_section_ice_liquid_coeffs(rᵢ)
    collision_cross_section_ice_liquid_coeffs(state, Dᵢ)

Monomial coefficients `(k₀, k₁, k₂)` of the ice-liquid collision cross-section as
a polynomial in the liquid diameter `Dₗ`,

```math
σ(Dᵢ, Dₗ) = π (rᵢ + Dₗ/2)² = k₀ + k₁ Dₗ + k₂ Dₗ²,
```

with `k₀ = π rᵢ²`, `k₁ = π rᵢ`, `k₂ = π/4`, where the ice effective radius
is `rᵢ = √(ice_area(state, Dᵢ)/π)`; see [`ice_area`](@ref).

Used in [`collision_cross_section_ice_liquid`](@ref)
"""
@inline collision_cross_section_ice_liquid_coeffs(rᵢ::FT) where {FT} =
    (π * rᵢ^2, π * rᵢ, FT(π / 4))
@inline collision_cross_section_ice_liquid_coeffs(state, Dᵢ) =
    collision_cross_section_ice_liquid_coeffs(√(ice_area(state, Dᵢ) / π))

"""
    collision_cross_section_ice_liquid(state, Dᵢ, Dₗ)

Ice-liquid collision cross-section [m²], `π (rᵢ(Dᵢ) + Dₗ/2)²`, evaluated by
Horner from the shared [`collision_cross_section_ice_liquid_coeffs`](@ref).
"""
collision_cross_section_ice_liquid(state, Dᵢ, Dₗ) =
    evalpoly(Dₗ, collision_cross_section_ice_liquid_coeffs(state, Dᵢ))

"""
    volumetric_collision_rate_integrand(velocity_params, ρₐ, state)

Returns a function that computes the volumetric collision rate integrand for ice-liquid collisions [m³/s].
The returned function takes ice and liquid particle diameters as arguments.

# Arguments
- `velocity_params`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density
- `state`: [`P3State`](@ref)

# Returns
A function `(D_ice, D_liq) -> E * K * |vᵢ - vₗ|` where:
- `D_ice` and `D_liq` are the (maximum) diameters of the ice and liquid particles
- `E` is the collision efficiency
- `K` is the collision cross section
- `vᵢ` and `vₗ` are the terminal velocities of ice and liquid particles

Note that `E`, `K`, `vᵢ` and `vₗ` are all, in general, functions of `D_ice` and `D_liq`.

This function is a component of integrals like

```math
∫ ∫ E * K * |vᵢ - vₗ| * N'_i * N'_l dD_i dD_l
```
"""
function volumetric_collision_rate_integrand(velocity_params, ρₐ, state)
    v_ice = ice_particle_terminal_velocity(velocity_params, ρₐ, state)
    v_liq = CO.particle_terminal_velocity(velocity_params.rain, ρₐ)
    function integrand(D_ice::FT, D_liq::FT) where {FT}
        E = FT(1)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
        K = collision_cross_section_ice_liquid(state, D_ice, D_liq)
        return E * K * abs(v_ice(D_ice) - v_liq(D_liq))
    end

    return integrand
end

"""
    compute_max_freeze_rate(aps, tps, velocity_params, ρₐ, Tₐ, state)

Returns a function `max_freeze_rate(Dᵢ)` that returns the maximum possible freezing rate [kg/s]
    for an ice particle of diameter `Dᵢ` [m]. Evaluates to `0` if `T ≥ T_freeze`.

# Arguments
- `aps`: [`CMP.AirProperties`](@ref)
- `tps`: `TDP.ThermodynamicsParameters`
- `velocity_params`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]
- `Tₐ`: air temperature [K]
- `state`: [`P3State`](@ref)

This rate represents the thermodynamic upper limit to collisional freezing,
which occurs when the heat transfer from the ice particle to the environment is
balanced by the latent heat of fusion.

From Eq (A7) in Musil (1970), [Musil1970](@cite).
"""
function compute_max_freeze_rate(aps, tps, velocity_params, ρₐ, Tₐ, state)
    (; D_vapor, K_therm) = aps
    cp_l = TDI.cp_l(tps)
    T_frz = TDI.T_freeze(tps)
    Lᵥ = TDI.Lᵥ(tps, Tₐ)
    L_f = TDI.Lf(tps, Tₐ)
    Tₛ = T_frz  # the surface of the ice particle is assumed to be at the freezing temperature
    ΔT = Tₛ - Tₐ  # temperature difference between the surface of the ice particle and the air
    Δρᵥ_sat =
        ρₐ * (  # saturation vapor density difference between the surface of the ice particle and the air
            TDI.p2q(tps, Tₛ, ρₐ, TDI.saturation_vapor_pressure_over_ice(tps, Tₛ)) -
            TDI.p2q(tps, Tₐ, ρₐ, TDI.saturation_vapor_pressure_over_ice(tps, Tₐ))
        )
    v_term = ice_particle_terminal_velocity(velocity_params, ρₐ, state)
    F_v = CO.ventilation_factor(state.params.vent, aps, v_term)
    # Musil (1970) dry-growth formula: the denominator `(L_f - cp_l·ΔT)`
    # represents the *net* latent heat per unit mass available to freeze a
    # colliding droplet. At Tₐ ≲ 220 K (ΔT ≳ L_f/cp_l ≈ 53 K with
    # T-dependent L_f, see Eq. A7 in Musil 1970), the denominator flips
    # sign, making `max_freeze_rate < 0` — which is unphysical. Cold air
    # is *further from* the dry/wet-growth transition, not closer to it:
    # the physical answer is `f_frz → 1` (every colliding droplet
    # freezes). We enforce that by returning `floatmax(FT)` when the
    # denominator is non-positive, so `min(∂ₜM_col, ∂ₜM_max) = ∂ₜM_col` and
    # `f_frz = 1`.
    denom = L_f - cp_l * ΔT
    function max_freeze_rate(Dᵢ)
        # fallback values typed by the promotion of the node and the captured state
        # (mixed plain/Dual under differentiation)
        FT = UT.promote_typeof(Dᵢ, ΔT, Δρᵥ_sat, denom)
        Tₐ ≥ T_frz && return zero(FT)  # No collisional freezing above the freezing temperature
        denom > 0 || return floatmax(FT)
        return 2 * (π * Dᵢ) * F_v(Dᵢ) * (K_therm * ΔT + Lᵥ * D_vapor * Δρᵥ_sat) / denom
    end
    return max_freeze_rate
end

"""
    compute_local_rime_density(velocity_params, ρₐ, T, state)

Provides a function `ρ′_rim(Dᵢ, Dₗ)` that computes the local rime density [kg/m³]
    for a given ice particle diameter `Dᵢ` [m] and liquid particle diameter `Dₗ` [m].

# Arguments
- `velocity_params`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]
- `T`: temperature [K]
- `state`: [`P3State`](@ref)

# Returns
A function that computes the local rime density [kg/m³] using the equation:

```math
ρ'_{rim} = a + b R_i + c R_i^2
```
where
```math
R_i = \\frac{ 10^6 ⋅ D_{liq} ⋅ |v_{liq} - v_{ice}| }{ 2 T_{sfc} }
```
and ``T_{sfc}`` is the surface temperature [°C], ``D_{liq}`` is the liquid particle
diameter [m], ``v_{liq/ice}`` is the particle terminal velocity [m/s].
With the ``10^6`` factor converting ``D_{liq}`` from [m] to [μm], the units of
``R_i`` are [μm m s⁻¹ °C⁻¹]. The units of ``ρ'_{rim}`` are [kg/m³].

We assume for simplicity that ``T_{sfc}`` equals ``T``, the ambient air temperature.
For real graupel, ``T_{sfc}`` is slightly higher than ``T`` due to latent heat release
of freezing liquid particles onto the ice particle. Morrison & Milbrandt (2013)
found little sensitivity to "realistic" increases in ``T_{sfc}``.

See also [`LocalRimeDensity`](@ref CloudMicrophysics.Parameters.LocalRimeDensity).

# Extended help

 Implementation follows Cober and List (1993), Eq. 16 and 17.
 See also the P3 fortran code, `microphy_p3.f90`, Line 3315-3323,
 which extends the range of the calculation to ``R_i ≤ 12``, the upper limit of which
 then equals the solid bulk ice density, ``ρ_ice = 916.7 kg/m^3``.

 Note that Morrison & Milbrandt (2015) [MorrisonMilbrandt2015](@cite) only uses this
 parameterization for collisions with cloud droplets.
 For rain drops, they use a value near the solid bulk ice density, ``ρ^* = 900 kg/m^3``.
 We do not consider this distinction, and use this parameterization for all liquid particles.
"""
function compute_local_rime_density(velocity_params, ρₐ, T, state)
    (; T_freeze, ρ_rim_local) = state.params
    T°C = T - T_freeze  # Convert to °C
    μm = 1_000_000  # Note: m to μm factor, c.f. units of rₘ in Eq. 16 in Cober and List (1993)

    v_ice = ice_particle_terminal_velocity(velocity_params, ρₐ, state)
    v_liq = CO.particle_terminal_velocity(velocity_params.rain, ρₐ)
    function ρ′_rim(Dᵢ, Dₗ)
        v_term = abs(v_ice(Dᵢ) - v_liq(Dₗ))
        Rᵢ = (Dₗ * μm * v_term) / (2 * T°C)  # Eq. 16 in Cober and List (1993). Note: no `-` due to absolute value in v_term
        return ρ_rim_local(Rᵢ)
    end
    return ρ′_rim
end

"""
    get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, liq_bounds; [quad])

Returns a function `liquid_integrals(Dᵢ)` that computes the liquid particle integrals
    for a given ice particle diameter `Dᵢ`.

# Arguments
- `n`: liquid particle size distribution function `n(D)`
- `∂ₜV`: volumetric collision rate integrand function `∂ₜV(Dᵢ, D)`
- `m_liq`: liquid particle mass function `m_liq(D)`
- `ρ′_rim`: local rime density function `ρ′_rim(Dᵢ, D)`
- `liq_bounds`: integration bounds for liquid particles

# Keyword arguments
- `quad`: quadrature rule, default is `ChebyshevGauss(100)`

# Notes
The function `liquid_integrals(Dᵢ)` returns a tuple `(∂ₜN_col, ∂ₜM_col, ∂ₜB_col)`
    of collision rates at `Dᵢ`, where:
- `∂ₜN_col`: number collision rate [1/s]
- `∂ₜM_col`: mass collision rate [kg/s]
- `∂ₜB_col`: rime volume collision rate [m³/s]
"""
@inline function get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, liq_bounds; quad = ChebyshevGauss(100))
    function liquid_integrals(Dᵢ)
        integrand = D -> begin
            V_val = ∂ₜV(Dᵢ, D)
            n_val = n(D)
            m_val = m_liq(D)
            term1 = V_val * n_val
            term2 = term1 * m_val
            term3 = term2 / ρ′_rim(Dᵢ, D)
            return SA.SVector((term1, term2, term3))
        end
        (∂ₜN_col, ∂ₜM_col, ∂ₜB_col) = integrate(integrand, liq_bounds, quad)
        return ∂ₜN_col, ∂ₜM_col, ∂ₜB_col
    end
    return liquid_integrals
end

"""
    crossover_diameter(v_target, v_l, D_min, D_max)

Find the diameter `D` in `[D_min, D_max]` where `v_l(D) = v_target`
"""
function crossover_diameter(v_target, v_l::F, D_min, D_max) where {F}
    FT = float(promote_type(typeof(v_target), typeof(D_min), typeof(D_max)))
    f(D) = v_l(D) - v_target
    maxiters = FT === Float32 ? 8 : 10
    sol = RS.find_zero(f,
        RS.BrentsMethod(FT(D_min), FT(D_max)), RS.CompactSolution(),
        FixedIterations{FT}(), maxiters,
    )
    return sol.root
end

"""
    closed_rain_inner_NM(
        Dᵢ, v_i_at_Dᵢ, v_l, rᵢ, ρw, ai, bi, ci, D_min, D_max, N₀r, Dr_mean,
    )

Closed-form `(∂ₜN_col, ∂ₜM_col)` for the rain inner integral at one outer `Dᵢ`.
"""
function closed_rain_inner_NM(Dᵢ, v_i_at_Dᵢ, v_l::F, rᵢ, ρw, ai, bi, ci, D_min, D_max, N₀r, Dr_mean) where {F}
    FT = float(eltype(ai))
    λ = inv(Dr_mean)  # rain PSD slope: n_r(D) ∝ e^{-λ D}
    Dstar = crossover_diameter(v_i_at_Dᵢ, v_l, D_min, D_max)

    # Compute rain PSD incomplete moments weighted by ice-liquid collision
    # cross-section `K`, and sedimentation velocity difference `|vᵢ - vₗ|`
    coeffs = SA.SVector(collision_cross_section_ice_liquid_coeffs(rᵢ))
    function Iᵖ(a, b, p, α)
        acc = @inbounds coeffs[1] * gamma_inc_moment(a, b, p, α)
        @inbounds for i in 2:lastindex(coeffs)
            acc += coeffs[i] * gamma_inc_moment(a, b, p + (i - 1), α)
        end
        return acc
    end
    function flux(a, b, p)  # ≡ ∫ₐᵇ K(Dᵢ, Dₗ) ⋅ (vᵢ(Dᵢ) - vₗ(Dₗ)) ⋅ n_r(Dₗ) dDₗ
        s = v_i_at_Dᵢ * Iᵖ(a, b, p, λ)  # vᵢ ⋅ ∫ₐᵇ K ⋅ n_r dDₗ
        @inbounds for j in eachindex(ai)  # - ∫ₐᵇ K ⋅ vₗ ⋅ n_r dDₗ
            s -= ai[j] * Iᵖ(a, b, p + bi[j], λ + ci[j])
        end
        return s
    end
    crossing(p) = flux(D_min, Dstar, p) - flux(Dstar, D_max, p)  # sign flip at Dstar
    mfac = ρw * CO.volume_sphere_D(one(FT))  # m_liq(D) = mfac Dₗ³
    return (N₀r * crossing(FT(0)), N₀r * mfac * crossing(FT(3)))  # number: D⁰, mass: D³
end

"""
    get_liquid_integrals_rain_closed(
        psd_r::RainParticlePDF_SB2006, vel::Chen2022VelType,
        n_r, ρₐ, L_r, N_r, state, ∂ₜV, m_liq, ρ′_rim, bounds_r; quad
    )

Returns a function `liquid_integrals(Dᵢ) -> (∂ₜN_col, ∂ₜM_col, ∂ₜB_col)` 
where N and M are the exact incomplete-gamma closed form and 
B_rim is computed by quadrature
"""
@inline function get_liquid_integrals_rain_closed(
    psd_r::CMP.RainParticlePDF_SB2006, vel::CMP.Chen2022VelType,
    n_r, ρₐ, L_r, N_r, state, ∂ₜV, m_liq, ρ′_rim, bounds_r; quad,
)
    FT = promote_type(eltype(state), UT.promote_typeof(ρₐ, L_r, N_r))
    ρw = psd_r.ρw
    (; N₀r, Dr_mean) = CM2.pdf_rain_parameters(psd_r, L_r / ρₐ, ρₐ, N_r)
    ai_t, bi_t, ci_t = CO.Chen2022_vel_coeffs(vel.rain, ρₐ)
    ai, bi, ci = SA.SVector(ai_t), SA.SVector(bi_t), SA.SVector(ci_t)
    v_l = CO.particle_terminal_velocity(vel.rain, ρₐ)
    v_i = ice_particle_terminal_velocity(vel, ρₐ, state)
    D_min, D_max = bounds_r
    zero_rates = (zero(FT), zero(FT), zero(FT))
    function liquid_integrals(Dᵢ)
        if iszero(N₀r) || !(D_max > D_min)
            return zero_rates
        end
        v_i_at_Dᵢ = v_i(Dᵢ)
        rᵢ = sqrt(ice_area(state, Dᵢ) / π)
        ∂ₜN_col, ∂ₜM_col = closed_rain_inner_NM(
            Dᵢ, v_i_at_Dᵢ, v_l, rᵢ, ρw, ai, bi, ci,
            D_min, D_max, N₀r, Dr_mean,
        )
        if !(isfinite(∂ₜN_col) && isfinite(∂ₜM_col))
            return zero_rates
        end
        ∂ₜB_col = integrate(
            D -> ∂ₜV(Dᵢ, D) * n_r(D) * m_liq(D) / ρ′_rim(Dᵢ, D),
            bounds_r,
            quad,
        )
        return (∂ₜN_col, ∂ₜM_col, ∂ₜB_col)
    end
    return liquid_integrals
end

@inline _rain_inner_integrals(
    psd_r::CMP.RainParticlePDF_SB2006, vel::CMP.Chen2022VelType,
    n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r, ρₐ, L_r, N_r, state; quad,
) = get_liquid_integrals_rain_closed(
    psd_r, vel, n_r, ρₐ, L_r, N_r, state, ∂ₜV, m_liq, ρ′_rim, bounds_r;
    quad,
)
@inline _rain_inner_integrals(
    ::Any, ::Any,
    n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r, ρₐ, L_r, N_r, state; quad,
) = get_liquid_integrals(n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r; quad)

"""
    ∫liquid_ice_collisions(
        n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; [quad]
    )

Computes the bulk collision rate integrands between ice and liquid particles.

# Arguments
- `n_i`: ice particle size distribution function n_i(D)
- `∂ₜM_max`: maximum freezing rate function ∂ₜM_max(Dᵢ)
- `cloud_integrals`: inner liquid integrals for cloud particles, e.g. from [`get_liquid_integrals`](@ref)
- `rain_integrals`: inner liquid integrals for rain particles, e.g. from [`get_liquid_integrals`](@ref) or `get_liquid_integrals_rain_closed`
- `ice_bounds`: integration bounds for ice particles, from [`integral_bounds`](@ref)

# Keyword arguments
- `quad`: quadrature rule, default is `ChebyshevGauss(100)`

# Returns
A 10-element vector of integrated rates,
`(QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL, ∫𝟙_wet_M_col)`;
see the method below for the definition of each entry.
"""
@inline function ∫liquid_ice_collisions(
    n_i,
    ∂ₜM_max,
    cloud_integrals,
    rain_integrals,
    ice_bounds;
    # Chebyshev–Gauss (not the Gauss–Legendre default used elsewhere):
    # this integrand has √-type endpoint/weight behavior that Gauss–Legendre resolves
    # poorly, so a higher node count is used here. TODO: run the same convergence check
    # used for the Gauss-Legendre default and lower this to the smallest node count within
    # tolerance.
    quad = ChebyshevGauss(100),
)
    function liquid_ice_collisions_integrands(Dᵢ)
        # Inner integrals over liquid particle diameters
        ∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col = cloud_integrals(Dᵢ)
        ∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col = rain_integrals(Dᵢ)

        # Partition the mass collisions between freezing and shedding
        ∂ₜM_col = ∂ₜM_c_col + ∂ₜM_r_col  # [kg / s]

        ∂ₜM_frz = min(∂ₜM_col, ∂ₜM_max(Dᵢ))
        f_frz = iszero(∂ₜM_col) ? zero(∂ₜM_frz) : ∂ₜM_frz / ∂ₜM_col
        𝟙_wet = ∂ₜM_col > ∂ₜM_frz  # Used for wet densification

        n = n_i(Dᵢ)
        # Integrating over `Dᵢ` gives another unit of `[m]`, so `[X / s / m]` --> `[X / s]`
        # ∂ₜX = ∫ ∂ₜX(Dᵢ) nᵢ(Dᵢ) dDᵢ
        return SA.SVector((
            n * ∂ₜM_c_col * f_frz,        # QCFRZ
            n * ∂ₜM_c_col * (1 - f_frz),  # QCSHD
            n * ∂ₜN_c_col,                # NCCOL
            n * ∂ₜM_r_col * f_frz,        # QRFRZ
            n * ∂ₜM_r_col * (1 - f_frz),  # QRSHD
            n * ∂ₜN_r_col,                # NRCOL
            n * ∂ₜM_col,                  # ∫M_col,      total collision rate
            n * ∂ₜB_c_col * f_frz,        # BCCOL,       ∂ₜB_rim source
            n * ∂ₜB_r_col * f_frz,        # BRCOL,       ∂ₜB_rim source
            n * 𝟙_wet * ∂ₜM_col,          # ∫𝟙_wet_M_col, wet growth indicator
        ))
    end
    return integrate(liquid_ice_collisions_integrands, ice_bounds, quad)
end

"""
    ∫liquid_ice_collisions(
        state, logλ, psd_c, psd_r, L_c, N_c, L_r, N_r,
        aps, tps, vel, ρₐ, T, m_liq; quad,
    )

Compute key liquid-ice collision rates and quantities. Used by [`bulk_liquid_ice_collision_sources`](@ref).

# Arguments
- `state`: [`P3State`](@ref)
- `logλ`: the log of the slope parameter [log(1/m)]
- `psd_c`: [`CMP.CloudParticlePDF_SB2006`](@ref)
- `psd_r`: [`CMP.RainParticlePDF_SB2006`](@ref)
- `L_c`: cloud liquid water content [kg/m³]
- `N_c`: cloud liquid water number concentration [1/m³]
- `L_r`: rain water content [kg/m³]
- `N_r`: rain number concentration [1/m³]
- `aps`: [`CMP.AirProperties`](@ref)
- `tps`: `TDP.ThermodynamicsParameters`
- `vel`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]
- `T`: temperature [K]
- `m_liq`: liquid particle mass function `m_liq(D)`

# Keyword arguments
- `quad`: A `QuadratureRule` instance (required)

# Returns
A tuple `(QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL, ∫𝟙_wet_M_col)`, where:
1. `QCFRZ` - Cloud mass collision rate due to freezing [kg/s]
2. `QCSHD` - Cloud mass collision rate due to shedding [kg/s]
3. `NCCOL` - Cloud number collision rate [1/s]
4. `QRFRZ` - Rain mass collision rate due to freezing [kg/s]
5. `QRSHD` - Rain mass collision rate due to shedding [kg/s]
6. `NRCOL` - Rain number collision rate [1/s]
7. `∫M_col` - Total collision rate [kg/s]
8. `BCCOL` - Cloud rime volume source [m³/m³/s]
9. `BRCOL` - Rain rime volume source [m³/m³/s]
10. `∫𝟙_wet_M_col` - Wet growth indicator [kg/s]
"""
@inline function ∫liquid_ice_collisions(
    state, logλ,
    psd_c, psd_r, L_c, N_c, L_r, N_r,
    aps, tps, vel, ρₐ, T, m_liq; quad,
)
    FT = eltype(state)

    # Particle size distributions
    n_c = DT.size_distribution(psd_c, L_c / ρₐ, ρₐ, N_c)  # n_c(Dₗ)
    n_r = DT.size_distribution(psd_r, L_r / ρₐ, ρₐ, N_r)  # n_r(Dₗ)
    n_i = DT.size_distribution(state, logλ)               # n_i(Dᵢ)

    # Initialize integration buffers by evaluating a representative integral
    p = FT(0.00001)
    ice_bounds = integral_bounds(state, logλ; p)
    bounds_c = CM2.get_size_distribution_bounds(psd_c, L_c / ρₐ, ρₐ, N_c, p)
    bounds_r = CM2.get_size_distribution_bounds(psd_r, L_r / ρₐ, ρₐ, N_r, p)

    # Integrand components
    # NOTE: We assume collision efficiency, shape (spherical), and terminal velocity is the
    #   same for cloud and precipitating liquid particles ⟹ same volumetric collision rate, ∂ₜV
    ∂ₜV = volumetric_collision_rate_integrand(vel, ρₐ, state)  # ∂ₜV(Dᵢ, Dₗ)
    ρ′_rim = compute_local_rime_density(vel, ρₐ, T, state)  # ρ′_rim(Dᵢ, Dₗ)
    ∂ₜM_max = compute_max_freeze_rate(aps, tps, vel, ρₐ, T, state)  # ∂ₜM_max(Dᵢ)

    cloud_integrals = get_liquid_integrals(n_c, ∂ₜV, m_liq, ρ′_rim, bounds_c; quad)  # (∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col)
    # Rain inner: exact closed form for the (SB2006-exp PSD, Chen-2022) pair
    # Numerical fallback for any other PSD/velocity type.
    rain_integrals = _rain_inner_integrals(
        psd_r, vel, n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r,
        ρₐ, L_r, N_r, state; quad,
    )  # (∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col)

    return ∫liquid_ice_collisions(n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; quad)
end

"""
    bulk_liquid_ice_collision_sources(
        state, logλ,
        psd_c, psd_r, L_c, N_c, L_r, N_r,
        aps, tps, vel, ρₐ, T; quad,
    )

Computes the bulk rates for ice and liquid particle collisions.

# Arguments
- `state`: [`P3State`](@ref)
- `logλ`: the log of the slope parameter [log(1/m)]
- `psd_c`: a [`CMP.CloudParticlePDF_SB2006`](@ref)
- `psd_r`: a [`CMP.RainParticlePDF_SB2006`](@ref)
- `L_c`: cloud liquid water content [kg/m³]
- `N_c`: cloud liquid water number concentration [1/m³]
- `L_r`: rain water content [kg/m³]
- `N_r`: rain number concentration [1/m³]
- `aps`: [`CMP.AirProperties`](@ref)
- `tps`: thermodynamics parameters
- `vel`: the velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]
- `T`: temperature [K]

# Keyword arguments
- `quad`: quadrature rule, default is `ChebyshevGauss(100)`

# Returns
A `NamedTuple` of `(; ∂ₜq_c, ∂ₜq_r, ∂ₜN_c, ∂ₜN_r, ∂ₜL_rim, ∂ₜL_ice, ∂ₜB_rim)`, where:
1. `∂ₜq_c`: cloud liquid water content tendency [kg/kg/s]
2. `∂ₜq_r`: rain water content tendency [kg/kg/s]
3. `∂ₜN_c`: cloud number concentration tendency [1/m³/s]
4. `∂ₜN_r`: rain number concentration tendency [1/m³/s]
5. `∂ₜL_rim`: riming mass tendency [kg/m³/s]
6. `∂ₜL_ice`: ice water content tendency [kg/m³/s]
7. `∂ₜB_rim`: rime volume tendency [m³/m³/s]
"""
@inline function bulk_liquid_ice_collision_sources(
    state, logλ,
    psd_c, psd_r, L_c, N_c, L_r, N_r,
    aps, tps, vel, ρₐ, T; quad = ChebyshevGauss(100),
)
    FT = promote_type(eltype(state), UT.promote_typeof(L_c, N_c, L_r, N_r, ρₐ, T))
    (; τ_wet, ρ_i) = state.params
    D_shd = FT(1e-3) # 1mm  # TODO: Externalize this parameter

    ρw = psd_c.ρw
    @assert ρw == psd_r.ρw "Cloud and rain should have the same liquid water density"
    m_liq(Dₗ) = ρw * CO.volume_sphere_D(Dₗ)

    rates = ∫liquid_ice_collisions(
        state, logλ,
        psd_c, psd_r, L_c, N_c, L_r, N_r,
        aps, tps, vel, ρₐ, T, m_liq; quad,
    )
    (QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫∂ₜM_col, BCCOL, BRCOL, ∫𝟙_wet_M_col) = rates

    # Bulk wet growth fraction
    f_wet = iszero(∫∂ₜM_col) ? zero(∫∂ₜM_col) : ∫𝟙_wet_M_col / ∫∂ₜM_col

    # Shedding of rain
    # QRSHD = ∫∂ₜM_col - (QCFRZ + QRFRZ)
    NRSHD = QRSHD / m_liq(D_shd)
    # NCSHD = QCSHD / m_liq(D_shd)

    # Densification of rime
    (; ρq_ice, F_rim, ρ_rim) = state
    B_rim = iszero(ρ_rim) ? zero(ρ_rim) : (ρq_ice * F_rim) / ρ_rim  # from: ρ_rim = L_rim / B_rim
    QIWET = f_wet * ρq_ice * (1 - F_rim) / τ_wet   # densification of rime mass
    BIWET = f_wet * (ρq_ice / ρ_i - B_rim) / τ_wet  # densification of rime volume

    # Bulk rates
    ## Liquid phase
    ∂ₜq_c = (-QCFRZ - QCSHD) / ρₐ
    ∂ₜq_r = (-QRFRZ + QCSHD) / ρₐ
    ∂ₜN_c = -NCCOL
    ∂ₜN_r = -NRCOL + NRSHD
    ## Ice phase
    ∂ₜL_rim = QCFRZ + QRFRZ + QIWET
    ∂ₜL_ice = QCFRZ + QRFRZ
    # ∂ₜN_ice = 0
    ∂ₜB_rim = BCCOL + BRCOL + BIWET

    return @NamedTuple{∂ₜq_c::FT, ∂ₜq_r::FT, ∂ₜN_c::FT, ∂ₜN_r::FT, ∂ₜL_rim::FT, ∂ₜL_ice::FT, ∂ₜB_rim::FT}(
        (∂ₜq_c, ∂ₜq_r, ∂ₜN_c, ∂ₜN_r, ∂ₜL_rim, ∂ₜL_ice, ∂ₜB_rim)
    )
end

"""
    ice_self_collection(state, logλ, vel, ρₐ; [quad])

Computes the ice self-collection (aggregation) rate, which decreases the ice number concentration
while leaving mass, rime mass, and rime volume unchanged.

# Arguments
- `state`: [`P3State`](@ref)
- `logλ`: the log of the slope parameter [log(1/m)]
- `vel`: the velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]

# Keyword arguments
- `quad`: quadrature rule, default is `ChebyshevGauss(100)`

# Returns
A `NamedTuple` of `(; dNdt)`, where:
1. `dNdt`: ice number concentration tendency due to self-collection `[1/m³/s]` (always positive or zero, represents a loss rate)
"""
@inline function ice_self_collection(state, logλ, vel, ρₐ; quad = ChebyshevGauss(100))
    n_i = DT.size_distribution(state, logλ)
    v_ice = ice_particle_terminal_velocity(vel, ρₐ, state)

    p = eps(one(ρₐ))
    ice_bounds = integral_bounds(state, logλ; p)

    function inner_integral(D_1)
        v1 = v_ice(D_1)
        r1 = sqrt(ice_area(state, D_1) / π)
        # Volumetric collision rate integrand: E · K · |v₁ − v₂| · n(D₂)
        # where E = 1 (collision efficiency),
        #       K = π (r₁ + r₂)² is the geometric collision cross-section,
        #       |v₁ − v₂| is the differential sedimentation speed,
        #       r = √(A/π) is the effective radius from projected ice area A.
        integrand = D_2 -> begin
            v2 = v_ice(D_2)
            r2 = sqrt(ice_area(state, D_2) / π)
            K = π * (r1 + r2)^2                    # collision cross section
            return K * abs(v1 - v2) * n_i(D_2)     # E = 1 implied
        end
        # Split the inner integral at the |v1 - v2| cusp (D_2 = D_1, where the relative
        # fall speed vanishes and the integrand has a kink). Each half is then smooth, so
        # the quadrature converges like the other (cusp-free) P3 integrals rather than
        # being cusp-limited, letting a much lower node count reach target accuracy.
        D_lo, D_hi = first(ice_bounds), last(ice_bounds)
        rate_at_D1 = integrate(integrand, (D_lo, D_1), quad) + integrate(integrand, (D_1, D_hi), quad)
        return rate_at_D1 * n_i(D_1)
    end

    total_rate = integrate(inner_integral, ice_bounds, quad)

    # The 0.5 factor accounts for double-counting in self-collection
    FT = eltype(state)
    dNdt = FT(0.5) * total_rate
    return (; dNdt)
end
