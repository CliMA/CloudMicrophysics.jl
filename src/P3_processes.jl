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
    VolumetricCollisionRate(state, v_i, v_l)

Volumetric collision rate for ice-liquid collisions [m³/s], evaluated as
`(D_ice, D_liq) -> E * K * |vᵢ - vₗ|` where:
- `D_ice` and `D_liq` are the (maximum) diameters of the ice and liquid particles
- `E` is the collision efficiency
- `K` is the collision cross section
- `vᵢ` and `vₗ` are the terminal velocities of ice and liquid particles

Note that `E`, `K`, `vᵢ` and `vₗ` are all, in general, functions of `D_ice` and `D_liq`.

This is a component of integrals like

```math
∫ ∫ E * K * |vᵢ - vₗ| * N'_i * N'_l dD_i dD_l
```

The terminal-velocity closures `v_i` and `v_l` are fields, so consumers of the
collision rate can query the velocities and their structure directly: the
fall-speed crossing (see [`crossing_integral_bounds`](@ref)), the
[`velocity_breakpoints`](@ref), and the coefficients of a
[`CO.Chen2022VelocityCurve`](@ref) used by the closed-form rain integrals.
"""
struct VolumetricCollisionRate{S, VI, VL} <: Function
    state::S
    v_i::VI
    v_l::VL
end
@inline function (∂ₜV::VolumetricCollisionRate)(D_ice::FT, D_liq::FT) where {FT}
    E = FT(1)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
    K = collision_cross_section_ice_liquid(∂ₜV.state, D_ice, D_liq)
    return E * K * abs(∂ₜV.v_i(D_ice) - ∂ₜV.v_l(D_liq))
end

"""
    volumetric_collision_rate_integrand(velocity_params, ρₐ, state)

Construct the [`VolumetricCollisionRate`](@ref) for the Chen 2022 ice and
rain terminal velocities at air density `ρₐ`.

!!! note
    We use the same terminal velocity parameterization for cloud and rain water.
"""
volumetric_collision_rate_integrand(velocity_params, ρₐ, state) = VolumetricCollisionRate(
    state,
    ice_particle_terminal_velocity(velocity_params, ρₐ, state),
    CO.particle_terminal_velocity(velocity_params.rain, ρₐ),
)

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
    # sign, making `max_freeze_rate < 0` - which is unphysical. Cold air
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
        # `rate` is non-finite when `denom ≤ 0`; the selection below discards it
        rate = 2 * (π * Dᵢ) * F_v(Dᵢ) * (K_therm * ΔT + Lᵥ * D_vapor * Δρᵥ_sat) / denom
        # zero above the freezing temperature; floatmax when denom ≤ 0 (see above)
        return ifelse(Tₐ ≥ T_frz, zero(FT), ifelse(denom > 0, FT(rate), floatmax(FT)))
    end
    return max_freeze_rate
end

"""
    RimeDensityRate{R, FT, VI, VL}

Callable `ρ′_rim(Dᵢ, Dₗ)` returned by [`compute_local_rime_density`](@ref); computes
the local rime density [kg/m³] for a given ice particle diameter `Dᵢ` [m] and liquid
particle diameter `Dₗ` [m].
"""
struct RimeDensityRate{R, FT, VI, VL} <: Function
    ρ_rim_local::R
    T°C::FT
    v_ice::VI
    v_liq::VL
end
@inline function (ρ′::RimeDensityRate)(Dᵢ, Dₗ)
    v_term = abs(ρ′.v_ice(Dᵢ) - ρ′.v_liq(Dₗ))
    return rime_density_at(ρ′, v_term, Dₗ)
end

"""
    rime_density_at(ρ′::RimeDensityRate, v_term, Dₗ)

Local rime density [kg/m³] at a precomputed sedimentation velocity difference
`v_term = |v_ice(Dᵢ) - v_liq(Dₗ)|`, shared between the collision-rate and rime-density
evaluations at the same liquid diameter node. See [`compute_local_rime_density`](@ref).
"""
@inline function rime_density_at(ρ′::RimeDensityRate, v_term, Dₗ)
    μm = 1_000_000  # m to μm factor, c.f. units of rₘ in Eq. 16 in Cober and List (1993)
    # Leading minus: Cober and List (1993), Eq. 16, and fortran `microphy_p3.f90`
    # (`Ri = -(0.5e6*D_c)*V_impact*iTc`); T°C < 0 then makes Rᵢ positive.
    Rᵢ = -(Dₗ * μm * v_term) / (2 * ρ′.T°C)
    # At and above the melting point the wet-growth limit is selected directly, which is the
    # same value the `Rᵢ → ∞` limit of the sub-zero branch reaches.
    Rᵢ = ifelse(ρ′.T°C < 0, Rᵢ, oftype(Rᵢ, CMP.RIME_DENSITY_Rᵢ_MAX))
    return ρ′.ρ_rim_local(Rᵢ)
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
R_i = -\\frac{ 10^6 ⋅ D_{liq} ⋅ |v_{liq} - v_{ice}| }{ 2 T_{sfc} }
```
and ``T_{sfc} < 0`` is the sub-zero surface temperature [°C], ``D_{liq}`` is the liquid
particle diameter [m], ``v_{liq/ice}`` is the particle terminal velocity [m/s].
With the ``10^6`` factor converting ``D_{liq}`` from [m] to [μm], the units of
``R_i`` are [μm m s⁻¹ °C⁻¹]. The units of ``ρ'_{rim}`` are [kg/m³].

We assume for simplicity that ``T_{sfc}`` equals ``T``, the ambient air temperature.
For real graupel, ``T_{sfc}`` is slightly higher than ``T`` due to latent heat release
of freezing liquid particles onto the ice particle. Morrison & Milbrandt (2013)
found little sensitivity to "realistic" increases in ``T_{sfc}``.

See also [`LocalRimeDensity`](@ref CloudMicrophysics.Parameters.LocalRimeDensity).

# Extended help

 Implementation follows Cober and List (1993), Eq. 16 and 17, and the P3 fortran code,
 `microphy_p3.f90`. The leading minus sign in ``R_i`` matches both: with ``T_{sfc} < 0``
 the minus makes ``R_i`` positive, and the fortran carries the explicit minus in
 `Ri = -(0.5e6 D_c) V_impact iTc` (its `0.5e6 D_c` equals ``10^6 D_{liq} / 2``).

 See also the P3 fortran code, Line 3315-3323, which extends the range of the calculation
 to ``R_i ≤ 12``, the upper limit of which then equals the solid bulk ice density,
 ``ρ_ice = 916.7 kg/m^3``.
 ``R_i`` diverges as ``T_{sfc} → 0^-``, and [`LocalRimeDensity`](@ref
 CloudMicrophysics.Parameters.LocalRimeDensity) clamps it there at `RIME_DENSITY_Rᵢ_MAX`,
 returning ``ρ_ice``. That is the wet-growth limit: rime deposited at vanishing supercooling soaks and
 freezes to solid ice. At and above ``T_{sfc} = 0`` the same limit is selected directly.
 For a 20 μm droplet meeting ice at 5 m/s the clamp is reached already at ``|T_{sfc}| = 4.2``°C,
 and it binds for every ``D_{liq} |v_{liq} - v_{ice}| ≥ 2.4 ⋅ 10^{-8}`` m²/s.

 Note that Morrison & Milbrandt (2015) [MorrisonMilbrandt2015](@cite) only uses this
 parameterization for collisions with cloud droplets.
 For rain drops, they use a value near the solid bulk ice density, ``ρ^* = 900 kg/m^3``.
 We do not consider this distinction, and use this parameterization for all liquid particles.
"""
function compute_local_rime_density(velocity_params, ρₐ, T, state)
    (; T_freeze, ρ_rim_local) = state.params
    T°C = T - T_freeze
    v_ice = ice_particle_terminal_velocity(velocity_params, ρₐ, state)
    v_liq = CO.particle_terminal_velocity(velocity_params.rain, ρₐ)
    return RimeDensityRate(ρ_rim_local, T°C, v_ice, v_liq)
end

"""
    get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, liq_bounds; quad)

Return a function `liquid_integrals(Dᵢ)` that computes the liquid particle integrals
    for a given ice particle diameter `Dᵢ`.

# Arguments
- `n`: liquid particle size distribution function `n(D)`
- `∂ₜV`: the [`VolumetricCollisionRate`](@ref) `∂ₜV(Dᵢ, D)`
- `m_liq`: liquid particle mass function `m_liq(D)`
- `ρ′_rim`: local rime density function `ρ′_rim(Dᵢ, D)`
- `liq_bounds`: integration bounds for liquid particles; the fall-speed
    crossing `∂ₜV.v_l(D) = ∂ₜV.v_i(Dᵢ)` is inserted as a subinterval boundary,
    see [`crossing_integral_bounds`](@ref)

# Keyword arguments
- `quad`: quadrature rule (a `Quadrature.QuadratureRule`)

# Notes
The function `liquid_integrals(Dᵢ)` returns a tuple `(∂ₜN_col, ∂ₜM_col, ∂ₜB_col)`
    of collision rates at `Dᵢ`, where:
- `∂ₜN_col`: number collision rate [1/s]
- `∂ₜM_col`: mass collision rate [kg/s]
- `∂ₜB_col`: rime volume collision rate [m³/s]
"""
@inline function get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, liq_bounds; quad)
    function liquid_integrals(Dᵢ)
        integrand = D -> begin
            ∂ₜV_D = ∂ₜV(Dᵢ, D)
            ∂ₜN = ∂ₜV_D * n(D)          # number collision rate
            ∂ₜM = ∂ₜN * m_liq(D)        # mass collision rate
            ∂ₜB = ∂ₜM / ρ′_rim(Dᵢ, D)   # rime volume collision rate
            return SA.SVector((∂ₜN, ∂ₜM, ∂ₜB))
        end
        bnds = crossing_integral_bounds(liq_bounds, ∂ₜV, Dᵢ)
        (∂ₜN_col, ∂ₜM_col, ∂ₜB_col) = integrate(integrand, bnds, quad)
        return ∂ₜN_col, ∂ₜM_col, ∂ₜB_col
    end
    return liquid_integrals
end

# `∂ₜV.v_i(Dᵢ)` and the cross-section coefficients derived from `ice_area(state, Dᵢ)`
# depend only on the outer diameter `Dᵢ`; hoist both once per outer node instead of
# reevaluating them at every inner quadrature node. `∂ₜV.v_l(D)` is likewise shared
# between the collision-rate and rime-density evaluations at each inner node.
@inline function get_liquid_integrals(
    n, ∂ₜV::VolumetricCollisionRate, m_liq, ρ′_rim::RimeDensityRate, liq_bounds; quad,
)
    function liquid_integrals(Dᵢ)
        v_i_at_Dᵢ = ∂ₜV.v_i(Dᵢ)
        coeffs = collision_cross_section_ice_liquid_coeffs(∂ₜV.state, Dᵢ)
        integrand = D -> begin
            v_l_at_D = ∂ₜV.v_l(D)
            v_term = abs(v_i_at_Dᵢ - v_l_at_D)
            E = one(v_term)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
            ∂ₜN = E * evalpoly(D, coeffs) * v_term * n(D)  # number collision rate
            ∂ₜM = ∂ₜN * m_liq(D)                           # mass collision rate
            ∂ₜB = ∂ₜM / rime_density_at(ρ′_rim, v_term, D)  # rime volume collision rate
            return SA.SVector(∂ₜN, ∂ₜM, ∂ₜB)
        end
        bnds = crossing_integral_bounds(liq_bounds, ∂ₜV, Dᵢ, v_i_at_Dᵢ)
        (∂ₜN_col, ∂ₜM_col, ∂ₜB_col) = integrate(integrand, bnds, quad)
        return ∂ₜN_col, ∂ₜM_col, ∂ₜB_col
    end
    return liquid_integrals
end

"""
    crossing_integral_bounds(liq_bounds, ∂ₜV, Dᵢ, [v_i_at_Dᵢ])

Insert the fall-speed crossing `∂ₜV.v_l(D) = ∂ₜV.v_i(Dᵢ)` into `liq_bounds`, so
that the derivative discontinuity of `|v_i(Dᵢ) - v_l(D)|` lies on a subinterval
boundary. `v_i_at_Dᵢ` defaults to `∂ₜV.v_i(Dᵢ)`; pass it explicitly to reuse an
already-computed value.

For a collision-rate integrand that does not carry the terminal-velocity
closures, the bounds are returned unchanged.

Called from [`get_liquid_integrals`](@ref).
"""
@inline function crossing_integral_bounds(
    liq_bounds::NTuple{2, Any}, ∂ₜV::VolumetricCollisionRate, Dᵢ, v_i_at_Dᵢ = ∂ₜV.v_i(Dᵢ),
)
    (D_min, D_max) = liq_bounds
    Dstar = crossover_diameter(v_i_at_Dᵢ, ∂ₜV.v_l, D_min, D_max)
    return (D_min, clamp(Dstar, D_min, D_max), D_max)
end
@inline crossing_integral_bounds(liq_bounds, ∂ₜV, Dᵢ, v_i_at_Dᵢ = nothing) = liq_bounds

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
    liquid_species_segment(n_liq, m_liq, ρ′_rim, ∂ₜV, v_i_at_Dᵢ, coeffs, quad, a, b)

Hand-written Gauss quadrature sum of the `(∂ₜN, ∂ₜM, ∂ₜB)` collision-rate
integrand over one segment `[a, b]`, equivalent to `integrate` applied to the
integrand built by [`get_liquid_integrals`](@ref)'s `VolumetricCollisionRate`-specialized
method. `v_i_at_Dᵢ` and `coeffs` are the outer-node quantities hoisted by the caller.
"""
@inline function liquid_species_segment(
    n_liq, m_liq, ρ′_rim::RimeDensityRate, ∂ₜV::VolumetricCollisionRate,
    v_i_at_Dᵢ, coeffs, quad, a::FT, b::FT,
) where {FT}
    zero3 = SA.SVector(zero(FT), zero(FT), zero(FT))
    a < b || return zero3
    nnodes = quad.n
    scale_factor = (b - a) / 2
    shift = (a + b) / 2
    acc = zero3
    @inbounds for i in 1:nnodes
        y = node(quad, FT(i), nnodes)
        D = scale_factor * y + shift
        w = inv_weight_fun(quad, y) * weight(quad, FT(i), nnodes)
        v_term = abs(v_i_at_Dᵢ - ∂ₜV.v_l(D))
        E = one(v_term)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
        ∂ₜN = E * evalpoly(D, coeffs) * v_term * n_liq(D)  # number collision rate
        ∂ₜM = ∂ₜN * m_liq(D)                               # mass collision rate
        ∂ₜB = ∂ₜM / rime_density_at(ρ′_rim, v_term, D)      # rime volume collision rate
        acc += SA.SVector(∂ₜN, ∂ₜM, ∂ₜB) * w
    end
    return scale_factor * acc
end

"""
    rain_B_segment(n_liq, m_liq, ρ′_rim, ∂ₜV, v_i_at_Dᵢ, coeffs, quad, a, b)

Hand-written Gauss quadrature sum of the rime-volume-collision-rate integrand
over one segment `[a, b]`, equivalent to `integrate` applied to the B_rim
integrand in [`get_liquid_integrals_rain_closed`](@ref). `v_i_at_Dᵢ` and
`coeffs` are the outer-node quantities hoisted by the caller.
"""
@inline function rain_B_segment(
    n_liq, m_liq, ρ′_rim::RimeDensityRate, ∂ₜV::VolumetricCollisionRate,
    v_i_at_Dᵢ, coeffs, quad, a::FT, b::FT,
) where {FT}
    a < b || return zero(FT)
    nnodes = quad.n
    scale_factor = (b - a) / 2
    shift = (a + b) / 2
    acc = zero(FT)
    @inbounds for i in 1:nnodes
        y = node(quad, FT(i), nnodes)
        D = scale_factor * y + shift
        w = inv_weight_fun(quad, y) * weight(quad, FT(i), nnodes)
        v_term = abs(v_i_at_Dᵢ - ∂ₜV.v_l(D))
        E = one(v_term)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
        acc += E * evalpoly(D, coeffs) * v_term * n_liq(D) * m_liq(D) / rime_density_at(ρ′_rim, v_term, D) * w
    end
    return scale_factor * acc
end

"""
    closed_rain_inner_NM_setup(ai, bi, ci, D_min, D_max, λ)

Precompute, once per point, the [`gamma_inc_moment_channel_setup`](@ref) for
the ice-velocity channel (rate `λ`, moment orders `0:5`) and for each
rain-velocity channel `j` (rate `λ + ci[j]`, moment orders `bi[j] .+ (0:5)`).
Pass the result to [`closed_rain_inner_NM`](@ref)'s `channel_setups` argument
to avoid rebuilding it at every outer ice diameter.
"""
@inline function closed_rain_inner_NM_setup(ai, bi, ci, D_min, D_max, λ)
    ice_setup = gamma_inc_moment_channel_setup(0, λ, D_min, D_max)
    rain_setups = map((b, c) -> gamma_inc_moment_channel_setup(b, λ + c, D_min, D_max), bi, ci)
    return (ice_setup, rain_setups)
end

"""
    closed_rain_inner_NM(
        v_i_at_Dᵢ, Dstar, rᵢ, ρw, ai, bi, ci, D_min, D_max, N₀r, Dr_mean;
        [channel_setups],
    )

Closed-form `(∂ₜN_col, ∂ₜM_col)` for the rain inner integral at one outer ice
diameter, where `v_i_at_Dᵢ` is the ice particle terminal velocity there and
`Dstar` the fall-speed crossing from [`crossover_diameter`](@ref).
`channel_setups` defaults to a fresh [`closed_rain_inner_NM_setup`](@ref);
pass the per-point value from [`get_liquid_integrals_rain_closed`](@ref) to
avoid rebuilding it at every outer ice diameter.
"""
function closed_rain_inner_NM(
    v_i_at_Dᵢ, Dstar, rᵢ, ρw, ai, bi, ci, D_min, D_max, N₀r, Dr_mean;
    channel_setups = closed_rain_inner_NM_setup(ai, bi, ci, D_min, D_max, inv(Dr_mean)),
)
    FT = float(eltype(ai))
    (ice_setup, rain_setups) = channel_setups
    coeffs = SA.SVector(collision_cross_section_ice_liquid_coeffs(rᵢ))

    # Each channel's `gamma_inc_moment_channel_finish` gives the six
    # consecutive-order crossing-split moments at once (cross-section orders
    # 0,1,2 for the N moment, indices 1:3; the same orders shifted by the M
    # moment's base order 3, indices 4:6), weighted here by `coeffs` and the
    # channel's velocity-curve weight (`v_i_at_Dᵢ` for the ice channel,
    # `-ai[j]` for rain channel `j`).
    function channel_NM(setup, weight)
        moments = gamma_inc_moment_channel_finish(setup, D_min, Dstar, D_max)
        N_lo = @inbounds coeffs[1] * moments[1][1]
        N_hi = @inbounds coeffs[1] * moments[1][2]
        M_lo = @inbounds coeffs[1] * moments[4][1]
        M_hi = @inbounds coeffs[1] * moments[4][2]
        @inbounds for i in 2:lastindex(coeffs)
            N_lo += coeffs[i] * moments[i][1]
            N_hi += coeffs[i] * moments[i][2]
            M_lo += coeffs[i] * moments[i + 3][1]
            M_hi += coeffs[i] * moments[i + 3][2]
        end
        return (weight * N_lo, weight * N_hi, weight * M_lo, weight * M_hi)
    end

    (s_lo_N, s_hi_N, s_lo_M, s_hi_M) = channel_NM(ice_setup, v_i_at_Dᵢ)
    @inbounds for j in eachindex(ai)
        (lo_N, hi_N, lo_M, hi_M) = channel_NM(rain_setups[j], -ai[j])
        s_lo_N += lo_N
        s_hi_N += hi_N
        s_lo_M += lo_M
        s_hi_M += hi_M
    end
    mfac = ρw * CO.volume_sphere_D(one(FT))  # m_liq(D) = mfac Dₗ³
    return (N₀r * (s_lo_N - s_hi_N), N₀r * mfac * (s_lo_M - s_hi_M))  # number: D⁰, mass: D³
end

"""
    get_liquid_integrals_rain_closed(
        psd_r::RainParticlePDF_SB2006,
        n_r, ρₐ, L_r, N_r, state, ∂ₜV, m_liq, ρ′_rim, bounds_r; quad
    )

Return a function `liquid_integrals(Dᵢ) -> (∂ₜN_col, ∂ₜM_col, ∂ₜB_col)`
where N and M are the exact incomplete-gamma closed form and
B_rim is computed by quadrature, split at the fall-speed crossing.
The velocities and the rain velocity-curve coefficients come from the
[`VolumetricCollisionRate`](@ref) `∂ₜV`.
"""
@inline function get_liquid_integrals_rain_closed(
    psd_r::CMP.RainParticlePDF_SB2006,
    n_r, ρₐ, L_r, N_r, state, ∂ₜV::VolumetricCollisionRate, m_liq, ρ′_rim::RimeDensityRate, bounds_r; quad,
)
    FT = promote_type(eltype(state), UT.promote_typeof(ρₐ, L_r, N_r))
    ρw = psd_r.ρw
    (; N₀r, Dr_mean) = CM2.pdf_rain_parameters(psd_r, L_r / ρₐ, ρₐ, N_r)
    (; v_i, v_l) = ∂ₜV
    ai, bi, ci = SA.SVector(v_l.ai), SA.SVector(v_l.bi), SA.SVector(v_l.ci)
    D_min, D_max = bounds_r
    zero_rates = (zero(FT), zero(FT), zero(FT))
    # `D_min`, `D_max`, and the rain PSD slope `λ = inv(Dr_mean)` do not depend
    # on the outer ice diameter; build the closed-form channels' incomplete-gamma
    # setup once per point instead of once per outer node.
    channel_setups = closed_rain_inner_NM_setup(ai, bi, ci, D_min, D_max, inv(Dr_mean))
    function liquid_integrals(Dᵢ)
        if iszero(FD.value(N₀r)) || !(D_max > D_min)
            return zero_rates
        end
        v_i_at_Dᵢ = v_i(Dᵢ)
        rᵢ = sqrt(ice_area(state, Dᵢ) / π)
        coeffs = collision_cross_section_ice_liquid_coeffs(rᵢ)
        Dstar = crossover_diameter(v_i_at_Dᵢ, v_l, D_min, D_max)
        ∂ₜN_col, ∂ₜM_col = closed_rain_inner_NM(
            v_i_at_Dᵢ, Dstar, rᵢ, ρw, ai, bi, ci,
            D_min, D_max, N₀r, Dr_mean; channel_setups,
        )
        if !(isfinite(∂ₜN_col) && isfinite(∂ₜM_col))
            return zero_rates
        end
        # Reuse `v_i_at_Dᵢ` and `coeffs` (hoisted above) and share `v_l(D)` between the
        # collision rate and the rime density at each inner node.
        ∂ₜB_col = integrate(
            D -> begin
                v_l_at_D = v_l(D)
                v_term = abs(v_i_at_Dᵢ - v_l_at_D)
                E = one(v_term)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
                E * evalpoly(D, coeffs) * v_term * n_r(D) * m_liq(D) / rime_density_at(ρ′_rim, v_term, D)
            end,
            (D_min, clamp(Dstar, D_min, D_max), D_max),
            quad,
        )
        return (∂ₜN_col, ∂ₜM_col, ∂ₜB_col)
    end
    return liquid_integrals
end

"""
    get_combined_liquid_integrals(
        psd_r::RainParticlePDF_SB2006,
        n_c, n_r, ρₐ, L_r, N_r, state, ∂ₜV, m_liq, ρ′_rim, bounds_c, bounds_r; quad,
    )

Return a function `combined_integrals(Dᵢ) -> (∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col, ∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col)`
evaluating the cloud and rain inner integrals at a shared outer ice diameter `Dᵢ`.
Rain N and M are the exact incomplete-gamma closed form
([`closed_rain_inner_NM`](@ref), unchanged); cloud's three integrals and rain's
B_rim are the [`liquid_species_segment`](@ref)/[`rain_B_segment`](@ref) hand-written
per-segment quadrature sums, called in interleaved cloud/rain order and sharing
the outer-node quantities `v_i_at_Dᵢ` and `coeffs` between the two species.
"""
@inline function get_combined_liquid_integrals(
    psd_r::CMP.RainParticlePDF_SB2006,
    n_c, n_r, ρₐ, L_r, N_r, state, ∂ₜV::VolumetricCollisionRate, m_liq, ρ′_rim::RimeDensityRate,
    bounds_c, bounds_r; quad,
)
    FT = promote_type(eltype(state), UT.promote_typeof(ρₐ, L_r, N_r))
    ρw = psd_r.ρw
    (; N₀r, Dr_mean) = CM2.pdf_rain_parameters(psd_r, L_r / ρₐ, ρₐ, N_r)
    (; v_l) = ∂ₜV
    ai, bi, ci = SA.SVector(v_l.ai), SA.SVector(v_l.bi), SA.SVector(v_l.ci)
    D_min_r, D_max_r = bounds_r
    # `D_min_r`, `D_max_r`, and the rain PSD slope `λ = inv(Dr_mean)` do not depend on
    # the outer ice diameter; build the closed-form channels' incomplete-gamma setup
    # once per point instead of once per outer node.
    channel_setups = closed_rain_inner_NM_setup(ai, bi, ci, D_min_r, D_max_r, inv(Dr_mean))
    function combined_integrals(Dᵢ)
        v_i_at_Dᵢ = ∂ₜV.v_i(Dᵢ)
        rᵢ = sqrt(ice_area(state, Dᵢ) / π)
        coeffs = collision_cross_section_ice_liquid_coeffs(rᵢ)

        cb = crossing_integral_bounds(bounds_c, ∂ₜV, Dᵢ, v_i_at_Dᵢ)
        c1 = liquid_species_segment(n_c, m_liq, ρ′_rim, ∂ₜV, v_i_at_Dᵢ, coeffs, quad, cb[1], cb[2])

        rain_ok = !iszero(FD.value(N₀r)) && D_max_r > D_min_r
        Dstar_r = rain_ok ? crossover_diameter(v_i_at_Dᵢ, v_l, D_min_r, D_max_r) : D_min_r
        (∂ₜN_r, ∂ₜM_r) =
            rain_ok ?
            closed_rain_inner_NM(
                v_i_at_Dᵢ, Dstar_r, rᵢ, ρw, ai, bi, ci,
                D_min_r, D_max_r, N₀r, Dr_mean; channel_setups,
            ) : (zero(FT), zero(FT))
        rain_finite = rain_ok && isfinite(∂ₜN_r) && isfinite(∂ₜM_r)
        Dstar_r_clamped = clamp(Dstar_r, D_min_r, D_max_r)

        r1 =
            rain_finite ? rain_B_segment(n_r, m_liq, ρ′_rim, ∂ₜV, v_i_at_Dᵢ, coeffs, quad, D_min_r, Dstar_r_clamped) :
            zero(FT)
        c2 = liquid_species_segment(n_c, m_liq, ρ′_rim, ∂ₜV, v_i_at_Dᵢ, coeffs, quad, cb[2], cb[3])
        r2 =
            rain_finite ? rain_B_segment(n_r, m_liq, ρ′_rim, ∂ₜV, v_i_at_Dᵢ, coeffs, quad, Dstar_r_clamped, D_max_r) :
            zero(FT)

        (∂ₜN_c, ∂ₜM_c, ∂ₜB_c) = Tuple(c1 + c2)
        if !rain_finite
            ∂ₜN_r, ∂ₜM_r = zero(FT), zero(FT)
        end
        ∂ₜB_r = r1 + r2
        return (∂ₜN_c, ∂ₜM_c, ∂ₜB_c, ∂ₜN_r, ∂ₜM_r, ∂ₜB_r)
    end
    return combined_integrals
end

@inline _rain_inner_integrals(
    psd_r::CMP.RainParticlePDF_SB2006,
    n_r, ∂ₜV::VolumetricCollisionRate{<:Any, <:Any, <:CO.Chen2022VelocityCurve},
    m_liq, ρ′_rim, bounds_r, ρₐ, L_r, N_r, state; quad,
) = get_liquid_integrals_rain_closed(
    psd_r, n_r, ρₐ, L_r, N_r, state, ∂ₜV, m_liq, ρ′_rim, bounds_r;
    quad,
)
@inline _rain_inner_integrals(
    psd_r, n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r, ρₐ, L_r, N_r, state; quad,
) = get_liquid_integrals(n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r; quad)

"""
    ∫liquid_ice_collisions(
        n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; [quad]
    )

Computes the bulk collision rate integrands between ice and liquid particles.

# Arguments
- `n_i`: ice particle size distribution function n_i(D)
- `∂ₜM_max`: maximum freezing rate function ∂ₜM_max(Dᵢ)
- `cloud_integrals`: an instance of [`get_liquid_integrals`](@ref) for cloud particles
- `rain_integrals`: an instance of [`get_liquid_integrals`](@ref) for rain particles
- `ice_bounds`: integration bounds for ice particles, from [`velocity_integral_bounds`](@ref)

# Keyword arguments
- `quad`: quadrature rule (a `Quadrature.QuadratureRule`)

# Returns
A 9-element vector of integrated rates,
`(QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL)`;
see the method below for the definition of each entry.
"""
@inline function liquid_ice_collisions_partition(
    n, ∂ₜM_max_at_Dᵢ, ∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col, ∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col,
)
    # Partition the mass collisions between freezing and shedding
    ∂ₜM_col = ∂ₜM_c_col + ∂ₜM_r_col  # [kg / s]

    ∂ₜM_frz = min(∂ₜM_col, ∂ₜM_max_at_Dᵢ)
    f_frz = iszero(FD.value(∂ₜM_col)) ? zero(∂ₜM_frz) : ∂ₜM_frz / ∂ₜM_col

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
    ))
end

@inline function ∫liquid_ice_collisions(n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; quad)
    function liquid_ice_collisions_integrands(Dᵢ)
        # Inner integrals over liquid particle diameters
        ∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col = cloud_integrals(Dᵢ)
        ∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col = rain_integrals(Dᵢ)
        return liquid_ice_collisions_partition(
            n_i(Dᵢ), ∂ₜM_max(Dᵢ), ∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col, ∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col,
        )
    end
    return integrate(liquid_ice_collisions_integrands, ice_bounds, quad)
end

"""
    ∫liquid_ice_collisions(n_i, ∂ₜM_max, combined_integrals, ice_bounds; [quad])

Same as [`∫liquid_ice_collisions`](@ref)`(n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; quad)`,
but with the cloud and rain inner integrals evaluated together by a single
`combined_integrals(Dᵢ) -> (∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col, ∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col)`
closure, e.g. [`get_combined_liquid_integrals`](@ref).
"""
@inline function ∫liquid_ice_collisions_combined(n_i, ∂ₜM_max, combined_integrals, ice_bounds; quad)
    function liquid_ice_collisions_integrands(Dᵢ)
        ∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col, ∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col = combined_integrals(Dᵢ)
        return liquid_ice_collisions_partition(
            n_i(Dᵢ), ∂ₜM_max(Dᵢ), ∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col, ∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col,
        )
    end
    return integrate(liquid_ice_collisions_integrands, ice_bounds, quad)
end

"""
    PartitionedOuter()
    SplitCorrection(; n_scan = 4, n_bisect = 4)
    SplitCorrection(quad_corr; n_scan = 4, n_bisect = 4)

How [`∫liquid_ice_collisions`](@ref) assembles its nine channels.

`PartitionedOuter` is the form the entry has always had: the freeze/shed partition is evaluated
inside the outer integral, and the two wet-growth onset diameters are inserted as subinterval
boundaries so the quadrature does not straddle the kink.

`SplitCorrection` removes the partition from the outer integral and restores it with a fixed-order
rule on the wet set, which is what [`∫liquid_ice_collisions_split`](@ref) computes. The two forms
agree to quadrature error and the second is the one a lookup table can hold, because its full-range
part carries no partition and therefore no dependence on the liquid magnitudes or the temperature.

`n_scan` and `n_bisect` set the bracket that locates the wet set.

`SplitCorrection` DEFAULTS ITS CORRECTION RULE TO THE RULE THE ENTRY IS ALREADY INTEGRATING WITH,
and does not build one. A quadrature rule is host data: `GaussLegendre(FT, n)` solves for its nodes
through `FastGaussQuadrature`, which allocates, so building one inside the assembly builds one per
evaluated cell. Measured, that is 26784 bytes per call against zero on the partitioned path, and on
a device the allocation is refused at compile time with `unsupported call to an unknown function
(call to julia.new_gc_frame)`. Nor can the assembly build it once and hold it, because the assembly
is itself a keyword default and is therefore constructed per call as well; only a rule that is
already in hand costs nothing.

The order that the wet-set correction was measured at is the order production integrates with, so
the default coincides with what was measured, and a coarser production rule coarsens the correction
with it rather than the other way round. A caller wanting a different correction rule passes one:
`SplitCorrection(rule)`.
"""
struct PartitionedOuter end

struct SplitCorrection{Q}
    n_scan::Int
    n_bisect::Int
    "the correction rule, or `nothing` to use the rule the entry integrates with"
    quad_corr::Q
end
SplitCorrection(; n_scan::Int = 4, n_bisect::Int = 4, quad_corr = nothing) =
    SplitCorrection(n_scan, n_bisect, quad_corr)
SplitCorrection(quad_corr::QuadratureRule; kwargs...) =
    SplitCorrection(; quad_corr, kwargs...)

"""
    _plain_rule(quad)

The quadrature RULE inside `quad`, which is `quad` itself unless it is a carrier of something else.

The correction rule is taken through this rather than from `quad` directly, so that a carrier
delegates to the rule it holds instead of being handed to a second integral as though it were one.
"""
@inline _plain_rule(quad) = quad

@inline _correction_rule(a::SplitCorrection, quad) =
    a.quad_corr === nothing ? _plain_rule(quad) : a.quad_corr

@inline _assemble(::PartitionedOuter, n_i, ∂ₜM_max, comb, ice_bounds, ice_bounds_plain,
    psd_c, psd_r, ∂ₜV, state, L_c, N_c, ρₐ, bounds_r, L_r, N_r; quad) =
    ∫liquid_ice_collisions_combined(n_i, ∂ₜM_max, comb, ice_bounds; quad)

@inline function _assemble(a::SplitCorrection, n_i, ∂ₜM_max, comb, ice_bounds, ice_bounds_plain,
    psd_c, psd_r, ∂ₜV, state, L_c, N_c, ρₐ, bounds_r, L_r, N_r; quad)
    balance = hybrid_wet_balance(psd_c, psd_r, ∂ₜV, ∂ₜM_max, state, L_c, N_c, ρₐ, bounds_r, L_r, N_r)
    # The full range keeps the velocity subintervals and drops only the wet onsets: the split does
    # not need them, and leaving them in would make the surrogate integral depend on a locator it no
    # longer uses.
    return _split_positional(
        n_i, ∂ₜM_max, comb, ice_bounds_plain, balance, quad, _correction_rule(a, quad),
        a.n_scan, a.n_bisect, nothing,
    )
end

"""
    _split_positional(n_i, ∂ₜM_max, comb, ice_bounds, balance, quad, quad_corr,
                      n_scan, n_bisect, full_range)

[`∫liquid_ice_collisions_split`](@ref) reached POSITIONALLY, which is not a stylistic preference.

Called with its five keywords from inside the assembly, the keyword `NamedTuple` was left with a
non-concrete tuple type and the call survived optimization as a dynamic `Core.kwcall`, which
allocates: 28128 bytes per evaluation, and a device compile refused for `gpu_gc_pool_alloc`. The
same arguments passed positionally through this wrapper resolve statically and allocate nothing.

The condition is a property of the CALL SITE rather than of any argument: every piece measured
alone allocates nothing, including the seam with the entry's own captured closure and this function
called directly with the same values. It appears only once the tabulated path makes the enclosing
method large enough, which is why it showed up when the fourth table completed the set.
"""
@inline function _split_positional(n_i, ∂ₜM_max, comb, ice_bounds, balance, quad, quad_corr,
    n_scan, n_bisect, full_range)
    return ∫liquid_ice_collisions_split(
        n_i, ∂ₜM_max, comb, ice_bounds, balance;
        quad, quad_corr, n_scan, n_bisect, full_range)
end

"""
    ∫liquid_ice_collisions_split(n_i, ∂ₜM_max, combined_integrals, ice_bounds, balance;
                                 quad, quad_corr, n_scan, n_bisect)

The nine channels of [`∫liquid_ice_collisions`](@ref), assembled as a FULL-RANGE integral that
carries no freeze/shed partition plus a correction supported on the wet set.

The partition enters the integrand only through `f_frz = min(1, ∂ₜM_max/∂ₜM_col)`, which is constant
wherever the balance `∂ₜM_col - ∂ₜM_max` does not change sign. Replacing `∂ₜM_max` by a constant
surrogate that reproduces that value therefore leaves every channel unchanged outside the wet set,
and the correction restores the true partition inside it:

```
    result = ∫_range  (surrogate partition)
           - ∫_W      (surrogate partition)
           + ∫_W      (true partition)
```

The surrogate is `floatmax` at a subfreezing state, which makes `f_frz = 1`, and zero at a state at
or above freezing, where `∂ₜM_max` vanishes identically and `f_frz = 0` is already exact. In the
second case the two correction terms have the same integrand and cancel to the last bit, so a warm
state needs no branch to skip a correction it does not want.

# Why this form

The full-range integral is what a lookup table can hold. It is linear in the liquid number
separately for cloud and for rain, no single channel mixes the two size distributions, and the
temperature and the liquid magnitudes reach the outputs only through the position of the wet set.
The correction needs no table: it is a fixed-order rule over an interval whose ends are located in
closed form, and a state with an empty wet set integrates a zero-width interval rather than taking a
different path.

`balance` is the sign function of the wet set, normally [`hybrid_wet_balance`](@ref).
`quad` integrates the full range and `quad_corr` the correction; bracket-budget (3) measured order 6
to be sufficient for the correction against a converged reference.
"""
@inline function ∫liquid_ice_collisions_split(
    n_i, ∂ₜM_max, combined_integrals, ice_bounds, balance;
    quad, quad_corr, n_scan::Int = 4, n_bisect::Int = 4,
)
    D_lo, D_hi = first(ice_bounds), last(ice_bounds)
    FT = typeof(FD.value(D_lo))
    # One scalar per state rather than a predicate per node. `compute_max_freeze_rate` returns
    # exactly zero at and above the freezing temperature, so a single evaluation decides which
    # constant partition the full-range integral must carry.
    is_warm = iszero(FD.value(∂ₜM_max(sqrt(D_lo * D_hi))))
    M̄ = ifelse(is_warm, zero(FT), floatmax(FT))
    ∂ₜM̄_max = _ -> M̄
    acc = ∫liquid_ice_collisions_combined(n_i, ∂ₜM̄_max, combined_integrals, ice_bounds; quad)
    segments = wet_set_bracket(balance, D_lo, D_hi; n_scan, n_bisect)
    for (a, b) in segments
        acc =
            acc -
            ∫liquid_ice_collisions_combined(n_i, ∂ₜM̄_max, combined_integrals, (a, b); quad = quad_corr) +
            ∫liquid_ice_collisions_combined(n_i, ∂ₜM_max, combined_integrals, (a, b); quad = quad_corr)
    end
    return acc
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
A tuple `(QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL)`, where:
1. `QCFRZ` - Cloud mass collision rate due to freezing [kg/s]
2. `QCSHD` - Cloud mass collision rate due to shedding [kg/s]
3. `NCCOL` - Cloud number collision rate [1/s]
4. `QRFRZ` - Rain mass collision rate due to freezing [kg/s]
5. `QRSHD` - Rain mass collision rate due to shedding [kg/s]
6. `NRCOL` - Rain number collision rate [1/s]
7. `∫M_col` - Total collision rate [kg/s]
8. `BCCOL` - Cloud rime volume source [m³/m³/s]
9. `BRCOL` - Rain rime volume source [m³/m³/s]
"""
@inline function ∫liquid_ice_collisions(
    state, logλ,
    psd_c, psd_r, L_c, N_c, L_r, N_r,
    aps, tps, vel, ρₐ, T, m_liq; quad, assembly = PartitionedOuter(),
)
    FT = eltype(state)

    # Particle size distributions
    n_c = DT.size_distribution(psd_c, L_c / ρₐ, ρₐ, N_c)  # n_c(Dₗ)
    n_r = DT.size_distribution(psd_r, L_r / ρₐ, ρₐ, N_r)  # n_r(Dₗ)
    n_i = DT.size_distribution(state, logλ)               # n_i(Dᵢ)

    # Integrand components; the collision rate carries the ice and liquid
    # terminal-velocity closures used below
    # NOTE: We assume collision efficiency, shape (spherical), and terminal velocity is the
    #   same for cloud and precipitating liquid particles ⟹ same volumetric collision rate, ∂ₜV
    ∂ₜV = volumetric_collision_rate_integrand(vel, ρₐ, state)  # ∂ₜV(Dᵢ, Dₗ)
    ρ′_rim = compute_local_rime_density(vel, ρₐ, T, state)  # ρ′_rim(Dᵢ, Dₗ)
    ∂ₜM_max = compute_max_freeze_rate(aps, tps, vel, ρₐ, T, state)  # ∂ₜM_max(Dᵢ)

    p = FT(0.00001)
    ice_bounds = velocity_integral_bounds(state, logλ, ∂ₜV.v_i; p)
    bounds_c = CM2.get_size_distribution_bounds(psd_c, L_c / ρₐ, ρₐ, N_c, p)
    bounds_r = CM2.get_size_distribution_bounds(psd_r, L_r / ρₐ, ρₐ, N_r, p)

    # The freeze/shed partition and the wet-growth indicator change branch at the
    # wet-growth onset diameter, so the onset is a subinterval boundary of the
    # outer integral
    # The velocity bounds BEFORE the wet onsets are inserted. `SplitCorrection` needs them: its
    # full-range integral carries no partition, so it must not be subdivided at a kink the partition
    # does not put there, and collapsing the velocity subintervals instead would give away the
    # subdivision that the `min(v, v_term_ice_max)` cap requires.
    ice_bounds_plain = ice_bounds
    (D_wet₁, D_wet₂) = wet_growth_onset_diameter(
        psd_c, psd_r, ∂ₜV, ∂ₜM_max, state,
        L_c, N_c, L_r, N_r, ρₐ, bounds_r,
        first(ice_bounds), last(ice_bounds),
    )
    ice_bounds = Tuple(
        SA.sort(
            SA.SVector(
                ice_bounds...,
                clamp(D_wet₁, first(ice_bounds), last(ice_bounds)),
                clamp(D_wet₂, first(ice_bounds), last(ice_bounds)),
            ),
        ),
    )

    return _∫liquid_ice_collisions_inner(
        psd_r, n_c, n_r, n_i, ∂ₜV, ρ′_rim, m_liq, ∂ₜM_max,
        bounds_c, bounds_r, ice_bounds, ρₐ, L_r, N_r, state; quad,
        assembly, psd_c, L_c, N_c, ice_bounds_plain,
    )
end

# Fast path: cloud and rain-B_rim inner integrals evaluated by the interleaved
# `get_combined_liquid_integrals` restructure; rain N/M stay the exact closed form
# ([`closed_rain_inner_NM`](@ref)), for the (SB2006-exp rain PSD, Chen-2022 velocity) pair.
@inline function _∫liquid_ice_collisions_inner(
    psd_r::CMP.RainParticlePDF_SB2006, n_c, n_r, n_i,
    ∂ₜV::VolumetricCollisionRate{<:Any, <:Any, <:CO.Chen2022VelocityCurve},
    ρ′_rim::RimeDensityRate, m_liq, ∂ₜM_max, bounds_c, bounds_r, ice_bounds, ρₐ, L_r, N_r, state; quad,
    assembly = PartitionedOuter(), psd_c = nothing, L_c = nothing, N_c = nothing,
    ice_bounds_plain = ice_bounds,
)
    combined_integrals = get_combined_liquid_integrals(
        psd_r, n_c, n_r, ρₐ, L_r, N_r, state, ∂ₜV, m_liq, ρ′_rim, bounds_c, bounds_r; quad,
    )
    return _assemble(assembly, n_i, ∂ₜM_max, combined_integrals, ice_bounds, ice_bounds_plain,
        psd_c, psd_r, ∂ₜV, state, L_c, N_c, ρₐ, bounds_r, L_r, N_r; quad)
end
# Numerical fallback for any other PSD/velocity type: cloud and rain inner integrals
# evaluated by two independent `get_liquid_integrals`/`_rain_inner_integrals` closures.
@inline function _∫liquid_ice_collisions_inner(
    psd_r, n_c, n_r, n_i, ∂ₜV, ρ′_rim, m_liq, ∂ₜM_max,
    bounds_c, bounds_r, ice_bounds, ρₐ, L_r, N_r, state; quad,
    assembly = PartitionedOuter(), psd_c = nothing, L_c = nothing, N_c = nothing,
    ice_bounds_plain = ice_bounds,
)
    # The split assembly is defined for the combined-integral path alone, because its correction
    # needs the closed-form rain term that only that path builds. A caller that asks for it on the
    # generic path is asking for something this method cannot supply, and saying so is better than
    # returning the partitioned result under the other name.
    assembly isa PartitionedOuter ||
        error("∫liquid_ice_collisions_split needs the SB2006 rain path; got $(typeof(psd_r))")
    cloud_integrals = get_liquid_integrals(n_c, ∂ₜV, m_liq, ρ′_rim, bounds_c; quad)
    rain_integrals = _rain_inner_integrals(
        psd_r, n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r, ρₐ, L_r, N_r, state; quad,
    )
    return ∫liquid_ice_collisions(n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; quad)
end

"""
    hybrid_wet_balance(psd_c, psd_r, ∂ₜV, ∂ₜM_max, state, L_c, N_c, ρₐ, bounds_r, L_r, N_r)

Return a function `g(Dᵢ)` of the ice diameter whose sign locates the wet set: it is positive where
the collected mass rate exceeds what the particle can freeze, and negative where it does not.

`g` costs no quadrature. The cloud term is production's own polynomial in the size-distribution
moments, with the droplet fall speed neglected against the ice fall speed, and the rain term is the
exact incomplete-gamma closed form [`closed_rain_inner_NM`](@ref) that the collision entry itself
uses. `∂ₜM_max` is closed form already.

# Why the rain term and not the cloud term

[`wet_growth_onset_diameter`](@ref) approximates BOTH terms, and evaluates the rain fall speed once
at the mean rain diameter. Measured against a 401-point scan of the true balance on 767 AMIP states,
that locator returns an EMPTY set at 4 of the 10 states whose wet set is non-empty, and those 4 hold
99.973 percent of the wet collision rate. Restoring the rain term alone and keeping the cloud
polynomial unchanged takes the miss to none of the 10, and adds no false detection at any of the 757
dry states. The whole of the defect is the rain approximation.

Scanned at 4 points and refined by 4 bisections, the resulting bracket reproduces the true-balance
bracket channel by channel at correction order 6, and costs 5.1 evaluations of `g` on average.
"""
function hybrid_wet_balance(
    psd_c::CMP.CloudParticlePDF_SB2006, psd_r::CMP.RainParticlePDF_SB2006,
    ∂ₜV::VolumetricCollisionRate{<:Any, <:Any, <:CO.Chen2022VelocityCurve},
    ∂ₜM_max, state, L_c, N_c, ρₐ, bounds_r, L_r, N_r,
)
    FT = promote_type(eltype(state), UT.promote_typeof(L_c, N_c, L_r, N_r, ρₐ))
    πFT = FT(π)
    (; v_i, v_l) = ∂ₜV
    # The cloud side: production's own three moments, formed once per state.
    (; λc, νcD, μcD) = CM2.pdf_cloud_parameters(psd_c, L_c / ρₐ, ρₐ, N_c)
    mfac = psd_c.ρw * CO.volume_sphere_D(one(FT))
    Mc₃ = mfac * DT.generalized_gamma_Mⁿ(νcD, μcD, λc, N_c, 3)
    Mc₄ = mfac * DT.generalized_gamma_Mⁿ(νcD, μcD, λc, N_c, 4)
    Mc₅ = mfac * DT.generalized_gamma_Mⁿ(νcD, μcD, λc, N_c, 5)
    # The rain side: the same closed form the entry uses, with its setup built once per state.
    (; N₀r, Dr_mean) = CM2.pdf_rain_parameters(psd_r, L_r / ρₐ, ρₐ, N_r)
    ai, bi, ci = SA.SVector(v_l.ai), SA.SVector(v_l.bi), SA.SVector(v_l.ci)
    D_min_r, D_max_r = bounds_r
    rain_live = !iszero(FD.value(N₀r)) && D_max_r > D_min_r
    rain_setup = closed_rain_inner_NM_setup(ai, bi, ci, D_min_r, D_max_r, inv(Dr_mean))
    function balance(Dᵢ)
        v = v_i(Dᵢ)
        rᵢ = sqrt(ice_area(state, Dᵢ) / πFT)
        (k₀, k₁, k₂) = collision_cross_section_ice_liquid_coeffs(rᵢ)
        cloud_rate = v * (k₀ * Mc₃ + k₁ * Mc₄ + k₂ * Mc₅)
        rain_rate = zero(cloud_rate)
        if rain_live
            Dstar = crossover_diameter(v, v_l, D_min_r, D_max_r)
            (_, Mr) = closed_rain_inner_NM(
                v, Dstar, rᵢ, psd_r.ρw, ai, bi, ci, D_min_r, D_max_r, N₀r, Dr_mean;
                channel_setups = rain_setup,
            )
            isfinite(FD.value(Mr)) && (rain_rate = Mr)
        end
        return cloud_rate + rain_rate - ∂ₜM_max(Dᵢ)
    end
    return balance
end

"""
    wet_set_bracket(g, D_lo, D_hi; n_scan, n_bisect)

The subintervals of `[D_lo, D_hi]` on which the balance `g` is positive, as an `SVector` of
`(lo, hi)` pairs padded to two entries with degenerate ones.

The balance is scanned at `n_scan` points in log diameter and each sign change is refined by
`n_bisect` bisections, so the cost is `n_scan + 1 + n_bisect * (crossings)` evaluations of `g`. A
segment of zero width contributes exactly zero to any integral over it, so a state with an empty wet
set runs the same instruction sequence as a state with a full one and no branch separates them.

The default of 4 points and 4 bisections is what bracket-budget (3) measured to reproduce a
401-point scan on every channel of the family at correction order 6.
"""
@inline function wet_set_bracket(g, D_lo, D_hi; n_scan::Int = 4, n_bisect::Int = 4)
    llo, lhi = log(D_lo), log(D_hi)
    Δl = (lhi - llo) / n_scan
    n_found = 0
    l_prev = llo
    g_prev = g(D_lo)
    roots = (D_lo, D_lo)
    for i in 1:n_scan
        l = llo + i * Δl
        D = exp(l)
        gᵢ = g(D)
        if FD.value(gᵢ) * FD.value(g_prev) < 0
            a, b = exp(l_prev), D
            sa = sign(FD.value(g_prev))
            for _ in 1:n_bisect
                m = sqrt(a * b)
                (sign(FD.value(g(m))) == sa) ? (a = m) : (b = m)
            end
            root = sqrt(a * b)
            n_found += 1
            roots = n_found == 1 ? (root, roots[2]) : (roots[1], root)
        end
        l_prev = l
        g_prev = gᵢ
    end
    # The endpoints of the candidate segments, in order, then the segments whose midpoint has a
    # positive balance. Two crossings at most is what the scan can return, and bracket-budget (3)
    # found no state with more.
    #
    # The three are written out rather than formed by `ntuple(3) do k ... end`, which captures `g`,
    # `pts` and `D_lo` in a closure that `ntuple` applies recursively. That is a simplification and
    # NOT a fix for anything: it was written to clear a device fault and, measured, it does not.
    #
    # THE FAULT IS STILL OPEN. `SplitCorrection` reached from a thin per-lane wrapper faults with
    # `ERROR_MISALIGNED_ADDRESS` on an A100 at 205312 lanes, on a plain `GaussLegendre(6)` rule and
    # on a table carrier alike, each on its own clean context; the same assembly reached through
    # `bulk_microphysics_tendencies` completes, so it depends on the enclosing kernel. Rewriting
    # this function flat changed none of that. No CPU test sees it, and `PartitionedOuter` and
    # `BulkPartition` never reach this code.
    #
    # What remains untested is the depth below: `wet_set_bracket` calls `g` about `n_scan +
    # n_scan*n_bisect + 3` times per lane, and each call runs a Brent solve inside
    # `crossover_diameter` and an incomplete-gamma chain inside `closed_rain_inner_NM`. Comparing
    # `n_scan = n_bisect = 1` against `4` separates a local-memory depth problem from the balance
    # closure itself.
    (r₁, r₂) = roots
    rlo, rhi = min(r₁, r₂), max(r₁, r₂)
    return (
        _wet_segment(g, D_lo, D_lo, rlo),
        _wet_segment(g, D_lo, rlo, rhi),
        _wet_segment(g, D_lo, rhi, D_hi),
    )
end

"""
    _wet_segment(g, D_lo, a, b)

The subinterval `(a, b)` where it is non-degenerate and the balance `g` is positive at its
geometric midpoint, and the degenerate `(D_lo, D_lo)` where it is not.

Written as a function of its arguments rather than as a closure inside `wet_set_bracket`; see the
comment at its call site for what that costs on a device.
"""
@inline function _wet_segment(g, D_lo, a, b)
    wide = FD.value(b) > FD.value(a)
    pos = wide && FD.value(g(sqrt(a * b))) > 0
    return pos ? (a, b) : (D_lo, D_lo)
end

"""
    wet_growth_onset_diameter(
        psd_c, psd_r, ∂ₜV, ∂ₜM_max, state,
        L_c, N_c, L_r, N_r, ρₐ, bounds_r, D_lo, D_hi,
    )

Return up to two ice diameters in `[D_lo, D_hi]` where the collected liquid
mass rate balances the freeze limit `∂ₜM_max` (the boundaries of the wet-growth
window; `D_lo` stands in for absent crossings). The balance is evaluated with
the liquid fall speeds simplified: the cloud collection term neglects the
droplet fall speed relative to the ice fall speed, and the rain term evaluates
the rain fall speed once at the mean rain diameter, so both reduce to polynomial
moments of their size distributions. Crossings are located on a log-spaced scan
of the interval and refined by fixed-iteration bisection.

This applies to the (`CMP.CloudParticlePDF_SB2006`,
`CMP.RainParticlePDF_SB2006`) distributions with a
[`CO.Chen2022VelocityCurve`](@ref) liquid velocity; for any other combination
`(D_lo, D_lo)` is returned.
"""
function wet_growth_onset_diameter(
    psd_c::CMP.CloudParticlePDF_SB2006, psd_r::CMP.RainParticlePDF_SB2006,
    ∂ₜV::VolumetricCollisionRate{<:Any, <:Any, <:CO.Chen2022VelocityCurve},
    ∂ₜM_max, state,
    L_c, N_c, L_r, N_r, ρₐ, bounds_r, D_lo, D_hi,
)
    (; v_i, v_l) = ∂ₜV
    FT = promote_type(eltype(state), UT.promote_typeof(L_c, N_c, L_r, N_r, ρₐ))
    πFT = FT(π)
    # Locate the onset window with the cloud droplet fall speed neglected (cloud
    # droplets fall far slower than the ice sizes that matter here) and the rain
    # fall speed evaluated once at the mean rain diameter (rain and ice fall
    # speeds are comparable, so dropping it entirely biases the search):
    #   ∫ K(D, Dₗ) n(Dₗ) m_liq(Dₗ) dDₗ = ∑ⱼ Kⱼ(rᵢ) (ρw π/6) M⁽ʲ⁺³⁾,
    # with K quadratic in Dₗ and M⁽ᵏ⁾ the (untruncated) size-distribution
    # moments. The onset diameters only bound the outer-integral subintervals;
    # the collision rates entering the physics keep the full fall-speed difference.
    (; λc, νcD, μcD) = CM2.pdf_cloud_parameters(psd_c, L_c / ρₐ, ρₐ, N_c)
    ρw = psd_c.ρw
    mfac = ρw * CO.volume_sphere_D(one(FT))
    M₃ = mfac * DT.generalized_gamma_Mⁿ(νcD, μcD, λc, N_c, 3)
    M₄ = mfac * DT.generalized_gamma_Mⁿ(νcD, μcD, λc, N_c, 4)
    M₅ = mfac * DT.generalized_gamma_Mⁿ(νcD, μcD, λc, N_c, 5)

    (; N₀r, Dr_mean) = CM2.pdf_rain_parameters(psd_r, L_r / ρₐ, ρₐ, N_r)
    mfac_r = psd_r.ρw * CO.volume_sphere_D(one(FT))
    M₃r = mfac_r * N₀r * FT(6) * Dr_mean^4    # k=3: k! = 6
    M₄r = mfac_r * N₀r * FT(24) * Dr_mean^5   # k=4: k! = 24
    M₅r = mfac_r * N₀r * FT(120) * Dr_mean^6  # k=5: k! = 120
    v_l_r = v_l(Dr_mean)

    function excess_mass_rate(D)
        v = v_i(D)
        rᵢ = sqrt(ice_area(state, D) / πFT)
        (k₀, k₁, k₂) = collision_cross_section_ice_liquid_coeffs(rᵢ)
        cloud_rate = v * (k₀ * M₃ + k₁ * M₄ + k₂ * M₅)
        rain_rate = abs(v - v_l_r) * (k₀ * M₃r + k₁ * M₄r + k₂ * M₅r)
        return cloud_rate + rain_rate - ∂ₜM_max(D)
    end
    # The balance can cross twice (a wet-growth window: collection outgrows the
    # freeze limit at intermediate sizes and falls behind again at large sizes),
    # so locate sign changes on a log-spaced grid, then refine each crossing in
    # log diameter within its grid bracket.
    llo, lhi = log(FT(D_lo)), log(FT(D_hi))
    n_scan = 16
    Δl = (lhi - llo) / n_scan
    maxiters = FT === Float32 ? 8 : 10
    tol = FixedIterations{FT}()
    # The crossing is returned with its derivative stripped, so under `ForwardDiff`
    # it carries zero partials. Brent's inverse-quadratic step divides by residual
    # differences that fall to the rounding floor near convergence, which leaves the
    # root's value finite but its partials not, and those partials would otherwise
    # poison the collision quadrature bounds. The crossing only bounds quadrature
    # subintervals, and the integrands are continuous across it apart from the
    # freezing partition, so the Leibniz boundary terms cancel and dropping the
    # derivative costs no accuracy in the rates themselves.
    #
    # The residual itself stays in the dual lane. Solving on a stripped `l` instead
    # would hand the captured maximum-freeze-rate closure a plain diameter while its
    # state is dual, and that closure types its result from its arguments, so it
    # could not represent the dual rate it computes.
    function refine_crossing(l₁, l₂)
        sol = RS.find_zero(l -> excess_mass_rate(exp(l)),
            RS.BrentsMethod(l₁, l₂), RS.CompactSolution(),
            tol, maxiters,
        )
        return FT(FD.value(exp(sol.root)))
    end
    onset₁ = FT(D_lo)
    onset₂ = FT(D_lo)
    l_prev = llo
    g_prev = excess_mass_rate(FT(D_lo))
    for i in 1:n_scan
        l = llo + i * Δl
        g = excess_mass_rate(exp(l))
        if g * g_prev < 0
            root = refine_crossing(l_prev, l)
            if onset₁ == FT(D_lo)
                onset₁ = root
            elseif onset₂ == FT(D_lo)
                onset₂ = root
            end
        end
        l_prev = l
        g_prev = g
    end
    return onset₁, onset₂
end
wet_growth_onset_diameter(
    psd_c, psd_r, ∂ₜV, ∂ₜM_max, state,
    L_c, N_c, L_r, N_r, ρₐ, bounds_r, D_lo, D_hi,
) = (D_lo, D_lo)

"""
    bulk_liquid_ice_collision_sources(
        state, logλ,
        psd_c, psd_r, L_c, N_c, L_r, N_r,
        aps, tps, vel, ρₐ, T; B_rim, quad,
    )

Computes the bulk rates for ice and liquid particle collisions.

# Arguments
- `state`: a [`P3State`](@ref)
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
- `B_rim`: the PROGNOSTIC rime volume concentration [m³/m³]. Wet-growth densification relaxes
  the `(L_rim, B_rim)` pair toward the fully-soaked solid endpoint, so it needs the prognostic
  volume rather than `L_rim/ρ_rim` reconstructed from the state's clamped and tapered quotient.
  On a consistent state the two are identical; on an inconsistent one only the prognostic form
  relaxes an above-`ρ_i` quotient back down.
- `quad`: the quadrature rule for the collision integrals

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

followed by the same transfer split by which liquid donor supplies it. What appears in that split
is decided by what a consumer needs, not by what this function happens to compute: a quantity that
is formed here and discarded cannot be reached by a caller that needs it, and one that a caller
cannot reconstruct from the other outputs must therefore be returned. `f_shd` and `∂ₜq_r_shd` are
both of that kind.

8. `∂ₜq_c_frz`: cloud mass frozen onto ice [kg/kg/s], `≤ 0`
9. `∂ₜq_c_shd`: cloud mass collected and shed as rain [kg/kg/s], `≤ 0`
10. `∂ₜq_r_frz`: rain mass frozen onto ice [kg/kg/s], `≤ 0`
11. `∂ₜq_r_shd`: rain mass collected and shed back to rain [kg/kg/s], `≤ 0`. It cancels out of
    `∂ₜq_r`, which is why it needs its own field: a consumer that wants the shed water cannot
    recover it from the net tendency.
12. `∂ₜB_rim_c`: rime volume from collected cloud [m³/m³/s]
13. `∂ₜB_rim_r`: rime volume from collected rain [m³/m³/s]
14. `f_shd`: the bulk shed fraction, dimensionless in `[0, 1]`, and zero at and above
    `T_freeze`. Exposed so a linearization can differentiate the `QIWET` and `BIWET`
    self-terms (`-f_shd/τ_wet` on both rime rows) without recomputing the collision
    quadrature; see
    [`_jacobian_2mp3_manual`](@ref
    CloudMicrophysics.BulkMicrophysicsTendencies._jacobian_2mp3_manual)'s `brim_brim`.

The bulk rates 1, 2 and 6 are linear combinations of rates 8-10, so the split
cannot be recovered from them: `∂ₜq_c + ∂ₜq_r + ∂ₜL_ice / ρₐ = 0` leaves only two
independent numbers. A linearization that needs the transfer routed to the donor
that supplies it, such as [`_jacobian_2mp3_manual`](@ref
CloudMicrophysics.BulkMicrophysicsTendencies._jacobian_2mp3_manual), takes 8-12.
"""
@inline function bulk_liquid_ice_collision_sources(
    state, logλ,
    psd_c, psd_r, L_c, N_c, L_r, N_r,
    aps, tps, vel, ρₐ, T; B_rim, quad,
)
    # `B_rim` is included in the promotion: it is a keyword, easy to miss, and
    # differentiating w.r.t. it alone (holding every other argument at plain `Float64`) is
    # exactly the scenario a single-component AD probe uses. Without it, `FT` silently narrows
    # to `Float64` and the `@NamedTuple{...::FT}` return annotation forces a
    # `Float64(::ForwardDiff.Dual)` conversion on `∂ₜB_rim`, which throws. Production's own
    # 8-wide Jacobian seeds every component together and never hits this narrowing.
    FT = promote_type(eltype(state), UT.promote_typeof(L_c, N_c, L_r, N_r, ρₐ, T, B_rim))
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
    (QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫∂ₜM_col, BCCOL, BRCOL) = rates

    # The densification driver: the shed fraction, gated to the subfreezing regime.
    #
    # `f_shd` is the fraction of the collected mass that is not frozen, and it drives the rime
    # densification below.
    #
    # Wet growth is a subfreezing process. `compute_max_freeze_rate` returns exactly zero at and
    # above `T_freeze`, so every collision is unfrozen there and the fraction would saturate at
    # one, running the densification at its maximum rate on particles that are melting. The
    # threshold is read from `tps`, the same source `compute_max_freeze_rate` reads, so the gate is
    # the exact complement of the region where the maximum freezing rate vanishes rather than an
    # approximation of it.
    #
    # The gate is applied to `f_shd` rather than to `QIWET` and `BIWET` below because `f_shd` also
    # reaches the manual Jacobian through `riming_split`, so one predicate in one place keeps the
    # rate and its two diagonals from disagreeing.
    is_subfreezing = T < TDI.T_freeze(tps)
    f_shd = (iszero(FD.value(∫∂ₜM_col)) || !is_subfreezing) ? zero(∫∂ₜM_col) :
            (QCSHD + QRSHD) / ∫∂ₜM_col

    # Shedding of rain
    # QRSHD = ∫∂ₜM_col - (QCFRZ + QRFRZ)
    NRSHD = QRSHD / m_liq(D_shd)
    # NCSHD = QCSHD / m_liq(D_shd)

    # Densification of rime: relax the `(L_rim, B_rim)` pair toward the fully-soaked solid
    # endpoint `(ρq_ice, ρq_ice/ρ_i)`, whose implied density is `ρ_i`. For `f_shd·h/τ_wet < 1`
    # the one-step update is a convex combination of the current pair and that endpoint, so the
    # implied density of the increment is `ρ_i` and the bulk quotient can only move toward it.
    #
    # `B_rim` is the PROGNOSTIC rime volume, not `ρq_ice·F_rim/ρ_rim` reconstructed from the
    # state's clamped and tapered quotient. On a consistent state the two agree; on a state
    # above `ρ_i` the reconstruction returns the volume the CLAMP implies rather than the one
    # the prognostics carry, which is exactly the volume the endpoint already is - so the term
    # sees no excess to remove and its restoring action is cancelled at precisely the states
    # that need it. With the prognostic volume, wet growth is the one process that actively
    # pulls an above-`ρ_i` quotient back down whenever it fires, and the convexity argument
    # holds without reference to the constructor's bounds.
    (; ρq_ice, F_rim) = state
    QIWET = f_shd * ρq_ice * (1 - F_rim) / τ_wet   # densification of rime mass
    BIWET = f_shd * (ρq_ice / ρ_i - B_rim) / τ_wet  # densification of rime volume

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
    ## The same transfer, split by donor. `∂ₜq_c`, `∂ₜq_r` and `∂ₜL_ice` sum to zero, so
    ## two of them determine the third and the split is not recoverable downstream.
    ∂ₜq_c_frz = -QCFRZ / ρₐ
    ∂ₜq_c_shd = -QCSHD / ρₐ
    ∂ₜq_r_frz = -QRFRZ / ρₐ
    # Collected rain that does not freeze returns to rain, so this rate cancels out of `∂ₜq_r` and
    # leaves no trace in the mass budget. It is returned anyway, because the entry's outputs are the
    # interface its consumers read and a quantity that is computed and discarded cannot be reached
    # by one that needs it. `f_shd` is the present instance: a consumer cannot form it without
    # `QRSHD`, which is why the entry forms it rather than exposing the parts. A prognostic liquid
    # fraction would need this rate directly, since the shed water is what it carries.
    ∂ₜq_r_shd = -QRSHD / ρₐ
    ∂ₜB_rim_c = BCCOL
    ∂ₜB_rim_r = BRCOL

    return @NamedTuple{
        ∂ₜq_c::FT, ∂ₜq_r::FT, ∂ₜN_c::FT, ∂ₜN_r::FT, ∂ₜL_rim::FT, ∂ₜL_ice::FT, ∂ₜB_rim::FT,
        ∂ₜq_c_frz::FT, ∂ₜq_c_shd::FT, ∂ₜq_r_frz::FT, ∂ₜq_r_shd::FT,
        ∂ₜB_rim_c::FT, ∂ₜB_rim_r::FT, f_shd::FT,
    }((
        ∂ₜq_c, ∂ₜq_r, ∂ₜN_c, ∂ₜN_r, ∂ₜL_rim, ∂ₜL_ice, ∂ₜB_rim,
        ∂ₜq_c_frz, ∂ₜq_c_shd, ∂ₜq_r_frz, ∂ₜq_r_shd, ∂ₜB_rim_c, ∂ₜB_rim_r, f_shd,
    ))
end


"""
    collision_cross_section_ice_ice(state, D_1, D_2)

Ice-ice collision cross-section [m²], `π (r(D_1) + r(D_2))²`, where the ice
effective radius is `r(D) = √(ice_area(state, D) / π)`; see [`ice_area`](@ref).
Used in [`ice_self_collection`](@ref).
"""
function collision_cross_section_ice_ice(state, D_1, D_2)
    r_eff(D) = √(ice_area(state, D) / π)
    return π * (r_eff(D_1) + r_eff(D_2))^2  # collision cross section
end

@inline r_eff_ice(state, D) = √(ice_area(state, D) / π)

"""
    ice_self_collection_inner_segment(state, c0, c1, v_1, v_ice, n_i, quad, a, b)

Hand-written Gauss quadrature sum of the ice self-collection integrand over
one inner segment `[a, b]`, equivalent to `integrate` applied to the inner
integrand in [`ice_self_collection`](@ref). `c0 = π r(D_1)²` and `c1 = 2π r(D_1)`
are the outer-node-only terms of [`collision_cross_section_ice_ice`](@ref)'s
expansion `π(r(D_1) + r(D_2))² = c0 + c1 r(D_2) + π r(D_2)²`, hoisted by the caller.
"""
@inline function ice_self_collection_inner_segment(
    state, c0, c1, v_1, v_ice, n_i, quad, a::FT, b::FT,
) where {FT}
    a < b || return zero(FT)
    nnodes = quad.n
    scale_factor = (b - a) / 2
    shift = (a + b) / 2
    acc = zero(FT)
    @inbounds for i in 1:nnodes
        y = node(quad, FT(i), nnodes)
        D_2 = scale_factor * y + shift
        w = inv_weight_fun(quad, y) * weight(quad, FT(i), nnodes)
        r_2 = r_eff_ice(state, D_2)
        cs = c0 + c1 * r_2 + π * r_2^2
        acc += cs * abs(v_1 - v_ice(D_2)) * n_i(D_2) * w
    end
    return scale_factor * acc
end

"""
    ice_self_collection_outer_segment(state, ice_bounds, D_min, D_max, v_ice, n_i, quad, a, b)

Hand-written Gauss quadrature sum of the ice self-collection outer integrand
over one segment `[a, b]`. At each outer node `D_1`, the inner integral over
`D_2 ∈ [D_1, D_max]` is evaluated by [`ice_self_collection_inner_segment`](@ref)
over each of `ice_bounds`'s subintervals, clamped up to `D_1`.
"""
@inline function ice_self_collection_outer_segment(
    state, ice_bounds, D_min, D_max, v_ice, n_i, quad, a::FT, b::FT,
) where {FT}
    a < b || return zero(FT)
    nnodes = quad.n
    scale_factor = (b - a) / 2
    shift = (a + b) / 2
    acc = zero(FT)
    @inbounds for i in 1:nnodes
        y = node(quad, FT(i), nnodes)
        D_1 = scale_factor * y + shift
        w = inv_weight_fun(quad, y) * weight(quad, FT(i), nnodes)
        v_1 = v_ice(D_1)
        n_1 = n_i(D_1)
        r_1 = r_eff_ice(state, D_1)
        c0 = π * r_1^2
        c1 = π * r_1 * 2
        D_lo = clamp(D_1, D_min, D_max)
        inner_val = zero(FT)
        for (lo, hi) in subintervals(ice_bounds)
            inner_val += ice_self_collection_inner_segment(
                state, c0, c1, v_1, v_ice, n_i, quad, max(lo, D_lo), max(hi, D_lo),
            )
        end
        acc += n_1 * inner_val * w
    end
    return scale_factor * acc
end

"""
    ice_sticking_efficiency(p3, T, F_rim)

The dimensionless efficiency with which colliding ice particles stick, `eii(T) · Eii_fact(F_rim)`.

**Applied by the CALLER as an exact scalar on the rate, never folded into the integrand.** Both
factors are per-cell scalars independent of particle size, so they come out of the collision integral
exactly; multiplying afterwards keeps [`ice_self_collection`](@ref) a pure geometric kernel. Three
things follow, and the third is the one that would otherwise bite:

  - the piecewise-linear kinks at `F_rim_lo` and `F_rim_hi` are evaluated EXACTLY per cell, where
    folding them into a tabulated kernel would interpolate across them;
  - a tabulated self-collection kernel needs no temperature axis and no regeneration when these
    values change, since neither factor is inside what it tabulates;
  - **the quadrature and tabulated paths cannot disagree**, because the factor is applied above the
    point where they diverge. Folding it in would require applying it in both, and the tabulated one
    has no temperature to apply it with.

The reference implementation applies it at its own call site, outside its lookup table, for the same
reason. See [`CMP.IceStickingEfficiency`](@ref) for the values and their provenance.
"""
@inline function ice_sticking_efficiency(p3, T, F_rim)
    (; e_cold, e_warm, T_cold, F_rim_lo, F_rim_hi) = p3.sticking
    FT = typeof(F_rim)
    # linear in T between the two saturating ends; the slope is derived from the endpoints rather
    # than hard-coded, so moving `T_cold` cannot silently leave the ramp inconsistent
    t = clamp((T - T_cold) / (p3.T_freeze - T_cold), zero(FT), one(FT))
    eii = e_cold + (e_warm - e_cold) * t
    # one below `F_rim_lo`, zero at and above `F_rim_hi`, linear between
    fact = clamp((F_rim_hi - F_rim) / (F_rim_hi - F_rim_lo), zero(FT), one(FT))
    return eii * fact
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
- `quad`: quadrature rule (a `Quadrature.QuadratureRule`)

# Returns
A `NamedTuple` of `(; dNdt)`, where:
1. `dNdt`: ice number concentration tendency due to self-collection `[1/m³/s]` (always positive or zero, represents a loss rate)
"""
@inline ice_self_collection(state, logλ, vel, ρₐ; quad) =
    _ice_self_collection(state, logλ, vel, ρₐ, quad)

# The body dispatches on `quad` POSITIONALLY, because Julia keywords do not participate in dispatch.
# A tabulated variant of this method lives in `P3_lut.jl`.
@inline function _ice_self_collection(state, logλ, vel, ρₐ, quad)
    n_i = DT.size_distribution(state, logλ)
    v_ice = ice_particle_terminal_velocity(vel, ρₐ, state)

    p = eps(one(ρₐ))
    ice_bounds = velocity_integral_bounds(state, logλ, v_ice; p)
    D_min, D_max = ice_bounds[1], ice_bounds[end]

    # Integrate the upper triangle D_1 ≤ D_2, counting each unordered particle
    # pair once - the self-collection rate.
    dNdt = zero(D_min)
    for (lo, hi) in subintervals(ice_bounds)
        dNdt += ice_self_collection_outer_segment(state, ice_bounds, D_min, D_max, v_ice, n_i, quad, lo, hi)
    end
    return (; dNdt)
end
