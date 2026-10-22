"""
    het_ice_nucleation(aerosol, tps, q_lcl, N_lcl, RH, T, ρₐ)

Calculate the ice nucleation rate from heterogeneous freezing due to some `aerosol`.

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
    ice_melt(velocity_params, aps, tps, state, ρₐ, Tₐ)

Compute the melting rates of ice content and number concentration.

The rates are QIMLT in [MorrisonMilbrandt2015](@cite).

# Arguments
 - `velocity_params`: [`CMP.Chen2022VelType`](@ref)
 - `aps`: [`CMP.AirProperties`](@ref)
 - `tps`: thermodynamics parameters
 - `state`: a [`P3State`](@ref) object
 - `ρₐ`: air density
 - `Tₐ`: temperature (K)

# Returns
A `NamedTuple` `(; dNdt, dLdt)` with the melting rates of ice number
  concentration [1/m³/s] and ice content [kg/m³/s].

A particle melts by heat conduction at `dm/dt = 2π D K_therm (Tₐ - T_freeze) F_v(D) / L_f`,
for capacitance ``C = D/2`` and without the vapor diffusion contribution in subsaturated air.
Mass and number melt at the same fractional rate, that of the size distribution,

```math
f = \\frac{∫ \\frac{dm}{dt} N'(D) \\, dD}{∫ m(D) N'(D) \\, dD}, \\qquad
\\frac{dL}{dt} = L f, \\qquad \\frac{dN}{dt} = N f,
```
with `N` the prognostic ice number.
"""
@inline function ice_melt(velocity_params, aps, tps, state, ρₐ, Tₐ)
    FT = eltype(state)
    (; particles, integration, dist) = state
    # Note: process not dependent on `F_liq`
    # (we want ice core shape params)
    # Get constants
    (; K_therm) = aps
    L_f = TDI.Lf(tps, Tₐ)

    (; ρq_ice) = state
    (; T_freeze, vent) = particles.params

    v_term = ice_particle_terminal_velocity(velocity_params, particles, ρₐ)
    F_v = CO.ventilation_factor(vent, aps, v_term)
    N′ = size_distribution(dist)

    # Integrate
    bnds = velocity_integral_bounds(particles, dist, v_term; p = 1e-6)
    ∫DFvN = integrate(D -> D * F_v(D) * N′(D), bnds, integration)

    # only consider melting (not fusion)
    ΔT = max(Tₐ - T_freeze, 0)
    distribution_melt_rate = 2 * FT(π) * K_therm * ΔT / L_f * ∫DFvN  # ∫ dm/dt N′ dD
    distribution_mass = dist.ρn * exp(logLdivN(particles, dist.logλ))  # ∫ m N′ dD
    melt_frac = UT.guarded_quotient(distribution_melt_rate, distribution_mass)
    return (; dNdt = state.ρn_ice * melt_frac, dLdt = ρq_ice * melt_frac)
end
