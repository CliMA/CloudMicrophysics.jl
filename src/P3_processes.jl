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
    ice_ventilation_integral(velocity_params, aps, state, ρₐ)

Return the ventilation integral `∫ D F_v(D) N′(D) dD` of the ice size distribution [m⁻²].

`F_v` is the ventilation factor of [`CO.ventilation_factor`](@ref) with the coefficients of
`state.particles.params.vent`. The method of

    ice_ventilation_integral(velocity_params, aps, integration, state, ρₐ)

for the type of `state.integration` evaluates the integral. For [`P3Tables`](@ref), the integral
is `a M_1 + b M_v`, with the two moments of [`ice_ventilation_moments`](@ref) and the coefficients
`(a, b)` of [`CO.ventilation_coefficients`](@ref).

# Arguments
- `velocity_params`: [`CMP.Chen2022VelType`](@ref)
- `aps`: [`CMP.AirProperties`](@ref)
- `state`: [`P3State`](@ref)
- `ρₐ`: air density [kg/m³]
"""
@inline ice_ventilation_integral(velocity_params, aps, state, ρₐ) =
    ice_ventilation_integral(velocity_params, aps, state.integration, state, ρₐ)

@inline function ice_ventilation_integral(velocity_params, aps, integration::QuadratureRule, state, ρₐ)
    (; particles, dist) = state
    v_term = ice_particle_terminal_velocity(velocity_params, particles, ρₐ)
    F_v = CO.ventilation_factor(particles.params.vent, aps, v_term)
    N′ = size_distribution(dist)
    bnds = velocity_integral_bounds(particles, dist, v_term; p = 1e-6)
    return integrate(D -> D * F_v(D) * N′(D), bnds, integration)
end

"""
    ice_ventilation_moments(velocity_params, state, ρₐ)

Return the two moments of the ice size distribution that make up the ventilation integral.

The ventilation factor `F_v(D) = aᵥ + bᵥ Sc^(1/3) (D v(D) / ν)^(1/2)` of
[`CO.ventilation_factor`](@ref) separates the integral of [`ice_ventilation_integral`](@ref) into

```math
∫ D F_v(D) N'(D) \\, dD = a_v M_1 + b_v \\, Sc^{1/3} ν^{-1/2} M_v, \\qquad
M_1 = ∫ D N'(D) \\, dD, \\qquad M_v = ∫ D^{3/2} v(D)^{1/2} N'(D) \\, dD,
```

so the moments do not depend on the air properties, and `M_1` does not depend on the air
density. The method of

    ice_ventilation_moments(velocity_params, integration, state, ρₐ)

for the type of `state.integration` evaluates the moments. For [`P3Tables`](@ref), `M_1` is the
closed form `ρn (μ + 1) / λ` of the distribution, and `M_v` per unit number is read from the
tables, which hold it for the terminal velocity the tables were built with, so
`velocity_params` is not used.

# Arguments
- `velocity_params`: [`CMP.Chen2022VelType`](@ref)
- `state`: [`P3State`](@ref)
- `ρₐ`: air density [kg/m³]

# Returns
- `(; M_1, M_v)`: the moments `M_1` [m⁻²] and `M_v` [m⁻¹ s^(-1/2)]
"""
@inline ice_ventilation_moments(velocity_params, state, ρₐ) =
    ice_ventilation_moments(velocity_params, state.integration, state, ρₐ)

@inline function ice_ventilation_moments(velocity_params, integration::QuadratureRule, state, ρₐ)
    (; particles, dist) = state
    v_term = ice_particle_terminal_velocity(velocity_params, particles, ρₐ)
    N′ = size_distribution(dist)
    bnds = velocity_integral_bounds(particles, dist, v_term; p = 1e-6)
    (M_1, M_v) = integrate(VentilationMomentsIntegrand(N′, v_term), bnds, integration)
    return (; M_1, M_v)
end

# Integrand of `ice_ventilation_moments`: `(D N′, D^(3/2) v^(1/2) N′)` at diameter `D`
struct VentilationMomentsIntegrand{N, V} <: Function
    N′::N
    v_term::V
end
@inline function (f::VentilationMomentsIntegrand)(D)
    DN′ = D * f.N′(D)
    return SA.SVector(DN′, DN′ * sqrt(D * f.v_term(D)))
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
for capacitance ``C = D/2`` and without the vapor diffusion contribution in subsaturated air,
with `L_f` at `T_freeze`, the temperature of the particle surface.
Mass and number melt at the same fractional rate, that of the size distribution,

```math
f = \\frac{∫ \\frac{dm}{dt} N'(D) \\, dD}{∫ m(D) N'(D) \\, dD}, \\qquad
\\frac{dL}{dt} = L f, \\qquad \\frac{dN}{dt} = N f,
```
with `N` the prognostic ice number.
"""
@inline function ice_melt(velocity_params, aps, tps, state, ρₐ, Tₐ)
    FT = eltype(state)
    (; particles, dist) = state
    # Note: process not dependent on `F_liq`
    # (we want ice core shape params)
    # Get constants
    (; K_therm) = aps
    (; ρq_ice) = state
    (; T_freeze) = particles.params
    L_f = TDI.Lf(tps, T_freeze)  # at the particle surface

    ∫DFvN = ice_ventilation_integral(velocity_params, aps, state, ρₐ)

    # only consider melting (not fusion)
    ΔT = max(Tₐ - T_freeze, 0)
    distribution_melt_rate = 2 * FT(π) * K_therm * ΔT / L_f * ∫DFvN  # ∫ dm/dt N′ dD
    distribution_mass = dist.ρn * exp(logLdivN(particles, dist.logλ))  # ∫ m N′ dD
    melt_frac = UT.guarded_quotient(distribution_melt_rate, distribution_mass)
    return (; dNdt = state.ρn_ice * melt_frac, dLdt = ρq_ice * melt_frac)
end
