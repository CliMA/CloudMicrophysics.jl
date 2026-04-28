export AbstractActivationScheme,
    NoActivation,
    DiagnosticNc,
    TwomeyActivation,
    FixedARGActivation

"""
    AbstractActivationScheme

Abstract super-type for aerosol-activation closures used inside
`BulkMicrophysicsTendencies`.

Every concrete subtype is a pure, immutable parameter struct that fully
specifies one tier of activation physics. BMT dispatches on the scheme
type to compute the cloud-droplet number source `dn_lcl_activation_dt`
[kg⁻¹ s⁻¹] given the local thermodynamic state.

The three intended tiers — in order of increasing fidelity — are:

1. `DiagnosticNc` — relax `n_lcl` toward a prescribed target on a fixed
   timescale whenever cloud mass is present. This is the RCEMIP-I
   default and the pre-ARG GCM standard.
2. `TwomeyActivation` — Twomey's (1959) empirical `N_CCN = C · S^k`
   closure, gated by a minimum updraft speed.
3. `FixedARGActivation` — the full Abdul-Razzak-Ghan (2000) kernel over
   a prescribed aerosol distribution. Requires updraft, pressure, and
   a positive `S > 0`.

For the null case, use `NoActivation`, which unconditionally returns a
zero tendency.

The activation tier is chosen at model configuration time; BMT does not
attempt to auto-promote between tiers.
"""
abstract type AbstractActivationScheme <: ParametersType end

"""
    NoActivation()

Singleton scheme that suppresses aerosol activation. `activation_source`
returns zero regardless of state. Useful for tests that isolate other
microphysical processes, for warm-rain-only setups without prognostic
CCN, and for profiling.
"""
struct NoActivation <: AbstractActivationScheme end

"""
    DiagnosticNc(; N_c, q_thresh = 1e-7, τ_relax = 60)

Diagnostic-`N_c` activation scheme (Tier 1). Relaxes the droplet number
toward a prescribed target.

- While `q_lcl > q_thresh`, the droplet target is `N_c` [kg⁻¹ air] and
  the tendency is `(N_c − n_lcl) / τ_relax`.
- Otherwise the target is zero (cloud has evaporated), and the tendency
  is `-n_lcl / τ_relax`.

# Design notes
- `τ_relax` is a FIXED physical timescale, NOT the integrator timestep.
  Relaxation by `1/dt` produces an infinitely-stiff source as adaptive
  ODE solvers refine the step (documented KiD-debug lesson, issue 006
  Round 7). A fixed `τ_relax = 60 s` keeps `n_lcl` tracking cloud mass
  within one macro-step but remains integrable.
- The threshold `q_thresh` avoids enabling activation into essentially
  vapor-only cells where a tiny amount of numerical noise in `q_lcl`
  would otherwise create phantom droplets.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct DiagnosticNc{FT} <: AbstractActivationScheme
    "Droplet number target [kg⁻¹ air] while cloud mass is present"
    N_c::FT
    "Cloud-existence threshold on `q_lcl` [kg/kg]. Default: 1e-7."
    q_thresh::FT = FT(1e-7)
    "Fixed relaxation timescale [s]. Default: 60."
    τ_relax::FT = FT(60)
end

"""
    TwomeyActivation(; C, k = 0.4, w_min = 0.01, q_thresh = 1e-7,
                     τ_relax = 60)

Twomey (1959) empirical activation scheme (Tier 2). Targets a droplet
number per mass of air set by the local liquid supersaturation `S`:

    N_CCN = C · S^k    [kg⁻¹ air]

and relaxes `n_lcl` toward this target on `τ_relax` whenever both an
updraft `w > w_min` and cloud mass `q_lcl > q_thresh` are present.

When either gate is closed (subsidence, stratiform saturation, cloud
evaporation) the target is zero and existing droplets are relaxed out.

`C` is the aerosol-content parameter (typical continental ~1e8 kg⁻¹,
maritime ~5e7 kg⁻¹); `k` is the Twomey exponent (typical 0.3–0.6).

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct TwomeyActivation{FT} <: AbstractActivationScheme
    "Twomey coefficient C [kg⁻¹ air] in `N_CCN = C · S^k`"
    C::FT
    "Twomey exponent k, dimensionless. Default: 0.4."
    k::FT = FT(0.4)
    "Minimum updraft gate [m/s] below which activation is inhibited. Default: 0.01."
    w_min::FT = FT(0.01)
    "Cloud-existence threshold on `q_lcl` [kg/kg]. Default: 1e-7."
    q_thresh::FT = FT(1e-7)
    "Fixed relaxation timescale [s]. Default: 60."
    τ_relax::FT = FT(60)
end

"""
    FixedARGActivation(; act_params, distribution, q_thresh = 1e-7,
                       τ_relax = 60)

Abdul-Razzak-Ghan 2000 activation kernel (Tier 3) over a fixed aerosol
`distribution`. Computes `N_act` by calling
`CM.AerosolActivation.total_N_activated`, then relaxes `n_lcl` toward
`N_act` on `τ_relax` when the computed number exceeds the current droplet
concentration (i.e. only droplets are added, never removed by this
closure — evaporation handles the reverse).

Requires the per-cell ambient inputs `w, p` to be passed positionally to
BMT (after `inpc_log_shift`); the host can pass `zero(ρ), zero(ρ)` for
schemes that don't need them, but for ARG those zeros silently zero out
activation. Falls back to zero tendency when any of the following holds:
- supersaturation `S ≤ 0`
- updraft `w ≤ 0`
- `N_act` is not finite
- `N_act ≤ n_lcl` (kernel implies no new activation)
- cloud-mass gate `q_lcl ≤ q_thresh`

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct FixedARGActivation{FT, AP, AIP, AD} <: AbstractActivationScheme
    "ARG activation parameters (`AerosolActivationParameters`)"
    act_params::AP
    "Air properties parameters (`AirProperties`)"
    air_properties::AIP
    "Aerosol distribution (e.g. `AerosolModel.AerosolDistribution`)"
    distribution::AD
    "Cloud-existence threshold on `q_lcl` [kg/kg]. Default: 1e-7."
    q_thresh::FT = FT(1e-7)
    "Fixed relaxation timescale [s]. Default: 60."
    τ_relax::FT = FT(60)
end

# Allow broadcasting as a scalar
Base.broadcastable(x::AbstractActivationScheme) = tuple(x)
