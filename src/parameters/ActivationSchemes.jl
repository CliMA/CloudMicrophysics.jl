export AbstractActivationScheme, NoActivation, DiagnosticNc

"""
    AbstractActivationScheme

Abstract super-type for aerosol-activation closures used inside
`BulkMicrophysicsTendencies`.

Every concrete subtype is a pure, immutable parameter struct that fully
specifies one activation closure. BMT dispatches on the scheme type to
compute the cloud-droplet number source `dn_lcl_activation_dt` [kg⁻¹ s⁻¹]
given the local thermodynamic state — the same options-pattern the 1-moment
scheme uses for its prescribed-droplet-number autoconversion
([`PrescribedNd`](@ref)).

Currently provided: [`NoActivation`](@ref) (the null source, default) and
[`DiagnosticNc`](@ref) (relaxation toward a prescribed droplet number).
Supersaturation-driven closures (Twomey, Abdul-Razzak–Ghan) are deferred to
a broader activation API.
"""
abstract type AbstractActivationScheme <: ParametersType end

"""
    NoActivation()

Singleton scheme that suppresses aerosol activation: the activation source
is zero regardless of state. The default — used when the host model supplies
its own activation, and for tests that isolate other processes.
"""
struct NoActivation <: AbstractActivationScheme end

"""
    DiagnosticNc(; N_c, q_thresh = 1e-7, τ_relax = 60)

Diagnostic-`N_c` activation: relax the droplet number toward a prescribed
target. The RCEMIP-I default and the pre-ARG GCM standard.

- While `q_lcl > q_thresh`, the droplet target is `N_c` [kg⁻¹ air] and the
  tendency is `(N_c − n_lcl) / τ_relax`.
- Otherwise the target is zero (cloud has evaporated) and the tendency is
  `-n_lcl / τ_relax`.

# Design notes
- `τ_relax` is a FIXED physical timescale, NOT the integrator timestep:
  relaxation by `1/dt` produces an infinitely-stiff source as adaptive
  solvers refine the step. The 60 s default keeps `n_lcl` tracking cloud
  mass within one macro-step while remaining integrable.
- `q_thresh` avoids activating into essentially vapor-only cells, where
  numerical noise in `q_lcl` would otherwise create phantom droplets.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct DiagnosticNc{FT} <: AbstractActivationScheme
    "Droplet number target [kg⁻¹ air] while cloud mass is present"
    N_c::FT
    "Cloud-existence threshold on `q_lcl` [kg/kg]. Default: 1e-7."
    q_thresh::FT = oftype(N_c, 1e-7)
    "Fixed relaxation timescale [s]. Default: 60."
    τ_relax::FT = oftype(N_c, 60)
end

Base.broadcastable(x::AbstractActivationScheme) = tuple(x)
