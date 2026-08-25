"""
    BulkMicrophysicsTendencies

Fused bulk microphysics source terms for atmospheric models.

Provides a dispatch-based API to compute all microphysics tendencies
in a single function call, enabling:
- Simplified integration in atmospheric models
- Point-evaluation suitable for quadrature over subgrid-scale fluctuations
- GPU-friendly pure function design

# Usage

```julia
using CloudMicrophysics.BulkMicrophysicsTendencies

# For 1-moment microphysics
tendencies = bulk_microphysics_tendencies(
    Instantaneous(), Microphysics1Moment(), mp, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno
)
(; dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt) = tendencies
```
"""
module BulkMicrophysicsTendencies

import ..Parameters as CMP
import ..Utilities as UT
import ..Microphysics0M as CM0
import ..Microphysics1M as CM1
import ..Microphysics2M as CM2
import ..MicrophysicsNonEq as CMNonEq
import ..P3Scheme as CMP3
import ..HetIceNucleation as CM_HetIce
import ...ThermodynamicsInterface as TDI
import ..Common as CO
import ForwardDiff as FD
import StaticArrays as SA

export MicrophysicsScheme,
    Microphysics0Moment,
    Microphysics1Moment,
    Microphysics2Moment,
    TendencyMode,
    Instantaneous,
    InstantaneousVerbose,
    LinearizedAverage,
    RosenbrockAverage,
    Verbose,
    Jacobian,
    DonorJacobian,
    CoupledDonorJacobian,
    ExactJacobian,
    GrowthTreatment,
    ImplicitGrowth,
    ExplicitGrowthDiagonal,
    TendencyLimiter,
    NoLimiter,
    EndStateSaturationAdjustment,
    rosenbrock_coupled,
    rosenbrock_exact,
    bulk_microphysics_tendencies

#####
##### Singleton types for dispatch
#####

"""
    MicrophysicsScheme

Abstract type for microphysics scheme dispatch.
"""
abstract type MicrophysicsScheme end

"""
    Microphysics0Moment <: MicrophysicsScheme

Singleton type for 0-moment microphysics scheme dispatch.
"""
struct Microphysics0Moment <: MicrophysicsScheme end

"""
    Microphysics1Moment <: MicrophysicsScheme

Singleton type for 1-moment microphysics scheme dispatch.
"""
struct Microphysics1Moment <: MicrophysicsScheme end

"""
    Microphysics2Moment <: MicrophysicsScheme

Singleton type for 2-moment microphysics scheme dispatch.

This unified scheme handles:
- Warm rain only (Seifert-Beheng 2006) when ice parameters are not provided
- Warm rain + P3 ice when ice state is provided
"""
struct Microphysics2Moment <: MicrophysicsScheme end

# --- Tendency output mode dispatch ---

"""
    TendencyMode

Abstract type for selecting the output mode of `bulk_microphysics_tendencies`.
"""
abstract type TendencyMode end

# A mode can carry configuration ([`RosenbrockAverage`](@ref) does), so it is not a
# singleton in general and a host broadcasting a tendency call over its fields has to be
# told to treat it as a scalar.
Base.broadcastable(m::TendencyMode) = tuple(m)

"""
    Instantaneous <: TendencyMode

Return raw nonlinear tendencies from a single evaluation of all microphysical
processes (no linearization, no time-averaging).
"""
struct Instantaneous <: TendencyMode end

"""
    InstantaneousVerbose <: TendencyMode

Return all individual source terms alongside aggregated `dq_*_dt` tendencies.
Useful for model diagnostics. Only works with instantaneous (nonlinear) evaluation.
"""
struct InstantaneousVerbose <: TendencyMode end

"""
    LinearizedAverage <: TendencyMode

Return time-averaged tendencies computed via repeated linearized implicit substeps.
This is the mode used operationally by ClimaAtmos.
"""
struct LinearizedAverage <: TendencyMode end

# --- The Rosenbrock-Euler substepping framework ---

"""
    Jacobian

Abstract type selecting the matrix used in each linearized-implicit substep of
[`RosenbrockAverage`](@ref). The supported types and how to add another are described in
the microphysics numerics documentation. The 1-moment scheme supports the donor-based
matrices and the exact derivative; the 2M+P3 scheme supports the exact derivative and the
hand-built matrices.
"""
abstract type Jacobian end

"""
    DonorJacobian <: Jacobian

The donor-based linearization of the tendency: each transfer is linearized in its donor
species and rate-floored. The matrix used by [`LinearizedAverage`](@ref).
"""
struct DonorJacobian <: Jacobian end

"""
    CoupledDonorJacobian <: Jacobian

The donor-based linearization with the vapor-competition and collector couplings of the
exact derivative restored.
"""
struct CoupledDonorJacobian <: Jacobian end

"""
    ExactJacobian <: Jacobian

The exact derivative of the tendency, formed with `ForwardDiff`.
"""
struct ExactJacobian <: Jacobian end

"""
    ManualJacobian <: Jacobian

A hand-built 2M+P3 substep Jacobian: the stiff condensation/deposition and
number-adjustment couplings carried as closed-form analytic derivatives, the warm-rain and
freezing transfers as donor-based linearizations, and the mixed-phase quadrature transfers
(ice melt, liquid-ice collision) donor-linearized on their own donor species with the
quadrature rates held frozen. Avoids the `ForwardDiff` pass and its `gamma_inc`
shape-derivative block.
"""
struct ManualJacobian <: Jacobian end

"""
    TemperatureCoupledJacobian <: Jacobian

The [`ManualJacobian`](@ref) extended to the temperature-coupled substep state
`(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, T)`. The phase-change rates carry
their bare relaxation `−1/τ` (the psychrometric damping emerges from the coupled `(q, T)`
block instead of the folded `−1/(τ·Γ)`), the temperature column carries the saturation
shift `−∂q_sat/∂T / τ` and the melting rate's linear temperature dependence, and the
temperature row is the latent-heating combination of the species rows.
"""
struct TemperatureCoupledJacobian <: Jacobian end

"""
    GrowthTreatment

Abstract type selecting how the positive (growth) diagonal of the Jacobian enters the
implicit operator.
"""
abstract type GrowthTreatment end

"""
    ImplicitGrowth <: GrowthTreatment

Use the Jacobian unchanged.
"""
struct ImplicitGrowth <: GrowthTreatment end

"""
    ExplicitGrowthDiagonal <: GrowthTreatment

Zero the positive diagonal of the Jacobian, so a growth mode is taken explicitly and only
the decay diagonal remains in the implicit operator.
"""
struct ExplicitGrowthDiagonal <: GrowthTreatment end

"""
    TendencyLimiter

Abstract type selecting a limiter applied to the realized substep increment.
"""
abstract type TendencyLimiter end

"""
    NoLimiter <: TendencyLimiter

Apply no limiter to the increment.
"""
struct NoLimiter <: TendencyLimiter end

"""
    EndStateSaturationAdjustment <: TendencyLimiter

Scale a substep increment so the latent-heated end state stays at or above saturation over
its more-supersaturated phase (ice when ice dominates, liquid when liquid dominates), for
cells that begin at or above saturation. Derived and analyzed in the Rosenbrock substepping
documentation.
"""
struct EndStateSaturationAdjustment <: TendencyLimiter end

"""
    RosenbrockAverage(jacobian, growth, limiter) <: TendencyMode
    RosenbrockAverage(; jacobian = DonorJacobian(), growth = ImplicitGrowth(),
        limiter = NoLimiter())

Time-averaged tendencies from repeated linearized-implicit (Rosenbrock-Euler) substeps. The
[`Jacobian`](@ref), [`GrowthTreatment`](@ref) and [`TendencyLimiter`](@ref) options select
the substep matrix, the growth-diagonal treatment and the increment limiter. See
[`rosenbrock_coupled`](@ref), [`rosenbrock_exact`](@ref), [`rosenbrock_manual`](@ref) and
[`rosenbrock_manual_temperature`](@ref) for the supported configurations.

The positivity treatment of the substep increment is not an option here: it follows from
the substep state, so a state whose components are not independent, such as the 2M+P3 rime
mass/volume pair, carries its own floor.
"""
struct RosenbrockAverage{
    J <: Jacobian, G <: GrowthTreatment, L <: TendencyLimiter,
} <: TendencyMode
    jacobian::J
    growth::G
    limiter::L
end
RosenbrockAverage(;
    jacobian = DonorJacobian(),
    growth = ImplicitGrowth(),
    limiter = NoLimiter(),
) = RosenbrockAverage(jacobian, growth, limiter)

"""
    rosenbrock_coupled()

[`RosenbrockAverage`](@ref) with the coupled donor-based Jacobian.
"""
rosenbrock_coupled() = RosenbrockAverage(CoupledDonorJacobian(), ImplicitGrowth(), NoLimiter())

"""
    rosenbrock_exact()

[`RosenbrockAverage`](@ref) with the exact Jacobian, an explicit growth diagonal, and the
end-state saturation adjustment.
"""
rosenbrock_exact() =
    RosenbrockAverage(ExactJacobian(), ExplicitGrowthDiagonal(), EndStateSaturationAdjustment())

"""
    rosenbrock_manual()

[`RosenbrockAverage`](@ref) with the hand-built 2M+P3 [`ManualJacobian`](@ref), the
explicit growth diagonal, and the end-state saturation adjustment.
"""
rosenbrock_manual() = RosenbrockAverage(
    ManualJacobian(), ExplicitGrowthDiagonal(), EndStateSaturationAdjustment(),
)

"""
    rosenbrock_manual_temperature()

[`RosenbrockAverage`](@ref) on the temperature-coupled substep state with the
[`TemperatureCoupledJacobian`](@ref), the explicit growth diagonal, and no increment
limiter: temperature evolves implicitly inside each substep, so the saturation-overshoot
adjustment is not required.
"""
rosenbrock_manual_temperature() =
    RosenbrockAverage(TemperatureCoupledJacobian(), ExplicitGrowthDiagonal(), NoLimiter())

"""
    Verbose(mode) <: TendencyMode

Diagnostic wrapper returning, alongside the net tendencies, the per-process tendencies
realized by the implicit solve of `mode`. Each process is attributed through the same
substep factorization, so the per-process tendencies sum to the net of the unlimited solve;
for a `mode` with a [`TendencyLimiter`](@ref), the wrapped net excludes the limiter. This
is a diagnostic path, separate from the model time step.
"""
struct Verbose{M <: TendencyMode} <: TendencyMode
    mode::M
end

"""
    Trace(mode) <: TendencyMode

Diagnostic wrapper recording the full substep sequence of `mode`'s implicit solve: the
state and the per-process rate breakdown evaluated at that same state, both taken AFTER
each substep, one record per substep. The per-process rates answer "what is each process
doing at the state this substep just reached", not an attribution of the increment that
reached it (that question is [`Verbose`](@ref)'s). Calls the same substep physics
production runs (`_rosenbrock_substep`) rather than a separate implementation, so a traced
run cannot drift from what the model actually does. This is a diagnostic path, separate
from the model time step; the number of substeps traced is a type parameter (`Val`), so
tracing costs nothing on any entry that does not construct a `Trace`.
"""
struct Trace{M <: TendencyMode} <: TendencyMode
    mode::M
end

# --- 1-Moment Microphysics ---

# --- Internal helpers ---

"""
Compute all individual 1-moment microphysics source terms in a single pass.

This is the **single source of truth** for which microphysical processes are
called and with what arguments. Both the raw tendency aggregation and the
linearized operator construction consume this output.

Constructs two `NamedTuple`s that are passed to all process functions
(see `Microphysics1M` module docs for the full convention):
- `micro = (; q_tot, q_lcl, q_icl, q_rai, q_sno)` — specific humidities (kg/kg)
- `thermo = (; ρ, T, w)` — air density (kg/m³), temperature (K) and vertical velocity [m/s]

Naming convention: `S_process_species1_species2`
 - process: physical mechanism (phase_change, acnv, accr, melt, accr_melt, accr_freeze)
 - species1, species2: interacting pair (not from/to)
 - `_cold` / `_warm` suffix for two-sided collision arms (inactive arm = zero)

Returns a `NamedTuple` of 18 scalar source terms plus the two relaxation timescales
`τ_phase_change_vap_lcl` / `τ_phase_change_vap_icl` used by the linearized solver.  All two-sided collision
processes are pre-routed by temperature, so consumers never need `is_warm`.
"""
@inline function _microphysics_source_terms(
    ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    # Clamp negative inputs to zero (robustness against numerical errors)
    ρ = UT.clamp_to_nonneg(ρ)
    q_tot = UT.clamp_to_nonneg(q_tot)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_icl = UT.clamp_to_nonneg(q_icl)
    q_rai = UT.clamp_to_nonneg(q_rai)
    q_sno = UT.clamp_to_nonneg(q_sno)

    FT = UT.promote_typeof(ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno)
    procs = mp.processes

    # Construct state tuples (reused across all process calls)
    micro = (; q_tot, q_lcl, q_icl, q_rai, q_sno)
    thermo = (; ρ, T, w)

    # Size-distribution / fall-speed quantities (λ⁻¹, n₀, v₀ for rain/snow/ice) are
    # pow/exp-heavy and shared by several processes per species: compute once and pass
    # to each process via its `sd` argument so they are not recomputed per process.
    sd = CM1.size_distr_parameters(mp, micro, thermo)

    # --- Phase change: vapor ↔ cloud condensate (bidirectional, ±) ---
    S_phase_change_vap_lcl = CMNonEq.conv_q_vap_to_q_lcl(procs.cloud_liquid_formation, mp, tps, micro, thermo)
    S_phase_change_vap_icl = CMNonEq.conv_q_vap_to_q_icl(procs.cloud_ice_formation, mp, tps, micro, thermo)
    # Relaxation timescales (without Γ) of the two phase changes above; `Inf` when disabled.
    # Used by the linearized solver to integrate the relaxation implicitly.
    τ_phase_change_vap_lcl = CMNonEq.τ_vap_to_q_lcl(procs.cloud_liquid_formation, mp, tps, micro, thermo)
    τ_phase_change_vap_icl = CMNonEq.τ_vap_to_q_icl(procs.cloud_ice_formation, mp, tps, micro, thermo)

    # --- Autoconversion (cloud → precipitation, ≥ 0) ---
    S_acnv_lcl_rai = CM1.conv_q_lcl_to_q_rai(procs.rain_autoconversion, mp, tps, micro, thermo)
    S_acnv_icl_sno = CM1.conv_q_icl_to_q_sno(procs.snow_autoconversion, mp, tps, micro, thermo, sd)

    # --- Accretion (collisions between species) ---
    is_warm = T >= TDI.T_freeze(tps)

    # Cloud liquid + rain → rain
    S_accr_lcl_rai = CM1.accretion(procs.cloud_liquid_rain_accretion, mp, tps, micro, thermo, sd)

    # Cloud liquid + snow: product goes to sno (cold) or rai (warm), plus thermal melt.
    # The previous accretion versions return a single number, this one has to return a tuple.
    (; S_accr, S_melt) = if isnothing(procs.cloud_liquid_snow_accretion)
        (; S_accr = zero(FT), S_melt = zero(FT))
    else
        CM1.accretion(procs.cloud_liquid_snow_accretion, mp, tps, micro, thermo, sd)
    end
    S_accr_lcl_sno_cold = ifelse(is_warm, zero(FT), S_accr)    # lcl → sno (cold)
    S_accr_lcl_sno_warm = ifelse(is_warm, S_accr, zero(FT))    # lcl → rai (warm)
    S_accr_melt_lcl_sno = S_melt                                # thermal melt of sno from warm lcl (already zero when cold)

    # Cloud ice + rain → snow (ice-side sink)
    S_accr_icl_rai = CM1.accretion(procs.cloud_ice_rain_accretion, mp, tps, micro, thermo, sd)

    # Rain frozen in cloud ice + rain collision → snow (rain sink)
    S_accr_freeze_icl_rai = CM1.accretion_rain_sink(procs.cloud_ice_rain_accretion, mp, tps, micro, thermo, sd)

    # Cloud ice + snow → snow
    S_accr_icl_sno = CM1.accretion(procs.cloud_ice_snow_accretion, mp, tps, micro, thermo, sd)

    # Rain-snow collisions: split into cold/warm arms (inactive arm = zero)
    (; S_rai_sno, S_sno_rai, S_melt) = CM1.accretion_snow_rain(procs.rain_snow_accretion, mp, tps, micro, thermo, sd)
    S_accr_rai_sno_cold = ifelse(is_warm, zero(FT), S_rai_sno) # cold arm: rai freezes → sno
    S_accr_rai_sno_warm = ifelse(is_warm, S_sno_rai, zero(FT)) # warm arm: sno melts → rai
    S_accr_melt_rai_sno = ifelse(is_warm, S_melt, zero(FT))    # thermal melt of sno from warm rai

    # --- Phase change: precipitation ↔ vapor ---
    S_phase_change_vap_rai = CM1.conv_q_rai_to_q_vap(procs.rain_condensation_evaporation, mp, tps, micro, thermo, sd)
    S_phase_change_vap_sno = CM1.conv_q_sno_to_q_vap(procs.snow_deposition_sublimation, mp, tps, micro, thermo, sd)

    # --- Melting ---
    S_melt_icl_lcl = CM1.conv_q_icl_to_q_lcl(procs.cloud_ice_melt, mp, tps, micro, thermo, sd)
    S_melt_sno_rai = CM1.conv_q_sno_to_q_rai(procs.snow_melt, mp, tps, micro, thermo, sd)

    # --- Freezing ---
    S_freeze_lcl_icl = CM1.conv_q_lcl_to_q_icl(procs.cloud_liquid_freezing, mp, tps, micro, thermo)

    return (;
        S_phase_change_vap_lcl, S_phase_change_vap_icl,
        τ_phase_change_vap_lcl, τ_phase_change_vap_icl,
        S_acnv_lcl_rai, S_acnv_icl_sno,
        S_accr_lcl_rai, S_accr_lcl_sno_cold, S_accr_lcl_sno_warm, S_accr_melt_lcl_sno,
        S_accr_icl_rai, S_accr_freeze_icl_rai, S_accr_icl_sno,
        S_accr_rai_sno_cold, S_accr_rai_sno_warm, S_accr_melt_rai_sno,
        S_phase_change_vap_rai, S_phase_change_vap_sno,
        S_melt_icl_lcl, S_melt_sno_rai,
        S_freeze_lcl_icl,
    )
end

"""
Aggregate individual source terms into the four hydrometeor tendency totals.

This is the **single location** where the sign convention of source terms
to tendency accumulators is defined.  All temperature-dependent routing is
pre-applied in `_microphysics_source_terms` (cold/warm arms), so every term
here appears with a fixed sign — no `ifelse` branching.
"""
@inline function _aggregate_tendencies(src)
    dq_lcl_dt =
        src.S_phase_change_vap_lcl - src.S_acnv_lcl_rai - src.S_accr_lcl_rai -
        src.S_accr_lcl_sno_cold - src.S_accr_lcl_sno_warm + src.S_melt_icl_lcl -
        src.S_freeze_lcl_icl

    dq_icl_dt =
        src.S_phase_change_vap_icl - src.S_acnv_icl_sno - src.S_accr_icl_rai -
        src.S_accr_icl_sno - src.S_melt_icl_lcl +
        src.S_freeze_lcl_icl

    dq_rai_dt =
        src.S_acnv_lcl_rai + src.S_accr_lcl_rai +
        src.S_accr_lcl_sno_warm + src.S_accr_melt_lcl_sno -
        src.S_accr_freeze_icl_rai -
        src.S_accr_rai_sno_cold + src.S_accr_rai_sno_warm + src.S_accr_melt_rai_sno +
        src.S_phase_change_vap_rai + src.S_melt_sno_rai

    dq_sno_dt =
        src.S_acnv_icl_sno +
        src.S_accr_lcl_sno_cold - src.S_accr_melt_lcl_sno +
        src.S_accr_icl_rai + src.S_accr_freeze_icl_rai +
        src.S_accr_icl_sno +
        src.S_accr_rai_sno_cold - src.S_accr_rai_sno_warm - src.S_accr_melt_rai_sno +
        src.S_phase_change_vap_sno - src.S_melt_sno_rai

    return (; dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt)
end

"""
    _relaxation_transfer(S, τ, Δt)

Transfer over a substep `Δt` of a linear relaxation with instantaneous rate `S`
and timescale `τ` (Morrison & Milbrandt 2015, Appendix C, eqs. C6-C7):
`S τ (1 - exp(-Δt/τ))`, i.e. `S Δt φ(Δt/τ)` with `φ(x) = (1 - e⁻ˣ)/x`. It tends
to `S Δt` for `Δt ≪ τ` (the instantaneous rate is recovered) and to `S τ` for
`Δt ≫ τ` (the transfer saturates at the equilibrium amount, so the substep never
overshoots the equilibrium). Only `x φ(x) = 1 - e⁻ˣ` is needed, so there is no
division by `x`; `τ` is clamped to `[eps, floatmax]` so that a disabled process
(`τ = Inf`, `S = 0`) yields `0` rather than `0 ⋅ Inf`.
"""
@inline function _relaxation_transfer(S, τ, Δt)
    τ_c = clamp(τ, eps(typeof(τ)), floatmax(typeof(τ)))
    return S * τ_c * -expm1(-Δt / τ_c)
end

"""
Construct a local linear approximation of 1-moment microphysics tendencies
from pre-computed source terms:

    dq/dt ≈ M * q + e

using a donor-based linearization:
- donor → receiver transfers are represented as `D * q_donor`
- vapor ↔ cloud condensate phase changes are relaxations toward their
  equilibrium, `S = (q* - q)/τ`. Their transfer over the substep is the time
  average of the exact relaxation, `Δq = S Δt φ(Δt/τ)` with
  `φ(x) = (1 - exp(-x))/x` (Morrison & Milbrandt 2015, Appendix C), so the
  substep never crosses `q*` for any `Δt/τ` and the instantaneous rate is
  recovered for `Δt ≪ τ`. A source (`Δq > 0`) enters `e` as a non-negative
  constant. A sink (`Δq < 0`) enters `M` as an implicit decay `-D q` whose
  coefficient, `D = |Δq| / ((q + Δq) Δt)`, removes exactly `|Δq|` when acting
  alone; combined with other sinks it keeps `q ≥ 0`, and, unlike the plain
  `S/q` decay, it does not over-sublimate when `τ ≪ Δt`. This matters when τ is
  a few seconds (e.g. `PrescribedIceNumber` with a large prescribed `N_0`),
  where feeding the instantaneous rate to the substep produced a
  deposition/sublimation flip-flop.
- vapor → snow deposition is treated as a constant source (`e`)
- other condensate sinks are treated as linear sinks (`-D * q`)

The other `D` coefficients use `D = S / max(q_min, q_donor)` for robustness.

Returns a `NamedTuple` containing the nonzero entries of `M` and `e`.
"""
@inline function _linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min, Δt)
    FT = typeof(src.S_phase_change_vap_lcl)

    M11 = zero(FT)
    M12 = zero(FT)
    M21 = zero(FT)
    M22 = zero(FT)
    M31 = zero(FT)
    M33 = zero(FT)
    M34 = zero(FT)
    M41 = zero(FT)
    M42 = zero(FT)
    M43 = zero(FT)
    M44 = zero(FT)
    e1 = zero(FT)
    e2 = zero(FT)
    e4 = zero(FT)

    # --- Phase change: vapor ↔ cloud condensate (time-averaged relaxation) ---
    # The non-equilibrium schemes relax the condensate toward q* = q + S τ at rate
    # 1/τ (τ excludes Γ, which is already inside S; `Inf` when the process is
    # disabled). Δq is the transfer over the substep (see `_relaxation_transfer`).
    # A source enters `e` as a non-negative constant; a sink enters `M` as an
    # implicit decay whose coefficient removes exactly |Δq| on its own
    # (|Δq| ≤ |S| τ ≤ q by the tendency bound, so q + Δq ≥ 0 up to round-off; the
    # `q_min` floor keeps the coefficient finite when the whole pool sublimates).
    Δq = _relaxation_transfer(src.S_phase_change_vap_lcl, src.τ_phase_change_vap_lcl, Δt)
    e1 += max(zero(Δq), Δq) / Δt
    M11 -= ifelse(Δq < zero(Δq), -Δq / (max(q_lcl + Δq, q_min) * Δt), zero(Δq))

    Δq = _relaxation_transfer(src.S_phase_change_vap_icl, src.τ_phase_change_vap_icl, Δt)
    e2 += max(zero(Δq), Δq) / Δt
    M22 -= ifelse(Δq < zero(Δq), -Δq / (max(q_icl + Δq, q_min) * Δt), zero(Δq))

    # --- Melt: ice cloud → liquid cloud ---
    D = src.S_melt_icl_lcl / max(q_min, q_icl)
    M22 -= D
    M12 += D

    # --- Freeze: liquid cloud → ice cloud ---
    D = src.S_freeze_lcl_icl / max(q_min, q_lcl)
    M11 -= D
    M21 += D

    # --- Autoconversion: donor-based transfer ---
    D = src.S_acnv_lcl_rai / max(q_min, q_lcl)
    M11 -= D
    M31 += D

    D = src.S_acnv_icl_sno / max(q_min, q_icl)
    M22 -= D
    M42 += D

    # --- Accretion: donor-based transfer ---
    D = src.S_accr_lcl_rai / max(q_min, q_lcl)
    M11 -= D
    M31 += D

    # lcl + sno accretion (cold/warm arms already zeroed)
    D_cold = src.S_accr_lcl_sno_cold / max(q_min, q_lcl)
    D_warm = src.S_accr_lcl_sno_warm / max(q_min, q_lcl)
    M11 -= D_cold + D_warm
    M31 += D_warm           # warm: lcl → rai
    M41 += D_cold           # cold: lcl → sno

    # thermal melt of sno from warm lcl
    D = src.S_accr_melt_lcl_sno / max(q_min, q_sno)
    M44 -= D
    M34 += D

    D = src.S_accr_icl_rai / max(q_min, q_icl)
    M22 -= D
    M42 += D

    D = src.S_accr_icl_sno / max(q_min, q_icl)
    M22 -= D
    M42 += D

    # rain frozen in icl + rai collision
    D = src.S_accr_freeze_icl_rai / max(q_min, q_rai)
    M33 -= D
    M43 += D

    # warm arm: sno melts → rai (already zero when cold)
    D = src.S_accr_rai_sno_warm / max(q_min, q_sno)
    M44 -= D
    M34 += D

    # thermal melt of sno from warm rai (already zero when cold)
    D = src.S_accr_melt_rai_sno / max(q_min, q_sno)
    M44 -= D
    M34 += D

    # cold arm: rai freezes → sno (already zero when warm)
    D = src.S_accr_rai_sno_cold / max(q_min, q_rai)
    M33 -= D
    M43 += D

    # --- Rain phase change: sink to vapor (always zero or negative) ---
    D = (-src.S_phase_change_vap_rai) / max(q_min, q_rai)
    M33 -= D

    # --- Snow phase change: deposition/sublimation ---
    D = src.S_phase_change_vap_sno / max(q_min, q_sno)
    is_source = src.S_phase_change_vap_sno >= zero(FT)
    e4 += ifelse(is_source, src.S_phase_change_vap_sno, zero(FT))
    M44 += ifelse(is_source, zero(FT), D)

    # --- Snow melt: snow → rain ---
    D = src.S_melt_sno_rai / max(q_min, q_sno)
    M44 -= D
    M34 += D

    return (
        M11 = M11, M12 = M12, M21 = M21, M22 = M22,
        M31 = M31, M33 = M33, M34 = M34,
        M41 = M41, M42 = M42, M43 = M43, M44 = M44,
        e1 = e1, e2 = e2, e4 = e4,
    )
end

"""
Compute time-averaged 1-moment microphysics tendencies over a single linearized substep.

Solves the linearized implicit system

    (q* - q⁰) / Δt = M q* + e

and returns the average tendency

    dq/dt = (q* - q⁰) / Δt.

The system uses a sparse structure specific to the 1-moment microphysics model.
`q_lcl` and `q_icl` as well as `q_rai` and `q_sno` are solved from a coupled 2×2 system.

Because sinks are linearized as `-D q`, they are effectively integrated as
exponential decays over the substep. The vapor → cloud condensate sources in `e`
are the time-averaged transfers over the substep (see `_linearize`); evaporation and
sublimation are part of `M`.
"""
@inline function _linearized_implicit_step(
    ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt::AbstractFloat,
)

    FT = typeof(q_tot)

    src = _microphysics_source_terms(
        Microphysics1Moment(), mp, tps,
        ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )
    q_min = TDI.TD.Parameters.q_min(tps)
    lin = _linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min, Δt)

    invΔt = one(FT) / Δt

    # Cap vap→condensate sources jointly so the substep cannot drive
    # `q_v` below `min(q_sat_liq, q_sat_ice)`. Preserves relative rates.
    q_sat_min = min(
        TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ),
        TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ),
    )
    q_v = q_tot - q_lcl - q_icl - q_rai - q_sno
    α = min(
        one(FT),
        max(zero(FT), q_v - q_sat_min) * invΔt /
        max(lin.e1 + lin.e2 + lin.e4, eps(FT)),
    )

    # A = I/Δt - M
    a11 = invΔt - lin.M11
    a12 = -lin.M12
    a21 = -lin.M21
    a22 = invΔt - lin.M22
    a31 = -lin.M31
    a33 = invΔt - lin.M33
    a34 = -lin.M34
    a41 = -lin.M41
    a42 = -lin.M42
    a43 = -lin.M43
    a44 = invΔt - lin.M44

    # rhs = e + q_0/Δt (vap→condensate `e` terms scaled by `α` above)
    # e3 = 0 by the 1m model
    b1 = α * lin.e1 + invΔt * q_lcl
    b2 = α * lin.e2 + invΔt * q_icl
    b3 = invΔt * q_rai
    b4 = α * lin.e4 + invΔt * q_sno

    # Solve 2×2 system for q_lcl, q_icl (coupled via ice melt M12 and liquid freezing M21)
    det12 = muladd(-a12, a21, a11 * a22)
    q_lcl_new = (b1 * a22 - a12 * b2) / det12
    q_icl_new = (a11 * b2 - a21 * b1) / det12

    # Reduced 2x2 system for q_rai_new, q_sno_new
    r3 = muladd(-a31, q_lcl_new, b3)
    r4 = muladd(-a41, q_lcl_new, muladd(-a42, q_icl_new, b4))

    det = muladd(-a34, a43, a33 * a44)
    # det is a positive because a33a44 is guaranteed to be larger than a34a43
    q_rai_new = (r3 * a44 - a34 * r4) / det
    q_sno_new = (a33 * r4 - r3 * a43) / det

    dq_lcl_dt = (q_lcl_new - q_lcl) * invΔt
    dq_icl_dt = (q_icl_new - q_icl) * invΔt
    dq_rai_dt = (q_rai_new - q_rai) * invΔt
    dq_sno_dt = (q_sno_new - q_sno) * invΔt

    return (; dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt)
end

# --- Public API: bulk_microphysics_tendencies with TendencyMode dispatch ---

"""
    bulk_microphysics_tendencies(
        ::Instantaneous, ::Microphysics1Moment, mp, tps,
        ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )

Compute all 1-moment microphysics tendencies in one fused call.

Returns a NamedTuple with the aggregated tendency for each water category.
This is a pure function of local thermodynamic state, suitable for:
- Point quadrature over subgrid-scale (T, q_tot) distributions
- GPU kernel evaluation
- Unit testing in isolation

# Arguments
- `mp`: Microphysics1MParams parameter container
- `tps`: Thermodynamics parameters
- `ρ`: Air density [kg/m³]
- `T`: Temperature [K]
- `w`: air vertical velocity [m/s]
- `q_tot`: Total water specific content [kg/kg]
- `q_lcl`: Cloud liquid water specific content [kg/kg]
- `q_icl`: Cloud ice specific content [kg/kg]
- `q_rai`: Rain specific content [kg/kg]
- `q_sno`: Snow specific content [kg/kg]

# Returns
`NamedTuple` with fields:
- `dq_lcl_dt`: Cloud liquid tendency [kg/kg/s]
- `dq_icl_dt`: Cloud ice tendency [kg/kg/s]
- `dq_rai_dt`: Rain tendency [kg/kg/s]
- `dq_sno_dt`: Snow tendency [kg/kg/s]

# Notes
- Negative specific contents are clamped to zero for robustness.
- Does NOT apply timestep-dependent limiters.
"""
@inline function bulk_microphysics_tendencies(
    ::Instantaneous, ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    src = _microphysics_source_terms(
        Microphysics1Moment(), mp, tps,
        ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )
    return _aggregate_tendencies(src)
end

"""
    bulk_microphysics_tendencies(
        ::InstantaneousVerbose, ::Microphysics1Moment, mp, tps,
        ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )

Compute all 1-moment microphysics tendencies and return both aggregated
tendencies (`dq_*_dt`), all individual source terms (`S_*`) and the two
`τ_phase_change_vap_*` relaxation timescales.

Useful for model diagnostics. The `dq_*_dt` fields are identical to those
returned by `Instantaneous()`.

# Returns
`NamedTuple` with all fields from `Instantaneous()` plus individual source
terms: `S_phase_change_vap_lcl`, `S_phase_change_vap_icl`, `S_acnv_lcl_rai`,
`S_acnv_icl_sno`, etc.
"""
@inline function bulk_microphysics_tendencies(
    ::InstantaneousVerbose, ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    src = _microphysics_source_terms(
        Microphysics1Moment(), mp, tps,
        ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )
    agg = _aggregate_tendencies(src)
    return merge(agg, src)
end

"""
    bulk_microphysics_tendencies(
        ::LinearizedAverage, ::Microphysics1Moment, mp, tps,
        ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1,
    )

Compute average 1-moment microphysics tendencies over `Δt` using repeated
linearized implicit substeps.

The interval `Δt` is divided into `nsub` equal substeps. At each substep, a local
linearized microphysics system is rebuilt from the current state and solved
implicitly for cloud liquid, cloud ice, rain, and snow. Temperature is then
updated from the latent heating implied by the substep tendencies.

The returned tendencies are the net change in the hydrometeor species over the
full interval divided by `Δt`.

Increasing `nsub` improves how well the method captures nonlinear changes in the
active microphysical processes, including regime changes near freezing.

# Returns
`NamedTuple` with fields:
- `dq_lcl_dt`: Cloud liquid tendency [kg/kg/s]
- `dq_icl_dt`: Cloud ice tendency [kg/kg/s]
- `dq_rai_dt`: Rain tendency [kg/kg/s]
- `dq_sno_dt`: Snow tendency [kg/kg/s]
"""
@inline function bulk_microphysics_tendencies(
    ::LinearizedAverage,
    cm::Microphysics1Moment,
    mp::CMP.Microphysics1MParams,
    tps,
    ρ,
    T,
    w,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    Δt::AbstractFloat,
    nsub::Integer = 1,
)
    FT = typeof(q_tot)

    q_lcl_0 = q_lcl
    q_icl_0 = q_icl
    q_rai_0 = q_rai
    q_sno_0 = q_sno

    Δt_sub = Δt / FT(nsub)

    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)

    for _ in 1:nsub
        rates = _linearized_implicit_step(
            cm,
            mp,
            tps,
            ρ,
            T,
            w,
            q_tot,
            q_lcl,
            q_icl,
            q_rai,
            q_sno,
            Δt_sub,
        )

        q_lcl += rates.dq_lcl_dt * Δt_sub
        q_icl += rates.dq_icl_dt * Δt_sub
        q_rai += rates.dq_rai_dt * Δt_sub
        q_sno += rates.dq_sno_dt * Δt_sub

        T +=
            (
                Lv_over_cp * (rates.dq_lcl_dt + rates.dq_rai_dt) +
                Ls_over_cp * (rates.dq_icl_dt + rates.dq_sno_dt)
            ) * Δt_sub
    end

    dq_lcl_dt = (q_lcl - q_lcl_0) / Δt
    dq_icl_dt = (q_icl - q_icl_0) / Δt
    dq_rai_dt = (q_rai - q_rai_0) / Δt
    dq_sno_dt = (q_sno - q_sno_0) / Δt

    return (; dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt)
end


# --- 0-Moment Microphysics ---
"""
    bulk_microphysics_tendencies(::Microphysics0Moment, mp, tps, T, q_lcl, q_icl)
    bulk_microphysics_tendencies(::Microphysics0Moment, mp, tps, T, q_lcl, q_icl, q_vap_sat)

Compute 0-moment microphysics tendencies in one fused call.

Returns the total water tendency `dq_tot_dt` (a scalar, in kg/kg/s) from precipitation removal.

The first form uses the fixed condensate threshold `qc_0`;
the second form uses the supersaturation threshold `S_0 * q_vap_sat`.

# Arguments
- `mp`: Microphysics0MParams (contains τ_precip, qc_0, S_0)
- `tps`: Thermodynamics parameters
- `T`: Temperature [K]
- `q_lcl`: Cloud liquid specific content [kg/kg]
- `q_icl`: Cloud ice specific content [kg/kg]
- `q_vap_sat`: (second method only) Saturation specific humidity [kg/kg]

# Notes
- Does NOT apply limiters (caller applies based on timestep)
"""
@inline function bulk_microphysics_tendencies(
    ::Microphysics0Moment, mp::CMP.Microphysics0MParams, tps,
    T, q_lcl, q_icl,
)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_icl = UT.clamp_to_nonneg(q_icl)
    dq_tot_dt = CM0.remove_precipitation(mp.precip, q_lcl, q_icl)
    return dq_tot_dt
end
@inline function bulk_microphysics_tendencies(
    ::Microphysics0Moment,
    mp::CMP.Microphysics0MParams,
    tps,
    T,
    q_lcl,
    q_icl,
    q_vap_sat,
)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_icl = UT.clamp_to_nonneg(q_icl)
    dq_tot_dt = CM0.remove_precipitation(mp.precip, q_lcl, q_icl, q_vap_sat)
    return dq_tot_dt
end

# --- 2-Moment Microphysics Helper Functions ---

"""
    _bare_rate_conv_q_vap_to_q_lcl(τ, tps, micro, thermo)

The phase-change rate `sat_excess / τ` toward saturation over liquid, with the evaporation
branch bounded by the donor content.

`τ` here is a DERIVED diffusional growth timescale, [`CM2.cloud_condensation_timescale`](@ref),
rather than a phenomenological closure parameter, and that is what decides the form. The
cell-scale relaxation `ds/dt = -Γ s / τ` comes out of the coupled mass and energy bookkeeping on
its own, so folding the psychrometric factor `Γ` into `τ` at this call site would apply that
correction a second time. The shared [`CMNonEq.conv_q_vap_to_q_lcl`](@ref) closes over a
phenomenological relaxation timescale instead and keeps the fold, which is the right pairing
there: a prescribed constant is given its cleanest meaning as a uniform supersaturation-relaxation
time. The boundary is which `τ`, not which call site.
"""
@inline function _bare_rate_conv_q_vap_to_q_lcl(τ, tps, micro, thermo)
    (; q_lcl) = micro
    (; ρ, T) = thermo
    qᵥ = TDI.q_vap(micro.q_tot, micro.q_lcl + micro.q_rai, micro.q_icl + micro.q_sno)
    qᵥ_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
    sat_excess = qᵥ - qᵥ_sat_liq
    return ifelse(
        sat_excess < 0,
        -min(-sat_excess, max(0, q_lcl)) / τ,
        sat_excess / τ,
    )
end

"""
    warm_rain_tendencies_2m(sb, q_lcl, q_rai, ρ, n_lcl, n_rai)

Internal helper function that computes 2M warm rain processes:
cloud condensation/evaporation, autoconversion, self-collection, accretion,
rain breakup, rain evaporation, and number adjustment for mass limits.

Used by both warm-only and warm+ice dispatch methods to reduce code duplication.

# Arguments
- `warm_rain`: warm-rain parameters (contains the SB2006 scheme, air
  properties, and condensation/evaporation options)
- `tps`: Thermodynamics parameters
- `T`: Temperature (K)
- `q_tot`: Total water specific content (kg/kg)
- `q_lcl`: Cloud liquid specific content (kg/kg)
- `q_rai`: Rain specific content (kg/kg)
- `q_ice`: Ice specific content (kg/kg)
- `ρ`: Air density (kg/m³)
- `n_lcl`: Cloud droplet number per kg air (1/kg)
- `n_rai`: Rain number per kg air (1/kg)
- `w`: Vertical velocity (m/s), default `0`
- `p`: Air pressure (Pa), default `0`

# Returns
`NamedTuple` with warm rain tendencies:
- `dq_lcl_dt`: Cloud liquid mass tendency (kg/kg/s)
- `dq_rai_dt`: Rain mass tendency (kg/kg/s)
- `dn_lcl_dt`: Cloud number tendency (1/kg/s)
- `dn_rai_dt`: Rain number tendency (1/kg/s)
- `dn_lcl_activation_dt`: Cloud number activation tendency (1/kg/s)
"""
@inline function warm_rain_tendencies_2m(
    warm_rain, tps, T, q_tot, q_lcl, q_rai, q_ice, ρ, n_lcl, n_rai,
    w = zero(ρ), p = zero(ρ),
)

    # Unpack parameters
    sb = warm_rain.seifert_beheng
    aps = warm_rain.air_properties

    # Convert to number densities for CM2 functions
    N_lcl = ρ * n_lcl
    N_rai = ρ * n_rai

    # Initialize tendencies
    FT = typeof(ρ)
    dq_lcl_dt = zero(FT)
    dq_rai_dt = zero(FT)
    dn_lcl_dt = zero(FT)
    dn_rai_dt = zero(FT)

    # --- Aerosol activation ---
    dn_lcl_activation_dt = zero(FT)

    # --- Condensation of vapor / evaporation of cloud liquid water ---
    # The relaxation timescale follows the droplet population's capacitance integral, so the
    # unbounded timescale diverges as the population vanishes. `cloud_condensation_timescale`
    # caps it at `CLOUD_COND_TIMESCALE_MAX` to stay finite, and `sat_excess / τ_max` is not zero,
    # so an existence threshold is required here. The same rate and the same gate serve the
    # per-process decomposition, whose sum has to equal what this entry returns.
    micro_mock = (; q_tot, q_lcl, q_icl = q_ice, q_rai, q_sno = zero(q_ice))
    thermo_mock = (; ρ, T)
    τ_cond = CM2.cloud_condensation_timescale(sb.pdf_c, aps, tps, T, ρ, q_lcl, N_lcl)
    ∂ₜq_lcl_cond = _bare_rate_conv_q_vap_to_q_lcl(τ_cond, tps, micro_mock, thermo_mock)
    # No droplets, no surfaces: condensation AND evaporation are both exactly zero.
    ∂ₜq_lcl_cond = ifelse(
        CM2.cloud_condensation_is_degenerate(τ_cond), zero(∂ₜq_lcl_cond), ∂ₜq_lcl_cond)
    ∂ₜn_lcl_cond = zero(∂ₜq_lcl_cond)  # neglect number change from condensation/evaporation
    dq_lcl_dt += ∂ₜq_lcl_cond
    dn_lcl_dt += ∂ₜn_lcl_cond

    # --- Evaporation of rain ---
    evap = CM2.rain_evaporation(sb, aps, tps, q_tot, q_lcl, q_ice, q_rai, zero(q_ice), ρ, N_rai, T)
    dq_rai_dt += evap.∂ₜq_rai
    dn_rai_dt += evap.∂ₜρn_rai / ρ

    # --- Autoconversion ---
    acnv = CM2.autoconversion(sb.acnv, sb.pdf_c, q_lcl, q_rai, ρ, N_lcl)
    dq_lcl_dt += acnv.dq_lcl_dt
    dq_rai_dt += acnv.dq_rai_dt
    dn_lcl_dt += acnv.dN_lcl_dt / ρ
    dn_rai_dt += acnv.dN_rai_dt / ρ

    # --- Cloud liquid self-collection ---
    ∂ₜN_lcl_sc = CM2.cloud_liquid_self_collection(sb.acnv, sb.pdf_c, q_lcl, ρ, N_lcl, acnv.dN_lcl_dt)
    dn_lcl_dt += ∂ₜN_lcl_sc / ρ

    # --- Accretion ---
    accr = CM2.accretion(sb, q_lcl, q_rai, ρ, N_lcl)
    dq_lcl_dt += accr.dq_lcl_dt
    dq_rai_dt += accr.dq_rai_dt
    dn_lcl_dt += accr.dN_lcl_dt / ρ

    # --- Rain self-collection ---
    ∂ρn_rai_sc_∂t = CM2.rain_self_collection(sb.pdf_r, sb.self, q_rai, ρ, N_rai)
    dn_rai_dt += ∂ρn_rai_sc_∂t / ρ

    # --- Rain breakup ---
    ∂ρn_rai_br_∂t = CM2.rain_breakup(sb.pdf_r, sb.brek, q_rai, ρ, N_rai, ∂ρn_rai_sc_∂t)
    dn_rai_dt += ∂ρn_rai_br_∂t / ρ

    # --- Number adjustment for mass limits ---
    # Cloud liquid
    numadj_lcl = (; sb.numadj.τ, x_min = sb.pdf_c.xc_min, x_max = sb.pdf_c.xc_max)
    ∂ₜn_lcl_numadj = CM2.number_tendency_from_mass_limits(numadj_lcl, q_lcl, n_lcl)
    dn_lcl_dt += ∂ₜn_lcl_numadj
    # Rain
    numadj_rai = (; sb.numadj.τ, x_min = sb.pdf_r.xr_min, x_max = sb.pdf_r.xr_max)
    ∂ₜn_rai_numadj = CM2.number_tendency_from_mass_limits(numadj_rai, q_rai, n_rai)
    dn_rai_dt += ∂ₜn_rai_numadj

    return (; dq_lcl_dt, dq_rai_dt, dn_lcl_dt, dn_rai_dt, dn_lcl_activation_dt)
end

# --- 2-Moment Microphysics (Unified Warm + Optional Ice) ---


"""
    bulk_microphysics_tendencies(
        ::Microphysics2Moment,
        mp::Microphysics2MParams{WR, Nothing}, tps,
        ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai,
        q_ice = 0, n_ice = 0, q_rim = 0, b_rim = 0, logλ = 0,
        inpc_log_shift = 0, w = 0, p = 0,
    )

Compute 2-moment **warm rain only** microphysics tendencies (Seifert-Beheng 2006).

This method is type-stable and GPU-optimized for warm rain processes only.
For warm rain + P3 ice, see the method that accepts `Microphysics2MParams{WR, <:P3IceParams}`.

# Arguments
- `mp`: Microphysics2MParams with `mp.ice == nothing` (warm rain only)
- `tps`: Thermodynamics parameters
- `ρ`: Air density (kg/m³)
- `T`: Temperature (K)
- `q_tot`: Total water specific content (kg/kg)
- `q_lcl`: Cloud liquid specific content (kg/kg)
- `n_lcl`: Cloud droplet number per kg air (1/kg)
- `q_rai`: Rain specific content (kg/kg)
- `n_rai`: Rain number per kg air (1/kg)
- `q_ice`, `n_ice`, `q_rim`, `b_rim`, `logλ`: (optional, default `0`)
  ice-state placeholders, accepted for interface uniformity and unused here
- `inpc_log_shift`, `w`, `p`: (optional, default `0`) INP-concentration log
  shift, vertical velocity (m/s), and air pressure (Pa)

# Returns
`NamedTuple` with warm rain tendency fields:
- `dq_lcl_dt`: Cloud liquid tendency (kg/kg/s)
- `dn_lcl_dt`: Cloud number tendency (1/kg/s)
- `dq_rai_dt`: Rain tendency (kg/kg/s)
- `dn_rai_dt`: Rain number tendency (1/kg/s)
- `dq_ice_dt`: Ice tendency (always zero for warm-only)
- `dq_rim_dt`: Rime mass tendency (always zero for warm-only)
- `db_rim_dt`: Rime volume tendency (always zero for warm-only)
- `dn_lcl_activation_dt`: Cloud number activation tendency (1/kg/s)
"""
@inline function bulk_microphysics_tendencies(  # TODO: Delete this function
    ::Microphysics2Moment, mp::CMP.Microphysics2MParams{WR, Nothing}, tps,
    ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai,
    q_ice = zero(ρ), n_ice = zero(ρ), q_rim = zero(ρ), b_rim = zero(ρ), logλ = zero(ρ),
    inpc_log_shift = zero(ρ),
    w = zero(ρ), p = zero(ρ),
) where {WR}
    # Clamp negative inputs to zero (robustness against numerical errors)
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

    # Initialize ice-related tendencies (always zero for warm-only)
    dq_ice_dt = zero(ρ)
    dq_rim_dt = zero(ρ)
    db_rim_dt = zero(ρ)

    # --- Warm Rain Processes
    warm = warm_rain_tendencies_2m(mp.warm_rain, tps, T, q_tot, q_lcl, q_rai, q_ice, ρ, n_lcl, n_rai, w, p)
    dq_lcl_dt = warm.dq_lcl_dt
    dn_lcl_dt = warm.dn_lcl_dt
    dq_rai_dt = warm.dq_rai_dt
    dn_rai_dt = warm.dn_rai_dt
    dn_lcl_activation_dt = warm.dn_lcl_activation_dt

    return (; dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt,
        dq_ice_dt, dq_rim_dt, db_rim_dt, dn_lcl_activation_dt)
end


include("BMT_rosenbrock_core.jl")
include("BMT_2mp3.jl")
include("BMT_2mp3_march.jl")

end # module BulkMicrophysicsTendencies
