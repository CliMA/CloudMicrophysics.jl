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
import StaticArrays as SA
import UnrolledUtilities as UU

export MicrophysicsScheme,
    Microphysics0Moment,
    Microphysics1Moment,
    Microphysics2Moment,
    TendencyMode,
    Instantaneous,
    InstantaneousVerbose,
    LinearizedAverage,
    LinearizedAverageVerbose,
    bulk_microphysics_tendencies

#####
##### Singleton types for dispatch
#####

"""
    MicrophysicsScheme

Abstract type for microphysics scheme dispatch.
"""
abstract type MicrophysicsScheme end
Base.broadcastable(x::MicrophysicsScheme) = tuple(x)

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
Base.broadcastable(x::TendencyMode) = tuple(x)

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

"""
    LinearizedAverageVerbose <: TendencyMode

Return the time-averaged tendencies of [`LinearizedAverage`](@ref) and the time-averaged
rate of each process, with the `S_*` names of [`InstantaneousVerbose`](@ref).
"""
struct LinearizedAverageVerbose <: TendencyMode end

include("bulk_process_terms.jl")

"""
    Condensates1M{FT} <: StaticArrays.FieldVector{4, FT}

Specific contents of the four 1-moment condensate species [kg/kg].

The fields are `q_lcl`, `q_icl`, `q_rai` and `q_sno`. The state is indexed by position or
by species name, for example `q[:q_lcl]`.
"""
struct Condensates1M{FT} <: SA.FieldVector{4, FT}
    q_lcl::FT
    q_icl::FT
    q_rai::FT
    q_sno::FT
end
SA.similar_type(::Type{<:Condensates1M}, ::Type{FT}, ::SA.Size{(4,)}) where {FT} = Condensates1M{FT}

@inline Base.getindex(q::Condensates1M, species::Symbol) = getproperty(q, species)
@inline Base.setindex(q::Q, v, species::Symbol) where {Q <: Condensates1M} =
    Base.setindex(q, v, Base.fieldindex(Q, species))

include("bulk_donor_step.jl")

# --- 1-Moment Microphysics ---

# --- Internal helpers ---

"""
    _microphysics_source_terms(
        ::Microphysics1Moment, mp, tps, ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )

Compute the term of each 1-moment microphysics process.

Each term is a [`Transfer`](@ref), [`VaporExchange`](@ref) or [`VaporRelaxation`](@ref), and
its type gives the species of the process and the direction of its rate. A collision with a
cold and a warm branch has one term for each branch, with the suffix `_cold` or `_warm`, and
the term of the inactive branch is zero.

# Arguments
- `mp`: 1-moment microphysics parameters.
- `tps`: thermodynamics parameters.
- `ρ`: air density [kg/m³].
- `T`: temperature [K].
- `w`: vertical velocity [m/s].
- `q_tot`, `q_lcl`, `q_icl`, `q_rai`, `q_sno`: specific contents of total water, cloud
  liquid, cloud ice, rain and snow [kg/kg].

# Returns
- `NamedTuple` of process terms. The terms of the two cloud phase changes store their
  relaxation timescale `τ` [s], which is `Inf` if the process is disabled.
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
    # Relaxation timescales (without Γ); `Inf` when disabled
    S_phase_change_vap_lcl = VaporRelaxation(:q_lcl,
        CMNonEq.conv_q_vap_to_q_lcl(procs.cloud_liquid_formation, mp, tps, micro, thermo),
        CMNonEq.τ_vap_to_q_lcl(procs.cloud_liquid_formation, mp, tps, micro, thermo),
    )

    S_phase_change_vap_icl = VaporRelaxation(:q_icl,
        CMNonEq.conv_q_vap_to_q_icl(procs.cloud_ice_formation, mp, tps, micro, thermo),
        CMNonEq.τ_vap_to_q_icl(procs.cloud_ice_formation, mp, tps, micro, thermo),
    )

    # --- Autoconversion (cloud → precipitation, ≥ 0) ---
    S_acnv_lcl_rai =
        Transfer(:q_lcl => :q_rai, CM1.conv_q_lcl_to_q_rai(procs.rain_autoconversion, mp, tps, micro, thermo))

    S_acnv_icl_sno =
        Transfer(:q_icl => :q_sno, CM1.conv_q_icl_to_q_sno(procs.snow_autoconversion, mp, tps, micro, thermo, sd))

    # --- Accretion (collisions between species) ---
    is_warm = T >= TDI.T_freeze(tps)

    # Cloud liquid + rain → rain
    S_accr_lcl_rai =
        Transfer(:q_lcl => :q_rai, CM1.accretion(procs.cloud_liquid_rain_accretion, mp, tps, micro, thermo, sd))

    # Cloud liquid + snow: product goes to sno (cold) or rai (warm), plus thermal melt.
    # The previous accretion versions return a single number, this one has to return a tuple.
    (; S_accr, S_melt) = if isnothing(procs.cloud_liquid_snow_accretion)
        (; S_accr = zero(FT), S_melt = zero(FT))
    else
        CM1.accretion(procs.cloud_liquid_snow_accretion, mp, tps, micro, thermo, sd)
    end
    S_accr_lcl_sno_cold = Transfer(:q_lcl => :q_sno, ifelse(is_warm, zero(FT), S_accr)) # cold
    S_accr_lcl_sno_warm = Transfer(:q_lcl => :q_rai, ifelse(is_warm, S_accr, zero(FT))) # warm
    S_accr_melt_lcl_sno = Transfer(:q_sno => :q_rai, S_melt) # thermal melt, zero when cold

    # Cloud ice + rain → snow (ice-side sink)
    S_accr_icl_rai =
        Transfer(:q_icl => :q_sno, CM1.accretion(procs.cloud_ice_rain_accretion, mp, tps, micro, thermo, sd))

    # Rain frozen in cloud ice + rain collision → snow (rain sink)
    S_accr_freeze_icl_rai =
        Transfer(:q_rai => :q_sno, CM1.accretion_rain_sink(procs.cloud_ice_rain_accretion, mp, tps, micro, thermo, sd))

    # Cloud ice + snow → snow
    S_accr_icl_sno =
        Transfer(:q_icl => :q_sno, CM1.accretion(procs.cloud_ice_snow_accretion, mp, tps, micro, thermo, sd))

    # Rain-snow collisions: split into cold/warm arms (inactive arm = zero)
    (; S_rai_sno, S_sno_rai, S_melt) = CM1.accretion_snow_rain(procs.rain_snow_accretion, mp, tps, micro, thermo, sd)
    S_accr_rai_sno_cold = Transfer(:q_rai => :q_sno, ifelse(is_warm, zero(FT), S_rai_sno)) # cold
    S_accr_rai_sno_warm = Transfer(:q_sno => :q_rai, ifelse(is_warm, S_sno_rai, zero(FT))) # warm
    S_accr_melt_rai_sno = Transfer(:q_sno => :q_rai, ifelse(is_warm, S_melt, zero(FT)))    # thermal melt

    # --- Phase change: precipitation ↔ vapor ---
    S_phase_change_vap_rai =
        VaporExchange(:q_rai, CM1.conv_q_rai_to_q_vap(procs.rain_condensation_evaporation, mp, tps, micro, thermo, sd))
    S_phase_change_vap_sno =
        VaporExchange(:q_sno, CM1.conv_q_sno_to_q_vap(procs.snow_deposition_sublimation, mp, tps, micro, thermo, sd))

    # --- Melting ---
    S_melt_icl_lcl =
        Transfer(:q_icl => :q_lcl, CM1.conv_q_icl_to_q_lcl(procs.cloud_ice_melt, mp, tps, micro, thermo, sd))
    S_melt_sno_rai =
        Transfer(:q_sno => :q_rai, CM1.conv_q_sno_to_q_rai(procs.snow_melt, mp, tps, micro, thermo, sd))

    # --- Freezing ---
    S_freeze_lcl_icl =
        Transfer(:q_lcl => :q_icl, CM1.conv_q_lcl_to_q_icl(procs.cloud_liquid_freezing, mp, tps, micro, thermo))

    return (;
        S_phase_change_vap_lcl, S_phase_change_vap_icl,
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
    linearized_step_1m(mp, tps, ρ, T, w, q_tot, q, Δt)

Compute one linearized substep of width `Δt` for the 1-moment species `q`.

The substep returns the time-averaged tendencies of the species and the rate of each
process over the substep.

The vapor sources are scaled by a common factor `α ≤ 1`, so that the substep does not reduce
the vapor below the smaller of its saturation contents over liquid and over ice.

# Arguments
- `q`: species, as a [`Condensates1M`](@ref) [kg/kg].
- `Δt`: width of the substep [s].
- The other arguments are those of `_microphysics_source_terms`.

# Returns
- `(; dq_dt, rates)`: the tendencies, as a [`Condensates1M`](@ref) [kg/kg/s], and the rate
  of each process from [`donor_rates`](@ref) [kg/kg/s].
"""
@inline function linearized_step_1m(mp::CMP.Microphysics1MParams, tps, ρ, T, w, q_tot, q, Δt)
    FT = eltype(q)
    terms = _microphysics_source_terms(Microphysics1Moment(), mp, tps, ρ, T, w, q_tot, q...)
    q_min = TDI.TD.Parameters.q_min(tps)
    lin = donor_linearization(terms, q, q_min, Δt)

    invΔt = one(FT) / Δt
    # Cap vap→condensate sources jointly so the substep cannot drive
    # `q_v` below `min(q_sat_liq, q_sat_ice)`. Preserves relative rates.
    q_sat_min = min(
        TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ),
        TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ),
    )
    q_v = q_tot - q.q_lcl - q.q_icl - q.q_rai - q.q_sno
    α = min(one(FT), max(zero(FT), q_v - q_sat_min) * invΔt / max(sum(lin.e), eps(FT)))

    q_new = backward_euler_solve(lin, q, α, Δt)
    return (; dq_dt = (q_new - q) * invΔt, rates = donor_rates(terms, q, q_new, α, q_min, Δt))
end

"""
    _linearized_implicit_step(
        ::Microphysics1Moment, mp, tps, ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
    )

Compute the tendencies of `linearized_step_1m` as the `NamedTuple`
`(; dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt)` [kg/kg/s].
"""
@inline function _linearized_implicit_step(
    ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt::AbstractFloat,
)
    q = Condensates1M(q_lcl, q_icl, q_rai, q_sno)
    (; dq_dt, rates) = linearized_step_1m(mp, tps, ρ, T, w, q_tot, q, Δt)
    return _output_1m(LinearizedAverage(), dq_dt, rates)
end

"""
    _substep_1m(mp, tps, ρ, T, w, q_tot, q, Δt_sub)

Advance the species `q` and the temperature `T` over one substep of width `Δt_sub`, and
return them as `(; q, T, increments)` with the increment of each process [kg/kg].

The temperature changes by the latent heat of the species increments.
"""
@inline function _substep_1m(mp, tps, ρ, T, w, q_tot, q, Δt_sub)
    cp_d = TDI.TD.Parameters.cp_d(tps)
    Lv = TDI.TD.Parameters.LH_v0(tps)
    Ls = TDI.TD.Parameters.LH_s0(tps)
    (; dq_dt, rates) = linearized_step_1m(mp, tps, ρ, T, w, q_tot, q, Δt_sub)
    (; q_lcl, q_rai, q_icl, q_sno) = dq_dt
    q += dq_dt * Δt_sub
    T += (Lv * (q_lcl + q_rai) + Ls * (q_icl + q_sno)) / cp_d * Δt_sub
    return (; q, T, increments = UU.unrolled_map(r -> r * Δt_sub, rates))
end

"""
    _zero_increments_1m(mp, tps, ρ, T, w, q_tot, q, Δt_sub)

Return zero increments of the processes, of the type that `_substep_1m` returns for the same arguments.
"""
@inline function _zero_increments_1m(args...)
    Substep = Core.Compiler.return_type(_substep_1m, typeof(args))
    Increments = fieldtype(Substep, :increments)
    return Increments(ntuple(_ -> zero(eltype(Increments)), Val(fieldcount(Increments))))
end

"""
    _output_1m(mode, dq_dt, rates)

Return the output of `bulk_microphysics_tendencies` for the 1-moment tendencies `dq_dt`, a
[`Condensates1M`](@ref) [kg/kg/s], and the rate of each process `rates` [kg/kg/s].

The output has the tendencies `dq_lcl_dt`, `dq_icl_dt`, `dq_rai_dt` and `dq_sno_dt`, and,
with `InstantaneousVerbose` and `LinearizedAverageVerbose`, also the rates.
"""
@inline function _output_1m(mode, dq_dt::Condensates1M, rates)
    tendencies = (;
        dq_lcl_dt = dq_dt.q_lcl, dq_icl_dt = dq_dt.q_icl, dq_rai_dt = dq_dt.q_rai, dq_sno_dt = dq_dt.q_sno,
    )
    return mode isa Union{InstantaneousVerbose, LinearizedAverageVerbose} ? merge(tendencies, rates) : tendencies
end

# --- Public API: bulk_microphysics_tendencies with TendencyMode dispatch ---

"""
    bulk_microphysics_tendencies(
        ::Union{Instantaneous, InstantaneousVerbose}, ::Microphysics1Moment, mp, tps,
        ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )

Compute all 1-moment microphysics tendencies in one fused call.

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
- With `InstantaneousVerbose`, also the rate of each process, such as `S_acnv_lcl_rai`
  [kg/kg/s], positive in the direction of the process.

# Notes
- Negative specific contents are clamped to zero for robustness.
- Does NOT apply timestep-dependent limiters.
"""
@inline function bulk_microphysics_tendencies(
    mode::Union{Instantaneous, InstantaneousVerbose}, ::Microphysics1Moment,
    mp::CMP.Microphysics1MParams, tps, ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    terms = _microphysics_source_terms(
        Microphysics1Moment(), mp, tps,
        ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )
    dq_dt = species_tendency(terms, Condensates1M{typeof(q_tot)})
    return _output_1m(mode, dq_dt, UU.unrolled_map(t -> t.S, terms))
end

"""
    bulk_microphysics_tendencies(
        ::Union{LinearizedAverage, LinearizedAverageVerbose}, ::Microphysics1Moment, mp, tps,
        ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1,
    )

Compute average 1-moment microphysics tendencies over `Δt` using repeated linearized
implicit substeps.

The interval `Δt` is divided into `nsub` equal substeps. At each substep, the
linearized tendency is rebuilt from the current state and solved implicitly for cloud
liquid, cloud ice, rain, and snow. Temperature is then updated from the latent
heating implied by the substep tendencies. Increasing `nsub` improves how well the method
captures nonlinear changes in the active microphysical processes, including regime changes
near freezing.

# Returns
- `NamedTuple` with the net change of each species over `Δt` divided by `Δt` [kg/kg/s]:
  `dq_lcl_dt`, `dq_icl_dt`, `dq_rai_dt` and `dq_sno_dt`.
- With `LinearizedAverageVerbose`, also the average rate of each process over `Δt`, with
  the names of `InstantaneousVerbose` [kg/kg/s].
"""
@inline function bulk_microphysics_tendencies(
    mode::Union{LinearizedAverage, LinearizedAverageVerbose}, ::Microphysics1Moment,
    mp::CMP.Microphysics1MParams, tps, ρ, T, w,
    q_tot, q_lcl, q_icl, q_rai, q_sno,
    Δt::AbstractFloat, nsub::Integer = 1,
)
    FT = typeof(q_tot)
    q₀ = Condensates1M(q_lcl, q_icl, q_rai, q_sno)
    Δt_sub = Δt / FT(nsub)

    q = q₀
    increments = _zero_increments_1m(mp, tps, ρ, T, w, q_tot, q₀, Δt_sub)
    for _ in 1:nsub
        substep = _substep_1m(mp, tps, ρ, T, w, q_tot, q, Δt_sub)
        (; q, T) = substep
        increments = UU.unrolled_map(+, increments, substep.increments)
    end
    dq_dt = (q - q₀) / Δt
    return _output_1m(mode, dq_dt, UU.unrolled_map(Δ -> Δ / Δt, increments))
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
    warm_rain_tendencies_2m(
        warm_rain, tps, T, q_tot, q_lcl, q_rai, q_ice, ρ, n_lcl, n_rai,
        w = zero(ρ), p = zero(ρ),
    )

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
    condevap = warm_rain.condevap

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
    micro_mock = (; q_tot, q_lcl, q_icl = q_ice, q_rai, q_sno = zero(q_ice))
    thermo_mock = (; ρ, T)
    ∂ₜq_lcl_cond = CMNonEq._conv_q_vap_to_q_lcl_const(
        condevap.τ_relax, tps, micro_mock, thermo_mock,
    )
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
    ∂ₜN_lcl_sc = CM2.cloud_liquid_self_collection(sb.acnv, sb.pdf_c, q_lcl, ρ, acnv.dN_lcl_dt)
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

"""
    bulk_microphysics_tendencies(
        ::Microphysics2Moment,
        mp::Microphysics2MParams{WR, <:P3IceParams}, tps,
        ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai,
        q_ice, n_ice, q_rim, b_rim, logλ,
        inpc_log_shift = 0, w = 0, p = 0,
    )

Compute 2-moment **warm rain + P3 ice** microphysics tendencies.

This method is type-stable and GPU-optimized. The P3 ice parameters are guaranteed
to be non-Nothing, eliminating runtime type checks and dynamic dispatch.

# Arguments
## Required
- `mp`: Microphysics2MParams with P3 ice parameters present
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
- `logλ`: Log of P3 distribution slope parameter, log(1/m)

## Optional
- `inpc_log_shift`: Additive shift to log(INPC) (default `0`)
- `w`: Vertical velocity (m/s), default `0`
- `p`: Air pressure (Pa), default `0`

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
- `dn_lcl_activation_dt`: Cloud number activation tendency (1/kg/s)
"""
@inline function bulk_microphysics_tendencies(
    ::Microphysics2Moment, mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai,
    q_ice, n_ice, q_rim, b_rim, logλ,
    inpc_log_shift = zero(ρ),
    w = zero(ρ), p = zero(ρ),
) where {WR, ICE <: CMP.P3IceParams}
    FT = eltype(ρ)
    ϵₘ = UT.ϵ_numerics_2M_M(FT)
    ϵₙ = UT.ϵ_numerics_2M_N(FT)
    ϵB = UT.ϵ_numerics_P3_B(FT)
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

    # Convert to volumetric quantities for P3 functions
    L_lcl = q_lcl * ρ  # [kg lcl / m³ air]
    L_rai = q_rai * ρ  # [kg rai / m³ air]
    N_lcl = n_lcl * ρ  # [1 / m³ air]
    N_rai = n_rai * ρ  # [1 / m³ air]
    L_ice = q_ice * ρ  # [kg ice / m³ air]
    N_ice = n_ice * ρ  # [1 / m³ air]
    L_rim = q_rim * ρ  # [kg rim / m³ air]
    B_rim = b_rim * ρ  # [m³ rim / m³ air]
    state = CMP3.state_from_prognostic(mp.ice.scheme, L_ice, N_ice, L_rim, B_rim)

    # Unpack warm rain parameters
    aps = mp.warm_rain.air_properties
    subdep = mp.warm_rain.subdep

    # Initialize ice-related tendencies
    dq_ice_dt = zero(ρ)
    dn_ice_dt = zero(ρ)
    dq_rim_dt = zero(ρ)
    db_rim_dt = zero(ρ)

    # --- Warm Rain Processes
    warm = warm_rain_tendencies_2m(mp.warm_rain, tps, T, q_tot, q_lcl, q_rai, q_ice, ρ, n_lcl, n_rai, w, p)
    dq_lcl_dt = warm.dq_lcl_dt
    dn_lcl_dt = warm.dn_lcl_dt
    dq_rai_dt = warm.dq_rai_dt
    dn_rai_dt = warm.dn_rai_dt
    dn_lcl_activation_dt = warm.dn_lcl_activation_dt

    # --- P3 Ice Processes
    p3 = mp.ice.scheme
    vel = mp.ice.terminal_velocity
    pdf_c = mp.ice.cloud_pdf
    pdf_r = mp.ice.rain_pdf
    ice_nucleation = mp.ice.ice_nucleation
    inp_depletion_model = mp.ice.inp_depletion_model
    quad = mp.ice.quad

    # Only compute ice processes if there is ice mass/number present
    if q_ice > ϵₘ && n_ice > ϵₙ

        # --- Liquid-ice collisions
        coll = CMP3.bulk_liquid_ice_collision_sources(
            state, logλ, pdf_c, pdf_r, L_lcl, N_lcl, L_rai, N_rai, aps, tps, vel, ρ, T;
            quad,
        )
        dq_lcl_dt += coll.∂ₜq_c
        dq_rai_dt += coll.∂ₜq_r
        dn_lcl_dt += coll.∂ₜN_c / ρ
        dn_rai_dt += coll.∂ₜN_r / ρ
        dq_ice_dt += coll.∂ₜL_ice / ρ
        dq_rim_dt += coll.∂ₜL_rim / ρ
        db_rim_dt += coll.∂ₜB_rim / ρ

        # --- Ice self-collection (aggregation)
        S_ice_agg = CMP3.ice_self_collection(state, logλ, vel, ρ; quad)
        dn_ice_dt -= S_ice_agg.dNdt / ρ

        # Ice melting (above freezing temperature)
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        melt = ifelse(T > T_freeze,
            CMP3.ice_melt(vel, aps, tps, T, ρ, state, logλ; quad),
            (; dNdt = zero(ρ), dLdt = zero(ρ)),
        )
        # Specific (per-kg-air) ice-mass melt rate.
        ∂ₜq_ice_melt = melt.dLdt / ρ
        ∂ₜn_ice_melt = melt.dNdt / ρ
        # Melting converts ice to rain.
        dq_rai_dt += ∂ₜq_ice_melt
        dn_rai_dt += ∂ₜn_ice_melt  # Melted ice becomes rain drops
        dq_ice_dt -= ∂ₜq_ice_melt
        dn_ice_dt -= ∂ₜn_ice_melt  # Ice particles consumed by melting
        # Rim mass and rim volume drain proportionally to ice mass during melting
        dq_rim_dt -= ∂ₜq_ice_melt * state.F_rim
        db_rim_dt -= ifelse(state.ρ_rim > 0, ∂ₜq_ice_melt * state.F_rim / state.ρ_rim, zero(FT))
    end

    # --- Ice nucleation (F23 + Bigg)
    τ_act = inp_depletion_model.τ_act
    # Vapor deposition nucleation size. TODO: put into ClimaParams.
    D_nuc = FT(10e-6)  # 10 μm nascent crystal - small-D tail of the P3
    m_nuc = p3.ρ_i * CO.volume_sphere_D(D_nuc)

    # F23 INP-activation depletion proxy.
    n_active = CM_HetIce.n_active(inp_depletion_model, n_ice)

    # --- deposition nucleation (vapor → pristine ice)
    dep = CM_HetIce.deposition_rate(
        ice_nucleation, tps, T, ρ, q_tot, q_lcl + q_rai, q_ice, n_active;
        m_nuc, τ_act, inpc_log_shift,
    )

    dn_ice_dt += dep.∂ₜn_frz
    dq_ice_dt += dep.∂ₜq_frz
    # No contribution to q_rim, b_rim — pristine deposition crystals have F_rim = 0.

    # --- F23-bounded Bigg immersion freezing of cloud drops
    cld_bigg = CM_HetIce.liquid_freezing_rate(
        mp.ice.rain_freezing, pdf_c, tps, q_lcl, ρ, N_lcl, T,
    )
    cld_cap = CM_HetIce.immersion_limit_rate(
        ice_nucleation, T, ρ; τ = τ_act, inpc_log_shift, n_active,
    )
    ∂ₜn_imm = min(cld_bigg.∂ₜn_frz, cld_cap.∂ₜn_frz)
    ∂ₜq_imm = ifelse(cld_bigg.∂ₜn_frz > 0, cld_bigg.∂ₜq_frz * ∂ₜn_imm / cld_bigg.∂ₜn_frz, zero(FT))

    # Drain liquid:
    dq_lcl_dt -= ∂ₜq_imm
    dn_lcl_dt -= ∂ₜn_imm
    # Add to ice as fully-rimed embryo graupel:
    dq_ice_dt += ∂ₜq_imm
    dn_ice_dt += ∂ₜn_imm
    dq_rim_dt += ∂ₜq_imm           # F_rim = 1 (frozen drop)
    db_rim_dt += ∂ₜq_imm / p3.ρ_i  # solid-ice rime volume

    # --- Ice Sublimation / Deposition
    n_per_q_ice = ifelse(q_ice > ϵₘ, n_ice / q_ice, zero(n_ice))
    # Deposition/sublimation of cloud ice
    micro_mock = (; q_tot, q_lcl, q_icl = q_ice, q_rai, q_sno = zero(q_ice))
    thermo_mock = (; ρ, T)
    ∂ₜq_ice_dep = CMNonEq._conv_q_vap_to_q_icl_const(
        subdep.τ_relax, tps, micro_mock, thermo_mock,
    )
    # No ice deposition above freezing (lack of INPs)
    ∂ₜq_ice_dep = ifelse(T > tps.T_freeze, min(∂ₜq_ice_dep, zero(T)), ∂ₜq_ice_dep)
    # During sublimation, the number of ice particles decreases in proportion to the mean ice mass
    # During deposition, the number of ice particles remain unchanged
    ∂ₜn_ice_dep = ifelse(∂ₜq_ice_dep < 0, n_per_q_ice * ∂ₜq_ice_dep, zero(∂ₜq_ice_dep))
    dq_ice_dt += ∂ₜq_ice_dep
    dn_ice_dt += ∂ₜn_ice_dep
    ∂ₜq_ice_sub = min(∂ₜq_ice_dep, 0)   # ≤ 0; zero on the deposition branch
    dq_rim_dt += ∂ₜq_ice_sub * state.F_rim
    db_rim_dt += ifelse(state.ρ_rim > 0, ∂ₜq_ice_sub * state.F_rim / state.ρ_rim, zero(FT))

    # --- Ice number adjustment for mass limits
    # Nudges n_ice toward [q_ice / x_max, q_ice / x_min] over timescale τ.
    numadj = (;  # TODO: put into ClimaParams
        τ = FT(100),
        x_min = FT(1e-12),  # min mean ice particle mass [kg] (~10 μm crystal)
        x_max = FT(1e-5),   # max mean ice particle mass [kg] (~5 mm aggregate)
    )
    ∂ₜn_ice_numadj = CM2.number_tendency_from_mass_limits(numadj, q_ice, n_ice)
    dn_ice_dt += ∂ₜn_ice_numadj

    # --- Rain Heterogeneous Freezing (Bigg 1953)
    rain_frz = CM_HetIce.liquid_freezing_rate(mp.ice.rain_freezing, pdf_r, tps, q_rai, ρ, N_rai, T)

    # Rain → ice (frozen rain is fully rimed, per MM15)
    dq_rai_dt -= rain_frz.∂ₜq_frz
    dn_rai_dt -= rain_frz.∂ₜn_frz
    dq_ice_dt += rain_frz.∂ₜq_frz
    dn_ice_dt += rain_frz.∂ₜn_frz
    dq_rim_dt += rain_frz.∂ₜq_frz
    db_rim_dt += rain_frz.∂ₜq_frz / p3.ρ_i  # ρ_i = 916.7 kg m⁻³, the density of solid bulk ice

    # Aerosol activation is folded into `warm_rain_tendencies_2m` above —
    # `dn_lcl_activation_dt` from `warm` is already included in `dn_lcl_dt`.

    return (; dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt,
        dq_ice_dt, dn_ice_dt, dq_rim_dt, db_rim_dt,
        dn_lcl_activation_dt)
end

end # module BulkMicrophysicsTendencies
