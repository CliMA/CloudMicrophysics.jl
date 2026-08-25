#####
##### 1M Rosenbrock-Euler substepping (`RosenbrockAverage`)
#####

"""
    MicroState1M{FT}

The four prognostic 1M species as a `StaticArrays.FieldVector`. Internal to the
[`RosenbrockAverage`](@ref) implementation, mirroring [`MicroState2MP3`](@ref)
at dimension 4 (1M carries no number species).
"""
struct MicroState1M{FT} <: SA.FieldVector{4, FT}
    q_lcl::FT
    q_icl::FT
    q_rai::FT
    q_sno::FT
end
SA.similar_type(::Type{<:MicroState1M}, ::Type{FT}, ::SA.Size{(4,)}) where {FT} =
    MicroState1M{FT}

"""
    _condensate_phases(x::MicroState1M)

The condensed water of a [`MicroState1M`](@ref) split by phase: cloud liquid and rain
together, then cloud ice and snow together.
"""
@inline _condensate_phases(x::MicroState1M) = (x.q_lcl + x.q_rai, x.q_icl + x.q_sno)

"""
    _instantaneous_1m_tendency(mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno)

The raw instantaneous 1M tendency projected onto the four prognostic species:
the unlimited process rates of the `Microphysics1Moment` `Instantaneous` entry
(`_microphysics_source_terms` aggregated by `_aggregate_tendencies`), without
timestep-dependent clipping.
"""
@inline function _instantaneous_1m_tendency(mp, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    src = _microphysics_source_terms(Microphysics1Moment(), mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )
    return _aggregate_tendencies(src)
end

"""
    Raw1MTendency(mp, tps, ρ, T, q_tot)

Callable bundling the frozen per-substep context; applying it to the species
vector evaluates [`_instantaneous_1m_tendency`](@ref), mirroring
[`Instantaneous2MP3Tendency`](@ref). `q_tot` is promoted to the state's element
type at the call; `T` and `ρ` stay plain.
"""
struct Raw1MTendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    T::F
    q_tot::F
end
@inline function (g::Raw1MTendency)(x::SA.StaticVector{4})
    (q_lcl, q_icl, q_rai, q_sno) = x
    tend = _instantaneous_1m_tendency(g.mp, g.tps,
        g.ρ, g.T, eltype(x)(g.q_tot),
        q_lcl, q_icl, q_rai, q_sno,
    )
    return MicroState1M(tend.dq_lcl_dt, tend.dq_icl_dt, tend.dq_rai_dt, tend.dq_sno_dt)
end

"""
    _per_process_1m(src)

Project each individual 1M source term onto the four prognostic species as a
[`MicroState1M`](@ref) `(q_lcl, q_icl, q_rai, q_sno)`, returning a `NamedTuple`
of these per-process contribution vectors. The signs match
[`_aggregate_tendencies`](@ref), so the sum over the returned processes equals
the aggregated raw tendency the [`Raw1MTendency`](@ref) functor evaluates. Used
by the verbose post-solve attribution to supply the per-process right-hand sides
`f_p`.
"""
@inline function _per_process_1m(src)
    FT = typeof(src.S_phase_change_vap_lcl)
    o = zero(FT)
    return (;
        phase_change_vap_lcl = MicroState1M(src.S_phase_change_vap_lcl, o, o, o),
        phase_change_vap_icl = MicroState1M(o, src.S_phase_change_vap_icl, o, o),
        acnv_lcl_rai = MicroState1M(-src.S_acnv_lcl_rai, o, src.S_acnv_lcl_rai, o),
        acnv_icl_sno = MicroState1M(o, -src.S_acnv_icl_sno, o, src.S_acnv_icl_sno),
        accr_lcl_rai = MicroState1M(-src.S_accr_lcl_rai, o, src.S_accr_lcl_rai, o),
        accr_lcl_sno_cold = MicroState1M(-src.S_accr_lcl_sno_cold, o, o, src.S_accr_lcl_sno_cold),
        accr_lcl_sno_warm = MicroState1M(-src.S_accr_lcl_sno_warm, o, src.S_accr_lcl_sno_warm, o),
        accr_melt_lcl_sno = MicroState1M(o, o, src.S_accr_melt_lcl_sno, -src.S_accr_melt_lcl_sno),
        accr_icl_rai = MicroState1M(o, -src.S_accr_icl_rai, o, src.S_accr_icl_rai),
        accr_freeze_icl_rai = MicroState1M(o, o, -src.S_accr_freeze_icl_rai, src.S_accr_freeze_icl_rai),
        accr_icl_sno = MicroState1M(o, -src.S_accr_icl_sno, o, src.S_accr_icl_sno),
        accr_rai_sno_cold = MicroState1M(o, o, -src.S_accr_rai_sno_cold, src.S_accr_rai_sno_cold),
        accr_rai_sno_warm = MicroState1M(o, o, src.S_accr_rai_sno_warm, -src.S_accr_rai_sno_warm),
        accr_melt_rai_sno = MicroState1M(o, o, src.S_accr_melt_rai_sno, -src.S_accr_melt_rai_sno),
        phase_change_vap_rai = MicroState1M(o, o, src.S_phase_change_vap_rai, o),
        phase_change_vap_sno = MicroState1M(o, o, o, src.S_phase_change_vap_sno),
        melt_icl_lcl = MicroState1M(src.S_melt_icl_lcl, -src.S_melt_icl_lcl, o, o),
        melt_sno_rai = MicroState1M(o, o, src.S_melt_sno_rai, -src.S_melt_sno_rai),
        freeze_lcl_icl = MicroState1M(-src.S_freeze_lcl_icl, src.S_freeze_lcl_icl, o, o),
    )
end

"""
    Verbose1MTendency(mp, tps, ρ, T, q_tot)

Per-process companion to [`Raw1MTendency`](@ref): applying it to the species
vector returns a `NamedTuple` of per-process tendency contributions (each a
[`MicroState1M`](@ref)) via [`_per_process_1m`](@ref), instead of only their
sum. Evaluated at the primal state only, it supplies the right-hand sides `f_p`
for the linear post-solve attribution.
"""
struct Verbose1MTendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    T::F
    q_tot::F
end
@inline function (g::Verbose1MTendency)(x::SA.StaticVector{4, FT}) where {FT}
    (q_lcl, q_icl, q_rai, q_sno) = x
    src = _microphysics_source_terms(Microphysics1Moment(), g.mp, g.tps,
        g.ρ, g.T, FT(g.q_tot),
        q_lcl, q_icl, q_rai, q_sno,
    )
    return _per_process_1m(src)
end

"""
    _rosenbrock_species_mask(x::MicroState1M)

Diagonal of the species projection `P` for the 1M state: 1 for active species,
0 for near-empty ones (mass below `1e-10`, per species: liquid, ice, rain,
snow). See the [`MicroState2MP3`](@ref) method for the role of `P` in
[`_rosenbrock_update_diag`](@ref).
"""
@inline function _rosenbrock_species_mask(x::MicroState1M{FT}) where {FT}
    ϵ_empty = FT(1e-10)
    lcl = ifelse(x.q_lcl < ϵ_empty, zero(FT), one(FT))
    icl = ifelse(x.q_icl < ϵ_empty, zero(FT), one(FT))
    rai = ifelse(x.q_rai < ϵ_empty, zero(FT), one(FT))
    sno = ifelse(x.q_sno < ϵ_empty, zero(FT), one(FT))
    return MicroState1M(lcl, icl, rai, sno)
end

"""
    _apply_limiter(::EndStateSaturationAdjustment, x::MicroState1M, d, ρ, Tsub,
        q_tot, Lv_over_cp, Ls_over_cp, tps)

[`EndStateSaturationAdjustment`](@ref) on the 1M state: the supersaturation
`Ssat` reads the cloud and rain mass as liquid and the cloud ice and snow mass
as ice, and the latent heating of an increment weights the same two groups by
`Lv_over_cp` and `Ls_over_cp`. The scaling itself is
[`_saturation_bisection`](@ref)'s.
"""
@inline function _apply_limiter(::EndStateSaturationAdjustment,
    x::MicroState1M{FT}, d::MicroState1M{FT},
    ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps,
) where {FT}
    Ssat(xx, TT) = max(
        TDI.supersaturation_over_ice(tps, q_tot, xx.q_lcl + xx.q_rai, xx.q_icl + xx.q_sno, ρ, TT),
        TDI.supersaturation_over_liquid(tps, q_tot, xx.q_lcl + xx.q_rai, xx.q_icl + xx.q_sno, ρ, TT),
    )
    latent(dd) = Lv_over_cp * (dd.q_lcl + dd.q_rai) + Ls_over_cp * (dd.q_icl + dd.q_sno)
    return _saturation_bisection(Ssat, latent, x, d, Tsub)
end

"""
    _jacobian_1m(FT; lcl_lcl, lcl_icl, ..., sno_sno)

Assemble the 4×4 Jacobian of the 1M tendency over the state
`(q_lcl, q_icl, q_rai, q_sno)` from named species-pair entries: keyword
`<receiver>_<donor>` is `∂(dq_receiver/dt)/∂q_donor`. Entries not supplied are
zero. Centralizes the index layout so a hand-built Jacobian is filled by physical
coupling name rather than numeric `[i, j]` position.
"""
@inline _jacobian_1m(::Type{FT};
    lcl_lcl = zero(FT), lcl_icl = zero(FT), lcl_rai = zero(FT), lcl_sno = zero(FT),
    icl_lcl = zero(FT), icl_icl = zero(FT), icl_rai = zero(FT), icl_sno = zero(FT),
    rai_lcl = zero(FT), rai_icl = zero(FT), rai_rai = zero(FT), rai_sno = zero(FT),
    sno_lcl = zero(FT), sno_icl = zero(FT), sno_rai = zero(FT), sno_sno = zero(FT),
) where {FT} = SA.SMatrix{4, 4, FT}(
    # column-major: each column is a donor q_*, each row a receiver dq_*/dt
    lcl_lcl, icl_lcl, rai_lcl, sno_lcl,  # ∂/∂q_lcl
    lcl_icl, icl_icl, rai_icl, sno_icl,  # ∂/∂q_icl
    lcl_rai, icl_rai, rai_rai, sno_rai,  # ∂/∂q_rai
    lcl_sno, icl_sno, rai_sno, sno_sno,  # ∂/∂q_sno
)

"""
    _jacobian_1m_linearized(g::Raw1MTendency, x::MicroState1M)

Donor-based Jacobian provider for [`_rosenbrock_average_1m`](@ref): the donor-based
linearized system matrix `M` that [`LinearizedAverage`](@ref) uses (built by
`_linearize`), assembled by name via [`_jacobian_1m`](@ref).

`M` is the donor-based linearization, not the exact derivative of the raw
tendency: each donor to receiver transfer is linearized only in its donor species,
so collector-species couplings are absent. It is evaluated at the same
`(ρ, Tsub, q_tot, x)` as the raw tendency `f = g(x)`.
"""
@inline function _jacobian_1m_linearized(g::Raw1MTendency, x::MicroState1M{FT}) where {FT}
    (; q_lcl, q_icl, q_rai, q_sno) = x
    src = _microphysics_source_terms(
        Microphysics1Moment(), g.mp, g.tps, g.ρ, g.T, FT(g.q_tot),
        q_lcl, q_icl, q_rai, q_sno,
    )
    q_min = TDI.TD.Parameters.q_min(g.tps)
    M = _linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min)
    return _jacobian_1m(FT;
        lcl_lcl = M.M11, lcl_icl = M.M12,
        icl_lcl = M.M21, icl_icl = M.M22,
        rai_lcl = M.M31, rai_rai = M.M33, rai_sno = M.M34,
        sno_lcl = M.M41, sno_icl = M.M42, sno_rai = M.M43, sno_sno = M.M44,
    )
end

"""
    _vapor_exchange_rates(g::Raw1MTendency, x::MicroState1M, q_tot)

The four vapor-to-species phase-change rates `(S_phase_change_vap_lcl,
S_phase_change_vap_icl, S_phase_change_vap_rai, S_phase_change_vap_sno)` of the
1M source terms at state `x` and total water `q_tot`, as a `MicroState1M` in the
`(lcl, icl, rai, sno)` receiver order. Used to differentiate the vapor-exchange
coupling with respect to `q_tot`.
"""
@inline function _vapor_exchange_rates(g::Raw1MTendency, x::MicroState1M, q_tot)
    (; q_lcl, q_icl, q_rai, q_sno) = x
    src = _microphysics_source_terms(
        Microphysics1Moment(), g.mp, g.tps, g.ρ, g.T, q_tot,
        q_lcl, q_icl, q_rai, q_sno,
    )
    return MicroState1M(
        src.S_phase_change_vap_lcl, src.S_phase_change_vap_icl,
        src.S_phase_change_vap_rai, src.S_phase_change_vap_sno,
    )
end

"""
    _jacobian_1m_coupled(g::Raw1MTendency, x::MicroState1M)

Coupled donor-based Jacobian provider for [`_rosenbrock_average_1m`](@ref): the
donor-based matrix [`_jacobian_1m_linearized`](@ref) with the vapor-competition
(Wegener-Bergeron-Findeisen) coupling restored.

The coupling is approximated by its vapor part: each rate depends on the
vapor specific content `q_vap = q_tot - q_lcl - q_icl - q_rai - q_sno`, giving a
condensate derivative `-∂(rate)/∂q_vap`, obtained from the derivative of the
vapor-to-species rates ([`_vapor_exchange_rates`](@ref)) with respect to `q_tot`
and added to each receiver row. The direct condensate dependence of the rates
(for example rain ventilation and the condensate availability terms) is not
recovered; use [`ExactJacobian`](@ref) for the full derivative.
"""
@inline function _jacobian_1m_coupled(g::Raw1MTendency, x::MicroState1M{FT}) where {FT}
    Jdonor = _jacobian_1m_linearized(g, x)
    dS_dq_tot = FD.derivative(qt -> _vapor_exchange_rates(g, x, qt), FT(g.q_tot))
    wbf = -dS_dq_tot
    coupling = SA.SMatrix{4, 4, FT}(
        wbf[1], wbf[2], wbf[3], wbf[4],
        wbf[1], wbf[2], wbf[3], wbf[4],
        wbf[1], wbf[2], wbf[3], wbf[4],
        wbf[1], wbf[2], wbf[3], wbf[4],
    )
    return Jdonor + coupling
end


"""
    _tendency_and_jacobian(::DonorJacobian, g, x)
    _tendency_and_jacobian(::CoupledDonorJacobian, g, x)

The raw substep tendency and the donor-based substep Jacobian at `x`, as
`(f, J)`. The donor-based matrices produce no tendency by-product, so `f = g(x)`
is evaluated separately.
"""
@inline _tendency_and_jacobian(::DonorJacobian, g, x) = (g(x), _jacobian_1m_linearized(g, x))
@inline _tendency_and_jacobian(::CoupledDonorJacobian, g, x) =
    (g(x), _jacobian_1m_coupled(g, x))

"""
    _rosenbrock_average_1m(mode, mp, tps, ρ, T, q_tot,
                           q_lcl, q_icl, q_rai, q_sno, Δt, nsub)

Shared 1M Rosenbrock-average substep driver. The substep loop, equilibration,
positivity clamp, and explicit between-substep `T` update are fixed; the
[`RosenbrockAverage`](@ref) `mode` selects the Jacobian, the growth treatment,
and the increment limiter through [`_tendency_and_jacobian`](@ref),
[`_species_mask`](@ref), [`_apply_growth`](@ref), and [`_apply_limiter`](@ref).
"""
@inline function _rosenbrock_average_1m(
    mode::RosenbrockAverage, mp::CMP.Microphysics1MParams, tps, ρ, T, q_tot,
    q_lcl, q_icl, q_rai, q_sno, Δt, nsub,
)
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)
    mask = _species_mask(mode.jacobian, mode.growth)

    x = MicroState1M{FT}(q_lcl, q_icl, q_rai, q_sno)
    x₀ = x
    Tsub = T
    for _ in 1:nsub_eff
        g = Raw1MTendency(mp, tps, ρ, Tsub, q_tot)
        x_prev = x
        if all(isfinite, x)
            f, J_raw = _tendency_and_jacobian(mode.jacobian, g, x)
            J = _apply_growth(mode.growth, J_raw)
            z = mask(x)
            d = if all(isfinite, J)
                # No rime pair in a one-moment state, so the floor here is the generic
                # `max.(x .+ Δx, 0)` and the density bounds it would take are genuinely
                # unused; the substep is entered directly rather than through a wrapper
                # that defaults them, which the 2M+P3 state cannot use.
                first(_rosenbrock_update_diag(x, f, J, z, h)) - x
            else
                _euler_update(x, f, h) - x
            end
            d = _apply_limiter(mode.limiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps)
            x = max.(x .+ d, 0)
        else
            f = g(x)
            x = _euler_update(x, f, h)
        end
        Δ = x - x_prev
        Tsub += Lv_over_cp * (Δ.q_lcl + Δ.q_rai) + Ls_over_cp * (Δ.q_icl + Δ.q_sno)
    end

    rates = (x - x₀) / Δt
    return (;
        dq_lcl_dt = rates.q_lcl, dq_icl_dt = rates.q_icl,
        dq_rai_dt = rates.q_rai, dq_sno_dt = rates.q_sno,
    )
end

"""
    bulk_microphysics_tendencies(::RosenbrockAverage, ::Microphysics1Moment,
        mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1)

Compute average 1M microphysics tendencies over `Δt` using `nsub`
linearized-implicit (Rosenbrock-Euler) substeps of the raw instantaneous
tendency. The [`Jacobian`](@ref), [`GrowthTreatment`](@ref), and
[`TendencyLimiter`](@ref) options of `mode` select the substep matrix, the
growth-diagonal treatment, and the increment limiter.

# Algorithm

For each substep of `h = Δt / nsub`:

1. Build the substep Jacobian `J` from the [`Jacobian`](@ref) provider and apply
   the [`GrowthTreatment`](@ref) to it.
2. Advance with [`_rosenbrock_update_diag`](@ref): solve `(I/h - P J P) Δx = f` in
   equilibrated variables, where the projection `P` is selected by
   [`_species_mask`](@ref).
3. Apply the [`TendencyLimiter`](@ref) to the increment.
4. Update the local temperature from the latent heating of the realized
   increments (constant latent heats, liquid+rain on `L_v`, ice+snow on `L_s`).

A non-finite state or Jacobian falls back to a forward-Euler substep of the raw
tendency; `q_tot` is held fixed across substeps.

Returns the net change in the species over `Δt` divided by `Δt`, in the same
fields as the `Instantaneous` entry.
"""
@inline bulk_microphysics_tendencies(mode::RosenbrockAverage, ::Microphysics1Moment,
    mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1) =
    _rosenbrock_average_1m(mode, mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub)

"""
    bulk_microphysics_tendencies(
        ::LinearizedAverage, ::Microphysics1Moment, mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1,
    )

Compute average 1-moment microphysics tendencies over `Δt` using repeated
linearized implicit substeps. Forwards to the donor-based configuration of
[`RosenbrockAverage`](@ref): [`DonorJacobian`](@ref) with
[`ImplicitGrowth`](@ref) and [`NoLimiter`](@ref).

# Returns
`NamedTuple` with fields:
- `dq_lcl_dt`: Cloud liquid tendency [kg/kg/s]
- `dq_icl_dt`: Cloud ice tendency [kg/kg/s]
- `dq_rai_dt`: Rain tendency [kg/kg/s]
- `dq_sno_dt`: Snow tendency [kg/kg/s]
"""
@inline bulk_microphysics_tendencies(
    ::LinearizedAverage, cm::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1,
) = bulk_microphysics_tendencies(
    RosenbrockAverage(DonorJacobian(), ImplicitGrowth(), NoLimiter()), cm, mp, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub,
)
