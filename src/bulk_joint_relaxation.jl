#####
##### Joint Γ-consistent relaxation of the vapor-driven phase changes
#####

"""
    _relaxation_coefficient(S, τ, Γ, δ, tol)

Coefficient `c` of a vapor ↔ cloud-condensate phase change written as `rate = c δ`, with `δ`
the vapor excess over the saturation of that phase: `c = 1/(τ Γ)`, the non-equilibrium
scheme's rate per unit excess when its pool is not limiting (`_joint_vapor_transfers` lowers it
to the scheme's pool-bounded rate for a phase that stays a sink over the substep), or `0` when the process is disabled (`τ = Inf`) or its rate was zeroed by a switch
or limiter (`S = 0` at `|δ| > tol`).
"""
@inline function _relaxation_coefficient(S, τ, Γ, δ, tol)
    FT = typeof(S)
    τ_c = clamp(τ, eps(FT), floatmax(FT))
    c = ifelse(τ < floatmax(FT), one(FT) / (τ_c * Γ), zero(FT))
    off = (S == zero(FT)) & (abs(δ) > tol)
    return ifelse(off, zero(FT), c)
end

"""
    _ratio_coefficient(S, δ, tol)

Coefficient `c = S/δ` of a precipitation vapor exchange (rain evaporation, snow deposition
or sublimation) whose 1-moment rate is linear in the vapor excess `δ`; `0` when the excess
vanishes (`|δ| ≤ tol`) or the rate was zeroed by a switch. The 1-moment rain rate is a sink
only (`conv_q_rai_to_q_vap` returns `min(0, ⋅)`); a positive rain rate would be accepted.
"""
@inline function _ratio_coefficient(S, δ, tol)
    FT = typeof(S)
    excess = abs(δ) > tol
    return ifelse(excess, max(S / ifelse(excess, δ, one(FT)), zero(FT)), zero(FT))
end

"""
    _averaged_excess(δ₀, A, invτ, Δt)

Time average over `Δt` of the exact solution of `dδ/dt = A - δ/τ` (Morrison &
Milbrandt 2015, eq. C5): `δ̄ = A τ + (δ₀ - A τ) (τ/Δt) (1 - e^{-Δt/τ})`. Reduces to
`δ₀` when `invτ = 1/τ = 0` (no active process).
"""
@inline function _averaged_excess(δ₀, A, invτ, Δt)
    FT = typeof(δ₀)
    x = invτ * Δt                       # Δt/τ
    active = x > eps(FT)                # for x ≤ eps the average is δ₀ to round-off
    x_c = max(x, eps(FT))
    τ = Δt / x_c
    φ = -expm1(-x_c) / x_c
    return ifelse(active, A * τ + (δ₀ - A * τ) * φ, δ₀)
end

"""
    _joint_averaged_excesses(c_l, c_r, c_i, c_s, Γₗ, Γᵢ, κ_il, κ_li, δₗ, δᵢ, Δs, Δt)

Substep-averaged vapor excesses over liquid and over ice saturation of the joint
relaxation (see `_joint_vapor_transfers`), written in the excess of the primary phase.
"""
@inline function _joint_averaged_excesses(c_l, c_r, c_i, c_s, Γₗ, Γᵢ, κ_il, κ_li, δₗ, δᵢ, Δs, Δt)
    FT = typeof(δₗ)
    a = c_l + c_r
    b = c_i + c_s
    liquid_primary = Γₗ * a >= Γᵢ * b
    Γp = ifelse(liquid_primary, Γₗ, Γᵢ)
    κ = ifelse(liquid_primary, κ_il, κ_li)
    Cp = ifelse(liquid_primary, a, b)
    Cs = ifelse(liquid_primary, b, a)
    δp = ifelse(liquid_primary, δₗ, δᵢ)
    σ = ifelse(liquid_primary, one(FT), -one(FT))
    invτ = Γp * Cp + κ * Cs
    A = -σ * Δs * κ * Cs
    δ̄p = _averaged_excess(δp, A, invτ, Δt)
    δ̄s = δ̄p + σ * Δs
    return (ifelse(liquid_primary, δ̄p, δ̄s), ifelse(liquid_primary, δ̄s, δ̄p))
end

"""
    _saturation_state(tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno)

Vapor content (with the same expression as the rate functions, so that a rate of exactly zero
corresponds to a zero excess), saturations over liquid and ice, the latent-heat factors `Γ`
(own phase) and `κ` (other phase acting on this saturation), and the quantities of the vapor
check of `_linearized_implicit_step`: the lower saturation `q_smin`, its temperature derivative
`λ_min`, the corresponding `Γ_min`, and the round-off tolerance `q_tol` (one ulp of the larger saturation content) for a vanishing excess.
`Γ_min` uses the reference latent heat and the dry-air heat capacity of the substep temperature
update (`_latent_heating`), as the vapor limiter needs; `Γ` and `κ` use `L(T)` and the moist heat
capacity of the rate functions.
"""
@inline function _saturation_state(tps::TDI.PS, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno)
    FT = typeof(q_tot)
    qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_icl + q_sno)
    q_sl = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
    q_si = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
    Lᵥ = TDI.Lᵥ(tps, T)
    Lₛ = TDI.Lₛ(tps, T)
    Rᵥ = TDI.Rᵥ(tps)
    cₚ = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)
    dqsl_dT = CMNonEq.dqcld_dT(q_sl, Lᵥ, Rᵥ, T)
    dqsi_dT = CMNonEq.dqcld_dT(q_si, Lₛ, Rᵥ, T)
    Γₗ = CMNonEq.gamma_helper(Lᵥ, cₚ, dqsl_dT)
    Γᵢ = CMNonEq.gamma_helper(Lₛ, cₚ, dqsi_dT)
    κ_il = CMNonEq.gamma_helper(Lₛ, cₚ, dqsl_dT)   # ice-phase latent heat acting on the liquid saturation
    κ_li = CMNonEq.gamma_helper(Lᵥ, cₚ, dqsi_dT)   # liquid-phase latent heat acting on the ice saturation
    # the vapor limiter's factor uses the reference latent heat and the dry-air heat capacity of the
    # substep temperature update (`_latent_heating`), so that the scaled sources leave the vapor on the floor
    cp_d = TDI.TD.Parameters.cp_d(tps)
    Γₗ_min = CMNonEq.gamma_helper(TDI.TD.Parameters.LH_v0(tps), cp_d, dqsl_dT)
    Γᵢ_min = CMNonEq.gamma_helper(TDI.TD.Parameters.LH_s0(tps), cp_d, dqsi_dT)
    (q_smin, λ_min, Γ_min) = ifelse(q_si <= q_sl, (q_si, dqsi_dT, Γᵢ_min), (q_sl, dqsl_dT, Γₗ_min))
    # one ulp of the saturation content: excesses and source sums below it cannot be told from zero
    q_tol = eps(FT) * max(q_sl, q_si)
    return (; qᵥ, q_sl, q_si, Γₗ, Γᵢ, κ_il, κ_li, q_smin, λ_min, Γ_min, q_tol)
end

"""
    _joint_vapor_transfers(src, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt)

Transfers over the substep of the four vapor-driven phase changes (cloud liquid
condensation/evaporation, rain evaporation, cloud ice deposition/sublimation, snow
deposition/sublimation), treated as ONE relaxation of
the shared vapor excess with the latent-heat feedback of all of them (Morrison &
Milbrandt 2015, Appendix C, eqs. C1-C7, extended to snow and rain).

Each process has an instantaneous rate `c_k δ_k`, where `δ_k` is the excess over the
saturation of its phase and `c_k` its rate per unit excess (`_relaxation_coefficient`,
`_ratio_coefficient`). Written in the excess `δ_p` of the primary phase (the phase
with the larger excess-decay rate, so that the single-process limits are exact), with
`δ_s = δ_p + σ (q*_liq - q*_ice)` for the other phase,

    dδ_p/dt = A - δ_p/τ,   1/τ = Γ_p C_p + κ_sp C_s,   A = -σ (q*_liq - q*_ice) κ_sp C_s,

where `C_p`, `C_s` sum the coefficients of the processes of each phase,
`Γ_p = 1 + (L_p/c_p) dq*_p/dT`, and `κ_sp = 1 + (L_s/c_p) dq*_p/dT` is the effect of
the secondary phase's latent heat on the primary saturation. `A` is the
Wegener-Bergeron-Findeisen drive (liquid evaporating while ice deposits). `δ_p` is
time-averaged exactly over the substep (`_averaged_excess`) and the transfers are
`Δq_k = c_k δ̄_k Δt`, with sinks clamped to their pools. For a phase that is a sink at the
start of the substep and on average, the coefficient is lowered to the scheme's pool-bounded
rate per unit excess, `min(c_k, S_k/δ_k)`, and the averages are evaluated again, so a thin
cloud in dry air evaporates at the scheme's rate `S_k` for `Δt ≪ τ` instead of at `|δ_k|/(τ Γ)`. The saturation difference and
the `Γ`, `κ` factors are held at their start-of-substep values (as in MM15); as in P3,
pools are not tracked within the substep: the implicit solve shares each pool among its
sinks and a Γ-consistent vapor check on the solved state (`_linearized_implicit_step`)
removes any deposition that an exhausted pool could not feed.

Returns the transfers `(; Δq_lcl, Δq_rai, Δq_icl, Δq_sno)` [kg/kg over the substep,
positive = vapor → condensate] and, for the vapor check, the lower saturation
`q_smin = min(q*_liq, q*_ice)`, its temperature derivative `λ_min`, the corresponding
`Γ_min`, and the round-off tolerance `q_tol` used for the excesses.
"""
@inline function _joint_vapor_transfers(src, tps::TDI.PS, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt)
    (; qᵥ, q_sl, q_si, Γₗ, Γᵢ, κ_il, κ_li, q_smin, λ_min, Γ_min, q_tol) =
        _saturation_state(tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno)
    δₗ = qᵥ - q_sl
    δᵢ = qᵥ - q_si
    Δs = q_sl - q_si
    S_l = src.S_phase_change_vap_lcl.S
    S_i = src.S_phase_change_vap_icl.S
    S_r = src.S_phase_change_vap_rai.S
    S_s = src.S_phase_change_vap_sno.S
    c_l = _relaxation_coefficient(S_l, src.S_phase_change_vap_lcl.τ, Γₗ, δₗ, q_tol)
    c_i = _relaxation_coefficient(S_i, src.S_phase_change_vap_icl.τ, Γᵢ, δᵢ, q_tol)
    c_r = _ratio_coefficient(S_r, δₗ, q_tol)
    c_s = _ratio_coefficient(S_s, δᵢ, q_tol)
    (δ̄ₗ, δ̄ᵢ) = _joint_averaged_excesses(c_l, c_r, c_i, c_s, Γₗ, Γᵢ, κ_il, κ_li, δₗ, δᵢ, Δs, Δt)
    # pool-bounded evaporation and sublimation: when the scheme bounds the rate by the pool (S = -q/τ
    # instead of δ/(τΓ)), a phase that is a sink at the start of the substep and on average keeps the
    # scheme's rate per unit excess, min(c, S/δ), and the averages are evaluated again; a phase whose
    # averaged excess turns positive (pool exhausted, deposition starting) is not bounded
    FT = typeof(δₗ)
    sink_l = (δₗ < -q_tol) & (δ̄ₗ < zero(FT))
    sink_i = (δᵢ < -q_tol) & (δ̄ᵢ < zero(FT))
    c_l = ifelse(sink_l, min(c_l, max(S_l / min(δₗ, -q_tol), zero(FT))), c_l)
    c_i = ifelse(sink_i, min(c_i, max(S_i / min(δᵢ, -q_tol), zero(FT))), c_i)
    (δ̄ₗ, δ̄ᵢ) = _joint_averaged_excesses(c_l, c_r, c_i, c_s, Γₗ, Γᵢ, κ_il, κ_li, δₗ, δᵢ, Δs, Δt)
    Δq_lcl = max(c_l * δ̄ₗ * Δt, -q_lcl)
    Δq_rai = max(c_r * δ̄ₗ * Δt, -q_rai)
    Δq_icl = max(c_i * δ̄ᵢ * Δt, -q_icl)
    Δq_sno = max(c_s * δ̄ᵢ * Δt, -q_sno)
    return (; Δq_lcl, Δq_rai, Δq_icl, Δq_sno, q_smin, λ_min, Γ_min, q_tol)
end
