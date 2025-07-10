import CloudMicrophysics.Common as CMO
import CloudMicrophysics.HetIceNucleation as CMI_het
import CloudMicrophysics.HomIceNucleation as CMI_hom
import CloudMicrophysics.MicrophysicsNonEq as MNE
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI
import Distributions as DS
import SpecialFunctions as SF

# helper function to limit the tendency for noneq
function limit(q, dt, n::Int)
    return q / dt / n
end

function aerosol_activation(::Empty, state)
    FT = eltype(state)
    return FT(0)
end

function aerosol_activation(params::AeroAct, state)
    (; aap, aerosol, aero_σ_g, r_nuc, const_dt) = params
    (; T, Sₗ, Nₐ) = state
    FT = eltype(state)

    ad = AM.Mode_κ(
        r_nuc,
        aero_σ_g,
        params.Nₐ,
        (FT(1.0),),
        (FT(1.0),),
        (aerosol.M,),
        (aerosol.κ,),
    )
    all_ad = AM.AerosolDistribution((ad,))

    smax = (Sₗ - 1) < 0 ? FT(0) : (Sₗ - 1)
    sm = AA.critical_supersaturation(aap, all_ad, T)
    u = 2 * log(sm[1] / smax) / 3 / sqrt(2) / log(ad.stdev)

    return max(FT(0), ad.N * (1 / 2) * (1 - SF.erf(u)) - (params.Nₐ - Nₐ)) / FT(1)
end

function deposition_nucleation(::Empty, state, dY)
    FT = eltype(state)
    return FT(0)
end

function deposition_nucleation(params::MohlerAF, state, dY)
    (; ips, aerosol, const_dt, tps) = params
    (; Sₗ, T, Nₐ, Nᵢ) = state
    Sᵢ = ξ(tps, T) * Sₗ
    FT = eltype(state)
    if Sᵢ >= ips.deposition.Sᵢ_max
        @warn("Supersaturation exceeds Sᵢ_max. No dust will be activated.")
    end

    AF =
        Sᵢ >= ips.deposition.Sᵢ_max ? FT(0) :
        AF = CMI_het.dust_activated_number_fraction(
            aerosol,
            ips.deposition,
            Sᵢ,
            T,
        )
    return max(FT(0), AF * Nₐ - Nᵢ) / const_dt
end

function deposition_nucleation(params::MohlerRate, state, dY)
    (; ips, aerosol, tps, const_dt) = params
    (; T, Nₐ, Sₗ) = state
    FT = eltype(state)
    Sᵢ = ξ(tps, T) * Sₗ
    dSᵢdt = ξ(tps, T) * dY[1]

    if Sᵢ >= ips.deposition.Sᵢ_max
        @warn("Supersaturation exceeds Sᵢ_max. No dust will be activated.")
        dNi_dt = FT(0)
    else
        dNi_dt = CMI_het.MohlerDepositionRate(
            aerosol,
            ips.deposition,
            Sᵢ,
            T,
            dSᵢdt,
            Nₐ,
        )
    end

    return min(max(dNi_dt, FT(0)), Nₐ / const_dt)
end

function deposition_nucleation(params::ABDINM, state, dY)
    FT = eltype(state)
    (; tps, aerosol, r_nuc, const_dt) = params
    (; T, p_air, qᵥ, qₗ, qᵢ, Nₐ) = state

    Rᵥ = TDI.Rᵥ(tps)
    R_air = TDI.Rₘ(tps, qᵥ + qₗ + qᵢ, qₗ, qᵢ)
    e = eᵥ(qᵥ, p_air, R_air, Rᵥ)

    Δa_w = CMO.a_w_eT(tps, e, T) - CMO.a_w_ice(tps, T)
    J = CMI_het.deposition_J(aerosol, Δa_w)
    A = 4 * FT(π) * r_nuc^2
    return min(J * Nₐ * A, Nₐ / const_dt)
end

function deposition_nucleation(params::P3_dep, state, dY)
    FT = eltype(state)
    (; ips, const_dt) = params
    (; T, Nᵢ, Nₐ) = state
    Nᵢ_depo = CMI_het.P3_deposition_N_i(ips.p3, T)
    return min(max(FT(0), (Nᵢ_depo - Nᵢ) / const_dt), Nₐ / const_dt)
end

function immersion_freezing(::Empty, PSD_liq, state)
    FT = eltype(state)
    return FT(0)
end

function immersion_freezing(params::ABIFM, PSD_liq, state)
    (; T, p_air, qᵥ, qₗ, qᵢ, Nₗ) = state
    (; tps, aerosol, A_aer, const_dt) = params
    FT = eltype(state)

    Rᵥ = TDI.Rᵥ(tps)
    R_air = TDI.Rₘ(tps, qᵥ + qₗ + qᵢ, qₗ, qᵢ)
    e = eᵥ(qᵥ, p_air, R_air, Rᵥ)

    Δa_w = CMO.a_w_eT(tps, e, T) - CMO.a_w_ice(tps, T)

    J = CMI_het.ABIFM_J(aerosol, Δa_w)
    return min(J * Nₗ * A_aer, (Nₗ / const_dt))
end

function immersion_freezing(params::P3_het, PSD_liq, state)
    FT = eltype(state)
    (; const_dt, ips) = params
    (; T, Nₗ, Nᵢ) = state
    Nᵢ_het = CMI_het.P3_het_N_i(ips.p3, T, Nₗ, PSD_liq.V, const_dt)
    return min(max(FT(0), (Nᵢ_het - Nᵢ) / const_dt), (Nₗ / const_dt))
end

function immersion_freezing(params::Frostenberg_random, PSD_liq, state)
    FT = eltype(state)
    (; ip, sampling_interval, const_dt) = params
    (; T, Nₗ, Nᵢ, t) = state
    if mod(t, sampling_interval) == 0
        μ = CMI_het.INP_concentration_mean(T)
        INPC = exp(rand(DS.Normal(μ, ip.σ)))
    else
        INPC = 0
    end
    return min(Nₗ, max(FT(0), INPC - Nᵢ)) / const_dt
end

function immersion_freezing(params::Frostenberg_mean, PSD_liq, state)
    FT = eltype(state)
    (; ip, const_dt) = params
    (; T, Nₗ, Nᵢ) = state
    INPC = exp(CMI_het.INP_concentration_mean(T))
    return min(Nₗ, max(FT(0), INPC - Nᵢ)) / const_dt
end

function INPC_model(params, state)
    FT = eltype(state)
    return FT(0)
end

function INPC_model(params::Frostenberg_stochastic, state)
    FT = eltype(state)
    (; ip, γ, const_dt) = params
    (; T, ln_INPC, t) = state

    μ = CMI_het.INP_concentration_mean(T)
    g = ip.σ * sqrt(2 * γ)

    dln_INPC = -γ * (ln_INPC - μ) * const_dt + g * sqrt(const_dt) * randn()

    return dln_INPC / const_dt
end

function immersion_freezing(params::Frostenberg_stochastic, PSD_liq, state)
    FT = eltype(state)
    (; ip, γ, const_dt) = params
    (; T, Nₗ, Nᵢ, ln_INPC, t) = state
    return min(Nₗ, max(FT(0), exp(ln_INPC) - Nᵢ)) / const_dt
end

function homogeneous_freezing(::Empty, PSD_liq, state)
    FT = eltype(state)
    return FT(0)
end

function homogeneous_freezing(params::ABHOM, PSD_liq, state)
    FT = eltype(state)
    (; tps, ips, const_dt) = params
    (; T, p_air, qᵥ, qₗ, qᵢ, Nₗ, Sₗ) = state

    e = TDI.saturation_vapor_pressure_over_liquid(tps, T) * Sₗ

    Δa_w = CMO.a_w_eT(tps, e, T) - CMO.a_w_ice(tps, T)
    J = CMI_hom.homogeneous_J_linear(ips.homogeneous, Δa_w)

    return min(J * Nₗ * PSD_liq.V, Nₗ / const_dt)
end

function homogeneous_freezing(params::P3_hom, PSD_liq, state)
    FT = eltype(state)
    (; Nₗ, T) = state
    (; const_dt) = params
    #TODO - use ClimaParams
    return T < FT(233.15) && Nₗ > FT(0) ? Nₗ / const_dt : FT(0)
end

function condensation(::Empty, PSD_liq, state, ρ_air)
    FT = eltype(state)
    return FT(0)
end

function condensation(params::CondParams, PSD_liq, state, ρ_air)
    FT = eltype(state)
    (; aps, tps, const_dt) = params
    (; Sₗ, T, Nₗ, qᵥ, qₗ) = state
    Gₗ = CMO.G_func_liquid(aps, tps, T)
    dqₗ_dt = 4 * FT(π) / ρ_air * (Sₗ - 1) * Gₗ * PSD_liq.r * Nₗ

    return ifelse(
        dqₗ_dt > FT(0),
        min(dqₗ_dt, limit(qᵥ, const_dt, 1)),
        -min(abs(dqₗ_dt), limit(qₗ, const_dt, 1)),
    )
end

function condensation(params::NonEqCondParams_simple, PSD, state, ρ_air)

    FT = eltype(state)
    (; Sₗ, qₗ, qᵥ) = state
    (; tps, liquid) = params

    q_sat_liq = max(Sₗ * qᵥ - qᵥ, 0)

    new_q = MNE.conv_q_vap_to_q_liq_ice(liquid, q_sat_liq, qₗ)

    return new_q
end

function condensation(params::NonEqCondParams, PSD, state, ρ_air)
    FT = eltype(state)
    (; T, qₗ, qᵥ, qᵢ) = state
    (; tps, liquid, dt) = params

    if qᵥ + qₗ > FT(0)

        qₜ = qᵥ + qₗ + qᵢ

        cond_rate = MNE.conv_q_vap_to_q_liq_ice_MM2015(liquid, tps, qₜ, qₗ, qᵢ, FT(0), FT(0), ρ_air, T)

        # Using same limiter as ClimaAtmos for now
        # Not sure why, but without intermediate storing of the tendencies for the
        # if/else branch this code segfaults on julia v1.11 (works fine on v1.10)
        cond_limit = min(cond_rate, limit(qᵥ, dt, 1))
        evap_limit = min(abs(cond_rate), limit(qₗ, dt, 1))
        ret = ifelse(cond_rate > FT(0), cond_limit, -1 * evap_limit)
        return ret
    else
        return FT(0)
    end
end

function deposition(::Empty, PSD_ice, state, ρ_air)
    FT = eltype(state)
    return FT(0)
end

function deposition(params::DepParams, PSD_ice, state, ρ_air)
    (; aps, tps, const_dt) = params
    (; T, Sₗ, Nᵢ, qᵥ, qᵢ) = state
    FT = eltype(state)
    Sᵢ = ξ(tps, T) * Sₗ
    Gᵢ = CMO.G_func_ice(aps, tps, T)
    dqᵢ_dt = 4 * FT(π) / ρ_air * (Sᵢ - 1) * Gᵢ * PSD_ice.r * Nᵢ

    return ifelse(
        dqᵢ_dt > FT(0),
        min(dqᵢ_dt, limit(qᵥ, const_dt, 1)),
        -min(abs(dqᵢ_dt), limit(qᵢ, const_dt, 1)),
    )
end

function deposition(params::NonEqDepParams_simple, PSD, state, ρ_air)
    FT = eltype(state)
    (; T, Sₗ, qᵥ, qᵢ) = state

    (; tps, ice) = params

    Sᵢ = S_i(tps, T, Sₗ)
    q_sat_ice = max(Sᵢ * qᵥ - qᵥ, 0)

    new_q = MNE.conv_q_vap_to_q_liq_ice(ice, q_sat_ice, qᵢ)

    return new_q
end

function deposition(params::NonEqDepParams, PSD, state, ρ_air)
    FT = eltype(state)
    (; T, qₗ, qᵥ, qᵢ) = state

    (; tps, ice, dt) = params

    if qᵥ + qᵢ > FT(0)
        qₜ = qᵥ + qₗ + qᵢ

        dep_rate = MNE.conv_q_vap_to_q_liq_ice_MM2015(ice, tps, qₜ, qₗ, qᵢ, FT(0), FT(0), ρ_air, T)

        # Using same limiter as ClimaAtmos for now
        # Not sure why, but without intermediate storing of the tendencies for the
        # if/else branch this code segfaults on julia v1.11 (works fine on v1.10)
        dep_limit = min(dep_rate, limit(qᵥ, dt, 1))
        sub_limit = min(abs(dep_rate), limit(qᵢ, dt, 1))
        ret = ifelse(dep_rate > FT(0), dep_limit, -1 * sub_limit)
        return ret
    else
        return FT(0)
    end
end
