import CUDA
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Utilities as UT
using Printf
const FT = Float64; const N = 1024
tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
mp0 = CMP.Microphysics1MParams(CP.create_toml_dict(FT))
out = CUDA.zeros(FT, N); arr(v) = CUDA.fill(FT(v), N)
ρ, T, qt, ql, qi, qr, qs = arr(1.1), arr(268), arr(0.01), arr(1e-4), arr(1e-4), arr(1e-4), arr(1e-4)
idx() = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
pr(label, k, args...) = (kern = CUDA.@cuda launch=false k(args...); @printf("%-40s regs=%3d local=%d B\n", label, CUDA.registers(kern), CUDA.memory(kern).local); flush(stdout))
function k_step!(out, mp, tps, ρ, T, qt, ql, qi, qr, qs, Δt)
    i = idx()
    if i <= length(out)
        @inbounds begin
            r = BMT._linearized_implicit_step(BMT.Microphysics1Moment(), mp, tps, ρ[i], T[i], qt[i], ql[i], qi[i], qr[i], qs[i], Δt)
            out[i] = r.dq_lcl_dt + r.dq_icl_dt + r.dq_rai_dt + r.dq_sno_dt
        end
    end
    return nothing
end
function k_lin!(out, mp, tps, ρ, T, qt, ql, qi, qr, qs, Δt, nsub)
    i = idx()
    if i <= length(out)
        @inbounds begin
            r = BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(), BMT.Microphysics1Moment(), mp, tps, ρ[i], T[i], qt[i], ql[i], qi[i], qr[i], qs[i], Δt, nsub)
            out[i] = r.dq_lcl_dt + r.dq_icl_dt + r.dq_rai_dt + r.dq_sno_dt
        end
    end
    return nothing
end
cpu_args_cold = (FT(1.1), FT(268), FT(0.01), FT(2e-4), FT(1e-4), FT(3e-4), FT(2e-4), FT(60))
cpu_args_warm = (FT(1.1), FT(276), FT(0.012), FT(2e-4), FT(1e-4), FT(3e-4), FT(2e-4), FT(60))
ref_cold = BMT._linearized_implicit_step(BMT.Microphysics1Moment(), mp0, tps, cpu_args_cold...)
ref_warm = BMT._linearized_implicit_step(BMT.Microphysics1Moment(), mp0, tps, cpu_args_warm...)
println("== baseline =="); pr("step", k_step!, out, mp0, tps, ρ, T, qt, ql, qi, qr, qs, FT(60)); pr("LinAvg rt nsub=3", k_lin!, out, mp0, tps, ρ, T, qt, ql, qi, qr, qs, FT(60), 3)

# ================= compute-once-and-pass prototype =================
@eval CM1 begin
    # λ⁻¹ with a caller-supplied n0 (avoids recomputing snow n0 = μ (ρ q)^ν)
    @inline function lambda_inverse_n0(mass::CMP.ParticleMass, q, ρ, n0)
        FT = UT.promote_typeof(q, ρ)
        (; r0, m0, me, Δm, χm, gamma_coeff) = mass
        qp = UT.clamp_to_nonneg(q); ρp = UT.clamp_to_nonneg(ρ)
        denom = χm * m0 * max(n0, UT.ϵ_numerics(n0)) * gamma_coeff
        λ_inv = (ρp * qp * r0^(me + Δm) / denom)^(1 / (me + Δm + 1))
        return max(r0 * FT(1e-5), λ_inv)
    end
    @inline function size_distr_parameters(mp, micro, thermo)
        (; q_rai, q_sno, q_icl) = micro
        ρ = thermo.ρ
        rai = mp.precip.rain; sno = mp.precip.snow; ice = mp.cloud.ice
        vr = mp.terminal_velocity.rain; vs = mp.terminal_velocity.snow
        n0_rai = get_n0(rai.pdf, q_rai, ρ)
        n0_sno = get_n0(sno.pdf, q_sno, ρ)
        n0_icl = get_n0(ice.pdf)
        λ_inv_rai = lambda_inverse_n0(rai.mass, q_rai, ρ, n0_rai)
        λ_inv_sno = lambda_inverse_n0(sno.mass, q_sno, ρ, n0_sno)
        λ_inv_icl = lambda_inverse_n0(ice.mass, q_icl, ρ, n0_icl)
        v0_rai = get_v0(vr, ρ); v0_sno = get_v0(vs, ρ)
        # per-species power laws shared by several processes (one pow each)
        accr_pow_rai = (rai.mass.r0 / λ_inv_rai)^(rai.area.ae + vr.ve + rai.area.Δa + vr.Δv)
        accr_pow_sno = (sno.mass.r0 / λ_inv_sno)^(sno.area.ae + vs.ve + sno.area.Δa + vs.Δv)
        vent_pow_rai = (rai.mass.r0 / λ_inv_rai)^((vr.ve + vr.Δv) / 2)
        vent_pow_sno = (sno.mass.r0 / λ_inv_sno)^((vs.ve + vs.Δv) / 2)
        vt_rai = terminal_velocity(rai, vr, ρ, q_rai, v0_rai, λ_inv_rai)
        vt_sno = terminal_velocity(sno, vs, ρ, q_sno, v0_sno, λ_inv_sno)
        return (; λ_inv_rai, n0_rai, v0_rai, λ_inv_sno, n0_sno, v0_sno, λ_inv_icl, n0_icl,
                  accr_pow_rai, accr_pow_sno, vent_pow_rai, vent_pow_sno, vt_rai, vt_sno)
    end
    @inline function _accretion_k(precip, vel, E, q_clo, q_pre, n0, v0, λ_inv, accr_pow)
        (; χv, gamma_accr) = vel
        (; a0, χa) = precip.area
        accr_rate = q_clo * E * n0 * a0 * v0 * χa * χv * λ_inv * gamma_accr / accr_pow
        ϵ = UT.ϵ_numerics(UT.promote_typeof(q_clo, q_pre))
        return ifelse(q_clo > ϵ && q_pre > ϵ, accr_rate, zero(accr_rate))
    end
    @inline accretion(::CMP.CloudLiquidRainAccretion, mp, tps, micro, thermo, sd = size_distr_parameters(mp, micro, thermo)) =
        _accretion_k(mp.precip.rain, mp.terminal_velocity.rain, mp.process_params.cloud_liquid_rain_accretion.e,
            micro.q_lcl, micro.q_rai, sd.n0_rai, sd.v0_rai, sd.λ_inv_rai, sd.accr_pow_rai)
    @inline function accretion(::CMP.CloudLiquidSnowAccretion, mp, tps, micro, thermo, sd = size_distr_parameters(mp, micro, thermo))
        S = _accretion_k(mp.precip.snow, mp.terminal_velocity.snow, mp.process_params.cloud_liquid_snow_accretion.e,
            micro.q_lcl, micro.q_sno, sd.n0_sno, sd.v0_sno, sd.λ_inv_sno, sd.accr_pow_sno)
        α = warm_accretion_melt_factor(tps, thermo.T)
        return (; S_accr = S, S_melt = α * S)
    end
    @inline accretion(::CMP.CloudIceRainAccretion, mp, tps, micro, thermo, sd = size_distr_parameters(mp, micro, thermo)) =
        _accretion_k(mp.precip.rain, mp.terminal_velocity.rain, mp.process_params.cloud_ice_rain_accretion.e,
            micro.q_icl, micro.q_rai, sd.n0_rai, sd.v0_rai, sd.λ_inv_rai, sd.accr_pow_rai)
    @inline accretion(::CMP.CloudIceSnowAccretion, mp, tps, micro, thermo, sd = size_distr_parameters(mp, micro, thermo)) =
        _accretion_k(mp.precip.snow, mp.terminal_velocity.snow, mp.process_params.cloud_ice_snow_accretion.e,
            micro.q_icl, micro.q_sno, sd.n0_sno, sd.v0_sno, sd.λ_inv_sno, sd.accr_pow_sno)

    # rain–snow: terminal velocities passed in, one pow per arm
    @inline function _accretion_snow_rain_k(type_j, E_ij, coeff_disp, q_i, q_j, ρ, n0_i, n0_j, λ_i_inv, λ_j_inv, v_ti, v_tj)
        (; r0, m0, me, Δm, χm, gamma_coeff) = type_j.mass
        δ = me + Δm
        Δv_eff = sqrt((v_ti - v_tj)^2 + coeff_disp * (v_ti^2 + v_tj^2))
        accr_rate = π / ρ * n0_i * n0_j * m0 * χm * E_ij * Δv_eff * gamma_coeff *
            (λ_j_inv / r0)^δ * λ_j_inv * λ_i_inv *
            (2 * λ_i_inv^2 + 2 * (δ + 1) * λ_i_inv * λ_j_inv + (δ + 2) * (δ + 1) * λ_j_inv^2)
        ϵ = UT.ϵ_numerics(UT.promote_typeof(q_i, q_j))
        return ifelse(q_i > ϵ && q_j > ϵ, accr_rate, zero(accr_rate))
    end
    @inline function accretion_snow_rain(::CMP.RainSnowAccretion, mp, tps, micro, thermo, sd = size_distr_parameters(mp, micro, thermo))
        (; q_rai, q_sno) = micro; (; ρ, T) = thermo
        (; e, coeff_disp) = mp.process_params.rain_snow_accretion
        S_rai_sno = _accretion_snow_rain_k(mp.precip.rain, e, coeff_disp, q_sno, q_rai, ρ, sd.n0_sno, sd.n0_rai, sd.λ_inv_sno, sd.λ_inv_rai, sd.vt_sno, sd.vt_rai)
        S_sno_rai = _accretion_snow_rain_k(mp.precip.snow, e, coeff_disp, q_rai, q_sno, ρ, sd.n0_rai, sd.n0_sno, sd.λ_inv_rai, sd.λ_inv_sno, sd.vt_rai, sd.vt_sno)
        α = warm_accretion_melt_factor(tps, T)
        return (; S_rai_sno, S_sno_rai, S_melt = α * S_rai_sno)
    end

    # ventilation factor shared by rain evap / snow subl / snow melt
    @inline function _vent_factor(vent, vel, aps, v0, λ_inv, vent_pow)
        (; ν_air, D_vapor) = aps
        (; χv, gamma_vent) = vel
        Sc = ν_air / max(D_vapor, UT.ϵ_numerics(D_vapor))
        return vent.a + vent.b * cbrt(Sc) / vent_pow * sqrt(2 * v0 * χv / ν_air * λ_inv) * gamma_vent
    end
    @inline function conv_q_rai_to_q_vap(::CMP.RainEvaporation, mp, tps, micro, thermo, sd = size_distr_parameters(mp, micro, thermo))
        (; q_tot, q_lcl, q_icl, q_rai, q_sno) = micro; (; ρ, T) = thermo
        aps = mp.air_properties; FT = eltype(ρ)
        S = TDI.supersaturation_over_liquid(tps, q_tot, q_lcl + q_rai, q_icl + q_sno, ρ, T)
        G = CO.G_func_liquid(aps, tps, T)
        F = _vent_factor(mp.precip.rain.vent, mp.terminal_velocity.rain, aps, sd.v0_rai, sd.λ_inv_rai, sd.vent_pow_rai)
        evap_rate = 4 * FT(π) * sd.n0_rai / ρ * S * G * sd.λ_inv_rai^2 * F
        cond = q_rai > UT.ϵ_numerics(q_rai) && S < 0
        return min(0, ifelse(cond, evap_rate, zero(evap_rate)))
    end
    @inline function _snow_subl_dep_rate(mp, tps, micro, thermo, sd = size_distr_parameters(mp, micro, thermo))
        (; q_tot, q_lcl, q_icl, q_rai, q_sno) = micro; (; ρ, T) = thermo
        aps = mp.air_properties; FT = eltype(ρ)
        S = TDI.supersaturation_over_ice(tps, q_tot, q_lcl + q_rai, q_icl + q_sno, ρ, T)
        G = CO.G_func_ice(aps, tps, T)
        F = _vent_factor(mp.precip.snow.vent, mp.terminal_velocity.snow, aps, sd.v0_sno, sd.λ_inv_sno, sd.vent_pow_sno)
        subl_rate = 4 * FT(π) * sd.n0_sno / ρ * S * G * sd.λ_inv_sno^2 * F
        return ifelse(q_sno > UT.ϵ_numerics(q_sno), subl_rate, zero(subl_rate))
    end
    @inline function conv_q_sno_to_q_rai(::CMP.SnowMelt, mp, tps, micro, thermo, sd = size_distr_parameters(mp, micro, thermo))
        q_sno = micro.q_sno; (; ρ, T) = thermo
        aps = mp.air_properties; FT = eltype(ρ)
        T_freeze = TDI.T_freeze(tps)
        L = TDI.Lf(tps, T)
        F = _vent_factor(mp.precip.snow.vent, mp.terminal_velocity.snow, aps, sd.v0_sno, sd.λ_inv_sno, sd.vent_pow_sno)
        snow_melt_rate = 4 * FT(π) * sd.n0_sno / ρ * aps.K_therm / L * (T - T_freeze) * sd.λ_inv_sno^2 * F
        cond = q_sno > UT.ϵ_numerics(q_sno) && T > T_freeze
        return ifelse(cond, snow_melt_rate, zero(snow_melt_rate))
    end
end

println("== compute-once-and-pass prototype (28 -> 16 pow) ==")
pr("step", k_step!, out, mp0, tps, ρ, T, qt, ql, qi, qr, qs, FT(60))
pr("LinAvg rt nsub=3", k_lin!, out, mp0, tps, ρ, T, qt, ql, qi, qr, qs, FT(60), 3)
new_cold = BMT._linearized_implicit_step(BMT.Microphysics1Moment(), mp0, tps, cpu_args_cold...)
new_warm = BMT._linearized_implicit_step(BMT.Microphysics1Moment(), mp0, tps, cpu_args_warm...)
rd(a, b) = maximum(abs.(values(a) .- values(b)) ./ max.(abs.(values(a)), 1e-30))
println("max rel diff vs baseline: cold=", rd(ref_cold, new_cold), " warm=", rd(ref_warm, new_warm))
# timing
threads = 256; blocks = cld(N, threads)
kern = CUDA.@cuda launch=false k_lin!(out, mp0, tps, ρ, T, qt, ql, qi, qr, qs, FT(60), 3)
kern(out, mp0, tps, ρ, T, qt, ql, qi, qr, qs, FT(60), 3; threads, blocks); CUDA.synchronize()
t = minimum(CUDA.@elapsed(kern(out, mp0, tps, ρ, T, qt, ql, qi, qr, qs, FT(60), 3; threads, blocks)) for _ in 1:20)
@printf("kernel time (N=%d, nsub=3): %.1f µs   [baseline earlier run: 182.3 µs]\n", N, 1e6 * t)
