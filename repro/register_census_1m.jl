# Measure register usage / local-memory spills of the 1M bulk tendency kernels.
import CUDA
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.ThermodynamicsInterface as TDI
using Printf

const FT = parse(Bool, get(ENV, "F32", "false")) ? Float32 : Float64
const N = 1024

mp = CMP.Microphysics1MParams(FT)
tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

# ---- kernels ---------------------------------------------------------------
function k_inst!(out, mp, tps, ρ, T, qt, ql, qi, qr, qs)
    i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if i <= length(out)
        @inbounds begin
            r = BMT.bulk_microphysics_tendencies(BMT.Instantaneous(), BMT.Microphysics1Moment(),
                mp, tps, ρ[i], T[i], qt[i], ql[i], qi[i], qr[i], qs[i])
            out[i] = r.dq_lcl_dt + r.dq_icl_dt + r.dq_rai_dt + r.dq_sno_dt
        end
    end
    return nothing
end

# runtime nsub (Int argument) -> real loop
function k_lin!(out, mp, tps, ρ, T, qt, ql, qi, qr, qs, Δt, nsub)
    i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if i <= length(out)
        @inbounds begin
            r = BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(), BMT.Microphysics1Moment(),
                mp, tps, ρ[i], T[i], qt[i], ql[i], qi[i], qr[i], qs[i], Δt, nsub)
            out[i] = r.dq_lcl_dt + r.dq_icl_dt + r.dq_rai_dt + r.dq_sno_dt
        end
    end
    return nothing
end

# compile-time nsub = 1 (loop fully unrolled / removed)
function k_lin1!(out, mp, tps, ρ, T, qt, ql, qi, qr, qs, Δt)
    i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if i <= length(out)
        @inbounds begin
            r = BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(), BMT.Microphysics1Moment(),
                mp, tps, ρ[i], T[i], qt[i], ql[i], qi[i], qr[i], qs[i], Δt, 1)
            out[i] = r.dq_lcl_dt + r.dq_icl_dt + r.dq_rai_dt + r.dq_sno_dt
        end
    end
    return nothing
end

# one substep only, no outer loop (the body of the loop)
function k_step!(out, mp, tps, ρ, T, qt, ql, qi, qr, qs, Δt)
    i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if i <= length(out)
        @inbounds begin
            r = BMT._linearized_implicit_step(BMT.Microphysics1Moment(),
                mp, tps, ρ[i], T[i], qt[i], ql[i], qi[i], qr[i], qs[i], Δt)
            out[i] = r.dq_lcl_dt + r.dq_icl_dt + r.dq_rai_dt + r.dq_sno_dt
        end
    end
    return nothing
end

# just the source terms (18 outputs summed)
function k_src!(out, mp, tps, ρ, T, qt, ql, qi, qr, qs)
    i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if i <= length(out)
        @inbounds begin
            r = BMT._microphysics_source_terms(BMT.Microphysics1Moment(),
                mp, tps, ρ[i], T[i], qt[i], ql[i], qi[i], qr[i], qs[i])
            out[i] = sum(values(r))
        end
    end
    return nothing
end

# ---- driver ----------------------------------------------------------------
out = CUDA.zeros(FT, N)
arr(v) = CUDA.fill(FT(v), N)
ρ, T, qt, ql, qi, qr, qs = arr(1.1), arr(268), arr(0.01), arr(1e-4), arr(1e-4), arr(1e-4), arr(1e-4)
Δt = FT(60)

function report(label, k, args...)
    kern = CUDA.@cuda launch=false k(args...)
    regs = CUDA.registers(kern)
    mem = CUDA.memory(kern)
    # time it too
    threads = 256; blocks = cld(N, threads)
    kern(args...; threads, blocks); CUDA.synchronize()
    t = CUDA.@elapsed kern(args...; threads, blocks)
    @printf("%-28s regs=%3d  local(spill)=%5d B  shared=%d  const=%d  t=%.1f µs\n",
        label, regs, mem.local, mem.shared, mem.constant, 1e6 * t)
    return kern
end

println("FT = $FT, device = ", CUDA.name(CUDA.device()))
report("Instantaneous", k_inst!, out, mp, tps, ρ, T, qt, ql, qi, qr, qs)
report("source_terms only", k_src!, out, mp, tps, ρ, T, qt, ql, qi, qr, qs)
report("one implicit step (no loop)", k_step!, out, mp, tps, ρ, T, qt, ql, qi, qr, qs, Δt)
report("LinearizedAverage nsub=1 lit", k_lin1!, out, mp, tps, ρ, T, qt, ql, qi, qr, qs, Δt)
report("LinearizedAverage nsub=1 rt", k_lin!, out, mp, tps, ρ, T, qt, ql, qi, qr, qs, Δt, 1)
report("LinearizedAverage nsub=3 rt", k_lin!, out, mp, tps, ρ, T, qt, ql, qi, qr, qs, Δt, 3)
