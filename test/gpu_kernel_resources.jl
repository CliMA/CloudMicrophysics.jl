#=
Static GPU kernel-resource regression tests.

Unlike `gpu_performance.jl`, this file measures nothing with a clock. It compiles
the hot microphysics kernels and reads the resources the compiler assigned:
register count and local-memory (spill) bytes. Those numbers are a deterministic
function of the source and the CUDA toolkit, so they can gate CI with hard
thresholds without being flaky under GPU contention.

Why registers: the per-thread register count sets a hard ceiling on how many
threads can be resident on an SM at once. `maxthreads` below is that ceiling as
computed by the driver, and it moves inversely with register usage. A change that
inflates register usage throttles occupancy long before it ever spills, so
registers are the earliest available warning.

Measured on an A100 (CUDA 12, Julia 1.11) comparing the 1-moment source before
and after the register-pressure work:

    1M Float64:  199 regs / 256 threads  ->  158 regs / 384 threads
    1M Float32:  115 regs / 512 threads  ->   71 regs / 896 threads
    2M (both FT): unchanged, as expected -- only the 1M path was touched

Note what did *not* move: spill bytes were byte-identical across that change
(240/744/136/384). Spill is therefore a regression guard here, not the signal
that separates the two implementations -- do not tighten it to zero.

Caveat: these standalone kernels are a *proxy*. The register-bound kernel seen in
full AMIP Nsight profiles is a much larger fused ClimaCore broadcast that inlines
far more than this isolated call, and it does hit the 255-register limit. This
test catches register-pressure regressions in the 1M path early and cheaply; it
does not reproduce the AMIP kernel's numbers.

Run with:
    CLIMACOMMS_DEVICE=CUDA julia --project=test test/gpu_kernel_resources.jl
=#

import Test as TT
using ClimaComms
ClimaComms.@import_required_backends

import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3

ClimaComms.device() isa ClimaComms.CUDADevice || error("No GPU found")
using CUDA
CUDA.allowscalar(false)

# The architectural per-thread register limit. Above this the compiler must spill.
const MAX_REGS = 255

"""
    kernel_resources(kernel, args...)

Compile `kernel(args...)` for the GPU *without launching it* and return the
static resource usage the compiler assigned.
"""
function kernel_resources(kernel, args...)
    k = CUDA.@cuda launch = false kernel(args...)
    mem = CUDA.memory(k)
    return (;
        registers = CUDA.registers(k),
        spill_bytes = mem.local,
        shared_bytes = mem.shared,
        maxthreads = CUDA.maxthreads(k),
    )
end

function report(name, res; max_registers, max_spill_bytes, min_maxthreads)
    @info """
    $name
      registers      = $(res.registers) / $MAX_REGS  (budget $max_registers)
      spill (local)  = $(res.spill_bytes) bytes      (budget $max_spill_bytes)
      max threads/blk= $(res.maxthreads)  (floor $min_maxthreads)
    """
    TT.@test res.registers <= max_registers
    TT.@test res.spill_bytes <= max_spill_bytes
    TT.@test res.maxthreads >= min_maxthreads
end

# --- Kernels -----------------------------------------------------------------
# Plain CUDA.jl kernels, one thread per grid point, mirroring how these functions
# are reached from a ClimaCore `@.` broadcast in ClimaAtmos.

function bulk_1m_kernel!(out, mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    @inbounds if i <= length(out)
        out[i] = BMT.bulk_microphysics_tendencies(
            BMT.Instantaneous(), BMT.Microphysics1Moment(), mp, tps,
            ρ[i], T[i], q_tot[i], q_lcl[i], q_icl[i], q_rai[i], q_sno[i],
        )
    end
    return nothing
end

function bulk_2m_kernel!(out, mp, tps, ρ, T, q_tot, q_lcl, q_rai, n_lcl, n_rai)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    @inbounds if i <= length(out)
        out[i] = BMT.bulk_microphysics_tendencies(
            BMT.Microphysics2Moment(), mp, tps,
            ρ[i], T[i], q_tot[i], q_lcl[i], n_lcl[i], q_rai[i], n_rai[i],
        )
    end
    return nothing
end

# --- Tests -------------------------------------------------------------------

# Budgets are set from measured values with headroom, tight enough that a revert
# to the pre-optimization 1-moment source (199 regs F64 / 115 regs F32) fails.
# Spill budgets sit just above the currently observed values -- they guard against
# new spilling, they do not encode a target of zero.
const BUDGETS = Dict(
    (:m1, Float64) => (max_registers = 175, max_spill_bytes = 256, min_maxthreads = 320),
    (:m1, Float32) => (max_registers = 90, max_spill_bytes = 160, min_maxthreads = 768),
    (:m2, Float64) => (max_registers = 135, max_spill_bytes = 768, min_maxthreads = 448),
    (:m2, Float32) => (max_registers = 90, max_spill_bytes = 416, min_maxthreads = 640),
)

function run_kernel_resource_tests(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    N = 1024
    dev(x) = CUDA.fill(FT(x), N)

    ρ, T = dev(1.0), dev(280.0)
    q_tot, q_lcl, q_icl = dev(1e-2), dev(1e-4), dev(1e-5)
    q_rai, q_sno = dev(1e-5), dev(1e-6)
    n_lcl, n_rai = dev(1e8), dev(1e4)

    TT.@testset "1-moment bulk tendencies kernel resources ($FT)" begin
        mp = CMP.Microphysics1MParams(FT)
        DT = @NamedTuple{dq_lcl_dt::FT, dq_icl_dt::FT, dq_rai_dt::FT, dq_sno_dt::FT}
        # `CUDA.zeros` would need `zero(::Type{<:NamedTuple})`, which is not
        # defined; the contents are irrelevant since we never launch.
        out = CUDA.CuArray{DT}(undef, N)
        res = kernel_resources(
            bulk_1m_kernel!, out, mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
        )
        report("Microphysics1Moment / Instantaneous ($FT)", res; BUDGETS[(:m1, FT)]...)
    end

    TT.@testset "2-moment bulk tendencies kernel resources ($FT)" begin
        mp = CMP.Microphysics2MParams(FT)
        DT = @NamedTuple{
            dq_lcl_dt::FT, dn_lcl_dt::FT, dq_rai_dt::FT, dn_rai_dt::FT,
            dq_ice_dt::FT, dq_rim_dt::FT, db_rim_dt::FT, dn_lcl_activation_dt::FT,
        }
        # `CUDA.zeros` would need `zero(::Type{<:NamedTuple})`, which is not
        # defined; the contents are irrelevant since we never launch.
        out = CUDA.CuArray{DT}(undef, N)
        res = kernel_resources(
            bulk_2m_kernel!, out, mp, tps, ρ, T, q_tot, q_lcl, q_rai, n_lcl, n_rai,
        )
        report("Microphysics2Moment ($FT)", res; BUDGETS[(:m2, FT)]...)
    end
end

TT.@testset "GPU kernel resource regressions" begin
    for FT in (Float64, Float32)
        run_kernel_resource_tests(FT)
    end
end
nothing
