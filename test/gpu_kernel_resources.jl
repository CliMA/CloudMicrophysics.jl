#=
Static GPU kernel-resource regression tests.

Unlike `gpu_performance.jl`, this file measures nothing with a clock. It compiles
the hot microphysics kernels and reads the resources the compiler assigned:
register count, local-memory (spill) bytes, and the resulting occupancy ceiling.
Nothing is timed, so it cannot flake under GPU contention.

Why registers: the per-thread register count caps how many threads can be
resident on an SM. `maxthreads` is that ceiling as computed by the driver, and it
moves inversely with register usage. A change that inflates register usage
throttles occupancy long before it ever spills, so registers are the earliest
available warning.

## Baselines are per environment, not absolute

Register allocation is deterministic for a fixed (GPU architecture, CUDA toolkit,
Julia/LLVM) triple, but differs substantially across them. Measured evidence: the
2-moment kernel, whose source did not change, compiles to 118 registers on one
A100 machine and 175 on another. Absolute thresholds calibrated on one machine
therefore produce false failures everywhere else, so this test compares against
baselines recorded per environment key in

    gpu_kernel_resource_baselines.toml

An environment with no recorded baseline reports its numbers and skips gating
rather than failing. To calibrate a machine:

    CM_UPDATE_KERNEL_BASELINES=1 julia --project=test test/gpu_kernel_resources.jl

then commit the updated file. Updates merge, so calibrating one machine does not
discard another's numbers.

## What moves and what doesn't

Comparing the 1-moment source before and after the register-pressure work on one
A100 (CUDA 13, Julia 1.11.5):

    1M Float64:  199 regs / 256 threads  ->  158 regs / 384 threads
    1M Float32:  115 regs / 512 threads  ->   71 regs / 896 threads
    2M (both FT): unchanged, as expected -- only the 1M path was touched

Spill bytes were byte-identical across that change (240/744/136/384). Spill is a
regression guard here, not the signal separating the two implementations -- do
not tighten it to zero.

## CI hardware is not production hardware

Measured on the same (mod) source, same test, different machines:

    kernel        A100 / sm_80          P100 / sm_60 (CI)
    1M Float64    158 regs, 240 B, 384  255 regs, 288 B, 256
    2M Float64    118 regs, 744 B, 512  175 regs, 744 B, 256
    1M Float32     71 regs, 136 B, 896  118 regs, 136 B, 512
    2M Float32     74 regs, 384 B, 768   92 regs, 384 B, 640

Two things to take from this. Spill bytes are largely architecture-stable (three
of four identical); register counts are not, differing by up to 60%. And on the
P100 the 1M Float64 kernel sits at the 255-register ceiling even on the optimized
source -- so on that architecture this case may not be able to distinguish an
optimized 1-moment path from an unoptimized one, and its gate is close to
vacuous. Check the recorded P100 baseline before reading meaning into it.

The CliMA GPU CI agents are Pascal-era P100s while AMIP production runs on A100s,
so numbers recorded here do not predict production behavior. This test is a
*regression detector per environment*, not a cross-machine performance model.

Caveat: these standalone kernels are also a *proxy* for the AMIP kernel itself.
The register-bound kernel in full AMIP Nsight profiles is a much larger fused
ClimaCore broadcast that inlines far more than this isolated call.

Run with:
    CLIMACOMMS_DEVICE=CUDA julia --project=test test/gpu_kernel_resources.jl
=#

import Test as TT
import TOML
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

# Absolute register counts are NOT portable. Register allocation is deterministic
# for a fixed (GPU architecture, CUDA toolkit, Julia/LLVM) triple, but varies a
# lot across them -- the untouched 2-moment kernel compiles to 118 registers on
# one machine and 175 on another from identical source. So baselines are recorded
# per environment and compared against, rather than hardcoded as constants.
const BASELINE_PATH = joinpath(@__DIR__, "gpu_kernel_resource_baselines.toml")

# Regeneration: run this file with CM_UPDATE_KERNEL_BASELINES=1 on the machine in
# question and commit the resulting file.
const UPDATING = get(ENV, "CM_UPDATE_KERNEL_BASELINES", "false") in
                 ("1", "true", "yes")

"""
    env_key()

Fingerprint the parts of the toolchain that determine register allocation. Two
machines sharing this key should produce identical numbers; two that differ may
legitimately produce different ones.
"""
function env_key()
    cap = CUDA.capability(CUDA.device())
    return string(
        "sm_", cap.major, cap.minor,
        "-cuda", replace(string(CUDA.runtime_version()), "." => "_"),
        "-julia", replace(string(VERSION), "." => "_"),
    )
end

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

# Slack absorbs run-to-run jitter while still catching the tens-of-registers moves
# a real regression produces. The jitter is not hypothetical: the same source, on
# the same GPU, in back-to-back processes, has been observed compiling the 2M
# Float32 kernel to both 70 and 74 registers -- Julia's inference results are not
# bit-identical across processes, and register allocation follows inference.
const REG_SLACK = 8
const SPILL_SLACK = 64

load_baselines() =
    isfile(BASELINE_PATH) ? TOML.parsefile(BASELINE_PATH) : Dict{String, Any}()

function report(name, res, baselines, recorded)
    key = env_key()
    env_baselines = get(baselines, key, nothing)
    ref = env_baselines === nothing ? nothing : get(env_baselines, name, nothing)

    @info """
    $name  [$key]
      registers      = $(res.registers) / $MAX_REGS$(ref === nothing ? "" : "  (baseline $(ref["registers"]))")
      spill (local)  = $(res.spill_bytes) bytes$(ref === nothing ? "" : "  (baseline $(ref["spill_bytes"]))")
      max threads/blk= $(res.maxthreads)$(ref === nothing ? "" : "  (baseline $(ref["maxthreads"]))")
    """

    # Always stash the measurement so an update run writes a complete file.
    get!(recorded, key, Dict{String, Any}())[name] = Dict(
        "registers" => Int(res.registers),
        "spill_bytes" => Int(res.spill_bytes),
        "maxthreads" => Int(res.maxthreads),
    )

    UPDATING && return nothing

    if ref === nothing
        # An unknown toolchain is not a regression. Report and move on rather
        # than failing CI on a machine nobody has calibrated yet.
        @warn """
        No recorded baseline for this environment; skipping resource gating.
        To calibrate this machine:
            CM_UPDATE_KERNEL_BASELINES=1 julia --project=test test/gpu_kernel_resources.jl
        then commit $(basename(BASELINE_PATH)).
        """ env = key kernel = name
        return nothing
    end

    TT.@test res.registers <= ref["registers"] + REG_SLACK
    TT.@test res.spill_bytes <= ref["spill_bytes"] + SPILL_SLACK
    # `maxthreads` is deliberately reported but not gated. It is a derived
    # quantity -- the driver computes it from the register count -- so asserting
    # on it double-counts the register check, and because it moves in quantized
    # steps, a few registers of harmless jitter can drop it a whole step and fail
    # a run that the register budget correctly accepted.
    return nothing
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

function run_kernel_resource_tests(FT, baselines, recorded)
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
        report("Microphysics1Moment / Instantaneous ($FT)", res, baselines, recorded)
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
        report("Microphysics2Moment ($FT)", res, baselines, recorded)
    end
end

let
    baselines = load_baselines()
    recorded = Dict{String, Any}()

    TT.@testset "GPU kernel resource regressions" begin
        for FT in (Float64, Float32)
            run_kernel_resource_tests(FT, baselines, recorded)
        end
    end

    if UPDATING
        # Merge rather than overwrite: other machines' baselines must survive an
        # update run performed on this one.
        merged = merge(baselines, recorded)
        open(BASELINE_PATH, "w") do io
            TOML.print(io, merged)
        end
        @info "Wrote kernel resource baselines" path = BASELINE_PATH env = env_key()
    end
end
nothing
