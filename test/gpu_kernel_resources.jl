#=
Static GPU kernel-resource regression tests.

Unlike `gpu_performance.jl`, this file measures nothing with a clock. It compiles
the hot microphysics kernels and reads the resources the compiler assigned:
register count, stack frame, and register spilling. Nothing is timed, so it
cannot flake under GPU contention.

Why registers: the per-thread register count caps how many threads can be
resident on an SM, so a change that inflates register usage throttles occupancy
long before it ever spills. Registers are the earliest available warning.

## Kernels are compiled for a FIXED architecture, not the local GPU

Register allocation depends on the target architecture, so numbers measured on
one GPU do not carry to another: the untouched 2-moment kernel compiles to 118
registers for sm_80 and 175 for sm_60 from identical source. That makes results
from a CI agent's GPU meaningless as a gate on the hardware production runs on.

So this test does not measure the local GPU at all. It compiles each kernel for
`TARGET_ARCH` (sm_80, A100 -- what CliMA actually runs science on) via
`code_ptx(...; cap)` and reads the resource usage back from `ptxas -v`. Any
CUDA-capable GPU can host the test; the numbers describe an A100 regardless.
Verified on an A100 (CUDA 13.0): this path reproduces `CUDA.registers` exactly
for all four kernels below.

Two details make it work, both easy to get wrong:

  * `dump_module = true` is REQUIRED. Without it `code_ptx` emits only the entry
    function, ptxas cannot resolve the called Julia device functions, and the
    only way to assemble is `--compile-only` (relocatable). That forces ABI
    calling conventions at the unresolved call sites and inflates the 1M Float64
    kernel from 158 registers to 250 -- a plausible-looking number that is wrong.
  * Do NOT pass `--compile-only`. See above.

## Stack frame is not spilling

`CUDA.memory(k).local` reports *local memory*, which is dominated by the stack
frame and is not the same thing as register spilling. ptxas separates them, and
for these kernels the split is: nonzero stack frame (the returned NamedTuple),
zero spill. Do not read the stack frame numbers as spill.

## Baselines

The numbers below are for `TARGET_ARCH` and were measured with the CUDA version
recorded in `BASELINE_CUDA`. Pinning the architecture removes GPU-to-GPU
variation but not toolkit variation: a different ptxas can still allocate
differently, so on a mismatched CUDA version the test reports and skips rather
than failing on a difference it cannot attribute.

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

ClimaComms.device() isa ClimaComms.CUDADevice || error("No GPU found")
using CUDA
CUDA.allowscalar(false)

# The architecture we gate on: A100. Independent of the GPU running the test.
const TARGET_CAP = v"8.0"
const TARGET_ARCH = "sm_$(TARGET_CAP.major)$(TARGET_CAP.minor)"

# CUDA toolkit the baselines were measured with. A different major version may
# allocate differently; see the note above.
const BASELINE_CUDA = v"13.0"

# The architectural per-thread register limit; above it the compiler must spill.
const MAX_REGS = 255

# Measured for sm_80 with CUDA 13.0. Registers are the gate; stack frame is
# recorded to catch a blow-up in returned-value size; spill should stay zero.
const BASELINES = Dict(
    ("1M", Float64) => (registers = 158, stack_frame = 240, spill = 0),
    ("2M", Float64) => (registers = 118, stack_frame = 744, spill = 0),
    ("1M", Float32) => (registers = 71, stack_frame = 136, spill = 0),
    ("2M", Float32) => (registers = 74, stack_frame = 384, spill = 0),
)

# Slack absorbs run-to-run jitter: the same source on the same GPU in
# back-to-back processes has produced both 70 and 74 registers for the 2M
# Float32 kernel, because Julia's inference is not bit-identical across
# processes and register allocation follows it. Still far tighter than the
# tens-of-registers move a real regression produces (199 vs 158 for 1M Float64).
const REG_SLACK = 8
const STACK_SLACK = 64
const SPILL_TOLERANCE = 0

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

# --- Measurement -------------------------------------------------------------

find_ptxas() = Sys.which("ptxas")

grab(re, s) = (m = match(re, s); m === nothing ? -1 : parse(Int, m[1]))

"""
    target_resources(kernel, args...)

Compile `kernel(args...)` for `TARGET_ARCH` and return the resources ptxas
assigned. The kernel is never launched, and never loaded onto the local device --
which is what allows a foreign architecture to be targeted.
"""
function target_resources(kernel, args...)
    ptxas = find_ptxas()
    ptxas === nothing && return nothing

    gargs = map(CUDA.cudaconvert, args)
    tt = Tuple{map(Core.Typeof, gargs)...}

    io = IOBuffer()
    # `dump_module = true` is load-bearing; see the header note.
    CUDA.code_ptx(io, kernel, tt; cap = TARGET_CAP, kernel = true,
                  dump_module = true)
    ptx_file = tempname() * ".ptx"
    write(ptx_file, String(take!(io)))

    out = IOBuffer()
    err = IOBuffer()
    cmd = `$ptxas -arch=$TARGET_ARCH -v -o $(tempname()).cubin $ptx_file`
    proc = run(pipeline(ignorestatus(cmd); stdout = out, stderr = err))
    log = String(take!(err)) * String(take!(out))
    if !success(proc)
        @warn "ptxas failed; cannot measure $TARGET_ARCH resources" log
        return nothing
    end

    return (
        registers = grab(r"Used (\d+) registers", log),
        stack_frame = max(grab(r"(\d+) bytes stack frame", log), 0),
        spill_stores = max(grab(r"(\d+) bytes spill stores", log), 0),
        spill_loads = max(grab(r"(\d+) bytes spill loads", log), 0),
    )
end

toolkit_matches() = CUDA.runtime_version().major == BASELINE_CUDA.major

function check(name, key, kernel, args...)
    res = target_resources(kernel, args...)
    if res === nothing
        @warn "Skipping $name: could not compile for $TARGET_ARCH (ptxas missing?)"
        return nothing
    end
    ref = BASELINES[key]

    @info """
    $name  [compiled for $TARGET_ARCH]
      registers   = $(res.registers) / $MAX_REGS   (baseline $(ref.registers))
      stack frame = $(res.stack_frame) bytes       (baseline $(ref.stack_frame))
      spill       = $(res.spill_stores) st / $(res.spill_loads) ld bytes
    """

    if !toolkit_matches()
        @warn """
        CUDA $(CUDA.runtime_version()) differs from the CUDA $(BASELINE_CUDA) the
        baselines were measured with; reporting only. Re-measure and update
        BASELINES (and BASELINE_CUDA) if this toolkit is the new normal.
        """ kernel = name
        return nothing
    end

    TT.@test res.registers <= ref.registers + REG_SLACK
    TT.@test res.stack_frame <= ref.stack_frame + STACK_SLACK
    # These kernels do not spill at sm_80; any spilling is a real regression.
    TT.@test res.spill_stores <= ref.spill + SPILL_TOLERANCE
    return nothing
end

# --- Tests -------------------------------------------------------------------

function run_kernel_resource_tests(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    N = 1024
    dev(x) = CUDA.fill(FT(x), N)

    ρ, T = dev(1.0), dev(280.0)
    q_tot, q_lcl, q_icl = dev(1e-2), dev(1e-4), dev(1e-5)
    q_rai, q_sno = dev(1e-5), dev(1e-6)
    n_lcl, n_rai = dev(1e8), dev(1e4)

    TT.@testset "1-moment bulk tendencies kernel resources ($FT)" begin
        DT = @NamedTuple{dq_lcl_dt::FT, dq_icl_dt::FT, dq_rai_dt::FT, dq_sno_dt::FT}
        # `CUDA.zeros` would need `zero(::Type{<:NamedTuple})`, which is not
        # defined; contents are irrelevant since the kernel is never launched.
        out = CUDA.CuArray{DT}(undef, N)
        check("Microphysics1Moment / Instantaneous ($FT)", ("1M", FT),
              bulk_1m_kernel!, out, CMP.Microphysics1MParams(FT), tps,
              ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno)
    end

    TT.@testset "2-moment bulk tendencies kernel resources ($FT)" begin
        DT = @NamedTuple{
            dq_lcl_dt::FT, dn_lcl_dt::FT, dq_rai_dt::FT, dn_rai_dt::FT,
            dq_ice_dt::FT, dq_rim_dt::FT, db_rim_dt::FT, dn_lcl_activation_dt::FT,
        }
        out = CUDA.CuArray{DT}(undef, N)
        check("Microphysics2Moment ($FT)", ("2M", FT),
              bulk_2m_kernel!, out, CMP.Microphysics2MParams(FT), tps,
              ρ, T, q_tot, q_lcl, q_rai, n_lcl, n_rai)
    end
end

@info "GPU kernel resource test" host_gpu = CUDA.name(CUDA.device()) target =
    TARGET_ARCH cuda = CUDA.runtime_version()

TT.@testset "GPU kernel resource regressions ($TARGET_ARCH)" begin
    for FT in (Float64, Float32)
        run_kernel_resource_tests(FT)
    end
end
nothing
