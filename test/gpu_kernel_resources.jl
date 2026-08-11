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
`TARGET_ARCH` (sm_80, A100, etc.) via
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

## Spilling only shows up under a register budget

Unconstrained, none of these kernels spill -- so `spill == 0` asserts something
that holds trivially and cannot catch a regression. The kernel this gate exists
for behaves quite differently in production: reached through a ClimaCore
broadcast in ClimaAtmos it sits at the 255-register cap and spilled ~25% of its
memory traffic before being fixed. A standalone one-thread-per-point kernel does
not reproduce that pressure.

So each kernel is also assembled repeatedly under `ptxas -maxrregcount=N`, walking
N down until it spills. That onset is a direct measure of how much register
pressure the function itself wants -- the property that pushed the real kernel
over -- and it moves long before the unconstrained register count does. It is a
proxy for the broadcast context, not a reproduction of it; reproducing it would
need ClimaCore in this package's test environment, and CloudMicrophysics sits
upstream of ClimaAtmos.

## Baselines

Baselines are indexed by CUDA toolkit version, because pinning the architecture
removes GPU-to-GPU variation but not toolkit variation: identical source gives
158 registers on 13.0 and 156 on 12.8 for the 1M Float64 kernel, and moves in
the other direction for 2M. When no set matches the running toolkit the test
warns and gates against the nearest one -- skipping is how this file once ran
zero assertions on CI while reporting success. Prefer adding a set over widening
REG_SLACK: slack spent on toolkit scatter is slack unavailable for detection.

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

# The architectural per-thread register limit; above it the compiler must spill.
const MAX_REGS = 255

# Baselines are indexed by CUDA toolkit version. Pinning the architecture removes
# GPU-to-GPU variation but not toolkit variation, and the two sets below differ by
# up to 9 registers from identical source -- in BOTH directions, which is what
# rules out a code cause: the 2M kernels are untouched by any change here, so they
# act as a control.
#
#   13.0  clima's A100s, where the science runs happen
#   12.8  the central cluster's CI agents (P100 hosts; the target is still sm_80)
#
# `onset` is the highest register budget at which ptxas still spills the kernel:
# a direct measure of how much register pressure the function wants. Gating on it
# is what makes spilling a live signal -- unconstrained, none of these kernels
# spill at all, so `spill == 0` alone asserts nothing.
#
# It is recorded PER TOOLKIT, and must be: the untouched 2M kernels spill under
# 12.8 at budgets where 13.0 does not, purely because ptxas allocates differently
# (122 vs 118 registers unconstrained). A shared budget would report that as a
# regression, which is exactly the confusion the per-toolkit split exists to
# prevent. `onset = nothing` means not yet measured on that toolkit: the test
# reports the measured value and tells you to record it, and keeps gating
# everything else.
const TOOLKIT_BASELINES = Dict(
    # Measured on clima's A100s.
    v"13.0" => Dict(
        ("1M", Float64) => (registers = 158, stack_frame = 240, spill = 0, onset = 160),
        ("2M", Float64) => (registers = 118, stack_frame = 744, spill = 0, onset = 112),
        ("1M", Float32) => (registers = 71, stack_frame = 136, spill = 0, onset = 67),
        ("2M", Float32) => (registers = 74, stack_frame = 384, spill = 0, onset = 62),
    ),
    # Registers and stack from a CI run; onsets not yet measured on this toolkit.
    v"12.8" => Dict(
        ("1M", Float64) => (registers = 156, stack_frame = 240, spill = 0, onset = nothing),
        ("2M", Float64) => (registers = 122, stack_frame = 744, spill = 0, onset = nothing),
        ("1M", Float32) => (registers = 62, stack_frame = 136, spill = 0, onset = nothing),
        ("2M", Float32) => (registers = 80, stack_frame = 384, spill = 0, onset = nothing),
    ),
)

"""
    baseline_set()

Return `(version, baselines, exact)` for the running toolkit. An unknown minor
version falls back to the newest set sharing its major version, and an unknown
major to the newest set overall -- warning, but still gating. Skipping is how
this file once ran zero assertions on CI while reporting success.

Spill onsets ARE per toolkit; see the note above TOOLKIT_BASELINES for why
sharing them across toolkits produces false regressions.
"""
function baseline_set()
    v = CUDA.runtime_version()
    key = VersionNumber(v.major, v.minor)
    haskey(TOOLKIT_BASELINES, key) && return key, TOOLKIT_BASELINES[key], true
    same_major = sort([k for k in keys(TOOLKIT_BASELINES) if k.major == v.major])
    fallback = isempty(same_major) ? maximum(keys(TOOLKIT_BASELINES)) : last(same_major)
    return fallback, TOOLKIT_BASELINES[fallback], false
end

# Slack absorbs run-to-run jitter: the same source on the same GPU in
# back-to-back processes has produced both 70 and 74 registers for the 2M
# Float32 kernel, because Julia's inference is not bit-identical across
# processes and register allocation follows it. Still far tighter than the
# tens-of-registers move a real regression produces (199 vs 158 for 1M Float64).
const REG_SLACK = 8
const STACK_SLACK = 64
const SPILL_TOLERANCE = 0
# Onsets jitter a little with inference, like registers do.
const ONSET_SLACK = 8

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

# PTX ISA versions to fall back through, newest first, when the default choice is
# rejected. Only used on toolchain mismatch; see `target_resources`. Floor is 7.0:
# sm_80 was introduced in PTX ISA 7.0, and ptxas rejects `-arch=sm_80` outright
# for anything older (verified: 6.5 and 6.2 both fail).
const PTX_ISA_FALLBACKS = [v"8.5", v"8.3", v"8.2", v"8.0", v"7.8", v"7.5", v"7.0"]

"""
    find_ptxas()

Prefer the `ptxas` CUDA.jl itself uses. It must match the PTX ISA that
`code_ptx` emits, and a `ptxas` found on PATH need not: on a CI agent whose
system CUDA was older than CUDA.jl's runtime, `Sys.which("ptxas")` produced

    fatal : Unsupported .version 8.7; current version is '8.2'

for every kernel, which silently reduced this whole file to a no-op.
"""
function find_ptxas()
    for f in (() -> only(CUDA.ptxas().exec), () -> Sys.which("ptxas"))
        p = try
            f()
        catch
            nothing
        end
        p !== nothing && isfile(p) && return p
    end
    return nothing
end

grab(re, s) = (m = match(re, s); m === nothing ? -1 : parse(Int, m[1]))

"""
    emit_ptx(kernel, tt; ptx)

Emit PTX for `TARGET_CAP` and return the file path. Emission is separated from
assembly because finding a kernel's spill onset means assembling the SAME PTX at
many register budgets; re-emitting each time would dominate the runtime.
"""
function emit_ptx(kernel, tt; ptx = nothing)
    io = IOBuffer()
    kwargs = (; cap = TARGET_CAP, kernel = true, dump_module = true)
    # `dump_module = true` is load-bearing; see the header note.
    CUDA.code_ptx(io, kernel, tt; (ptx === nothing ? kwargs : (; kwargs..., ptx))...)
    path = tempname() * ".ptx"
    write(path, String(take!(io)))
    return path
end

"""
    assemble(ptxas, ptx_file; maxregs)

Run `ptxas`, returning `(ok, log)`. `maxregs` caps registers per thread via
`-maxrregcount`, which is what makes spilling measurable at all.
"""
function assemble(ptxas, ptx_file; maxregs = nothing)
    out, err = IOBuffer(), IOBuffer()
    budget = maxregs === nothing ? `` : `-maxrregcount=$maxregs`
    cmd = `$ptxas -arch=$TARGET_ARCH $budget -v -o $(tempname()).cubin $ptx_file`
    proc = run(pipeline(ignorestatus(cmd); stdout = out, stderr = err))
    return success(proc), String(take!(err)) * String(take!(out))
end

"""
    resolve_ptx(ptxas, kernel, tt)

Emit PTX that this `ptxas` accepts, returning `(path, log)` from the
unconstrained assembly, or `nothing`.

CUDA.jl stamps the PTX ISA from the *runtime*, which can outrun the `ptxas`
binary actually available: a CI agent with a 12.8 runtime reported
"Unsupported .version 8.7; current version is '8.2'". Naming the exact ISA ptxas
reported is not enough -- `code_ptx` rejects any ISA newer than *LLVM* supports
("Requested PTX ISA 8.2 is not supported by LLVM 16.0.6"), and neither ceiling is
queryable without reaching into CUDA.jl internals. So walk candidates downwards
and take the first that clears both.

This does not change the generated code. CUDA.jl always has LLVM emit at LLVM's
own highest ISA (`llvm_ptx = maximum(llvm_ptxs)`, independent of the request) and
only rewrites the `.version` directive afterwards, and `-arch` still pins the
architecture. Verified on an A100 with CUDA 13.0: the 1M Float64 kernel gives 158
registers / 240 B stack at the default ISA and identically at 7.8, 7.5 and 7.0.
"""
function resolve_ptx(ptxas, kernel, tt)
    path = emit_ptx(kernel, tt)
    ok, log = assemble(ptxas, path)
    ok && return path, log

    m = match(r"current version is '(\d+\.\d+)'", log)
    ceiling = m === nothing ? nothing : VersionNumber(m.captures[1])
    # NB: not named `isa` -- that is an infix operator, so `@info "..." isa`
    # followed by `break` parses as the single expression `@info "..." isa break`,
    # swallowing the loop exit into the log call.
    for isa_ver in PTX_ISA_FALLBACKS
        ceiling !== nothing && isa_ver > ceiling && continue
        local candidate
        try
            candidate = emit_ptx(kernel, tt; ptx = isa_ver)
        catch err
            # LLVM cannot target this ISA; try an older one.
            @debug "PTX ISA $isa_ver unavailable" err
            continue
        end
        ok, log = assemble(ptxas, candidate)
        if ok
            @info "Re-emitted at a lower PTX ISA (codegen unchanged)" isa_ver
            return candidate, log
        end
    end

    @warn "ptxas failed; cannot measure $TARGET_ARCH resources" log
    return nothing
end

parse_resources(log) = (
    registers = grab(r"Used (\d+) registers", log),
    stack_frame = max(grab(r"(\d+) bytes stack frame", log), 0),
    spill_stores = max(grab(r"(\d+) bytes spill stores", log), 0),
    spill_loads = max(grab(r"(\d+) bytes spill loads", log), 0),
)

"""
    spill_onset(ptxas, ptx_file, start)

The highest register budget at which `ptxas` still spills this kernel -- a direct
measure of how much register pressure the function wants, and the property that
makes the real kernel spill inside a ClimaCore broadcast.

Scans downward rather than bisecting: spill is not monotonic in the budget. The
2M Float64 kernel spills 4 B at 112, nothing at 110 and 108, then 4 B again at
106, so a bisection would land on whichever side it happened to probe.
"""
function spill_onset(ptxas, ptx_file, start)
    for b in start:-2:32
        ok, log = assemble(ptxas, ptx_file; maxregs = b)
        ok || continue
        max(grab(r"(\d+) bytes spill stores", log), 0) > 0 && return b
    end
    return nothing
end

function check(name, key, kernel, args...)
    baseline_version, baselines, exact = baseline_set()
    ref = baselines[key]

    ptxas = find_ptxas()
    if ptxas === nothing
        # Deliberately a failure, not a skip. An earlier version warned and
        # returned, so a broken toolchain produced "0 tests, all passed" -- a
        # gate that silently protects nothing is worse than one that complains.
        @error "No ptxas found; cannot measure $TARGET_ARCH resources for $name"
        TT.@test false
        return nothing
    end

    gargs = map(CUDA.cudaconvert, args)
    tt = Tuple{map(Core.Typeof, gargs)...}
    resolved = resolve_ptx(ptxas, kernel, tt)
    if resolved === nothing
        @error "Could not measure $TARGET_ARCH resources for $name; see ptxas log above"
        TT.@test false
        return nothing
    end
    ptx_file, log = resolved
    res = parse_resources(log)

    # Start the scan above the unconstrained allocation: the onset can sit ABOVE
    # it, because under -maxrregcount ptxas changes strategy and can spill at a
    # budget it would not otherwise need (1M Float64: 158 unconstrained, onset 160).
    onset = spill_onset(ptxas, ptx_file, res.registers + 40)

    @info """
    $name  [compiled for $TARGET_ARCH, baselines from CUDA $baseline_version]
      registers   = $(res.registers) / $MAX_REGS   (baseline $(ref.registers))
      stack frame = $(res.stack_frame) bytes       (baseline $(ref.stack_frame))
      spill       = $(res.spill_stores) st / $(res.spill_loads) ld bytes
      spill onset = $(something(onset, "none")) registers   (baseline $(something(ref.onset, "unrecorded")))
    """

    # An unmatched toolkit is noted but does NOT disable the gate. Skipping here
    # is how this file previously ran zero assertions on CI while reporting
    # success. If a toolkit really does allocate differently, that surfaces as a
    # failure carrying both numbers, which is diagnosable; silence is not.
    exact || @warn """
    No baselines recorded for CUDA $(CUDA.runtime_version()); gating against the
    CUDA $baseline_version set instead. If this fails on registers only and the
    delta is small, re-measure on this toolkit and add a set to
    TOOLKIT_BASELINES rather than widening REG_SLACK.
    """ kernel = name

    TT.@test res.registers <= ref.registers + REG_SLACK
    TT.@test res.stack_frame <= ref.stack_frame + STACK_SLACK
    # Unconstrained these never spill, so this is a floor, not a live assertion.
    TT.@test res.spill_stores <= ref.spill + SPILL_TOLERANCE

    # The live one: how hard the kernel has to be squeezed before it spills.
    if ref.onset === nothing
        # Gating an onset against another toolkit's number is precisely the bug
        # this replaced -- the untouched 2M kernels spilled at CUDA 13.0's budgets
        # under 12.8 purely because ptxas allocates differently. Report instead of
        # asserting against a baseline we do not have; everything else still gates.
        @warn """
        No spill-onset baseline for CUDA $(CUDA.runtime_version()). Measured
        $(something(onset, "none")) for this kernel -- record `onset = $(something(onset, "nothing"))`
        in its TOOLKIT_BASELINES entry to gate it.
        """ kernel = name
    elseif onset === nothing
        # Never spilled in the window: strictly better than baseline.
        @info "No spill found down to a 32-register budget" kernel = name
    else
        TT.@test onset <= ref.onset + ONSET_SLACK
    end
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
