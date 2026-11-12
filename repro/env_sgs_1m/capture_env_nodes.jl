# Usage: julia --project=.buildkite capture_env_nodes.jl <out dir>
#
# Spins up the `set_CM_cache.jl` state on the GPU, profiles
# `set_microphysics_tendency_cache!` for reference, then copies the inputs of
# the environment SGS-quadrature 1M broadcast to the CPU and re-evaluates that
# broadcast there with `_call_bmt_1m(::LinearizedAverage, ...)` instrumented,
# recording the inputs and outputs of every CloudMicrophysics call (one per
# quadrature node that does not exit early), tagged with its cell.
using CUDA, Serialization, Printf
const OUT = abspath(ARGS[1])
mkpath(OUT)
const MODELS = "/home/treddy/ClimaAtmos.jl/perf/run-and-time/models"

# `case_setup` from the harness, without the trailing profiling lines.
src = read(joinpath(MODELS, "set_CM_cache.jl"), String)
src = first(split(src, "\nimport ClimaCore\nClimaCoreCudaExt"))
const ATMOS_CONFIG = joinpath(MODELS, "atmos_config.jl")
src = replace(
    src,
    "@__DIR__" => repr(MODELS),
    "include(\"atmos_config.jl\")" => "include($(repr(ATMOS_CONFIG)))",
)
cd(() -> include_string(Main, src, joinpath(MODELS, "set_CM_cache.jl")), MODELS)

import ClimaAtmos as CA
import ClimaCore as CC
const BMT = CA.BMT
const CAP = CA.Parameters

state = case_setup()
(; Y, p) = state
mm, tm = p.atmos.microphysics_model, p.atmos.turbconv_model
CA.set_microphysics_tendency_cache!(Y, p, mm, tm)
CUDA.synchronize()
CA.set_microphysics_tendency_cache!(Y, p, mm, tm)
CUDA.synchronize()

# Reference: the routine's kernels on this state, longest first.
prof = CUDA.@profile trace = false CA.set_microphysics_tendency_cache!(Y, p, mm, tm)
d = prof.device
rows = sort(
    [(1e3 * (d.stop[i] - d.start[i]), d.registers[i], d.local_mem[i], d.block[i], d.name[i])
     for i in eachindex(d.name) if !ismissing(d.registers[i])];
    rev = true,
)
open(joinpath(OUT, "reference_kernels.txt"), "w") do io
    for r in rows
        println(io, @sprintf("%8.3f ms  regs=%3d  local=%s  block=%s  %s", r[1], r[2], r[3], r[4], first(r[5], 200)))
    end
end
println("REFERENCE longest kernel: ", @sprintf("%.3f ms regs=%d local=%s block=%s", rows[1][1:4]...))

# Inputs of the environment broadcast (`set_microphysics_tendency_cache!` for
# `NonEquilibriumMicrophysics1M`, `PrognosticEDMFX`): the scratch fields it
# just filled and the precomputed quantities it reads.
pc, s = p.precomputed, p.scratch
ᶜw⁰ = @. CA.w_component(CA.WVec(pc.ᶜu⁰))
cpu = CC.to_cpu
ρ, T, w, qt = cpu(s.ᶜtemp_scalar), cpu(pc.ᶜT⁰), cpu(ᶜw⁰), cpu(pc.ᶜq_tot_nonneg⁰)
ql, qi, qr, qs = cpu(s.ᶜtemp_scalar_2), cpu(s.ᶜtemp_scalar_3), cpu(s.ᶜtemp_scalar_4), cpu(s.ᶜtemp_scalar_5)
λ, mu_S = cpu(s.ᶜtemp_scalar_6), cpu(s.ᶜtemp_scalar_7)
TT, qq, corr = cpu(pc.ᶜT′T′), cpu(pc.ᶜq′q′), cpu(pc.ᶜcorr_Tq)
mom, α, pf = cpu(pc.ᶜsgs_moments), cpu(pc.ᶜα), cpu(pc.ᶜprecip_frac)
gpu_out = cpu(pc.ᶜmp_tendency⁰)

dt = CA.mp_averaging_dt(p, p.mp_averaging)
nsubs_quad = mm.n_substeps_quad
sgs_quad = p.atmos.sgs_quadrature
opts = CA.sgs_microphysics_options(p.params)
cmp = CAP.microphysics_1m_params(p.params)
thp = CAP.thermodynamics_params(p.params)
println("dt = $dt ($(typeof(dt)))  nsubs_quad = $nsubs_quad  mp_averaging = $(p.mp_averaging)")
open(joinpath(OUT, "parameters.txt"), "w") do io
    println(io, "dt = $dt  nsubs_quad = $nsubs_quad  mp_averaging = $(p.mp_averaging)")
    println(io, "\nmicrophysics_model:"); show(io, MIME"text/plain"(), mm)
    println(io, "\n\nmicrophysics_1m_options:"); show(io, MIME"text/plain"(), CA.microphysics_1m_options(mm))
    println(io, "\n\nsgs_quad:"); show(io, MIME"text/plain"(), sgs_quad)
    println(io, "\n\nopts:"); show(io, MIME"text/plain"(), opts)
    println(io, "\n\ncmp:"); show(io, MIME"text/plain"(), cmp)
    println(io, "\n\nthp:"); show(io, MIME"text/plain"(), thp)
end

# Instrument the CloudMicrophysics call of every quadrature node. GPU work is
# done; this only affects code compiled from here on (the CPU re-evaluation).
const REC_CELL = Int32[]
const REC_IN = NTuple{8, Float32}[]
const REC_OUT = NTuple{4, Float32}[]
const DT_NSUBS = Set{Tuple{Float32, Int}}()
const CELL = Ref{Int32}(0)
@noinline function record!(inputs, out, dt, nsubs)
    push!(REC_CELL, CELL[])
    push!(REC_IN, Float32.(inputs))
    push!(REC_OUT, Float32.(values(out)))
    push!(DT_NSUBS, (Float32(dt), Int(nsubs)))
    return nothing
end
@eval CA @inline function _call_bmt_1m(
    mode::BMT.LinearizedAverage, scheme, mp, tps, ρ, T, w,
    q_tot, q_lcl, q_icl, q_rai, q_sno, dt, nsubs, args,
)
    out = BMT.bulk_microphysics_tendencies(
        mode, scheme, mp, tps, ρ, T, w,
        q_tot, q_lcl, q_icl, q_rai, q_sno, dt, nsubs, args...,
    )
    Main.record!((ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno), out, dt, nsubs)
    return out
end

ᶜidx = similar(ρ)
parent(ᶜidx) .= reshape(Float32.(1:length(parent(ᶜidx))), size(parent(ᶜidx)))
with_cell(i, args...) = (CELL[] = Int32(i); CA.microphysics_tendencies_1m(args...))
cpu_out = similar(gpu_out)
scheme = BMT.Microphysics1Moment()
t_cpu = @elapsed @. cpu_out = with_cell(
    ᶜidx, scheme, sgs_quad, cmp, thp, ρ, T, w, qt, ql, qi, qr, qs,
    TT, qq, corr, mom, α, dt, nsubs_quad, λ, mu_S, pf, opts,
)
a, b = parent(cpu_out), parent(gpu_out)
println("CPU re-evaluation: $(round(t_cpu; digits = 1)) s; max |cpu - gpu| / max |gpu| per tendency: ",
    [maximum(abs, a[:, :, :, f, :] .- b[:, :, :, f, :]) / max(maximum(abs, b[:, :, :, f, :]), floatmin(Float32)) for f in axes(a, 4)])
@assert length(DT_NSUBS) == 1 "dt/nsubs vary: $DT_NSUBS"

ncell = length(parent(ᶜidx))
counts = zeros(Int32, ncell)
foreach(c -> counts[c] += 1, REC_CELL)
perm = sortperm(REC_CELL; alg = Base.Sort.DEFAULT_STABLE)
offsets = Int32.(cumsum(vcat(1, counts)))
inputs = Matrix{Float32}(undef, 8, length(perm))
outputs = Matrix{Float32}(undef, 4, length(perm))
for (j, i) in enumerate(perm)
    inputs[:, j] .= REC_IN[i]
    outputs[:, j] .= REC_OUT[i]
end
(dtr, nsr) = only(DT_NSUBS)
println("cells = $ncell (parent size $(size(parent(ᶜidx)))), CM calls = $(length(REC_CELL)), ",
    "calls/cell mean = $(round(length(REC_CELL) / ncell; digits = 2)), max = $(maximum(counts)), ",
    "cells with 0 calls = $(round(100 * count(iszero, counts) / ncell; digits = 1))%")
serialize(joinpath(OUT, "env_nodes.jls"), (;
    parent_size = size(parent(ᶜidx)), offsets, inputs, outputs,
    input_names = (:ρ, :T, :w, :q_tot, :q_lcl, :q_icl, :q_rai, :q_sno),
    output_names = (:dq_lcl_dt, :dq_icl_dt, :dq_rai_dt, :dq_sno_dt),
    dt = dtr, nsubs = nsr,
    reference_ms = rows[1][1], reference_regs = rows[1][2],
))
println("wrote ", joinpath(OUT, "env_nodes.jls"), " (", round(filesize(joinpath(OUT, "env_nodes.jls")) / 2^20; digits = 1), " MiB)")
