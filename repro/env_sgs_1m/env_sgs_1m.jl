# Minimal reproducer of the slowest GPU kernel of ClimaAtmos'
# `perf/run-and-time/models/set_CM_cache.jl`: the environment SGS-quadrature 1M
# broadcast of `set_microphysics_tendency_cache!(Y, p, ::NonEquilibriumMicrophysics1M,
# ::PrognosticEDMFX)`.
#
# In ClimaAtmos each thread handles one cell. It builds a conditional quadrature
# rule over the subgrid PDF (up to 8 nodes) and, at each node that does not exit
# early, calls a `@noinline` evaluator whose body is
#
#     BMT.bulk_microphysics_tendencies(BMT.LinearizedAverage(), BMT.Microphysics1Moment(),
#         mp, tps, ρ, T_node, w, q_tot_node, q_lcl_node, q_icl_node, q_rai_node, q_sno_node,
#         dt, nsubs)
#
# with the parameters folded into the kernel as constants (`Val(mp)`, `Val(tps)`).
# This script keeps only that CloudMicrophysics part: one thread per cell, a loop
# over the cell's recorded node inputs, and the same `@noinline` barrier around the
# same call. The inputs of every such call were recorded from the ClimaAtmos state
# after the 40-step spin-up of the ClimaCoupler benchmark `amip_progedmf_1m_land_he16`
# (dt = 30 s, Float32), with ClimaAtmos' outputs for validation, by
# `capture_env_nodes.jl` (run in ClimaAtmos' `.buildkite` environment). Cells are
# ordered as the parent array, `(v, i, j, h)` with `v` fastest. `env_nodes.jls` here
# keeps every 8th horizontal element of the 1536 (242688 of the 1941504 cells);
# `parameters.txt` holds the parameters ClimaAtmos used.
#
# Usage (needs a GPU):
#     julia --project=repro/env_sgs_1m repro/env_sgs_1m/env_sgs_1m.jl
# Environment variables: DATA sets the input file (default `env_nodes.jls` here);
# ALWAYS_INLINE=false compiles without the `always_inline = true` that ClimaCore
# uses for its kernels; FOLD=false passes the parameters as kernel arguments
# instead of compile-time constants; NSAMPLES sets the number of timed launches.
import CUDA
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.ThermodynamicsInterface as TDI
using Printf, Serialization

const FT = Float32
const FOLD = parse(Bool, get(ENV, "FOLD", "true"))
const ALWAYS_INLINE = parse(Bool, get(ENV, "ALWAYS_INLINE", "true"))
const NSAMPLES = parse(Int, get(ENV, "NSAMPLES", "20"))

data = deserialize(get(ENV, "DATA", joinpath(@__DIR__, "env_nodes.jls")))
(; offsets, inputs, outputs, dt, nsubs) = data
ncell = length(offsets) - 1

# Parameters as ClimaAtmos builds them for this configuration
# (`parameters.txt` holds the values recorded from the ClimaAtmos run). Cloud
# liquid and ice formation are off: ClimaAtmos treats them in a separate
# implicit call, so this bulk set leaves them out.
toml_dict = CP.create_toml_dict(FT; override_file = joinpath(@__DIR__, "amip_progedmf_1m.toml"))
mp = CMP.Microphysics1MParams(
    toml_dict;
    cloud_liquid_formation = nothing,
    cloud_ice_formation = nothing,
    cloud_ice_melt = CMP.CloudIceMelt(),
    cloud_liquid_freezing = CMP.HomogeneousAndHeterogeneous(),
    rain_autoconversion = CMP.Kessler1M(),
    snow_autoconversion = CMP.NoSupersaturation(),
    rain_condensation_evaporation = CMP.RainEvaporation(),
    snow_deposition_sublimation = CMP.DepositionAndSublimation(),
    snow_melt = CMP.SnowMelt(),
    cloud_liquid_rain_accretion = CMP.CloudLiquidRainAccretion(),
    cloud_liquid_snow_accretion = CMP.CloudLiquidSnowAccretion(),
    cloud_ice_rain_accretion = CMP.CloudIceRainAccretion(),
    cloud_ice_snow_accretion = CMP.CloudIceSnowAccretion(),
    rain_snow_accretion = CMP.RainSnowAccretion(),
    vapor_relaxation = CMP.JointVaporRelaxation(),
)
tps = TDI.TD.Parameters.ThermodynamicsParameters(toml_dict)

@inline unval(::Val{x}) where {x} = x
@inline unval(x) = x

# The per-node function barrier of ClimaAtmos' `Microphysics1MEvaluator`. As
# there, the parameters cross the barrier still wrapped in `Val` and are
# unwrapped inside it, so the call is compiled with them as constants.
@noinline node_tendencies(mp, tps, ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno, dt, nsubs) =
    BMT.bulk_microphysics_tendencies(
        BMT.LinearizedAverage(), BMT.Microphysics1Moment(), unval(mp), unval(tps),
        ρ, T, w, q_tot, q_lcl, q_icl, q_rai, q_sno, dt, nsubs,
    )

@inline function node_at(mp, tps, inputs, k, dt, nsubs)
    @inbounds node_tendencies(
        mp, tps,
        inputs[1, k], inputs[2, k], inputs[3, k], inputs[4, k],
        inputs[5, k], inputs[6, k], inputs[7, k], inputs[8, k], dt, nsubs,
    )
end

# One thread per cell, looping over the cell's nodes (the reproducer).
function cell_kernel!(out, offsets, inputs, mp, tps, dt, nsubs)
    i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if i < length(offsets)
        a = b = c = d = zero(eltype(out))
        @inbounds for k in offsets[i]:(offsets[i + 1] - 1)
            r = node_at(mp, tps, inputs, k, dt, nsubs)
            a += r.dq_lcl_dt
            b += r.dq_icl_dt
            c += r.dq_rai_dt
            d += r.dq_sno_dt
        end
        @inbounds out[1, i], out[2, i], out[3, i], out[4, i] = a, b, c, d
    end
    return nothing
end

# One thread per node, for validating each call against ClimaAtmos' output.
function node_kernel!(out, inputs, mp, tps, dt, nsubs)
    k = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if k <= size(inputs, 2)
        r = node_at(mp, tps, inputs, k, dt, nsubs)
        @inbounds out[1, k], out[2, k], out[3, k], out[4, k] =
            r.dq_lcl_dt, r.dq_icl_dt, r.dq_rai_dt, r.dq_sno_dt
    end
    return nothing
end

function launch(kern, n, args...)
    threads = min(n, CUDA.launch_configuration(kern.fun).threads)
    kern(args...; threads, blocks = cld(n, threads))
end

function report(name, f, n, args)
    kern = CUDA.@cuda launch = false always_inline = ALWAYS_INLINE f(args...)
    launch(kern, n, args...)
    CUDA.synchronize()
    times = [CUDA.@elapsed(launch(kern, n, args...)) for _ in 1:NSAMPLES]
    sort!(times)
    threads = min(n, CUDA.launch_configuration(kern.fun).threads)
    @printf("%-12s median %8.3f ms  min %8.3f ms  regs %3d  local %4d B  threads %d  blocks %d\n",
        name, 1e3 * times[cld(end, 2)], 1e3 * times[1], CUDA.registers(kern),
        CUDA.memory(kern).local, threads, cld(n, threads))
    return kern
end

fmp, ftps = FOLD ? (Val(mp), Val(tps)) : (mp, tps)
d_offsets, d_inputs = CUDA.CuArray(offsets), CUDA.CuArray(inputs)
nnode = size(inputs, 2)
println(CUDA.device(), "; parameters folded into the kernel: ", FOLD, "; always_inline: ", ALWAYS_INLINE)
@printf("%d cells, %d CloudMicrophysics calls (%.2f per cell, %.1f%% of cells with none), dt = %g s, nsubs = %d\n",
    ncell, nnode, nnode / ncell, 100 * count(i -> offsets[i + 1] == offsets[i], 1:ncell) / ncell, dt, nsubs)
haskey(data, :subsample) && println("Subsample: ", data.subsample)
@printf("ClimaAtmos kernel on the full state: %.3f ms, %d registers\n", data.reference_ms, data.reference_regs)

d_cell_out = CUDA.zeros(FT, 4, ncell)
report("per cell", cell_kernel!, ncell, (d_cell_out, d_offsets, d_inputs, fmp, ftps, FT(dt), nsubs))
d_node_out = CUDA.zeros(FT, 4, nnode)
report("per node", node_kernel!, nnode, (d_node_out, d_inputs, fmp, ftps, FT(dt), nsubs))

# Each tendency against ClimaAtmos' value of the same call (computed on the CPU),
# relative to the largest magnitude of that tendency. The GPU and CPU round
# differently, which matters at a few nodes near process thresholds.
node_out = Array(d_node_out)
for (j, name) in enumerate(data.output_names)
    scale = max(maximum(abs, view(outputs, j, :)), floatmin(FT))
    rel = abs.(view(node_out, j, :) .- view(outputs, j, :)) ./ scale
    @printf("%-10s max |Δ| / max |ClimaAtmos| = %.2e; calls above 1e-3: %d, above 1e-6: %d\n",
        name, maximum(rel), count(>(1e-3), rel), count(>(1e-6), rel))
end
