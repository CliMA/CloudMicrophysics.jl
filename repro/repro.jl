# MVP: GPU kernel-compile time, raw CuArray broadcast vs ClimaCore Field broadcast.
#
# The same scalar per-point function is broadcast two ways on the GPU:
#   (1) over plain CUDA.CuArrays            -- "raw" GPU broadcast
#   (2) over a ClimaCore CenterFiniteDifferenceSpace Field
# and the compile time of each is measured. Run with the GPU device:
#   CLIMACOMMS_DEVICE=CUDA julia --project=<env-with ClimaCore,CUDA,SpecialFunctions,CloudMicrophysics> mvp_compile.jl
#
# Measurement rigor:
#   * warmup()  -- a trivial throwaway broadcast over both backends compiles all the
#                  one-time infrastructure (CUDA runtime/context, GPUCompiler, ClimaCore
#                  Field-broadcast machinery, LLVM/ptxas) so the timed calls pay only the
#                  target function's kernel compile, not "regular" package compilation.
#   * run-twice -- compile = first_call - second_call, so kernel exec is separated out too
#                  (it is microseconds at this size; the array length does not affect compile).

import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC
import CUDA
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.ThermodynamicsInterface as TDI
import SpecialFunctions as SF
using Printf

const FT = Float64
const N = 1000
make_column() =
    CC.CommonGrids.ColumnGrid(FT; z_elem = N, z_min = 0, z_max = 1000) |>
    CC.Spaces.CenterFiniteDifferenceSpace

# (1) light synthetic, no CloudMicrophysics
@inline syn_light(ρ, T, q1, q2) = ρ * q1 + T * q2 + exp(-q1) + log1p(q2)

# (2) gamma_inc in a runtime LOOP -- one call site
@inline function syn_gammaloop(ρ, T, q1, q2)
    acc = zero(ρ)
    @inbounds for i in 1:12
        a = q1 * FT(i) + one(ρ) + T
        x = q2 * FT(10) * FT(i) + FT(i)
        (p, qq) = SF.gamma_inc(a, x)
        acc += SF.gamma(a) * (p - qq) * exp(-x) + SF.loggamma(a)
    end
    acc + log1p(abs(acc))
end

# (3) gamma_inc unrolled -- 16 distinct call sites (same work as the loop)
@inline function syn_gammaunroll(ρ, T, q1, q2)
    acc = zero(ρ)
    Base.Cartesian.@nexprs 16 i -> begin
        a_i = q1 * FT(i) + one(ρ) + T
        x_i = q2 * FT(10) * FT(i) + FT(i)
        pq_i = SF.gamma_inc(a_i, x_i)
        acc += SF.gamma(a_i) * (pq_i[1] - pq_i[2]) * exp(-x_i) + SF.loggamma(a_i)
    end
    acc + log1p(abs(acc))
end

# (4) real CloudMicrophysics one-moment (1M) bulk tendency, reduced to a scalar
const MP1 = CMP.Microphysics1MParams(FT)
const TPS = TDI.TD.Parameters.ThermodynamicsParameters(FT)
const RHO = FT(1.2)
const QT = FT(0.015)
const TEMP = TDI.TD.Parameters.T_freeze(TPS) - FT(5)
@inline g1m(q_liq, q_ice, q_rai, q_sno) = sum(
    values(
        BMT.bulk_microphysics_tendencies(
            BMT.Instantaneous(), BMT.Microphysics1Moment(), MP1, TPS,
            RHO, TEMP, QT, q_liq, q_ice, q_rai, q_sno,
        ),
    ),
)

# trivial function for warmup -- distinct from all targets so it pre-compiles only infra
@inline warm(a, b, c, d) = a + b + c + d

function warmup()
    ca = ntuple(_ -> CUDA.fill(FT(1e-4), N), 4)
    Base.broadcast(warm, ca...)
    CUDA.synchronize()
    sp = make_column()
    fl = ntuple(_ -> (x = zeros(sp); fill!(parent(x), FT(1e-4)); x), 4)
    Base.broadcast(warm, fl...)
    CUDA.synchronize()
    return nothing
end

# infra warm => the first call below pays only f's kernel compile; run-twice strips exec
function measure(label, f)
    sp = make_column()
    fl = ntuple(_ -> (x = zeros(sp); fill!(parent(x), FT(1e-4)); x), 4)
    GC.gc()
    t = time()
    Base.broadcast(f, fl...)
    CUDA.synchronize()
    f1 = time() - t
    GC.gc()
    t = time()
    Base.broadcast(f, fl...)
    CUDA.synchronize()
    f2 = time() - t
    ca = ntuple(_ -> CUDA.fill(FT(1e-4), N), 4)
    GC.gc()
    t = time()
    Base.broadcast(f, ca...)
    CUDA.synchronize()
    cu1 = time() - t
    GC.gc()
    t = time()
    Base.broadcast(f, ca...)
    CUDA.synchronize()
    cu2 = time() - t
    @printf("  %-13s  CuArray %6.1f s (exec %.4f)   Field %6.1f s (exec %.4f)   ratio %5.1fx\n",
        label, cu1 - cu2, cu2, f1 - f2, f2, (f1 - f2) / max(cu1 - cu2, FT(1e-6)))
    flush(stdout)
end

println("device = ", ClimaComms.device())
println("warming up infrastructure (CUDA / GPUCompiler / ClimaCore) with a trivial kernel...")
warmup()
println("=== kernel compile, infra-separated & exec-separated (compile = 1st - 2nd call) ===")
measure("1M", g1m)
measure("gamma_loop", syn_gammaloop)
measure("syn_light", syn_light)
# measure("gamma_unroll", syn_gammaunroll)
println("DONE")
