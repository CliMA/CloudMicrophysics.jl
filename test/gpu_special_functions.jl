# GPU harness for `ADSpecialFunctions`.
#
# Backend-agnostic via KernelAbstractions. The point of this module is that
# `SpecialFunctions` compiles slowly on GPU (and `gamma_inc_inv` not at all),
# so this harness (1) confirms every function gives the host answer inside a
# GPU kernel, including through a `ForwardDiff.Dual`, and (2) reports compile
# (first launch) and steady-state time.
#
# Because the `Float32` path is genuine `Float32` arithmetic (no internal
# promotion to `Float64`), it runs on Float32-only GPUs such as Apple Metal:
#
#     julia --project=test -e 'using Metal; include("test/gpu_special_functions.jl"); \
#         run_gpu_special_functions(MetalBackend(), MtlArray)'
#     julia --project=test -e 'using CUDA;  include("test/gpu_special_functions.jl"); \
#         run_gpu_special_functions(CUDABackend(), CuArray)'
#
# Included directly (no GPU), it runs a CPU-backend smoke test.

import Test as TT
using KernelAbstractions
import CloudMicrophysics.ADSpecialFunctions as ADSF
import SpecialFunctions as SF
import ForwardDiff as FD

# One kernel per function. Scalar, branch-light bodies keep the GPU compile small.
@kernel inbounds = true function _gamma_inc_kernel!(out, a, x)
    i = @index(Global, Linear)
    out[i] = ADSF.gamma_inc(a[i], x[i])[1]
end
@kernel inbounds = true function _gamma_inc_inv_kernel!(out, a, p)
    i = @index(Global, Linear)
    out[i] = ADSF.gamma_inc_inv(a[i], p[i])
end
@kernel inbounds = true function _beta_inc_kernel!(out, a, b, x)
    i = @index(Global, Linear)
    out[i] = ADSF.beta_inc(a[i], b[i], x[i])[1]
end
@kernel inbounds = true function _erf_kernel!(out, x)
    i = @index(Global, Linear)
    out[i] = ADSF.erf(x[i])
end
# AD inside a kernel: d/dx P(a,·) carried by a Dual argument.
@kernel inbounds = true function _dPdx_kernel!(out, a, x)
    i = @index(Global, Linear)
    out[i] = FD.derivative(xx -> ADSF.gamma_inc(a[i], xx)[1], x[i])
end

function _launch!(backend, kern!, out, args...)
    kern!(backend)(out, args...; ndrange = length(out))
    KernelAbstractions.synchronize(backend)
    out
end

_line(name, tc, tr, err, tol) = println(
    "  $(rpad(name, 15)) compile=$(round(tc; digits = 2))s  " *
    "run≈$(round(tr * 1e3; digits = 3))ms  maxrelerr=$(round(err; sigdigits = 2))  (tol $tol)",
)

"""
    run_gpu_special_functions(backend = CPU(), ArrayType = Array; FT = Float32, n = 4096)

Launch each `ADSpecialFunctions` kernel on `backend`, check against a host
reference, and print compile/run timings. `ArrayType` is the device array
constructor (`Array`, `CuArray`, `MtlArray`).
"""
function run_gpu_special_functions(backend = CPU(), ArrayType = Array; FT = Float32, n = 4096)
    println("backend = $(backend), ArrayType = $ArrayType, FT = $FT, n = $n")
    dev(v) = ArrayType(FT.(v))

    # ranges representative of the gamma-PSD use in CloudMicrophysics
    ah = collect(range(0.3, 25.0, length = n))
    xh = collect(range(0.05, 40.0, length = n))
    bh = collect(range(0.4, 12.0, length = n))
    ph = collect(range(0.01, 0.99, length = n))
    xeh = collect(range(-5.0, 5.0, length = n))

    a, x, b, p, xe = dev(ah), dev(xh), dev(bh), dev(ph), dev(xeh)
    xbeta = dev(xh ./ 50)
    out = ArrayType(zeros(FT, n))

    relerr(v, r) = abs(r) < eps(FT) ? abs(v) : abs(v - r) / abs(r)
    maxrelerr(o, ref) = maximum(relerr.(Float64.(Array(o)), ref))

    cases = (
        ("gamma_inc P", () -> _launch!(backend, _gamma_inc_kernel!, out, a, x),
            [SF.gamma_inc(ah[i], xh[i])[1] for i in 1:n]),
        ("gamma_inc_inv", () -> _launch!(backend, _gamma_inc_inv_kernel!, out, a, p),
            [SF.gamma_inc_inv(ah[i], ph[i], 1 - ph[i]) for i in 1:n]),
        ("beta_inc I", () -> _launch!(backend, _beta_inc_kernel!, out, a, b, xbeta),
            [SF.beta_inc(ah[i], bh[i], xh[i] / 50)[1] for i in 1:n]),
        ("erf", () -> _launch!(backend, _erf_kernel!, out, xe),
            [SF.erf(xeh[i]) for i in 1:n]),
        ("dP/dx (AD)", () -> _launch!(backend, _dPdx_kernel!, out, a, x),
            [xh[i]^(ah[i] - 1) * exp(-xh[i]) / SF.gamma(ah[i]) for i in 1:n]),
    )

    TT.@testset "GPU special functions ($backend, $FT)" begin
        for (name, run, ref) in cases
            t_compile = @elapsed run()                 # first launch: includes compile
            t_run = @elapsed for _ in 1:10
                run()
            end
            t_run /= 10
            err = maxrelerr(out, ref)
            # F32 normal-range floor (~1e-5) plus headroom for tail/underflow points
            tol = FT == Float32 ? 5e-3 : 1e-9
            _line(name, t_compile, t_run, err, tol)
            TT.@test err < tol
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_gpu_special_functions(CPU(), Array)
end
