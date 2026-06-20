import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI
using CairoMakie

const SVector = BMT.SA.SVector

# Convergence of the 2M+P3 averaging modes under substep refinement.
#
# The reference is rosenbrock_exact with a large substep count. Each mode is
# evaluated over the same Δt with a range of substep counts and the relative
# end-state error against the reference is reported.
#
# Both modes use the same consistent end-state saturation adjustment, a
# uniform-increment bisection that activates only when a substep would cross
# saturation and so vanishes as the substep length h -> 0, plus a per-substep
# positivity floor that also vanishes under refinement. Both therefore converge
# in every regime. They differ in how many substeps they need: the
# linearized-implicit rosenbrock_exact damps the stiff modes and converges at
# few substeps, while the explicit forward-Euler SubsteppedAverage needs h small
# relative to the fastest process time scale. The warm-rain panel shows both
# converging; the mixed-phase panel shows rosenbrock_exact converging at few
# substeps while SubsteppedAverage converges more slowly through the stiff regime.

const FT = Float64

tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
mp2 = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
p3 = mp2.ice.scheme

consistent_logλ(ρ, x) = P3.get_distribution_logλ(
    P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]),
)

Δt = FT(30)

# per-species scale floors so a collapsing species cannot dominate the metric
# x = [q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim]
scales = FT[4e-4, 8e7, 2.1e-3, 5e4, 8e-4, 5e5, 5e-4, 9e-7]
relerr(x, xref) = maximum(abs.(x .- xref) ./ scales)

step(mode, ρ, T, qtot, x0, logλ, n) =
    SVector{8, FT}(x0...) .+
    Δt .* SVector(
        values(
            BMT.bulk_microphysics_tendencies(mode, BMT.Microphysics2Moment(), mp2, tps,
                ρ, T, qtot, x0..., logλ, Δt, n))...,
    )

# warm-rain state: liquid and rain only, no ice
ρw, Tw, qtotw = FT(1.05), FT(288), FT(0.015)
x0w = FT[4e-4, 8e7, 2.1e-3, 5e4, 0, 0, 0, 0]
logλw = FT(-Inf)

# mixed-phase state: liquid, rain, P3 ice, and rime present near freezing
ρm, Tm, qtotm = FT(0.78), FT(273.5), FT(0.009)
x0m = FT[2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8]
logλm = consistent_logλ(ρm, x0m)

n_ref = 16384
nsubs = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]

ros_mode(n) = BMT.rosenbrock_exact(; n_substeps = n)
sub_mode = BMT.SubsteppedAverage(; limiter = BMT.EndStateSaturationAdjustment())

function errors(ρ, T, qtot, x0, logλ)
    xref = step(ros_mode(n_ref), ρ, T, qtot, x0, logλ, n_ref)
    err_sub = [relerr(step(sub_mode, ρ, T, qtot, x0, logλ, n), xref) for n in nsubs]
    err_ros = [relerr(step(ros_mode(n), ρ, T, qtot, x0, logλ, n), xref) for n in nsubs]
    return err_sub, err_ros
end

err_sub_w, err_ros_w = errors(ρw, Tw, qtotw, x0w, logλw)
err_sub_m, err_ros_m = errors(ρm, Tm, qtotm, x0m, logλm)

# ---- figure ----
fig = Figure(size = (900, 400))

for (col, (title, err_sub, err_ros)) in enumerate((
    ("2M + P3 warm rain (Δt = 30 s)", err_sub_w, err_ros_w),
    ("2M + P3 mixed phase (Δt = 30 s)", err_sub_m, err_ros_m),
))
    ax = Axis(
        fig[1, col];
        title,
        xlabel = "number of substeps",
        ylabel = "relative end-state error",
        xscale = log2,
        yscale = log10,
    )
    scatterlines!(ax, nsubs, max.(err_ros, eps(FT));
        label = "rosenbrock_exact", color = :black)
    scatterlines!(ax, nsubs, max.(err_sub, eps(FT));
        label = "SubsteppedAverage (satadj)", color = :firebrick)
    axislegend(ax; position = :lb)
end

save(joinpath(@__DIR__, "..", "averaging_modes_convergence.svg"), fig)
fig
