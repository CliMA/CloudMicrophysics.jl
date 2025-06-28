import CairoMakie as PL
PL.activate!(type = "svg")

import SpecialFunctions as SF

import Thermodynamics as TD
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

#! format: off

FT = Float64
res = 100 # plot resolution

# Plot Fig 3 from Blahak 2010
# ---------------------------

a_range = range(1e-1, stop = 30, length = res)
c1 = [P3.c₁(a) for a in a_range]
c2 = [P3.c₂(a) for a in a_range]
c3 = [P3.c₃(a) for a in a_range]
c4 = [P3.c₄(a) for a in a_range]

fig = PL.Figure(size = (1500, 1000), fontsize=22, linewidth=3)

ax1 = PL.Axis(fig[1, 1])
ax2 = PL.Axis(fig[1, 2])
ax3 = PL.Axis(fig[2, 1])
ax4 = PL.Axis(fig[2, 2])

ax1.ylabel = "c₁"
ax2.ylabel = "c₂"
ax3.ylabel = "c₃"
ax4.ylabel = "c₄"

ax3.xlabel = "a"
ax4.xlabel = "a"

PL.ylims!(ax2, 0, 2)
PL.ylims!(ax4, 0, 10)

c1l = PL.lines!(ax1, a_range,  c1)
c2l = PL.lines!(ax2, a_range,  c2)
c3l = PL.lines!(ax3, a_range,  c3)
c4l = PL.lines!(ax4, a_range,  c4)

PL.save("P3_GammaInc_c_coeffs.svg", fig)


# Plot figure 1 from Blahak 2010
# ------------------------------

fig2 = PL.Figure(size = (1000, 500), fontsize=22)
ax = PL.Axis(fig2[1,1])
ax.xlabel = "z"
ax.ylabel = "P"

z_range = range(1e-1, stop = 25, length=res)

P_a1  = [P3.Γ_lower(FT(1), z)  / SF.gamma(FT(1))  for z in z_range]
P_a3  = [P3.Γ_lower(FT(3), z)  / SF.gamma(FT(3))  for z in z_range]
P_a5  = [P3.Γ_lower(FT(5), z)  / SF.gamma(FT(5))  for z in z_range]
P_a10 = [P3.Γ_lower(FT(10), z) / SF.gamma(FT(10)) for z in z_range]
P_a15 = [P3.Γ_lower(FT(15), z) / SF.gamma(FT(15)) for z in z_range]
P_a30 = [P3.Γ_lower(FT(30), z) / SF.gamma(FT(30)) for z in z_range]

G_a1   = [SF.gamma_inc(FT(1), z)[1]  for z in z_range]
G_a3   = [SF.gamma_inc(FT(3), z)[1]  for z in z_range]
G_a5   = [SF.gamma_inc(FT(5), z)[1]  for z in z_range]
G_a10  = [SF.gamma_inc(FT(10), z)[1] for z in z_range]
G_a15  = [SF.gamma_inc(FT(15), z)[1] for z in z_range]
G_a30  = [SF.gamma_inc(FT(30), z)[1] for z in z_range]

l1  = PL.lines!(ax, z_range, P_a1,  label="a=1",  color="midnightblue", linewidth=3)
l3  = PL.lines!(ax, z_range, P_a3,  label="a=3",  color="slateblue",    linewidth=3)
l5  = PL.lines!(ax, z_range, P_a5,  label="a=5",  color="steelblue1",   linewidth=3)
l10 = PL.lines!(ax, z_range, P_a10, label="a=10", color="orange",       linewidth=3)
l15 = PL.lines!(ax, z_range, P_a15, label="a=15", color="orangered2",   linewidth=3)
l30 = PL.lines!(ax, z_range, P_a30, label="a=30", color="brown",        linewidth=3)

g1  = PL.lines!(ax, z_range, G_a1,  label="a=1", color="midnightblue", linewidth=3, linestyle=:dash)
g3  = PL.lines!(ax, z_range, G_a3,  label="a=1", color="slateblue",    linewidth=3, linestyle=:dash)
g5  = PL.lines!(ax, z_range, G_a5,  label="a=1", color="steelblue1",   linewidth=3, linestyle=:dash)
g10 = PL.lines!(ax, z_range, G_a10, label="a=1", color="orange",       linewidth=3, linestyle=:dash)
g15 = PL.lines!(ax, z_range, G_a15, label="a=1", color="orangered2",   linewidth=3, linestyle=:dash)
g30 = PL.lines!(ax, z_range, G_a30, label="a=1", color="brown",        linewidth=3, linestyle=:dash)

PL.Legend(
  fig2[1, 2],
  [l1, l3, l5, l10, l15, l30],
  ["a=1", "a=3", "a=5", "a=10", "a=15", "a=30"],
)
PL.save("P3_GammaInc_P_check.svg", fig2)

# Plot figure 4 from Blahak 2010
# ------------------------------
# The relative error blows up for cases where the Γ_lower is very small.
# The reference value we obtain from Julia Special Functions is off
# by several orders of magnitude when compared with Wolfram Mathematica

y_range = range(1e-1, stop = 2.5, length = res)
z_range = y_range .* (a_range .+ 1)

val     = zeros(res, res)
err_rel = zeros(res, res)
err_abs = zeros(res, res)
for it in range(1, stop=res, step=1)      # a
    for jt in range(1, stop=res, step=1)  # z
        a = FT(a_range[it])
        z = FT(z_range[jt])
        val[it, jt] = P3.Γ_lower(a, z) / SF.gamma(a)
        err_abs[it, jt] = (P3.Γ_lower(a, z) / SF.gamma(a)) - SF.gamma_inc(a, z)[1]
        err_rel[it, jt] = (P3.Γ_lower(a, z) / SF.gamma(a)) / SF.gamma_inc(a, z)[1] - 1
    end
end

fig3 = PL.Figure(size = (1800, 600), fontsize=22)
ax1 = PL.Axis(fig3[1, 1], xlabel = "a", ylabel = "z  / (a + 1)", title = "Relative error")
ax2 = PL.Axis(fig3[1, 3], xlabel = "a", title = "Absolute error")
ax3 = PL.Axis(fig3[1, 5], xlabel = "a", title = "Γ_lower(a, z) / Γ(a)")

PL.ylims!(ax1, 0, 2.5)
PL.ylims!(ax2, 0, 2.5)
PL.ylims!(ax3, 0, 2.5)

cplot1 = PL.heatmap!(
    ax1, a_range, y_range, err_rel,
    colormap = PL.cgrad(:viridis, 20, categorical=true),
    colorrange = (-0.1, 10), highclip=:red, lowclip=:orange,
)
PL.Colorbar(fig3[1,2], cplot1)

cplot2 = PL.contourf!(
    ax2, a_range, y_range, err_abs,
    levels=-0.03:0.005:0.03,
    colormap="RdBu"
)
PL.Colorbar(fig3[1,4], cplot2)

cplot3 = PL.contourf!(
    ax3, a_range, y_range, val,
)
PL.Colorbar(fig3[1,6], cplot3)

PL.save("P3_GammaInc_error.svg", fig3)
