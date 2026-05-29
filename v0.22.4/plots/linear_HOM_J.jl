import CairoMakie as MK
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.HomIceNucleation as CMH

FT = Float32
const ip = CMP.IceNucleationParameters(FT)

# Koop2000 parameterization
Koop_Δa = collect(0.26:0.0025:0.34)
KoopJ = @. CMH.homogeneous_J_cubic(ip.homogeneous, Koop_Δa) # SI units
log10KoopJ = @. log10(KoopJ .* 1e-6)

# Linear fit on Koop 2000
M = [ones(length(Koop_Δa)) Koop_Δa]
linear_coeffs = M \ log10KoopJ
new_log10J =
    [linear_coeffs[2] * Delta_a + linear_coeffs[1] for Delta_a in Koop_Δa]

# Plotting J vs Δa
fig = MK.Figure(resolution = (800, 600), fontsize = 18)
ax1 = MK.Axis(
    fig[1, 1],
    xlabel = "Δa_w [-]",
    ylabel = "log10(J), J = [cm^-3 s^-1]",
    title = "Linear fit to Koop2000",
)

MK.lines!(
    ax1,
    Koop_Δa,
    new_log10J,
    label = "log10(J) = $(linear_coeffs[2]) Δa + $(linear_coeffs[1])",
)
MK.lines!(ax1, Koop_Δa, log10KoopJ, label = "Koop 2000 Parameterization")

MK.axislegend(ax1, position = :rb)

#! format: on

MK.save("linear_HOM_J.svg", fig)
