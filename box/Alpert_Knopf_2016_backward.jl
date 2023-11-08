import CairoMakie as MK
import Dierckx as DX
import Random as RD
import Distributions as DS
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CMO
import CloudMicrophysics.HetIceNucleation as CMI_het

FT = Float64

# AK 2016 data and the frozen fraction fit to data
include(joinpath(pkgdir(CM), "box", "Alpert_Knopf_2016_data.jl"))
AK16_ff = DX.Spline1D(AK16_T_frozen_fraction, AK16_frozen_fraction, k = 1)

# Eq 11
# n_frz  - number of freezing events that occur in a given time interval
# N_ufrz - number of unforzen droplets at temperature T
function apparent_J(n_frz, N_ufrz, Ag, dt)
    # Need some way of preventing divisions by 0
    return (N_ufrz < 1 || n_frz == 0) ? 1e-10 : n_frz / N_ufrz / Ag / dt
end
# Eq 12
# A_sum - total surface area from the droplets that are still liquid
function actual_J(n_frz, A_sum, dt)
    # Need some way of preventing divisions by 0
    return A_sum < 1e-11 ? 1 : n_frz / A_sum / dt
end

# temperature range and paper based frozen fraction
T = range(260, 230, 50)
ff = AK16_ff(T)
# initial number of liquid droplets
N₀ = 1000
# cooling rate and fake model time step
cooling_rate = FT(0.5 / 60) # K s^-1
dt = (T[1] - T[2]) / cooling_rate

# compute n_frz from frozen fraction
n_frz = similar(T)
n_frz[1] = FT(0)
for it in range(2, length(n_frz) - 1)
    n_frz[it] = (ff[it] - ff[it - 1]) * N₀
end

# compute the N_ufrz from frozen fraction
N_ufrz = @. (1 - ff) * N₀

# assumed distribution of available surface area
Ag = 1e-9
σ = 10
A_distr = DS.LogNormal(log(Ag), log(σ))
# generate a sample of available areas
function create_Aj(A_distr, N₀)
    A_sample = DS.rand(A_distr, N₀)
    Aj = zeros(N₀)
    for it in range(N₀, stop = 1, step = -1)
        Aj[it] = A_sample[it]
    end
    return sort(Aj, rev = true)
end
Aj_sorted = zeros(N₀)
Aj_sorted = create_Aj(A_distr, N₀)

# Compute the sum of Aj from the unfrozen droplets
# assuming they will be freezing from the largest to the smallest.
# There is probably a better way to do it, here I'm just zeroing out
# areas going from largest to smallest based on the numer of freezing events
# that should occurr.
function compute_Asum!(A_sum, Aj_sorted, T)
    track = 0
    for it in range(start = 1, stop = length(T), step = 1)
        nf = Int64(floor(n_frz[it]))
        if nf > 0
            for j in range(start = 1, stop = nf, step = 1)
                idx = min(track + j, length(Aj_sorted))
                Aj_sorted[idx] = FT(0)
            end
        end
        track += nf
        A_sum[it] = sum(Aj_sorted)
    end
end
A_sum = similar(T)
compute_Asum!(A_sum, Aj_sorted, T)

# Compute the immersion freezing rate from CloudMicrophysics
tps = CMP.ThermodynamicsParameters(FT)
aerosol = CMP.Illite(FT)
Δa = @. FT(1) - CMO.a_w_ice(tps, T)
J_immer = @. CMI_het.ABIFM_J(aerosol, Δa) # m^-2 s^-1

#! format: off
# Plot results
fig = MK.Figure(resolution = (1200, 400))
ax1 = MK.Axis(fig[1, 1], ylabel = "J_immer [cm^-2 s^-1]", xlabel = "Temperature [K]", yscale = log10)
ax2 = MK.Axis(fig[1, 2], ylabel = "frozen fraction",      xlabel = "Temperature [K]")
ax3 = MK.Axis(fig[1, 3], ylabel = "total A [cm^2]",       xlabel = "Temperature [K]", yscale = log10)
MK.ylims!(ax1, 1e-4, 1e10)

MK.lines!(ax1, T, apparent_J.(n_frz, N_ufrz, Ag, dt) .* 1e-4, color = :orange, label = "apparent J", linewidth = 3)
MK.lines!(ax1, T, actual_J.(n_frz, A_sum, dt)        .* 1e-4, color = :green3, label = "actual J",   linewidth = 3)
MK.lines!(ax1, T, J_immer                            .* 1e-4, color = :red,    label = "ABIFM",      linewidth = 3)

MK.lines!(ax2, T, ff, color = :black, label = "prescribed", linewidth = 3)

MK.lines!(ax3, T[N_ufrz .> 0], N_ufrz[N_ufrz .> 0] .* Ag .* 1e4, color = :orange, label = "const A", linewidth = 3)
MK.lines!(ax3, T[A_sum .> 0], A_sum[A_sum .> 0] .* 1e4, color = :green3, label = "variable A", linewidth = 3)

MK.axislegend(ax1, framevisible = false, labelsize = 14, orientation = :horizontal, nbanks = 2, position = :rt)
MK.axislegend(ax2, framevisible = false, labelsize = 14, orientation = :horizontal, nbanks = 2, position = :rt)
MK.axislegend(ax3, framevisible = false, labelsize = 14, orientation = :horizontal, nbanks = 2, position = :lt)
MK.save("Alpert_Knopf_2016_backward.svg", fig)
#! format: on
