import CairoMakie: CairoMakie, Makie
CairoMakie.activate!(type = "svg")

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Common as CO
import CloudMicrophysics.DistributionTools as DT
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.ThermodynamicsInterface as TDI
FT = Float64

ice = CMP.Microphysics2MParams(FT; with_ice = true).ice
aps = CMP.AirProperties(FT)
tps = TDI.PS(FT)

# The collected cloud and rain mass rates and the maximum freezing rate of one ice particle, against its
# diameter, and the ends of the wet set, for a rimed ice population collecting cloud and some rain at -10 °C
function wet_set_curves(ice, aps, tps)
    vel = ice.terminal_velocity
    ρₐ, T = FT(1), FT(263.15)
    L_ice, N_ice, F_rim, ρ_rim = FT(1e-3), FT(1e3), FT(0.5), FT(500)
    L_c, N_c, L_r, N_r = FT(1e-3), FT(1e8), FT(1e-4), FT(1e3)

    state = P3.P3State(ice.scheme, ice.integration, L_ice, N_ice, F_rim * L_ice, F_rim * L_ice / ρ_rim)
    (; particles, dist) = state
    ∂ₜV = P3.volumetric_collision_rate_integrand(vel, particles, ρₐ)
    ρ′_rim = P3.compute_local_rime_density(vel, particles, ρₐ, T)
    ∂ₜM_max = P3.compute_max_freeze_rate(aps, tps, vel, particles, ρₐ, T)
    m_liq(D) = ice.cloud_pdf.ρw * CO.volume_sphere_D(D)
    p = P3.liquid_ice_collision_quantile(FT)
    function collected(pdf, L, N)
        n = DT.size_distribution(pdf, L / ρₐ, ρₐ, N)
        bounds = CM2.get_size_distribution_bounds(pdf, L / ρₐ, ρₐ, N, p)
        return P3.get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, bounds, ice.integration)
    end
    cloud = collected(ice.cloud_pdf, L_c, N_c)
    rain = collected(ice.rain_pdf, L_r, N_r)

    bounds = P3.velocity_integral_bounds(particles, dist, P3.ice_particle_terminal_velocity(vel, particles, ρₐ); p)
    D_lo, D_hi = first(bounds), last(bounds)
    D_wet = P3.wet_growth_onset_diameters(cloud, rain, ∂ₜM_max, D_lo, D_hi)
    Ds = exp.(range(log(D_lo), log(D_hi); length = 400))
    return (; Ds, u_c = [cloud(D)[2] for D in Ds], u_r = [rain(D)[2] for D in Ds], m = ∂ₜM_max.(Ds), D_wet)
end
(; Ds, u_c, u_r, m, D_wet) = wet_set_curves(ice, aps, tps)

fig = Makie.with_theme(Makie.theme_minimal()) do
    fig = Makie.Figure(size = (600, 380), figure_padding = 20)
    ax = Makie.Axis(fig[1, 1];
        title = "Wet set of one ice population", xscale = log10, yscale = log10,
        xlabel = "ice particle diameter D (m)", ylabel = "mass rate per particle (kg/s)",
    )
    Makie.vspan!(ax, extrema(D_wet)...; color = (:royalblue, 0.12), label = "wet set")
    Makie.lines!(ax, Ds, u_c; color = :gray, linestyle = :dash, label = "collected cloud")
    Makie.lines!(ax, Ds, u_r; color = :gray, linestyle = :dot, label = "collected rain")
    Makie.lines!(ax, Ds, u_c .+ u_r; color = :black, label = "collected liquid u + v")
    Makie.lines!(ax, Ds, m; color = :darkorange, label = "maximum freezing rate m")
    Makie.axislegend(ax; position = :lt)
    fig
end
Makie.save("P3WetSet.svg", fig)
