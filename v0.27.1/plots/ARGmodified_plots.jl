import OrdinaryDiffEq as ODE
import CairoMakie as MK

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))

function run_parcel_model(Nₐ, Nₗ, Nᵢ, rₗ, rᵢ, w, FT)

    # Get free parameters
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    wps = CMP.WaterProperties(FT)

    # Initial conditions
    T₀ = FT(230)
    cᵥ₀ = FT(5 * 1e-5)
    ln_INPC = FT(0)

    # Constants
    ρₗ = wps.ρw
    ρᵢ = wps.ρi
    ϵₘ = TDI.Rd_over_Rv(tps)
    eₛ = TDI.saturation_vapor_pressure_over_liquid(tps, T₀)
    qᵥ = ϵₘ / (ϵₘ - 1 + 1 / cᵥ₀)
    Sₗ = FT(1.0)
    e = Sₗ * eₛ
    p₀ = e / cᵥ₀
    ρ_air = TDI.air_density(tps, T₀, p₀, qᵥ, FT(0), FT(0))
    qₗ = Nₗ * FT(4 / 3 * π) * (rₗ)^3 * ρₗ / ρ_air
    qᵢ = Nᵢ * FT(4 / 3 * π) * (rᵢ)^3 * ρᵢ / ρ_air
    IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC, qₗ, Nₗ]

    # Simulation parameters passed into ODE solver
    const_dt = FT(1e-2)    # model timestep
    t_max = FT(60)     # total time
    aerosol = CMP.Sulfate(FT)

    condensation_growth = "Condensation"
    deposition_growth = "Deposition"
    aerosol_act = "AeroAct"     # turn on aerosol activation
    aero_σ_g = FT(2.3)
    r_nuc = FT(4e-8)

    params = parcel_params{FT}(
        w = w,
        const_dt = const_dt,
        aerosol_act = aerosol_act,
        aerosol = aerosol,
        aero_σ_g = aero_σ_g,
        r_nuc = r_nuc,
        condensation_growth = condensation_growth,
        deposition_growth = deposition_growth,
        Nₐ = Nₐ,
        liq_size_distribution = "MonodisperseMix",
    )

    # solve ODE
    sol = run_parcel(IC, FT(0), t_max, params)

    S_max = maximum(sol[1, :]) - FT(1)

    # ARG results
    ad = AM.Mode_κ(
        params.r_nuc,
        params.aero_σ_g,
        params.Nₐ,
        (FT(1.0),),
        (FT(1.0),),
        (params.aerosol.M,),
        (params.aerosol.κ,),
    )
    all_ad = AM.AerosolDistribution((ad,))
    S_max_ARG = AA.max_supersaturation(params.aap, all_ad, params.aps, params.tps, T₀, p₀, w, qᵥ + qₗ + qᵢ, qₗ, qᵢ)
    S_max_mod =
        AA.max_supersaturation(params.aap, all_ad, params.aps, params.tps, T₀, p₀, w, qᵥ + qₗ + qᵢ, qₗ, qᵢ, Nₗ, Nᵢ)

    return S_max / S_max_ARG, S_max_mod / S_max_ARG
end


FT = Float32
Nₐ = FT(5e8)
Nₗ = FT(1e8)
Nᵢ = FT(1e6)
rₗ = FT(20e-6)
rᵢ = FT(20e-6)
w = FT(1.2)   # updraft speed

n_points = 10
Nₗ_range = range(FT(0), stop = FT(1e8), length = n_points)
Nᵢ_range = range(FT(0), stop = FT(1e6), length = n_points)
rₗ_range = range(FT(0e-6), stop = FT(25e-6), length = n_points)
rᵢ_range = range(FT(0e-6), stop = FT(25e-6), length = n_points)

S_max_parcel_Nₗ = zeros(n_points)
S_max_ARGmod_Nₗ = zeros(n_points)
S_max_parcel_Nᵢ = zeros(n_points)
S_max_ARGmod_Nᵢ = zeros(n_points)
S_max_parcel_rₗ = zeros(n_points)
S_max_ARGmod_rₗ = zeros(n_points)
S_max_parcel_rᵢ = zeros(n_points)
S_max_ARGmod_rᵢ = zeros(n_points)
for i in 1:n_points
    S_max_parcel_Nₗ[i], S_max_ARGmod_Nₗ[i] = run_parcel_model.(Nₐ, Nₗ_range[i], FT(0), rₗ, rᵢ, w, FT)
    S_max_parcel_Nᵢ[i], S_max_ARGmod_Nᵢ[i] = run_parcel_model.(Nₐ, FT(0), Nᵢ_range[i], rₗ, rᵢ, w, FT)
    S_max_parcel_rₗ[i], S_max_ARGmod_rₗ[i] = run_parcel_model.(Nₐ, Nₗ, FT(0), rₗ_range[i], rᵢ, w, FT)
    S_max_parcel_rᵢ[i], S_max_ARGmod_rᵢ[i] = run_parcel_model.(Nₐ, FT(0), Nᵢ, rₗ, rᵢ_range[i], w, FT)
end

# Plot results
fig = MK.Figure(size = (800, 600), fontsize = 20)
ax1 = MK.Axis(fig[1, 1], ylabel = "S_max / S_max_ARG", xlabel = "Nₗ [cm⁻³]", title = "No ice particle, rₗ=20 μm")
ax2 = MK.Axis(fig[1, 2], ylabel = "S_max / S_max_ARG", xlabel = "Nᵢ [cm⁻³]", title = "No liquid droplets, rᵢ=20 μm")
ax3 = MK.Axis(fig[2, 1], ylabel = "S_max / S_max_ARG", xlabel = "rₗ [μm]", title = "No ice particle, Nₗ=100 cm⁻³")
ax4 = MK.Axis(fig[2, 2], ylabel = "S_max / S_max_ARG", xlabel = "rᵢ [μm]", title = "No liquid droplets, Nᵢ=1 cm⁻³")

MK.lines!(ax1, Nₗ_range * 1e-6, S_max_parcel_Nₗ, label = "Parcel", linewidth = 2, color = :blue)
MK.lines!(ax1, Nₗ_range * 1e-6, S_max_ARGmod_Nₗ, label = "Modified ARG", linewidth = 2, color = :red)
MK.lines!(ax2, Nᵢ_range * 1e-6, S_max_parcel_Nᵢ, label = "Parcel", linewidth = 2, color = :blue)
MK.lines!(ax2, Nᵢ_range * 1e-6, S_max_ARGmod_Nᵢ, label = "Modified ARG", linewidth = 2, color = :red)
MK.lines!(ax3, rₗ_range * 1e6, S_max_parcel_rₗ, label = "Parcel", linewidth = 2, color = :blue)
MK.lines!(ax3, rₗ_range * 1e6, S_max_ARGmod_rₗ, label = "Modified ARG", linewidth = 2, color = :red)
MK.lines!(ax4, rᵢ_range * 1e6, S_max_parcel_rᵢ, label = "Parcel", linewidth = 2, color = :blue)
MK.lines!(ax4, rᵢ_range * 1e6, S_max_ARGmod_rᵢ, label = "Modified ARG", linewidth = 2, color = :red)

MK.ylims!(ax1, -0.05, 1.1)
MK.ylims!(ax3, -0.05, 1.1)
MK.axislegend(ax1, framevisible = false, labelsize = 16, position = :rc)
MK.axislegend(ax2, framevisible = false, labelsize = 16, position = :rc)
MK.axislegend(ax3, framevisible = false, labelsize = 16, position = :rc)
MK.axislegend(ax4, framevisible = false, labelsize = 16, position = :rc)

MK.save("parcel_vs_modifiedARG_aerosol_activation.svg", fig)
