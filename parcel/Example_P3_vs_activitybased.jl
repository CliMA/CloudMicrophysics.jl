import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaParams as CP

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float32
# get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
aps = CMP.AirProperties(FT)
wps = CMP.WaterProperties(FT)
ip = CMP.IceNucleationParameters(FT)

# Constants
ρₗ = wps.ρw
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)

# Initial conditions
Nₐ = FT(2000)
Nₗ = FT(2000)
Nᵢ = FT(0)
rₗ = FT(1.25e-6)
p₀ = FT(20000)
T₀_dep = FT(238)
T₀_het = FT(239)
T₀_hom = FT(236.5)
qᵥ = FT(8.3e-4)
qₗ = FT(Nₗ * 4 / 3 * π * rₗ^3 * ρₗ / 1.2)
qᵢ = FT(0)
ln_INPC = FT(0)

# Moisture dependent initial conditions
q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
ts_dep = TD.PhaseNonEquil_pTq(tps, p₀, T₀_dep, q)
ts_het = TD.PhaseNonEquil_pTq(tps, p₀, T₀_het, q)
ts_hom = TD.PhaseNonEquil_pTq(tps, p₀, T₀_hom, q)
ρₐ_dep = TD.air_density(tps, ts_dep)
ρₐ_het = TD.air_density(tps, ts_het)
ρₐ_hom = TD.air_density(tps, ts_hom)
Rₐ = TD.gas_constant_air(tps, q)
eₛ_dep = TD.saturation_vapor_pressure(tps, T₀_dep, TD.Liquid())
eₛ_het = TD.saturation_vapor_pressure(tps, T₀_het, TD.Liquid())
eₛ_hom = TD.saturation_vapor_pressure(tps, T₀_hom, TD.Liquid())
e = eᵥ(qᵥ, p₀, Rₐ, R_v)

Sₗ_dep = FT(e / eₛ_dep)
Sₗ_het = FT(e / eₛ_het)
Sₗ_hom = FT(e / eₛ_hom)
IC_dep = [Sₗ_dep, p₀, T₀_dep, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]
IC_het = [Sₗ_het, p₀, T₀_het, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]
IC_hom = [Sₗ_hom, p₀, T₀_hom, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]

# Simulation parameters passed into ODE solver
r_nuc = FT(1.25e-6)                     # assumed size of nucleated particles
w = FT(0.5)                             # updraft speed
α_m = FT(0.5)                           # accomodation coefficient
const_dt = FT(0.1)                      # model timestep
t_max = FT(50)
aerosol_none = []
aerosol_deposition = [CMP.Feldspar(FT), CMP.Ferrihydrite(FT), CMP.Kaolinite(FT)]
aerosol_heterogeneous = [CMP.DesertDust(FT), CMP.Illite(FT), CMP.Kaolinite(FT)]
deposition_modes = ["P3_dep", "ABDINM"]
heterogeneous_modes = ["P3_het", "ABIFM"]
homogeneous_modes = ["P3_hom", "ABHOM"]
deposition_growth = "Deposition"
size_distribution = "Monodisperse"

# Plotting
fig = MK.Figure(resolution = (1000, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Ice Saturation [-]")
ax2 = MK.Axis(fig[1, 2], ylabel = "ICNC [cm^-3]", yscale = log10)
ax3 = MK.Axis(fig[1, 3], ylabel = "Temperature [K]")
ax4 = MK.Axis(fig[2, 1], ylabel = "Ice Saturation [-]")
ax5 = MK.Axis(fig[2, 2], ylabel = "ICNC [cm^-3]", yscale = log10)
ax6 = MK.Axis(fig[2, 3], ylabel = "Temperature [K]")
ax7 = MK.Axis(fig[3, 1], ylabel = "Ice Saturation [-]", xlabel = "Time [s]")
ax8 = MK.Axis(fig[3, 2], ylabel = "ICNC [cm^-3]", xlabel = "Time [s]")
ax9 = MK.Axis(fig[3, 3], ylabel = "Temperature [K]", xlabel = "Time [s]")

MK.ylims!(ax2, 1e-9, 1)
MK.ylims!(ax5, 1e-16, 1e-5)

# solve ODE
#! format: off
for deposition in deposition_modes
    if deposition == "P3_dep"
        local params = parcel_params{FT}(
            const_dt = const_dt,
            r_nuc = r_nuc,
            w = w,
            aerosol = aerosol_none,
            deposition = deposition,
            deposition_growth = deposition_growth,
            liq_size_distribution = size_distribution,
        )
        # solve ODE
        local sol = run_parcel(IC_dep, FT(0), t_max, params)
        # Plot
        MK.lines!(ax1, sol.t, S_i.(tps, sol[3, :], (sol[1, :])), label = deposition, color = :blue) # saturation
        MK.lines!(ax2, sol.t, sol[9, :] * 1e-6, color = :blue)                                      # ICNC
        MK.lines!(ax3, sol.t, sol[3, :], color = :blue)                                             # Temperature
    elseif deposition == "ABDINM"
        for aerosol in aerosol_deposition
            local params = parcel_params{FT}(
                const_dt = const_dt,
                r_nuc = r_nuc,
                w = w,
                aerosol = aerosol,
                deposition = deposition,
                deposition_growth = deposition_growth,
                liq_size_distribution = size_distribution,
            )
            # solve ODE
            local sol = run_parcel(IC_dep, FT(0), t_max, params)
            # Plot
            if aerosol == CMP.Feldspar(FT)
                line_color = :green
            elseif aerosol == CMP.Ferrihydrite(FT)
                line_color = :orange
            elseif aerosol == CMP.Kaolinite(FT)
                line_color = :magenta
            end
            MK.lines!(ax1, sol.t, S_i.(tps, sol[3, :], (sol[1, :])), label = deposition, color = line_color) # saturation
            MK.lines!(ax2, sol.t, sol[9, :] * 1e-6, color = line_color)                                      # ICNC
            MK.lines!(ax3, sol.t, sol[3, :], color = line_color)                                             # Temperature
        end
    end
end

for heterogeneous in heterogeneous_modes
    if heterogeneous == "P3_het"
        local params = parcel_params{FT}(
            const_dt = const_dt,
            r_nuc = r_nuc,
            w = w,
            aerosol = aerosol_none,
            heterogeneous = heterogeneous,
            deposition_growth = deposition_growth,
            liq_size_distribution = size_distribution,
        )
        # solve ODE
        local sol = run_parcel(IC_het, FT(0), t_max, params)
        # Plot
        MK.lines!(ax4, sol.t, S_i.(tps, sol[3, :], (sol[1, :])), label = heterogeneous, color = :blue)  # saturation
        MK.lines!(ax5, sol.t, sol[9, :] .* 1e-6, color = :blue)                                         # ICNC
        MK.lines!(ax6, sol.t, sol[3, :], color = :blue)                                                 # Temperature
    elseif heterogeneous == "ABIFM"
        for aerosol in aerosol_heterogeneous
            local params = parcel_params{FT}(
                const_dt = const_dt,
                r_nuc = r_nuc,
                w = w,
                aerosol = aerosol,
                heterogeneous = heterogeneous,
                deposition_growth = deposition_growth,
                liq_size_distribution = size_distribution,
            )
            # solve ODE
            local sol = run_parcel(IC_het, FT(0), t_max, params)
            # Plot
            if aerosol == CMP.DesertDust(FT)
                line_color = :green
                MK.lines!(ax4, sol.t, S_i.(tps, sol[3, :], (sol[1, :])), label = heterogeneous, color = line_color) # saturation
                MK.lines!(ax5, sol.t, sol[9, :] .* 1e-6, color = line_color)                                        # ICNC
                MK.lines!(ax6, sol.t, sol[3, :], color = line_color)                                                # Temperature
            elseif aerosol == CMP.Illite(FT)
                line_color = :orange
                MK.lines!(ax4, sol.t, S_i.(tps, sol[3, :], (sol[1, :])), label = heterogeneous, color = line_color) # saturation
                MK.lines!(ax5, sol.t, sol[9, :] .* 1e-6, color = line_color)                                        # ICNC
                MK.lines!(ax6, sol.t, sol[3, :], color = line_color)                                                # Temperature
            elseif aerosol == CMP.Kaolinite(FT)
                line_color = :magenta
                MK.lines!(ax4, sol.t, S_i.(tps, sol[3, :], (sol[1, :])), label = heterogeneous, color = line_color) # saturation
                MK.lines!(ax5, sol.t, sol[9, :] .* 1e-6, color = line_color)                                        # ICNC
                MK.lines!(ax6, sol.t, sol[3, :], color = line_color)                                                # Temperature
            end
        end
    end
end

for homogeneous in homogeneous_modes
    local params = parcel_params{FT}(
        const_dt = const_dt,
        r_nuc = r_nuc,
        w = w,
        aerosol = aerosol_none,
        homogeneous = homogeneous,
        deposition_growth = deposition_growth,
        liq_size_distribution = size_distribution,
    )
    # solve ODE
    local sol = run_parcel(IC_hom, FT(0), t_max, params)
    if homogeneous == "P3_hom"
        line_color = :blue
    elseif homogeneous == "ABHOM"
        line_color = :orange
    end
    MK.lines!(ax7, sol.t, S_i.(tps, sol[3, :], (sol[1, :])), label = homogeneous, color = line_color)# saturation
    MK.lines!(ax8, sol.t, sol[9, :] .* 1e-6, color = line_color)                                     # ICNC
    MK.lines!(ax9, sol.t, sol[3, :], color = line_color)                                             # Temperature
end

MK.axislegend(ax1, framevisible = false, labelsize = 10, nbanks = 3, orientation = :horizontal, position = :lt)
MK.axislegend(ax4, framevisible = false, labelsize = 10, nbanks = 3, orientation = :horizontal, position = :lt)
MK.axislegend(ax7, framevisible = false, labelsize = 10, nbanks = 2, orientation = :horizontal, position = :lt)

#! format: on

MK.save("P3_vs_activitybased.svg", fig)
