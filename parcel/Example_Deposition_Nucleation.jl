import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float32

# get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)

# Initial conditions
Nₐ = FT(2000 * 1e3)
Nₗ = FT(0)
Nᵢ = FT(0)
p₀ = FT(20000)
T₀ = FT(230)
qᵥ = FT(3.3e-4)
qₗ = FT(0)
qᵢ = FT(0)
x_sulph = FT(0)

# Moisture dependent initial conditions
q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
ts = TD.PhaseNonEquil_pTq(tps, p₀, T₀, q)
ρₐ = TD.air_density(tps, ts)
Rₐ = TD.gas_constant_air(tps, q)
Rᵥ = TD.Parameters.R_v(tps)
eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
e = eᵥ(qᵥ, p₀, Rₐ, Rᵥ)

Sₗ = FT(e / eₛ)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]

# Simulation parameters passed into ODE solver
r_nuc = FT(1.25e-6)                     # assumed size of nucleated particles
w = FT(3.5 * 1e-2)                      # updraft speed
const_dt = FT(0.1)                      # model timestep
t_max = FT(100)
aerosol_1 = [CMP.DesertDust(FT), CMP.ArizonaTestDust(FT)]               # aersol type for DustDeposition
aerosol_2 = [CMP.Feldspar(FT), CMP.Ferrihydrite(FT), CMP.Kaolinite(FT)] # aersol type for DepositionNucleation
deposition_modes = ["MohlerRate", "ABDINM"]
deposition_growth = "Deposition"
size_distribution = "Monodisperse"

# Plotting
fig = MK.Figure(size = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Ice Saturation [-]")
ax2 = MK.Axis(fig[1, 2], ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "q_ice [g/kg]", xlabel = "time")
ax4 = MK.Axis(fig[2, 2], ylabel = "N_ice [m^-3]", xlabel = "time")

for deposition in deposition_modes
    if deposition == "MohlerRate"
        for aerosol in aerosol_1
            local params = parcel_params{FT}(
                const_dt = const_dt,
                r_nuc = r_nuc,
                w = w,
                aerosol = aerosol,
                deposition = deposition,
                deposition_growth = deposition_growth,
                size_distribution = size_distribution,
            )
            # solve ODE
            local sol = run_parcel(IC, FT(0), t_max, params)

            # Plot results
            if aerosol == CMP.DesertDust(FT)
                aero_label = "DesertDust"
            elseif aerosol == CMP.ArizonaTestDust(FT)
                aero_label = "ArizonaTestDust"
            end
            MK.lines!(                               # saturation
                ax1,
                sol.t,
                S_i.(tps, sol[3, :], sol[1, :]),
                label = aero_label,
            )
            MK.lines!(ax2, sol.t, sol[3, :])         # temperature
            MK.lines!(ax3, sol.t, sol[6, :] * 1e3)   # q_ice
            MK.lines!(ax4, sol.t, sol[9, :])         # N_ice
        end

    elseif deposition == "ABDINM"
        for aerosol in aerosol_2
            local params = parcel_params{FT}(
                const_dt = const_dt,
                r_nuc = r_nuc,
                w = w,
                aerosol = aerosol,
                deposition = deposition,
                deposition_growth = deposition_growth,
                size_distribution = size_distribution,
            )
            # solve ODE
            local sol = run_parcel(IC, FT(0), t_max, params)

            # Plot results
            if aerosol == CMP.Feldspar(FT)
                aero_label = "Feldspar"
            elseif aerosol == CMP.Ferrihydrite(FT)
                aero_label = "Ferrihydrite"
            elseif aerosol == CMP.Kaolinite(FT)
                aero_label = "Kaolinite"
            end
            MK.lines!(                                                      # saturation
                ax1,
                sol.t,
                S_i.(tps, sol[3, :], sol[1, :]),
                label = aero_label,
                linestyle = :dash,
            )
            MK.lines!(ax2, sol.t, sol[3, :], linestyle = :dash)             # temperature
            MK.lines!(ax3, sol.t, sol[6, :] * 1e3, linestyle = :dash)       # q_ice
            MK.lines!(ax4, sol.t, sol[9, :], linestyle = :dash)             # N_ice
        end
    end
end

MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 12,
    orientation = :horizontal,
    nbanks = 3,
    position = :lt,
)

MK.save("deposition_nucleation.svg", fig)
