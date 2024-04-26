import DelimitedFiles

import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions
import Random
import Distributions
import CairoMakie as MK
import LinearAlgebra
import OrdinaryDiffEq as ODE

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import Thermodynamics as TD

FT = Float64
include(
    joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "calibration.jl"),
)

# Check for data file in AIDA_data folder
fpath = joinpath(
    pkgdir(CM),
    "papers",
    "ice_nucleation_2024",
    "AIDA_data",
)
data_file = "in05_17_aida.edf"

if isfile(joinpath(fpath, data_file))
    data_path = joinpath(fpath, data_file)
else
    # if data does not exist, download from dropbox
    !isfile(joinpath(fpath, data_file))
    url = "https://www.dropbox.com/scl/fi/zm38lvhdu4sm8um2fvz9y/in05_17_aida.edf?rlkey=156zoi14y4xix0e569ja0ftzv&dl=0"
    download(url, data_file)
    mkdir(joinpath(fpath, data_file))
    mv(
        joinpath(pkgdir(CM), data_file),
        joinpath(fpath, data_file),
        force = true,
    )
    data_path = joinpath(fpath, data_file)
end

IN0517 = DelimitedFiles.readdlm(data_path, skipstart = 125)
IN0517_ICNC = IN0517[100:end, 6] # Nᵢ [cm^-3]

# Plotting AIDA data
# IN05_17
IN0517_fig = MK.Figure(size = (800, 600))
ax1 = MK.Axis(IN0517_fig[1, 1], ylabel = "ICNC [cm^-3]", xlabel = "time [s]", title = "AIDA IN05_17")
MK.lines!(ax1, IN0517[100:end, 1], IN0517[100:end, 6], label = "AIDA")
#MK.axislegend(ax1, framevisible = true, labelsize = 12, position = :rc)
#IN0517_fig

# Define parameters
# IN05_17

# Define initial conditions
# IN05_17

# Calibration
# IN05_17
IN0517_output = calibrate_J_parameters(
    FT,
    "ABHOM",
    IN0517_params,
    IN0517_IC,
    IN0517_ICNC,
    IN0517_Γ,
)
IN0517_calibrated_parameters = [IN0517_output[1], IN0517_output[2]]
IN0517_ensemble_means = ensemble_means(
    IN0517_output[3],
    size(IN0517_output[3])[1],
    size(IN0517_output[3][1])[2],
)

# Plots
# IN05_17
IN0517_calibrated_fig = MK.Figure(size = (800, 600))
ax3 = MK.Axis(IN0517_calibrated_fig[1, 1], ylabel = "m coefficient [-]", xlabel = "iteration number", title = "IN05_17 (HOM)")
ax4 = MK.Axis(IN0517_calibrated_fig[2, 1], ylabel = "c coefficient [-]", xlabel = "iteration number")
MK.lines!(ax3, iterations, IN0517_ensemble_means[1], label = "ensemble mean", color = :orange)
MK.lines!(ax4, iterations, IN0517_ensemble_means[2], label = "ensemble mean", color = :orange)
