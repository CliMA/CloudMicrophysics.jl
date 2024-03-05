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

include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float64

# Check for data file in AIDA_data folder
fpath = joinpath(
    pkgdir(CM),
    "calibration",
    "parcel_ice_nucleation",
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