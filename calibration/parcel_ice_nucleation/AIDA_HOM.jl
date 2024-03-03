import EDF

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

include(joinpath(
    pkgdir(CM),
    "calibration",
    "parcel_ice_nucleation",
    "AIDA Chamber Data/homogeneous freezing/IN05_17/in05_17_aida.edf"
))

