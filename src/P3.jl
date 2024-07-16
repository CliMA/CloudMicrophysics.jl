"""
Predicted particle properties scheme P3 for ice.
Implementation of Morrison and Milbrandt 2015 doi: 10.1175/JAS-D-14-0065.1

Note: Particle size is defined as its maximum length (i.e. max dimesion).
"""
module P3Scheme

import SpecialFunctions as SF
import QuadGK as QGK
import RootSolvers as RS
import HCubature as HC

import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics as CM
import CloudMicrophysics.TerminalVelocity as TV
import CloudMicrophysics.Microphysics2M as CM2

const PSP3 = CMP.ParametersP3

export thresholds, distribution_parameter_solver

include("P3_helpers.jl")
include("P3_particle_properties.jl")
include("P3_size_distribution.jl")
include("P3_terminal_velocity.jl")
include("P3_processes.jl")

end