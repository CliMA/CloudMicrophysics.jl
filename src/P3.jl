"""
Predicted particle properties scheme P3 for ice [MorrisonMilbrandt2015](@cite).

See online docs for more information.
"""
module P3Scheme

using DocStringExtensions

import SpecialFunctions as SF
import QuadGK as QGK
import RootSolvers as RS
import HCubature as HC
import LogExpFunctions

import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics.HetIceNucleation as CM_HetIce
import CloudMicrophysics.Microphysics2M as CM2

import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

const PSP3 = CMP.ParametersP3

include("P3_particle_properties.jl")
include("P3_size_distribution.jl")
include("P3_integral_properties.jl")
include("P3_terminal_velocity.jl")
include("P3_processes.jl")

end
