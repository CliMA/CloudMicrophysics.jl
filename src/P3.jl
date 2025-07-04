"""
Predicted particle properties scheme P3 for ice [MorrisonMilbrandt2015](@cite).

See online docs for more information.
"""
module P3Scheme

using DocStringExtensions

import SpecialFunctions as SF
import QuadGK as QGK
import RootSolvers as RS
import LogExpFunctions
import StaticArrays as SA

import ClimaParams as CP

import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics.DistributionTools as DT
import CloudMicrophysics.HetIceNucleation as CM_HetIce
import CloudMicrophysics.Microphysics2M as CM2

include("P3_particle_properties.jl")
include("P3_size_distribution.jl")
include("P3_integral_properties.jl")
include("P3_terminal_velocity.jl")
include("P3_processes.jl")

end
