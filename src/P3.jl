"""
Predicted particle properties scheme P3 for ice [MorrisonMilbrandt2015](@cite).

See online docs for more information.
"""
module P3Scheme

using DocStringExtensions

import SpecialFunctions as SF
import RootSolvers as RS
import LogExpFunctions
import StaticArrays as SA
import UnrolledUtilities as UU

import ClimaParams as CP

import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics.DistributionTools as DT
import CloudMicrophysics.HetIceNucleation as CM_HetIce
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Utilities as UT
import CloudMicrophysics: ShowMethods

# Quadrature rules live in their own module (included before `Parameters` so a
# constructed rule can be stored on a parameter struct). Bring the names into
# `P3Scheme` so the `integrate(...; quad = ...)` API and the `P3.GaussLegendre`
# / `P3.ChebyshevGauss` references used in the scheme and its tests are unchanged.
import CloudMicrophysics.Quadrature:
    QuadratureRule, ChebyshevGauss, GaussLegendre, integrate, node, weight,
    inv_weight_fun, subintervals

include("P3_particle_properties.jl")
include("P3_size_distribution.jl")
include("P3_integral_properties.jl")
include("P3_terminal_velocity.jl")
include("P3_processes.jl")

end
