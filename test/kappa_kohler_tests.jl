using Test

using Thermodynamics

using CloudMicrophysics.AerosolModel
include("//home/idularaz/CloudMicrophysics.jl/src/KappaKohlerAerosolActivation.jl")
include("/home/idularaz/CloudMicrophysics.jl/src/KappaKohlerAerosolModel.jl")
using CloudMicrophysics.AerosolActivation

using CLIMAParameters
using CLIMAParameters: gas_constant
using CLIMAParameters.Planet:
    molmass_water, ρ_cloud_liq, grav, cp_d, surface_tension_coeff
using CLIMAParameters.Atmos.Microphysics

struct EarthParameterSet <: AbstractEarthParameterSet end
const EPS = EarthParameterSet
const param_set = EarthParameterSet()

# Atmospheric conditions
T = 298.15       # air temperature
p = 100000.0   # air pressure
w = 0.5        # vertical velocity

# Accumulation mode
r_dry_accum = 0.243 * 1e-6 # μm
stdev_accum = 1.4          # -
N_accum = 100.0 * 1e6      # 1/m3

# Sea Salt - universal parameters
M_seasalt = 0.058443
ρ_seasalt = 2170.0
ϕ_seasalt = 0.9
ν_seasalt = 2.0
ϵ_seasalt = 1.0

# Sea Salt - Kappa-Kohler value 
kappa_seasalt = 1.12

# Test case 1: Aerosol Model without Kappa parameter
accum_seasalt = Mode(
    r_dry_accum,
    stdev_accum,
    N_accum,
    (1.0,),
    (ϵ_seasalt,),
    (ϕ_seasalt,),
    (M_seasalt,),
    (ν_seasalt,),
    (ρ_seasalt,),
    1
)
AM_1 = AerosolDistribution((accum_seasalt,))

# Test case 2: Aerosol Model with Kappa parameter
kappa_accum_seasalt = KappaKohlerMode(
    r_dry_accum,
    stdev_accum, 
    N_accum,
    (1.0,),
    (M_seasalt,),
    kappa_seasalt,
    1,
)

AM_2 = KappaKohlerAerosolDistribution((kappa_accum_seasalt,))


println(total_N_activated(param_set, AM_1, T, p, w))
println(KK_total_N_activated(param_set, AM_2, T, p, w))
# println(KK_critical_supersaturation(param_set, AM_2, T))
println(critical_supersaturation(param_set, AM_1, T))
