import Test
import CLIMAParameters

import CloudMicrophysics

const TT = Test
const AM = CloudMicrophysics.AerosolModel
const CG = CloudMicrophysics.Coagulation

include("../src/Coagulation.jl")
using .Coagulation
local_exp_file = joinpath(@__DIR__,  "temp_params.toml")
FT = Float64
toml_dict =
    CLIMAParameters.create_toml_dict(FT; override_file = local_exp_file)

param_names = ["MSLP", "T_surf_ref", "k_Boltzmann"]
params = CLIMAParameters.get_parameter_values!(toml_dict, param_names)
params = (; params...)

# Adapted from microphysics_tests.jl

# Accumulation mode
r_dry_accum = 0.243 * 1e-6 # m
stdev_accum = 1.8          # -
N_accum = 100.0 * 1e6      # 1/m3
# Aitken mode
r_dry_paper = 0.05 * 1e-6  # m
stdev_paper = 1.6          # -
N_1_paper = 100.0 * 1e6    # 1/m3

# TODO - move areosol properties to CLIMAParameters
# Sea Salt - universal parameters
M_seasalt = 0.058443
ρ_seasalt = 2170.0
ϕ_seasalt = 0.9
ν_seasalt = 2.0
ϵ_seasalt = 1.0
κ_seasalt = 1.12   # TODO - which value to take?
# Sulfate - universal parameters
M_sulfate = 0.132
ρ_sulfate = 1770.0
ϕ_sulfate = 1.0
ν_sulfate = 3.0
ϵ_sulfate = 1.0
κ_sulfate = 0.53   # TODO - which value to take?

# Aerosol size distribution modes

aitken_sulfate_κ = AM.Mode_κ(
    r_dry_paper,
    stdev_paper,
    N_1_paper,
    (1.0,),
    (1.0,),
    (M_sulfate,),
    (κ_sulfate,),
    1,
)

accum_sulfate_κ = AM.Mode_κ(
    r_dry_accum,
    stdev_accum,
    N_accum,
    (1.0,),
    (1.0,),
    (M_sulfate,),
    (κ_sulfate,),
    1,
)

ad = AM.AerosolDistribution((aitken_sulfate_κ,accum_sulfate_κ))

air_pressure = params.MSLP        #  standard surface pressure (pa)
air_temp = params.T_surf_ref  #  standard surface temperature (K)

# TODO: Add better values
particle_density_acc = 0.1
particle_density_ait = 0.1
gas_viscosity = 1e-5

Coagulation.coagulation_quadrature(ad, particle_density_ait, particle_density_acc, gas_viscosity, air_pressure, air_temp, params)
