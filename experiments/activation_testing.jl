import Plots

import CloudMicrophysics
import CLIMAParameters
import Thermodynamics

const PL = Plots
const AM = CloudMicrophysics.AerosolModel
const AA = CloudMicrophysics.AerosolActivation
const CMP = CloudMicrophysics.Parameters
const CP =  CLIMAParameters
const TD = Thermodynamics

include(joinpath(pkgdir(CloudMicrophysics), "test", "create_parameters.jl"))
FT = Float64
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
const param_set = cloud_microphysics_parameters(toml_dict)
thermo_params = CMP.thermodynamics_params(param_set)

# Atmospheric conditions
T = 294.0         # air temperature    # TODO - change 290, 295, 300, 305
p = 1000.0 *1e2   # air pressure       # TODO - change 980, 1000, 1100
w = 0.5           # vertical velocity  # TODO - 0.1, 0.5, 1.0, 2.5

# We need the phase partition here only so that we can compute the
# moist R_m and cp_m in aerosol activation module.
# We are assuming here saturated conditions and no liquid water or ice.
# This is consistent with the assumptions of the aerosol activation scheme.
p_vs = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
q_vs = 1 / (1 - CMP.molmass_ratio(param_set) * (p_vs - p) / p_vs)
q = TD.PhasePartition(q_vs, 0.0, 0.0)

# Aitken mode
r_atk = 0.05 * 1e-6 # um   #TODO: test for  0.01 - 10 um
s_atk = 1.6         # -
N_atk = 100.0 * 1e6 # 1/m3

# Accumulation mode
r_acc = 0.2 * 1e-6  # um   #TODO: test for
s_acc = 1.5         # -
N_acc = 100.0 * 1e6 # 1/m3

# Coarse mode
r_crs = 0.6 * 1e-6  # um   #TODO: test for
s_crs = 1.8         # -
N_crs = 10.0 * 1e6  # 1/m3

# Sulfate - universal parameters
M_sulfate = 0.132
ρ_sulfate = 1770.0
κ_sulfate = 0.7 # TODO - check

# Sea salt - universal parameters
M_seasalt = 0.5844
ρ_seasalt = 2160.0 # TODO - check
κ_seasalt = 1.2 # TODO - check

# For now assuming one spiece per mode
n_components = 1
vol_mix_ratio = (1.0,)
mass_mix_ratio = (1.0,)

mode_atk = AM.Mode_κ(
    r_atk,
    s_atk,
    N_atk,
    vol_mix_ratio,
    mass_mix_ratio,
    (M_sulfate,),
    (κ_sulfate,),
    n_components,
)
mode_acc = AM.Mode_κ(
    r_acc,
    s_acc,
    N_acc,
    vol_mix_ratio,
    mass_mix_ratio,
    (M_sulfate,),
    (κ_sulfate,),
    n_components,
)
mode_crs = AM.Mode_κ(
    r_crs,
    s_crs,
    N_crs,
    vol_mix_ratio,
    mass_mix_ratio,
    (M_sulfate,),
    (κ_sulfate,),
    n_components,
)

AD =  AM.AerosolDistribution((mode_atk, mode_acc, mode_crs))
N_act = AA.total_N_activated(param_set, AD, T, p, w, q)
S_max = AA.max_supersaturation(param_set, AD, T, p, w, q)

print("N_act [1/cm3] = ", N_act * 1e-6, ", fraction [%] = ", N_act / (N_atk + N_acc + N_crs) * 100, ", S_max [%] = ", S_max * 100 )
