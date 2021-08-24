using Test
# import Pkg
# Pkg.add("Plots")
using Plots
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
T = [282.8054379700115, 282.82011913510763, 282.8182541470769, 282.816547286481, 
     282.8152302772072, 282.83225081516434, 282.83107153716537, 282.8300594435075, 
     282.82918785057166, 282.8284346224201, 282.8277816799783, 282.8272141498667, 
     282.8267193163115, 282.8262868328548, 282.8259080527273, 282.825575586196, 
     282.8252831734087, 282.8250255282393, 282.82479813646944, 282.8245971038684]     # air temperature
p = [99352.28143153602, 99381.95622299137, 99381.94656495028, 99381.85770219186, 
     99381.90580284428, 99411.64002911752, 99411.58811377818, 99411.54174606412, 
     99411.50019030878, 99411.46280238577, 99411.42915075587, 99411.39885950244, 
     99411.37147577599, 99411.34671755912, 99411.32435993387, 99411.30416091887, 
     99411.28589558789, 99411.26938254002, 99411.25445547458, 99411.24096135839] # air pressure
w = .5        # vertical velocity

r_factors = [0.26666667, 0.29333333, 0.32, 0.34666667, 0.37333333, 0.4, 
             0.42666667, 0.45333333, 0.48, 0.50666667, 0.53333333, 0.56,
             0.58666667, 0.61333333, 0.64, 0.66666667, 0.69333333, 0.72,
             0.74666667, 0.77333333]

r_dry = 0.04 * 1e-6
s1 = 1.4
s2 = 1.6
kapp = .53
N = 1000.
M_as = .12314  

mode1 = KappaKohlerMode(
    r_dry,
    s1, 
    N,
    (1.0,),
    (M_as,),
    (kapp,),
    1
)

N_frac = zeros(length(r_factors))
for i in 1:length(r_factors)
    mode2 = KappaKohlerMode(
    r_dry/r_factors[i],
    s2, 
    N,
    (1.0,),
    (M_as,),
    (kapp,),
    1
    )
    model = KappaKohlerAerosolDistribution((mode1, mode2))
    N_total = model.KK_Modes[1].N + model.KK_Modes[2].N 
    N_act = KK_total_N_activated(param_set, model, T[i], p[i], w)
    N_frac[i] = N_act/N_total
end
print(N_frac)