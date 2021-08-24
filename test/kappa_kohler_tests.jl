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
T = 298.15       # air temperature
p = 100000.0   # air pressure
w = .05        # vertical velocity

# Accumulation mode
r_dry_accum = 0.243 * 1e-6 # μm
stdev_accum = 1.4          # -
N_accum = 100.0 *10^6      # 1/m3

# Sea Salt - universal parameters
M_seasalt = 0.058443
ρ_seasalt = 2170.0
ϕ_seasalt = 0.9
ν_seasalt = 2.0
ϵ_seasalt = 1.0

# Sea Salt - Kappa-Kohler value 
kappa_seasalt = 1.12

# Ammonium Sulfate - universal parameters
# Ammonum Sulfate (AS) - universal parameters
M_as = .12314      # molar mass
ρ_as = 1130.       # density
ϕ_as = 0.9         # osmotic coeff
ν_as = 3.0         # n ions 

# Ammonium Sulfate - Kappa-Kohler value 
κ_as = .51         


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
    (kappa_seasalt,),
    1,
)

AM_2 = KappaKohlerAerosolDistribution((kappa_accum_seasalt,))
println(AM_2.KK_Modes)

accum_ammonium = Mode(
    r_dry_accum,
    stdev_accum,
    N_accum,
    (1.0,),
    (1.0,),
    (ϕ_as,),
    (M_as,),
    (ν_as,),
    (ρ_as,),
    1
)

AM_3 = AerosolDistribution((accum_ammonium,))

kappa_accum_ammonium = KappaKohlerMode(
    r_dry_accum,
    stdev_accum, 
    N_accum,
    (1.0,),
    (M_as,),
    (κ_as,),
    1,
)

AM_4 = KappaKohlerAerosolDistribution((kappa_accum_ammonium,))

accum_seasalt_ammonium = Mode(
    r_dry_accum,
    stdev_accum,
    N_accum,
    (.75, .25),
    (ϵ_seasalt, ϵ_seasalt),
    (ϕ_seasalt,ϕ_as),
    (M_seasalt, M_as),
    (ν_seasalt, ν_as),
    (ρ_seasalt, ρ_as),
    2    
)
AM_5 = AerosolDistribution((accum_seasalt_ammonium,))

kappa_accum_seasalt_ammonium = KappaKohlerMode(
    r_dry_accum,
    stdev_accum, 
    N_accum,
    ((.75 / ρ_seasalt) /  (.75 / ρ_seasalt + .25 / ρ_as), 1-(.75 / ρ_seasalt) /  (.75 / ρ_seasalt + .25 / ρ_as) ),
    (M_seasalt, M_as),
    (kappa_seasalt,κ_as ),
    2,
)

AM_6 = KappaKohlerAerosolDistribution((kappa_accum_seasalt_ammonium,))


# println(total_N_activated(param_set, AM_1, T, p, w))
# println(KK_total_N_activated(param_set, AM_2, T, p, w))
# println(KK_critical_supersaturation(param_set, AM_2, T))
# println(critical_supersaturation(param_set, AM_1, T))
ws = range(0.001, stop = 1, length = 1000)
Ts = range(230, stop = 320, length = 1000)

frac_activated_nom = zeros(length(ws))
frac_activated_kk = zeros(length(ws))
# for i in 1:length(ws)
#     total_N_act_nom = total_N_activated(param_set, AM_5, T, p, ws[i])
#     total_N_act_kk = KK_total_N_activated(param_set, AM_6, T, p, ws[i])
#     frac_activated_nom[i] = total_N_act_nom / N_accum
#     frac_activated_kk[i] = total_N_act_kk / N_accum
# end

for i in 1:length(Ts)
    total_N_act_nom = total_N_activated(param_set, AM_1, Ts[i], p, w)
    total_N_act_kk = KK_total_N_activated(param_set, AM_2, Ts[i], p, w)
    frac_activated_nom[i] = total_N_act_nom / N_accum
    frac_activated_kk[i] = total_N_act_kk / N_accum
end

println("ws")
println(Ts)
println("Nominal")
println(frac_activated_nom)
println("Kappa")
println(frac_activated_kk)

# display(plot(Ts, frac_activated_nom, label="Kelvin-Kohler"))
# display(plot!(Ts, frac_activated_kk, label="Kappa-Kohler"))
# savefig("test.png")