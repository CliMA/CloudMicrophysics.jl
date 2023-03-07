using Test
using CLIMAParameters

using CloudMicrophysics
using CairoMakie

const TT = Test
const AM = CloudMicrophysics.AerosolModel
const CG = CloudMicrophysics.Coagulation


FT = Float64
toml_dict = CLIMAParameters.create_toml_dict(FT)

param_names = ["MSLP", "T_surf_ref", "k_Boltzmann"]
params = CLIMAParameters.get_parameter_values!(toml_dict, param_names)
params = (; params...)

# Adapted from microphysics_tests.jl

# Stdev, radius, number concentration from Whitby 1978 - Table 2.
# Aitken mode
r_dry_paper = 0.013  # um
stdev_paper = 1.7          # -
N_1_paper = 7.7e3 # no / cm3
# Accumulation mode
r_dry_accum = 0.069  # um
stdev_accum = 2.03          # -
N_accum = 1.3e3


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

air_pressure = params.MSLP        #  standard surface pressure (pa)
air_temp = params.T_surf_ref  #  standard surface temperature (K)

# TODO: Add better values
particle_density_acc = 0.1
particle_density_ait = 0.1

ad = AM.AerosolDistribution((aitken_sulfate_κ, accum_sulfate_κ))



function plot(ad::AM.AerosolDistribution)
    step=1e-3
    stop=2
    f = Figure()
    ax = Axis(f[1, 1],
        title = "",
        xlabel = "Particle Diameter",
        ylabel = "Number"
    )
    sizes = range(eps(), stop=stop, step=step)
    for mode in ad.Modes
        lognormal = CG.lognormal_dist(mode)
        distribution = lognormal.(sizes)
        # Normalize
        distribution *= mode.N / sum(distribution)
        lines!(sizes, distribution)
    end
    save("test.png", f)
    return f
end
plot(ad)

# Test that cubature returns 1 for integrating over the distribution
function test_cubature(am)
    dist = CG.lognormal_dist(am)
    integrand = dp -> dist(dp[1]) * dist(dp[2])
    result = CG.cubature(integrand) / am.N^2
    @test result ≈ 1 rtol=1e-9
end

test_cubature(ad.Modes[1])

function timestep_coagulation_changes(
    ad::AM.AerosolDistribution,
    particle_density_ait,
    particle_density_acc,
    air_pressure,
    air_temp,
    params,
    t
)
    for i in 1:t
        rates = CG.quadrature(
            ad,
            particle_density_ait,
            particle_density_acc,
            air_pressure,
            air_temp,
            params
        )
        (intra_rates, inter_rates) = rates
        for rate in rates
            for (i, mode) in enumerate(ad.Modes)
                modal_rates = rate[i]
                mode.N += modal_rates[1]
                mode.r_dry += modal_rates[2]/2
                # Need to add surface area, volume?
            end
        end
        if i%100 == 0
            println(rates)
            println(ad)
        end
    end
    return ad
end

function binkowski_roselle_rates(
    ad::AM.AerosolDistribution,
    particle_density_ait,
    particle_density_acc,
    air_pressure,
    air_temp,
    params)

    ad = timestep_coagulation_changes(
        ad::AM.AerosolDistribution,
        particle_density_ait,
        particle_density_acc,
        air_pressure,
        air_temp,
        params,
        3600
    )

end

binkowski_roselle_rates(
    ad,
    particle_density_ait,
    particle_density_acc,
    air_pressure,
    air_temp,
    params
)