"""
Flexible N-moment microphysics representation, including: 
  - Generalized collisional-coalescence described by a kernel K(x,y), with
  default parameterization based on Long collision kernel 
  - Power-series representation of fall-speed 
  - Power-series representation of condensational growth
  TODO: no representation of ventilation effects
"""
module CloudyExt

import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

import CloudMicrophysics.Common as CO
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsFlexible:
    CLSetup, coalescence, condensation, weighted_vt

import Cloudy as CL
import Cloudy.ParticleDistributions as CPD
import Cloudy.KernelFunctions as CLK

"""
A structure containing the subdistributions, their moments, and additional
dynamical parameters corresponding to rates of collision, sedimentation, and
condensation/evaporation
"""
function CLSetup{FT}(;
    pdists::Vector{<:CPD.PrimitiveParticleDistribution{FT}} = Vector([
        CPD.ExponentialPrimitiveParticleDistribution(
            FT(100 * 1e6),
            FT(1e5 * 1e-18 * 1e3),
        ), # 100/cm^3; 10^5 µm^3 = 1e-10 kg
        CPD.GammaPrimitiveParticleDistribution(
            FT(1 * 1e6),
            FT(1e6 * 1e-18 * 1e3),
            FT(1),
        ),   # 1/cm^3; 10^6 µm^3 = 1e-9 kg; k=1
    ]),
    mom::Vector{FT} = FT.([100.0 * 1e6, 1e-2, 1.0 * 1e6, 1e-3, 2e-12]),
    NProgMoms::Vector{Int} = [Integer(CPD.nparams(dist)) for dist in pdists],
    KernelFunc::CLK.KernelFunction{FT} = CLK.LongKernelFunction(
        FT(5.236e-10),  # 5.236e-10 kg;  
        FT(9.44e9),     # 9.44e9 m^3/kg^2/s;
        FT(5.78),       # 5.78 m^3/kg/s
    ),
    mass_thresholds::Vector{FT} = [FT(1e-9), FT(Inf)],
    kernel_order::Int = 1,
    kernel_limit::FT = FT(1 * 1e-9 * 1e3),   # ~1mm Diam ~ 1 mm^3 volume
    coal_data = nothing,
    vel::Vector{Tuple{FT, FT}} = [(FT(50.0), FT(1.0 / 6))], # 50 m/s/kg^(1/6),
    norms::Vector{FT} = [1e0, 1e0], # 1e6 / m^3; 1e-9 kg / m^3
) where {FT <: AbstractFloat}

    CLSetup{FT}(
        pdists,
        mom,
        NProgMoms,
        KernelFunc,
        mass_thresholds,
        kernel_order,
        kernel_limit,
        coal_data,
        vel,
        norms,
    )
end

"""
    coalescence(clinfo)

 - `clinfo` - kwarg structure containing pdists, moments, and coalescence parameters
 TODO: currently implemented only for analytical coalescence style

Returns a vector of moment tendencies due to collisional coalescence
"""
function coalescence(clinfo::CLSetup{FT}) where {FT}
    kernel_tensor = CL.KernelTensors.CoalescenceTensor(
        clinfo.KernelFunc,
        clinfo.kernel_order,
        clinfo.kernel_limit,
    )
    if isnothing(clinfo.coal_data)
        clinfo.coal_data = CL.Coalescence.initialize_coalescence_data(
            CL.EquationTypes.AnalyticalCoalStyle(),
            kernel_tensor,
            clinfo.NProgMoms,
            norms = clinfo.norms,
            dist_thresholds = clinfo.mass_thresholds,
        )
    end
    mom_norms =
        CL.get_moments_normalizing_factors(clinfo.NProgMoms, clinfo.norms)
    mom_normalized = clinfo.mom ./ mom_norms
    # first, update the particle distributions
    for (i, dist) in enumerate(clinfo.pdists)
        ind_rng = CL.get_dist_moments_ind_range(clinfo.NProgMoms, i)
        CPD.update_dist_from_moments!(dist, mom_normalized[ind_rng]) #clinfo.mom[ind_rng])
    end
    CL.Coalescence.update_coal_ints!(
        CL.EquationTypes.AnalyticalCoalStyle(),
        clinfo.pdists,
        clinfo.coal_data,
    )
    return clinfo.coal_data.coal_ints .* mom_norms
end

"""
    condensation(clinfo, aps, tps, q, ρ, T)

 - `clinfo` - kwarg structure containing pdists, moments, and coalescence parameters
 - `aps` - air properties
 - `tps` - thermodynamics parameters
 - `T` - air temperature
 - `S` - saturation ratio (supersaturation = S - 1)
Returns a vector of moment tendencies due to condensation/evaporation
"""
function condensation(
    clinfo::CLSetup{FT},
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    T::FT,
    S::FT,
) where {FT}
    ξ = CO.G_func(aps, tps, T, TD.Liquid())
    # first, update the particle distributions
    for (i, dist) in enumerate(clinfo.pdists)
        ind_rng = CL.get_dist_moments_ind_range(clinfo.NProgMoms, i)
        CPD.update_dist_from_moments!(dist, clinfo.mom[ind_rng])
    end
    return CL.Condensation.get_cond_evap(
        S - 1,
        (; ξ = ξ, pdists = clinfo.pdists),
    ) * (4 * pi / 3)^(2 / 3) / (1000.0^(1 / 3)) # TODO: correct this scaling in Cloudy.jl
end

"""
    weighted_vt(clinfo)

 - `clinfo` - kwarg structure containing pdists, moments, and coalescence parameters
Returns the integrated fall speeds corresponding to the rate of change of prognostic moments
"""
function weighted_vt(clinfo::CLSetup{FT}) where {FT}
    for (i, dist) in enumerate(clinfo.pdists)
        ind_rng = CL.get_dist_moments_ind_range(clinfo.NProgMoms, i)
        CPD.update_dist_from_moments!(dist, clinfo.mom[ind_rng])
    end
    sed_flux =
        CL.Sedimentation.get_sedimentation_flux(clinfo.pdists, clinfo.vel)
    return -sed_flux ./ clinfo.mom
end


end #module
