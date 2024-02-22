"""
Flexible N-moment microphysics representation, including: 
  - Generalized collisional-coalescence described by a kernel K(x,y), with
  default parameterization based on Long collision kernel 
  - Power-series representation of fall-speed 
  - Power-series representation of condensational growth
  TODO: no representation of ventilation effects
"""
module MicrophysicsFlexible

import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

import ..Common as CO
import ..Parameters as CMP

import Cloudy as CL
import Cloudy.ParticleDistributions as CPD
import Cloudy.KernelFunctions as CLK

export coalescence

"""
A structure containing the subdistributions, their moments, and additional
dynamical parameters corresponding to rates of collision, sedimentation, and
condensation/evaporation
"""
Base.@kwdef mutable struct CLSetup{FT} 
    "Subdistributions of the hydrometeor size distribution; defaults to
    exponential cloud mode and gamma rain mode"
    pdists::Vector{CPD.PrimitiveParticleDistribution{FT}} = Vector([
        CPD.ExponentialPrimitiveParticleDistribution(FT(0), FT(1.0)),
        CPD.GammaPrimitiveParticleDistribution(FT(0), FT(1.0), FT(1.0))
    ])
    "Moments of the subdistributions; should correspond with pdists"
    mom::Vector{FT} = [FT(0) for i in 1:5]
    "Total number of prognostic moments for each subdistribution"
    NProgMoms::Vector{Int} = [Integer(CPD.nparams(dist)) for dist in pdists]
    "Kernel function for collisional coalescence"
    KernelFunc::CLK.KernelFunction{FT} = CLK.LongKernelFunction(FT(0.5236), FT(9.44e-3), FT(5.78e-3))
    "Particle mass thresholds for analytical integration style"
    mass_thresholds::Vector{FT} = [10.0, Inf]
    "Polynomial order of the kernel approximation for coalescence"
    kernel_order::Int = 1
    "Upper limit for evaluation of the polynomial kernel approximation"
    kernel_limit::FT = FT(500)
    "Coalescence data"
    coal_data = nothing
    "Sedimentation rate parameters"
    vel::Vector{Tuple{FT,FT}} = [(FT(2.0), FT(1.0 / 6))]
end

"""
    coalescence(clinfo)

 - `clinfo` - kwarg structure containing pdists, moments, and coalescence parameters
 TODO: currently implemented only for analytical coalescence style

Returns a vector of moment tendencies due to collisional coalescence
"""
function coalescence(clinfo::CLSetup{FT}) where {FT}
    kernel_tensor = CL.KernelTensors.CoalescenceTensor(clinfo.KernelFunc, clinfo.kernel_order, clinfo.kernel_limit)
    if isnothing(clinfo.coal_data)
        clinfo.coal_data = CL.Coalescence.initialize_coalescence_data(
            CL.EquationTypes.AnalyticalCoalStyle(),
            kernel_tensor,
            clinfo.NProgMoms,
            dist_thresholds = clinfo.mass_thresholds
        )
    end
    # first, update the particle distributions
    for (i, dist) in enumerate(clinfo.pdists)
        ind_rng = CL.get_dist_moments_ind_range(clinfo.NProgMoms, i)
        CPD.update_dist_from_moments!(dist, clinfo.mom[ind_rng])
    end
    CL.Coalescence.update_coal_ints!(
        CL.EquationTypes.AnalyticalCoalStyle(),
        clinfo.pdists,
        clinfo.coal_data
    )
    return clinfo.coal_data.coal_ints
end

"""
    condensation(clinfo, aps, tps, q, ρ, T)

 - `clinfo` - kwarg structure containing pdists, moments, and coalescence parameters
 - `aps` - air properties
 - `tps` - thermodynamics parameters
 - `q` - phase partition
 - `ρ` - air density
 - `T` - air temperature
Returns a vector of moment tendencies due to condensation/evaporation
"""
function condensation(
    clinfo::CLSetup{FT},
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q::TD.PhasePartition{FT},
    ρ::FT,
    T::FT,
) where {FT}
    S = TD.supersaturation(tps, q, ρ, T, TD.Liquid())
    ξ = CO.G_func(aps, tps, T, TD.Liquid())

    # first, update the particle distributions
    for (i, dist) in enumerate(clinfo.pdists)
        ind_rng = CL.get_dist_moments_ind_range(clinfo.NProgMoms, i)
        CPD.update_dist_from_moments!(dist, clinfo.mom[ind_rng])
    end
    return CL.Condensation.get_cond_evap(S, (; ξ=ξ, pdists=clinfo.pdists))
end

"""
    sedimentation(clinfo)

 - `clinfo` - kwarg structure containing pdists, moments, and coalescence parameters
Returns the integrated fall speeds corresponding to the rate of change of prognostic moments
"""
function sedimentation(clinfo::CLSetup{FT}) where {FT}
    return CL.Sedimentation.get_sedimentation_flux((; pdists=clinfo.pdists, vel=clinfo.vel))
end


end #module