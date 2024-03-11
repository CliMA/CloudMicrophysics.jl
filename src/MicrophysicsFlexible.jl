"""
Flexible N-moment microphysics representation, including: 
  - Generalized collisional-coalescence described by a kernel K(x,y), with
  default parameterization based on Long collision kernel 
  - Power-series representation of fall-speed 
  - Power-series representation of condensational growth
  TODO: no representation of ventilation effects
  TODO: conversion back to N_rai, N_liq, q_rai, q_liq
"""
module MicrophysicsFlexible

"""
A structure containing the subdistributions, their moments, and additional
dynamical parameters corresponding to rates of collision, sedimentation, and
condensation/evaporation
"""
mutable struct CLSetup{FT}
    "Subdistributions of the hydrometeor size distribution; defaults to
    exponential cloud mode and gamma rain mode"
    pdists::Vector
    "Moments of the subdistributions; should correspond with pdists"
    mom::Vector{FT}
    "Total number of prognostic moments for each subdistribution"
    NProgMoms::Vector{Int}
    "Kernel function for collisional coalescence"
    KernelFunc::Any
    "Particle mass thresholds for analytical integration style"
    mass_thresholds::Vector{FT}
    "Polynomial order of the kernel approximation for coalescence"
    kernel_order::Int
    "Upper limit for evaluation of the polynomial kernel approximation"
    kernel_limit::FT
    "Coalescence data"
    coal_data::Any
    "Sedimentation rate parameters"
    vel::Vector{Tuple{FT, FT}}
    "Normalizing factors for number density and particle mass"
    norms::Vector{FT}
end

"""
    coalescence(clinfo)

 - `clinfo` - kwarg structure containing pdists, moments, and coalescence parameters
 TODO: currently implemented only for analytical coalescence style

Returns a vector of moment tendencies due to collisional coalescence
"""
function coalescence end

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
function condensation end

"""
    weighted_vt(clinfo)

 - `clinfo` - kwarg structure containing pdists, moments, and coalescence parameters
Returns the integrated fall speeds corresponding to the rate of change of prognostic moments
"""
function weighted_vt end

end
