"""
Flexible N-moment microphysics representation, including: 
  - Generalized collisional-coalescence described by a kernel K(x,y), with
  default parameterization based on Long collision kernel 
  - Power-series representation of fall-speed 
  - Power-series representation of condensational growth
  TODO: no representation of ventilation effects
"""
module MicrophysicsFlexible

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
end

function coalescence end

function condensation end

function sedimentation end

end
