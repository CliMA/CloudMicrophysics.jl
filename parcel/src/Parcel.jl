"""
    Parcel

An adiabatic parcel model for testing and illustrating the CloudMicrophysics
aerosol activation, ice nucleation, and cloud condensate growth parameterizations.
"""
module Parcel

import CloudMicrophysics as CM

export parcel_params, run_parcel, distribution_moments, S_i, eᵥ, ξ

include("ParcelParameters.jl")
include("ParcelDistributions.jl")
include("ParcelCommon.jl")
include("ParcelTendencies.jl")
include("ParcelModel.jl")

end
