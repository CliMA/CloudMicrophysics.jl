include("Parameters.jl")
include("Common.jl")
include("CommonTypes.jl")
include("AerosolModel.jl")

Mode_κ(
    1,
    2,
    400,
    (1,2,3),
    (4,5,6),
    (7,8,9),
    (5,6,7),
    3,
)

function test(
    ad::AerosolModel.AerosolDistribution{NTuple{N, T}},
    ) where {N, T <: AM.Mode_κ}
    
end