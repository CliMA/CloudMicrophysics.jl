#ENV["CLIMACOMMS_DEVICE"] = "CUDA"

import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore: Fields, Domains, Geometry, Meshes, Spaces

import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMN

"""
    A helper function to create a ClimaCore space
"""
function make_space(
    ::Type{FT};
    zlim = (FT(0), FT(10000)),
    nelements = 1000,
) where {FT}
    column = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(zlim[1]),
        Geometry.ZPoint{FT}(zlim[2]);
        boundary_names = (:bottom, :top),
    )
    mesh = Meshes.IntervalMesh(column; nelems = nelements)
    space = Spaces.CenterFiniteDifferenceSpace(mesh)
    return space
end

"""
    Try to reproduce the setup of how terminal velocity is used in Atmos
"""
function set_sedimentation_precomputed_quantities(Y, p, t)

    (; w) = p
    (; params) = p

    @. w = CMN.terminal_velocity(
        params.liquid,
        params.Ch2022.rain,
        Y.ρ,
        max(0, Y.ρq / Y.ρ),
    )
    return nothing
end

function main(::Type{FT}) where {FT}

    Ch2022 = CMP.Chen2022VelType(FT)
    liquid = CMP.CloudLiquid(FT)
    ice = CMP.CloudIce(FT)
    rain = CMP.Rain(FT)
    snow = CMP.Snow(FT)

    params = (; liquid, ice, Ch2022)

    space_3d_ρq = make_space(FT)
    space_3d_ρ = make_space(FT)
    space_3d_w = make_space(FT)

    ρq = Fields.ones(space_3d_ρq) .* 1e-3
    ρ = Fields.ones(space_3d_ρ)
    w = Fields.zeros(space_3d_w)

    Y = (; ρq, ρ)
    p = (; w, params)

    t = 1

    set_sedimentation_precomputed_quantities(Y, p, t)
    return nothing
end

using Test
@testset "GPU inference failure" begin
    main(Float64)
end
