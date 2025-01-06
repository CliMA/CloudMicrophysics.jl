#ENV["CLIMACOMMS_DEVICE"] = "CUDA"

import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC

import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMN

"""
    A helper function to create a ClimaCore 1d column space
"""
function make_column(::Type{FT}) where {FT}

    context = ClimaComms.SingletonCommsContext(ClimaComms.CUDADevice())
    #context = ClimaComms.context()

    vert_domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(FT(0)),
        CC.Geometry.ZPoint{FT}(FT(1000));
        boundary_names = (:bottom, :top),
    )
    vert_mesh = CC.Meshes.IntervalMesh(vert_domain; nelems = 1000)
    vert_topology = CC.Topologies.IntervalTopology(context, vert_mesh)
    vert_center_space = CC.Spaces.CenterFiniteDifferenceSpace(vert_topology)
    return vert_center_space
end

"""
    A helper function to create a ClimaCore extruded sphere space
"""
function make_extruded_sphere(::Type{FT}) where {FT}

    context = ClimaComms.SingletonCommsContext(ClimaComms.CUDADevice())
    #context = ClimaComms.context()

    # Define vertical
    # domain
    vert_domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(FT(0)),
        CC.Geometry.ZPoint{FT}(FT(1000));
        boundary_names = (:bottom, :top),
    )
    # mesh
    vert_mesh = CC.Meshes.IntervalMesh(vert_domain; nelems = 1000)
    # topology
    vert_topology = CC.Topologies.IntervalTopology(context, vert_mesh)
    # grid
    vert_grid = CC.Grids.FiniteDifferenceGrid(vert_topology)
    #vert_center_space = CC.Spaces.CenterFiniteDifferenceSpace(vert_topology)

    # Define horizontal:
    # domain
    horz_domain = CC.Domains.SphereDomain(FT(30))
    # mesh
    horz_mesh = CC.Meshes.EquiangularCubedSphere(horz_domain, 4)
    # topology
    horz_topology = CC.Topologies.Topology2D(
        context,
        horz_mesh,
        CC.Topologies.spacefillingcurve(horz_mesh),
    )
    # space
    horz_space = CC.Spaces.SpectralElementSpace2D(
        horz_topology,
        CC.Quadratures.GLL{3 + 1}();
        enable_bubble = true,
    )
    # grid
    horz_grid = CC.Spaces.grid(horz_space)

    # Define surface
    z_surface = zeros(horz_space)
    hypsography = CC.Hypsography.Flat()

    # Define grid
    deep = false
    grid = CC.Grids.ExtrudedFiniteDifferenceGrid(
        horz_grid,
        vert_grid,
        hypsography;
        deep,
    )

    # Define 3D space
    center_extruded_space = CC.Spaces.CenterExtrudedFiniteDifferenceSpace(grid)
    return center_extruded_space
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

function main_1d(::Type{FT}) where {FT}

    Ch2022 = CMP.Chen2022VelType(FT)
    liquid = CMP.CloudLiquid(FT)
    ice = CMP.CloudIce(FT)
    rain = CMP.Rain(FT)
    snow = CMP.Snow(FT)

    params = (; liquid, ice, Ch2022)

    space_1d_ρq = make_column(FT)
    space_1d_ρ = make_column(FT)
    space_1d_w = make_column(FT)

    ρq = CC.Fields.ones(space_1d_ρq) .* FT(1e-3)
    ρ = CC.Fields.ones(space_1d_ρ)
    w = CC.Fields.zeros(space_1d_w)

    Y = (; ρq, ρ)
    p = (; w, params)

    t = 1

    set_sedimentation_precomputed_quantities(Y, p, t)
    return nothing
end

function main_3d(::Type{FT}) where {FT}

    Ch2022 = CMP.Chen2022VelType(FT)
    liquid = CMP.CloudLiquid(FT)
    ice = CMP.CloudIce(FT)
    rain = CMP.Rain(FT)
    snow = CMP.Snow(FT)

    params = (; liquid, ice, Ch2022)

    space_3d_ρq = make_extruded_sphere(FT)
    space_3d_ρ = make_extruded_sphere(FT)
    space_3d_w = make_extruded_sphere(FT)

    ρq = CC.Fields.ones(space_3d_ρq) .* FT(1e-3)

    ρ = CC.Fields.ones(space_3d_ρ)
    w = CC.Fields.zeros(space_3d_w)

    Y = (; ρq, ρ)
    p = (; w, params)

    t = 1

    set_sedimentation_precomputed_quantities(Y, p, t)
    return nothing
end

using Test
@testset "GPU inference failure 1D Float64" begin
    main_1d(Float64)
end
@testset "GPU inference failure 3D Float64" begin
    main_3d(Float64)
end
@testset "GPU inference failure 1D Float32" begin
    main_1d(Float32)
end
@testset "GPU inference failure 3D Float32" begin
    main_3d(Float32)
end
