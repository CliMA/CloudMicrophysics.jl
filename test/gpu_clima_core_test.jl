#ENV["CLIMACOMMS_DEVICE"] = "CUDA"

import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC

import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMN
import CloudMicrophysics.Microphysics1M as CM1

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

    (; wₗ, wᵢ, wᵣ, wₛ) = p
    (; params) = p

    @. wₗ = CMN.terminal_velocity(
        params.liquid,
        params.STVel,
        Y.ρ,
        max(0, Y.ρq / Y.ρ),
    )
    @. wᵢ = CMN.terminal_velocity(
        params.ice,
        params.Ch2022.small_ice,
        Y.ρ,
        max(0, Y.ρq / Y.ρ),
    )
    @. wᵣ = CM1.terminal_velocity(
        params.rain,
        params.Ch2022.rain,
        Y.ρ,
        max(0, Y.ρq / Y.ρ),
    )
    @. wₛ = CM1.terminal_velocity(
        params.snow,
        params.Ch2022.large_ice,
        Y.ρ,
        max(0, Y.ρq / Y.ρ),
    )
    return nothing
end

function main_1d(::Type{FT}) where {FT}

    Ch2022 = CMP.Chen2022VelType(FT)
    STVel = CMP.StokesRegimeVelType(FT)
    liquid = CMP.CloudLiquid(FT)
    ice = CMP.CloudIce(FT)
    rain = CMP.Rain(FT)
    snow = CMP.Snow(FT)

    params = (; liquid, ice, rain, snow, Ch2022, STVel)

    space_1d_ρq = make_column(FT)
    space_1d_ρ = make_column(FT)
    space_1d_w = make_column(FT)

    ρq = CC.Fields.ones(space_1d_ρq) .* FT(1e-3)
    ρ = CC.Fields.ones(space_1d_ρ)
    wₗ = CC.Fields.zeros(space_1d_w)
    wᵢ = CC.Fields.zeros(space_1d_w)
    wᵣ = CC.Fields.zeros(space_1d_w)
    wₛ = CC.Fields.zeros(space_1d_w)

    Y = (; ρq, ρ)
    p = (; wₗ, wᵢ, wᵣ, wₛ, params)

    t = 1

    set_sedimentation_precomputed_quantities(Y, p, t)
    return nothing
end

function main_3d(::Type{FT}) where {FT}

    Ch2022 = CMP.Chen2022VelType(FT)
    STVel = CMP.StokesRegimeVelType(FT)
    liquid = CMP.CloudLiquid(FT)
    ice = CMP.CloudIce(FT)
    rain = CMP.Rain(FT)
    snow = CMP.Snow(FT)

    params = (; liquid, ice, rain, snow, Ch2022, STVel)

    space_3d_ρq = make_extruded_sphere(FT)
    space_3d_ρ = make_extruded_sphere(FT)
    space_3d_w = make_extruded_sphere(FT)

    ρq = CC.Fields.ones(space_3d_ρq) .* FT(1e-3)

    ρ = CC.Fields.ones(space_3d_ρ)
    wₗ = CC.Fields.zeros(space_3d_w)
    wᵢ = CC.Fields.zeros(space_3d_w)
    wᵣ = CC.Fields.zeros(space_3d_w)
    wₛ = CC.Fields.zeros(space_3d_w)

    Y = (; ρq, ρ)
    p = (; wₗ, wᵢ, wᵣ, wₛ, params)

    t = 1

    set_sedimentation_precomputed_quantities(Y, p, t)
    return nothing
end

import Test as TT

TT.@testset "ClimaCore GPU inference failure $nD $FT" for nD in ("1D", "3D"), FT in (Float64, Float32)
    main_nd = nD == "1D" ? main_1d : main_3d
    main_nd(FT)
end
nothing
