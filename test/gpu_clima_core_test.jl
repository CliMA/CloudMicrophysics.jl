#ENV["CLIMACOMMS_DEVICE"] = "CUDA"

import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC

import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMN
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.P3Scheme as P3

"""
    A helper function to create a ClimaCore 1d column space
"""
make_column(::Type{FT}) where {FT} =
    CC.CommonGrids.ColumnGrid(FT;
        z_elem = 1000, z_min = 0, z_max = 1000,
        # device = ClimaComms.CUDADevice(),
    ) |> CC.Spaces.CenterFiniteDifferenceSpace

"""
    A helper function to create a ClimaCore extruded sphere space
"""
make_extruded_sphere(::Type{FT}) where {FT} =
    CC.CommonGrids.ExtrudedCubedSphereGrid(FT;
        z_elem = 1000, z_min = 0, z_max = 1000, radius = 30,
        h_elem = 4, n_quad_points = 3 + 1, enable_bubble = true,
        # device = ClimaComms.CUDADevice(),
    ) |> CC.Spaces.CenterExtrudedFiniteDifferenceSpace

"""
    Try to reproduce the setup of how terminal velocity is used in Atmos
"""
function set_sedimentation_precomputed_quantities(Y, p, t)
    (; wₗ, wᵢ, wᵣ, wₛ) = p
    (; liquid, STVel, ice, Ch2022, rain, snow) = p.params

    @. wₗ = CMN.terminal_velocity(liquid, STVel, Y.ρ, max(0, Y.ρq / Y.ρ))
    @. wᵢ = CMN.terminal_velocity(ice, Ch2022.small_ice, Y.ρ, max(0, Y.ρq / Y.ρ))
    @. wᵣ = CM1.terminal_velocity(rain, Ch2022.rain, Y.ρ, max(0, Y.ρq / Y.ρ))
    @. wₛ = CM1.terminal_velocity(snow, Ch2022.large_ice, Y.ρ, max(0, Y.ρq / Y.ρ))
    return nothing
end

get_params(::Type{FT}) where {FT} = (;
    liquid = CMP.CloudLiquid(FT),
    ice = CMP.CloudIce(FT),
    rain = CMP.Rain(FT),
    snow = CMP.Snow(FT),
    Ch2022 = CMP.Chen2022VelType(FT),
    STVel = CMP.StokesRegimeVelType(FT),
)

get_precomputed_quantities(::Type{FT}, space) where {FT} = (;
    params = get_params(FT),
    wₗ = zeros(space),
    wᵢ = zeros(space),
    wᵣ = zeros(space),
    wₛ = zeros(space),
)

get_state(::Type{FT}, space) where {FT} = (;
    ρq = ones(space) .* FT(1e-3),
    ρ = ones(space),
)

function main_1d(::Type{FT}) where {FT}
    space = make_column(FT)
    Y = get_state(FT, space)
    p = get_precomputed_quantities(FT, space)
    t = 1
    set_sedimentation_precomputed_quantities(Y, p, t)
    return nothing
end

function main_3d(::Type{FT}) where {FT}
    space = make_extruded_sphere(FT)
    Y = get_state(FT, space)
    p = get_precomputed_quantities(FT, space)
    t = 1
    set_sedimentation_precomputed_quantities(Y, p, t)
    return nothing
end

function rcemipii_z_mesh(::Type{FT}) where {FT}
    z_max = 33_000
    z_elem = 74
    boundary_layer =
        FT[0, 37, 112, 194, 288, 395, 520, 667, 843, 1062, 1331, 1664, 2055, 2505]
    n_bl = length(boundary_layer) - 1
    free_atmosphere = range(3000, FT(z_max), length = z_elem - n_bl)  # z_elem=74, z_max=33_000m --> 500m spacing
    CT = CC.Geometry.ZPoint{FT}
    faces = CT.([boundary_layer; free_atmosphere])
    z_domain = CC.Domains.IntervalDomain(CT(0), CT(z_max); boundary_names = (:bottom, :top))
    z_mesh = CC.Meshes.IntervalMesh(z_domain, faces)
    return z_mesh
end

function get_rcemipii_grid(::Type{FT}) where {FT}
    z_mesh = rcemipii_z_mesh(FT)
    x_max = y_max = 96_000
    grid = CC.CommonGrids.Box3DGrid(FT;
        z_elem = z_elem = length(z_mesh.faces) - 1, z_mesh,
        x_min = 0, x_max, y_min = 0, y_max, z_min = 0, z_max = 33_000,
        periodic_x = true, periodic_y = true,
        n_quad_points = 3 + 1, x_elem = 4, y_elem = 4,
        global_geometry = CC.Geometry.CartesianGlobalGeometry(),
        enable_bubble = true,
    )
    return grid
end

get_rcemipii_center_space(::Type{FT}) where {FT} = CC.Spaces.CenterExtrudedFiniteDifferenceSpace(get_rcemipii_grid(FT))

get_p3_fields(::Type{FT}, space) where {FT} = (;
    params = CMP.ParametersP3(FT),
    L_ice = ones(space) .* FT(1e-4),   # [kg/m³] ice mass content
    N_ice = ones(space) .* FT(1e4),    # [1/m³]  ice number concentration
    F_rim = zeros(space),              # [-]     rime mass fraction (no rime)
    ρ_rim = zeros(space),              # [kg/m³] rime density      (no rime)
    logλ = zeros(space),               # cache field for the result
)

"""
    Test that `get_distribution_logλ` can be broadcast over ClimaCore fields.
    Fields are initialised to physically-plausible non-zero values so that the
    root-finder has a well-posed problem.
"""
p3_logλ_1d(::Type{FT}) where {FT} = p3_logλ(make_column(FT))

p3_logλ_3d(::Type{FT}) where {FT} = p3_logλ(get_rcemipii_center_space(FT))

function p3_logλ(space)
    (; params, L_ice, N_ice, F_rim, ρ_rim, logλ) = get_p3_fields(FT, space)
    @. logλ = P3.get_distribution_logλ(params, L_ice, N_ice, F_rim, ρ_rim)

    @. logλ = P3.get_distribution_logλ(
        P3.P3State(params, L_ice, N_ice, F_rim, ρ_rim)
    )
    return logλ
end


import Test as TT

TT.@testset "ClimaCore GPU inference failure $nD $FT" for nD in ("1D", "3D"), FT in (Float64, Float32)
    TT.@testset "sedimentation velocities" begin
        main_nd = nD == "1D" ? main_1d : main_3d
        main_nd(FT)
    end
    TT.@testset "P3 logλ" begin
        p3_logλ_nd = nD == "1D" ? p3_logλ_1d : p3_logλ_3d
        p3_logλ_nd(FT)
    end
end

nothing
