using CairoMakie
CairoMakie.activate!(type = "svg")
import Thermodynamics as TD
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMNe

include(joinpath(@__DIR__, "plotting_utilities.jl"))  # spliced_cmap

FT = Float64
thp = TD.Parameters.ThermodynamicsParameters(FT)
g_kg⁻¹ = 1e-3

qₗ = -10g_kg⁻¹:(0.01g_kg⁻¹):10g_kg⁻¹; xlabel = "liquid water humidity, q_liq (g/kg)"; x_suf = "ql"
x_sc = qₗ / g_kg⁻¹

function calc_S(thp, T, ρ, q_tot, q_liq, q_ice)
    qᵥ = q_tot - q_liq - q_ice
    pₛᵥ = TD.saturation_vapor_pressure(thp, T, TD.Liquid())
    qₛₗ = TD.q_vap_from_p_vap(thp, T, ρ, pₛᵥ)
    return qᵥ / qₛₗ
end

struct HydrostaticBalance_qₗ_z end
function generate_cond_evap_rate(::HydrostaticBalance_qₗ_z)
    title = "Condensation/evaporation rate (g/kg/s) \n"

    # x-axis, qₗ
    qₗ = -10g_kg⁻¹:(0.01g_kg⁻¹):10g_kg⁻¹; xlabel = "liquid water humidity, q_liq (g/kg)"; x_suf = "ql"

    # y-axis, z (parameterizing T, ρ)
    title *= "assuming hydrostatic balance ("
    Γ = 6.5; title *= "Γ=$(Γ)K/km, "  # K/km, lapse rate
    R = 287  # J/(kg·K), specific gas constant for dry air
    g = 9.81  # m/s², acceleration due to gravity
    T₀ = 295; title *= "T₀=$(T₀)K, "  # K, at z=0
    p₀ = 100_000; title *= "p₀=$(Int(p₀/100))hPa)\n"  # Pa, at z=0
    ρ₀ = p₀ / (R * T₀) # kg/m³, at z=0
    
    T_fn(z) = T₀ - Γ * z
    ρ_fn(z) = ρ₀ * (T_fn(z) / T₀)^(g/R/Γ - 1)

    z = 0:0.1:15; ylabel = "height, z (km)"; y_suf = "z"

    T = T_fn.(z)'
    ρ = ρ_fn.(z)'

    # Fixed values
    qₜ_gkg = 10; qₜ = qₜ_gkg * g_kg⁻¹; title *= "qₜ=$(qₜ_gkg)g/kg, "
    qᵢ = qᵣ = qₛ = 0g_kg⁻¹; title *= "qᵢ=qᵣ=qₛ=0, "
    ρ = 1; title *= "ρ=$(ρ)kg/m³, "
    dt = 1; title *= "dt=$(dt)s, "
    τ_relax = 1; title *= "τ_relax=$(τ_relax)s, "
    cm_params = CM.Parameters.CloudLiquid(FT)
    clf = CMP.CloudLiquidFormation()

    data = broadcast(qₗ, T) do q_lcl, temp
        CMNe.conv_q_vap_to_q_lcl(
            clf,
            (; cloud = (; liquid = cm_params),
                process_params = (; cloud_liquid_formation = (; τ_relax = FT(τ_relax)))),
            thp,
            (; q_tot = qₜ, q_lcl, q_icl = qᵢ, q_rai = qᵣ, q_sno = qₛ),
            (; ρ, T = temp)
        )
    end
    colorrange = extrema(data)
    @. data[iszero(data)] = NaN  # set zero values to NaN, then use `nan_color=:gray` to show them as gray. These are clipped values.
    colormap = spliced_cmap(:blues, :reds, colorrange...; mid = 0, categorical = true, symmetrize_color_ranges = true)
    S = @. calc_S(thp, T, ρ, qₜ, qₗ, qᵢ)

    file_name = "condensation_evaporation_$(x_suf)_$(y_suf).svg"

    return (; title, xlabel, ylabel, colorrange, colormap, file_name,
        x = qₗ / g_kg⁻¹, y = z, data, S,
    )
end

struct Range_qₗ_T end
function generate_cond_evap_rate(::Range_qₗ_T)
    title = "Condensation/evaporation rate (g/kg/s) \n"

    # x-axis, qₗ
    qₗ = -10g_kg⁻¹:(0.01g_kg⁻¹):10g_kg⁻¹; xlabel = "liquid water humidity, q_liq (g/kg)"; x_suf = "ql"

    # y-axis, T
    T = 240:0.01:273.15; ylabel = "Temperature (K)"; y_suf = "T"

    # Fixed values
    qₜ_gkg = 10; qₜ = qₜ_gkg * g_kg⁻¹; title *= "qₜ=$(qₜ_gkg)g/kg, "
    qᵢ = qᵣ = qₛ = 0g_kg⁻¹; title *= "qᵢ=qᵣ=qₛ=0, "
    ρ = 1; title *= "ρ=$(ρ)kg/m³, "
    dt = 1; title *= "dt=$(dt)s, "
    τ_relax = 1; title *= "τ_relax=$(τ_relax)s, "
    cm_params = CM.Parameters.CloudLiquid(FT)
    clf = CMP.CloudLiquidFormation()

    data = broadcast(qₗ, T') do q_lcl, temp
        CMNe.conv_q_vap_to_q_lcl(
            clf,
            (; cloud = (; liquid = cm_params),
                process_params = (; cloud_liquid_formation = (; τ_relax = FT(τ_relax)))),
            thp,
            (; q_tot = qₜ, q_lcl, q_icl = qᵢ, q_rai = qᵣ, q_sno = qₛ),
            (; ρ, T = temp)
        )
    end
    colorrange = extrema(data)
    @. data[iszero(data)] = NaN  # set zero values to NaN, then use `nan_color=:gray` to show them as gray. These are clipped values.
    colormap = spliced_cmap(:blues, :reds, colorrange...; mid = 0, categorical = true, symmetrize_color_ranges = true)
    S = @. calc_S(thp, T', ρ, qₜ, qₗ, qᵢ)

    file_name = "condensation_evaporation_$(x_suf)_$(y_suf).svg"

    return (; title, xlabel, ylabel, colorrange, colormap, file_name,
        x = qₗ / g_kg⁻¹, y = T, data, S,
    )
end


minussign(s) = replace(s, "-" => "−")

### Generate figures ###
function make_figure(case)
    (; title, xlabel, ylabel, colorrange, colormap, file_name,
        x, y, data, S,
    ) = generate_cond_evap_rate(case)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel, ylabel, title)
    hm = heatmap!(ax, x, y, data; colorrange, colormap, nan_color = :gray, rasterize = 2)
    contour!(ax, x, y, S; labels = true, color = :black,
        levels = [0.5, 1, 1.5],
        labelformatter = l -> "S = " * minussign(Makie.contour_label_formatter(l)),
    )
    label =
        rich("Condensation(+)", color = :red) * " / " * rich("Evaporation(−)", color = :blue) * " rate (g/kg/s) " *
        rich("(gray≡0)", color = :gray)
    Colorbar(fig[1, 2], hm; colorrange, label)
    save(file_name, fig)
end

make_figure(HydrostaticBalance_qₗ_z())
make_figure(Range_qₗ_T())

