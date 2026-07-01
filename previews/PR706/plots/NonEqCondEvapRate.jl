using CairoMakie
CairoMakie.activate!(type = "svg")
import Thermodynamics as TD
import CloudMicrophysics as CM
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
    cm_params = CM.Parameters.CloudLiquid(FT(τ_relax), cm_params.ρw, cm_params.r_eff, cm_params.N_0) # overwrite τ_relax

    data = @. CMNe.conv_q_vap_to_q_lcl_icl_MM2015(cm_params, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T)
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
    cm_params = CM.Parameters.CloudLiquid(FT(τ_relax), cm_params.ρw, cm_params.r_eff, cm_params.N_0) # overwrite τ_relax

    data = @. CMNe.conv_q_vap_to_q_lcl_icl_MM2015(cm_params, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T')
    colorrange = extrema(data)
    @. data[iszero(data)] = NaN  # set zero values to NaN, then use `nan_color=:gray` to show them as gray. These are clipped values.
    colormap = spliced_cmap(:blues, :reds, colorrange...; mid = 0, categorical = true, symmetrize_color_ranges = true)
    S = @. calc_S(thp, T', ρ, qₜ, qₗ, qᵢ)

    file_name = "condensation_evaporation_$(x_suf)_$(y_suf).svg"

    return (; title, xlabel, ylabel, colorrange, colormap, file_name,
        x = qₗ / g_kg⁻¹, y = T, data, S,
    )
end

# --- Derivative line plots ---

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

function make_deriv_line_plot()
    # Temperature range
    T_range = collect(220:0.5:310)

    # Fixed parameters
    qₜ = FT(10e-3)    # 10 g/kg total water
    qₗ = FT(1e-3)     # 1 g/kg cloud liquid
    qᵢ = FT(1e-3)     # 1 g/kg cloud ice
    qᵣ = FT(0)
    qₛ = FT(0)
    ρ  = FT(1)        # kg/m³

    τ_relax = FT(1)    # 1 s relaxation time
    _liq = CM.Parameters.CloudLiquid(FT)
    cm_liq = CM.Parameters.CloudLiquid(τ_relax, _liq.ρw, _liq.r_eff, _liq.N_0)
    _ice = CM.Parameters.CloudIce(FT)
    cm_ice = CM.Parameters.CloudIce(_ice.pdf, _ice.mass, _ice.r0, _ice.r_ice_snow, τ_relax, _ice.ρᵢ, _ice.r_eff, _ice.N_0)

    # Compute tendency and derivatives for each temperature
    S_liq       = zeros(FT, length(T_range))
    dS_liq_simple = zeros(FT, length(T_range))
    dS_liq_full = zeros(FT, length(T_range))

    S_ice       = zeros(FT, length(T_range))
    dS_ice_simple = zeros(FT, length(T_range))
    dS_ice_full = zeros(FT, length(T_range))

    for (i, T) in enumerate(T_range)
        T_ft = FT(T)
        s = CMNe.conv_q_vap_to_q_lcl_icl_MM2015(
            cm_liq, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T_ft,
        )
        ds_s = CMNe.∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld(
            cm_liq, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T_ft,
        )
        ds_f = CMNe.∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld(
            cm_liq, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T_ft;
            simplified = false,
        )
        S_liq[i] = s
        dS_liq_simple[i] = ds_s
        dS_liq_full[i] = ds_f

        s = CMNe.conv_q_vap_to_q_lcl_icl_MM2015(
            cm_ice, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T_ft,
        )
        ds_s = CMNe.∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld(
            cm_ice, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T_ft,
        )
        ds_f = CMNe.∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld(
            cm_ice, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T_ft;
            simplified = false,
        )
        S_ice[i] = s
        dS_ice_simple[i] = ds_s
        dS_ice_full[i] = ds_f
    end

    T_freeze = FT(TD.Parameters.T_freeze(thp))

    # --- Figure: 2 panels ---
    fig = Figure(size = (900, 450))

    # Left panel: tendency vs T
    ax1 = Axis(fig[1, 1];
        xlabel = "Temperature (K)",
        ylabel = "Tendency (1/s)",
        title = "Condensation/Evaporation & Deposition/Sublimation\n(qₜ=10, qₗ=1, qᵢ=1 g/kg, ρ=1 kg/m³, τ=1 s)",
    )
    lines!(ax1, T_range, S_liq; color = :blue, label = "Liquid (cond/evap)")
    lines!(ax1, T_range, S_ice; color = :magenta, label = "Ice (dep/subl)")
    vlines!(ax1, [T_freeze]; color = :gray, linestyle = :dash, label = "T_freeze")
    hlines!(ax1, [0]; color = :black, linewidth = 0.5)
    text!(ax1, T_freeze + 2, minimum(S_ice) * 0.5;
        text = "INP limiter\n(ice dep → 0)",
        fontsize = 10, color = :gray,
    )
    axislegend(ax1; position = :lb, framevisible = false)

    # Right panel: derivative vs T
    ax2 = Axis(fig[1, 2];
        xlabel = "Temperature (K)",
        ylabel = "∂(tendency)/∂q (1/s)",
        title = "Derivative of tendency w.r.t. cloud condensate",
    )
    lines!(ax2, T_range, dS_liq_full;   color = :blue, label = "Liquid (full)")
    lines!(ax2, T_range, dS_ice_full;   color = :magenta, label = "Ice (full)")
    # Simple derivative is -1/τ for both phases (identical when τ_liq == τ_ice)
    lines!(ax2, T_range, dS_liq_simple; color = :black, linestyle = :dash, label = "Simple = −1/τ (both)")
    vlines!(ax2, [T_freeze]; color = :gray, linestyle = :dash, label = "T_freeze")
    text!(ax2, T_freeze + 2, minimum(dS_ice_full) * 0.5;
        text = "INP limiter\n(∂ice → 0)",
        fontsize = 10, color = :gray,
    )
    axislegend(ax2; position = :lt, framevisible = false)

    save("condensation_evaporation_deriv.svg", fig)
end

function make_deriv_vs_qvap_plot()
    # Available water vapor range (q_vap = q_tot - q_lcl - q_icl - q_rai - q_sno)
    # We vary q_tot while keeping q_lcl, q_icl fixed, so q_vap changes
    qₗ = FT(1e-3)     # 1 g/kg cloud liquid (fixed)
    qᵢ = FT(1e-3)     # 1 g/kg cloud ice (fixed)
    qᵣ = FT(0)
    qₛ = FT(0)
    ρ  = FT(1)        # kg/m³

    # q_tot range: from 2 g/kg (all condensate, no vapor) to 20 g/kg
    qₜ_range = collect(range(FT(2e-3), stop = FT(20e-3), length = 181))
    qᵥ_range = qₜ_range .- qₗ .- qᵢ  # available water vapor

    τ_relax = FT(1)
    _liq = CM.Parameters.CloudLiquid(FT)
    cm_liq = CM.Parameters.CloudLiquid(τ_relax, _liq.ρw, _liq.r_eff, _liq.N_0)
    _ice = CM.Parameters.CloudIce(FT)
    cm_ice = CM.Parameters.CloudIce(_ice.pdf, _ice.mass, _ice.r0, _ice.r_ice_snow, τ_relax, _ice.ρᵢ, _ice.r_eff, _ice.N_0)

    # Evaluate at two representative temperatures
    T_vals = [FT(250), FT(280)]
    T_labels = ["T=250K", "T=280K"]
    T_colors_liq = [:royalblue, :cornflowerblue]
    T_colors_ice = [:magenta, :hotpink]

    fig = Figure(size = (900, 450))

    ax1 = Axis(fig[1, 1];
        xlabel = "Water vapor, qᵥ (g/kg)",
        ylabel = "Tendency (1/s)",
        title = "Tendency vs available water vapor\n(qₗ=1, qᵢ=1 g/kg, ρ=1 kg/m³, τ=1 s)",
    )
    ax2 = Axis(fig[1, 2];
        xlabel = "Water vapor, qᵥ (g/kg)",
        ylabel = "∂(tendency)/∂q (1/s)",
        title = "Derivative vs available water vapor",
    )

    for (j, T) in enumerate(T_vals)
        S_liq     = zeros(FT, length(qₜ_range))
        dS_liq_f  = zeros(FT, length(qₜ_range))
        dS_liq_s  = zeros(FT, length(qₜ_range))
        S_ice     = zeros(FT, length(qₜ_range))
        dS_ice_f  = zeros(FT, length(qₜ_range))

        for (i, qₜ) in enumerate(qₜ_range)
            s = CMNe.conv_q_vap_to_q_lcl_icl_MM2015(
                cm_liq, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T,
            )
            ds_s = CMNe.∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld(
                cm_liq, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T,
            )
            ds_f = CMNe.∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld(
                cm_liq, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T;
                simplified = false,
            )
            S_liq[i] = s
            dS_liq_f[i] = ds_f
            dS_liq_s[i] = ds_s

            s = CMNe.conv_q_vap_to_q_lcl_icl_MM2015(
                cm_ice, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T,
            )
            ds_f = CMNe.∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld(
                cm_ice, thp, qₜ, qₗ, qᵢ, qᵣ, qₛ, ρ, T;
                simplified = false,
            )
            S_ice[i] = s
            dS_ice_f[i] = ds_f
        end

        qᵥ_gkg = qᵥ_range * 1e3

        # Left panel: tendency
        lines!(ax1, qᵥ_gkg, S_liq; color = T_colors_liq[j], label = "Liquid $(T_labels[j])")
        lines!(ax1, qᵥ_gkg, S_ice; color = T_colors_ice[j], label = "Ice $(T_labels[j])")

        # Right panel: derivative
        lines!(ax2, qᵥ_gkg, dS_liq_f; color = T_colors_liq[j], label = "Liquid $(T_labels[j])")
        lines!(ax2, qᵥ_gkg, dS_ice_f; color = T_colors_ice[j], label = "Ice $(T_labels[j])")

        if j == 1
            lines!(ax2, qᵥ_gkg, dS_liq_s; color = :black, linestyle = :dash, label = "Simple = −1/τ")
        end
    end

    hlines!(ax1, [0]; color = :black, linewidth = 0.5)
    axislegend(ax1; position = :lt, framevisible = false)
    axislegend(ax2; position = :lb, framevisible = false)

    save("condensation_evaporation_deriv_vs_qvap.svg", fig)
end

make_figure(HydrostaticBalance_qₗ_z())
make_figure(Range_qₗ_T())
make_deriv_line_plot()
make_deriv_vs_qvap_plot()

