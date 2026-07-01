import CloudMicrophysics.ThermodynamicsInterface as TDI
using CairoMakie

# Thermodynamics handles (a representative air density is fixed throughout).
const FT = Float64
const tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
const Tf = TDI.T_freeze(tps)
const cp_d = TDI.TD.Parameters.cp_d(tps)
const Lv = TDI.TD.Parameters.LH_v0(tps)
const Ls = TDI.TD.Parameters.LH_s0(tps)
const RHO = FT(0.8)

# Saturation specific humidity over each phase, recovered from the supersaturation interface
# (S = q_vap/q_sat - 1 evaluated at zero condensate, so q_sat = q_ref / (1 + S)).
const QREF = FT(1e-2)
q_sat_liq(T) = QREF / (1 + TDI.supersaturation_over_liquid(tps, QREF, FT(0), FT(0), RHO, T))
q_sat_ice(T) = QREF / (1 + TDI.supersaturation_over_ice(tps, QREF, FT(0), FT(0), RHO, T))

# Simplified 0-D vapor-exchange-only model: cloud liquid and cloud ice exchange mass with the shared vapor by
# relaxation toward their own saturation, with liquid kinetically fast and ice slow. Collection and precipitation
# are excluded so the Wegener-Bergeron-Findeisen transfer is isolated.
function wbf_box(; T0, q_liq0, q_ice0, tau_liq, tau_ice, t_end, dt)
    q_tot = q_sat_liq(T0) + q_liq0 + q_ice0           # start at liquid saturation
    q_liq, q_ice, T = q_liq0, q_ice0, T0
    n = Int(round(t_end / dt))
    ts = zeros(n + 1)
    QL = zeros(n + 1)
    QI = zeros(n + 1)
    SL = zeros(n + 1)
    SI = zeros(n + 1)
    for k in 0:n
        q_vap = q_tot - q_liq - q_ice
        ts[k + 1] = k * dt
        QL[k + 1] = q_liq
        QI[k + 1] = q_ice
        SL[k + 1] = q_vap / q_sat_liq(T) - 1
        SI[k + 1] = q_vap / q_sat_ice(T) - 1
        cond = (q_liq > 0 || q_vap > q_sat_liq(T)) ? (q_vap - q_sat_liq(T)) / tau_liq : 0.0
        dep = (q_ice > 0 || q_vap > q_sat_ice(T)) ? (q_vap - q_sat_ice(T)) / tau_ice : 0.0
        dq_liq = max(cond * dt, -q_liq)               # cannot evaporate more liquid than present
        dq_ice = max(dep * dt, -q_ice)
        q_liq += dq_liq
        q_ice += dq_ice
        T += (Lv * dq_liq + Ls * dq_ice) / cp_d
    end
    (; ts, QL, QI, SL, SI)
end

colors = Makie.wong_colors()
fig = Figure(size = (1200, 360))

r = wbf_box(; T0 = Tf - 15, q_liq0 = 1.5e-3, q_ice0 = 1e-5, tau_liq = 3.0, tau_ice = 50.0, t_end = 700.0, dt = 0.5)

# Panel (a): supersaturation -- the WBF signature. While liquid is present the fast liquid pins the vapor at water
# saturation (S_liq ~ 0), so the parcel stays supersaturated over ice (S_ice > 0) and ice deposits; once the liquid
# is gone the vapor relaxes to ice saturation (S_ice -> 0, S_liq < 0).
axa = Axis(fig[1, 1], xlabel = "time [s]", ylabel = "supersaturation",
    title = "(a) WBF supersaturation, T = T_f − 15 K")
lines!(axa, r.ts, r.SL, color = colors[1], linewidth = 2, label = "S_liq")
lines!(axa, r.ts, r.SI, color = colors[2], linewidth = 2, label = "S_ice")
hlines!(axa, [0.0], color = :gray, linestyle = :dot)
axislegend(axa, position = :rt, framevisible = false)

# Panel (b): the mass transfer -- liquid evaporates and ice grows at sustained S_ice > 0 (the WBF transfer).
axb = Axis(fig[1, 2], xlabel = "time [s]", ylabel = "specific content [g/kg]",
    title = "(b) WBF mass transfer")
lines!(axb, r.ts, r.QL .* 1e3, color = colors[1], linewidth = 2, label = "cloud liquid")
lines!(axb, r.ts, r.QI .* 1e3, color = colors[2], linewidth = 2, label = "cloud ice")
axislegend(axb, position = :rc, framevisible = false)

# Panel (c): the thermodynamic basis of the max(S_ice, S_liq) criterion. The two saturation curves cross at the
# freezing point; the limiter holds the vapor at the saturation of the more-supersaturated phase = the lower of the
# two q_sat curves (smaller q_sat is larger S). The binding phase changes from ice to liquid at T_f.
Ts = collect(range(Tf - 25, Tf + 8, length = 200))
qsl = q_sat_liq.(Ts) .* 1e3
qsi = q_sat_ice.(Ts) .* 1e3
axc = Axis(fig[1, 3], xlabel = "temperature [K]", ylabel = "saturation specific humidity [g/kg]",
    title = "(c) max-criterion floor")
lines!(axc, Ts, qsl, color = colors[1], linewidth = 2, label = "q_sat over liquid")
lines!(axc, Ts, qsi, color = colors[2], linewidth = 2, label = "q_sat over ice")
lines!(axc, Ts, min.(qsl, qsi), color = colors[6], linewidth = 4, label = "vapor floor = min(q_sat)")
vlines!(axc, [Tf], color = :gray, linestyle = :dash)
axislegend(axc, position = :lt, framevisible = false)

save("SatAdjustmentWBF.svg", fig)
