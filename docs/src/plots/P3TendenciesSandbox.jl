import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CairoMakie as Plt
import CloudMicrophysics.Microphysics1M as CM1
import Thermodynamics as TD

import Plots as PL

const PSP3 = CMP.ParametersP3

FT = Float64

p3 = CMP.ParametersP3(FT)
Chen2022 = CMP.Chen2022VelType(FT)

const liquid = CMP.CloudLiquid(FT)
const ice = CMP.CloudIce(FT)
const rain = CMP.Rain(FT)
const snow = CMP.Snow(FT)

const Blk1MVel = CMP.Blk1MVelType(FT)
const ce = CMP.CollisionEff(FT)

# Get expected N from n_0 as provided in the 1M scheme and q 
function getN(q, n_0)
    ρ_w = FT(1000)
    λ = (8 * n_0 * π * ρ_w / q)^(1 / 4)
    return n_0 / λ
end

qᵢ = FT(0.001)
qᵣ = FT(0.1)
q_c = FT(0.005)
Nᵢ = FT(1e8)
Nᵣ = FT(1e8)
N_c = FT(1e8)
ρ_r = FT(500)
F_r = FT(0.5)
ρ_a = FT(1.2)
T = FT(300)
ρ = ρ_a

q_rain_range = range(1e-8, stop = 0.005, length = 100)
q_snow_range = q_rain_range
n_0_rain = 16 * 1e6
n_0_ice = 2 * 1e7

PL.plot(
    q_rain_range * 1e3,
    [
        P3.ice_collisions(
            "rain",
            p3,
            Chen2022,
            qᵢ,
            Nᵢ,
            q,
            N_c,
            ρ_a,
            F_r,
            ρ_r,
            T,
        ) for q in q_rain_range
    ],
    linewidth = 3,
    xlabel = "q_precipitation [g/kg]",
    ylabel = "collision/accretion rate [1/s]",
    label = "P3 rain and ice",
    linestyle = :dash,
    color = :blue,
)

PL.plot!(
    q_rain_range * 1e3,
    [
        CM1.accretion(ice, rain, Blk1MVel.rain, ce, qᵢ, q_rai, ρ_a) for
        q_rai in q_rain_range
    ],
    linewidth = 3,
    label = "1 Moment rain and ice",
    color = :blue,
)

PL.plot!(
    q_snow_range * 1e3,
    [
        P3.ice_collisions(
            "cloud",
            p3,
            Chen2022,
            q,
            Nᵢ,
            q_c,
            N_c,
            ρ_a,
            F_r,
            ρ_r,
            T,
            FT(0.1),
        ) for q in q_snow_range
    ],
    label = "P3 cloudwater and ice",
    xlabel = "q_precipitation [g/kg]",
    ylabel = "collision/accretion rate [1/s]",
    linewidth = 3,
    linestyle = :dash,
    color = :red,
)

PL.plot!(
    q_snow_range * 1e3,
    [
        CM1.accretion(liquid, snow, Blk1MVel.snow, ce, q_c, q_sno, ρ_a) for
        q_sno in q_snow_range
    ],
    linewidth = 3,
    label = "1 Moment liquid and snow",
    color = :red,
)

PL.savefig("CollisionComparisons.svg")

const tps = TD.Parameters.ThermodynamicsParameters(FT)
const aps = CMP.AirProperties(FT)

# 1 Moment Parameters to match graph
T, p = 273.15 + 15, 90000.0
ϵ = 1.0 / TD.Parameters.molmass_ratio(tps)
p_sat = TD.saturation_vapor_pressure(tps, T, TD.Ice())
q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
q_tot = 15e-3
q_vap = 0.15 * q_sat
q_liq = 0.0
q_ice = q_tot - q_vap - q_liq
q = TD.PhasePartition(q_tot, q_liq, q_ice)
R = TD.gas_constant_air(tps, q)
ρ = p / R / T
T = 273.15

PL.plot(
    q_snow_range * 1e3,
    [
        CM1.snow_melt(snow, Blk1MVel.snow, aps, tps, q_sno, ρ, T + 2) for
        q_sno in q_snow_range
    ],
    xlabel = "q_snow [g/kg]",
    linewidth = 3,
    ylabel = "snow melt rate [1/s]",
    label = "1M: T=2C",
    color = :blue,
)
PL.plot!(
    q_snow_range * 1e3,
    [
        CM1.snow_melt(snow, Blk1MVel.snow, aps, tps, q_sno, ρ, T + 4) for
        q_sno in q_snow_range
    ],
    xlabel = "q_snow [g/kg]",
    linewidth = 3,
    label = "1M: T=4C",
    color = :red,
)
PL.plot!(
    q_snow_range * 1e3,
    [
        CM1.snow_melt(snow, Blk1MVel.snow, aps, tps, q_sno, ρ, T + 6) for
        q_sno in q_snow_range
    ],
    xlabel = "q_snow [g/kg]",
    linewidth = 3,
    label = "1M: T=6C",
    color = :green,
)

PL.plot!(
    q_snow_range * 1e3,
    [
        P3.p3_melt(p3, Chen2022, aps, tps, q_sno, Nᵢ, T + 2, ρ_a, F_r, ρ_r) for
        q_sno in q_snow_range
    ],
    xlabel = "q_snow [g/kg]",
    linewidth = 3,
    label = "P3: T=2C",
    linestyle = :dash,
    color = :blue,
)

PL.plot!(
    q_snow_range * 1e3,
    [
        P3.p3_melt(p3, Chen2022, aps, tps, q_sno, Nᵢ, T + 4, ρ_a, F_r, ρ_r) for
        q_sno in q_snow_range
    ],
    xlabel = "q_snow [g/kg]",
    linewidth = 3,
    label = "P3: T=4C",
    linestyle = :dash,
    color = :red,
)

PL.plot!(
    q_snow_range * 1e3,
    [
        P3.p3_melt(p3, Chen2022, aps, tps, q_sno, Nᵢ, T + 6, ρ_a, F_r, ρ_r) for
        q_sno in q_snow_range
    ],
    xlabel = "q_snow [g/kg]",
    linewidth = 3,
    label = "P3: T=6C",
    linestyle = :dash,
    color = :green,
)

PL.savefig("MeltRateComparisons.svg")
