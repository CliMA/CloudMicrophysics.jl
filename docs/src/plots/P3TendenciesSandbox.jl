import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CairoMakie as Plt
import CloudMicrophysics.Microphysics1M as CM1

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
q_c = FT(0.0005)
Nᵢ = FT(1e8)
Nᵣ = FT(1e8)
N_c = FT(1e8)
ρ_r = FT(500)
F_r = FT(0.5)
ρ_a = FT(1.2)
T = FT(300)
ρ = ρ_a

q_rain_range = range(1e-8, stop = 0.005, length = 15)
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
            N_c,
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
