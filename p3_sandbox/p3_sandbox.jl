import OrdinaryDiffEq as ODE

import CLIMAParameters
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.HetIceNucleation as CMI_het
import CloudMicrophysics.Common as CMO
import Thermodynamics as TD

"""
    ODE problem definitions
"""
function p3_sandbox(dY, Y, p, t)

    # Numerical precision used in the simulation
    FT = eltype(Y)
    (; const_dt, tps, T, pₐ, qᵥ, qₗ, Nₗ, rₗ, aero_type) = p
    p3 = CMP.ParametersP3(FT)

    # Our state vector
    Nᵢ = Y[1]      # Ice number concentration
    qᵢ = Y[2]      # Ice mixing ratio
    qᵣ = Y[3]      # Rime mass mixing ratio
    Bᵣ = Y[4]      # Rime volume

    F_r = qᵣ / qᵢ
    ρ_r = ifelse(Bᵣ == 0, FT(0), qᵣ / Bᵣ)
    q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)

    Rₐ = TD.gas_constant_air(tps, TD.PhasePartition(qᵥ))
    R_v = TD.Parameters.R_v(tps)
    e = qᵥ * pₐ * R_v / Rₐ

    a_w = CMO.a_w_eT(tps, e, T)
    a_w_ice = CMO.a_w_ice(tps, T)
    Δa_w = a_w - a_w_ice
    J_immersion = CMI_het.ABIFM_J(aero_type, Δa_w)
    dNᵢ_dt = J_immersion * Nₗ * 4 * π * rₗ^2
    println("J = ", J_immersion)

    sol = P3.thresholds(p3, ρ_r, F_r)
    println(" ")
    println("D_cr = ", sol.D_cr)
    println("D_gr = ", sol.D_gr)
    println("ρ_g =  ", sol.ρ_g)
    println("ρ_dr = ", sol.ρ_d)

    # TODO - compute the corresponding N0 and λ
    # TODO - add ice nucleation source terms
    # TODO - allow for zero rimed mass and volume

    # Set tendencies
    dY[1] = dNᵢ_dt       # ice number concentration
    dY[2] = FT(0)        # ice mixing ratio
    dY[3] = FT(0)        # rime mixing ratio
    dY[4] = FT(0)        # rime volume
end

FT = Float64

const_dt = 1.0
t_0 = 0
t_end = 2

tps = TD.Parameters.ThermodynamicsParameters(FT) # thermodynamics free parameters
wps = CMP.WaterProperties(FT)

T = FT(251)
pₐ = FT(800 * 1e2)   # air pressure
ρₗ = wps.ρw
Nₗ = FT(500 * 1e3)   # number of cloud droplets
rₗ = FT(1e-6)        # radius of dropletes
qᵥ = FT(8.1e-4)      # mixing ratio of water vapor
qₗ = Nₗ * 4 / 3 * π * rₗ^3 * ρₗ / 1.2 # 1.2 should be ρₐ
aero_type = CMP.Illite(FT)

p = (; const_dt, tps, T, pₐ, qᵥ, qₗ, Nₗ, rₗ, aero_type)

Nᵢ_0 = 100 * 1e6
qᵢ_0 = 1e-3
qᵣ_0 = 1e-4
Bᵣ_0 = 1 / 200 * 1e-4

IC = [FT(Nᵢ_0), FT(qᵢ_0), FT(qᵣ_0), FT(Bᵣ_0)]

problem = ODE.ODEProblem(p3_sandbox, IC, (FT(t_0), FT(t_end)), p)
sol = ODE.solve(
    problem,
    ODE.Euler(),
    dt = const_dt,
    reltol = 10 * eps(FT),
    abstol = 10 * eps(FT),
)

println("number of ice crystals = ", sol[1, :])
