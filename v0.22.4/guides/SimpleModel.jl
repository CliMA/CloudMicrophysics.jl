# # Simple Model

# This guide shows how to build a simple model using `CloudMicrophysics.jl` parameterizations
# and `OrdinaryDiffEq.jl` solver.

# We start by importing the needed external packages and `CloudMicrophysics.jl` modules.
import OrdinaryDiffEq as ODE
import UnicodePlots as UP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics1M as CM1

nothing #hide

# We define the problem for the ordinary differential equations solver.
# The `dY` and `Y` are the tendencies and the state vector of the solved problem,
# `p` stores the additional parameters of the simulation and `t` is the simulation time.
# See the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) for more details
# about the solver interface.
function rain_formation(dY, Y, p, t)
    FT = eltype(Y) # Floating point precision type
    (; ρₐ, rain, liquid, ce, v_term) = p # Additional parameters passed through p

    qₗ = Y[1] # Cloud water specific humidity
    qᵣ = Y[2] # Rain water specific humidity

    acnv = CM1.conv_q_liq_to_q_rai(rain.acnv1M, qₗ) # Rain autoconversion rate
    accr = CM1.accretion(liquid, rain, v_term.rain, ce, qₗ, qᵣ, ρₐ) # Rain accretion rate

    dY[1] = -acnv - accr # Add the tendecies for cloud water
    dY[2] = acnv + accr  # and rain
end

nothing #hide

# We choose the simulation precision type and grab the default values of simulation parameters.
# We store the simulation parameters in the named tuple `p`.
FT = Float32

rain = CMP.Rain(FT) # Rain drop parameters for the 1-moment scheme
liquid = CMP.CloudLiquid(FT) # Cloud droplet parameters for the 1-moment scheme
ce = CMP.CollisionEff(FT) # Collision efficiencies
v_term = CMP.Blk1MVelType(FT) # Terminal velocity parameters
ρₐ = FT(1) # Air density
p = (; ρₐ, rain, liquid, ce, v_term)

nothing #hide

# Finally we define the simulation time and initial conditions for cloud and rain water specific humidities.
# We define the ODE problem, pass it to the solver, and visualize the results.
t₀ = FT(0)
t_end = FT(10 * 60)
TS = (t₀, t_end)

qₗ0 = FT(5e-3)
qᵣ0 = FT(0)
IC = [FT(qₗ0), FT(qᵣ0)]

problem = ODE.ODEProblem(rain_formation, IC, TS, p)
sol = ODE.solve(problem, ODE.Tsit5(), reltol = eps(FT), abstol = eps(FT))

plt = UP.lineplot(
    sol.t,
    sol[1, :] .* 1e3,
    name = "cloud",
    xlabel = "time [s]",
    ylabel = "q [g/kg]",
)
UP.lineplot!(plt, sol.t, sol[2, :] .* 1e3, name = "rain")
