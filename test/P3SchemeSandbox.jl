import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CLIMAParameters as CP

FT = Float64

const PSP3 = CMP.ParametersP3
p3 = CMP.ParametersP3(FT)

N = FT(100.0)
q = FT(1e-1)
F_r = FT(0.5)
ρ_r = FT(1e-4)
D = FT(1e-5)

#th = P3.thresholds(p3, ρ_r, F_r)
#mass = P3.p3_mass(p3, D, F_r, th)
(λ, N_0) = P3.distribution_parameter_solver(p3, q, N, ρ_r, F_r)
println("λ = ", λ)
println("N_0 = ", N_0)