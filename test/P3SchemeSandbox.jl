import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CLIMAParameters as CP
import Integrals as IN
import SpecialFunctions as SF
import RootSolvers as RS

FT = Float64

const PSP3 = CMP.ParametersP3
p3 = CMP.ParametersP3(FT)

"""
    N′(D, p)
    
 - D - maximum particle dimension
 - p - a tuple containing N_0, λ, μ (intrcept, slope, and shape parameters for N′ respectively)
 
 Returns the value of N′ 
 Eq. 2 in Morrison and Milbrandt (2015).   
"""
N′(D, p) = p.N_0 * exp(-p.λ * D)


"""
    q_helper(N_0, λ)

 - p3 - a struct with P3 scheme parameters
 - N_0 - intercept parameter of N′ gamma distribution
 - λ - slope parameter of N′ gamma distribution 
 - F_r - rime mass fraction (q_rim/q_i)
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]

Returns the prognostic mass mixing ratio
Eq. 5 in Morrison and Milbrandt (2015).
"""
function q_helper(p3::PSP3{FT}, N_0::FT, λ::FT, F_r::FT, ρ_r::FT) where {FT}
    th = P3.thresholds(p3, ρ_r, F_r)
    q′(D, p) = P3.p3_mass(p3, D, F_r, th) * N′(D, p)
    problem = IN.IntegralProblem(q′, 0, Inf, (N_0 = N_0, λ = λ))
    sol = IN.solve(problem, IN.HCubatureJL(), reltol = 1e-3, abstol = 1e-3)
    return FT(sol.u)
end

"""
    lambda_expected()

Returns the expected lambda to be compared to the non linear solver
"""
function λ_expected(q_i::FT, n_0::FT) where{FT}
    m_e = FT(3)
    r_0 = FT(1e-5)
    ρ = FT(1)       # density of air 
    ρ_i = FT(917)   # kg/m^3
    m_0 = FT(4/3 * π * ρ_i * r_0^3)
    X_m = FT(1)
    Δ_m = FT(0)

    return ((SF.gamma(m_e + Δ_m + 1) * X_m * m_0 * n_0)/(q_i * ρ * r_0 ^ (m_e + Δ_m)))^(1/(m_e + Δ_m + 1))

end

function test_1M(q_i, n_0, ρ_w)
    shape_problem(λ) = q_i - 4/3 * ρ_w * n_0 * (6* λ^(-4))
    λ = RS.find_zero(shape_problem, RS.SecantMethod(FT(1000), FT(100000)), RS.CompactSolution(), RS.RelativeSolutionTolerance(1e-3), 10,).root 

    println("True λ = ", λ_expected(q_i, n_0))  
    println("Modeled λ = ", λ)
end

function test_solver(p3, q, N, ρ_r, F_r) 
    (λ, N_0) = P3.distribution_parameter_solver(p3, q, N, ρ_r, F_r)
    println("λ solved = ", λ) 
    println("N_0 solved = ", N_0)
end

q_i = FT(1e-6)
n_0 = FT(16 * 1e6)
ρ_w = FT(1e3)

N = FT(1e8)
q = FT(1e-3)
q_r = FT(0.5 * 1e-6) 
B_r = FT(1/200 * 1e-4)
n_0 = FT(16 * 1e6)
ρ_w = FT(1e3)

F_r = q_r/q_i
ρ_r = ifelse(B_r == 0, FT(0), q_r/B_r)

#th = P3.thresholds(p3, ρ_r, F_r)
#λ_ex = λ_expected(q_i, n_0)
#q_calculated1 = P3.q_gamma(p3, F_r, n_0, λ_ex, 0.00191*λ_ex^0.8 - 2, th)

#N = n_0/λ_ex
#println("N_expected = ", N)
    # N = n_0 * λ^(-0.00191*λ_ex^(0.8) + 2 - 1) * SF.gamma(1 + 0.00191*λ_ex^0.8 - 2)

#q_calculated2 = P3.q_helper(p3, n_0, λ_ex, F_r, ρ_r)
#N_calculated = P3.N_helper(n_0, λ_ex)

#println("with gamma functions, q = ", q_calculated1)
#println("with integral, q = ", q_calculated2)
#println()

#println("n_0 entered = ", n_0)
#println("λ expected = ", λ_ex)

test_1M(q_i, n_0, ρ_w)
#test_solver(p3, q_calculated1, N, ρ_r, F_r)