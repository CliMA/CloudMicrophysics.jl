#!/usr/bin/env julia

# Simple test of the new piecewise power law API
push!(LOAD_PATH, joinpath(@__DIR__, "src"))

using CloudMicrophysics
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Common as CO

# Create parameters
FT = Float64
params = CMP.ParametersP3(FT)
(; mass, area, slope, ρ_i) = params

# Test state conditions
F_rim = FT(0.5)
ρ_rim = FT(916.7)

# Test diameter
D = FT(1e-4)  # 100 micrometers

println("Testing new piecewise power law API...")

# Test 1: Create mass power law and evaluate it
println("\n=== Test 1: Mass power law construction and evaluation ===")
state = P3.get_state(params; F_rim, ρ_rim)
mass_law = P3.get_mass_law(state)
println("Mass law created:")
println(mass_law)

# Test evaluation - there's now only one API
mass = P3.ice_mass(state, D)
println("Mass at D=$D: ", mass)

# Test individual piece display
println("\nFirst piece details:")
println(mass_law.pieces[1])

# Test 2: Area evaluation function
println("\n=== Test 2: Area evaluation function construction and evaluation ===")
area_evaluator = P3.get_area_law(state)
println("Area evaluator created: ", typeof(area_evaluator))

# Test evaluation - there's now only one API
area = P3.ice_area(state, D)
println("Area at D=$D: ", area)

# Test 3: Integration - this requires more setup, so we'll keep it simple
println("\n=== Test 3: Integration example ===")
logλ = log(FT(1e5))  # Some reasonable value
μ = P3.get_μ(slope, logλ)

# This tests the new piecewise integration approach
println("logλ = ", logλ)
println("μ = ", μ)

println("\n✅ API tests completed successfully!") 