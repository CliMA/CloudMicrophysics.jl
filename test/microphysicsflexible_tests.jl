import Test as TT

import Thermodynamics as TD
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsFlexible as CMF
import Cloudy.ParticleDistributions as CPD
import Cloudy.KernelFunctions as CPK

@info "Microphysics Tests"

function test_microphysics_flexible(FT)
    # Thermodynamics and air properties parameters
    aps = CMP.AirProperties(FT)
    tps = TD.Parameters.ThermodynamicsParameters(FT)

    TT.@testset "Flexible microphysics - unit tests" begin
        # Create Cloudy struct with different distributions 
        clinfo = CMF.CLSetup{FT}(;
            pdists = [
                CPD.LognormalPrimitiveParticleDistribution(
                    FT(10.0),
                    FT(1.0),
                    FT(1.0),
                ),
            ],
            mom = CPD.get_moments(
                CPD.LognormalPrimitiveParticleDistribution(
                    FT(10.0),
                    FT(1.0),
                    FT(1.0),
                ),
            ),
            NProgMoms = [3],
            KernelFunc = CPK.ConstantKernelFunction(FT(1)),
            mass_thresholds = [FT(100.0)],
            kernel_order = 4,
            kernel_limit = FT(100),
        )

        # Create Cloudy struct with defaults
        clinfo = CMF.CLSetup{FT}()
        TT.@test length(clinfo.pdists) == 2
        TT.@test length(clinfo.mom) == 5
        TT.@test length(clinfo.NProgMoms) == 2

        # Test coalescence
        TT.@test isnothing(clinfo.coal_data)
        dmom = CMF.coalescence(clinfo)
        TT.@test ~isnothing(clinfo.coal_data)
        TT.@test all(dmom[1:clinfo.NProgMoms[1]] .< FT(0))
        TT.@test all(dmom[(end - 1):end] .> FT(0))

        # Test evaporation
        T = FT(288.15)
        S = 0.99
        cond_evap = CMF.condensation(clinfo, aps, tps, T, S)
        # conserve number density
        TT.@test cond_evap[1] == FT(0)
        TT.@test cond_evap[clinfo.NProgMoms[1] + 1] == FT(0)
        # decrease mass density
        TT.@test cond_evap[2] < FT(0)
        TT.@test cond_evap[clinfo.NProgMoms[1] + 2] < FT(0)

        # Test condensation 
        S = 1.01
        cond_evap = CMF.condensation(clinfo, aps, tps, T, S)
        # conserve number density
        TT.@test cond_evap[1] == FT(0)
        TT.@test cond_evap[clinfo.NProgMoms[1] + 1] == FT(0)
        # increase cloud mass density
        TT.@test cond_evap[2] > FT(0)

        # Test sedimentation 
        weighted_vt = CMF.weighted_vt(clinfo)
        TT.@test all(weighted_vt .> FT(0))
    end
end

println("Testing Float64")
test_microphysics_flexible(Float64)

# println("Testing Float32")
# test_microphysics_flexible(Float32)
