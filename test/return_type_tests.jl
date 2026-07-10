using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.MicrophysicsNonEq as CMNonEq
import CloudMicrophysics.DistributionTools as DT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Utilities as UT
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ForwardDiff as FD

# Functions must return CONCRETE types for every mix of plain-float and Dual
# arguments: fallback returns typed from a single argument produce silent,
# heap-boxed unions otherwise (see Utilities.promote_typeof). Each entry below
# is checked with every single-Dual signature and the all-Dual signature.

const FT64 = Float64
const D8 = FD.Dual{Nothing, Float64, 8}

function concrete_for_all_mixes(f, pre::Tuple, n::Int; post::Tuple = ())
    ok = true
    for i in 0:n  # 0 = all-Dual; i >= 1 = Dual only in slot i
        numeric = ntuple(j -> (i == 0 || i == j) ? D8 : FT64, n)
        rts = Base.return_types(f, Tuple{pre..., numeric..., post...})
        ok &= length(rts) >= 1 && all(isconcretetype, rts)
    end
    return ok
end

@testset "early returns are concretely typed under mixed Dual/plain args" begin
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT64)
    mp = CMP.Microphysics2MParams(FT64; with_ice = true, is_limited = true)
    mpnl = CMP.Microphysics2MParams(FT64; with_ice = true, is_limited = false)
    sb = mp.warm_rain.seifert_beheng
    sbnl = mpnl.warm_rain.seifert_beheng
    aps = mp.warm_rain.air_properties

    P = typeof
    @test concrete_for_all_mixes(CM2.pdf_rain_parameters, (P(sb.pdf_r),), 3)
    @test concrete_for_all_mixes(CM2.pdf_rain_parameters, (P(sbnl.pdf_r),), 3)
    @test concrete_for_all_mixes(CM2.pdf_rain_parameters_mass, (P(sb.pdf_r),), 3)
    @test concrete_for_all_mixes(CM2.log_pdf_cloud_parameters_mass, (P(sb.pdf_c),), 3)
    @test concrete_for_all_mixes(CM2.pdf_cloud_parameters_mass, (P(sb.pdf_c),), 3)
    @test concrete_for_all_mixes(CM2.pdf_cloud_parameters, (P(sb.pdf_c),), 3)
    @test concrete_for_all_mixes(CM2.get_size_distribution_bounds, (P(sb.pdf_r),), 3)
    @test concrete_for_all_mixes(CM2.get_size_distribution_bounds, (P(sb.pdf_c),), 3)
    @test concrete_for_all_mixes(CM2.autoconversion, (P(sb.acnv), P(sb.pdf_c)), 4)
    @test concrete_for_all_mixes(CM2.accretion, (P(sb),), 4)
    @test concrete_for_all_mixes(CM2.cloud_liquid_self_collection, (P(sb.acnv), P(sb.pdf_c)), 3)
    @test concrete_for_all_mixes(CM2.rain_self_collection, (P(sb.pdf_r), P(sb.self)), 3)
    @test concrete_for_all_mixes(CM2.rain_breakup, (P(sb.pdf_r), P(sb.brek)), 4)
    @test concrete_for_all_mixes(CM2.rain_self_collection_and_breakup, (P(sb),), 3)
    @test concrete_for_all_mixes(CM2.rain_evaporation, (P(sb), P(aps), P(tps)), 8)
    numadj = (; sb.numadj.τ, x_min = sb.pdf_c.xc_min, x_max = sb.pdf_c.xc_max)
    @test concrete_for_all_mixes(CM2.number_tendency_from_mass_limits, (P(numadj),), 2)
    @test concrete_for_all_mixes(
        CM2.rain_terminal_velocity, (P(sb), P(CMP.SB2006VelType(FT64))), 3,
    )
    @test concrete_for_all_mixes(
        CM2.rain_terminal_velocity, (P(sb), P(mp.ice.terminal_velocity.rain)), 3,
    )
    @test concrete_for_all_mixes(
        CM2.cloud_terminal_velocity, (P(sb.pdf_c), P(CMP.StokesRegimeVelType(FT64))), 3,
    )
    @test concrete_for_all_mixes(
        CO.Chen2022_exponential_pdf, (), 4; post = (Int,),
    )
    @test concrete_for_all_mixes(CO.logistic_function_integral, (), 3)
    @test concrete_for_all_mixes(UT._regularised_ratio, (), 2)
    @test concrete_for_all_mixes(
        CMNonEq.τ_relax,
        (P(CMP.CloudIce(FT64)), P(CMP.AirProperties(FT64)), P(mp.ice.ice_nucleation)), 2,
    )
    @test concrete_for_all_mixes(P3.loggamma_inc_moment, (), 4)
    # P3 aspect ratio: plain state with Dual D and vice versa
    st = P3.state_from_prognostic(mp.ice.scheme, FT64(1e-4), FT64(1e4), FT64(2e-5), FT64(4e-8))
    @test concrete_for_all_mixes(P3.ϕᵢ, (P(st),), 1)
    # ice terminal velocities: plain state, Dual ρₐ (fallback previously typed
    # off the state alone). Keyword wrapper -> probe the positional core.
    _glq = P3.GaussLegendre(FT64, 12)
    for f in (P3.ice_terminal_velocity_number_weighted, P3.ice_terminal_velocity_mass_weighted)
        g = (vp, ρ, s, l) -> f(vp, ρ, s, l; quad = _glq)
        rts = Base.return_types(g, Tuple{P(mp.ice.terminal_velocity), D8, P(st), FT64})
        @test all(isconcretetype, rts)
    end

    # LclRaiRates: mixed-type construction promotes instead of MethodError
    r = CM2.LclRaiRates(D8(1.0), 0.0, D8(2.0), 0.0)
    @test r isa CM2.LclRaiRates{D8}

    # size-distribution closures: zero branch must match the main branch when
    # the node D and the captured parameters mix plain/Dual
    for (pdf, q) in ((sb.pdf_r, FT64(1e-4)), (sb.pdf_c, FT64(1e-3)))
        n_plain = DT.size_distribution(pdf, q, FT64(1.1), FT64(1e7))
        @test isconcretetype(only(Base.return_types(n_plain, Tuple{D8})))
        n_dual = DT.size_distribution(pdf, D8(q), FT64(1.1), D8(1e7))
        @test isconcretetype(only(Base.return_types(n_dual, Tuple{FT64})))
        n_zero = DT.size_distribution(pdf, zero(FT64), FT64(1.1), zero(FT64))
        @test n_zero(D8(1e-4)) isa D8  # zero branch, Dual node
    end

    # NOTE: the 1M Instantaneous entry is NOT asserted here: its process
    # leaves (Microphysics1M.jl) carry their own single-argument-typed
    # returns. The 1M leaves share this pattern and are out of scope for
    # this 2M/P3 sweep (the BMT-level source-terms FT is fixed).
end
