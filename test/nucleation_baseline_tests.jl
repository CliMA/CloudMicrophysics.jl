import Test as TT

import ClimaParams as CP

import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.HetIceNucleation as CMI_het

"""
    micro_thermo_at(tps, FT; T, ρ, S_i, n_ice)

Build the `(micro, thermo)` pair the process functions take, at a prescribed ice
supersaturation. The vapor content is solved from `S_i` rather than prescribed, so the
same construction works at any temperature.
"""
function micro_thermo_at(tps, FT; T, ρ, S_i, n_ice, q_lcl = FT(0), q_rai = FT(0), q_ice = FT(0))
    q_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
    q_vap = (1 + S_i) * q_sat_ice
    micro = (;
        q_tot = q_vap + q_lcl + q_rai + q_ice,
        q_lcl,
        q_rai,
        q_ice,
        n_ice,
    )
    thermo = (; ρ, T)
    return (micro, thermo)
end

function test_nucleation_baseline(FT)

    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true)
    p3 = mp.ice.scheme
    inp = mp.ice.ice_nucleation

    # A cold cirrus state, where the deposition slot lives: 215 K at about 215 hPa.
    T_cold = FT(215)
    ρ_cold = FT(0.348)

    TT.@testset "deposition closure defaults" begin
        # The default deposition closure is the exponential-in-supercooling spectrum.
        TT.@test inp isa CMP.ExponentialSupercoolingINP

        # The shipped values are Cooper (1986) in the Thompson et al. (2004) form, with the
        # reference implementation's activation window and ceiling.
        TT.@test inp.a ≈ FT(5)
        TT.@test inp.b ≈ FT(0.304)
        TT.@test inp.N_max == FT(1.0e5)
        TT.@test inp.T_thr == FT(258.15)
        TT.@test inp.S_thr == FT(0.05)
    end

    TT.@testset "nascent crystal single source" begin
        seed = CMP.ice_seed(p3)

        # The nascent diameter is 2 μm, so the starter mass is that of a 2 μm solid ice
        # sphere. Both consumers of the diameter read it from here.
        TT.@test p3.D_nuc == FT(2e-6)
        TT.@test seed.r_nuc == p3.D_nuc / 2
        TT.@test seed.ρ_i == p3.ρ_i
        TT.@test seed.m_nuc ≈ FT(3.8399e-15) rtol = FT(1e-4)
        TT.@test seed.m_nuc ≈ p3.ρ_i * FT(π) / 6 * p3.D_nuc^3 rtol = 2 * eps(FT)
    end

    TT.@testset "target spectrum and its ceiling" begin
        # The target never exceeds the ceiling, anywhere.
        for T in FT(180):FT(1):FT(280)
            TT.@test inp(T) <= inp.N_max
        end
        # It is at the ceiling well inside the cirrus range, and below it at the warm edge
        # of the activation window.
        TT.@test inp(FT(213)) == inp.N_max
        TT.@test inp(T_cold) == inp.N_max
        TT.@test inp(inp.T_thr) < inp.N_max
        # Below the ceiling the target increases as the air cools.
        TT.@test inp(FT(250)) > inp(FT(255)) > inp(FT(260))
    end

    TT.@testset "activation window: temperature" begin
        # No nucleation at or above the temperature threshold, with the supersaturation
        # condition satisfied, so that only the temperature gate can be closing the slot.
        for T in (inp.T_thr, inp.T_thr + FT(1), FT(270))
            (micro, thermo) =
                micro_thermo_at(tps, FT; T, ρ = FT(0.8), S_i = FT(0.2), n_ice = FT(0))
            rate = CMI_het.deposition_rate(inp, mp, tps, micro, thermo)
            TT.@test !CMI_het.is_active(inp, T, FT(0.2))
            TT.@test rate.∂ₜn_frz == FT(0)
            TT.@test rate.∂ₜq_frz == FT(0)
        end
        # Just colder than the threshold, the slot is open.
        T_open = inp.T_thr - FT(0.01)
        (micro, thermo) =
            micro_thermo_at(tps, FT; T = T_open, ρ = FT(0.8), S_i = FT(0.2), n_ice = FT(0))
        rate = CMI_het.deposition_rate(inp, mp, tps, micro, thermo)
        TT.@test CMI_het.is_active(inp, T_open, FT(0.2))
        TT.@test rate.∂ₜn_frz > FT(0)
    end

    TT.@testset "activation window: ice supersaturation" begin
        # No nucleation below the supersaturation threshold, at a temperature well inside
        # the open half of the temperature gate.
        for S_i in (FT(-0.1), FT(0), FT(0.049))
            (micro, thermo) =
                micro_thermo_at(tps, FT; T = T_cold, ρ = ρ_cold, S_i, n_ice = FT(0))
            rate = CMI_het.deposition_rate(inp, mp, tps, micro, thermo)
            TT.@test !CMI_het.is_active(inp, T_cold, S_i)
            TT.@test rate.∂ₜn_frz == FT(0)
            TT.@test rate.∂ₜq_frz == FT(0)
        end
        # At the threshold itself the GATE is open: the condition is `S_i ≥ S_thr`.
        TT.@test CMI_het.is_active(inp, T_cold, inp.S_thr)
        # The RATE is asserted a hair above it, not at it. `deposition_rate` does not receive
        # `S_i`; it recomputes it as `q_vap / q_sat_ice - 1` from a `q_vap` the fixture builds
        # as `(1 + S_i) * q_sat_ice`. That multiply-then-divide round trip lands a hair ABOVE
        # the threshold at Float64 and a hair BELOW it at Float32, so an exact-boundary
        # assertion on the rate tests the rounding of the fixture rather than the closure.
        S_above = inp.S_thr + eps(FT) * 8
        (micro, thermo) =
            micro_thermo_at(tps, FT; T = T_cold, ρ = ρ_cold, S_i = S_above, n_ice = FT(0))
        rate = CMI_het.deposition_rate(inp, mp, tps, micro, thermo)
        TT.@test rate.∂ₜn_frz > FT(0)
    end

    TT.@testset "seed delivery time" begin
        # The delivery time is the diffusional growth time of a nascent crystal. At the
        # activation threshold in cold cirrus it is of order 10 s; the assertion is an
        # order-of-magnitude bound, because the value follows from the growth coefficient
        # and moves with the thermodynamics parameters.
        inv_τ = CMI_het.delivery_rate(inp, mp, tps, T_cold, inp.S_thr)
        τ_act = 1 / inv_τ
        TT.@test FT(1) < τ_act < FT(100)

        # It is the fundamental form: linear in the growth coefficient and in the
        # supersaturation, and inverse in the squared nascent radius.
        (; r_nuc, ρ_i) = CMP.ice_seed(p3)
        G_ice = CO.G_func_ice(mp.warm_rain.air_properties, tps, T_cold)
        TT.@test inv_τ ≈ 2 * G_ice * inp.S_thr / (ρ_i * r_nuc^2) rtol = 4 * eps(FT)

        # Delivery is faster where the air is more supersaturated, and it does not go
        # negative at or below ice saturation.
        TT.@test CMI_het.delivery_rate(inp, mp, tps, T_cold, FT(0.5)) >
                 CMI_het.delivery_rate(inp, mp, tps, T_cold, FT(0.1))
        TT.@test CMI_het.delivery_rate(inp, mp, tps, T_cold, FT(-0.3)) == FT(0)
    end

    TT.@testset "rate pair, depletion and disabling" begin
        (micro, thermo) =
            micro_thermo_at(tps, FT; T = T_cold, ρ = ρ_cold, S_i = FT(0.3), n_ice = FT(0))
        rate = CMI_het.deposition_rate(inp, mp, tps, micro, thermo)
        (; m_nuc) = CMP.ice_seed(p3)

        # Every crystal is created at the starter mass: there is no vapor-excess branch to
        # scale one moment without the other.
        TT.@test rate.∂ₜq_frz == m_nuc * rate.∂ₜn_frz
        TT.@test rate.∂ₜn_frz > FT(0)

        # The rate relaxes toward the target: it falls as the already-activated number
        # rises, and reaches zero once the target is met.
        n_target = inp(T_cold) / ρ_cold
        (micro_half, _) = micro_thermo_at(
            tps, FT; T = T_cold, ρ = ρ_cold, S_i = FT(0.3), n_ice = n_target / 2,
        )
        (micro_full, _) = micro_thermo_at(
            tps, FT; T = T_cold, ρ = ρ_cold, S_i = FT(0.3), n_ice = 2 * n_target,
        )
        rate_half = CMI_het.deposition_rate(inp, mp, tps, micro_half, thermo)
        rate_full = CMI_het.deposition_rate(inp, mp, tps, micro_full, thermo)
        TT.@test FT(0) < rate_half.∂ₜn_frz < rate.∂ₜn_frz
        TT.@test rate_full.∂ₜn_frz == FT(0)
        TT.@test rate_full.∂ₜq_frz == FT(0)

        # `nothing` disables the slot with an exactly zero pair.
        off = CMI_het.deposition_rate(nothing, mp, tps, micro, thermo)
        TT.@test off.∂ₜn_frz == FT(0)
        TT.@test off.∂ₜq_frz == FT(0)

        # Type stability.
        TT.@test eltype(rate.∂ₜn_frz) == FT
        TT.@test eltype(rate.∂ₜq_frz) == FT
        TT.@test eltype(off.∂ₜn_frz) == FT
    end

    TT.@testset "cloud droplet immersion is Bigg alone" begin
        hom = CMP.Koop2000(FT)
        rf = mp.ice.rain_freezing
        aps = mp.warm_rain.air_properties
        pdf_c = mp.ice.cloud_pdf
        T_freeze = TDI.T_freeze(tps)

        T = FT(255)
        ρ = FT(0.9)
        q_lcl = FT(5e-4)
        N_lcl = FT(1e8)
        qᵥ = FT(2e-3)
        r = CMI_het.cloud_freezing_rate(
            rf, hom, p3.vent, aps, tps, pdf_c, q_lcl, ρ, N_lcl, T, qᵥ,
        )

        # The heterogeneous coefficient is the Barklie-Gokhale form of Bigg (1953),
        # unmodified: no ice nucleating particle budget stands above it.
        TT.@test r.J_het == rf(T, T_freeze)
        # And no bound is carried alongside it.
        TT.@test !hasproperty(r, :J_cap)
        TT.@test r.∂ₜn_frz > FT(0)
        TT.@test isfinite(r.∂ₜq_frz)
    end
end

TT.@testset "P3 ice nucleation baseline ($FT)" for FT in (Float64, Float32)
    test_nucleation_baseline(FT)
end
nothing
