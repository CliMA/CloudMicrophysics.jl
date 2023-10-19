import Test as TT

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CMC
import CloudMicrophysics.Microphysics2M as CM2

@info "Microphysics Tests"

function test_microphysics2M(FT)

    # Different 2-moment autoconversion and accretion parameters
    KK2000 = CMP.KK2000(FT)
    B1994 = CMP.B1994(FT)
    TC1980 = CMP.TC1980(FT)
    LD2004 = CMP.LD2004(FT)
    VarTSc = CMP.VarTimescaleAcnv(FT)

    # Seifert and Beheng 2006 parameters
    SB2006 = CMP.SB2006(FT)

    # Thermodynamics and air properties parameters
    aps = CMP.AirProperties(FT)
    tps = CMP.ThermodynamicsParameters(FT)

    # Terminal velocity parameters
    SB2006Vel = CMP.SB2006VelType(FT)
    Chen2022Vel = CMP.Chen2022VelTypeRain(FT)

    TT.@testset "2M_microphysics - unit tests" begin

        ρ = FT(1)

        # no reference data available - checking if callable and not NaN
        q_liq = FT(0.5e-3)
        q_rai = FT(1e-6)

        TT.@test CM2.accretion(KK2000, q_liq, q_rai, ρ) != NaN
        TT.@test CM2.accretion(B1994, q_liq, q_rai, ρ) != NaN
        TT.@test CM2.accretion(TC1980, q_liq, q_rai) != NaN
        TT.@test CM2.conv_q_liq_to_q_rai(VarTSc, q_liq, ρ) != NaN

        # output should be zero if either q_liq or q_rai are zero
        q_liq = FT(0)
        q_rai = FT(1e-6)

        TT.@test CM2.conv_q_liq_to_q_rai(VarTSc, q_liq, ρ) == FT(0)
        TT.@test CM2.conv_q_liq_to_q_rai(KK2000, q_liq, ρ) == FT(0)
        TT.@test CM2.conv_q_liq_to_q_rai(B1994, q_liq, ρ) == FT(0)
        TT.@test CM2.conv_q_liq_to_q_rai(TC1980, q_liq, ρ) == FT(0)
        TT.@test CM2.conv_q_liq_to_q_rai(LD2004, q_liq, ρ) == FT(0)
        TT.@test CM2.accretion(KK2000, q_liq, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(B1994, q_liq, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(TC1980, q_liq, q_rai) == FT(0)

        q_liq = FT(0.5e-3)
        q_rai = FT(0)
        TT.@test CM2.accretion(KK2000, q_liq, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(B1994, q_liq, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(TC1980, q_liq, q_rai) == FT(0)

        TT.@test CM2.conv_q_liq_to_q_rai(VarTSc, q_liq, ρ, N_d = FT(1e8)) >
                 CM2.conv_q_liq_to_q_rai(VarTSc, q_liq, ρ, N_d = FT(1e9))

        # far from threshold points, autoconversion with and without smooth transition should
        # be approximately equal
        q_liq = FT(0.5e-3)
        TT.@test CM2.conv_q_liq_to_q_rai(
            B1994,
            q_liq,
            ρ,
            smooth_transition = true,
        ) ≈ CM2.conv_q_liq_to_q_rai(
            B1994,
            q_liq,
            ρ,
            smooth_transition = false,
        ) rtol = 0.2
        TT.@test CM2.conv_q_liq_to_q_rai(
            TC1980,
            q_liq,
            ρ,
            smooth_transition = true,
        ) ≈ CM2.conv_q_liq_to_q_rai(
            TC1980,
            q_liq,
            ρ,
            smooth_transition = false,
        ) rtol = 0.2
        TT.@test CM2.conv_q_liq_to_q_rai(
            LD2004,
            q_liq,
            ρ,
            smooth_transition = true,
        ) ≈ CM2.conv_q_liq_to_q_rai(
            LD2004,
            q_liq,
            ρ,
            smooth_transition = false,
        ) rtol = 0.2

    end

    TT.@testset "2M_microphysics - compare with Wood_2005" begin

        ρ = FT(1)
        q_liq = FT(0.5e-3)

        # compare with Wood 2005 Fig 1 panel a
        function compare(scheme, input, output; eps = 0.1)
            TT.@test CM2.conv_q_liq_to_q_rai(scheme, input * FT(1e-3), ρ) ≈
                     output atol = eps * output
        end
        compare(KK2000, FT(0.03138461538461537), FT(2.636846054348105e-12))
        compare(KK2000, FT(0.8738461538461537), FT(9.491665962977648e-9))
        compare(
            B1994,
            FT(0.13999999999999999),
            FT(4.584323122458155e-12),
            eps = 1,
        )
        compare(
            B1994,
            FT(0.9000000000000006),
            FT(5.4940586176564715e-8),
            eps = 1,
        )
        compare(TC1980, FT(0.2700000000000001), FT(3.2768635256661366e-8))
        compare(TC1980, FT(0.9000000000000006), FT(5.340418612468997e-7))
        compare(LD2004, FT(0.3700000000000002), FT(8.697439193234471e-9))
        compare(LD2004, FT(0.9000000000000006), FT(1.1325570516983242e-7))

        # compare with Wood 2005 Fig 1 panel b
        function compare_Nd(scheme, input, output; eps = 0.1)
            TT.@test CM2.conv_q_liq_to_q_rai(
                scheme,
                q_liq,
                N_d = input * FT(1e6),
                ρ,
            ) ≈ output atol = eps * output
        end
        compare_Nd(KK2000, FT(16.13564081404141), FT(6.457285532394289e-8))
        compare_Nd(KK2000, FT(652.093931356625), FT(8.604011482409198e-11))
        compare_Nd(B1994, FT(14.47851799831075), FT(4.2829062386778675e-7))
        compare_Nd(B1994, FT(693.0425211336465), FT(6.076294746898778e-12))
        compare_Nd(TC1980, FT(13.658073017575544), FT(2.7110779872658386e-7))
        compare_Nd(TC1980, FT(205.0970632305975), FT(1.0928660431622176e-7))
        compare_Nd(LD2004, FT(15.122629721719655), FT(1.1647783461546477e-7))
        compare_Nd(
            LD2004,
            FT(149.01220754857331),
            FT(1.3917890403908125e-8),
            eps = 1,
        )

    end

    # 2M_microphysics - Seifert and Beheng 2006 double moment scheme tests
    TT.@testset "limiting lambda_r and x_r - Seifert and Beheng 2006" begin
        #setup
        q_rai = [FT(0), FT(1e-4), FT(1e-2)]
        N_rai = [FT(1e1), FT(1e3), FT(1e5)]
        ρ = FT(1)

        xr_min = FT(2.6e-10)
        xr_max = FT(5e-6)
        λ_min = FT(1e3)
        λ_max = FT(1e4)

        for Nr in N_rai
            for qr in q_rai
                #action
                λ = CM2.raindrops_limited_vars(SB2006.pdf, qr, ρ, Nr).λr
                xr = CM2.raindrops_limited_vars(SB2006.pdf, qr, ρ, Nr).xr

                #test
                TT.@test λ_min <= λ <= λ_max
                TT.@test xr_min <= xr <= xr_max
            end
        end

    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 autoconversion and liquid self-collection" begin
        #setup
        ρ = FT(1)
        q_liq = FT(0.5e-3)
        N_liq = FT(1e8)
        q_rai = FT(1e-6)

        kcc = FT(4.44e9)
        xstar = FT(2.6e-10)
        νc = FT(2.0)
        ρ0 = FT(1.225)

        #action
        au = CM2.autoconversion(SB2006.acnv, q_liq, q_rai, ρ, N_liq)
        sc = CM2.liquid_self_collection(SB2006.acnv, q_liq, ρ, au.dN_liq_dt)
        au_sc = CM2.autoconversion_and_liquid_self_collection(
            SB2006,
            q_liq,
            q_rai,
            ρ,
            N_liq,
        )

        Lc = ρ * q_liq
        Lr = ρ * q_rai
        xc = min(xstar, Lc / N_liq)
        τ = 1 - Lc / (Lc + Lr)
        ϕ_au = 400 * τ^0.7 * (1 - τ^0.7)^3
        dqrdt_au =
            kcc / 20 / xstar * (νc + 2) * (νc + 4) / (νc + 1)^2 *
            Lc^2 *
            xc^2 *
            (1 + ϕ_au / (1 - τ)^2) *
            (ρ0 / ρ) / ρ
        dqcdt_au = -dqrdt_au
        dNcdt_au = 2 / xstar * ρ * dqcdt_au
        dNrdt_au = -0.5 * dNcdt_au
        dNcdt_sc = -kcc * (νc + 2) / (νc + 1) * (ρ0 / ρ) * Lc^2 - au.dN_liq_dt

        #test
        TT.@test au isa CM2.LiqRaiRates
        TT.@test au.dq_liq_dt ≈ dqcdt_au rtol = 1e-6
        TT.@test au.dq_rai_dt ≈ dqrdt_au rtol = 1e-6
        TT.@test au.dN_liq_dt ≈ dNcdt_au rtol = 1e-6
        TT.@test au.dN_rai_dt ≈ dNrdt_au rtol = 1e-6
        TT.@test sc ≈ dNcdt_sc rtol = 1e-6
        TT.@test au_sc isa NamedTuple
        TT.@test au_sc.au.dq_liq_dt ≈ dqcdt_au rtol = 1e-6
        TT.@test au_sc.au.dq_rai_dt ≈ dqrdt_au rtol = 1e-6
        TT.@test au_sc.au.dN_liq_dt ≈ dNcdt_au rtol = 1e-6
        TT.@test au_sc.au.dN_rai_dt ≈ dNrdt_au rtol = 1e-6
        TT.@test au_sc.sc ≈ dNcdt_sc rtol = 1e-6

        #action
        au = CM2.autoconversion(SB2006.acnv, FT(0), FT(0), ρ, N_liq)
        sc = CM2.liquid_self_collection(SB2006.acnv, FT(0), ρ, au.dN_liq_dt)
        au_sc = CM2.autoconversion_and_liquid_self_collection(
            SB2006,
            FT(0),
            FT(0),
            ρ,
            N_liq,
        )

        #test
        TT.@test au.dq_liq_dt ≈ FT(0) atol = eps(FT)
        TT.@test au.dq_rai_dt ≈ FT(0) atol = eps(FT)
        TT.@test au.dN_liq_dt ≈ FT(0) atol = eps(FT)
        TT.@test au.dN_rai_dt ≈ FT(0) atol = eps(FT)
        TT.@test sc ≈ FT(0) atol = eps(FT)
        TT.@test au_sc.au.dq_liq_dt ≈ FT(0) atol = eps(FT)
        TT.@test au_sc.au.dq_rai_dt ≈ FT(0) atol = eps(FT)
        TT.@test au_sc.au.dN_liq_dt ≈ FT(0) atol = eps(FT)
        TT.@test au_sc.au.dN_rai_dt ≈ FT(0) atol = eps(FT)
        TT.@test au_sc.sc ≈ FT(0) atol = eps(FT)
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 accretion" begin
        #setup
        ρ = FT(1.1)
        q_liq = FT(0.5e-3)
        N_liq = FT(1e8)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)

        kcr = FT(5.25)
        ρ0 = FT(1.225)

        #action
        ac = CM2.accretion(SB2006, q_liq, q_rai, ρ, N_liq)

        Lc = ρ * q_liq
        Lr = ρ * q_rai
        xc = Lc / N_liq
        τ = 1 - Lc / (Lc + Lr)
        ϕ_ac = (τ / (τ + 5e-5))^4

        dqrdt_ac = kcr * Lc * Lr * ϕ_ac * sqrt(ρ0 / ρ) / ρ
        dqcdt_ac = -dqrdt_ac
        dNcdt_ac = 1 / xc * ρ * dqcdt_ac
        dNrdt_ac = FT(0)

        #test
        TT.@test ac isa CM2.LiqRaiRates
        TT.@test ac.dq_liq_dt ≈ dqcdt_ac rtol = FT(1e-6)
        TT.@test ac.dq_rai_dt ≈ dqrdt_ac rtol = FT(1e-6)
        TT.@test ac.dN_liq_dt ≈ dNcdt_ac rtol = FT(1e-6)
        TT.@test ac.dN_rai_dt ≈ dNrdt_ac rtol = FT(1e-6)

        #action
        ac = CM2.accretion(SB2006, FT(0), FT(0), ρ, N_liq)

        #test
        TT.@test ac.dq_liq_dt ≈ FT(0) atol = eps(FT)
        TT.@test ac.dq_rai_dt ≈ FT(0) atol = eps(FT)
        TT.@test ac.dN_liq_dt ≈ FT(0) atol = eps(FT)
        TT.@test ac.dN_rai_dt ≈ FT(0) atol = eps(FT)
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 rain self-collection and breakup" begin
        #setup
        ρ = FT(1.1)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)

        krr = FT(7.12)
        κrr = FT(60.7)
        Deq = FT(9e-4)
        Dr_th = FT(3.5e-4)
        kbr = FT(1000)
        κbr = FT(2300)
        ρ0 = FT(1.225)

        #action
        sc_rai =
            CM2.rain_self_collection(SB2006.pdf, SB2006.self, q_rai, ρ, N_rai)
        br_rai =
            CM2.rain_breakup(SB2006.pdf, SB2006.brek, q_rai, ρ, N_rai, sc_rai)
        sc_br_rai =
            CM2.rain_self_collection_and_breakup(SB2006, q_rai, ρ, N_rai)

        λr =
            CM2.raindrops_limited_vars(SB2006.pdf, q_rai, ρ, N_rai).λr *
            FT(6 / π / 1000)^FT(1 / 3)
        dNrdt_sc = -krr * N_rai * ρ * q_rai * (1 + κrr / λr)^-5 * sqrt(ρ0 / ρ)

        Dr =
            (
                CM2.raindrops_limited_vars(SB2006.pdf, q_rai, ρ, N_rai).xr /
                1000 / FT(π) * 6
            )^FT(1 / 3)
        ΔDr = Dr - Deq
        ϕ_br =
            Dr < 0.35e-3 ? FT(-1) :
            ((Dr < 0.9e-3) ? kbr * ΔDr : 2 * (exp(κbr * ΔDr) - 1))

        dNrdt_br = -(ϕ_br + 1) * sc_rai

        #test
        TT.@test sc_rai ≈ dNrdt_sc rtol = 1e-6
        TT.@test CM2.rain_self_collection(
            SB2006.pdf,
            SB2006.self,
            FT(0),
            ρ,
            N_rai,
        ) ≈ FT(0) atol = eps(FT)
        TT.@test br_rai ≈ dNrdt_br rtol = 1e-6
        TT.@test sc_br_rai isa NamedTuple
        TT.@test sc_br_rai.sc ≈ dNrdt_sc rtol = 1e-6
        TT.@test sc_br_rai.br ≈ dNrdt_br rtol = 1e-6

        #setup
        q_rai = FT(0)

        #action
        sc_rai =
            CM2.rain_self_collection(SB2006.pdf, SB2006.self, q_rai, ρ, N_rai)
        br_rai =
            CM2.rain_breakup(SB2006.pdf, SB2006.brek, q_rai, ρ, N_rai, sc_rai)
        sc_br_rai =
            CM2.rain_self_collection_and_breakup(SB2006, q_rai, ρ, N_rai)

        #test
        TT.@test sc_rai ≈ FT(0) atol = eps(FT)
        TT.@test br_rai ≈ FT(0) atol = eps(FT)
        TT.@test sc_br_rai.sc ≈ FT(0) atol = eps(FT)
        TT.@test sc_br_rai.br ≈ FT(0) atol = eps(FT)
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 rain terminal velocity" begin
        #setup
        ρ = FT(1.1)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)

        ρ0 = FT(1.225)
        aR = FT(9.65)
        bR = FT(10.3)
        cR = FT(600)

        #action
        vt_rai = CM2.rain_terminal_velocity(SB2006, SB2006Vel, q_rai, ρ, N_rai)

        λr = CM2.raindrops_limited_vars(SB2006.pdf, q_rai, ρ, N_rai).λr
        vt0 = max(0, sqrt(ρ0 / ρ) * (aR - bR / (1 + cR / λr)))
        vt1 = max(0, sqrt(ρ0 / ρ) * (aR - bR / (1 + cR / λr)^4))

        #test
        TT.@test vt_rai isa Tuple
        TT.@test vt_rai[1] ≈ vt0 rtol = 1e-6
        TT.@test vt_rai[2] ≈ vt1 rtol = 1e-6
        TT.@test CM2.rain_terminal_velocity(
            SB2006,
            SB2006Vel,
            FT(0),
            ρ,
            N_rai,
        )[1] ≈ 0 atol = eps(FT)
        TT.@test CM2.rain_terminal_velocity(
            SB2006,
            SB2006Vel,
            FT(0),
            ρ,
            N_rai,
        )[2] ≈ 0 atol = eps(FT)
    end

    TT.@testset "2M_microphysics - Chen 2022 rain terminal velocity" begin
        #setup
        ρ = FT(1.1)
        q_rai = FT(5e-4)
        N_rai = FT(1e4)

        #action
        vt_rai =
            CM2.rain_terminal_velocity(SB2006, Chen2022Vel, q_rai, ρ, N_rai)
        v_bigger =
            CM2.rain_terminal_velocity(SB2006, Chen2022Vel, q_rai * 2, ρ, N_rai)

        #test
        TT.@test vt_rai isa Tuple
        TT.@test vt_rai[1] ≈ 1.2591475834547752 rtol = 1e-6
        TT.@test vt_rai[2] ≈ 4.552478635185714 rtol = 1e-6

        TT.@test CM2.rain_terminal_velocity(
            SB2006,
            Chen2022Vel,
            FT(0),
            ρ,
            N_rai,
        )[1] ≈ 0 atol = eps(FT)
        TT.@test CM2.rain_terminal_velocity(
            SB2006,
            Chen2022Vel,
            FT(0),
            ρ,
            N_rai,
        )[2] ≈ 0 atol = eps(FT)

        TT.@test v_bigger[1] > vt_rai[1]
        TT.@test v_bigger[2] > vt_rai[2]
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 rain evaporation" begin
        #setup
        ρ = FT(1.1)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)
        T = FT(288.15)
        q_tot = FT(1e-3)
        q = TD.PhasePartition(q_tot)

        av = FT(0.78)
        bv = FT(0.308)
        α = FT(159.0)
        β = FT(0.266)
        (; ν_air, D_vapor) = aps
        ρ0 = FT(1.225)

        #action
        evap = CM2.rain_evaporation(SB2006, aps, tps, q, q_rai, ρ, N_rai, T)

        G = CMC.G_func(aps, tps, T, TD.Liquid())
        S = TD.supersaturation(tps, q, ρ, T, TD.Liquid())

        xr = CM2.raindrops_limited_vars(SB2006.pdf, q_rai, ρ, N_rai).xr
        Dr = FT(6 / π / 1000.0)^FT(1 / 3) * xr^FT(1 / 3)
        N_Re = α * xr^β * sqrt(ρ0 / ρ) * Dr / ν_air

        a_vent_0 = av * FT(0.1915222379058504)
        b_vent_0 = bv * FT(0.2040123897555518)
        Fv0 = a_vent_0 + b_vent_0 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)
        a_vent_1 = av * FT(0.5503212081491045)
        b_vent_1 = bv * FT(0.5873135598802672)
        Fv1 = a_vent_1 + b_vent_1 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)

        evap0 = 2 * FT(π) * G * S * N_rai * Dr * Fv0 / xr
        evap1 = 2 * FT(π) * G * S * N_rai * Dr * Fv1 / ρ

        #test
        TT.@test evap isa Tuple
        TT.@test evap[1] ≈ (evap0 - 2.5) rtol = 1e-4
        TT.@test evap[2] ≈ evap1 rtol = 1e-5
        TT.@test CM2.rain_evaporation(
            SB2006,
            aps,
            tps,
            q,
            FT(0),
            ρ,
            N_rai,
            T,
        )[1] ≈ 0 atol = eps(FT)
        TT.@test CM2.rain_evaporation(
            SB2006,
            aps,
            tps,
            q,
            FT(0),
            ρ,
            N_rai,
            T,
        )[2] ≈ 0 atol = eps(FT)
    end
end

println("Testing Float64")
test_microphysics2M(Float64)

println("Testing Float32")
test_microphysics2M(Float32)
