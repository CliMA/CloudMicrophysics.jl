import Test as TT

import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaParams as CP

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CMC
import CloudMicrophysics.Microphysics2M as CM2

import QuadGK as QGK
import SpecialFunctions as SF

@info "Microphysics 2M Tests"

function test_microphysics2M(FT)

    # Different 2-moment autoconversion and accretion parameters
    KK2000 = CMP.KK2000(FT)
    B1994 = CMP.B1994(FT)
    TC1980 = CMP.TC1980(FT)
    LD2004 = CMP.LD2004(FT)
    VarTSc = CMP.VarTimescaleAcnv(FT)

    # Seifert and Beheng 2006 parameters
    override_file = joinpath(
        pkgdir(CM),
        "src",
        "parameters",
        "toml",
        "SB2006_limiters.toml",
    )
    toml_dict = CP.create_toml_dict(FT; override_file)
    SB2006 = CMP.SB2006(toml_dict)
    SB2006_no_limiters = CMP.SB2006(toml_dict, false)

    # Thermodynamics and air properties parameters
    aps = CMP.AirProperties(FT)
    tps = TD.Parameters.ThermodynamicsParameters(FT)

    # Terminal velocity parameters
    SB2006Vel = CMP.SB2006VelType(FT)
    Chen2022Vel = CMP.Chen2022VelTypeRain(FT)

    TT.@testset "2M_microphysics - unit tests" begin

        ρ = FT(1)

        # no reference data available - checking if callable and not NaN
        q_liq = FT(0.5e-3)
        q_rai = FT(1e-6)
        N_d = FT(1e8)

        TT.@test CM2.accretion(KK2000, q_liq, q_rai, ρ) != NaN
        TT.@test CM2.accretion(B1994, q_liq, q_rai, ρ) != NaN
        TT.@test CM2.accretion(TC1980, q_liq, q_rai) != NaN
        TT.@test CM2.conv_q_liq_to_q_rai(VarTSc, q_liq, ρ, N_d) != NaN

        # output should be zero if either q_liq or q_rai are zero
        q_liq = FT(0)
        q_rai = FT(1e-6)

        TT.@test CM2.conv_q_liq_to_q_rai(VarTSc, q_liq, ρ, N_d) == FT(0)
        TT.@test CM2.conv_q_liq_to_q_rai(KK2000, q_liq, ρ, N_d) == FT(0)
        TT.@test CM2.conv_q_liq_to_q_rai(B1994, q_liq, ρ, N_d) == FT(0)
        TT.@test CM2.conv_q_liq_to_q_rai(TC1980, q_liq, ρ, N_d) == FT(0)
        TT.@test CM2.conv_q_liq_to_q_rai(LD2004, q_liq, ρ, N_d) == FT(0)
        TT.@test CM2.accretion(KK2000, q_liq, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(B1994, q_liq, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(TC1980, q_liq, q_rai) == FT(0)

        q_liq = FT(0.5e-3)
        q_rai = FT(0)
        TT.@test CM2.accretion(KK2000, q_liq, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(B1994, q_liq, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(TC1980, q_liq, q_rai) == FT(0)

        TT.@test CM2.conv_q_liq_to_q_rai(VarTSc, q_liq, ρ, N_d) >
                 CM2.conv_q_liq_to_q_rai(VarTSc, q_liq, ρ, 10 * N_d)

        # far from threshold points, autoconversion with and without smooth transition should
        # be approximately equal
        q_liq = FT(0.5e-3)
        TT.@test CM2.conv_q_liq_to_q_rai(B1994, q_liq, ρ, N_d, true) ≈
                 CM2.conv_q_liq_to_q_rai(B1994, q_liq, ρ, N_d, false) rtol = 0.2
        TT.@test CM2.conv_q_liq_to_q_rai(TC1980, q_liq, ρ, N_d, true) ≈
                 CM2.conv_q_liq_to_q_rai(TC1980, q_liq, ρ, N_d, false) rtol =
            0.2
        TT.@test CM2.conv_q_liq_to_q_rai(LD2004, q_liq, ρ, N_d, true) ≈
                 CM2.conv_q_liq_to_q_rai(LD2004, q_liq, ρ, N_d, false) rtol =
            0.2

    end

    TT.@testset "2M_microphysics - compare with Wood_2005" begin

        ρ = FT(1)
        q_liq = FT(0.5e-3)
        N_d = FT(1e8)

        # compare with Wood 2005 Fig 1 panel a
        function compare(scheme, input, output; eps = 0.1)
            TT.@test CM2.conv_q_liq_to_q_rai(scheme, input * FT(1e-3), ρ, N_d) ≈
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
                ρ,
                input * FT(1e6),
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
        q_rai = [FT(0), FT(1e-3), FT(1e-4), FT(1e-2)]
        N_rai = [FT(1e1), FT(1e1), FT(1e3), FT(1e5)]
        ρ = FT(1)

        (; xr_min, xr_max, λ_min, λ_max) = SB2006.pdf_r

        for Nr in N_rai
            for qr in q_rai
                #action
                λ = CM2.pdf_rain(SB2006.pdf_r, qr, ρ, Nr).λr
                xr = CM2.pdf_rain(SB2006.pdf_r, qr, ρ, Nr).xr

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

        for SB in [SB2006, SB2006_no_limiters]
            (; kcc, x_star, ρ0) = SB.acnv
            (; νc) = SB.pdf_c

            #action
            au = CM2.autoconversion(SB.acnv, SB.pdf_c, q_liq, q_rai, ρ, N_liq)
            sc = CM2.liquid_self_collection(
                SB.acnv,
                SB.pdf_c,
                q_liq,
                ρ,
                au.dN_liq_dt,
            )
            au_sc = CM2.autoconversion_and_liquid_self_collection(
                SB,
                q_liq,
                q_rai,
                ρ,
                N_liq,
            )

            Lc = ρ * q_liq
            Lr = ρ * q_rai
            xc = min(x_star, Lc / N_liq)
            τ = 1 - Lc / (Lc + Lr)
            ϕ_au = 400 * τ^0.7 * (1 - τ^0.7)^3
            dqrdt_au =
                kcc / 20 / x_star * (νc + 2) * (νc + 4) / (νc + 1)^2 *
                Lc^2 *
                xc^2 *
                (1 + ϕ_au / (1 - τ)^2) *
                (ρ0 / ρ) / ρ
            dqcdt_au = -dqrdt_au
            dNcdt_au = 2 / x_star * ρ * dqcdt_au
            dNrdt_au = -0.5 * dNcdt_au
            dNcdt_sc =
                -kcc * (νc + 2) / (νc + 1) * (ρ0 / ρ) * Lc^2 - au.dN_liq_dt

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
            au = CM2.autoconversion(SB.acnv, SB.pdf_c, FT(0), FT(0), ρ, N_liq)
            sc = CM2.liquid_self_collection(
                SB.acnv,
                SB.pdf_c,
                FT(0),
                ρ,
                au.dN_liq_dt,
            )
            au_sc = CM2.autoconversion_and_liquid_self_collection(
                SB,
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
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 accretion" begin
        #setup
        ρ = FT(1.1)
        q_liq = FT(0.5e-3)
        N_liq = FT(1e8)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)

        for SB in [SB2006, SB2006_no_limiters]
            (; kcr, ρ0) = SB.accr

            #action
            ac = CM2.accretion(SB, q_liq, q_rai, ρ, N_liq)

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
            ac = CM2.accretion(SB, FT(0), FT(0), ρ, N_liq)

            #test
            TT.@test ac.dq_liq_dt ≈ FT(0) atol = eps(FT)
            TT.@test ac.dq_rai_dt ≈ FT(0) atol = eps(FT)
            TT.@test ac.dN_liq_dt ≈ FT(0) atol = eps(FT)
            TT.@test ac.dN_rai_dt ≈ FT(0) atol = eps(FT)
        end
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 rain self-collection and breakup" begin
        #setup
        ρ = FT(1.1)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)

        for SB in [SB2006, SB2006_no_limiters]
            (; krr, κrr) = SB.self
            (; Deq, kbr, κbr) = SB.brek
            ρ0 = SB.pdf_r.ρ0

            #action
            sc_rai =
                CM2.rain_self_collection(SB.pdf_r, SB.self, q_rai, ρ, N_rai)
            br_rai =
                CM2.rain_breakup(SB.pdf_r, SB.brek, q_rai, ρ, N_rai, sc_rai)
            sc_br_rai =
                CM2.rain_self_collection_and_breakup(SB, q_rai, ρ, N_rai)

            λr = CM2.pdf_rain(SB.pdf_r, q_rai, ρ, N_rai).Br

            dNrdt_sc =
                -krr * N_rai * ρ * q_rai * (1 + κrr / λr)^-5 * sqrt(ρ0 / ρ)

            Dr =
                (
                    CM2.pdf_rain(SB.pdf_r, q_rai, ρ, N_rai).xr / 1000 / FT(π) *
                    6
                )^FT(1 / 3)
            ΔDr = Dr - Deq
            ϕ_br =
                Dr < 0.35e-3 ? FT(-1) :
                ((Dr < 0.9e-3) ? kbr * ΔDr : 2 * (exp(κbr * ΔDr) - 1))

            dNrdt_br = -(ϕ_br + 1) * sc_rai

            #test
            TT.@test sc_rai ≈ dNrdt_sc rtol = 1e-6
            TT.@test CM2.rain_self_collection(
                SB.pdf_r,
                SB.self,
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
                CM2.rain_self_collection(SB.pdf_r, SB.self, q_rai, ρ, N_rai)
            br_rai =
                CM2.rain_breakup(SB.pdf_r, SB.brek, q_rai, ρ, N_rai, sc_rai)
            sc_br_rai =
                CM2.rain_self_collection_and_breakup(SB, q_rai, ρ, N_rai)

            #test
            TT.@test sc_rai ≈ FT(0) atol = eps(FT)
            TT.@test br_rai ≈ FT(0) atol = eps(FT)
            TT.@test sc_br_rai.sc ≈ FT(0) atol = eps(FT)
            TT.@test sc_br_rai.br ≈ FT(0) atol = eps(FT)
        end
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 rain terminal velocity with limiters" begin
        #setup
        ρ = FT(1.1)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)

        (; ρ0, aR, bR, cR) = SB2006Vel

        #action
        vt_rai = CM2.rain_terminal_velocity(SB2006, SB2006Vel, q_rai, ρ, N_rai)

        λr = CM2.pdf_rain(SB2006.pdf_r, q_rai, ρ, N_rai).λr
        vt0 = max(0, sqrt(ρ0 / ρ) * (aR - bR / (1 + cR / λr)))
        vt1 = max(0, sqrt(ρ0 / ρ) * (aR - bR / (1 + cR / λr)^4))

        #test
        TT.@test vt_rai isa Tuple
        TT.@test vt_rai[1] ≈ vt0 rtol = 1e-6
        TT.@test vt_rai[2] ≈ vt1 rtol = 1e-6

        TT.@test CM2.rain_terminal_velocity(
            SB2006,
            SB2006Vel,
            q_rai,
            ρ,
            FT(0),
        )[1] ≈ 0 atol = eps(FT)
        TT.@test CM2.rain_terminal_velocity(
            SB2006,
            SB2006Vel,
            FT(0),
            ρ,
            N_rai,
        )[2] ≈ 0 atol = eps(FT)
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 modified rain terminal velocity without limiters" begin
        #setup
        ρ = FT(1.1)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)

        (; ρ0, aR, bR, cR) = SB2006Vel

        #action
        vt_rai = CM2.rain_terminal_velocity(
            SB2006_no_limiters,
            SB2006Vel,
            q_rai,
            ρ,
            N_rai,
        )

        λr = CM2.pdf_rain(SB2006_no_limiters.pdf_r, q_rai, ρ, N_rai).λr
        _rc = -1 / (2 * cR) * log(aR / bR)
        _Γ_1(t) = exp(-t)
        _Γ_4(t) = (t^3 + 3 * t^2 + 6 * t + 6) * exp(-t)
        _pa0 = _Γ_1(2 * _rc * λr)
        _pb0 = _Γ_1(2 * _rc * (λr + cR))
        _pa1 = _Γ_4(2 * _rc * λr) / FT(6)
        _pb1 = _Γ_4(2 * _rc * (λr + cR)) / FT(6)
        vt0 = max(0, sqrt(ρ0 / ρ) * (aR * _pa0 - bR * _pb0 / (1 + cR / λr)))
        vt1 = max(0, sqrt(ρ0 / ρ) * (aR * _pa1 - bR * _pb1 / (1 + cR / λr)^4))

        #test
        TT.@test vt_rai isa Tuple
        TT.@test vt_rai[1] ≈ vt0 rtol = 1e-6
        TT.@test vt_rai[2] ≈ vt1 rtol = 1e-6

        TT.@test CM2.rain_terminal_velocity(
            SB2006_no_limiters,
            SB2006Vel,
            q_rai,
            ρ,
            FT(0),
        )[1] ≈ 0 atol = eps(FT)
        TT.@test CM2.rain_terminal_velocity(
            SB2006_no_limiters,
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

        for SB in [SB2006, SB2006_no_limiters]
            #action
            vt_rai =
                CM2.rain_terminal_velocity(SB, Chen2022Vel, q_rai, ρ, N_rai)
            v_bigger =
                CM2.rain_terminal_velocity(SB, Chen2022Vel, q_rai * 2, ρ, N_rai)

            #test
            TT.@test vt_rai isa Tuple
            TT.@test vt_rai[1] ≈ 1.0738503635546666 rtol = 1e-6
            TT.@test vt_rai[2] ≈ 4.00592218028957 rtol = 1e-6

            TT.@test CM2.rain_terminal_velocity(
                SB,
                Chen2022Vel,
                q_rai,
                ρ,
                FT(0),
            )[1] ≈ 0 atol = eps(FT)
            TT.@test CM2.rain_terminal_velocity(
                SB,
                Chen2022Vel,
                FT(0),
                ρ,
                N_rai,
            )[2] ≈ 0 atol = eps(FT)

            TT.@test v_bigger[1] > vt_rai[1]
            TT.@test v_bigger[2] > vt_rai[2]
        end
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 rain evaporation" begin
        #setup
        ρ = FT(1.1)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)
        T = FT(288.15)
        q_tot = FT(1e-3)
        q = TD.PhasePartition(q_tot)

        for SB in [SB2006, SB2006_no_limiters]

            (; av, bv, α, β, ρ0) = SB.evap
            (; ν_air, D_vapor) = aps

            #action
            evap = CM2.rain_evaporation(SB, aps, tps, q, q_rai, ρ, N_rai, T)

            G = CMC.G_func(aps, tps, T, TD.Liquid())
            S = TD.supersaturation(tps, q, ρ, T, TD.Liquid())

            xr = CM2.pdf_rain(SB.pdf_r, q_rai, ρ, N_rai).xr
            Dr = FT(6 / π / 1000.0)^FT(1 / 3) * xr^FT(1 / 3)
            N_Re = α * xr^β * sqrt(ρ0 / ρ) * Dr / ν_air

            a_vent_0 = av * FT(0.15344374450453543)
            b_vent_0 = bv * FT(0.17380986321413017)
            Fv0 = a_vent_0 + b_vent_0 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)
            a_vent_1 = av * FT(0.5503212081491045)
            b_vent_1 = bv * FT(0.5873135598802672)
            Fv1 = a_vent_1 + b_vent_1 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)

            evap0 = 2 * FT(π) * G * S * N_rai * Dr * Fv0 / xr
            evap1 = 2 * FT(π) * G * S * N_rai * Dr * Fv1 / ρ

            #test
            TT.@test evap isa NamedTuple
            TT.@test evap.evap_rate_0 ≈ evap0 rtol = 1e-4
            TT.@test evap.evap_rate_1 ≈ evap1 rtol = 1e-5
            TT.@test CM2.rain_evaporation(
                SB,
                aps,
                tps,
                q,
                q_rai,
                ρ,
                FT(0),
                T,
            ).evap_rate_0 ≈ 0 atol = eps(FT)
            TT.@test CM2.rain_evaporation(
                SB,
                aps,
                tps,
                q,
                FT(0),
                ρ,
                N_rai,
                T,
            ).evap_rate_1 ≈ 0 atol = eps(FT)
        end

        # test limit case: xr = 0 for SB with no limiters
        TT.@test CM2.rain_evaporation(
            SB2006_no_limiters,
            aps,
            tps,
            q,
            FT(0),
            ρ,
            N_rai,
            T,
        ).evap_rate_0 ≈ 0 atol = eps(FT)

    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 effective radius and reflectivity" begin
        #setup
        ρₐ = FT(1)

        q_liq = [FT(2.128e-4), FT(2.128e-20), FT(1.6e-12), FT(0), FT(1.037e-25)]
        N_liq = [FT(15053529), FT(3), FT(5512), FT(0), FT(5.225e-12)]
        q_rai = [FT(1.573e-4), FT(1.573e-4), FT(1.9e-15), FT(0), FT(2.448e-27)]
        N_rai = [FT(510859), FT(510859), FT(0), FT(0), FT(5.136e-18)]

        # reference values
        rr = [FT(-12.561951), FT(-12.579899), FT(-150), FT(-150), FT(-150)]
        reff = [FT(2.319383e-5), FT(6.91594e-5), FT(0), FT(0), FT(0)]

        for (qₗ, Nₗ, qᵣ, Nᵣ, rₑ, Z) in zip(q_liq, N_liq, q_rai, N_rai, reff, rr)
            for SB in [SB2006, SB2006_no_limiters]

                #action
                Z_val = CM2.radar_reflectivity(SB, qₗ, qᵣ, Nₗ, Nᵣ, ρₐ)
                rₑ_val = CM2.effective_radius(SB, qₗ, qᵣ, Nₗ, Nᵣ, ρₐ)

                #test
                TT.@test rₑ_val ≈ rₑ atol = FT(1e-6)
                TT.@test Z_val ≈ Z atol = FT(1e-4)
            end
        end
    end

    TT.@testset "2M_microphysics - '1/3' power law from Liu and Hallett (1997)" begin
        #setup
        ρ_air = FT(1)
        ρ_w = FT(1000)
        q_liq = FT(2.128e-4)
        N_liq = FT(15053529)
        q_rai = FT(1.573e-4)
        N_rai = FT(510859)

        #action
        reff = CM2.effective_radius_Liu_Hallet_97(
            q_liq,
            q_rai,
            N_liq,
            N_rai,
            ρ_air,
            ρ_w,
        )

        #test
        TT.@test reff ≈ FT(2.66e-05) atol = FT(8e-6)

    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 rain distribution sanity checks" begin

        # air and liquid water densities
        ρₐ = FT(1.2)
        ρₗ = SB2006.pdf_r.ρw

        # example number concentration and specific humidity
        Nᵣ = FT(5e5)
        qᵣ = FT(5e-4)

        # limited distribution parameters for rain
        Ar_l = CM2.pdf_rain(SB2006.pdf_r, qᵣ, ρₐ, Nᵣ).Ar
        Br_l = CM2.pdf_rain(SB2006.pdf_r, qᵣ, ρₐ, Nᵣ).Br
        λr_l = CM2.pdf_rain(SB2006.pdf_r, qᵣ, ρₐ, Nᵣ).λr
        αr_l = CM2.pdf_rain(SB2006.pdf_r, qᵣ, ρₐ, Nᵣ).αr
        # not limited distribution parameters for rain
        (; λr, αr, Ar, Br) = CM2.pdf_rain(SB2006_no_limiters.pdf_r, qᵣ, ρₐ, Nᵣ)

        # mass of liquid droplet as a function of its diameter
        m(D) = FT(π / 6) * ρₗ * D^3

        # rain drop diameter distribution (eq.(3) from 2M docs)
        f_D(D) = αr * exp(-λr * D)
        # rain drop diameter distribution, but using SB2006 limiters
        f_D_limited(D) = αr_l * exp(-λr_l * D)

        # rain drop mass distribution (eq.(4) from 2M docs)
        f_x(x) = Ar * x^(-2 / 3) * exp(-Br * x^(1 / 3))
        # rain drop mass distribution, but using the SB2006 limiters
        f_x_limited(x) = Ar_l * x^(-2 / 3) * exp(-Br_l * x^(1 / 3))

        # integral bounds
        D₀ = 1e-7
        D∞ = 1e-2
        m₀ = m(D₀)
        m∞ = m(D∞)

        # Sanity checks for number concentrations for rain
        ND = QGK.quadgk(x -> f_D(x), D₀, D∞)[1]
        Nx = QGK.quadgk(x -> f_x(x), m₀, m∞)[1]
        ND_lim = QGK.quadgk(x -> f_D_limited(x), D₀, D∞)[1]
        Nx_lim = QGK.quadgk(x -> f_x_limited(x), m₀, m∞)[1]
        TT.@test ND ≈ Nᵣ rtol = FT(1e-2)
        TT.@test Nx ≈ Nᵣ rtol = FT(1e-2)
        TT.@test ND_lim ≈ Nᵣ rtol = FT(1e-2)
        TT.@test Nx_lim ≈ Nᵣ rtol = FT(1e-2)

        # Sanity checks for specific humidities for rain
        qD = QGK.quadgk(x -> m(x) * f_D(x), D₀, D∞)[1] / ρₐ
        qx = QGK.quadgk(x -> x * f_x(x), m₀, m∞)[1] / ρₐ
        qD_lim = QGK.quadgk(x -> m(x) * f_D_limited(x), D₀, D∞)[1] / ρₐ
        qx_lim = QGK.quadgk(x -> x * f_x_limited(x), m₀, m∞)[1] / ρₐ
        TT.@test qD ≈ qᵣ atol = FT(1e-6)
        TT.@test qx ≈ qᵣ atol = FT(1e-6)

        # The mass integrals don't work with limiters
        TT.@test qx_lim ≈ qᵣ atol = FT(1e-6)
        TT.@test qD_lim ≈ qx_lim atol = FT(1e-6)
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 cloud distribution sanity checks" begin

        # example number concentration and specific humidity
        Nₗ = FT(1e9)
        qₗ = FT(1e-3)

        # air and liquid water densities in μg/m3
        ρₐ = FT(1.2)
        ρₗ = SB2006.pdf_r.ρw

        # distribution parameters for cloud, units: [Bc] = 1/μg, [Ac] = 1/m3 1/μg3
        (; χ, Ac, Bc, Cc, Ec, ϕc, ψc) =
            CM2.pdf_cloud(SB2006_no_limiters.pdf_c, qₗ, ρₐ, Nₗ)

        # mass of liquid droplet as a function of its diameter in μg
        m(D) = FT(π / 6) * ρₗ * FT(10)^χ * D^3

        # cloud droplet mass distribution (eq.(2) from 2M docs)
        f_x(x) = Ac * x^(2) * exp(-Bc * x^(1))

        # cloud droplet diameter distribution (eq.(7) from 2M docs)
        f_D(D) = Cc * D^ϕc * exp(-Ec * D^ψc)

        # integral bounds
        D₀ = 1e-8
        D∞ = 1e-4
        m₀ = m(D₀)
        m∞ = m(D∞)

        # Sanity checks specific humidity and number concentration with mass distribution
        # Sanity checks for number concentrations for cloud
        qx = QGK.quadgk(x -> x * f_x(x), m₀, m∞)[1] / (ρₐ * FT(10)^χ)
        Nx = QGK.quadgk(x -> f_x(x), m₀, m∞)[1]
        TT.@test qx ≈ qₗ rtol = FT(1e-6)
        TT.@test Nx ≈ Nₗ rtol = FT(1e-6)

        # mass of liquid droplets as a function of its diameter in mm and μg
        _m(D) = FT(π / 6) * ρₗ * FT(10)^(χ - 9) * D^3

        # integral bounds in millimiters
        _D₀ = 1e-8 * FT(1e3)
        _D∞ = 1e-4 * FT(1e3)

        # Sanity checks specific humidity and number concentration with diameter distribution
        qD =
            QGK.quadgk(x -> _m(x) * f_D(x), _D₀, _D∞)[1] / (ρₐ * FT(10)^(χ - 9))
        ND = QGK.quadgk(x -> f_D(x), _D₀, _D∞)[1] * FT(1e9)
        TT.@test ND ≈ Nₗ rtol = FT(1e-6)
        TT.@test qD ≈ qₗ rtol = FT(1e-6)

    end
end

println("Testing Float64")
test_microphysics2M(Float64)

println("Testing Float32")
test_microphysics2M(Float32)
