import Test as TT

import ClimaParams as CP

import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CMC
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.DistributionTools as DT

import QuadGK as QGK
import SpecialFunctions as SF

function test_microphysics2M(FT)

    # Different 2-moment autoconversion and accretion parameters
    KK2000 = CMP.KK2000(FT)
    B1994 = CMP.B1994(FT)
    TC1980 = CMP.TC1980(FT)
    LD2004 = CMP.LD2004(FT)
    VarTSc = CMP.VarTimescaleAcnv(FT)

    # Seifert and Beheng 2006 parameters
    override_file = joinpath(
        pkgdir(CM), "src", "parameters", "toml", "SB2006_limiters.toml",
    )
    toml_dict = CP.create_toml_dict(FT; override_file)
    SB2006 = CMP.SB2006(toml_dict)
    SB2006_no_limiters = CMP.SB2006(toml_dict, false)

    # Thermodynamics and air properties parameters
    aps = CMP.AirProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    # Terminal velocity parameters
    STVel = CMP.StokesRegimeVelType(FT)
    SB2006Vel = CMP.SB2006VelType(FT)
    Chen2022Vel = CMP.Chen2022VelTypeRain(FT)

    TT.@testset "2M_microphysics - unit tests" begin

        ρ = FT(1)

        # no reference data available - checking if callable and not NaN
        q_lcl = FT(0.5e-3)
        q_rai = FT(1e-6)
        N_d = FT(1e8)

        TT.@test CM2.accretion(KK2000, q_lcl, q_rai, ρ) != NaN
        TT.@test CM2.accretion(B1994, q_lcl, q_rai, ρ) != NaN
        TT.@test CM2.accretion(TC1980, q_lcl, q_rai) != NaN
        TT.@test CM2.conv_q_lcl_to_q_rai(VarTSc, q_lcl, ρ, N_d) != NaN

        # output should be zero if either q_liq or q_rai are zero
        q_lcl = FT(0)
        q_rai = FT(1e-6)

        TT.@test CM2.conv_q_lcl_to_q_rai(VarTSc, q_lcl, ρ, N_d) == FT(0)
        TT.@test CM2.conv_q_lcl_to_q_rai(KK2000, q_lcl, ρ, N_d) == FT(0)
        TT.@test CM2.conv_q_lcl_to_q_rai(B1994, q_lcl, ρ, N_d) == FT(0)
        TT.@test CM2.conv_q_lcl_to_q_rai(TC1980, q_lcl, ρ, N_d) == FT(0)
        TT.@test CM2.conv_q_lcl_to_q_rai(LD2004, q_lcl, ρ, N_d) == FT(0)
        TT.@test CM2.accretion(KK2000, q_lcl, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(B1994, q_lcl, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(TC1980, q_lcl, q_rai) == FT(0)

        q_lcl = FT(0.5e-3)
        q_rai = FT(0)
        TT.@test CM2.accretion(KK2000, q_lcl, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(B1994, q_lcl, q_rai, ρ) == FT(0)
        TT.@test CM2.accretion(TC1980, q_lcl, q_rai) == FT(0)

        TT.@test CM2.conv_q_lcl_to_q_rai(VarTSc, q_lcl, ρ, N_d) >
                 CM2.conv_q_lcl_to_q_rai(VarTSc, q_lcl, ρ, 10 * N_d)

        # far from threshold points, autoconversion with and without smooth transition should
        # be approximately equal
        q_lcl = FT(0.5e-3)
        TT.@test CM2.conv_q_lcl_to_q_rai(B1994, q_lcl, ρ, N_d, true) ≈
                 CM2.conv_q_lcl_to_q_rai(B1994, q_lcl, ρ, N_d, false) rtol = 0.2
        TT.@test CM2.conv_q_lcl_to_q_rai(TC1980, q_lcl, ρ, N_d, true) ≈
                 CM2.conv_q_lcl_to_q_rai(TC1980, q_lcl, ρ, N_d, false) rtol =
            0.2
        TT.@test CM2.conv_q_lcl_to_q_rai(LD2004, q_lcl, ρ, N_d, true) ≈
                 CM2.conv_q_lcl_to_q_rai(LD2004, q_lcl, ρ, N_d, false) rtol =
            0.2

    end

    TT.@testset "2M_microphysics - compare with Wood_2005" begin

        ρ = FT(1)
        q_lcl = FT(0.5e-3)
        N_d = FT(1e8)

        # compare with Wood 2005 Fig 1 panel a
        function compare(scheme, input, output; eps = 0.1)
            TT.@test CM2.conv_q_lcl_to_q_rai(scheme, input * FT(1e-3), ρ, N_d) ≈
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
            TT.@test CM2.conv_q_lcl_to_q_rai(
                scheme,
                q_lcl,
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
    TT.@testset "Seifert and Beheng 2006 - PDF parameters limiting behavior" begin
        N = 0.0
        q = 0.0
        ρₐ = 1.2
        # limited rain drop size distribution
        params = CM2.pdf_rain_parameters(SB2006.pdf_r, q, ρₐ, N)
        TT.@test all(iszero, params)
        n = CM2.size_distribution(SB2006.pdf_r, q, ρₐ, N)
        TT.@test all(iszero, (n(0), n(0.1), n(Inf)))
        bnds = CM2.get_size_distribution_bounds(SB2006.pdf_r, q, ρₐ, N)
        TT.@test all(iszero, bnds)
        # not limited rain drop size distribution
        params = CM2.pdf_rain_parameters(SB2006_no_limiters.pdf_r, q, ρₐ, N)
        TT.@test all(iszero, params)
        n = CM2.size_distribution(SB2006_no_limiters.pdf_r, q, ρₐ, N)
        TT.@test all(iszero, (n(0), n(0.1), n(Inf)))
        bnds = CM2.get_size_distribution_bounds(SB2006_no_limiters.pdf_r, q, ρₐ, N)
        TT.@test all(iszero, bnds)
        # cloud drop size distribution
        logA, logB = CM2.log_pdf_cloud_parameters_mass(SB2006.pdf_c, q, ρₐ, N)
        TT.@test logA == -Inf
        TT.@test logB == Inf
        A, B = CM2.pdf_cloud_parameters_mass(SB2006.pdf_c, q, ρₐ, N)
        TT.@test A == 0
        TT.@test B == Inf
        n = CM2.size_distribution(SB2006.pdf_c, q, ρₐ, N)
        TT.@test all(iszero, (n(0), n(0.1), n(Inf)))
    end
    TT.@testset "limiting lambda_r and x_r - Seifert and Beheng 2006" begin
        #setup
        q_rai = [FT(0), FT(1e-3), FT(1e-4), FT(1e-2)]
        N_rai = [FT(1e1), FT(1e1), FT(1e3), FT(1e5)]
        ρ = FT(1)

        (; xr_min, xr_max, λ_min, λ_max) = SB2006.pdf_r

        for Nr in N_rai
            for qr in q_rai
                #action
                (; Dr_mean, xr_mean) = CM2.pdf_rain_parameters(SB2006.pdf_r, qr, ρ, Nr)
                λ = 1 / Dr_mean

                # Test limits, with tolerance 1e-5
                tol = eps(λ)
                TT.@test λ_min - tol <= λ <= λ_max + tol
                TT.@test xr_min - tol <= xr_mean <= xr_max + tol
            end
        end

    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 autoconversion and cloud liquid self-collection" begin
        #setup
        ρ = FT(1)
        q_lcl = FT(0.5e-3)
        N_lcl = FT(1e8)
        q_rai = FT(1e-6)

        for SB in [SB2006, SB2006_no_limiters]
            (; kcc, x_star, ρ0) = SB.acnv
            (; νc) = SB.pdf_c

            #action
            au = CM2.autoconversion(SB.acnv, SB.pdf_c, q_lcl, q_rai, ρ, N_lcl)
            sc = CM2.cloud_liquid_self_collection(
                SB.acnv,
                SB.pdf_c,
                q_lcl,
                ρ,
                au.dN_lcl_dt,
            )
            au_sc = CM2.autoconversion_and_cloud_liquid_self_collection(
                SB,
                q_lcl,
                q_rai,
                ρ,
                N_lcl,
            )

            Lc = ρ * q_lcl
            Lr = ρ * q_rai
            xc = min(x_star, Lc / N_lcl)
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
                -kcc * (νc + 2) / (νc + 1) * (ρ0 / ρ) * Lc^2 - au.dN_lcl_dt

            #test
            TT.@test au isa CM2.LclRaiRates
            TT.@test au.dq_lcl_dt ≈ dqcdt_au rtol = 1e-6
            TT.@test au.dq_rai_dt ≈ dqrdt_au rtol = 1e-6
            TT.@test au.dN_lcl_dt ≈ dNcdt_au rtol = 1e-6
            TT.@test au.dN_rai_dt ≈ dNrdt_au rtol = 1e-6
            TT.@test sc ≈ dNcdt_sc rtol = 1e-6
            TT.@test au_sc isa NamedTuple
            TT.@test au_sc.au.dq_lcl_dt ≈ dqcdt_au rtol = 1e-6
            TT.@test au_sc.au.dq_rai_dt ≈ dqrdt_au rtol = 1e-6
            TT.@test au_sc.au.dN_lcl_dt ≈ dNcdt_au rtol = 1e-6
            TT.@test au_sc.au.dN_rai_dt ≈ dNrdt_au rtol = 1e-6
            TT.@test au_sc.sc ≈ dNcdt_sc rtol = 1e-6

            #action
            au = CM2.autoconversion(SB.acnv, SB.pdf_c, FT(0), FT(0), ρ, N_lcl)
            sc = CM2.cloud_liquid_self_collection(
                SB.acnv,
                SB.pdf_c,
                FT(0),
                ρ,
                au.dN_lcl_dt,
            )
            au_sc = CM2.autoconversion_and_cloud_liquid_self_collection(
                SB,
                FT(0),
                FT(0),
                ρ,
                N_lcl,
            )

            #test
            TT.@test au.dq_lcl_dt ≈ FT(0) atol = eps(FT)
            TT.@test au.dq_rai_dt ≈ FT(0) atol = eps(FT)
            TT.@test au.dN_lcl_dt ≈ FT(0) atol = eps(FT)
            TT.@test au.dN_rai_dt ≈ FT(0) atol = eps(FT)
            TT.@test sc ≈ FT(0) atol = eps(FT)
            TT.@test au_sc.au.dq_lcl_dt ≈ FT(0) atol = eps(FT)
            TT.@test au_sc.au.dq_rai_dt ≈ FT(0) atol = eps(FT)
            TT.@test au_sc.au.dN_lcl_dt ≈ FT(0) atol = eps(FT)
            TT.@test au_sc.au.dN_rai_dt ≈ FT(0) atol = eps(FT)
            TT.@test au_sc.sc ≈ FT(0) atol = eps(FT)
        end
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 accretion" begin
        #setup
        ρ = FT(1.1)
        q_lcl = FT(0.5e-3)
        N_lcl = FT(1e8)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)

        for SB in [SB2006, SB2006_no_limiters]
            (; kcr, ρ0) = SB.accr

            #action
            ac = CM2.accretion(SB, q_lcl, q_rai, ρ, N_lcl)

            Lc = ρ * q_lcl
            Lr = ρ * q_rai
            xc = Lc / N_lcl
            τ = 1 - Lc / (Lc + Lr)
            ϕ_ac = (τ / (τ + 5e-5))^4

            dqrdt_ac = kcr * Lc * Lr * ϕ_ac * sqrt(ρ0 / ρ) / ρ
            dqcdt_ac = -dqrdt_ac
            dNcdt_ac = 1 / xc * ρ * dqcdt_ac
            dNrdt_ac = FT(0)

            #test
            TT.@test ac isa CM2.LclRaiRates
            TT.@test ac.dq_lcl_dt ≈ dqcdt_ac rtol = FT(1e-6)
            TT.@test ac.dq_rai_dt ≈ dqrdt_ac rtol = FT(1e-6)
            TT.@test ac.dN_lcl_dt ≈ dNcdt_ac rtol = FT(1e-6)
            TT.@test ac.dN_rai_dt ≈ dNrdt_ac rtol = FT(1e-6)

            #action
            ac = CM2.accretion(SB, FT(0), FT(0), ρ, N_lcl)

            #test
            TT.@test ac.dq_lcl_dt ≈ FT(0) atol = eps(FT)
            TT.@test ac.dq_rai_dt ≈ FT(0) atol = eps(FT)
            TT.@test ac.dN_lcl_dt ≈ FT(0) atol = eps(FT)
            TT.@test ac.dN_rai_dt ≈ FT(0) atol = eps(FT)
        end
    end

    for SB in [SB2006, SB2006_no_limiters]
        sb_str = CMP.islimited(SB.pdf_r) ? "with limiters" : "without limiters"
        TT.@testset "2M_microphysics - Seifert and Beheng 2006 rain self-collection and breakup ($sb_str)" begin
            # Setup
            ρ = FT(1.1)
            q_rai = FT(1e-6)
            N_rai = FT(1e4)

            (; krr, κrr) = SB.self
            (; Deq, kbr, κbr) = SB.brek
            ρ0 = SB.pdf_r.ρ0

            # Action
            sc_rai = CM2.rain_self_collection(SB.pdf_r, SB.self, q_rai, ρ, N_rai)
            br_rai = CM2.rain_breakup(SB.pdf_r, SB.brek, q_rai, ρ, N_rai, sc_rai)
            sc_br_rai = CM2.rain_self_collection_and_breakup(SB, q_rai, ρ, N_rai)

            (; xr_mean) = CM2.pdf_rain_parameters(SB.pdf_r, q_rai, ρ, N_rai)
            (; Br) = CM2.pdf_rain_parameters_mass(SB.pdf_r, q_rai, ρ, N_rai)

            dNrdt_sc = -krr * N_rai * ρ * q_rai * (1 + κrr / Br)^-5 * √(ρ0 / ρ)
            Dr = cbrt(xr_mean / 1000 / FT(π) * 6)
            ΔDr = Dr - Deq
            ϕ_br =
                Dr < 0.35e-3 ? FT(-1) :
                ((Dr < 0.9e-3) ? kbr * ΔDr : 2 * (exp(κbr * ΔDr) - 1))

            dNrdt_br = -(ϕ_br + 1) * sc_rai

            # Test
            TT.@test sc_rai ≈ dNrdt_sc rtol = 1e-6
            TT.@test CM2.rain_self_collection(
                SB.pdf_r, SB.self, FT(0), ρ, N_rai,
            ) ≈ FT(0) atol = eps(FT)
            TT.@test br_rai ≈ dNrdt_br rtol = 1e-6
            TT.@test sc_br_rai isa NamedTuple
            TT.@test sc_br_rai.sc ≈ dNrdt_sc rtol = 1e-6
            TT.@test sc_br_rai.br ≈ dNrdt_br rtol = 1e-6

            #setup
            q_rai = FT(0)

            #action
            sc_rai = CM2.rain_self_collection(SB.pdf_r, SB.self, q_rai, ρ, N_rai)
            br_rai = CM2.rain_breakup(SB.pdf_r, SB.brek, q_rai, ρ, N_rai, sc_rai)
            sc_br_rai = CM2.rain_self_collection_and_breakup(SB, q_rai, ρ, N_rai)

            #test
            TT.@test sc_rai ≈ FT(0) atol = eps(FT)
            TT.@test br_rai ≈ FT(0) atol = eps(FT)
            TT.@test sc_br_rai.sc ≈ FT(0) atol = eps(FT)
            TT.@test sc_br_rai.br ≈ FT(0) atol = eps(FT)
        end
    end

    TT.@testset "2M_microphysics - cloud terminal velocity" begin
        #setup
        ρ = FT(1.1)
        q_liq = FT(1e-3)
        N_liq = FT(1e7)

        (; ρw, grav, ν_air) = STVel

        #action
        vt_liq = CM2.cloud_terminal_velocity(SB2006.pdf_c, STVel, q_liq, ρ, N_liq)

        (; νc, μc, ρw) = SB2006.pdf_c
        (; Bc) = CM2.pdf_cloud_parameters_mass(SB2006.pdf_c, q_liq, ρ, N_liq)
        terminal_velocity_prefactor = FT(2 / 9) * (3 / 4 / pi / ρw)^(2 / 3) * (ρw / ρ - 1) * grav / ν_air
        vt0 = terminal_velocity_prefactor * DT.generalized_gamma_Mⁿ(νc, μc, Bc, N_liq, FT(2 / 3)) / N_liq
        vt1 = terminal_velocity_prefactor * DT.generalized_gamma_Mⁿ(νc, μc, Bc, N_liq, FT(5 / 3)) / ρ / q_liq

        #test
        TT.@test vt_liq isa Tuple
        TT.@test vt_liq[1] ≈ vt0 rtol = 1e-6
        TT.@test vt_liq[2] ≈ vt1 rtol = 1e-6

        TT.@test CM2.cloud_terminal_velocity(
            SB2006.pdf_c, STVel, q_liq, ρ, FT(0),
        )[1] ≈ 0 atol = eps(FT)
        TT.@test CM2.cloud_terminal_velocity(
            SB2006.pdf_c, STVel, FT(0), ρ, N_liq,
        )[2] ≈ 0 atol = eps(FT)
        TT.@test CM2.cloud_terminal_velocity(
            SB2006.pdf_c, STVel, FT(0), ρ, FT(0),
        )[1] ≈ 0 atol = eps(FT)
        TT.@test CM2.cloud_terminal_velocity(
            SB2006.pdf_c, STVel, FT(0), ρ, FT(0),
        )[2] ≈ 0 atol = eps(FT)
    end

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 rain terminal velocity with limiters" begin
        #setup
        ρ = FT(1.1)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)

        (; ρ0, aR, bR, cR) = SB2006Vel

        #action
        vt_rai = CM2.rain_terminal_velocity(SB2006, SB2006Vel, q_rai, ρ, N_rai)

        (; Dr_mean) = CM2.pdf_rain_parameters(SB2006.pdf_r, q_rai, ρ, N_rai)
        vt0 = max(0, sqrt(ρ0 / ρ) * (aR - bR / (1 + cR * Dr_mean)))
        vt1 = max(0, sqrt(ρ0 / ρ) * (aR - bR / (1 + cR * Dr_mean)^4))

        #test
        TT.@test vt_rai isa Tuple
        TT.@test vt_rai[1] ≈ vt0 rtol = 1e-6
        TT.@test vt_rai[2] ≈ vt1 rtol = 1e-6

        TT.@test CM2.rain_terminal_velocity(
            SB2006, SB2006Vel, q_rai, ρ, FT(0),
        )[1] ≈ 0 atol = eps(FT)
        TT.@test CM2.rain_terminal_velocity(
            SB2006, SB2006Vel, FT(0), ρ, N_rai,
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
            SB2006_no_limiters, SB2006Vel, q_rai, ρ, N_rai,
        )

        (; Dr_mean) = CM2.pdf_rain_parameters(SB2006_no_limiters.pdf_r, q_rai, ρ, N_rai)
        λr = 1 / Dr_mean
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
            SB2006_no_limiters, SB2006Vel, q_rai, ρ, FT(0),
        )[1] ≈ 0 atol = eps(FT)
        TT.@test CM2.rain_terminal_velocity(
            SB2006_no_limiters, SB2006Vel, FT(0), ρ, N_rai,
        )[2] ≈ 0 atol = eps(FT)
    end

    TT.@testset "2M_microphysics - Chen 2022 rain terminal velocity" begin
        #setup
        ρ = FT(1.1)
        q_rai = FT(5e-4)
        N_rai = FT(1e4)

        for SB in [SB2006, SB2006_no_limiters]
            #action
            vt_rai = CM2.rain_terminal_velocity(SB, Chen2022Vel, q_rai, ρ, N_rai)
            v_bigger = CM2.rain_terminal_velocity(SB, Chen2022Vel, q_rai * 2, ρ, N_rai)

            #test
            TT.@test vt_rai isa Tuple
            TT.@test vt_rai[1] ≈ 1.0738503635546666
            TT.@test vt_rai[2] ≈ 4.00592218028957

            TT.@test CM2.rain_terminal_velocity(
                SB, Chen2022Vel, q_rai, ρ, FT(0),
            )[1] ≈ 0 atol = eps(FT)
            TT.@test CM2.rain_terminal_velocity(
                SB, Chen2022Vel, FT(0), ρ, N_rai,
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
        q_lcl = FT(0)
        q_icl = FT(0)
        q_sno = FT(0)

        for SB in [SB2006, SB2006_no_limiters]

            (; av, bv, α, β, ρ0) = SB.evap
            (; ν_air, D_vapor) = aps

            #action
            evap = CM2.rain_evaporation(SB, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, N_rai, T)

            G = CMC.G_func_liquid(aps, tps, T)
            S = TDI.supersaturation_over_liquid(tps, q_tot, q_lcl + q_rai, q_icl + q_sno, ρ, T)

            (; xr_mean) = CM2.pdf_rain_parameters(SB.pdf_r, q_rai, ρ, N_rai)
            Dr = FT(6 / π / 1000.0)^FT(1 / 3) * xr_mean^FT(1 / 3)
            N_Re = α * xr_mean^β * sqrt(ρ0 / ρ) * Dr / ν_air

            a_vent_0 = av * FT(0.15344374450453543)
            b_vent_0 = bv * FT(0.17380986321413017)
            Fv0 = a_vent_0 + b_vent_0 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)
            a_vent_1 = av * FT(0.5503212081491045)
            b_vent_1 = bv * FT(0.5873135598802672)
            Fv1 = a_vent_1 + b_vent_1 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)

            ∂ₜρn_rai = 2 * FT(π) * G * S * N_rai * Dr * Fv0 / xr_mean
            ∂ₜq_rai = 2 * FT(π) * G * S * N_rai * Dr * Fv1 / ρ

            #test
            TT.@test evap isa @NamedTuple{∂ₜρn_rai::FT, ∂ₜq_rai::FT}
            TT.@test evap.∂ₜρn_rai ≈ ∂ₜρn_rai rtol = 1e-4
            TT.@test evap.∂ₜq_rai ≈ ∂ₜq_rai rtol = 1e-5
            TT.@test CM2.rain_evaporation(
                SB, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, FT(0), T,
            ).∂ₜρn_rai ≈ 0 atol = eps(FT)
            TT.@test CM2.rain_evaporation(
                SB, aps, tps, q_tot, q_lcl, q_icl, FT(0), q_sno, ρ, N_rai, T,
            ).∂ₜq_rai ≈ 0 atol = eps(FT)
        end

        # test limit case: xr = 0 for SB with no limiters
        TT.@test CM2.rain_evaporation(
            SB2006_no_limiters, aps, tps, q_tot, q_lcl, q_icl, FT(0), q_sno, ρ, N_rai, T,
        ).∂ₜρn_rai ≈ 0 atol = eps(FT)

    end

    for pdf_r in [SB2006_no_limiters.pdf_r, SB2006.pdf_r]
        pdf_str = CMP.islimited(pdf_r) ? "with limiters" : "without limiters"
        TT.@testset "2M_microphysics - Seifert and Beheng 2006 rain distribution sanity checks ($pdf_str)" begin

            # air and liquid water densities
            ρₐ = FT(1.2)  # kg/m³
            (; νr, μr, ρw) = pdf_r

            # example number concentration and specific content
            Nᵣ = FT(0.5 * 1e6)   # 0.5 1/cm3
            qᵣ = FT(0.5 * 1e-3)  # 0.5 g/kg

            # distribution parameters for rain
            (; N₀r, Dr_mean) = CM2.pdf_rain_parameters(pdf_r, qᵣ, ρₐ, Nᵣ)
            (; Ar, Br) = CM2.pdf_rain_parameters_mass(pdf_r, qᵣ, ρₐ, Nᵣ)
            TT.@test all(x -> x isa FT, (N₀r, Dr_mean, Ar, Br))

            # mass of liquid droplet as a function of its diameter
            k_m = π * ρw / 6
            m(D) = k_m * D^3

            ### Write the size distribution functions manually
            # rain drop diameter distribution (eq.(3) from 2M docs)
            f_D(D) = N₀r * exp(-D / Dr_mean)
            # rain drop mass distribution (eq.(4) from 2M docs)
            f_x(x) = iszero(x) ? 0 : Ar * x^νr * exp(-Br * x^μr)

            ### Fetch the size distribution functions from the module
            psd = CM2.size_distribution(pdf_r, qᵣ, ρₐ, Nᵣ)

            Mⁿ(n, psd) = y -> y^n * psd(y)

            # integral bounds computed based on the size distribution
            p = FT(1e-6)
            D_min, D_max = CM2.get_size_distribution_bounds(pdf_r, ρₐ, qᵣ, Nᵣ, p)
            x_min = DT.generalized_gamma_quantile(νr, μr, Br, p)
            x_max = DT.generalized_gamma_quantile(νr, μr, Br, 1 - p)

            # Test that these bounds correspond to the correct probability levels
            TT.@test DT.generalized_gamma_cdf(νr, μr, Br, x_min) ≈ p
            TT.@test DT.generalized_gamma_cdf(νr, μr, Br, x_max) ≈ 1 - p
            TT.@test DT.exponential_cdf(Dr_mean, D_min) ≈ p
            TT.@test DT.exponential_cdf(Dr_mean, D_max) ≈ 1 - p

            # Sanity checks for number concentrations for rain
            ND = P3.integrate(f_D, D_min, D_max; quad = P3.ChebyshevGauss(1000))
            Nx = P3.integrate(f_x, x_min, x_max; quad = P3.ChebyshevGauss(100_000))
            ND_psd = P3.integrate(psd, D_min, D_max; quad = P3.ChebyshevGauss(1000))
            TT.@test ND ≈ Nᵣ rtol = 1e-6
            if FT == Float64
                TT.@test Nx ≈ Nᵣ rtol = 7e-3
            else
                TT.@test Nx ≈ Nᵣ rtol = 4e-2  # TODO: poor convergence for Float32
            end
            TT.@test ND_psd == ND

            # Sanity checks for specific contents for rain
            qD = P3.integrate(Mⁿ(3, f_D), D_min, D_max) * k_m / ρₐ
            qx = P3.integrate(Mⁿ(1, f_x), x_min, x_max) / ρₐ
            qD_psd = P3.integrate(Mⁿ(3, psd), D_min, D_max) * k_m / ρₐ
            TT.@test qD ≈ qᵣ rtol = 6e-4
            TT.@test qx ≈ qᵣ rtol = 5e-4
            TT.@test qD_psd == qD

            # Test relationship between exponential moments in diameter space and generalized gamma moments in mass space
            # For raindrops, we expect:
            # - 0th moment in D (number concentration) = 0th moment in mass
            # - 3rd moment in D (mass) = 1st moment in mass
            # - 6th moment in D (mass^2) = 2nd moment in mass
            M⁰_D = DT.exponential_Mⁿ(Dr_mean, Nᵣ, 0)
            M⁰_x = DT.generalized_gamma_Mⁿ(νr, μr, Br, Nᵣ, 0)
            TT.@test M⁰_D ≈ M⁰_x
            TT.@test M⁰_D ≈ Nᵣ
            TT.@test M⁰_x ≈ Nᵣ

            Lᵣ = qᵣ * ρₐ
            M³_D = DT.exponential_Mⁿ(Dr_mean, Nᵣ, 3) * k_m
            M¹_x = DT.generalized_gamma_Mⁿ(νr, μr, Br, Nᵣ, 1)
            TT.@test M³_D ≈ M¹_x
            TT.@test M³_D ≈ Lᵣ
            TT.@test M¹_x ≈ Lᵣ

            # Proportional to radar reflectivity
            M⁶_D = DT.exponential_Mⁿ(Dr_mean, Nᵣ, 6) * k_m^2
            M²_x = DT.generalized_gamma_Mⁿ(νr, μr, Br, Nᵣ, 2)
            TT.@test M⁶_D ≈ M²_x rtol = 1e-6

        end  # end of testset
    end  # end of loop over pdf_r

    TT.@testset "2M_microphysics - Seifert and Beheng 2006 cloud distribution sanity checks" begin

        # example number concentration and specific content
        Nₗ = FT(1e3 * 1e6) # 1000 1/cm3
        qₗ = FT(1e-3)      # 1 g/kg

        # air and liquid water densities in kg/m3
        ρₐ = FT(1.2)
        (; pdf_c) = SB2006
        (; νc, μc, ρw) = pdf_c
        # distribution parameters for cloud
        (; Ac, Bc) = CM2.pdf_cloud_parameters_mass(pdf_c, qₗ, ρₐ, Nₗ)
        (; logN₀c, λc, νcD, μcD) = CM2.pdf_cloud_parameters(pdf_c, qₗ, ρₐ, Nₗ)

        logAc, logBc = CM2.log_pdf_cloud_parameters_mass(pdf_c, qₗ, ρₐ, Nₗ)
        TT.@test all(x -> x isa FT, (Ac, Bc, logN₀c, λc, νcD, μcD, logAc, logBc))

        # mass of liquid droplet as a function of its diameter
        k_m = π * ρw / 6
        m(D) = k_m * D^3

        # cloud droplet mass distribution (Eq. (2) from 2M docs, but in log space)
        logf_x(x) = logAc + νc * log(x) - Bc * x^μc
        f_x(x) = exp(logf_x(x))

        # cloud droplet diameter distribution (Eq. (6) from 2M docs)
        logf_D(D) = logN₀c + (3νc + 2) * log(D) - λc * D^(3μc)
        f_D(D) = exp(logf_D(D))

        psd = CM2.size_distribution(pdf_c, qₗ, ρₐ, Nₗ)

        Mⁿ(n, psd) = y -> y^n * psd(y)

        # integral bounds guesstimated in meters for the mass distribution
        p = FT(1e-6)
        D_min, D_max = CM2.get_size_distribution_bounds(pdf_c, ρₐ, qₗ, Nₗ, p)
        x_min = DT.generalized_gamma_quantile(νc, μc, Bc, p)
        x_max = DT.generalized_gamma_quantile(νc, μc, Bc, 1 - p)

        # Test that these bounds correspond to the correct probability levels
        TT.@test DT.generalized_gamma_cdf(νc, μc, Bc, x_min) ≈ p
        TT.@test DT.generalized_gamma_cdf(νc, μc, Bc, x_max) ≈ 1 - p
        TT.@test DT.generalized_gamma_cdf(νcD, μcD, λc, D_min) ≈ p
        TT.@test DT.generalized_gamma_cdf(νcD, μcD, λc, D_max) ≈ 1 - p


        # Sanity checks of specific content and number concentration with mass distribution
        Nx = P3.integrate(Mⁿ(0, f_x), x_min, x_max)
        qx = P3.integrate(Mⁿ(1, f_x), x_min, x_max) / ρₐ
        TT.@test qx ≈ qₗ rtol = 2e-5
        TT.@test Nx ≈ Nₗ rtol = 1e-5

        # Sanity checks of specific content and number concentration with diameter distribution
        ND = P3.integrate(Mⁿ(0, f_D), D_min, D_max)
        ND_psd = P3.integrate(Mⁿ(0, psd), D_min, D_max)
        qD = P3.integrate(Mⁿ(3, f_D), D_min, D_max) * k_m / ρₐ
        qD_psd = P3.integrate(Mⁿ(3, psd), D_min, D_max) * k_m / ρₐ
        TT.@test ND ≈ Nₗ rtol = 1e-5
        TT.@test ND_psd ≈ Nₗ rtol = 1e-5
        TT.@test qD ≈ qₗ rtol = 2e-5
        TT.@test qD_psd ≈ qₗ rtol = 2e-5
    end

    TT.@testset "2M_microphysics - Horn 2012 number concentration adjustment" begin

        # Setup
        ρ = FT(1.2)       # kg/m³
        q = FT(1e-3)      # kg/kg
        x_min = FT(2.6e-10)   # kg
        x_max = FT(5e-6)      # kg
        NumAdj = SB2006.numadj
        (; τ) = NumAdj

        N_low = FT(1e2)        # 1/m³
        N_inrange = FT(1e4)    # 1/m³
        N_high = FT(1e7)       # 1/m³
        N_veryhigh = FT(1e15)  # 1/m³

        # Action
        dN_dt_low = (ρ * q / x_max - N_low) / τ
        dN_dt_high = (ρ * q / x_min - N_high) / τ

        # Test
        TT.@test CM2.number_increase_for_mass_limit(NumAdj, x_max, q, ρ, N_low) ≈ dN_dt_low
        TT.@test CM2.number_increase_for_mass_limit(NumAdj, x_max, q, ρ, N_inrange) ≈ FT(0)
        TT.@test CM2.number_increase_for_mass_limit(NumAdj, x_max, q, ρ, N_high) ≈ FT(0)
        TT.@test CM2.number_increase_for_mass_limit(NumAdj, x_max, FT(0), ρ, N_inrange) ≈ FT(0)

        TT.@test CM2.number_decrease_for_mass_limit(NumAdj, x_min, q, ρ, N_low) ≈ FT(0)
        TT.@test CM2.number_decrease_for_mass_limit(NumAdj, x_min, q, ρ, N_inrange) ≈ FT(0)
        TT.@test CM2.number_decrease_for_mass_limit(NumAdj, x_min, q, ρ, N_high) ≈ dN_dt_high
        TT.@test CM2.number_decrease_for_mass_limit(NumAdj, x_min, FT(0), ρ, N_inrange) ≈ -N_inrange / τ
        TT.@test CM2.number_decrease_for_mass_limit(NumAdj, FT(0), q, ρ, N_veryhigh) ≈ FT(0)
        TT.@test CM2.number_decrease_for_mass_limit(NumAdj, FT(0), FT(0), ρ, N_veryhigh) ≈ FT(0)
    end
end

TT.@testset "Microphysics 2M Tests ($FT)" for FT in (Float64, Float32)
    test_microphysics2M(FT)
end
nothing
