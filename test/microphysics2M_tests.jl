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
    SB2006_no_limiters = CMP.SB2006(toml_dict; is_limited = false)
    # `is_limited = true` now builds the single mean-mass window, so a test that characterises
    # the RETIRED SB2006 clamp cascade has to name it. Kept because the cascade is still what the
    # comparison arms and the sign-floor assertions are measured against.
    cascade_pdf = CMP.RainParticlePDF_SB2006_limited(toml_dict)
    SB2006_cascade = CMP.SB2006(toml_dict; rain_pdf = cascade_pdf)

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

        # output should be zero if either q_liq or q_rai are zero
        q_lcl = FT(0)
        q_rai = FT(1e-6)

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
    TT.@testset "single mean-mass window - rain PDF parameters" begin
        # `RainParticlePDF_SB2006_windowed` bounds one quantity, the mean drop mass, and lets
        # λ and N₀ inherit their ranges through the exponential-PSD identities. The two halves
        # asserted here are (a) the returned triple always describes ONE exponential PSD
        # consistent with the rain mass the state carries, and (b) the bounded quantity is in
        # fact bounded. The cascade satisfies (b) but not (a).
        win = CMP.RainParticlePDF_SB2006_windowed(toml_dict)
        (; xr_min, xr_max, ρw) = win
        ρₐ = FT(1.2)
        tol = FT === Float32 ? FT(1e-4) : FT(1e-10)

        # The grid straddles the presence gate at Float32, where `eps(FT) = 1.2e-7`; the gated
        # states return the zero triple and are asserted separately below, not here.
        for q in FT.((1e-9, 1e-6, 1.5e-4, 1e-3, 5e-3)),
            N in FT.((1e-2, 30, 1e3, 1e5, 1e7))

            (; N₀r, Dr_mean, xr_mean) = CM2.pdf_rain_parameters(win, q, ρₐ, N)
            if all(iszero, (N₀r, Dr_mean, xr_mean))
                TT.@test q < eps(FT) || N < eps(FT)
                continue
            end
            L = ρₐ * q
            # (a) one PSD: mean mass and mean diameter agree; the PSD's mass integral returns
            # the rain mass of the state rather than a clamped surrogate; and its number
            # integral returns the state's number bounded by the one window, not by a second.
            TT.@test π * ρw * Dr_mean^3 ≈ xr_mean rtol = tol
            TT.@test π * ρw * N₀r * Dr_mean^4 ≈ L rtol = tol
            TT.@test N₀r * Dr_mean ≈ clamp(N, L / xr_max, L / xr_min) rtol = tol
            # (b) the bounded quantity is bounded
            TT.@test xr_min * (1 - tol) <= xr_mean <= xr_max * (1 + tol)
        end

        # The seeded case from the design note: a state entirely inside the sanctioned mass
        # range, which the intercept window alone corrupts. The windowed inversion returns the
        # honest mean size and the state's own number; the cascade returns a PSD describing
        # several times as many drops as the state has.
        q, N = FT(1.5e-4), FT(30)
        w = CM2.pdf_rain_parameters(win, q, FT(1), N)
        l = CM2.pdf_rain_parameters(cascade_pdf, q, FT(1), N)
        TT.@test w.xr_mean ≈ q / N rtol = tol
        TT.@test w.N₀r * w.Dr_mean ≈ N rtol = tol
        TT.@test !isapprox(l.N₀r * l.Dr_mean, N; rtol = FT(0.1))

        # Empty on either moment, matching the presence gate the consuming rates apply.
        for (q, N) in ((FT(0), FT(0)), (FT(0), FT(1e5)), (FT(1e-4), FT(0)))
            TT.@test all(iszero, CM2.pdf_rain_parameters(win, q, ρₐ, N))
            TT.@test all(iszero, CM2.get_size_distribution_bounds(win, q, ρₐ, N))
            n = CM2.size_distribution(win, q, ρₐ, N)
            TT.@test all(iszero, (n(0), n(0.1), n(Inf)))
        end

        # Degenerate states stay bounded: mass without number is repaired by the window into
        # the largest sanctioned drop rather than described as a metre-scale one.
        for (q, N) in ((FT(1e-4), eps(FT)), (FT(5e-3), FT(1e-9)), (eps(FT), FT(1e5)))
            p = CM2.pdf_rain_parameters(win, q, ρₐ, N)
            TT.@test all(isfinite, p)
            TT.@test p.xr_mean <= xr_max * (1 + tol)
        end
    end

    TT.@testset "limiting lambda_r and x_r - Seifert and Beheng 2006 (retired cascade)" begin
        # Characterises the RETIRED cascade, so it names it: `λ_min`/`λ_max` are fields of that
        # variant alone. Under the mean-mass window λ has no window of its own - it inherits its
        # range from the mass bound through `λ = cbrt(π ρw / x̄)`, which the windowed testset
        # above asserts directly.
        q_rai = [FT(0), FT(1e-3), FT(1e-4), FT(1e-2)]
        N_rai = [FT(1e1), FT(1e1), FT(1e3), FT(1e5)]
        ρ = FT(1)

        (; xr_min, xr_max, λ_min, λ_max) = cascade_pdf

        for Nr in N_rai
            for qr in q_rai
                #action
                (; Dr_mean, xr_mean) = CM2.pdf_rain_parameters(cascade_pdf, qr, ρ, Nr)
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
                N_lcl,
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
            bound_factor_au = CM2.mean_mass_bound_factor(Lc / N_lcl, x_star)
            bound_factor_sc = CM2.mean_mass_bound_factor(Lc / N_lcl, x_star)
            τ = 1 - Lc / (Lc + Lr)
            ϕ_au = 400 * τ^0.7 * (1 - τ^0.7)^3
            dqrdt_au =
                kcc / 20 / x_star * (νc + 2) * (νc + 4) / (νc + 1)^2 *
                Lc^2 *
                xc^2 *
                (1 + ϕ_au / (1 - τ)^2) *
                (ρ0 / ρ) / ρ * bound_factor_au
            dqcdt_au = -dqrdt_au
            dNcdt_au = 2 / x_star * ρ * dqcdt_au
            dNrdt_au = -0.5 * dNcdt_au
            dNcdt_sc =
                -kcc * (νc + 2) / (νc + 1) * (ρ0 / ρ) * Lc^2 * bound_factor_sc - au.dN_lcl_dt

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
                N_lcl,
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

    TT.@testset "rain self-collection and breakup as a relaxation" begin
        # The pair relaxes the rain number toward the collisional equilibrium n_eq = L/x_eq.
        # Written that way its frozen-shape diagonal is -1/τ_eff ≤ 0 everywhere, where the
        # donor-linearization recipe (sc + br)/N goes POSITIVE wherever breakup dominates.
        # The recast needs a PSD whose parameters are honest functions of (L, N), so it is
        # exercised on the windowed variant; the cascade arm is here to show that it is not
        # merely a stylistic choice.
        win = CMP.RainParticlePDF_SB2006_windowed(toml_dict)
        (; self, brek) = SB2006
        ρ = FT(1.1)
        x_eq = FT(π) / 6 * win.ρw * brek.Deq^3
        tol = FT === Float32 ? FT(1e-4) : FT(1e-10)

        # equilibrium: n_eq is L/x_eq, and the pair vanishes there
        for q in FT.((1e-6, 1e-4, 2e-3))
            N_eq = CM2.rain_equilibrium_number(brek, win, q, ρ)
            TT.@test N_eq ≈ ρ * q / x_eq rtol = tol
            rel = CM2.rain_number_relaxation(win, self, brek, q, ρ, N_eq)
            TT.@test rel.N_eq ≈ N_eq rtol = tol
            # the rate is zero at the fixed point, measured against its size a factor 2 away
            scale = abs(CM2.rain_number_relaxation(win, self, brek, q, ρ, 2 * N_eq).∂ₜN_rai)
            TT.@test abs(rel.∂ₜN_rai) <= tol * scale
            # and the removable limit is finite and damping there
            TT.@test isfinite(rel.inv_τ_eff)
            TT.@test rel.inv_τ_eff > 0
        end

        # the tendency is the composition, to the bit, and the recast round trip reproduces it
        for q in FT.((1e-6, 1e-4, 2e-3)), r in FT.((0.05, 0.5, 0.9, 1.1, 3.0, 30.0))
            N = r * CM2.rain_equilibrium_number(brek, win, q, ρ)
            rel = CM2.rain_number_relaxation(win, self, brek, q, ρ, N)
            sc = CM2.rain_self_collection(win, self, q, ρ, N)
            br = CM2.rain_breakup(win, brek, q, ρ, N, sc)
            TT.@test rel.∂ₜN_rai === sc + br
            TT.@test -(N - rel.N_eq) * rel.inv_τ_eff ≈ sc + br rtol = tol
            # the sign structure, which is what the reformulation exists for
            TT.@test -rel.inv_τ_eff <= 0
            TT.@test sign(rel.∂ₜN_rai) == sign(rel.N_eq - N)
        end

        # the donor recipe this replaces IS anti-damping above the equilibrium diameter, so
        # the two are not interchangeable and the testset can fail in both directions
        q = FT(1e-4)
        N_small = FT(0.2) * CM2.rain_equilibrium_number(brek, win, q, ρ)   # Dr > Deq
        rel = CM2.rain_number_relaxation(win, self, brek, q, ρ, N_small)
        TT.@test rel.∂ₜN_rai / N_small > 0        # the donor recipe's entry
        TT.@test -rel.inv_τ_eff < 0               # the relaxation's entry

        # THE INEQUALITY THE SIGN GUARANTEE ACTUALLY RESTS ON, named here because nothing
        # else in the code ties these parameters together. Form A is clean on the windowed
        # variant for exactly one reason: clamping the mean mass into [x_min, x_max] cannot move
        # it ACROSS x_eq, because x_eq lies strictly inside the window. Measured (job 6940654,
        # 850 states, both precisions): push x_max below x_eq or x_min above it and the sign
        # correspondence breaks on 425 of 850 states, while a window only 10 percent wide that
        # still CONTAINS x_eq is clean. So it is containment that matters, not width.
        #
        # `Deq` comes from the SB2006 breakup fit and the window ends come from the cloud/rain
        # separation mass and the tail-mass criterion - three parameters set independently in
        # ClimaParams. If a future change to any of them breaks this, the sign floor keeps the
        # matrix safe and this assertion says why the floor started binding.
        TT.@test win.xr_min < x_eq < win.xr_max
        TT.@test CMP.RainParticlePDF_SB2006_notlimited(toml_dict).xr_min < x_eq
        TT.@test x_eq < CMP.RainParticlePDF_SB2006_notlimited(toml_dict).xr_max

        # THE PRECONDITION, asserted on the variant that VIOLATES it. Form A's sign guarantee
        # rests on the mean size tracking the number, which the clamp cascade breaks: there
        # `Φ_br`'s sign stops following `n_rai - n_eq` and the quotient can go negative. The
        # floor is what keeps `-1/τ_eff ≤ 0` structural on ANY variant a host builds, so this
        # loop runs the CASCADE deliberately - the windowed loops above cannot exercise it.
        cascade = cascade_pdf
        floor_bound = false
        for q in FT.((1e-5, 1e-4, 1e-3, 5e-3)),
            r in FT.((0.05, 0.5, 0.9, 1.19, 2.0, 10.0))

            N = r * CM2.rain_equilibrium_number(brek, cascade, q, ρ)
            rel = CM2.rain_number_relaxation(cascade, self, brek, q, ρ, N)
            TT.@test rel.inv_τ_eff >= 0
            TT.@test -rel.inv_τ_eff <= 0
            # the tendency is the composition either way: the floor touches the timescale only
            sc = CM2.rain_self_collection(cascade, self, q, ρ, N)
            TT.@test rel.∂ₜN_rai === sc + CM2.rain_breakup(cascade, brek, q, ρ, N, sc)
            # did the floor actually bind anywhere? if not, this testset proves nothing
            sign(rel.∂ₜN_rai) != sign(rel.N_eq - N) && !iszero(rel.∂ₜN_rai) &&
                (floor_bound = true)
        end
        TT.@test floor_bound

        # and it is INERT on the honest variant, which is the other half of the claim
        for q in FT.((1e-5, 1e-4, 1e-3, 5e-3)),
            r in FT.((0.05, 0.5, 0.9, 1.19, 2.0, 10.0))

            N = r * CM2.rain_equilibrium_number(brek, win, q, ρ)
            rel = CM2.rain_number_relaxation(win, self, brek, q, ρ, N)
            TT.@test sign(rel.∂ₜN_rai) == sign(rel.N_eq - N)
        end

        # empty and degenerate states: no relaxation, nothing manufactured
        for (q, N) in ((FT(0), FT(0)), (FT(0), FT(1e4)), (FT(1e-4), FT(0)))
            rel = CM2.rain_number_relaxation(win, self, brek, q, ρ, N)
            TT.@test rel.∂ₜN_rai == FT(0)
            TT.@test rel.inv_τ_eff == FT(0)
            TT.@test isfinite(rel.N_eq)
        end
        TT.@test CM2.rain_equilibrium_number(brek, win, FT(0), ρ) == FT(0)
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
        # The number consistent with the windowed mean mass, not the state's own. This test's
        # parameter set is `SB2006_limiters.toml`, whose `xc_max` is 6.54e-11 kg rather than the
        # 2.6e-10 default, so the ceiling BINDS at this otherwise ordinary cloud state (1 g/m³ on
        # 1e7 droplets/m³, a 59 μm mean droplet). Dividing a moment of the windowed distribution
        # by the state's mass while counting the state's droplets mixes two populations - the
        # inconsistency #74 removes - so the reference has to use the same pair production does,
        # or it asserts the defect.
        #
        # THIS STATE IS THE CHEAPEST IN-REPO CASE OF THE CEILING ACTUALLY ENGAGING, and it is
        # what decides the FORM of the fix rather than merely confirming its arithmetic. Clamping
        # the mean mass while leaving the number alone would divide the moment of a windowed
        # distribution by an unwindowed mass, so this perfectly ordinary cloud would sediment its
        # mass a factor `x̄_raw / xc_max` ≈ 1.68 too SLOWLY - a real warm-rain perturbation at a
        # shipped configuration, not a change confined to degenerate states. Rebuilding the
        # number holds it at the ceiling instead, which is the fall speed of the largest droplet
        # the category admits. That is why the pair is rebuilt together.
        #
        # Asserted rather than described, so the case cannot quietly go vacuous if the limiters
        # file or this state ever changes: the ceiling must genuinely bind here.
        TT.@test ρ * q_liq / N_liq > SB2006.pdf_c.xc_max
        x̄_ref, N_eff = CM2.cloud_mean_droplet_mass_and_number(SB2006.pdf_c, q_liq, ρ, N_liq)
        TT.@test N_eff > N_liq
        terminal_velocity_prefactor = FT(2 / 9) * (3 / 4 / pi / ρw)^(2 / 3) * (ρw / ρ - 1) * grav / ν_air
        vt0 = terminal_velocity_prefactor * DT.generalized_gamma_Mⁿ(νc, μc, Bc, N_eff, FT(2 / 3)) / N_eff
        # The distribution's own mass, `x̄ * N_eff`, is what the mass-weighted moment is
        # normalised by. At THIS state that equals the content, because the ceiling binds here
        # and the rebuild is mass-preserving; the identity is asserted rather than relied on, so
        # that a reference written as `/ ρ / q_liq` cannot look correct in general on the
        # strength of a state where the two happen to agree.
        TT.@test x̄_ref * N_eff ≈ ρ * q_liq rtol = 10 * eps(FT)
        vt1 = terminal_velocity_prefactor * DT.generalized_gamma_Mⁿ(νc, μc, Bc, N_eff, FT(5 / 3)) / (x̄_ref * N_eff)

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

    TT.@testset "2M_microphysics - the cloud mean-mass window is two-sided and consistent" begin
        # #74. The cloud PSD builder floored the mean droplet mass and left the ceiling open, so
        # an unbounded quotient reached the Stokes velocity law - 16,096 m/s at a per-stage
        # environment state in the r2 EDMF box record, from inputs clearing every absolute floor
        # in the call. The window closes it, and the number is rebuilt with it so that BOTH
        # velocity components are functions of the windowed mean mass alone.
        (; xc_min, xc_max) = SB2006.pdf_c
        ρ = FT(1.15308)
        q = FT(4.2275e-7)

        vel(x) = CM2.cloud_terminal_velocity(
            SB2006.pdf_c, STVel, q, ρ, FT(ρ * q / x),
        )

        # THE CENTRAL RUNG: above the ceiling the velocity is CONSTANT, because the windowed mean
        # mass is all either component depends on. This is the rung that discriminates - within a
        # single ray the velocity is invariant on EVERY variant (both moments are linear in N,
        # so the ratio is N-free whatever the builder did), and a testset asserting only that
        # would pass on unfixed code. Sampling ACROSS rays spanning decades above the ceiling is
        # what separates the forms: unfixed, the velocity climbs as the quotient does; with the
        # mean mass clamped but the number left alone, it FALLS as 1/x̄_raw. Only the consistent
        # pair holds it flat - the same flat-across-decades signature the rain sibling shows.
        above = [xc_max * FT(f) for f in (FT(2), FT(10), FT(1e2), FT(1e4), FT(1e6), FT(1.8e7))]
        v_ref = vel(above[1])
        for x in above
            v = vel(x)
            TT.@test v[1] ≈ v_ref[1] rtol = 100 * eps(FT)
            TT.@test v[2] ≈ v_ref[2] rtol = 100 * eps(FT)
        end

        # ... and BELOW the ceiling it must genuinely VARY. Without this the constancy rung could
        # be satisfied by a velocity that is constant everywhere, which would be a far worse bug
        # than the one being fixed.
        below = [xc_max * FT(f) for f in (FT(1e-3), FT(1e-2), FT(0.1), FT(0.5))]
        v_below = [vel(x) for x in below]
        for k in 2:length(below)
            TT.@test v_below[k][2] > v_below[k - 1][2]
        end
        TT.@test v_below[end][2] < v_ref[2]

        # The measured killing row itself: the state that produced the observed velocity, at the
        # ceiling on both components after the fix.
        v_kill = CM2.cloud_terminal_velocity(SB2006.pdf_c, STVel, q, ρ, FT(1.03679e-4))
        TT.@test v_kill[1] ≈ v_ref[1] rtol = 100 * eps(FT)
        TT.@test v_kill[2] ≈ v_ref[2] rtol = 100 * eps(FT)
        # and it is a physical fall speed for the largest droplet the category admits, not merely
        # a finite one - the check that still fails if the clamp is applied at the wrong end.
        TT.@test FT(0) < v_kill[2] < FT(1)

        # The window's own arithmetic, asserted rather than assumed: the state is outside the
        # window and both of the call's absolute presence floors are OPEN there, which is why
        # nothing else in the call can bound it.
        TT.@test ρ * q / FT(1.03679e-4) > xc_max
        TT.@test q > CM.Utilities.ϵ_numerics_2M_M(FT)
        TT.@test FT(1.03679e-4) > CM.Utilities.ϵ_numerics_2M_N(FT)

        # INERTNESS inside the window: the rebuilt number is the state's own number exactly
        # wherever the ceiling does not bind, so nothing inside the window moves at all.
        for x in below
            N = FT(ρ * q / x)
            TT.@test CM2.cloud_mean_droplet_mass_and_number(SB2006.pdf_c, q, ρ, N)[2] === N
        end
        # and the floor keeps the number, which is the asymmetry the helper's docstring argues:
        # a mass-preserving rebuild there would return zero droplets for a real population.
        N_hi = FT(ρ * q / (FT(1e-3) * xc_min))
        TT.@test CM2.cloud_mean_droplet_mass_and_number(SB2006.pdf_c, q, ρ, N_hi)[2] === N_hi
        TT.@test CM2.cloud_mean_droplet_mass_and_number(SB2006.pdf_c, q, ρ, N_hi)[1] ≈ xc_min

        # THE MIRROR RUNG, and it is the one the mass-weighted normalisation turns on. Everything
        # above characterises the ceiling; the floor was asserted on the helper alone, so the
        # VELOCITY was uncovered in the only regime the normalisation changes.
        #
        # Below the floor the mean mass is held at `xc_min` while the number is deliberately left
        # alone, so the distribution carries MORE mass than the state does. Normalising the
        # mass-weighted moment by the state's content therefore divides a numerator frozen at the
        # floor by a shrinking denominator, and the fall speed climbs in exact proportion to the
        # amplitude: measured on these rays the retired form reaches 1.5e4 m/s eight decades below
        # the floor, against 1.5e-4 m/s on it, and at a state taken from an EDMF box record it
        # returned 1.19e4 m/s where rain was 6.10 and ice 2.08. Normalised by the distribution's
        # own mass it is FLAT, for the same reason the ceiling rung is flat: both components are
        # functions of the windowed mean mass alone.
        #
        # Sampling ACROSS rays is again what discriminates. Within one ray every variant is
        # N-free, since both moments are linear in the amplitude.
        # The reference sits ON the floor, where the clamp is inclusive, and the rays sit strictly
        # below it. The two are kept apart because `ρ q / (ρ q / x)` is a round trip that need not
        # land back on `x`: a ray placed exactly at the boundary can come back a hair above it and
        # fail the strictness check below for a reason that is arithmetic rather than physics.
        below_floor = [xc_min * FT(f) for f in (FT(1e-1), FT(1e-2), FT(1e-4), FT(1e-6), FT(1e-8))]
        v_floor = vel(xc_min)
        for x in below_floor
            v = vel(x)
            TT.@test v[1] ≈ v_floor[1] rtol = 100 * eps(FT)
            TT.@test v[2] ≈ v_floor[2] rtol = 100 * eps(FT)
        end
        # and it is the fall speed of the SMALLEST droplet the category admits: physical, and
        # strictly slower than the ceiling value, which is the check that still fails if the two
        # bounds are ever applied at the wrong ends.
        TT.@test FT(0) < v_floor[2] < v_ref[2]
        # The rays are genuinely below the floor and genuinely inside both absolute presence
        # gates, so nothing else in the call can be doing the bounding.
        for x in below_floor
            TT.@test ρ * q / FT(ρ * q / x) < xc_min
        end
        TT.@test q > CM.Utilities.ϵ_numerics_2M_M(FT)
        TT.@test FT(ρ * q / below_floor[end]) > CM.Utilities.ϵ_numerics_2M_N(FT)
    end

    TT.@testset "2M_microphysics - SB2006 rain terminal velocity, retired cascade (untruncated)" begin
        # The reference expression here has NO moment factors, i.e. it is the untruncated
        # integral. That is a property of the retired cascade, not of limiting, so this testset
        # names the cascade. The windowed and unbounded variants both truncate the integral at the
        # diameter where the individual-drop fit turns positive; the testset below covers that
        # form and the one after it asserts the two agree.
        ρ = FT(1.1)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)

        (; ρ0, aR, bR, cR) = SB2006Vel

        #action
        vt_rai = CM2.rain_terminal_velocity(SB2006_cascade, SB2006Vel, q_rai, ρ, N_rai)

        (; Dr_mean) = CM2.pdf_rain_parameters(cascade_pdf, q_rai, ρ, N_rai)
        vt0 = max(0, sqrt(ρ0 / ρ) * (aR - bR / (1 + cR * Dr_mean)))
        vt1 = max(0, sqrt(ρ0 / ρ) * (aR - bR / (1 + cR * Dr_mean)^4))

        #test
        TT.@test vt_rai isa Tuple
        TT.@test vt_rai[1] ≈ vt0 rtol = 1e-6
        TT.@test vt_rai[2] ≈ vt1 rtol = 1e-6

        TT.@test CM2.rain_terminal_velocity(
            SB2006_cascade, SB2006Vel, q_rai, ρ, FT(0),
        )[1] ≈ 0 atol = eps(FT)
        TT.@test CM2.rain_terminal_velocity(
            SB2006_cascade, SB2006Vel, FT(0), ρ, N_rai,
        )[2] ≈ 0 atol = eps(FT)
    end

    TT.@testset "2M_microphysics - the window inherits the truncated velocity integral" begin
        # Rider (b) of D07, asserted rather than only documented: the windowed variant dispatches
        # to the same truncated helper as the unbounded one, so at a state INSIDE the window -
        # where the two inversions return the same PSD - the two must give the same fall speeds.
        # The cascade at the same state does not, and that difference is the truncation, not the
        # limiting.
        ρ = FT(1.1)
        q_rai = FT(1e-6)
        N_rai = FT(1e4)
        win_pdf = CMP.RainParticlePDF_SB2006_windowed(toml_dict)
        x̄ = ρ * q_rai / N_rai
        TT.@test win_pdf.xr_min < x̄ < win_pdf.xr_max      # the window does not bind here

        v_win = CM2.rain_terminal_velocity(SB2006, SB2006Vel, q_rai, ρ, N_rai)
        v_raw = CM2.rain_terminal_velocity(SB2006_no_limiters, SB2006Vel, q_rai, ρ, N_rai)
        v_cas = CM2.rain_terminal_velocity(SB2006_cascade, SB2006Vel, q_rai, ρ, N_rai)
        TT.@test v_win[1] ≈ v_raw[1] rtol = FT(1e-5)
        TT.@test v_win[2] ≈ v_raw[2] rtol = FT(1e-5)
        TT.@test !isapprox(v_win[2], v_cas[2]; rtol = FT(0.1))
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
            ND = P3.integrate(f_D, D_min, D_max, P3.GaussLegendre(FT, 16))
            Nx = P3.integrate(f_x, x_min, x_max, P3.GaussLegendre(FT, 45_000))
            ND_psd = P3.integrate(psd, D_min, D_max, P3.GaussLegendre(FT, 16))
            # D_min/D_max truncate 2p of the distribution, so the exact
            # truncated integral is Nᵣ(1 - 2p); ND recovers Nᵣ only to O(2p) = 2e-6.
            TT.@test ND ≈ Nᵣ rtol = 3e-6
            if FT == Float64
                TT.@test Nx ≈ Nᵣ rtol = 7e-3
            else
                TT.@test Nx ≈ Nᵣ rtol = 4e-2  # TODO: poor convergence for Float32
            end
            TT.@test ND_psd == ND

            # Sanity checks for specific contents for rain
            qD = P3.integrate(Mⁿ(3, f_D), D_min, D_max, P3.GaussLegendre(FT, 96)) * k_m / ρₐ
            qx = P3.integrate(Mⁿ(1, f_x), x_min, x_max, P3.GaussLegendre(FT, 96)) / ρₐ
            qD_psd = P3.integrate(Mⁿ(3, psd), D_min, D_max, P3.GaussLegendre(FT, 96)) * k_m / ρₐ
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
        Nx = P3.integrate(Mⁿ(0, f_x), x_min, x_max, P3.GaussLegendre(FT, 32))
        qx = P3.integrate(Mⁿ(1, f_x), x_min, x_max, P3.GaussLegendre(FT, 32)) / ρₐ
        TT.@test qx ≈ qₗ rtol = 2e-5
        TT.@test Nx ≈ Nₗ rtol = 1e-5

        # Sanity checks of specific content and number concentration with diameter distribution
        ND = P3.integrate(Mⁿ(0, f_D), D_min, D_max, P3.GaussLegendre(FT, 32))
        ND_psd = P3.integrate(Mⁿ(0, psd), D_min, D_max, P3.GaussLegendre(FT, 32))
        qD = P3.integrate(Mⁿ(3, f_D), D_min, D_max, P3.GaussLegendre(FT, 32)) * k_m / ρₐ
        qD_psd = P3.integrate(Mⁿ(3, psd), D_min, D_max, P3.GaussLegendre(FT, 32)) * k_m / ρₐ
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

        # number_tendency_from_mass_limits relaxes the mean mass x = q / n into
        # [x_min, x_max], covering both the upper (increase) and lower (decrease)
        # bounds. Work in specific quantities: q [kg/kg], n [1/kg].
        numadj_nt = (; x_min, x_max, τ)
        n_low, n_inrange, n_high = N_low / ρ, N_inrange / ρ, N_high / ρ
        # x > x_max: relax up towards q / x_max (positive tendency)
        TT.@test CM2.number_tendency_from_mass_limits(numadj_nt, q, n_low) ≈ (q / x_max - n_low) / τ
        # x within [x_min, x_max]: no adjustment
        TT.@test CM2.number_tendency_from_mass_limits(numadj_nt, q, n_inrange) ≈ FT(0)
        # x < x_min: relax down towards q / x_min (negative tendency)
        TT.@test CM2.number_tendency_from_mass_limits(numadj_nt, q, n_high) ≈ (q / x_min - n_high) / τ
        # q ≈ 0: target is zero number
        TT.@test CM2.number_tendency_from_mass_limits(numadj_nt, FT(0), n_inrange) ≈ -n_inrange / τ
        # x_min = 0: q / x_min is Inf, so the lower-bound relaxation never fires
        TT.@test CM2.number_tendency_from_mass_limits((; x_min = FT(0), x_max, τ), q, n_high) ≈ FT(0)
    end
end

TT.@testset "Microphysics 2M Tests ($FT)" for FT in (Float64, Float32)
    test_microphysics2M(FT)
end
nothing
