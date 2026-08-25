import Test as TT

import ClimaParams as CP

import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.Parameters as CMP

"""
    write_prescribed_aerosol_override(dir, entries)

Write a ClimaParams override file into `dir` carrying the prescribed-aerosol TOML names in
`entries` (a `NamedTuple` keyed by `PrescribedAerosol` field name), and return its path.
"""
function write_prescribed_aerosol_override(dir, entries)
    path = joinpath(dir, "prescribed_aerosol_override.toml")
    open(path, "w") do io
        for (field, value) in pairs(entries)
            name = getfield(CMP.PRESCRIBED_AEROSOL_TOML_NAMES, field)
            println(io, "[$name]")
            println(io, "value = $value")
            println(io, "type = \"float\"")
            println(io, "description = \"test override\"")
            println(io)
        end
    end
    return path
end

function test_prescribed_aerosol(FT)
    TT.@testset "PrescribedAerosol code defaults" begin
        pa = CMP.PrescribedAerosol(FT)
        TT.@test pa isa CMP.PrescribedAerosol{FT}
        # The defaults describe a physically admissible two-mode population, which is what the
        # activation scheme needs: positive radii and hygroscopicities, widths above one.
        TT.@test pa.r_dry_accum > 0
        TT.@test pa.r_dry_coarse > 0
        TT.@test pa.κ_accum > 0
        TT.@test pa.κ_coarse > 0
        TT.@test pa.stdev_accum > 1
        TT.@test pa.stdev_coarse > 1
        TT.@test pa.N_accum >= 0
        TT.@test pa.N_coarse >= 0
        # A stock parameter dictionary carries none of the TOML names, so it reproduces the
        # code defaults exactly.
        TT.@test CMP.PrescribedAerosol(CP.create_toml_dict(FT)) == pa
        # Keyword arguments name fields directly.
        TT.@test CMP.PrescribedAerosol(FT; N_accum = FT(1e8)).N_accum == FT(1e8)
    end

    TT.@testset "PrescribedAerosol TOML override path" begin
        mktempdir() do dir
            override_file = write_prescribed_aerosol_override(
                dir, (; N_accum = 1.2e8, κ_coarse = 0.9),
            )
            toml_dict = CP.create_toml_dict(FT; override_file)
            pa = CMP.PrescribedAerosol(toml_dict)
            default = CMP.PrescribedAerosol(FT)
            TT.@test pa.N_accum == FT(1.2e8)
            TT.@test pa.κ_coarse == FT(0.9)
            # Fields the override does not mention keep the code default.
            TT.@test pa.N_coarse == default.N_coarse
            TT.@test pa.r_dry_accum == default.r_dry_accum
            # An explicit keyword argument takes precedence over the override file.
            TT.@test CMP.PrescribedAerosol(toml_dict; N_accum = FT(7e6)).N_accum == FT(7e6)
        end
    end

    TT.@testset "Aerosol distribution built from PrescribedAerosol" begin
        pa = CMP.PrescribedAerosol(FT)
        ad = AM.aerosol_distribution(pa)
        TT.@test AM.n_modes(ad) == 2
        TT.@test isbits(ad)
        TT.@test ad.modes[1].N == pa.N_accum
        TT.@test ad.modes[2].N == pa.N_coarse
    end
end

function test_activation_rate(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    aip = CMP.AirProperties(FT)
    ap = CMP.AerosolActivationParameters(FT)
    pa = CMP.PrescribedAerosol(FT)
    ad = AM.aerosol_distribution(pa)

    ρₐ = FT(0.9)
    T = FT(285)
    p = FT(9e4)
    q_vs = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρₐ)
    x_seed = CMP.Microphysics2MParams(FT).warm_rain.seifert_beheng.pdf_c.xc_min

    # A cell held at a stated supersaturation over liquid, with no condensate.
    state_at(S, w) = (
        T, p, w, ρₐ, q_vs * (1 + S), zero(FT), zero(FT), zero(FT), FT(x_seed),
    )
    rate_at(S, w) = AA.cloud_droplet_activation_rate(ap, pa, aip, tps, state_at(S, w)...)

    TT.@testset "Activated number at a given supersaturation" begin
        # The activated number is monotone in the supersaturation and bounded by the aerosol
        # budget, since the sum over modes cannot exceed the particles that exist.
        N_a = pa.N_accum + pa.N_coarse
        Ns = [AA.total_N_activated(ap, ad, T, FT(S)) for S in (1e-4, 1e-3, 1e-2, 1e-1, 1)]
        TT.@test issorted(Ns)
        TT.@test all(N -> 0 <= N <= N_a * (1 + sqrt(eps(FT))), Ns)
        # At and below the supersaturation floor the activated number is a vanishing fraction of
        # the budget: the lognormal argument is deep in the tail where `erf` saturates.
        TT.@test AA.total_N_activated(ap, ad, T, zero(FT)) < FT(1e-6) * N_a
        # The closed-form derivative matches a centered difference on the smooth interior. The
        # step is wide enough that the difference of two nearly-saturated tails is resolved in
        # `Float32` as well.
        S = FT(2e-3)
        h = FT(1e-4)
        fd =
            (AA.total_N_activated(ap, ad, T, S + h) - AA.total_N_activated(ap, ad, T, S - h)) /
            (2 * h)
        TT.@test AA.∂N_activated_∂S(ap, ad, T, S) ≈ fd rtol = FT(2e-2)
    end

    TT.@testset "The ambient branch does not need an updraft" begin
        # Supersaturated air at rest still activates: the ambient supersaturation is the term
        # that fires, and the parcel branch is merely switched off.
        act = rate_at(FT(0.2), zero(FT))
        TT.@test act.∂ₜn_lcl > 0
        TT.@test act.∂ₜq_lcl ≈ x_seed * act.∂ₜn_lcl
        TT.@test act.inv_τ_act > 0
        TT.@test act.outside_parcel_regime
        # Subsaturated air does not activate at any updraft.
        for w in (FT(0), FT(1))
            quiet = rate_at(FT(-0.05), w)
            TT.@test quiet.∂ₜn_lcl == 0
            TT.@test quiet.∂ₜq_lcl == 0
            TT.@test quiet.inv_τ_act == 0
            TT.@test quiet.∂ₜn_∂S == 0
        end
    end

    TT.@testset "The parcel branch is selected inside its own regime" begin
        # At a marginally supersaturated state a resolved updraft generates more supersaturation
        # than the cell carries, so the parcel maximum is the term that governs and the cell is
        # reported as inside the regime the ARG closure is valid in.
        marginal = rate_at(FT(1e-6), FT(1))
        TT.@test !marginal.outside_parcel_regime
        TT.@test marginal.∂ₜn_lcl > 0
        # Every state stays finite, including the corners that make `max_supersaturation`
        # itself non-total: no updraft, subsiding air, no pressure.
        for (S, w, pp) in (
            (FT(0.2), FT(0), FT(9e4)), (FT(0.2), FT(-1), FT(9e4)), (FT(0.2), FT(1), FT(0)),
            (FT(0), FT(0), FT(0)),
        )
            act = AA.cloud_droplet_activation_rate(
                ap, pa, aip, tps, T, pp, w, ρₐ, q_vs * (1 + S),
                zero(FT), zero(FT), zero(FT), FT(x_seed),
            )
            TT.@test isfinite(act.∂ₜn_lcl)
            TT.@test isfinite(act.∂ₜq_lcl)
            TT.@test isfinite(act.inv_τ_act)
            TT.@test isfinite(act.∂ₜn_∂S)
        end
    end

    TT.@testset "An aerosol population with no particles is inert" begin
        empty_pa = CMP.PrescribedAerosol(FT; N_accum = zero(FT), N_coarse = zero(FT))
        act = AA.cloud_droplet_activation_rate(
            ap, empty_pa, aip, tps, T, p, FT(1), ρₐ, q_vs * FT(1.2),
            zero(FT), zero(FT), zero(FT), FT(x_seed),
        )
        TT.@test act.∂ₜn_lcl == 0
        TT.@test act.∂ₜq_lcl == 0
    end
end

function test_absent_aerosol(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    aip = CMP.AirProperties(FT)
    ap = CMP.AerosolActivationParameters(FT)

    TT.@testset "An aerosol-absent parameter bundle constructs" begin
        # The default parameter bundle carries no aerosol population, and building it must not
        # require any prescribed-aerosol parameter to exist.
        for with_ice in (false, true)
            mp = CMP.Microphysics2MParams(FT; with_ice)
            TT.@test isnothing(mp.warm_rain.aerosol)
            TT.@test mp.warm_rain.activation isa CMP.AerosolActivationParameters{FT}
        end
        # A configuration that states its aerosol population passes it in.
        pa = CMP.PrescribedAerosol(FT)
        TT.@test CMP.Microphysics2MParams(FT; aerosol = pa).warm_rain.aerosol == pa
    end

    TT.@testset "An aerosol-absent bundle computes zero activation" begin
        ρₐ = FT(0.9)
        T = FT(285)
        q_vs = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρₐ)
        # Strongly supersaturated, rising air: every branch of the populated rate would fire.
        act = AA.cloud_droplet_activation_rate(
            ap, nothing, aip, tps, T, FT(9e4), FT(1), ρₐ, q_vs * FT(1.2),
            zero(FT), zero(FT), zero(FT), FT(1e-15),
        )
        TT.@test act.∂ₜn_lcl == 0
        TT.@test act.∂ₜq_lcl == 0
        TT.@test act.inv_τ_act == 0
        TT.@test act.∂ₜn_∂S == 0
        TT.@test !act.outside_parcel_regime
        TT.@test act.qᵥ_sat ≈ q_vs
    end
end

TT.@testset "Prescribed aerosol and droplet activation rate" begin
    for FT in (Float32, Float64)
        test_prescribed_aerosol(FT)
        test_activation_rate(FT)
        test_absent_aerosol(FT)
    end
end
nothing
