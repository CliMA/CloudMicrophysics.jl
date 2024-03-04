import Thermodynamics as TD
import CloudMicrophysics.Common as CMO
import CloudMicrophysics.HetIceNucleation as CMI_het
import CloudMicrophysics.HomIceNucleation as CMI_hom
import CloudMicrophysics.Parameters as CMP

function deposition_nucleation(::Empty, state, dY)
    FT = eltype(state)
    return FT(0)
end

function deposition_nucleation(params::MohlerAF, state, dY)
    (; ips, aerosol, const_dt, tps) = params
    (; Sₗ, T, Nₐ, Nᵢ) = state
    Sᵢ = ξ(tps, T) * Sₗ
    FT = eltype(state)
    if Sᵢ >= ips.deposition.Sᵢ_max
        @warn("Supersaturation exceeds Sᵢ_max. No dust will be activated.")
    end

    AF =
        Sᵢ >= ips.deposition.Sᵢ_max ? FT(0) :
        AF = CMI_het.dust_activated_number_fraction(
            aerosol,
            ips.deposition,
            Sᵢ,
            T,
        )
    return max(FT(0), AF * Nₐ - Nᵢ) / const_dt
end

function deposition_nucleation(params::MohlerRate, state, dY)
    (; ips, aerosol, tps) = params
    (; T, Nₐ, Sₗ) = state
    Sᵢ = ξ(tps, T) * Sₗ
    dSᵢdt = ξ(tps, T) * dY[1]

    if Sᵢ >= ips.deposition.Sᵢ_max
        @warn("Supersaturation exceeds Sᵢ_max. No dust will be activated.")
    end

    return Sᵢ >= ips.deposition.Sᵢ_max ? FT(0) :
           CMI_het.MohlerDepositionRate(
        aerosol,
        ips.deposition,
        Sᵢ,
        T,
        dSᵢdt,
        Nₐ,
    )
end

function deposition_nucleation(params::ABDINM, state, dY)
    FT = eltype(state)
    (; tps, aerosol, r_nuc) = params
    (; T, p_air, qᵥ, qₗ, qᵢ) = state

    q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
    Rᵥ = TD.Parameters.R_v(tps)
    R_air = TD.gas_constant_air(tps, q)
    e = eᵥ(qᵥ, p_air, R_air, Rᵥ)

    Δa_w = CMO.a_w_eT(tps, e, T) - CMO.a_w_ice(tps, T)
    J = CMI_het.deposition_J(aerosol, Δa_w)
    A = 4 * FT(π) * r_nuc^2
    return max(FT(0), J * Nₐ * A)
end

function deposition_nucleation(params::P3_dep, state, dY)
    FT = eltype(state)
    (; ips, const_dt) = params
    (; T, Nᵢ) = state
    Nᵢ_depo = CMI_het.P3_deposition_N_i(ips.p3, T)
    return max(FT(0), (Nᵢ_depo - Nᵢ) / const_dt)
end

function immersion_freezing(::Empty, PSD, state)
    FT = eltype(state)
    return FT(0)
end

function immersion_freezing(params::ABIFM, PSD, state)
    (; T, xS, p_air, qᵥ, qₗ, qᵢ, Nₗ) = state
    (; H₂SO₄ps, tps, aerosol, A_aer) = params

    q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
    Rᵥ = TD.Parameters.R_v(tps)
    R_air = TD.gas_constant_air(tps, q)
    e = eᵥ(qᵥ, p_air, R_air, Rᵥ)

    # TODO - get rif of a_w_x option
    Δa_w =
        T > FT(185) && T < FT(235) ?
        CMO.a_w_xT(H₂SO₄ps, tps, xS, T) - CMO.a_w_ice(tps, T) :
        CMO.a_w_eT(tps, e, T) - CMO.a_w_ice(tps, T)

    J = CMI_het.ABIFM_J(aerosol, Δa_w)
    return max(FT(0), J * Nₗ * A_aer)
end

function immersion_freezing(params::P3_het, PSD, state)
    FT = eltype(state)
    (; const_dt, ips) = params
    (; T, Nₗ, Nᵢ) = state
    Nᵢ_het = CMI_het.P3_het_N_i(ips.p3, T, Nₗ, PSD.Vₗ, const_dt)
    return max(FT(0), (Nᵢ_het - Nᵢ) / const_dt)
end

function homogeneous_freezing(::Empty, PSD, state)
    FT = eltype(state)
    return FT(0)
end

function homogeneous_freezing(params::ABHOM, PSD, state)
    FT = eltype(state)
    (; tps, ips) = params
    (; T, p_air, qᵥ, qₗ, qᵢ, Nₐ, Nₗ) = state

    q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
    Rᵥ = TD.Parameters.R_v(tps)
    R_air = TD.gas_constant_air(tps, q)
    e = eᵥ(qᵥ, p_air, R_air, Rᵥ)

    Δa_w = CMO.a_w_eT(tps, e, T) - CMO.a_w_ice(tps, T)

    if Δa_w > ips.homogeneous.Δa_w_max || Δa_w < ips.homogeneous.Δa_w_min
        @warn("Clipping Δa_w for Homogeneous freezing")
    end

    Δa_w = max(min(ips.homogeneous.Δa_w_max, Δa_w), ips.homogeneous.Δa_w_min)
    J = CMI_hom.homogeneous_J(ips.homogeneous, Δa_w)

    return max(FT(0), J * Nₗ * PSD.Vₗ)
end

function homogeneous_freezing(params::P3_hom, PSD, state)
    FT = eltype(state)
    (; Nₗ, T) = state
    (; const_dt) = params
    #TODO - use ClimaParams
    return T < FT(233.15) && Nₗ > FT(0) ? Nₗ / const_dt : FT(0)
end

function condensation(::Empty, PSD, state, ρ_air)
    FT = eltype(state)
    return FT(0)
end

function condensation(params::CondParams, PSD, state, ρ_air)
    (; aps, tps) = params
    (; Sₗ, T, Nₗ) = state
    Gₗ = CMO.G_func(aps, tps, T, TD.Liquid())
    return 4 * FT(π) / ρ_air * (Sₗ - 1) * Gₗ * PSD.rₗ * Nₗ
end

function deposition(::Empty, PSD, state, ρ_air)
    FT = eltype(state)
    return FT(0)
end

function deposition(params::DepParams, PSD, state, ρ_air)
    (; aps, tps) = params
    (; T, Sₗ, Nᵢ) = state
    FT = eltype(state)
    Sᵢ = ξ(tps, T) * Sₗ
    Gᵢ = CMO.G_func(aps, tps, T, TD.Ice())
    return 4 * FT(π) / ρ_air * (Sᵢ - 1) * Gᵢ * PSD.rᵢ * Nᵢ
end
