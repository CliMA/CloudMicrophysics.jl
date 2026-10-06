#####
##### Process terms
#####

"""
    Transfer{Donor, Receiver, FT}
    Transfer(:Donor => :Receiver, S)

One-way transfer of mass from species `Donor` to species `Receiver`.

For example, autoconversion transfers mass from cloud liquid to rain.
`Donor` and `Receiver` are species names, such as `:q_lcl`.

# Fields
- `S`: transfer rate from `Donor` to `Receiver`, `S ≥ 0` [kg/kg/s].

# Examples
```julia
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
BMT.Transfer(:q_lcl => :q_rai, 1e-6)
```
"""
struct Transfer{Donor, Receiver, FT}
    S::FT
end
@inline Transfer(p::Pair{Symbol, Symbol}, S) = Transfer{p.first, p.second, typeof(S)}(S)

"""
    VaporExchange{Condensate, FT}
    VaporExchange(:Condensate, S)

Exchange of mass between the vapor and the condensate species `Condensate`.

The exchange goes in either direction, for example by deposition on snow and sublimation of
snow. `Condensate` is a species name, such as `:q_sno`.

# Fields
- `S`: exchange rate, positive from the vapor to `Condensate` [kg/kg/s].

# Examples
```julia
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
BMT.VaporExchange(:q_sno, 1e-7)  # deposition on snow
```
"""
struct VaporExchange{Condensate, FT}
    S::FT
end
@inline VaporExchange(species::Symbol, S) = VaporExchange{species, typeof(S)}(S)

"""
    VaporRelaxation{Condensate, FT}
    VaporRelaxation(:Condensate, S, τ)

Exchange of mass between the vapor and the condensate species `Condensate` that relaxes toward
equilibrium.

The exchange goes in either direction and relaxes with the timescale `τ`, for example by
condensation and evaporation of cloud liquid. `Condensate` is a species name, such as `:q_lcl`.

# Fields
- `S`: exchange rate, positive from the vapor to `Condensate` [kg/kg/s].
- `τ`: relaxation timescale [s].

# Examples
```julia
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
BMT.VaporRelaxation(:q_lcl, -2e-7, 30.0)  # evaporation of cloud liquid
```
"""
struct VaporRelaxation{Condensate, FT}
    S::FT
    τ::FT
end
@inline VaporRelaxation(species::Symbol, S, τ) = VaporRelaxation{species, typeof(S)}(S, τ)

"""
    species_tendency(terms, Q)

Sum the tendencies of `terms` into a state of the type `Q` [kg/kg/s].

`Q` is a state type with one field per species, for example `Condensates1M{FT}`. Each term adds its rate `S`:
- `Transfer(:Donor => :Receiver, S)`: `-S` to `Donor` and `S` to `Receiver`.
- `VaporExchange(:Condensate, S)` and `VaporRelaxation(:Condensate, S, τ)`: `S` to `Condensate`.

The vapor is not a species of `Q`; it changes by minus the sum of the tendencies.
"""
@inline species_tendency(terms, ::Type{Q}) where {Q} =
    UU.unrolled_reduce(_add_tendency, values(terms), zero(Q))

@inline function _add_tendency(dq_dt, t::Transfer{Donor, Receiver}) where {Donor, Receiver}
    dq_dt = Base.setindex(dq_dt, dq_dt[Donor] - t.S, Donor)
    return Base.setindex(dq_dt, dq_dt[Receiver] + t.S, Receiver)
end
@inline _add_tendency(
    dq_dt, t::Union{VaporExchange{Condensate}, VaporRelaxation{Condensate}},
) where {Condensate} = Base.setindex(dq_dt, dq_dt[Condensate] + t.S, Condensate)
