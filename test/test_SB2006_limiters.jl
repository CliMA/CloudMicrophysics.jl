import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics2M as CM2

import QuadGK as QGK
import SpecialFunctions as SF

FT = Float32

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
# old parameters
#SB2006 = CMP.SB2006(FT)

# air and liquid water densities
ρₐ = FT(1.2)
ρₗ = SB2006.pdf_r.ρw

# example precipitation number concentration and specific humidity
Nᵣ = FT(1e8)
qᵣ = FT(5e-4)

# mass of liquid droplet as a function of its diameter
function m(D)
    return π / 6 * ρₗ* D^3
end

# rain drop size distribution as a function of diameter (eq.(3) from 2M docs)
function f_D(D)
    β = (π * ρₗ * Nᵣ / (qᵣ * ρₐ))^(1 / 3)
    α = Nᵣ * β
    return α * exp(-β * D)
end

# rain drop mass distribution (eq.(4) from 2M docs)
function f_x(x)
    Br = (SF.gamma(1) / SF.gamma(4) * qᵣ * ρₐ / Nᵣ)^(-1/3)
    Ar = 1/3 * Nᵣ * Br
    return Ar * x^(-2/3) * exp(-Br * x^(1/3))
end

# rain drop mass distribution, but using the SB2006 limiters
function f_x_limited(x)
    (; Ar, Br) = CM2.raindrops_limited_vars(SB2006.pdf_r, qᵣ, ρₐ, Nᵣ)
    return Ar * x^(-2/3) * exp(-Br * x^(1/3))
end

# integral bounds
D₀ = 1e-7
D∞ = 1e-2
m₀ = m(D₀)
m∞ = m(D∞)

# Sanity checks for number concentrations
ND     = QGK.quadgk(x -> f_D(x),         D₀, D∞)[1]
Nx     = QGK.quadgk(x -> f_x(x),         m₀, m∞)[1]
Nx_lim = QGK.quadgk(x -> f_x_limited(x), m₀, m∞)[1]
@info(" ", Nᵣ, ND, Nx, Nx_lim)

# Sanity checks for specific humidities
qD =     QGK.quadgk(x -> m(x) * f_D(x),      D₀, m∞)[1] / ρₐ
qx =     QGK.quadgk(x -> x * f_x(x),         m₀, m∞)[1] / ρₐ
qx_lim = QGK.quadgk(x -> x * f_x_limited(x), m₀, m∞)[1] / ρₐ
@info(" ", qᵣ, qD, qx, qx_lim)

