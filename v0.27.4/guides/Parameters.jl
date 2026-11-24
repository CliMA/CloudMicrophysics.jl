# # Handling parameters

# This guide shows how to overwrite the default free parameters values and
# how to dispatch between different parameterization options based on parameters types.

# ## Overwriting parameters

# `CloudMicrophysics.jl` is designed to allow easy parameter calibrations.
# As a result, free parameters are not hard-coded in the source code but are instead
# passed as arguments to functions. The default values are stored in a separate
# repository [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl) in a `toml` file.

# We start by importing the `ClimaParams` package and the needed
# `CloudMicrophysics.jl` modules. We define the precision type.

import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.HetIceNucleation as CMI_het

FT = Float32

nothing #hide

# We define the default parameters struct using the default `toml` file via the constructor
# provided in the `CloudMicrophysics.jl` `Parameters` module.
# Additionally, we create a `toml` file in which we save the name, value and type
# of the parameter we are changing. In a typical application, this step would be done
# via the calibration algorithm searching for optimal parameters values.
# We create a second struct with changed parameters based on the `override_dict.toml` file.
const default = CMP.Rain(FT)

override_file = joinpath("override_dict.toml")
open(override_file, "w") do io
    println(io, "[rain_autoconversion_timescale]")
    println(io, "value = " * string(1))
    println(io, "type = \"float\"")
end
toml_dict = CP.create_toml_dict(FT; override_file)
isfile(override_file) && rm(override_file; force = true)
const overwrite = CMP.Rain(toml_dict)

nothing #hide

# Overwriting the parameters can also be done using dictionaries instead of toml files.
# As an example we create a dictionary were we define a different value of the
# rain autoconversion time scale and create another parameter struct based on it.
override_file = Dict(
    "rain_autoconversion_timescale" =>
        Dict("value" => 13, "type" => "float"),
)
toml_dict2 = CP.create_toml_dict(FT; override_file)
const overwrite2 = CMP.Rain(toml_dict2)

nothing #hide

# Finally we check the values of the rain autoconversion timescales
# and the corresponding rain formation rates.
qₗ = FT(1e-3)
default_acnv = CM1.conv_q_liq_to_q_rai(default.acnv1M, qₗ) # Rain autoconversion rate
overwrite_acnv = CM1.conv_q_liq_to_q_rai(overwrite.acnv1M, qₗ) # Rain autoconversion rate
overwrite_acnv2 = CM1.conv_q_liq_to_q_rai(overwrite2.acnv1M, qₗ) # Rain autoconversion rate

@info("Default:", default.acnv1M.τ, default_acnv)
@info("Overwrite:", overwrite.acnv1M.τ, overwrite_acnv)
@info("Overwrite from dict:", overwrite2.acnv1M.τ, overwrite_acnv2)

# ## Dispatching over parameter types

# `CloudMicrophysics.jl` `Parameters` module introduces type hierarchy that is used to
# dispatch over different parameterization options and inputs.
# For example `ABIFM_J` accepts `Illite` or `Kaolinite`
# as input type and will return immersion freezing nucleation rate coefficient
# that corresponds to the two different aerosol types.
const aerosol1 = CMP.Illite(FT)
const aerosol2 = CMP.Kaolinite(FT)

Δa_w = FT(0.28)
J1 = CMI_het.ABIFM_J(aerosol1, Δa_w)
J2 = CMI_het.ABIFM_J(aerosol2, Δa_w)

@info("ABIFM derived J: ", J1, J2)

# Similarily, the 2 moment microphysics offers different rain autoconversion
# formulations based on parameterizations from the literature.
# We can chose between them based on the free parameters type that is passed in as argument.

const KK2000 = CMP.KK2000(FT)  # Khairoutdinov and Kogan (2000)
const B1994 = CMP.B1994(FT)    # Beheng (1994)
const TC1980 = CMP.TC1980(FT)  # Tripoli and Cotton (1980)
const LD2004 = CMP.LD2004(FT)  # Liu and Daum (2004)

qₗ = FT(1e-3)
ρₐ = FT(1)
N = FT(1e8)
KK2000_rate = CM2.conv_q_liq_to_q_rai(KK2000, qₗ, ρₐ, N)
TC1980_rate = CM2.conv_q_liq_to_q_rai(TC1980, qₗ, ρₐ, N)
LD2004_rate = CM2.conv_q_liq_to_q_rai(LD2004, qₗ, ρₐ, N)
B1994_rate = CM2.conv_q_liq_to_q_rai(B1994, qₗ, ρₐ, N)

@info("Autoconversion: ", KK2000_rate, B1994_rate, TC1980_rate, LD2004_rate)
