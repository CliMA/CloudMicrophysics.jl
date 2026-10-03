# # Getting started

# This guide shows how to call a function from `CloudMicrophysics.jl` package.
# Please consult our [README](https://github.com/CliMA/CloudMicrophysics.jl?tab=readme-ov-file#installation-and-running-instructions)
# for the CloudMicrophysics.jl installation instructions.

# In this guide will call the `accretion` function that parameterizes the growth of rain drops
# through collisions with cloud droples.
# Check the [API documentation](https://clima.github.io/CloudMicrophysics.jl/dev/API/#CloudMicrophysics.Microphysics2M.accretion)
# and the [parameterization documentation](https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics2M/#Accretion)
# for more details.

# We start by defining the single precision floating point type
# that will be used in the computations.
# We import the `Microphysics2M` module in which the `accretion` function is defined
# and the `Parameters` module in which we store the default values of free parameters.
FT = Float32

import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Parameters as CMP

nothing #hide

# We grab the parameters needed by the accretion function from the parameters module
# and define the example input values.
# Note that both the free parameters and the input values are of the same floating point type.
# All values are defined in base SI units.
const SB2006 = CMP.SB2006(FT)
qₗ = FT(1e-3)  # Cloud liquid water specific humidity
qᵣ = FT(5e-4)  # Rain water specific humidity
ρₐ = FT(1)     # Air density
Nₗ = FT(1e8)   # Cloud droplet number concentration

nothing #hide

# Finally, we call `accretion`, which will return the accretion rates for
# cloud and rain water specific humidities, as well as cloud and rain water number concentrations.
(; dq_rai_dt, dq_liq_dt, dN_rai_dt, dN_liq_dt) =
    CM2.accretion(SB2006, qₗ, qᵣ, ρₐ, Nₗ)
@info("Accretion rates: ", dq_rai_dt, dq_liq_dt, dN_rai_dt, dN_liq_dt)
