export AerosolActivationParameters, PrescribedAerosol

"""
    AerosolActivationParameters{FT}

Parameters for Abdul-Razzak and Ghan 2000 aerosol activation scheme
DOI: 10.1029/1999JD901161

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AerosolActivationParameters{FT} <: ParametersType
    "molar mass of water [kg/mol]"
    M_w::FT
    "gas constant [J/mol/K]"
    R::FT
    "cloud water density [kg/m3]"
    ρ_w::FT
    "cloud ice density [kg/m3]"
    ρ_i::FT
    "surface tension of water [N/m]"
    σ::FT
    "gravitational acceleration [m/s2]"
    g::FT
    "scaling coefficient in Abdul-Razzak and Ghan 2000 [-]"
    f1::FT
    "scaling coefficient in Abdul-Razzak and Ghan 2000 [-]"
    f2::FT
    "scaling coefficient in Abdul-Razzak and Ghan 2000 [-]"
    g1::FT
    "scaling coefficient in Abdul-Razzak and Ghan 2000 [-]"
    g2::FT
    "power of (zeta / eta) in Abdul-Razzak and Ghan 2000 [-]"
    p1::FT
    "power of (S_m^2 / (zeta + 3 * eta)) in Abdul-Razzak and Ghan 2000 [-]"
    p2::FT
end

function AerosolActivationParameters(td::CP.ParamDict)
    name_map = (;
        :molar_mass_water => :M_w,
        :universal_gas_constant => :R,
        :density_liquid_water => :ρ_w,
        :density_ice_water => :ρ_i,
        :surface_tension_water => :σ,
        :gravitational_acceleration => :g,
        :ARG2000_f_coeff_1 => :f1,
        :ARG2000_f_coeff_2 => :f2,
        :ARG2000_g_coeff_1 => :g1,
        :ARG2000_g_coeff_2 => :g2,
        :ARG2000_pow_1 => :p1,
        :ARG2000_pow_2 => :p2,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return AerosolActivationParameters(; parameters...)
end

"""
    PrescribedAerosol{FT}

A prescribed, spatially uniform two-mode lognormal aerosol population, in the form the
Abdul-Razzak and Ghan activation scheme consumes: an accumulation mode and a coarse mode, each a
single-component `AerosolModel.Mode_κ`.

Prescribed rather than prognostic because the aerosol is a boundary condition on the
microphysics, not a state of it. Carrying it as a parameter is what lets the droplet number stay
fully prognostic while the SUPPLY of cloud condensation nuclei is a stated, per-configuration
choice: the RCEMIP-II protocol, for instance, asks schemes that specify cloud condensation nuclei
to set values consistent with `N_c = 1e8` per cubic metre, which is a statement about this struct.

The values below are code defaults, not ClimaParams keys, because the population is a property of
the configuration being run rather than of the parameterization. They are reachable three ways, in
increasing precedence: the defaults here, the TOML names listed at
[`PRESCRIBED_AEROSOL_TOML_NAMES`](@ref) if an override file supplies them, and keyword arguments to
the constructor.

Every field is a plain scalar so the struct is `isbits` and the distribution can be rebuilt inside
a GPU kernel; see
[`AerosolModel.aerosol_distribution`](@ref CloudMicrophysics.AerosolModel.aerosol_distribution).

# Constraints, not checked at construction because the struct is built inside a kernel
`r_dry` and `κ` must be strictly positive, or the critical supersaturation is infinite and the
activated number is undefined. `stdev` must be strictly greater than one, or the lognormal
collapses to a monodisperse population and `log(stdev)` divides by zero. `N` may be zero, which
switches the mode off cleanly.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct PrescribedAerosol{FT} <: ParametersType
    "geometric mean dry radius of the accumulation mode [m]. Sulfate-like, the same size as the
    MERRA2 sulfate aerosol radius."
    r_dry_accum::FT = 3.5e-7
    "geometric standard deviation of the accumulation mode [-]. The MAM3 accumulation-mode width."
    stdev_accum::FT = 1.8
    "number concentration of the accumulation mode [1/m3]. Together with `N_coarse` this sets the
    supply of cloud condensation nuclei, and it is the least constrained number in the struct.
    TODO: the value is a plausible continental-to-maritime midpoint rather than a measured or
    fitted one, and a configuration that cares about droplet number should state its own instead
    of inheriting this."
    N_accum::FT = 9.0e6
    "hygroscopicity parameter of the accumulation mode [-]. Sulfate. Source: Petters and
    Kreidenweis (2007), DOI: 10.5194/acp-7-1961-2007."
    κ_accum::FT = 0.53
    "geometric mean dry radius of the coarse mode [m]. Sea-salt-like."
    r_dry_coarse::FT = 1.0e-6
    "geometric standard deviation of the coarse mode [-]. The MAM3 coarse-mode width."
    stdev_coarse::FT = 1.8
    "number concentration of the coarse mode [1/m3]. TODO: as uncertain as `N_accum`, and stated
    per configuration for the same reason."
    N_coarse::FT = 1.0e6
    "hygroscopicity parameter of the coarse mode [-]. Sea salt. Source: Petters and Kreidenweis
    (2007), DOI: 10.5194/acp-7-1961-2007."
    κ_coarse::FT = 1.12
end

"""
    PRESCRIBED_AEROSOL_TOML_NAMES

The TOML parameter name each [`PrescribedAerosol`](@ref) field is read from when an override file
supplies it, as a `NamedTuple` mapping field name to TOML name.

None of these names is a ClimaParams default, so a stock parameter dictionary carries none of them
and the code defaults stand. Listing them here keeps the override path available without making the
values a central default: a configuration that states its own aerosol population adds any subset of
these names to its override file.
"""
const PRESCRIBED_AEROSOL_TOML_NAMES = (;
    r_dry_accum = "prescribed_aerosol_accumulation_radius",
    stdev_accum = "prescribed_aerosol_accumulation_stdev",
    N_accum = "prescribed_aerosol_accumulation_number",
    κ_accum = "prescribed_aerosol_accumulation_kappa",
    r_dry_coarse = "prescribed_aerosol_coarse_radius",
    stdev_coarse = "prescribed_aerosol_coarse_stdev",
    N_coarse = "prescribed_aerosol_coarse_number",
    κ_coarse = "prescribed_aerosol_coarse_kappa",
)

"""
    PrescribedAerosol(td::CP.ParamDict; kwargs...)

Build a [`PrescribedAerosol`](@ref) in the dictionary's float type, taking each field from the
TOML name in [`PRESCRIBED_AEROSOL_TOML_NAMES`](@ref) where the dictionary carries it and from the
struct's own default where it does not. Keyword arguments name fields directly and take precedence
over both.
"""
function PrescribedAerosol(td::CP.ParamDict; kwargs...)
    FT = CP.float_type(td)
    name_map = Dict(
        Symbol(toml_name) => field for
        (field, toml_name) in pairs(PRESCRIBED_AEROSOL_TOML_NAMES) if haskey(td.data, toml_name)
    )
    isempty(name_map) && return PrescribedAerosol{FT}(; kwargs...)
    overrides = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return PrescribedAerosol{FT}(; overrides..., kwargs...)
end
