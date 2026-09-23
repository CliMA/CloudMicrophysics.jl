#=
Single place for the mapping between ClimaParams parameter names and the fields of the CloudMicrophysics
parameter structs (and the parameters of the 1-moment process options).

`name_map(T)` is the `(; :climaparams_name => :local_name, ...)` NamedTuple that the `CP.ParamDict`
constructor of `T` reads. For structs whose constructor is nothing but "read the map and fill the
fields", the constructor is generated below from the map; constructors that validate, derive or
compose values are written by hand next to their struct and read their parameters with
`make_params(td, name_map(T))`, so that every ClimaParams name the library uses is listed here.
`name_map_summary()` collects all maps for the documentation and the tests.
=#

export name_map, make_params, name_map_summary, parameter_groups

"""
    name_map(::Type{T})
    name_map(::Type{ParticleMass}, ::Type{Rain})    # and the other type-dispatched constructors

Return the NamedTuple `(; :climaparams_name => :local_name, ...)` of the ClimaParams parameters
that the `CP.ParamDict` constructor of `T` reads. The local names are the struct fields, or the
inputs of derived fields (e.g. the power-law exponents that enter the gamma-function factors of
`Blk1MVelTypeRain`). For a 1-moment process option `T`, the map of the parameters returned by
`process_params_for(T(), td)`.
"""
function name_map end

"""
    make_params(td::CP.ParamDict, nm::NamedTuple)

Read the parameters of the name map `nm` (see [`name_map`](@ref)) from `td`, tag them as used by
CloudMicrophysics, and return them as a NamedTuple keyed by the local names, converted to the
float type of `td`. The field order of the result is not the order of `nm`, so construct structs
by keyword (`T(; make_params(td, name_map(T))...)`), never positionally.
"""
make_params(td::CP.ParamDict, nm::NamedTuple) = CP.get_parameter_values(td, nm, "CloudMicrophysics")

# ═══════════════════════════════════════════════════════════════════
# Air and water properties
# ═══════════════════════════════════════════════════════════════════

name_map(::Type{AirProperties}) = (;
    :thermal_conductivity_of_air => :K_therm,
    :diffusivity_of_water_vapor => :D_vapor,
    :kinematic_viscosity_of_air => :ν_air,
)

name_map(::Type{WaterProperties}) = (; :density_liquid_water => :ρw, :density_ice_water => :ρi)

# ═══════════════════════════════════════════════════════════════════
# 0-moment microphysics
# ═══════════════════════════════════════════════════════════════════

name_map(::Type{Parameters0M}) = (;
    :precipitation_timescale => :τ_precip,
    :specific_humidity_precipitation_threshold => :qc_0,
    :supersaturation_precipitation_threshold => :S_0,
)

# ═══════════════════════════════════════════════════════════════════
# 1-moment microphysics
# ═══════════════════════════════════════════════════════════════════

name_map(::Type{CloudLiquid}) = (;
    :density_liquid_water => :ρw,
    :liquid_cloud_effective_radius => :r_eff,
    :cloud_liquid_sedimentation_number_concentration => :N_0,
)
name_map(::Type{CloudIce}) = (;
    :cloud_ice_apparent_density => :ρᵢ,
    :cloud_ice_size_distribution_coefficient_n0 => :n0,
    :ice_cloud_effective_radius => :r_eff,
    :cloud_ice_sedimentation_number_concentration => :N_0,
)
name_map(::Type{IceNumberTemperatureFit}) = (;
    :cloud_ice_number_temperature_fit_prefactor => :N_ref,
    :cloud_ice_number_temperature_fit_intercept => :a,
    :cloud_ice_number_temperature_fit_slope => :b,
    :cloud_ice_number_max => :N_max,
    :temperature_water_freeze => :T_freeze,
)
name_map(::Type{ParticleMass}, ::Type{CloudIce}) = (;
    :cloud_ice_apparent_density => :ρᵢ,
    :cloud_ice_crystals_length_scale => :r0,
    :cloud_ice_mass_size_relation_coefficient_me => :me,
    :cloud_ice_mass_size_relation_coefficient_delm => :Δm,
    :cloud_ice_mass_size_relation_coefficient_chim => :χm,
)
name_map(::Type{Rain}) = (;
    :rain_drop_size_distribution_coefficient_n0 => :n0,
    :rain_ventilation_coefficient_a => :a,
    :rain_ventilation_coefficient_b => :b,
)
name_map(::Type{ParticleMass}, ::Type{Rain}) = (;
    :density_liquid_water => :ρ,
    :rain_drop_length_scale => :r0,
    :rain_mass_size_relation_coefficient_me => :me,
    :rain_mass_size_relation_coefficient_delm => :Δm,
    :rain_mass_size_relation_coefficient_chim => :χm,
)
name_map(::Type{ParticleArea}, ::Type{Rain}) = (;
    :rain_drop_length_scale => :r0,
    :rain_cross_section_size_relation_coefficient_ae => :ae,
    :rain_cross_section_size_relation_coefficient_dela => :Δa,
    :rain_cross_section_size_relation_coefficient_chia => :χa,
)
name_map(::Type{Snow}) = (;
    :snow_apparent_density => :ρᵢ,
    :density_ice_water => :ρᵢ_bulk,
    :snow_flake_size_distribution_coefficient_mu => :μ,
    :snow_flake_size_distribution_coefficient_nu => :ν,
    :snow_ventilation_coefficient_a => :a,
    :snow_ventilation_coefficient_b => :b,
    :snow_aspect_ratio => :ϕ,
    :snow_aspect_ratio_coefficient => :κ,
)
name_map(::Type{ParticleMass}, ::Type{Snow}) = (;
    :snow_flake_length_scale => :r0,
    :snow_mass_size_relation_coefficient_me => :me,
    :snow_mass_size_relation_coefficient_delm => :Δm,
    :snow_mass_size_relation_coefficient_chim => :χm,
)
name_map(::Type{ParticleArea}, ::Type{Snow}) = (;
    :snow_flake_length_scale => :r0,
    :snow_cross_section_size_relation_coefficient => :ae,
    :snow_cross_section_size_relation_coefficient_dela => :Δa,
    :snow_cross_section_size_relation_coefficient_chia => :χa,
)
name_map(::Type{VarTimescaleAcnv}) = (;
    :rain_autoconversion_timescale => :τ,
    :Variable_time_scale_autoconversion_coeff_alpha => :α,
    :prescribed_cloud_droplet_number_concentration => :Nc,
)
# The convective ("fast") values use the original Kessler keys; the quiescent ("slow") values
# have their own `_stratiform` keys (ClimaParams >= 1.1.12), equal to the fast ones by default.
name_map(::Type{KesslerAcnv}) = (;
    :rain_autoconversion_timescale_stratiform => :τ_slow,
    :rain_autoconversion_timescale => :τ_fast,
    :cloud_liquid_water_specific_humidity_autoconversion_threshold_stratiform => :q_threshold_slow,
    :cloud_liquid_water_specific_humidity_autoconversion_threshold => :q_threshold_fast,
    :rain_autoconversion_velocity_scale => :w_0,
    :threshold_smooth_transition_steepness => :k,
)

# ═══════════════════════════════════════════════════════════════════
# 1-moment process options
# ═══════════════════════════════════════════════════════════════════

name_map(::Type{TemperatureDependent}) = (; :sublimation_deposition_timescale => :τ_relax)
name_map(::Type{NoSupersaturation}) = (;
    :snow_autoconversion_timescale => :τ,
    :cloud_ice_specific_humidity_autoconversion_threshold => :q_threshold,
    :threshold_smooth_transition_steepness => :k,
)

name_map(::Type{CloudLiquidFormation}) = (;
    :condensation_evaporation_timescale => :τ_relax,
    :temperature_homogenous_nucleation => :T_hom,
)

name_map(::Type{ConstantTimescale}) = (; :sublimation_deposition_timescale => :τ_relax)

name_map(::Type{Homogeneous}) = (;
    :temperature_homogenous_nucleation => :T_hom,
    :homogeneous_freezing_timescale => :τ_hom,
)

name_map(::Type{Heterogeneous}) = (;
    :Reisner_et_al_A_parameter => :A,
    :Reisner_et_al_B_parameter => :B,
)

name_map(::Type{WithSupersaturation}) = (; :ice_snow_threshold_radius => :r_ice_snow)

name_map(::Type{CloudLiquidRainAccretion}) = (; :cloud_liquid_rain_collision_efficiency => :e)

name_map(::Type{CloudLiquidSnowAccretion}) = (; :cloud_liquid_snow_collision_efficiency => :e)

name_map(::Type{CloudIceRainAccretion}) = (; :cloud_ice_rain_collision_efficiency => :e)

name_map(::Type{CloudIceSnowAccretion}) = (; :cloud_ice_snow_collision_efficiency => :e)

name_map(::Type{RainSnowAccretion}) = (;
    :rain_snow_collision_efficiency => :e,
    :rain_snow_velocity_dispersion_coefficient => :coeff_disp,
)

# ═══════════════════════════════════════════════════════════════════
# 2-moment microphysics
# ═══════════════════════════════════════════════════════════════════

name_map(::Type{AcnvKK2000}) = (;
    :KK2000_autoconversion_coeff_A => :A,
    :KK2000_autoconversion_coeff_a => :a,
    :KK2000_autoconversion_coeff_b => :b,
    :KK2000_autoconversion_coeff_c => :c,
)
name_map(::Type{AccrKK2000}) = (;
    :KK2000_accretion_coeff_A => :A,
    :KK2000_accretion_coeff_a => :a,
    :KK2000_accretion_coeff_b => :b,
)
name_map(::Type{AcnvB1994}) = (;
    :B1994_autoconversion_coeff_C => :C,
    :B1994_autoconversion_coeff_a => :a,
    :B1994_autoconversion_coeff_b => :b,
    :B1994_autoconversion_coeff_c => :c,
    :B1994_autoconversion_coeff_N_0 => :N_0,
    :B1994_autoconversion_coeff_d_low => :d_low,
    :B1994_autoconversion_coeff_d_high => :d_high,
    :threshold_smooth_transition_steepness => :k,
)
name_map(::Type{AcnvTC1980}) = (;
    :TC1980_autoconversion_coeff_a => :a,
    :TC1980_autoconversion_coeff_b => :b,
    :TC1980_autoconversion_coeff_D => :D,
    :TC1980_autoconversion_coeff_r_0 => :r_0,
    :TC1980_autoconversion_coeff_me_liq => :me_liq,
    :threshold_smooth_transition_steepness => :k,
    :density_liquid_water => :m0_liq_coeff,
)
name_map(::Type{LD2004}) = (;
    :LD2004_R_6C_coeff => :R_6C_0,
    :LD2004_E_0_coeff => :E_0,
    :density_liquid_water => :ρ_w,
    :threshold_smooth_transition_steepness => :k,
)
name_map(::Type{RainParticlePDF_SB2006_limited}) = (;
    :SB2006_rain_distribution_coeff_nu => :νr,
    :SB2006_rain_distribution_coeff_mu => :μr,
    :SB2006_raindrops_min_mass => :xr_min,
    :SB2006_raindrops_max_mass => :xr_max,
    :SB2006_raindrops_size_distribution_coeff_N0_min => :N0_min,
    :SB2006_raindrops_size_distribution_coeff_N0_max => :N0_max,
    :SB2006_raindrops_size_distribution_coeff_lambda_min => :λ_min,
    :SB2006_raindrops_size_distribution_coeff_lambda_max => :λ_max,
    :density_liquid_water => :ρw,
    :SB2006_reference_air_density => :ρ0,
)
name_map(::Type{RainParticlePDF_SB2006_notlimited}) = (;
    :SB2006_rain_distribution_coeff_nu => :νr,
    :SB2006_rain_distribution_coeff_mu => :μr,
    :SB2006_raindrops_min_mass => :xr_min,
    :SB2006_raindrops_max_mass => :xr_max,
    :density_liquid_water => :ρw,
    :SB2006_reference_air_density => :ρ0,
)
name_map(::Type{CloudParticlePDF_SB2006}) = (;
    :SB2006_cloud_gamma_distribution_coeff_nu => :νc,
    :SB2006_cloud_gamma_distribution_coeff_mu => :μc,
    :SB2006_cloud_droplets_min_mass => :xc_min,
    :SB2006_raindrops_min_mass => :xc_max,
    :density_liquid_water => :ρw,
)
name_map(::Type{AcnvSB2006}) = (;
    :SB2006_collection_kernel_coeff_kcc => :kcc,
    :SB2006_raindrops_min_mass => :x_star,
    :SB2006_reference_air_density => :ρ0,
    :SB2006_autoconversion_correcting_function_coeff_A => :A,
    :SB2006_autoconversion_correcting_function_coeff_a => :a,
    :SB2006_autoconversion_correcting_function_coeff_b => :b,
)
name_map(::Type{AccrSB2006}) = (;
    :SB2006_collection_kernel_coeff_kcr => :kcr,
    :SB2006_accretion_correcting_function_coeff_tau0 => :τ0,
    :SB2006_reference_air_density => :ρ0,
    :SB2006_accretion_correcting_function_coeff_c => :c,
)
name_map(::Type{SelfColSB2006}) = (;
    :SB2006_collection_kernel_coeff_krr => :krr,
    :SB2006_collection_kernel_coeff_kapparr => :κrr,
    Symbol("SB2006_raindrops_self-collection_coeff_d") => :d,
)
name_map(::Type{BreakupSB2006}) = (;
    :SB2006_raindrops_equilibrium_mean_diameter => :Deq,
    :SB2006_raindrops_breakup_mean_diameter_threshold => :Dr_th,
    :SB2006_raindrops_breakup_coeff_kbr => :kbr,
    :SB2006_raindrops_breakup_coeff_kappabr => :κbr,
)
name_map(::Type{EvaporationSB2006}) = (;
    :SB2006_ventilation_factor_coeff_av => :av,
    :SB2006_ventilation_factor_coeff_bv => :bv,
    :SB2006_rain_evaporation_coeff_alpha => :α,
    :SB2006_rain_evaporation_coeff_beta => :β,
    :SB2006_reference_air_density => :ρ0,
)
name_map(::Type{NumberAdjustmentHorn2012}) = (;
    :Horn2012_number_concentration_adjustment_timescale => :τ,
)
name_map(::Type{CondEvap2M}) = (; :condensation_evaporation_timescale => :τ_relax)
name_map(::Type{SubDep2M}) = (; :sublimation_deposition_timescale => :τ_relax)

name_map(::Type{AccrB1994}) = (; :B1994_accretion_coeff_A => :A)

name_map(::Type{AccrTC1980}) = (; :TC1980_accretion_coeff_A => :A)

# ═══════════════════════════════════════════════════════════════════
# P3 scheme
# ═══════════════════════════════════════════════════════════════════

name_map(::Type{AreaPowerLaw}) = (; :M1996_area_coeff_gamma => :γ, :M1996_area_exponent_sigma => :σ)
name_map(::Type{SlopePowerLaw}) = (;
    :Heymsfield_mu_coeff1 => :a,
    :Heymsfield_mu_coeff2 => :b,
    :Heymsfield_mu_coeff3 => :c,
    :Heymsfield_mu_cutoff => :μ_max,
)
name_map(::Type{SlopeConstant}) = (; :P3_constant_slope_parameterization_value => :μ)
name_map(::Type{VentilationFactor}) = (;
    :SB2006_ventilation_factor_coeff_av => :aᵥ,
    :SB2006_ventilation_factor_coeff_bv => :bᵥ,
)
name_map(::Type{LocalRimeDensity}) = (;
    :CL1993_local_rime_density_constant_coeff => :a,
    :CL1993_local_rime_density_linear_coeff => :b,
    :CL1993_local_rime_density_quadratic_coeff => :c,
    :density_ice_water => :ρ_ice,
)

name_map(::Type{MassPowerLaw}) = (;
    :BF1995_mass_coeff_alpha => :α_va,
    :BF1995_mass_exponent_beta => :β_va,
)

name_map(::Type{ParametersP3}) = (;
    :density_ice_water => :ρ_i,  # TODO: Use `WaterProperties` struct for ice and liquid water density
    :density_liquid_water => :ρ_l,
    :temperature_water_freeze => :T_freeze,
    :P3_wet_growth_timescale => :τ_wet,
)

# ═══════════════════════════════════════════════════════════════════
# Terminal velocity
# ═══════════════════════════════════════════════════════════════════

name_map(::Type{StokesRegimeVelType}) = (;
    :density_liquid_water => :ρw,
    :kinematic_viscosity_of_air => :ν_air,
    :gravitational_acceleration => :grav,
)
name_map(::Type{SB2006VelType}) = (;
    :SB2006_reference_air_density => :ρ0,
    :SB2006_raindrops_terminal_velocity_coeff_aR => :aR,
    :SB2006_raindrops_terminal_velocity_coeff_bR => :bR,
    :SB2006_raindrops_terminal_velocity_coeff_cR => :cR,
    :density_liquid_water => :ρw,
    :kinematic_viscosity_of_air => :ν_air,
    :gravitational_acceleration => :grav,
)
# TODO: These should be array parameters.
name_map(::Type{Chen2022VelTypeSmallIce}) = (;
    :Chen2022_table_B3_As => :A,
    :Chen2022_table_B3_Bs => :B,
    :Chen2022_table_B3_Cs => :C,
    :Chen2022_table_B3_Es => :E,
    :Chen2022_table_B3_Fs => :F,
    :Chen2022_table_B3_Gs => :G,
    :Chen2022_ice_cutoff => :cutoff,
)
# TODO: These should be array parameters.
name_map(::Type{Chen2022VelTypeLargeIce}) = (;
    :Chen2022_table_B5_Al => :A,
    :Chen2022_table_B5_Bl => :B,
    :Chen2022_table_B5_Cl => :C,
    :Chen2022_table_B5_El => :E,
    :Chen2022_table_B5_Fl => :F,
    :Chen2022_table_B5_Gl => :G,
    :Chen2022_table_B5_Hl => :H,
    :Chen2022_ice_cutoff => :cutoff,
)
name_map(::Type{Chen2022VelTypeRain}) = (;
    :Chen2022_table_B1_q_coeff => :ρ0,
    :Chen2022_table_B1_ai => :a,
    :Chen2022_table_B1_a3_pow_coeff => :a3_pow,
    :Chen2022_table_B1_bi => :b,
    :Chen2022_table_B1_b_rho_coeff => :b_ρ,
    :Chen2022_table_B1_ci => :c,
)

name_map(::Type{Blk1MVelTypeRain}) = (;
    # TODO: rain's `r0` maps to `snow_flake_length_scale`, while the rain
    # mass and area power laws map to `rain_drop_length_scale`
    # (Microphysics1M.jl). Both are 1e-3 m in ClimaParams, so the two agree
    # numerically today; calibrating `rain_drop_length_scale` alone would
    # leave `get_v0` scaling rain fall speeds by the snow value. Confirm
    # which length scale is intended here and map to it.
    :snow_flake_length_scale => :r0,
    :rain_terminal_velocity_size_relation_coefficient_ve => :ve,
    :rain_terminal_velocity_size_relation_coefficient_delv => :Δv,
    :rain_terminal_velocity_size_relation_coefficient_chiv => :χv,
    :density_liquid_water => :ρw,
    :rain_drop_drag_coefficient => :C_drag,
    :gravitational_acceleration => :grav,
    # mass and area power-law exponents, used only in the gamma-function factors
    :rain_mass_size_relation_coefficient_me => :me,
    :rain_mass_size_relation_coefficient_delm => :Δm,
    :rain_cross_section_size_relation_coefficient_ae => :ae,
    :rain_cross_section_size_relation_coefficient_dela => :Δa,
)

name_map(::Type{Blk1MVelTypeSnow}) = (;
    :snow_flake_length_scale => :r0,
    :snow_terminal_velocity_size_relation_coefficient => :ve,
    :snow_terminal_velocity_size_relation_coefficient_delv => :Δv,
    :snow_terminal_velocity_size_relation_coefficient_chiv => :χv,
    # mass and area power-law exponents, used only in the gamma-function factors
    :snow_mass_size_relation_coefficient_me => :me,
    :snow_mass_size_relation_coefficient_delm => :Δm,
    :snow_cross_section_size_relation_coefficient => :ae,
    :snow_cross_section_size_relation_coefficient_dela => :Δa,
)

# ═══════════════════════════════════════════════════════════════════
# Aerosol activation and nucleation
# ═══════════════════════════════════════════════════════════════════

name_map(::Type{AerosolActivationParameters}) = (;
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

name_map(::Type{H2S04NucleationParameters}) = (;
    :mam3_nucleation_p_b_n_neutral => :p_b_n,
    :mam3_nucleation_p_b_i_ion_induced => :p_b_i,
    :mam3_nucleation_u_b_n_neutral => :u_b_n,
    :mam3_nucleation_u_b_i_ion_induced => :u_b_i,
    :mam3_nucleation_v_b_n_neutral => :v_b_n,
    :mam3_nucleation_v_b_i_ion_induced => :v_b_i,
    :mam3_nucleation_w_b_n_neutral => :w_b_n,
    :mam3_nucleation_w_b_i_ion_induced => :w_b_i,
    :mam3_nucleation_p_t_n_neutral => :p_t_n,
    :mam3_nucleation_p_t_i_ion_induced => :p_t_i,
    :mam3_nucleation_u_t_n_neutral => :u_t_n,
    :mam3_nucleation_u_t_i_ion_induced => :u_t_i,
    :mam3_nucleation_v_t_n_neutral => :v_t_n,
    :mam3_nucleation_v_t_i_ion_induced => :v_t_i,
    :mam3_nucleation_w_t_n_neutral => :w_t_n,
    :mam3_nucleation_w_t_i_ion_induced => :w_t_i,
    :mam3_nucleation_p_A_n_neutral => :p_A_n,
    :mam3_nucleation_p_A_i_ion_induced => :p_A_i,
    :mam3_nucleation_a_n_neutral => :a_n,
    :mam3_nucleation_a_i_ion_induced => :a_i,
)
name_map(::Type{OrganicNucleationParameters}) = (;
    :mam3_nucleation_a_1_neutral => :a_1,
    :mam3_nucleation_a_2_neutral => :a_2,
    :mam3_nucleation_a_3_ion_induced => :a_3,
    :mam3_nucleation_a_4_ion_induced => :a_4,
    :mam3_nucleation_a_5 => :a_5,
    :mam3_nucleation_Y_MTO3_percent => :Y_MTO3,
    :mam3_nucleation_Y_MTOH_percent => :Y_MTOH,
    :mam3_nucleation_k_MTO3_organic_factor => :k_MTO3,
    :mam3_nucleation_k_MTOH_organic_factor => :k_MTOH,
    :mam3_nucleation_exp_MTO3_organic_factor => :exp_MTO3,
    :mam3_nucleation_exp_MTOH_organic_factor => :exp_MTOH,
)
name_map(::Type{MixedNucleationParameters}) = (;
    :mam3_nucleation_k_H2SO4_mixed_organic_sulfuric_acid_factor =>
        :k_H2SO4org,
    :mam3_nucleation_k_MTOH_organic_factor => :k_MTOH,
    :mam3_nucleation_exp_MTOH_organic_factor => :exp_MTOH,
)

name_map(::Type{H2SO4SolutionParameters}) = (;
    :p_over_sulphuric_acid_solution_T_max => :T_max,
    :p_over_sulphuric_acid_solution_T_min => :T_min,
    :p_over_sulphuric_acid_solution_w_2 => :w_2,
    :p_over_sulphuric_acid_solution_c1 => :c1,
    :p_over_sulphuric_acid_solution_c2 => :c2,
    :p_over_sulphuric_acid_solution_c3 => :c3,
    :p_over_sulphuric_acid_solution_c4 => :c4,
    :p_over_sulphuric_acid_solution_c5 => :c5,
    :p_over_sulphuric_acid_solution_c6 => :c6,
    :p_over_sulphuric_acid_solution_c7 => :c7,
)

# ═══════════════════════════════════════════════════════════════════
# Ice nucleation
# ═══════════════════════════════════════════════════════════════════

name_map(::Type{Mohler2006}) = (;
    :Mohler2006_maximum_allowed_Si => :Sᵢ_max,
    :Mohler2006_threshold_T => :T_thr,
)
name_map(::Type{Koop2000}) = (;
    :Koop2000_min_delta_aw => :Δa_w_min,
    :Koop2000_max_delta_aw => :Δa_w_max,
    :Koop2000_J_hom_coeff1 => :c₁,
    :Koop2000_J_hom_coeff2 => :c₂,
    :Koop2000_J_hom_coeff3 => :c₃,
    :Koop2000_J_hom_coeff4 => :c₄,
    :Linear_J_hom_coeff1 => :linear_c₁,
    :Linear_J_hom_coeff2 => :linear_c₂,
)
name_map(::Type{MorrisonMilbrandt2014}) = (;
    :temperature_homogenous_nucleation => :T_dep_thres,
    :Thompson2004_c1_Cooper => :c₁,
    :Thompson2004_c2_Cooper => :c₂,
    :temperature_water_freeze => :T₀,
    :BarklieGokhale1959_a_parameter => :het_a,
    :BarklieGokhale1959_B_parameter => :het_B,
)
name_map(::Type{RainFreezing}) = (;
    :BarklieGokhale1959_a_parameter => :het_a,
    :BarklieGokhale1959_B_parameter => :het_B,
)
name_map(::Type{Frostenberg2023}) = (;
    :Frostenberg2023_standard_deviation => :σ,
    :Frostenberg2023_a_coefficient => :a,
    :Frostenberg2023_b_coefficient => :b,
    :temperature_water_freeze => :T_freeze,
)

# ═══════════════════════════════════════════════════════════════════
# Aerosol species
# ═══════════════════════════════════════════════════════════════════

name_map(::Type{ArizonaTestDust}) = (;
    :Mohler2006_S0_warm_ArizonaTestDust => :S₀_warm,
    :Mohler2006_S0_cold_ArizonaTestDust => :S₀_cold,
    :Mohler2006_a_warm_ArizonaTestDust => :a_warm,
    :Mohler2006_a_cold_ArizonaTestDust => :a_cold,
    :J_ABDINM_m_ArizonaTestDust => :deposition_m,
    :J_ABDINM_c_ArizonaTestDust => :deposition_c,
    :J_ABIFM_m_ArizonaTestDust => :ABIFM_m,
    :J_ABIFM_c_ArizonaTestDust => :ABIFM_c,
)

name_map(::Type{AsianDust}) = (;
    :J_ABDINM_m_AsianDust => :deposition_m,
    :J_ABDINM_c_AsianDust => :deposition_c,
    :J_ABIFM_m_AsianDust => :ABIFM_m,
    :J_ABIFM_c_AsianDust => :ABIFM_c,
)

name_map(::Type{DesertDust}) = (;
    :Mohler2006_S0_warm_DesertDust => :S₀_warm,
    :Mohler2006_S0_cold_DesertDust => :S₀_cold,
    :Mohler2006_a_warm_DesertDust => :a_warm,
    :Mohler2006_a_cold_DesertDust => :a_cold,
    :AlpertKnopf2016_J_ABIFM_m_DesertDust => :ABIFM_m,
    :AlpertKnopf2016_J_ABIFM_c_DesertDust => :ABIFM_c,
)

name_map(::Type{Dust}) = (;
    :J_ABDINM_m_Dust => :deposition_m,
    :J_ABDINM_c_Dust => :deposition_c,
    :J_ABIFM_m_Dust => :ABIFM_m,
    :J_ABIFM_c_Dust => :ABIFM_c,
)

name_map(::Type{Feldspar}) = (;
    :Alpert2022_J_deposition_m_Feldspar => :deposition_m,
    :Alpert2022_J_deposition_c_Feldspar => :deposition_c,
)

name_map(::Type{Ferrihydrite}) = (;
    :Alpert2022_J_deposition_m_Ferrihydrite => :deposition_m,
    :Alpert2022_J_deposition_c_Ferrihydrite => :deposition_c,
)

name_map(::Type{Illite}) = (;
    :J_ABDINM_m_Illite => :deposition_m,
    :J_ABDINM_c_Illite => :deposition_c,
    :KnopfAlpert2013_J_ABIFM_m_Illite => :ABIFM_m,
    :KnopfAlpert2013_J_ABIFM_c_Illite => :ABIFM_c,
)

name_map(::Type{Kaolinite}) = (;
    :China2017_J_deposition_m_Kaolinite => :deposition_m,
    :China2017_J_deposition_c_Kaolinite => :deposition_c,
    :KnopfAlpert2013_J_ABIFM_m_Kaolinite => :ABIFM_m,
    :KnopfAlpert2013_J_ABIFM_c_Kaolinite => :ABIFM_c,
)

name_map(::Type{MiddleEasternDust}) = (;
    :J_ABIFM_m_MiddleEasternDust => :ABIFM_m,
    :J_ABIFM_c_MiddleEasternDust => :ABIFM_c,
)

name_map(::Type{SaharanDust}) = (;
    :J_ABDINM_m_SaharanDust => :deposition_m,
    :J_ABDINM_c_SaharanDust => :deposition_c,
)

name_map(::Type{Seasalt}) = (;
    :seasalt_aerosol_molar_mass => :M,
    :seasalt_aerosol_density => :ρ,
    :seasalt_aerosol_osmotic_coefficient => :ϕ,
    :seasalt_aerosol_ion_number => :ν,
    :seasalt_aerosol_water_soluble_mass_fraction => :ϵ,
    :seasalt_aerosol_kappa => :κ,
)

name_map(::Type{Sulfate}) = (;
    :sulfate_aerosol_molar_mass => :M,
    :sulfate_aerosol_density => :ρ,
    :sulfate_aerosol_osmotic_coefficient => :ϕ,
    :sulfate_aerosol_ion_number => :ν,
    :sulfate_aerosol_water_soluble_mass_fraction => :ϵ,
    :sulfate_aerosol_kappa => :κ,
)

# ═══════════════════════════════════════════════════════════════════
# Generated constructors and summary
# ═══════════════════════════════════════════════════════════════════

# Structs whose `CP.ParamDict` constructor is exactly "read the name map, fill the fields by keyword".
const PLAIN_PARAMETER_TYPES = (
    ArizonaTestDust,
    AerosolActivationParameters,
    AsianDust,
    DesertDust,
    Dust,
    Feldspar,
    Ferrihydrite,
    Illite,
    Kaolinite,
    MiddleEasternDust,
    H2S04NucleationParameters,
    OrganicNucleationParameters,
    MixedNucleationParameters,
    SaharanDust,
    Seasalt,
    Sulfate,
    H2SO4SolutionParameters,
    AirProperties,
    Mohler2006,
    Koop2000,
    MorrisonMilbrandt2014,
    RainFreezing,
    Frostenberg2023,
    Parameters0M,
    CloudLiquid,
    IceNumberTemperatureFit,
    VarTimescaleAcnv,
    AcnvKK2000,
    AccrKK2000,
    AcnvB1994,
    LD2004,
    RainParticlePDF_SB2006_limited,
    RainParticlePDF_SB2006_notlimited,
    AcnvSB2006,
    AccrSB2006,
    SelfColSB2006,
    BreakupSB2006,
    NumberAdjustmentHorn2012,
    CondEvap2M,
    SubDep2M,
    SlopePowerLaw,
    SlopeConstant,
    VentilationFactor,
    LocalRimeDensity,
    StokesRegimeVelType,
    SB2006VelType,
    WaterProperties,
    AccrB1994,
    AccrTC1980,
)
for T in PLAIN_PARAMETER_TYPES
    @eval (::Type{$T})(td::CP.ParamDict) = $T(; make_params(td, name_map($T))...)
end

# 1-moment process options whose parameters are the NamedTuple of their name map.
const NAMEDTUPLE_OPTION_TYPES = (
    CloudLiquidFormation,
    ConstantTimescale,
    Homogeneous,
    Heterogeneous,
    WithSupersaturation,
    CloudLiquidRainAccretion,
    CloudLiquidSnowAccretion,
    CloudIceRainAccretion,
    CloudIceSnowAccretion,
    RainSnowAccretion,
)
for T in NAMEDTUPLE_OPTION_TYPES
    @eval process_params_for(::$T, td::CP.ParamDict) = make_params(td, name_map($T))
end

"""
    parameter_groups()

Ordered `title => keys` pairs listing every registered name map by documentation section. A key is a
parameter struct or option type, or a `(Constructor, Type)` tuple for the type-dispatched
constructors such as `(ParticleMass, Rain)`.
"""
parameter_groups() = (
    "Air and water properties" => (
        AirProperties,
        WaterProperties,
    ),
    "0-moment microphysics" => (
        Parameters0M,
    ),
    "1-moment microphysics" => (
        CloudLiquid,
        CloudIce,
        IceNumberTemperatureFit,
        (ParticleMass, CloudIce),
        Rain,
        (ParticleMass, Rain),
        (ParticleArea, Rain),
        Snow,
        (ParticleMass, Snow),
        (ParticleArea, Snow),
        VarTimescaleAcnv,
        KesslerAcnv,
    ),
    "1-moment process options" => (
        TemperatureDependent,
        NoSupersaturation,
        CloudLiquidFormation,
        ConstantTimescale,
        Homogeneous,
        Heterogeneous,
        WithSupersaturation,
        CloudLiquidRainAccretion,
        CloudLiquidSnowAccretion,
        CloudIceRainAccretion,
        CloudIceSnowAccretion,
        RainSnowAccretion,
    ),
    "2-moment microphysics" => (
        AcnvKK2000,
        AccrKK2000,
        AcnvB1994,
        AcnvTC1980,
        LD2004,
        RainParticlePDF_SB2006_limited,
        RainParticlePDF_SB2006_notlimited,
        CloudParticlePDF_SB2006,
        AcnvSB2006,
        AccrSB2006,
        SelfColSB2006,
        BreakupSB2006,
        EvaporationSB2006,
        NumberAdjustmentHorn2012,
        CondEvap2M,
        SubDep2M,
        AccrB1994,
        AccrTC1980,
    ),
    "P3 scheme" => (
        AreaPowerLaw,
        SlopePowerLaw,
        SlopeConstant,
        VentilationFactor,
        LocalRimeDensity,
        MassPowerLaw,
        ParametersP3,
    ),
    "Terminal velocity" => (
        StokesRegimeVelType,
        SB2006VelType,
        Chen2022VelTypeSmallIce,
        Chen2022VelTypeLargeIce,
        Chen2022VelTypeRain,
        Blk1MVelTypeRain,
        Blk1MVelTypeSnow,
    ),
    "Aerosol activation and nucleation" => (
        AerosolActivationParameters,
        H2S04NucleationParameters,
        OrganicNucleationParameters,
        MixedNucleationParameters,
        H2SO4SolutionParameters,
    ),
    "Ice nucleation" => (
        Mohler2006,
        Koop2000,
        MorrisonMilbrandt2014,
        RainFreezing,
        Frostenberg2023,
    ),
    "Aerosol species" => (
        ArizonaTestDust,
        AsianDust,
        DesertDust,
        Dust,
        Feldspar,
        Ferrihydrite,
        Illite,
        Kaolinite,
        MiddleEasternDust,
        SaharanDust,
        Seasalt,
        Sulfate,
    ),
)

name_map(key::Tuple) = name_map(key...)
name_map_label(T::Type) = Symbol(nameof(T))
name_map_label(key::Tuple) = Symbol(nameof(key[1]), "(", nameof(key[2]), ")")

"""
    name_map_summary()

`Dict(label => name_map(key))` over every registered name map (see [`parameter_groups`](@ref)); the
label is the struct or option name, or e.g. `Symbol("ParticleMass(Rain)")` for type-dispatched
constructors. Used by the documentation and by the parameter tests.
"""
name_map_summary() = Dict(name_map_label(k) => name_map(k) for (_, keys) in parameter_groups() for k in keys)
