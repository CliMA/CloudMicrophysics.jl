using Test
import CLIMAParameters as CP

include("../src/Nucleation.jl")
import .Nucleation


const FT = Float64
toml_dict = CP.create_toml_dict(FT)
param_names = (
    "u_b_n",
    "v_b_n",
    "w_b_n",
    "u_b_i",
    "v_b_i",
    "w_b_i",
    "u_t_n",
    "v_t_n",
    "w_t_n",
    "u_t_i",
    "v_t_i",
    "w_t_i",
    "p_t_n",
    "p_A_n",
    "a_n",
    "p_t_i",
    "p_A_i",
    "a_i",
    "p_b_n",
    "p_b_i",
    "a_1",
    "a_2",
    "a_3",
    "a_4",
    "a_5",
    "k_H2SO4org",
)
params = CP.get_parameter_values!(toml_dict, param_names)
params = (; params...)

@testset "Pure H2SO4 Nucleation Smoke Test" begin
    nh3_conc = 0.0
    negative_ion_conc = 0.0
    temp = 208
    test_h2so4(h2so4_c) =
        sum(
            Nucleation.h2so4_nucleation_rate(
                h2so4_c * 1e6,
                nh3_conc,
                negative_ion_conc,
                temp,
                params,
            ),
        ) * 1e-6

    h2so4_concentrations = 10 .^ (6:0.5:9)
    rates = [
        0.004530231195176104,
        0.4299078411040382,
        40.79720082267248,
        3871.5544026624943,
        367401.0272892414,
        3.4865457336815596e7,
        3.308646478955847e9,
    ]
    for (h2so4_conc, rate) in zip(h2so4_concentrations, rates)
        @test test_h2so4(h2so4_conc) ≈ rate rtol = 1e-5
    end
end


@testset "Pure Organic Nucleation Smoke Test" begin
    negative_ion_conc = 0.0
    test_organic(hom_c) =
        sum(
            Nucleation.organic_nucleation_rate(
                negative_ion_conc,
                hom_c,
                params,
            ),
        ) * 1e-6

    hom_concentrations = 10 .^ (6:0.5:8.5)
    rates = [
        7.778131746328608e-6,
        0.0024180450482058588,
        0.0400097,
        0.35954428147126755,
        2.944798579681883,
        24.176443961969532,
    ]
    for (hom_conc, rate) in zip(hom_concentrations, rates)
        @test test_organic(hom_conc) ≈ rate rtol = 1e-5
    end
end

@testset "Mixed Organic and H2SO4 Nucleation Smoke Test" begin
    h2so4_conc = 2.6e6
    test_mixed_nucleation(bioOxOrg) =
        Nucleation.organic_and_h2so4_nucleation_rate(
            h2so4_conc,
            bioOxOrg,
            params,
        ) * 1e6

    bioOxOrg_concentrations = 10 .^ (5.8:0.25:8.5)
    rates = [
        0.00697371914560178,
        0.012401221168017876,
        0.022052836262425032,
        0.03921610465843014,
        0.06973719145601781,
        0.12401221168017874,
        0.22052836262425027,
        0.39216104658430134,
        0.697371914560178,
        1.24012211680179,
        2.2052836262425073,
    ]
    for (bioOxOrg_conc, rate) in zip(bioOxOrg_concentrations, rates)
        @test test_mixed_nucleation(bioOxOrg_conc) ≈ rate rtol = 1e-5
    end
end
