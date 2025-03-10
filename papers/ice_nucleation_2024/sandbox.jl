import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import Thermodynamics as TD
import CairoMakie as MK

# To grab data
import DelimitedFiles
using LazyArtifacts
using ClimaUtilities.ClimaArtifacts
import NaNStatistics
import Interpolations as Intp

FT = Float64


function interpolated_data(raw_data, time_array)
    data_at_t = Intp.linear_interpolation(raw_data[:, 1], raw_data[:, 2])
    return data_at_t.(time_array)
end

function unpack_ACI04_07(FT)

    file_path = CM.ArtifactCalling.AIDA_ice_nucleation("ACI04_07_experiment_timeseries.edf")
    unpacked_data = DelimitedFiles.readdlm(file_path, ',',skipstart = 4362)

    N_ice1 = unpacked_data[:, 8]
    N_ice2 = unpacked_data[:, 10]
    for i in eachindex(N_ice1)
        if N_ice1[i] == -9999.99
            N_ice1[i] = FT(0)
        end
        if N_ice2[i] == -9999.990000000
            N_ice2[i] = FT(0)
        end
    end
    ACI04_07_t_profile = unpacked_data[:, 1]
    ACI04_07_T_profile = unpacked_data[:, 3]
    ACI04_07_P_profile = unpacked_data[:, 2] .* 1e2  # hPa to Pa
    ACI04_07_ICNC = (N_ice1 .+ N_ice2) .* 1e6      # Nᵢ [m^-3]
    ACI04_07_e = unpacked_data[:, 4] .* [TD.saturation_vapor_pressure(tps, T₀, TD.Ice()) for T₀ in ACI04_07_T_profile]
    
    data = (; ACI04_07_t_profile, ACI04_07_T_profile, ACI04_07_P_profile, ACI04_07_ICNC, ACI04_07_e)
    return data
end

function unpack_ACI04_22(FT)
    local data_file_name = "ACI04_22"
    local total_t = 400
    ACI04_22_t_profile = collect(0:1:total_t)  # TODO - start from 0 or 1?

    ACI04_22_T_profile = zeros(total_t)
    ACI04_22_P_profile = zeros(total_t)
    ACI04_22_ICNC = zeros(total_t)
    ACI04_22_e = zeros(total_t)

    ACI04_22_T_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_T.csv")
    ACI04_22_P_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_P.csv")
    ACI04_22_ICNC_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_N_ice.csv")
    ACI04_22_e_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_RH_water.csv")
    
    ACI04_22_raw_T_data = DelimitedFiles.readdlm(ACI04_22_T_file_path, ',' , header = false)
    ACI04_22_raw_P_data = DelimitedFiles.readdlm(ACI04_22_P_file_path, ',' , header = false)
    ACI04_22_raw_ICNC_data = zeros(11, 2)
    for i in eachindex(ACI04_22_raw_ICNC_data[:,1])
        ACI04_22_raw_ICNC_data[i,1] =  ACI04_22_raw_ICNC_data[i,1] .+ ((i-1) * 7)
    end
    ACI04_22_raw_ICNC_data = vcat(ACI04_22_raw_ICNC_data, DelimitedFiles.readdlm(ACI04_22_ICNC_file_path, ',' , header = false))
    ACI04_22_raw_e_data = DelimitedFiles.readdlm(ACI04_22_e_file_path, ',' , header = false)

    ACI04_22_T_profile = interpolated_data(ACI04_22_raw_T_data, ACI04_22_t_profile)
    ACI04_22_P_profile = interpolated_data(ACI04_22_raw_P_data, ACI04_22_t_profile) .* 100
    ACI04_22_ICNC = interpolated_data(ACI04_22_raw_ICNC_data, ACI04_22_t_profile) .* 1e6
    ACI04_22_e = interpolated_data(ACI04_22_raw_e_data, ACI04_22_t_profile)
    
    data = (; ACI04_22_t_profile, ACI04_22_T_profile, ACI04_22_P_profile, ACI04_22_ICNC, ACI04_22_e)
    return data
end

function unpack_EXP19(FT)
    local data_file_name = "EXP19"
    local total_t = 240
    EXP19_t_profile = collect(0:1:total_t)  # TODO - start from 0 or 1?

    EXP19_T_profile = zeros(total_t)
    EXP19_P_profile = zeros(total_t)
    EXP19_ICNC = zeros(total_t)
    EXP19_e = zeros(total_t)

    EXP19_T_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_T.csv")
    EXP19_P_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_P.csv")
    EXP19_ICNC_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_N_ice.csv")
    EXP19_e_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_RH_water.csv")
    
    EXP19_raw_T_data = vcat(DelimitedFiles.readdlm(EXP19_T_file_path, ',' , header = false), [205.7971 243.45; 228.9855 244.353; 246.37681 245.0297])
    EXP19_raw_P_data = DelimitedFiles.readdlm(EXP19_P_file_path, ',' , header = false)
    EXP19_raw_ICNC_data = DelimitedFiles.readdlm(EXP19_ICNC_file_path, ',' , header = false)
    EXP19_raw_e_data = vcat(DelimitedFiles.readdlm(EXP19_e_file_path, ',' , header = false), [242.8571 0.69655])

    EXP19_T_profile = interpolated_data(EXP19_raw_T_data, EXP19_t_profile)
    EXP19_P_profile = interpolated_data(EXP19_raw_P_data, EXP19_t_profile) .* 100
    EXP19_ICNC = interpolated_data(EXP19_raw_ICNC_data, EXP19_t_profile) .* 1e6
    EXP19_e = interpolated_data(EXP19_raw_e_data, EXP19_t_profile)
    
    data = (; EXP19_t_profile, EXP19_T_profile, EXP19_P_profile, EXP19_ICNC, EXP19_e)
    return data
end

function unpack_EXP45(FT)
    local data_file_name = "EXP45"
    local total_t = 200
    EXP45_t_profile = collect(0:1:total_t)  # TODO - start from 0 or 1?

    EXP45_T_profile = zeros(total_t)
    EXP45_P_profile = zeros(total_t)
    EXP45_ICNC = zeros(total_t)
    EXP45_e = zeros(total_t)

    EXP45_T_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_T.csv")
    EXP45_P_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_P.csv")
    EXP45_ICNC_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_N_ice.csv")
    EXP45_e_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_RH_water.csv")
    
    EXP45_raw_T_data = vcat(DelimitedFiles.readdlm(EXP45_T_file_path, ',' , header = false), [242.65 208.63])
    EXP45_raw_P_data = DelimitedFiles.readdlm(EXP45_P_file_path, ',' , header = false)
    EXP45_raw_ICNC_data = vcat(DelimitedFiles.readdlm(EXP45_ICNC_file_path, ',' , header = false), [184.59 9.1029; 212.87 8.7541])
    EXP45_raw_e_data = DelimitedFiles.readdlm(EXP45_e_file_path, ',' , header = false)

    EXP45_T_profile = interpolated_data(EXP45_raw_T_data, EXP45_t_profile)
    EXP45_P_profile = interpolated_data(EXP45_raw_P_data, EXP45_t_profile) .* 100
    EXP45_ICNC = interpolated_data(EXP45_raw_ICNC_data, EXP45_t_profile) .* 1e6
    EXP45_e = interpolated_data(EXP45_raw_e_data, EXP45_t_profile)
    
    data = (; EXP45_t_profile, EXP45_T_profile, EXP45_P_profile, EXP45_ICNC, EXP45_e)
    return data
end

function unpack_EXP28(FT)
    local data_file_name = "EXP28"
    local total_t = 300
    EXP28_t_profile = collect(0:1:total_t)  # TODO - start from 0 or 1?

    EXP28_T_profile = zeros(total_t)
    EXP28_P_profile = zeros(total_t)
    EXP28_ICNC = zeros(total_t)
    EXP28_e = zeros(total_t)

    EXP28_T_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_T.csv")
    EXP28_P_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_P.csv")
    EXP28_ICNC_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_N_ice.csv")
    EXP28_e_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_RH_water.csv")
    
    EXP28_raw_T_data = DelimitedFiles.readdlm(EXP28_T_file_path, ',' , header = false)
    EXP28_raw_P_data = DelimitedFiles.readdlm(EXP28_P_file_path, ',' , header = false)
    EXP28_raw_ICNC_data = DelimitedFiles.readdlm(EXP28_ICNC_file_path, ',' , header = false)
    EXP28_raw_e_data = DelimitedFiles.readdlm(EXP28_e_file_path, ',' , header = false)

    EXP28_T_profile = interpolated_data(EXP28_raw_T_data, EXP28_t_profile)
    EXP28_P_profile = interpolated_data(EXP28_raw_P_data, EXP28_t_profile) .* 100
    EXP28_ICNC = interpolated_data(EXP28_raw_ICNC_data, EXP28_t_profile) .* 1e6
    EXP28_e = interpolated_data(EXP28_raw_e_data, EXP28_t_profile)
    
    data = (; EXP28_t_profile, EXP28_T_profile, EXP28_P_profile, EXP28_ICNC, EXP28_e)
    return data
end

ACI04_07 = unpack_ACI04_07(FT)
ACI04_22 = unpack_ACI04_22(FT)
EXP19 = unpack_EXP19(FT)
EXP45 = unpack_EXP45(FT)
EXP28 = unpack_EXP28(FT)

(; ACI04_07_t_profile, ACI04_07_T_profile, ACI04_07_P_profile, ACI04_07_ICNC, ACI04_07_e) = ACI04_07
(; ACI04_22_t_profile, ACI04_22_T_profile, ACI04_22_P_profile, ACI04_22_ICNC, ACI04_22_e) = ACI04_22
(; EXP19_t_profile, EXP19_T_profile, EXP19_P_profile, EXP19_ICNC, EXP19_e) = EXP19
(; EXP45_t_profile, EXP45_T_profile, EXP45_P_profile, EXP45_ICNC, EXP45_e) = EXP45
(; EXP28_t_profile, EXP28_T_profile, EXP28_P_profile, EXP28_ICNC, EXP28_e) = EXP28

ACI04_07_fig = MK.Figure(size = (1000, 600), fontsize = 24)
ACI04_07_ax = MK.Axis(ACI04_07_fig[1, 1], ylabel = "ICNC [m^-3]", xlabel = "time [s]", title = "ACI04_07 data")
MK.lines!(ACI04_07_ax, ACI04_07_t_profile, ACI04_07_ICNC, label = "Raw AIDA", color =:blue, linestyle =:dash, linewidth = 2)
# ACI04_07_fig

ACI04_22_fig = MK.Figure(size = (1000, 600), fontsize = 24)
ACI04_22_ax = MK.Axis(ACI04_22_fig[1, 1], ylabel = "ICNC [m^-3]", xlabel = "time [s]", title = "ACI04_22 data")
MK.lines!(ACI04_22_ax, ACI04_22_t_profile, ACI04_22_ICNC, label = "Raw AIDA", color =:blue, linestyle =:dash, linewidth = 2)
# ACI04_22_fig

EXP19_fig = MK.Figure(size = (1000, 600), fontsize = 24)
EXP19_ax = MK.Axis(EXP19_fig[1, 1], ylabel = "ICNC [m^-3]", xlabel = "time [s]", title = "EXP19 data")
MK.lines!(EXP19_ax, EXP19_t_profile, EXP19_ICNC, label = "Raw AIDA", color =:blue, linestyle =:dash, linewidth = 2)
# EXP19_fig

EXP45_fig = MK.Figure(size = (1000, 600), fontsize = 24)
EXP45_ax = MK.Axis(EXP45_fig[1, 1], ylabel = "ICNC [m^-3]", xlabel = "time [s]", title = "EXP45 data")
MK.lines!(EXP45_ax, EXP45_t_profile, EXP45_ICNC, label = "Raw AIDA", color =:blue, linestyle =:dash, linewidth = 2)
# EXP45_fig

EXP28_fig = MK.Figure(size = (1000, 600), fontsize = 24)
EXP28_ax = MK.Axis(EXP28_fig[1, 1], ylabel = "ICNC [m^-3]", xlabel = "time [s]", title = "EXP28 data")
MK.lines!(EXP28_ax, EXP28_t_profile, EXP28_ICNC, label = "Raw AIDA", color =:blue, linestyle =:dash, linewidth = 2)
EXP28_fig
