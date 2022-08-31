# load simulated data, and fit the peaks, estimate the alignment factor, remove the noise and estimate the sensitivity matrices
## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)
rc("figure",max_open_warning=50)
using myPlot
color_array = ["tab:blue"; "tab:orange"; "tab:green"; "tab:red"; "tab:purple"; "tab:brown"; "tab:pink"; "tab:gray"; "tab:olive"; "tab:cyan"; "magenta"; "yellow"; "hotpink"; "darkmagenta"; "chartreuse"; "deepskyblue"; "navy"; "darkcyan"; "crimson"; "firebrick"]; 

# data manipulation (loading, writing, etc)
using Printf
# using XLSX # CSV does not deal with multiple sheets
# using DataValues
using XPSfile
using DataFrames
using Query

# scientific package from the official Julia repositories
using LinearAlgebra
using StatsBase
using Interpolations

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)
# using XPSsampling

PLOT_FIG  = true
SAVE_FIG  = true
SAVE_DATA = true
USING_GT_PROFILE = false

data_folder = "../data/cylinder_radius_10.0/peak_shift/eal_5_restricted_range/"

##
## load data and meta data
##

folder_content = readdir(data_folder);
list_match = match.(r"^offcenter_",folder_content)
data_folders = folder_content[list_match.!=nothing]
α_gt       = zeros(Cdouble,5*length(data_folders));
α_ratio    = zeros(Cdouble,5*length(data_folders));
xc_off     = zeros(Cdouble,5*length(data_folders));
α_al_noise = zeros(Cdouble,5*length(data_folders));
# τ_al_noise_gt = zeros(Cdouble,5*length(data_folders));
# unit conversion constant (some of the quantities are in μm, some in L and some in Mbarn)
NA = 6.022e23;
κ_simple_units = 1.0e-37*NA; # simplified model
κ_units        = 1.0e-25*NA; # original model

idx = 1
for folder_tag in data_folders
    # local dictAllData,df,symbolDict = dataAndMeta_xlsx2df(string(data_folder,folder_tag,"water/data.xlsx"),r"^hν_[0-9]*",r"meta");
    # local dictAllGeom,symbolDictInv = model_xlsx2df(string(data_folder,folder_tag,"water/model.xlsx"));
    # local dictAllData,df,symbolDict = dataAndMeta_xlsx2df(string(data_folder,folder_tag,"/new_bg/water/data.xlsx"),r"^hν_[0-9]*",r"meta");
    # local dictAllGeom,symbolDictInv = model_xlsx2df(string(data_folder,folder_tag,"/new_bg/water/model.xlsx"));
    local dictAllData,df,symbolDict = dataAndMeta_xlsx2df(string(data_folder,folder_tag,"/new_eal/water/data.xlsx"),r"^hν_[0-9]*",r"meta");
    local dictAllGeom,symbolDictInv = model_xlsx2df(string(data_folder,folder_tag,"/new_eal/water/model.xlsx"));
    local Ndata                     = length(dictAllData);

    global idx
    for sym_data in keys(symbolDict)
        α_gt[idx]    = dictAllGeom[symbolDict[sym_data]].α[1]
        α_ratio[idx] = dictAllData[sym_data].α[1]
        xc_off[idx]  = dictAllGeom[symbolDict[sym_data]].xc[1]
        local H_liq = dictAllGeom[symbolDict[sym_data]][!,:H]
        local Nr = length(H_liq);
        local S_liq_noisy = Array{Cdouble,1}(dictAllData[sym_data][!,:Snoisy]) - Array{Cdouble,1}(dictAllData[sym_data][!,:SpectrumA_1_gas]);
        local σ_cs_liq = dictAllData[sym_data][!,:σ_cs_dens]
        local Sbg = dictAllData[sym_data][!,:Sbg]
        if USING_GT_PROFILE
            local ρ = dictAllGeom[symbolDict[sym_data]][!,:ρ] # it works way too well with the true concentration profile
        else
            local ρ = 55.49*ones(Cdouble,Nr)
        end
        α_al_noise[idx],_    = noiseAndParameterEstimation(σ_cs_liq,H_liq,S_liq_noisy,Sbg,ρ)
        local Tj = dictAllData[sym_data][!,:T][1];
        local Fνj = dictAllData[sym_data][!,:F][1];
        local σ_tot = dictAllData[sym_data][!,:σ_tot][1];
        local Δt = dictAllData[sym_data][!,:Δt][1];
        α_al_noise[idx] = α_al_noise[idx]/(κ_units*Tj*Fνj*σ_tot*Δt)
        # α_al_noise[idx] = 1.0e12α_al_noise[idx] # convert to m^{-2} from μm^{-2}
        # τ_al_noise_gt[idx],_ = noiseAndParameterEstimation(dictAllData[symbol_h][1][:σ_cs_dens],dictAllGeom[simbol_λ][1][:H],
        #                     Array{Cdouble,1}(dictAllData[symbol_h][1][:Snoisy]),dictAllData[symbol_h][1][:Sbg],ρA_1)

        idx = idx + 1
    end

    # plot
    if PLOT_FIG
        local ph_sym_list = keys(symbolDict)
        figure(figsize=[12, 10])
        for (i,plot_sym) in zip(1:4,ph_sym_list)
            local Be = dictAllData[plot_sym].Be;
            local Sliq = dictAllData[plot_sym].SpectrumA_1;
            local Sgas = dictAllData[plot_sym].SpectrumA_1_gas;
            local Sbg = dictAllData[plot_sym].Sbg;
            local Snoisy = dictAllData[plot_sym].Snoisy;
            local Eph = dictAllData[plot_sym].hν[1];
            
            local ax = subplot(2,2,i)
            
            scatter(Be,Snoisy,label="Data")
            plot(Be,Sliq,label="Liquid signal")
            plot(Be,Sgas,label="Gas signal")
            plot(Be,Sbg,label="Background signal")
            ylim(0.0)
            xlabel("binding energy [eV]",fontsize=14); 
            ylabel("spectrum [count]",fontsize=14) 
            xticks(fontsize=14); yticks(fontsize=14); 
            ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
            ax.yaxis.offsetText.set_size(14)
            ax.xaxis.offsetText.set_size(14)
            legend(fontsize=14)
            ax.invert_xaxis()
            ax.text(0.1, 0.5, string("Eph = ",Eph," [eV]"), transform=ax.transAxes,fontsize=14)
        end
        tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
    end
end


# get the mean values for each estimated alignment so that it 
Ndata = 5;
α_gt_mean    = zeros(Cdouble,length(data_folders));
α_ratio_mean = zeros(Cdouble,length(data_folders));
α_al_mean    = zeros(Cdouble,length(data_folders));
for i in 1: length(data_folders)
    α_gt_mean[i]    = mean(α_gt[(i-1)*Ndata+1:i*Ndata])
    α_ratio_mean[i] = mean(α_ratio[(i-1)*Ndata+1:i*Ndata])
    α_al_mean[i]    = mean(α_al_noise[(i-1)*Ndata+1:i*Ndata])
end


if PLOT_FIG
    figure(figsize=[10, 5]); # scatter(1.0e9α_gt,α_ratio)
    ax = subplot(122)
    scatter(1.0e8α_gt_mean,1.0e-2α_ratio_mean,color="tab:orange",label="liquid/vapor ratio") # α_ratio_mean.^2.5
    scatter(1.0e8α_gt_mean,1.0e8α_al_mean,color="tab:green",label="noise estimation [cm\$^{-2}\$]") 
    xlim(-0.001,0.033); ylim(-0.001)
    xlabel("GT [cm\$^{-2}\$]",fontsize=14); 
    ylabel("estimation",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax.yaxis.offsetText.set_size(14)
    ax.xaxis.offsetText.set_size(14)
    legend(fontsize=14)

    ax = subplot(121)
    scatter(xc_off[1:5:end],1.0e8α_gt_mean,label="GT [cm\$^{-2}\$]")
    scatter(xc_off[1:5:end],1.0e-2α_ratio_mean,label="liquid/vapor area ratio")
    scatter(xc_off[1:5:end],1.0e8α_al_mean,label="noise estimation [cm\$^{-2}\$]")
    ylim(-0.001)
    xlabel("horizontal offset [\$\\mu\$m]",fontsize=14); 
    ylabel("alignment parameter",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax.yaxis.offsetText.set_size(14)
    ax.xaxis.offsetText.set_size(14)
    legend(fontsize=14)
    tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)

    if SAVE_FIG
        if USING_GT_PROFILE
            savefig(string(data_folder,"liquid_vapor_area_ratio_and_noise_estimation_vs_alignment_parameter_units_gt.png"))
            savefig(string(data_folder,"liquid_vapor_area_ratio_and_noise_estimation_vs_alignment_parameter_units_gt.pdf"))
        else
            savefig(string(data_folder,"liquid_vapor_area_ratio_and_noise_estimation_vs_alignment_parameter_units.png"))
            savefig(string(data_folder,"liquid_vapor_area_ratio_and_noise_estimation_vs_alignment_parameter_units.pdf"))
        end
    end
end


# figure(); scatter(collect(1:length(α_al_noise)),1.0e9α_al_noise)

df_α =  DataFrame( "α_gt"=> α_gt_mean, "α_ratio" =>  α_ratio_mean, "x_offset" => xc_off[1:5:end]);
if SAVE_DATA
    # XLSX.writetable(string(data_folder,"/alignment_list.xlsx"), collect(eachcol(df_α)), names(df_α));
    if USING_GT_PROFILE
        XLSX.writetable(string(data_folder,"/alignment_list_new_eal_gt.xlsx"), collect(eachcol(df_α)), names(df_α));
    else
        XLSX.writetable(string(data_folder,"/alignment_list_new_eal.xlsx"), collect(eachcol(df_α)), names(df_α));
    end
end
