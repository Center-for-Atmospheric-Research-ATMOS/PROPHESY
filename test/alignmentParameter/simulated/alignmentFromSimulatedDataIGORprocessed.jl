# load simulated data, and fit the peaks, estimate the alignment factor, remove the noise and estimate the sensitivity matrices
## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)
rc("figure",max_open_warning=50)
# using myPlot
color_array = ["tab:blue"; "tab:orange"; "tab:green"; "tab:red"; "tab:purple"; "tab:brown"; "tab:pink"; "tab:gray"; "tab:olive"; "tab:cyan"; "magenta"; "yellow"; "hotpink"; "darkmagenta"; "chartreuse"; "deepskyblue"; "navy"; "darkcyan"; "crimson"; "firebrick"]; 

# data manipulation (loading, writing, etc)
# using Printf
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
# using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)
# using XPSsampling

PLOT_FIG  = true
SAVE_FIG  = !true
SAVE_DATA = !true
USING_GT_PROFILE = false

data_folder = "../../../data/cylinder_radius_10.0/peak_shift/eal_5_restricted_range/"

##
## load data and meta data
##

folder_content = readdir(data_folder);
# list_match = match.(r"^offcenter_",folder_content)
# data_folders = folder_content[list_match.!=nothing]
# α_gt       = zeros(Cdouble,5*length(data_folders));
# α_ratio    = zeros(Cdouble,5*length(data_folders));
# xc_off     = zeros(Cdouble,5*length(data_folders));
# α_al_noise = zeros(Cdouble,5*length(data_folders));
# τ_al_noise_gt = zeros(Cdouble,5*length(data_folders));
# unit conversion constant (some of the quantities are in μm, some in L and some in Mbarn)
NA = 6.022e23;
κ_simple_units = 1.0e-37*NA; # simplified model
κ_units        = 1.0e-25*NA; # original model


filenames =["IGOR/data_c1s_newEAL_offcenter_0_0001.xlsx";
"IGOR/data_o1s_newEAL_offcenter_0.xlsx";
"IGOR/data_c1s_newEAL_offcenter_1_0001.xlsx";
"IGOR/data_o1s_newEAL_offcenter_1.xlsx";
"IGOR/data_c1s_newEAL_offcenter_2_0001.xlsx";
"IGOR/data_o1s_newEAL_offcenter_2.xlsx"];

filenamesModels = ["/offcenter_0/new_eal/0001/model.xlsx";
                    "/offcenter_0/new_eal/water/model.xlsx";
                    "/offcenter_1/new_eal/0001/model.xlsx";
                    "/offcenter_1/new_eal/water/model.xlsx";
                    "/offcenter_2/new_eal/0001/model.xlsx";
                    "/offcenter_2/new_eal/water/model.xlsx"];

idx_file = 6;
FLAG_O1S_VECT = [false; true; false; true; false; true];
FLAG_O1S = FLAG_O1S_VECT[idx_file];
regIGORsheet  = r"^hν_[0-9]*_IGOR"
regDataColumn = r"((\bBe\b)|(\bCurve_[0-9]*\b)|(\bAll_curves\b)|(\bBackground\b))"
regMetaColumn = r"((\b[Bb]inding_energy\b)|(\b[Ii]ntensity\b)|(\b[Aa]rea\b)|(\b[Ll]orentzian\b)|(\b[Gg]aussian\b)|(\b[Ff]ull_width\b))";

# # load data and model
# dictAllData,_,symbolDict = dataAndMeta_xlsx2df(string(data_folder,filenames[idx_file]),r"^hν_[0-9]*$",r"meta");
# dictAllGeom,_ = model_xlsx2df(string(data_folder,filenamesModels[idx_file]));

# # load IGOR results
# dataIGOR,metaIGOR = IGORcolumn_xlsx2df(string(data_folder,filenames[idx_file]),regIGORsheet,regDataColumn,regMetaColumn);


# automation of the APE
function APEsimu(dictAllData::Dict,dictAllGeom::Dict,dataIGOR::Dict,symbolDict::Dict,FLAG_O1S::Bool)
    dictAllAPE  = Dict()
    dictAllPlot = Dict()
    dictAllIGOR = Dict()
    for plot_sym in keys(symbolDict)
        local Beplot = dictAllData[plot_sym].Be;
        local dBe = abs(Beplot[2]-Beplot[1]);
        idx_1 = findfirst(abs.(Beplot.-dataIGOR[plot_sym].Be[1]).<0.001) 
        idx_2 = findfirst(abs.(Beplot.-dataIGOR[plot_sym].Be[end]).<0.001)
        local H_liq = dictAllGeom[symbolDict[plot_sym]][!,:H];
        if FLAG_O1S
            local S_noisy_GT = Array{Cdouble,1}(dictAllData[plot_sym][!,:Snoisy]) - Array{Cdouble,1}(dictAllData[plot_sym][!,:SpectrumA_1_gas]);
            local S_noisy = dictAllData[plot_sym].Snoisy[idx_1:idx_2] - dataIGOR[plot_sym].Curve_2 # or 2 # # 
            local σ_cs = dataIGOR[plot_sym].Curve_1/(dBe*sum(dataIGOR[plot_sym].Curve_1)) # dictAllData[plot_sym][!,:σ_cs_dens][idx_1:idx_2]
            local Spectrum_OI_GT = Array{Cdouble,1}(dictAllData[plot_sym][!,:SpectrumA_1])
            local Spectrum_OI_IGOR = dataIGOR[plot_sym].Curve_1
        else
            local S_noisy_GT = Array{Cdouble,1}(dictAllData[plot_sym][!,:Snoisy])
            local S_noisy = dictAllData[plot_sym].Snoisy[idx_1:idx_2]
            local σ_cs = (dataIGOR[plot_sym].Curve_1+dataIGOR[plot_sym].Curve_2+dataIGOR[plot_sym].Curve_3)/(dBe*sum(dataIGOR[plot_sym].Curve_1+dataIGOR[plot_sym].Curve_2+dataIGOR[plot_sym].Curve_3))
            local Spectrum_OI_GT = Array{Cdouble,1}(dictAllData[plot_sym][!,:SpectrumA_1])
            local Spectrum_OI_IGOR = dataIGOR[plot_sym].Curve_1+dataIGOR[plot_sym].Curve_2+dataIGOR[plot_sym].Curve_3
        end
        local Sbg = dataIGOR[plot_sym].Background;
        local ρ_gt = dictAllGeom[symbolDict[plot_sym]][!,:ρ] # it works way too well with the true concentration profile
        local ρ = ρ_gt[1]*ones(Cdouble,length(H_liq));
        global α_al_noise,_               = noiseAndParameterEstimation(σ_cs,H_liq,S_noisy,Sbg,ρ_gt)
        global α_al_approx,S_noise_approx = noiseAndParameterEstimation(σ_cs,H_liq,S_noisy,Sbg,ρ)
        local Fνj = dictAllData[plot_sym].F[1]
        local Δt = dictAllData[plot_sym].Δt[1]
        local σ_tot = dictAllData[plot_sym].σ_tot[1];
        local Tj = dictAllData[plot_sym].T[1];
        α_al_noise  = α_al_noise/(κ_units*Tj*Fνj*σ_tot*Δt)
        α_al_approx = α_al_approx/(κ_units*Tj*Fνj*σ_tot*Δt)
        global α_GT = dictAllGeom[symbolDict[plot_sym]].α[1]
        # create dataframe 
        dictData = Dict( "α_GT" => α_GT, "α_APE" => α_al_noise, "α_APE_approx" => α_al_approx, "name" => plot_sym, "other_name" => symbolDict[plot_sym], 
                            "hν" => dictAllData[plot_sym].hν[1], "λ" => dictAllGeom[symbolDict[plot_sym]].λ[1], "xc" => dictAllGeom[symbolDict[plot_sym]].xc[1])
        dictPlotData = Dict( "Be" => Beplot,                "S_noisy" => dictAllData[plot_sym].Snoisy,              "SOI_noisy" => S_noisy_GT, 
                            "Spectrum_OI" => Spectrum_OI_GT  , "bg" => dictAllData[plot_sym].Sbg, "noise estimation" => S_noise_approx)
        dictPlotIGOR = Dict( "Be" => dataIGOR[plot_sym].Be, "S_noisy" => dictAllData[plot_sym].Snoisy[idx_1:idx_2], "SOI_noisy" => S_noisy,    
                            "Spectrum_OI" => Spectrum_OI_IGOR, "bg" => dataIGOR[plot_sym].Background) 
        # push it the global 
        dictAllAPE[plot_sym] = dictData;
        dictAllPlot[plot_sym] = dictPlotData;
        dictAllIGOR[plot_sym] = dictPlotIGOR;
    end

    # create the data frame out of the disctionary (later, from each of these DataFrames, create a dictionary whose entries are the file names)
    # return
    dictAllAPE,dictAllPlot,dictAllIGOR
end

dictAllAPE = Dict();
dictAllPlot = Dict();
dictAllIGOR = Dict();
for (dataFileName,modelFileName,FLAG_O) in zip(filenames,filenamesModels,FLAG_O1S_VECT)
    # load data and model
    local dictAllData,_,symbolDict = dataAndMeta_xlsx2df(string(data_folder,dataFileName),r"^hν_[0-9]*$",r"meta");
    local dictAllGeom,_ = model_xlsx2df(string(data_folder,modelFileName));

    # load IGOR results
    local dataIGOR,metaIGOR = IGORcolumn_xlsx2df(string(data_folder,dataFileName),regIGORsheet,regDataColumn,regMetaColumn);

    # get the data, fits, background and alignments
    local dictAPE,dictPlot,dictIGOR = APEsimu(dictAllData,dictAllGeom,dataIGOR,symbolDict,FLAG_O)

    # push in global dictionary
    dictAllAPE[dataFileName] = dictAPE;
    dictAllPlot[dataFileName] = dictPlot;
    dictAllIGOR[dataFileName] = dictIGOR;
end


# plot the alignments
fileNamesC1sAll = filenames[match.(r"[cC]1[sS]",filenames).!=nothing]
fileNamesO1sAll = filenames[match.(r"[oO]1[sS]",filenames).!=nothing]
TWO_COLUMN = true
if TWO_COLUMN
    FONTSIZE = 16
else
    FONTSIZE = 14
end
if TWO_COLUMN
    LEGEND_FONTSIZE = 14
else
    LEGEND_FONTSIZE = 12
end
LINEWIDTH = 2.5
if TWO_COLUMN
    figure(figsize=[14, 5])
else
    figure(figsize=[12, 5])
end
ax1 = subplot(121)
ax1.plot(1.0e8*[1e-11; 3.0e-8],1.0e8*[1e-11; 3.0e-8],label="1:1",linewidth=LINEWIDTH)
xlim(1.0e8*1e-11,1.0e8*3.0e-8)
ylim(1.0e8*1e-11,1.0e8*3.0e-8)
yscale("log")
xscale("log")
xlabel("alignment parameter GT [cm\$^{-2}\$]",fontsize=FONTSIZE); 
ylabel("model estimation [cm\$^{-2}\$]",fontsize=FONTSIZE) 
xticks(fontsize=FONTSIZE); yticks(fontsize=FONTSIZE); 
# ax1.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax1.yaxis.offsetText.set_size(FONTSIZE)
ax1.xaxis.offsetText.set_size(FONTSIZE)
ax2 = subplot(122)
ax2.plot(1.0e8*[1e-11; 3.0e-8],1.0e8*[1e-11; 3.0e-8],label="1:1",linewidth=LINEWIDTH)
ax2.plot(1.0e8*[1e-11; 3.0e-8],1.0e8*0.5*[1e-11; 3.0e-8],label="2:1",linewidth=LINEWIDTH)
xlim(1.0e8*1e-11,1.0e8*3.0e-8)
ylim(1.0e8*1e-11,1.0e8*3.0e-8)
yscale("log")
xscale("log")
yscale("log")
xscale("log")
xlabel("alignment parameter GT [cm\$^{-2}\$]",fontsize=FONTSIZE); 
ylabel("model estimation [cm\$^{-2}\$]",fontsize=FONTSIZE) 
xticks(fontsize=FONTSIZE); yticks(fontsize=FONTSIZE); 
# ax2.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax2.yaxis.offsetText.set_size(FONTSIZE)
ax2.xaxis.offsetText.set_size(FONTSIZE)
for fileNameC1s in fileNamesC1sAll
    # get the data
    local α_GT         = zeros(Cdouble,length(dictAllAPE[fileNameC1s]))
    local α_APE        = zeros(Cdouble,length(dictAllAPE[fileNameC1s]))
    local α_APE_approx = zeros(Cdouble,length(dictAllAPE[fileNameC1s]))
    local idxArray = 1
    local xc = 0.0
    for dictTmp in values(dictAllAPE[fileNameC1s])
        α_GT[idxArray]         = 1.0e8dictTmp["α_GT"]
        α_APE[idxArray]        = 1.0e8dictTmp["α_APE"]
        α_APE_approx[idxArray] = 1.0e8dictTmp["α_APE_approx"]
        xc                     = dictTmp["xc"]
        idxArray               = idxArray + 1
    end
    # scatter plots
    ax1.scatter(α_GT,α_APE,label=string("C1s \$x_c\$=",xc,"\$\\mu\$m")) # replace(fileNameC1s[6:end-5],"_"=>" ")
    ax2.scatter(α_GT,α_APE_approx,label=string("C1s \$x_c\$=",xc,"\$\\mu\$m")) # replace(fileNameC1s[6:end-5],"_"=>" ")
end
for fileNameO1s in fileNamesO1sAll
    # get the data
    local α_GT         = zeros(Cdouble,length(dictAllAPE[fileNameO1s]))
    local α_APE        = zeros(Cdouble,length(dictAllAPE[fileNameO1s]))
    local α_APE_approx = zeros(Cdouble,length(dictAllAPE[fileNameO1s]))
    local idxArray = 1
    local xc = 0.0
    for dictTmp in values(dictAllAPE[fileNameO1s])
        α_GT[idxArray]         = 1.0e8dictTmp["α_GT"]
        α_APE[idxArray]        = 1.0e8dictTmp["α_APE"]
        α_APE_approx[idxArray] = 1.0e8dictTmp["α_APE_approx"]
        xc                     = dictTmp["xc"]
        idxArray               = idxArray + 1
    end
    # scatter plots
    ax1.scatter(α_GT,α_APE,label=string("O1s \$x_c\$=",xc,"\$\\mu\$m"))
    ax2.scatter(α_GT,α_APE_approx,label=string("O1s \$x_c\$=",xc,"\$\\mu\$m"))
end
ax1.legend(fontsize=LEGEND_FONTSIZE,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
ax2.legend(fontsize=LEGEND_FONTSIZE,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
ax1.text(-0.1, 0.95, "a)", transform=ax1.transAxes,fontsize=16)
ax2.text(-0.1, 0.95, "b)", transform=ax2.transAxes,fontsize=16)
ax1.text(0.5, 0.25, "\$\\rho\$ GT density", transform=ax1.transAxes,fontsize=16)
ax2.text(0.5, 0.25, "\$\\rho\$ constant density", transform=ax2.transAxes,fontsize=16)

# savefig(string(data_folder,"APE_vs_alignment_parameter_units_gt_IGOR.png"))
# savefig(string(data_folder,"APE_vs_alignment_parameter_units_gt_IGOR.pdf"))
# savefig(string(data_folder,"two_col_APE_vs_alignment_parameter_units_gt_IGOR.png"))
# savefig(string(data_folder,"two_col_APE_vs_alignment_parameter_units_gt_IGOR.pdf"))



# plot background and noise removal
sym_plot_bg = [:hν_650; :hν_958; :hν_1576; :hν_1884]
panel_label = ["a)"; "b)"; "c)"; "d)"]
panel_label_loc_x = [-0.09; -0.14; -0.09; -0.14]
AX_vect = Array{PyPlot.PyObject,1}(undef,4)
figure(figsize=[12, 10])
for i in 1:4
    AX_vect[i] = subplot(2,2,i)
    local Be     = dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["Be"];
    local BeIGOR = dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["Be"];
    AX_vect[i].plot(Be,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["bg"],label="background GT",color="darkblue")
    AX_vect[i].scatter(BeIGOR,dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["bg"],label="background IGOR",color="tab:blue")
    AX_vect[i].scatter(Be,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["S_noisy"],label="data",color="tab:red")
    AX_vect[i].plot(Be,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["Spectrum_OI"],label="SOI",color="darkgreen")
    AX_vect[i].scatter(BeIGOR,dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["Spectrum_OI"],label="SOI IGOR",color="tab:green")

    # plot the true noise, the estimated noise with the proposed method and the estimated noise with IGOR
    AX_vect[i].plot(Be,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["S_noisy"]-dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["bg"]-dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["Spectrum_OI"],label="noise GT",color=color_array[5])
    AX_vect[i].plot(BeIGOR,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["noise estimation"],label="noise estimation",color=color_array[6])
    AX_vect[i].scatter(BeIGOR,dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["S_noisy"]-dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["bg"]-dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["Spectrum_OI"],label="noise IGOR",color=color_array[7])

    
    # plot(Be,dictAllData[plot_sym].Sbg+dictAllData[plot_sym].SpectrumA_1,label="fits GT"); 
    # plot(Be,dictAllData[plot_sym].Snoisy-(dictAllData[plot_sym].Sbg+dictAllData[plot_sym].SpectrumA_1),label="noise GT"); 
    # scatter(Be,dictAllData[plot_sym].Snoisy,label="data")
    xlim(BeIGOR[end],BeIGOR[1])
    AX_vect[i].invert_xaxis();
    xlabel("binding energy [eV]",fontsize=14); 
    ylabel("spectrum [count]",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    AX_vect[i].ticklabel_format(axis="y",style="sci",scilimits=(-1,2),useOffset=true)
    AX_vect[i].yaxis.offsetText.set_size(14)
    AX_vect[i].xaxis.offsetText.set_size(14)
    AX_vect[i].text(0.33, 0.75, string("\$h\\nu\$=",convert(Int64,round(dictAllAPE[fileNamesC1sAll[3]][sym_plot_bg[i]]["hν"])),"[eV]"), transform=AX_vect[i].transAxes,fontsize=16)
    AX_vect[i].text(panel_label_loc_x[i], 0.95, panel_label[i] , transform=AX_vect[i].transAxes,fontsize=16)
    AX_vect[i].legend(fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
end
tight_layout(pad=1.0, w_pad=0.4, h_pad=0.2)





# plot background and noise removal

sym_plot_bg = [:hν_650; :hν_958; :hν_1576; :hν_1884]
panel_label = ["a)"; "b)"; "c)"; "d)"]
if TWO_COLUMN
    panel_label_loc_x = [-0.09; -0.14; -0.09; -0.14].-0.01
else
    panel_label_loc_x = [-0.09; -0.14; -0.09; -0.14]
end
AX_vect = Array{PyPlot.PyObject,1}(undef,4)
if TWO_COLUMN
    figure(figsize=[12, 8])
else
    figure(figsize=[12, 10])
end
if TWO_COLUMN
    FONTSIZE = 16
else
    FONTSIZE = 14
end
if TWO_COLUMN
    LEGEND_FONTSIZE = 14
else
    LEGEND_FONTSIZE = 12
end
LINEWIDTH = 2.5
for i in 1:4
    AX_vect[i] = subplot(2,2,i)
    local Be     = dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["Be"];
    local BeIGOR = dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["Be"];
    AX_vect[i].plot(Be,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["bg"],label="bg. GT",color="darkblue",linewidth=LINEWIDTH)
    AX_vect[i].scatter(BeIGOR,dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["bg"],label="bg. SPANCF",color="tab:blue")
    # AX_vect[i].scatter(Be,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["S_noisy"],label="data",color="tab:red")
    AX_vect[i].plot(Be,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["Spectrum_OI"],label="SOI GT",color="darkgreen",linewidth=LINEWIDTH)
    AX_vect[i].scatter(BeIGOR,dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["Spectrum_OI"],label="SOI SPANCF",color="tab:green")

    # plot the true noise, the estimated noise with the proposed method and the estimated noise with IGOR
    AX_vect[i].plot(Be,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["S_noisy"]-dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["bg"]-dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["Spectrum_OI"],label="noise GT",color="cyan",linewidth=LINEWIDTH) # color_array[5]
    AX_vect[i].plot(BeIGOR,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["noise estimation"],label="noise SVD-based",color=color_array[6],linewidth=LINEWIDTH-0.5)
    AX_vect[i].scatter(BeIGOR,dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["S_noisy"]-dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["bg"]-dictAllIGOR[fileNamesC1sAll[3]][sym_plot_bg[i]]["Spectrum_OI"],label="noise SPANCF",color=color_array[7])

    
    # plot(Be,dictAllData[plot_sym].Sbg+dictAllData[plot_sym].SpectrumA_1,label="fits GT"); 
    # plot(Be,dictAllData[plot_sym].Snoisy-(dictAllData[plot_sym].Sbg+dictAllData[plot_sym].SpectrumA_1),label="noise GT"); 
    # scatter(Be,dictAllData[plot_sym].Snoisy,label="data")
    xlim(BeIGOR[end],BeIGOR[1])
    AX_vect[i].invert_xaxis();
    xlabel("binding energy [eV]",fontsize=FONTSIZE); 
    ylabel("spectrum [count]",fontsize=FONTSIZE) 
    xticks(fontsize=FONTSIZE); yticks(fontsize=FONTSIZE); 
    AX_vect[i].ticklabel_format(axis="y",style="sci",scilimits=(-1,2),useOffset=true)
    AX_vect[i].yaxis.offsetText.set_size(FONTSIZE)
    AX_vect[i].xaxis.offsetText.set_size(FONTSIZE)
    if TWO_COLUMN
        AX_vect[i].text(0.38, 0.85, string("\$h\\nu\$=",convert(Int64,round(dictAllAPE[fileNamesC1sAll[3]][sym_plot_bg[i]]["hν"])),"[eV]"), transform=AX_vect[i].transAxes,fontsize=16)
    else
        AX_vect[i].text(0.33, 0.75, string("\$h\\nu\$=",convert(Int64,round(dictAllAPE[fileNamesC1sAll[3]][sym_plot_bg[i]]["hν"])),"[eV]"), transform=AX_vect[i].transAxes,fontsize=16)
    end
    AX_vect[i].text(panel_label_loc_x[i], 0.95, panel_label[i] , transform=AX_vect[i].transAxes,fontsize=16)
    AX_vect[i].legend(fontsize=LEGEND_FONTSIZE,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
end
tight_layout(pad=1.0, w_pad=0.4, h_pad=0.2)


# savefig(string(data_folder,"background_removal_peak_fit_and_noise_removal_IGOR.png"))
# savefig(string(data_folder,"background_removal_peak_fit_and_noise_removal_IGOR.pdf"))
# savefig(string(data_folder,"two_col_background_removal_peak_fit_and_noise_removal_IGOR.png"))
# savefig(string(data_folder,"two_col_background_removal_peak_fit_and_noise_removal_IGOR.pdf"))


# plot data

sym_plot_bg = [:hν_650; :hν_958; :hν_1576; :hν_1884]
panel_label = ["a)"; "b)"; "c)"; "d)"]
if TWO_COLUMN
    panel_label_loc_x = [-0.11; -0.13; -0.11; -0.13].-0.01
else
    panel_label_loc_x = [-0.11; -0.13; -0.11; -0.13]
end
AX_vect = Array{PyPlot.PyObject,1}(undef,4)
if TWO_COLUMN
    figure(figsize=[12, 8])
else
    figure(figsize=[12, 10])
end
if TWO_COLUMN
    FONTSIZE = 16
else
    FONTSIZE = 14
end
if TWO_COLUMN
    LEGEND_FONTSIZE = 14
else
    LEGEND_FONTSIZE = 12
end
LINEWIDTH = 2.5
for i in 1:4
    AX_vect[i] = subplot(2,2,i)
    local Be     = dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["Be"];
    AX_vect[i].plot(Be,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["bg"],label="background",color="tab:blue",linewidth=LINEWIDTH)
    AX_vect[i].scatter(Be,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["S_noisy"],label="data",color="tab:blue")
    AX_vect[i].plot(Be,dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["Spectrum_OI"]+dictAllPlot[fileNamesC1sAll[3]][sym_plot_bg[i]]["bg"],label="noise free spectrum",color="tab:orange",linewidth=LINEWIDTH)
    
    # plot(Be,dictAllData[plot_sym].Sbg+dictAllData[plot_sym].SpectrumA_1,label="fits GT"); 
    # plot(Be,dictAllData[plot_sym].Snoisy-(dictAllData[plot_sym].Sbg+dictAllData[plot_sym].SpectrumA_1),label="noise GT"); 
    # scatter(Be,dictAllData[plot_sym].Snoisy,label="data")
    # xlim(Be[end],Be[1])
    xlim(272.0,290.0)
    AX_vect[i].invert_xaxis();
    xlabel("binding energy [eV]",fontsize=FONTSIZE); 
    ylabel("spectrum [count]",fontsize=FONTSIZE) 
    xticks(fontsize=FONTSIZE); yticks(fontsize=FONTSIZE); 
    AX_vect[i].ticklabel_format(axis="y",style="sci",scilimits=(-1,2),useOffset=true)
    AX_vect[i].yaxis.offsetText.set_size(FONTSIZE)
    AX_vect[i].xaxis.offsetText.set_size(FONTSIZE)
    AX_vect[i].text(0.1, 0.5, string("\$h\\nu\$=",convert(Int64,round(dictAllAPE[fileNamesC1sAll[3]][sym_plot_bg[i]]["hν"])),"[eV]"), transform=AX_vect[i].transAxes,fontsize=16)
    AX_vect[i].text(panel_label_loc_x[i], 0.95, panel_label[i] , transform=AX_vect[i].transAxes,fontsize=16)
    AX_vect[i].legend(fontsize=LEGEND_FONTSIZE,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
end
tight_layout(pad=1.0, w_pad=0.4, h_pad=0.2)


# savefig(string(data_folder,replace(fileNamesC1sAll[3], "xlsx"=>"png", "/"=>"_")))
# savefig(string(data_folder,replace(fileNamesC1sAll[3], "xlsx"=>"pdf", "/"=>"_")))
# savefig(string(data_folder,"two_col_",replace(fileNamesC1sAll[3], "xlsx"=>"png", "/"=>"_")))
# savefig(string(data_folder,"two_col_",replace(fileNamesC1sAll[3], "xlsx"=>"pdf", "/"=>"_")))

