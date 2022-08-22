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
using XLSX # CSV does not deal with multiple sheets
using DataFrames
using Query
using DataValues

# scientific package from the official Julia repositories
using LinearAlgebra
using StatsBase
using Interpolations

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)
using XPSsampling

PLOT_FIG  = false
SAVE_FIG  = false
SAVE_DATA = false

data_folder = "../data/cylinder_radius_10.0/peak_shift/eal_5_restricted_range/"

# TODO: for each center offset, for the water case, load the data, the model and compare the estimated alignment parameter with the true parameter
# data selection
# FLAG_OFF_CENTER_0 = true;
# FLAG_OFF_CENTER_1 = false;
# FLAG_OFF_CENTER_2 = false;
# FLAG_OFF_CENTER_3 = false;

# # selection of the profile (one must be true and the others false)
# FLAG_0001  = false
# FLAG_0002  = false
# FLAG_0003  = false
# FLAG_0004  = false
# FLAG_WATER = true

# if FLAG_OFF_CENTER_0
#     data_folder = string(data_folder,"offcenter_0/")
# elseif FLAG_OFF_CENTER_1
#     data_folder = string(data_folder,"offcenter_1/")
# elseif FLAG_OFF_CENTER_2
#     data_folder = string(data_folder,"offcenter_2/")
# else
#     data_folder = string(data_folder,"offcenter_3/")
# end

# if FLAG_0001
#     data_folder = string(data_folder,"0001/")
# elseif FLAG_0002
#     data_folder = string(data_folder,"0002/")
# elseif FLAG_0003
#     data_folder = string(data_folder,"0003/")
# elseif FLAG_0004
#     data_folder = string(data_folder,"0004/")
# else
#     data_folder = string(data_folder,"water/")
# end


function dataAndMeta_xlsx2df(fileName::String,regData::Regex,regMeta::Regex)
    local xf_data = XLSX.readxlsx(fileName);
    local xf_sheet_names = XLSX.sheetnames(xf_data);
    local list_match_data = match.(regData,xf_sheet_names);
    local list_match_meta = match.(regMeta,xf_sheet_names);
    local xf_data_sheet_names = xf_sheet_names[list_match_data.!=nothing];
    local xf_meta_sheet_names = xf_sheet_names[list_match_meta.!=nothing];
    # the data
    dictAllData = Dict();
    symbolDict = Dict();
    for xf_name in xf_data_sheet_names
        local x = XLSX.getdata(xf_data[xf_name])
        local dataPairs = Array{Pair{String,Vector{Cdouble}}}(undef,size(x,2));
        for j in 1:size(x,2)
            dataPairs[j] = (x[1,j] => convert(Array{Cdouble,1},x[2:end,j]))
        end
        df = DataFrame(dataPairs);
        dictAllData[Symbol(xf_name)] = df 
        symbolDict[Symbol(xf_name)] = Symbol(string("λe_",df[!,:λ][1]))
    end
    # meta data
    local dictAllMeta = Dict();
    for xf_name in xf_meta_sheet_names
        local x = XLSX.getdata(xf_data[xf_name])
        local dataPairs = Array{Pair{String,Vector{Cdouble}}}(undef,size(x,2));
        for j in 1:size(x,2)
            dataPairs[j] = (x[1,j] => convert(Array{Cdouble,1},x[2:end,j]))
        end
        df = DataFrame(dataPairs);
        dictAllMeta[Symbol(xf_name)] = df 
    end

    local dfMeta = DataFrame();
    for (_,df) in dictAllMeta
        append!(dfMeta,df);
    end
    df = dfMeta |> @orderby(_[:hν]) |> DataFrame

    # return
    dictAllData,df,symbolDict
end

function model_xlsx2df(fileName::String)
    local xf_model = XLSX.readxlsx(fileName);
    local xf_sheet_names = XLSX.sheetnames(xf_model);
    # the data
    dictAllGeom = Dict();
    symbolDict = Dict();
    for xf_name in xf_sheet_names
        local x = XLSX.getdata(xf_model[xf_name])
        local dataPairs = Array{Pair{String,Union{Vector{Cdouble},Vector{String}}}}(undef,size(x,2));
        for j in 1:size(x,2)
            if (x[1,j]=="model")
                dataPairs[j] = (x[1,j] => string.(x[2:end,j]))
            else
                dataPairs[j] = (x[1,j] => convert(Array{Cdouble,1},x[2:end,j]))
            end
        end
        df = DataFrame(dataPairs);
        dictAllGeom[Symbol(xf_name)] = df 
        symbolDict[Symbol(xf_name)] = Symbol(string("hν_",round(Int64,df[!,:hν][1])))
    end

    # return
    dictAllGeom,symbolDict
end

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
idx = 1
for folder_tag in data_folders
    # local dictAllData,df,symbolDict = dataAndMeta_xlsx2df(string(data_folder,folder_tag,"water/data.xlsx"),r"^hν_[0-9]*",r"meta");
    # local dictAllGeom,symbolDictInv = model_xlsx2df(string(data_folder,folder_tag,"water/model.xlsx"));
    local dictAllData,df,symbolDict = dataAndMeta_xlsx2df(string(data_folder,folder_tag,"/water/data.xlsx"),r"^hν_[0-9]*",r"meta");
    local dictAllGeom,symbolDictInv = model_xlsx2df(string(data_folder,folder_tag,"/water/model.xlsx"));
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
        local ρ = 55.49*ones(Cdouble,Nr) # dictAllGeom[symbolDict[sym_data]][!,:ρ] # it works way too well with the true concentration profile
        α_al_noise[idx],_    = noiseAndParameterEstimation(σ_cs_liq,H_liq,S_liq_noisy,Sbg,ρ)
        local Tj = dictAllData[sym_data][!,:T][1];
        local Fνj = dictAllData[sym_data][!,:F][1];
        local σ_tot = dictAllData[sym_data][!,:σ_tot][1];
        α_al_noise[idx] = α_al_noise[idx]/(Tj*Fνj*σ_tot)
        # τ_al_noise_gt[idx],_ = noiseAndParameterEstimation(dictAllData[symbol_h][1][:σ_cs_dens],dictAllGeom[simbol_λ][1][:H],
        #                     Array{Cdouble,1}(dictAllData[symbol_h][1][:Snoisy]),dictAllData[symbol_h][1][:Sbg],ρA_1)

        idx = idx + 1
    end

    # plot
    if false
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

Ndata = 5;
α_gt_mean    = zeros(Cdouble,length(data_folders));
α_ratio_mean = zeros(Cdouble,length(data_folders));
α_al_mean    = zeros(Cdouble,length(data_folders));
for i in 1: length(data_folders)
    α_gt_mean[i]    = mean(α_gt[(i-1)*Ndata+1:i*Ndata])
    α_ratio_mean[i] = mean(α_ratio[(i-1)*Ndata+1:i*Ndata])
    α_al_mean[i]    = mean(α_al_noise[(i-1)*Ndata+1:i*Ndata])
end


# α_gt_mean    = [mean(α_gt[1:Ndata]); mean(α_gt[Ndata+1:2Ndata]); mean(α_gt[2Ndata+1:3Ndata]); mean(α_gt[3Ndata+1:4Ndata])];
# α_ratio_mean = [mean(α_ratio[1:Ndata]); mean(α_ratio[Ndata+1:2Ndata]); mean(α_ratio[2Ndata+1:3Ndata]); mean(α_ratio[3Ndata+1:4Ndata])];



figure(figsize=[10, 5]); # scatter(1.0e9α_gt,α_ratio)
ax = subplot(122)
scatter(1.0e10α_gt_mean,α_ratio_mean,color="tab:orange",label="liquid/vapor ratio") # .^2.5
scatter(1.0e10α_gt_mean,1.0e10α_al_mean,color="tab:green",label="noise estimation [x\$10^{10}\$]") 
xlabel("GT [x\$10^{9}\$]",fontsize=14); 
ylabel("estimation",fontsize=14) 
xticks(fontsize=14); yticks(fontsize=14); 
ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax.yaxis.offsetText.set_size(14)
ax.xaxis.offsetText.set_size(14)
legend(fontsize=14)

ax = subplot(121)
scatter(xc_off[1:5:end],1.0e10α_gt_mean,label="GT [x\$10^{10}\$]")
scatter(xc_off[1:5:end],α_ratio_mean,label="liquid/vapor area ratio")
scatter(xc_off[1:5:end],1.0e10α_al_mean,label="noise estimation [x\$10^{10}\$]")
xlabel("horizontal offset [\$\\mu\$m]",fontsize=14); 
ylabel("alignment parameter",fontsize=14) 
xticks(fontsize=14); yticks(fontsize=14); 
ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax.yaxis.offsetText.set_size(14)
ax.xaxis.offsetText.set_size(14)
legend(fontsize=14)
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)


savefig(string(data_folder,"liquid_vapor_area_ratio_and_noise_estimation_vs_alignment_parameter.png"))
savefig(string(data_folder,"liquid_vapor_area_ratio_and_noise_estimation_vs_alignment_parameter.pdf"))

# figure(); scatter(collect(1:length(α_al_noise)),1.0e9α_al_noise)

df_α =  DataFrame( "α_gt"=> α_gt_mean, "α_ratio" =>  α_ratio_mean, "x_offset" => xc_off[1:5:end]);
# XLSX.writetable(string(data_folder,"/alignment_list.xlsx"), collect(eachcol(df_α)), names(df_α));

#TODO: run the alignment parameter estimation form the variance of the noise


if false
    ph_sym_list = keys(symbolDict)

    # for i in df_Eph[!,Symbol("Photon energy")]
    # zip(1:4,Int64.(round.(collect(LinRange(1,Ndata,4)))))
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
    # savefig(string(fileName,"_plot.pdf"))
    # savefig(string(fileName,"_plot.png"))
end

if false
    # TODO: find the data and the fitted cross sections!
    

    
    # # sort sheets by increasing photon energy
    # idx_hν = findfirst(XLSX.getdata(xf_data[xf_data_sheet_names[1]])[1,:].=="hν");
    # hν_list = zeros(Cdouble,Ndata);
    # for i in 1:Ndata
    #     hν_list[i] = XLSX.getdata(xf_data[xf_data_sheet_names[i]])[2,idx_hν]
    # end
    # idx_perm_data = sortperm(hν_list)
    # xf_data_sheet_names = xf_data_sheet_names[idx_perm_data];

    


    x_fit = XLSX.getdata(xf_data[xf_fit_data_sheet_names[1]]);
    # remove missing columns
    x_fit = x_fit[:,broadcast(~,(ismissing.(x_fit[1, :])))]
    df_fit = DataFrame([col for col in eachcol(x_fit[2:end, :])], Symbol.(x_fit[1, :]))
    # remove missing rows
    df_fit = dropmissing(df_fit, disallowmissing=true)
    df_fit = df_fit |> @orderby(_[Symbol("Photon energy")]) |> @thenby(_[Symbol("Binding energy")]) |> DataFrame

    # disallowmissing!(dictAllData[Symbol("Eph=650eV")])

    # TODO: plot spectra  and signal of interest




    # some regex tests
    # str_test = "DataValue{Any}(\"570.103805\")"
    # str_test_int = "DataValue{Any}(\"1884\")"
    # match.(r"[0-9]*\.[0-9]*",str_test)
    # match.(r"[0-9]*",str_test_int)
    # play around with querries and conversions from DataValue{Any} to Cdouble using regular expressions


    # dictAllData[Symbol("Eph=650eV")] |> @filter(_.Spectrum<86.0)

    # import Base.convert
    function convertFloat(a::DataValue{Any})
        # parse(T,match(r"[0-9]*\.[0-9]*",string(a)).match)
        # parse.(match.(r"\"([0-9]*\.[0-9]*|[0-9]*)\"",string(a)).match[2:end-1])
        parse(Cdouble,match(r"[0-9]*\.[0-9]*",string(a)).match)
    end
    function convertInt(a::DataValue{Any})
        # println(a)
        # println(match(r"\([0-9]*\)",string(a)).match[2:end-1])
        parse(Int64,match(r"\([0-9]*\)",string(a)).match[2:end-1])
    end

    function isless(a::DataValue{Any}, b::Float64)
        convertFloat(a)<b
    end
    function isless(a::DataValue{Any}, b::Int64)
        convertInt(a)<b
    end
    function ismore(a::DataValue{Any}, b::Float64)
        convertFloat(a)>b
    end
    function ismore(a::DataValue{Any}, b::Int64)
        convertInt(a)>b
    end

    # df_fit |> @filter(isless(_.Area,400.0))
    # dff = df_fit |> @filter(isless(_.Area,400.0)) |> DataFrame
    df_fit |> @filter(isless(_[Symbol("Area")],400.0)) |> DataFrame
    df_fit |> @filter(convertInt(_[Symbol("Photon energy")])==1884)

    # list all the photon energies
    df_Eph = select(df_fit,Symbol("Photon energy")) |> @unique() |> @orderby(_[Symbol("Photon energy")]) |> DataFrame
    # put everything in nice boxes
    dictPeak = Dict();
    for i in df_Eph[!,Symbol("Photon energy")]
        local df = df_fit |> @filter(_[Symbol("Photon energy")]==i) |> DataFrame
        dictPeak[Symbol(string("hν_",i))] = df # (df,names(df))
    end


    # for i in df_Eph[!,Symbol("Photon energy")]
    figure(figsize=[12, 10])
    for i in 1:Ndata
        local plot_sym = Symbol(string("hν_",df_Eph[!,Symbol("Photon energy")][i]));
        local Be = dictAllData[plot_sym].Be;
        local dKe = median(abs.(Be[2:end]-Be[1:end-1]))
        local Be1 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Binding energy")][1])
        local Be2 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Binding energy")][2])
        local Be3 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Binding energy")][3])
        local σe1 = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("Total width")][1])
        local σe2 = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("Total width")][2])
        local σe3 = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("Total width")][3])
        local Ae1 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Area")][1])
        local Ae2 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Area")][2])
        local Ae3 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Area")][3])

        local σ_peak_1 = (1.0/sqrt(2π*σe1^2))*exp.(-0.5*((Be.-Be1)/σe1).^2)
        local σ_peak_2 = (1.0/sqrt(2π*σe2^2))*exp.(-0.5*((Be.-Be2)/σe2).^2)
        local σ_peak_3 = (1.0/sqrt(2π*σe3^2))*exp.(-0.5*((Be.-Be3)/σe3).^2)
        
        local ax = subplot(2,2,i)
        # title(string("Eph = ",df_Eph[!,Symbol("Photon energy")][i]," [eV]"),fontsize=14)
        plot(Be,dictAllData[plot_sym].Spectrum,label="Data")
        plot(Be,dKe*(Ae1*σ_peak_1+Ae2*σ_peak_2+Ae3*σ_peak_3),label="fitted spectra")
        # plot(Be,Ae1*σ_peak_1+Ae2*σ_peak_2+Ae3*σ_peak_3)
        ylim(0.0)
        xlabel("binding energy [eV]",fontsize=14); 
        ylabel("spectrum [count]",fontsize=14) 
        xticks(fontsize=14); yticks(fontsize=14); 
        legend(fontsize=14)
        ax.invert_xaxis()
        ax.text(0.1, 0.5, string("Eph = ",df_Eph[!,Symbol("Photon energy")][i]," [eV]"), transform=ax.transAxes,fontsize=14)
    end
    tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
    # savefig(string(fileName,"_plot.pdf"))
    # savefig(string(fileName,"_plot.png"))


end


if false

    toto = [5.40000000000000E+0002  8.26188850402832E+0000  8.25107097625732E+0000
            5.39950000000000E+0002  8.06315422058105E+0000  8.23270416259766E+0000
            5.39900000000000E+0002  8.24601459503174E+0000  8.41569614410400E+0000
            5.39850000000000E+0002  8.16820430755615E+0000  8.08174800872803E+0000
            5.39800000000000E+0002  8.25739955902100E+0000  7.39145565032959E+0000
            5.39750000000000E+0002  8.37970161437988E+0000  7.29849338531494E+0000
            5.39700000000000E+0002  7.69463205337524E+0000  7.31126785278320E+0000
            5.39650000000000E+0002  7.09918546676636E+0000  6.89618396759033E+0000
            5.39600000000000E+0002  7.46637010574341E+0000  6.69278287887573E+0000
            5.39550000000000E+0002  7.29692935943603E+0000  7.05577611923218E+0000
            5.39500000000000E+0002  7.19886445999145E+0000  7.27503299713135E+0000
            5.39450000000000E+0002  6.72597122192383E+0000  6.42251920700073E+0000
            5.39400000000000E+0002  6.45564842224121E+0000  6.31265640258789E+0000
            5.39350000000000E+0002  6.33080673217773E+0000  6.10994386672974E+0000
            5.39300000000000E+0002  6.03901004791260E+0000  5.97072029113769E+0000
            5.39250000000000E+0002  6.23915100097656E+0000  5.84074544906616E+0000
            5.39200000000000E+0002  6.26766777038574E+0000  6.30377483367920E+0000
            5.39150000000000E+0002  5.89980316162109E+0000  5.43377733230591E+0000
            5.39100000000000E+0002  5.30004501342773E+0000  5.16163587570190E+0000
            5.39050000000000E+0002  5.26628017425537E+0000  5.38411808013916E+0000
            5.39000000000000E+0002  5.15499877929687E+0000  4.61934995651245E+0000
            5.38950000000000E+0002  4.85021352767944E+0000  4.40010452270508E+0000
            5.38900000000000E+0002  4.54656648635864E+0000  4.13062858581543E+0000
            5.38850000000000E+0002  4.36090230941772E+0000  4.24173641204834E+0000
            5.38800000000000E+0002  3.81183791160583E+0000  3.84136915206909E+0000
            5.38750000000000E+0002  3.71233510971069E+0000  3.59217429161072E+0000
            5.38700000000000E+0002  3.44495677947998E+0000  3.52073574066162E+0000
            5.38650000000000E+0002  3.33858036994934E+0000  3.43424820899963E+0000
            5.38600000000000E+0002  3.31817722320557E+0000  3.17063045501709E+0000
            5.38550000000000E+0002  2.99812555313110E+0000  3.04732656478882E+0000
            5.38500000000000E+0002  2.89893174171448E+0000  3.08824372291565E+0000
            5.38450000000000E+0002  2.62877058982849E+0000  2.90934705734253E+0000
            5.38400000000000E+0002  2.57403302192688E+0000  2.68774771690369E+0000
            5.38350000000000E+0002  2.48610091209412E+0000  2.53858399391174E+0000
            5.38300000000000E+0002  2.52875971794128E+0000  2.46057200431824E+0000
            5.38250000000000E+0002  2.57645583152771E+0000  2.44990563392639E+0000
            5.38200000000000E+0002  2.32161808013916E+0000  2.31495642662048E+0000
            5.38150000000000E+0002  2.25885224342346E+0000  2.07874584197998E+0000
            5.38100000000000E+0002  2.08349418640137E+0000  2.22382712364197E+0000
            5.38050000000000E+0002  2.01342582702637E+0000  1.92315137386322E+0000
            5.38000000000000E+0002  1.84464025497437E+0000  1.86253595352173E+0000
            5.37950000000000E+0002  1.76427006721497E+0000  1.90676856040955E+0000
            5.37900000000000E+0002  1.98506546020508E+0000  1.64297091960907E+0000
            5.37850000000000E+0002  1.87575423717499E+0000  1.65947604179382E+0000
            5.37800000000000E+0002  1.69386613368988E+0000  1.72593808174133E+0000
            5.37750000000000E+0002  1.69524991512299E+0000  1.85407876968384E+0000
            5.37700000000000E+0002  1.50657320022583E+0000  1.72411561012268E+0000
            5.37650000000000E+0002  1.48299515247345E+0000  1.67904829978943E+0000
            5.37600000000000E+0002  1.73168599605560E+0000  1.70582818984985E+0000
            5.37550000000000E+0002  1.96001029014587E+0000  1.68904697895050E+0000
            5.37500000000000E+0002  1.71974897384644E+0000  1.54475772380829E+0000
            5.37450000000000E+0002  1.62870836257935E+0000  1.44112157821655E+0000
            5.37400000000000E+0002  1.63667619228363E+0000  1.57152974605560E+0000
            5.37350000000000E+0002  1.64352047443390E+0000  1.54455482959747E+0000
            5.37300000000000E+0002  1.79427194595337E+0000  1.66192364692688E+0000
            5.37250000000000E+0002  1.63601124286652E+0000  1.74302840232849E+0000
            5.37200000000000E+0002  1.53865563869476E+0000  1.69476950168610E+0000
            5.37150000000000E+0002  1.61197149753571E+0000  1.65572130680084E+0000
            5.37100000000000E+0002  1.64341950416565E+0000  1.72806203365326E+0000
            5.37050000000000E+0002  1.71643018722534E+0000  1.67598295211792E+0000
            5.37000000000000E+0002  1.80103266239166E+0000  1.63081550598145E+0000
            5.36950000000000E+0002  1.83639562129974E+0000  1.67951643466949E+0000
            5.36900000000000E+0002  1.74036419391632E+0000  1.62384200096130E+0000
            5.36850000000000E+0002  1.63200187683105E+0000  1.81717514991760E+0000
            5.36800000000000E+0002  1.73597514629364E+0000  1.66864573955536E+0000
            5.36750000000000E+0002  1.78348064422607E+0000  1.66303789615631E+0000
            5.36700000000000E+0002  1.88907384872437E+0000  1.61987578868866E+0000
            5.36650000000000E+0002  1.82160401344299E+0000  1.65241050720215E+0000
            5.36600000000000E+0002  1.79667985439301E+0000  1.91163110733032E+0000
            5.36550000000000E+0002  1.93657135963440E+0000  2.05275583267212E+0000
            5.36500000000000E+0002  1.78435480594635E+0000  2.16812539100647E+0000
            5.36450000000000E+0002  1.76684260368347E+0000  2.09257912635803E+0000
            5.36400000000000E+0002  2.08937358856201E+0000  2.04358768463135E+0000
            5.36350000000000E+0002  2.05039238929749E+0000  1.87699806690216E+0000
            5.36300000000000E+0002  1.92212402820587E+0000  2.15857744216919E+0000
            5.36250000000000E+0002  1.96474969387054E+0000  2.09885764122009E+0000
            5.36200000000000E+0002  2.18087291717529E+0000  2.16877031326294E+0000
            5.36150000000000E+0002  2.26230812072754E+0000  2.25971245765686E+0000
            5.36100000000000E+0002  2.14952015876770E+0000  2.42560601234436E+0000
            5.36050000000000E+0002  2.19217658042908E+0000  2.45444154739380E+0000
            5.36000000000000E+0002  1.98879206180573E+0000  2.47353911399841E+0000
            5.35950000000000E+0002  2.35838270187378E+0000  2.44890165328979E+0000
            5.35900000000000E+0002  2.47570061683655E+0000  2.70945644378662E+0000
            5.35850000000000E+0002  2.56155133247375E+0000  2.92019104957581E+0000
            5.35800000000000E+0002  2.70251059532166E+0000  2.64803600311279E+0000
            5.35750000000000E+0002  2.98957610130310E+0000  2.90046095848084E+0000
            5.35700000000000E+0002  3.03846669197082E+0000  3.25938820838928E+0000
            5.35650000000000E+0002  3.45357108116150E+0000  3.40209770202637E+0000
            5.35600000000000E+0002  3.40588188171387E+0000  3.64286351203918E+0000
            5.35550000000000E+0002  3.59390521049500E+0000  4.14451360702515E+0000
            5.35500000000000E+0002  4.13762044906616E+0000  4.41114044189453E+0000
            5.35450000000000E+0002  4.79208898544312E+0000  4.83177947998047E+0000
            5.35400000000000E+0002  5.22579002380371E+0000  5.65952539443970E+0000
            5.35350000000000E+0002  5.76665306091309E+0000  6.40190267562866E+0000
            5.35300000000000E+0002  6.78479814529419E+0000  7.10252904891968E+0000
            5.35250000000000E+0002  7.12573862075806E+0000  7.89675092697144E+0000
            5.35200000000000E+0002  8.12331485748291E+0000  9.33777236938477E+0000
            5.35150000000000E+0002  9.77906131744385E+0000  1.04814157485962E+0001
            5.35100000000000E+0002  1.14153251647949E+0001  1.22417268753052E+0001
            5.35050000000000E+0002  1.36651716232300E+0001  1.49302167892456E+0001
            5.35000000000000E+0002  1.60214881896973E+0001  1.75617885589600E+0001
            5.34950000000000E+0002  1.89288921356201E+0001  2.06696205139160E+0001
            5.34900000000000E+0002  2.21548442840576E+0001  2.46658515930176E+0001
            5.34850000000000E+0002  2.57045726776123E+0001  2.90967330932617E+0001
            5.34800000000000E+0002  2.94933910369873E+0001  3.29043579101562E+0001
            5.34750000000000E+0002  3.40315628051758E+0001  3.71354598999023E+0001
            5.34700000000000E+0002  3.97033500671387E+0001  4.21696434020996E+0001
            5.34650000000000E+0002  4.45704536437988E+0001  4.82840690612793E+0001
            5.34600000000000E+0002  4.97135887145996E+0001  5.45140762329102E+0001
            5.34550000000000E+0002  5.64068336486816E+0001  6.09708824157715E+0001
            5.34500000000000E+0002  6.20914154052734E+0001  6.62926025390625E+0001
            5.34450000000000E+0002  6.86242218017578E+0001  7.25934448242188E+0001
            5.34400000000000E+0002  7.49871520996094E+0001  7.80492401123047E+0001
            5.34350000000000E+0002  8.25041961669922E+0001  8.42156295776367E+0001
            5.34300000000000E+0002  8.87792816162109E+0001  9.08489456176758E+0001
            5.34250000000000E+0002  9.30383148193359E+0001  9.67926483154297E+0001
            5.34200000000000E+0002  9.84802474975586E+0001  1.02579277038574E+0002
            5.34150000000000E+0002  1.03385238647461E+0002  1.05608634948730E+0002
            5.34100000000000E+0002  1.08479759216309E+0002  1.08240440368652E+0002
            5.34050000000000E+0002  1.10838249206543E+0002  1.10717597961426E+0002
            5.34000000000000E+0002  1.10129074096680E+0002  1.10712158203125E+0002
            5.33950000000000E+0002  1.08838363647461E+0002  1.06278594970703E+0002
            5.33900000000000E+0002  1.06031082153320E+0002  1.04264564514160E+0002
            5.33850000000000E+0002  1.03303070068359E+0002  9.99715118408203E+0001
            5.33800000000000E+0002  9.73491744995117E+0001  9.38504333496094E+0001
            5.33750000000000E+0002  9.28102340698242E+0001  8.82932815551758E+0001
            5.33700000000000E+0002  8.71233673095703E+0001  8.17279052734375E+0001
            5.33650000000000E+0002  7.98241882324219E+0001  7.71325836181641E+0001
            5.33600000000000E+0002  7.58903198242188E+0001  7.40518264770508E+0001
            5.33550000000000E+0002  7.26585922241211E+0001  7.16459655761719E+0001
            5.33500000000000E+0002  7.12925109863281E+0001  6.92459640502930E+0001
            5.33450000000000E+0002  7.13373565673828E+0001  7.13161773681641E+0001
            5.33400000000000E+0002  6.99936904907227E+0001  7.39482955932617E+0001
            5.33350000000000E+0002  7.38861007690430E+0001  7.61587371826172E+0001
            5.33300000000000E+0002  7.73641891479492E+0001  8.07706604003906E+0001
            5.33250000000000E+0002  8.18825912475586E+0001  8.55679092407226E+0001
            5.33200000000000E+0002  8.82666625976563E+0001  9.13942794799805E+0001
            5.33150000000000E+0002  9.38216018676758E+0001  1.00043563842773E+0002
            5.33100000000000E+0002  1.00556289672852E+0002  1.06078575134277E+0002
            5.33050000000000E+0002  1.09885993957520E+0002  1.13100448608398E+0002
            5.33000000000000E+0002  1.16691192626953E+0002  1.22204360961914E+0002
            5.32950000000000E+0002  1.24850402832031E+0002  1.30415969848633E+0002
            5.32900000000000E+0002  1.33183898925781E+0002  1.37967453002930E+0002
            5.32850000000000E+0002  1.38799255371094E+0002  1.45694580078125E+0002
            5.32800000000000E+0002  1.45048431396484E+0002  1.54814041137695E+0002
            5.32750000000000E+0002  1.54416320800781E+0002  1.63508728027344E+0002
            5.32700000000000E+0002  1.63638290405273E+0002  1.72474838256836E+0002
            5.32650000000000E+0002  1.73162933349609E+0002  1.80452758789062E+0002
            5.32600000000000E+0002  1.81506759643555E+0002  1.85768417358398E+0002
            5.32550000000000E+0002  1.86642105102539E+0002  1.92832962036133E+0002
            5.32500000000000E+0002  1.94276977539062E+0002  1.98542648315430E+0002
            5.32450000000000E+0002  2.00454925537109E+0002  2.01730789184570E+0002
            5.32400000000000E+0002  2.04941482543945E+0002  2.06169128417969E+0002
            5.32350000000000E+0002  2.06338333129883E+0002  2.11825332641602E+0002
            5.32300000000000E+0002  2.12141311645508E+0002  2.14884597778320E+0002
            5.32250000000000E+0002  2.14302581787109E+0002  2.18497390747070E+0002
            5.32200000000000E+0002  2.14283874511719E+0002  2.19166366577148E+0002
            5.32150000000000E+0002  2.14637802124023E+0002  2.17274215698242E+0002
            5.32100000000000E+0002  2.15287719726562E+0002  2.16962265014648E+0002
            5.32050000000000E+0002  2.13017272949219E+0002  2.15331222534180E+0002
            5.32000000000000E+0002  2.10224334716797E+0002  2.10906646728516E+0002
            5.31950000000000E+0002  2.05123428344727E+0002  2.06114974975586E+0002
            5.31900000000000E+0002  2.04079666137695E+0002  1.99307159423828E+0002
            5.31850000000000E+0002  1.99370605468750E+0002  1.96053665161133E+0002
            5.31800000000000E+0002  1.90427856445313E+0002  1.88518630981445E+0002
            5.31750000000000E+0002  1.83916564941406E+0002  1.83971237182617E+0002
            5.31700000000000E+0002  1.77990737915039E+0002  1.74145507812500E+0002
            5.31650000000000E+0002  1.69175399780273E+0002  1.66644149780273E+0002
            5.31600000000000E+0002  1.61372009277344E+0002  1.56585906982422E+0002
            5.31550000000000E+0002  1.52214294433594E+0002  1.49176834106445E+0002
            5.31500000000000E+0002  1.42132492065430E+0002  1.38324569702148E+0002
            5.31450000000000E+0002  1.32864929199219E+0002  1.28521408081055E+0002
            5.31400000000000E+0002  1.23651878356934E+0002  1.18944931030273E+0002
            5.31350000000000E+0002  1.14578018188477E+0002  1.09080825805664E+0002
            5.31300000000000E+0002  1.05348838806152E+0002  1.00424880981445E+0002
            5.31250000000000E+0002  9.69440612792969E+0001  9.31253204345703E+0001
            5.31200000000000E+0002  9.01242065429688E+0001  8.54708633422852E+0001
            5.31150000000000E+0002  8.07712554931641E+0001  7.63154983520508E+0001
            5.31100000000000E+0002  7.36177215576172E+0001  6.89653549194336E+0001
            5.31050000000000E+0002  6.66551361083984E+0001  6.15665893554687E+0001
            5.31000000000000E+0002  5.93629531860351E+0001  5.46748161315918E+0001
            5.30950000000000E+0002  5.23209304809570E+0001  4.96983108520508E+0001
            5.30900000000000E+0002  4.64637527465820E+0001  4.29693641662598E+0001
            5.30850000000000E+0002  4.13088150024414E+0001  3.80141448974609E+0001
            5.30800000000000E+0002  3.65436096191406E+0001  3.35670127868652E+0001
            5.30750000000000E+0002  3.19180431365967E+0001  2.92353439331055E+0001
            5.30700000000000E+0002  2.75817604064941E+0001  2.51611785888672E+0001
            5.30650000000000E+0002  2.40366439819336E+0001  2.15736236572266E+0001
            5.30600000000000E+0002  2.14908638000488E+0001  1.85600357055664E+0001
            5.30550000000000E+0002  1.94003696441650E+0001  1.62392063140869E+0001
            5.30500000000000E+0002  1.69534606933594E+0001  1.40005016326904E+0001
            5.30450000000000E+0002  1.42176570892334E+0001  1.26679611206055E+0001
            5.30400000000000E+0002  1.18929929733276E+0001  1.09349575042725E+0001
            5.30350000000000E+0002  1.06784296035767E+0001  9.59691429138184E+0000
            5.30300000000000E+0002  9.88136863708496E+0000  8.39709377288818E+0000
            5.30250000000000E+0002  8.48249244689941E+0000  6.98136329650879E+0000
            5.30200000000000E+0002  7.50316381454468E+0000  6.67924261093140E+0000
            5.30150000000000E+0002  6.40300941467285E+0000  5.88806009292602E+0000
            5.30100000000000E+0002  5.48972511291504E+0000  5.20340824127197E+0000
            5.30050000000000E+0002  5.13008832931519E+0000  4.99361991882324E+0000
            5.30000000000000E+0002  4.48657512664795E+0000  4.58317041397095E+0000]


    figure()
    plot(toto[:,1],toto[:,2])
    plot(toto[:,1],toto[:,2]+sqrt.(toto[:,2]).*randn(length(toto[:,2])))

    figure()
    plot(toto[:,1],100toto[:,2])
    plot(toto[:,1],100toto[:,2]+10sqrt.(toto[:,2]).*randn(length(toto[:,2])))

    figure()
    plot(toto[:,1],37.642toto[:,2])
    plot(toto[:,1],37.642toto[:,2]+sqrt.(37.642toto[:,2]).*randn(length(toto[:,2])))

    figure()
    plot(toto[:,1],118toto[:,2])
    plot(toto[:,1],118toto[:,2]+sqrt.(118toto[:,2]).*randn(length(toto[:,2])))

    figure()
    plot(toto[:,1],1000toto[:,2])
    plot(toto[:,1],1000toto[:,2]+sqrt.(1000toto[:,2]).*randn(length(toto[:,2])))
end
