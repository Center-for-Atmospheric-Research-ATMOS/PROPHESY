##
## load data and meta data 
##

folder_content = readdir(data_folder);
list_match = match.(r"^offcenter_[0-3]$",folder_content)
data_folders = folder_content[list_match.!=nothing]
# folder_tag = data_folders[1]
folder_tag = data_folders[4]

# unit conversion constant (some of the quantities are in μm, some in L and some in Mbarn)
NA = 6.022e23;
κ_simple_units = 1.0e-37*NA; # simplified model
κ_units        = 1.0e-25*NA; # original model

# O1s data
dictAllData,df,symbolDict = dataAndMeta_xlsx2df(string(data_folder,folder_tag,"/new_eal/water/data.xlsx"),r"^hν_[0-9]*",r"meta");
dictAllGeom,symbolDictInv = model_xlsx2df(string(data_folder,folder_tag,"/new_eal/water/model.xlsx"));

# S2p data
if FLAG_0001
    dictAllDataS2p,dfS2p,symbolDictS2p = dataAndMeta_xlsx2df(string(data_folder,folder_tag,"/new_eal/0001/S2p/data.xlsx"),r"^hν_[0-9]*",r"meta")
    dictAllGeomS2p,symbolDictInvS2p = model_xlsx2df(string(data_folder,folder_tag,"/new_eal/0001/S2p/model.xlsx"));
end
if FLAG_0002
    dictAllDataS2p,dfS2p,symbolDictS2p = dataAndMeta_xlsx2df(string(data_folder,folder_tag,"/new_eal/0002/S2p/data.xlsx"),r"^hν_[0-9]*",r"meta")
    dictAllGeomS2p,symbolDictInvS2p = model_xlsx2df(string(data_folder,folder_tag,"/new_eal/0002/S2p/model.xlsx"));
end
if FLAG_0003
    dictAllDataS2p,dfS2p,symbolDictS2p = dataAndMeta_xlsx2df(string(data_folder,folder_tag,"/new_eal/0003/S2p/data.xlsx"),r"^hν_[0-9]*",r"meta")
    dictAllGeomS2p,symbolDictInvS2p = model_xlsx2df(string(data_folder,folder_tag,"/new_eal/0003/S2p/model.xlsx"));
end
if FLAG_0004
    dictAllDataS2p,dfS2p,symbolDictS2p = dataAndMeta_xlsx2df(string(data_folder,folder_tag,"/new_eal/0004/S2p/data.xlsx"),r"^hν_[0-9]*",r"meta")
    dictAllGeomS2p,symbolDictInvS2p = model_xlsx2df(string(data_folder,folder_tag,"/new_eal/0004/S2p/model.xlsx"));
end
if FLAG_0005
    dictAllDataS2p,dfS2p,symbolDictS2p = dataAndMeta_xlsx2df(string(data_folder,folder_tag,"/new_eal/0005/S2p/data.xlsx"),r"^hν_[0-9]*",r"meta")
    dictAllGeomS2p,symbolDictInvS2p = model_xlsx2df(string(data_folder,folder_tag,"/new_eal/0005/S2p/model.xlsx"));
end


# number of measurement
NdataO1s                     = length(dictAllData);
NdataS2p                     = length(dictAllDataS2p);

# sort essential entries (photon energy,attenuation length,peak probability) by increasing photon energy and increasing binding energy mode
sortBy = :hν
thenSortBy = :peak_mode
dfEssS2p = dfS2p |> @orderby(_[sortBy]) |> @thenby(_[thenSortBy]) |> @select(:hν,:λ,:peak_probability)  |> DataFrame;
EssS2p   = Matrix{Cdouble}(dfEssS2p);
dfhν     = dfEssS2p |> @select(:hν) |> @unique() |> DataFrame
dfλ      = dfEssS2p |> @select(:λ) |> @unique() |> DataFrame
hν       = dropdims(Array{Cdouble}(dfhν),dims=2);
λ        = dropdims(Array{Cdouble}(dfλ),dims=2);

y_peak_1 = zeros(Cdouble,NdataS2p);
y_peak_2 = zeros(Cdouble,NdataS2p);
A_tot    = zeros(Cdouble,NdataS2p); # total area under the spectrum (without background, e.g. from the fitted curve)
σnoise   = zeros(Cdouble,NdataS2p);
for i in 1:NdataS2p
    local symData = Symbol(string("hν_",round(Int64,EssS2p[2*(i-1)+1,1])));
    local spectrum = dictAllDataS2p[symData].SpectrumA_1;
    local Ke = dictAllDataS2p[symData].Ke;
    local p1 = EssS2p[2*(i-1)+1,3];
    local p2 = EssS2p[2*(i-1)+2,3];
    dKe = abs(Ke[2]-Ke[1])
    # total area (in the Reimann sense)
    A_tot[i] = dKe*sum(spectrum)
    # partial area (area under each peak)
    y_peak_1[i] = p1*A_tot[i]
    y_peak_2[i] = p2*A_tot[i]
    σnoise[i] = sqrt(var(dictAllDataS2p[symData].Snoisy-dictAllDataS2p[symData].SpectrumA_1-dictAllDataS2p[symData].Sbg))
end

##
## plotting
##
if PLOT_DATA
    ph_sym_list = keys(symbolDict)
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
        ax.text(0.05, 0.5, string("Eph = ",Eph," [eV]"), transform=ax.transAxes,fontsize=14)
    end
    tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)

    ph_sym_listS2p = keys(symbolDictS2p)
    figure(figsize=[12, 10])
    for (i,plot_sym) in zip(1:4,ph_sym_listS2p)
        local Be = dictAllDataS2p[plot_sym].Be;
        local Sliq = dictAllDataS2p[plot_sym].SpectrumA_1;
        local Sbg = dictAllDataS2p[plot_sym].Sbg;
        local Snoisy = dictAllDataS2p[plot_sym].Snoisy;
        local Eph = dictAllDataS2p[plot_sym].hν[1];
        
        local ax = subplot(2,2,i)
        
        scatter(Be,Snoisy,label="Data")
        plot(Be,Sliq,label="Liquid signal")
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
        ax.text(0.1, 0.3, string("Eph = ",Eph," [eV]"), transform=ax.transAxes,fontsize=14)
    end
    tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
end