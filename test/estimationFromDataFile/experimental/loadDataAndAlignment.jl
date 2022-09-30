# find the data files
folder_content = readdir(data_folderC1s);
list_match = match.(r"xlsx$",folder_content)
data_filesC1s = folder_content[list_match.!=nothing]

# C1s
# idx_file = 1; # not enough data
# idx_file = 2; # OK
# idx_file = 3; # not so good, probably problem in the data
# idx_file = 4; # OK

#S2p
# idx_file = 1;  # not enough data, but looks okish
# idx_file = 2; # OK
# idx_file = 3; # not OK
# idx_file = 4; # OK
idx_file = 5;

fileName = data_filesC1s[idx_file][1:end-5];

ρC1s_out  = 0.0                                                                  # [M] concentration of C 1s outside of the sharp edge volume (at the outter boundary)
ρC1s_bulk = 1.0e-3parse(Cdouble,match(r"SDS_[0-9]*mM",fileName).match[5:end-2]); # [M] bulk concentration of C 1s

dictAllData,df_fit,Ndata = dataAndFit_xlsx2df(string(data_folderC1s,fileName,".xlsx"); regData=regData, regFit=regFit, regEph=regEph, sortBy=sortBy, thenSortBy=thenSortBy);

# list all the photon energies
df_Eph = select(df_fit,photon_sym) |> @unique() |> @orderby(_[photon_sym]) |> DataFrame
# put everything in nice boxes
dictPeak = Dict();
for i in df_Eph[!,photon_sym]
    local df = df_fit |> @filter(_[photon_sym]==i) |> DataFrame
    dictPeak[Symbol(string("hν_",i))] = df # (df,names(df))
end


# alignement
α_al_noise = zeros(Cdouble,Ndata);
α_al_ratio = zeros(Cdouble,Ndata);
T          = zeros(Cdouble,Ndata);
Fν         = zeros(Cdouble,Ndata);
σ_C1s      = zeros(Cdouble,Ndata);
H_geom     = zeros(Cdouble,Ndata,Nr);
σ_data     = zeros(Cdouble,Ndata);
y_peak_1   = zeros(Cdouble,Ndata);
y_peak_2   = zeros(Cdouble,Ndata);
y_peak_3   = zeros(Cdouble,Ndata);
λ_all      = zeros(Cdouble,Ndata);

if FLAG_PLOT
    figure(figsize=[12, 10])
end
for i in 1:Ndata
    # get fitting data into arrays
    local plot_sym = Symbol(string("hν_",df_Eph[!,photon_sym][i]));
    local Be       = dictAllData[plot_sym].Wavelength;
    local dKe      = median(abs.(Be[2:end]-Be[1:end-1]))
    local Npeaks, BePeak, σePeak, σePeakL, AePeak, σ_peak = curveFromFit(dictPeak[plot_sym],Be,true;bind_sym=bind_sym, shift_sym=shift_sym, gauss_sym=gauss_sym, area_sym=area_sym, loren_sym=loren_sym);

    # get peak areas
    y_peak_1[i] = dKe*sum(dictAllData[plot_sym].Curve1); # dKe*AePeak[1]
    y_peak_2[i] = dKe*sum(dictAllData[plot_sym].Curve2); # dKe*AePeak[2]
    # y_peak_3[i] = dKe*sum(dictAllData[plot_sym].Curve3); # dKe*AePeak[3] # C1s

    # plot data and fits
    if FLAG_PLOT
        local ax = subplot(2,2,i)
        plot(Be,dictAllData[plot_sym].Raw_spectrum,label="Data")
        plot(Be,dictAllData[plot_sym].Background.+dKe*dropdims(AePeak'*σ_peak,dims=1),label="fitted spectra peaks")
        # plot(Be,dictAllData[plot_sym].Background.+dictAllData[plot_sym].Curve1.+dictAllData[plot_sym].Curve2.+dictAllData[plot_sym].Curve3,label="fitted spectra curves") # C1s
        plot(Be,dictAllData[plot_sym].Background.+dictAllData[plot_sym].Curve1.+dictAllData[plot_sym].Curve2,label="fitted spectra curves") # S2p
        ylim(0.0)
        xlabel("binding energy [eV]",fontsize=14); 
        ylabel("spectrum [count]",fontsize=14) 
        xticks(fontsize=14); yticks(fontsize=14); 
        legend(fontsize=14)
        ax.invert_xaxis()
        ax.text(0.1, 0.5, string("Eph = ",df_Eph[!,photon_sym][i]," [eV]"), transform=ax.transAxes,fontsize=14)
    end
    
    # measurement model
    local λe = 1.0e-3dictPeak[plot_sym][!,:IMFP][2]
    # local σ_all = dictAllData[plot_sym].Curve1.+dictAllData[plot_sym].Curve2.+dictAllData[plot_sym].Curve3; # C1s
    local σ_all = dictAllData[plot_sym].Curve1.+dictAllData[plot_sym].Curve2; # S2p
    σ_all = σ_all/(dKe*sum(σ_all));
    # compute the geomtry factor
    H_geom[i,:],_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λe);
    # estimate the alignment parameter
    local Sbg = dictAllData[plot_sym].Background;
    local S_noisy = dictAllData[plot_sym].Raw_spectrum;
    local ρ = ρC1s_bulk*ones(Cdouble,Nr)
    α_al_noise[i],noise_data = noiseAndParameterEstimation(σ_all,H_geom[i,:],S_noisy,Sbg,ρ)
    # α_al_noise[i]            = α_al_noise[i]/(κ_units*σ_C1s_exp(convert(Cdouble,df_Eph[!,photon_sym][i]))*dictPeak[plot_sym][!,ph_flu_sym][1])
    α_al_noise[i]            = α_al_noise[i]/(κ_units*σ_S2p_exp(convert(Cdouble,df_Eph[!,photon_sym][i]))*dictPeak[plot_sym][!,ph_flu_sym][1])
    α_al_ratio[i]            = dictPeak[plot_sym][!,α_al_sym][1]
    σ_data[i]                = sqrt.(var(noise_data))
    Fν[i]                    = dictPeak[plot_sym][!,ph_flu_sym][1];
    T[i]                     = 1.0;
    # σ_C1s[i]                 = σ_C1s_exp(convert(Cdouble,df_Eph[!,photon_sym][i]));
    σ_C1s[i]                 = σ_S2p_exp(convert(Cdouble,df_Eph[!,photon_sym][i]));
    λ_all[i]                 = dictPeak[plot_sym][!,:IMFP][2];
end

y_data_1 = y_peak_1./(α_al_noise.*T.*Fν.*σ_C1s.*ρC1s_bulk*κ_units) 
y_data_2 = y_peak_2./(α_al_noise.*T.*Fν.*σ_C1s.*ρC1s_bulk*κ_units)
y_data_3 = y_peak_3./(α_al_noise.*T.*Fν.*σ_C1s.*ρC1s_bulk*κ_units)
σ_noise  = σ_data./(α_al_noise.*T.*Fν.*σ_C1s.*ρC1s_bulk*κ_units)

SNR_1 = (y_data_1./σ_noise).^2;
SNR_2 = (y_data_2./σ_noise).^2;
SNR_3 = (y_data_3./σ_noise).^2;

if FLAG_PLOT
    tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
    if FLAG_SAVE_PLOT
        savefig(string(fileName,"_plot.pdf"))
        savefig(string(fileName,"_plot.png"))
    end
end
