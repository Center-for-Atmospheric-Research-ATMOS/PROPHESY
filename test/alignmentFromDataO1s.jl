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

#TODO: wait for the shifts
data_folderO1s = "../data/TK/O1s/"

# find the data files
folder_content = readdir(data_folderO1s);
list_match = match.(r"xlsx$",folder_content)
data_filesO1s = folder_content[list_match.!=nothing]
# fileName = data_filesO1s[end][1:end-5];
α_noiseO1s = Dict();
α_ratioO1s = Dict();
for idx_file in 1:length(data_filesO1s)
    # idx_file = 3
    fileName = data_filesO1s[idx_file][1:end-5];
    ρO1s_bulk = 1.0e-3parse(Cdouble,match(r"SDS_[0-9]*mM",fileName).match[5:end-2]); # extract concentration in mM

    xf_data = XLSX.readxlsx(string(data_folderO1s,fileName,".xlsx"));
    xf_data_sheet_names = XLSX.sheetnames(xf_data);

    # TODO: find the data and the fitted cross sections!
    list_match = match.(r"Eph ?= ?",xf_data_sheet_names);
    list_match_fit = match.(r"[Ff]itt?ing(_| )?results?",xf_data_sheet_names);
    xf_raw_data_sheet_names = xf_data_sheet_names[list_match.!=nothing];
    xf_fit_data_sheet_names = xf_data_sheet_names[list_match_fit.!=nothing];
    Ndata = length(xf_raw_data_sheet_names);
    dictAllData = Dict();
    for xf_name in xf_raw_data_sheet_names
        local x = XLSX.getdata(xf_data[xf_name])[2:end,:]
        local col_sym = string.(XLSX.gettable(xf_data[xf_name])[2])
        local dataPairs = Array{Pair{String,Vector{Cdouble}}}(undef,length(col_sym));
        for j in 1:length(col_sym)
            if (typeof(x[1,j])<:AbstractString)
                dataPairs[j] = (col_sym[j] => parse.(Cdouble,x[:,j]))
            else
                dataPairs[j] = (col_sym[j] => convert(Array{Cdouble,1},x[:,j]))
            end
        end
        df = DataFrame(dataPairs);
        dictAllData[Symbol(string("hν_",match.(r"=[0-9]*eV",xf_name).match[2:end-2]))] = df # (df,names(df))
        # figure()
        # scatter(df.Wavelength,df.Raw_spectrum)
        # plot(df.Wavelength,df.Background)
        # plot(df.Wavelength,df.Curve1)
        # plot(df.Wavelength,df.Curve2)
        # # plot(df.Wavelength,df.Curve3)
        # ylim(0.0)
    end


    x_fit = XLSX.getdata(xf_data[xf_fit_data_sheet_names[1]]);
    # remove missing columns
    x_fit = x_fit[:,broadcast(~,(ismissing.(x_fit[1, :])))]
    df_fit = DataFrame([col for col in eachcol(x_fit[2:end, :])], Symbol.(x_fit[1, :]))
    # remove missing rows
    filter!(x -> any(!ismissing, x), df_fit)
    # df_fit = dropmissing(df_fit, disallowmissing=true)
    df_fit = df_fit |> @orderby(_[Symbol("Photon energy")]) |> @thenby(_[Symbol("Binding energy")]) |> DataFrame

    # list all the photon energies
    df_Eph = select(df_fit,Symbol("Photon energy")) |> @unique() |> @orderby(_[Symbol("Photon energy")]) |> DataFrame
    # put everything in nice boxes
    dictPeak = Dict();
    for i in df_Eph[!,Symbol("Photon energy")]
        local df = df_fit |> @filter(_[Symbol("Photon energy")]==i) |> DataFrame
        dictPeak[Symbol(string("hν_",i))] = df # (df,names(df))
    end



    # geometry setup
    λe0 = 2.0e-3;         # reference penetration depth in μm
    δr = 2.0e-3           # transition to vacuum layer thickness (let's set about 1 nm)
    k0 = 10;              # compute the measurement model over a distance of k0*λe0
    Nr = 201;             # number of discretization points in the radial dimension
    Nθ = 256;             # number of discretization points in the polar angle dimension
    Ny = 256;             # number of discretization points in the cylinder axis dimension
    L  = 50.0;            # height of the irradiated sample (the vertical extent of the beam is more like 20μm instead of 100μm)
    μ0 = 10.0;            # radius of the cylinder
    x0 = sqrt(2.0)*5000.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
    y0 = 0.0;
    z0 = 5000.0;

    # spacial discretization 
    r = collect(range(μ0-k0*λe0,μ0+δr,length=Nr));
    θ0 = atan(x0,z0)
    θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
    y = collect(range(-L/2.0,L/2.0,length=Ny));

    #
    α_al_noise = zeros(Cdouble,Ndata);
    α_al_ratio = zeros(Cdouble,Ndata);
    global idx = 1

    # for i in df_Eph[!,Symbol("Photon energy")]
    if false
        figure(figsize=[12, 10])
    end
    for i in 1:Ndata
        local plot_sym = Symbol(string("hν_",df_Eph[!,Symbol("Photon energy")][i]));
        local Be = dictAllData[plot_sym].Wavelength;
        local dKe = median(abs.(Be[2:end]-Be[1:end-1]))
        if (typeof(dictPeak[plot_sym][!,Symbol("Binding energy")][1])<:AbstractString)
            local Be1 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Binding energy")][1])
            local Be2 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Binding energy")][2])
        else
            local Be1 = dictPeak[plot_sym][!,Symbol("Binding energy")][1]
            local Be2 = dictPeak[plot_sym][!,Symbol("Binding energy")][2]
        end
        # if (typeof(dictPeak[plot_sym][!,Symbol("Peak shift")][1])<:AbstractString)
        #     Be1 = Be1 + parse(Cdouble,dictPeak[plot_sym][!,Symbol("Peak shift")][1])
        #     Be2 = Be2 + parse(Cdouble,dictPeak[plot_sym][!,Symbol("Peak shift")][2])
        # else
        #     Be1 = Be1 + dictPeak[plot_sym][!,Symbol("Peak shift")][1]
        #     Be2 = Be2 + dictPeak[plot_sym][!,Symbol("Peak shift")][2]
        # end
        if (typeof(dictPeak[plot_sym][!,Symbol("FWHM(G)")][1])<:AbstractString)
            local σe1 = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("FWHM(G)")][1])
            local σe2 = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("FWHM(G)")][2])
        else
            local σe1 = 0.5*1.0e-3dictPeak[plot_sym][!,Symbol("FWHM(G)")][1];
            local σe2 = 0.5*1.0e-3dictPeak[plot_sym][!,Symbol("FWHM(G)")][2];
        end
        if (typeof(dictPeak[plot_sym][!,Symbol("FWHM(L)")][1])<:AbstractString)
            local σe1_L = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("FWHM(L)")][1])
            local σe2_L = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("FWHM(L)")][2])
        else
            local σe1_L = 0.5*1.0e-3dictPeak[plot_sym][!,Symbol("FWHM(L)")][1];
            local σe2_L = 0.5*1.0e-3dictPeak[plot_sym][!,Symbol("FWHM(L)")][2];
        end
        if (typeof(dictPeak[plot_sym][!,Symbol("Area")][1])<:AbstractString)
            local Ae1 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Area")][1])
            local Ae2 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Area")][2])
        else
            local Ae1 = dictPeak[plot_sym][!,Symbol("Area")][1]
            local Ae2 = dictPeak[plot_sym][!,Symbol("Area")][2]
        end

        local σ_peak_1 = (1.0/sqrt(2π*σe1^2))*exp.(-0.5*((Be.-Be1)/σe1).^2)
        # σ_peak_1 = σ_peak_1 + (1.0/(π*σe1_L))./(1.0 .+ ((Be.-Be1)/σe1_L).^2)
        local σ_peak_2 = (1.0/sqrt(2π*σe2^2))*exp.(-0.5*((Be.-Be2)/σe2).^2)
        # σ_peak_2 = σ_peak_2 + (1.0/(π*σe2_L))./(1.0 .+ ((Be.-Be2)/σe2_L).^2)

        
        if false
            ax = subplot(2,2,i)
            # title(string("Eph = ",df_Eph[!,Symbol("Photon energy")][i]," [eV]"),fontsize=14)
            plot(Be,dictAllData[plot_sym].Raw_spectrum,label="Data")
            plot(Be,dictAllData[plot_sym].Background.+dKe*(Ae1*σ_peak_1+Ae2*σ_peak_2),label="fitted spectra peaks")
            plot(Be,dictAllData[plot_sym].Background.+dictAllData[plot_sym].Curve1.+dictAllData[plot_sym].Curve2,label="fitted spectra curves")
            ylim(0.0)
            xlabel("binding energy [eV]",fontsize=14); 
            ylabel("spectrum [count]",fontsize=14) 
            xticks(fontsize=14); yticks(fontsize=14); 
            legend(fontsize=14)
            ax.invert_xaxis()
            ax.text(0.1, 0.5, string("Eph = ",df_Eph[!,Symbol("Photon energy")][i]," [eV]"), transform=ax.transAxes,fontsize=14)
        end
        # TODO: compute the measurement model
        local λe = 1.0e-3dictPeak[plot_sym][!,:IMFP][1]
        σ_all = dictAllData[plot_sym].Curve1 # .+dictAllData[plot_sym].Curve2.+dictAllData[plot_sym].Curve3;
        σ_all = σ_all/(dKe*sum(σ_all));
        # compute the geomtry factor
        H_geom,H_rθy,Arn,Aθj,Ayk = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λe);
        if false
            plot(r,H_geom)
        end
        # estimate the alignment parameter
        Sbg = dictAllData[plot_sym].Background;
        S_noisy = dictAllData[plot_sym].Raw_spectrum.-dictAllData[plot_sym].Curve2;
        ρ = ρO1s_bulk*ones(Cdouble,Nr)
        α_al_noise[idx],_    = noiseAndParameterEstimation(σ_all,H_geom,S_noisy,Sbg,ρ)
        α_al_noise[idx]      = α_al_noise[idx]/(dictPeak[plot_sym][!,Symbol("Sigma")][1]*dictPeak[plot_sym][!,Symbol("Photon flux")][1])
        α_al_ratio[idx]      = dictPeak[plot_sym][!,Symbol("Alignment")][1]
        global idx = idx + 1;
    end
    tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
    # savefig(string(fileName,"_plot.pdf"))
    # savefig(string(fileName,"_plot.png"))
    α_noiseO1s[Symbol(fileName)] = α_al_noise
    α_ratioO1s[Symbol(fileName)] = α_al_ratio
end


if false
    figure()
    for i in 1:length(α_noiseO1s)
        local plot_sym = Symbol(data_filesO1s[i][1:end-5]);
        scatter(1.0e9mean(α_noiseO1s[plot_sym]),mean(α_ratioO1s[plot_sym])^2.5)
    end
    for i in 1:length(α_noiseC1s)
        local plot_sym = Symbol(data_filesC1s[i][1:end-5]);
        scatter(1.0e9mean(α_noiseC1s[plot_sym]),mean(α_ratioC1s[plot_sym])^2.5)
    end
    xscale("log")
    yscale("log")


    figure(figsize=[10, 5])
    ax = subplot(121)
    for i in 1:length(α_noiseO1s)
        local plot_sym = Symbol(data_filesO1s[i][1:end-5]);
        # scatter(1.0e9mean(α_noise[plot_sym]),mean(α_ratio[plot_sym]))
        scatter(1.0e9α_noiseO1s[plot_sym],α_ratioO1s[plot_sym],label=replace(data_filesO1s[1][1:end-5],"_"=>" ")) #.^2.5
    end
    xlim(0.1,10.0)
    ylim(0.1,10.0)
    xlabel("model estimation [x\$10^{9}\$]",fontsize=14); 
    ylabel("liq O1s/gas O1s",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax.yaxis.offsetText.set_size(14)
    ax.xaxis.offsetText.set_size(14)
    xscale("log")
    yscale("log")
    legend(fontsize=14)

    ax = subplot(122)
    for i in 1:length(α_noiseC1s)
        local plot_sym = Symbol(data_filesC1s[i][1:end-5]);
        scatter(1.0e9α_noiseC1s[plot_sym],α_ratioC1s[plot_sym],label=replace(data_filesC1s[1][1:end-5],"_"=>" ")) # .^2.5
    end
    xlim(0.1,10.0)
    ylim(0.1,10.0)
    xlabel("model estimation [x\$10^{9}\$]",fontsize=14); 
    ylabel("liq O1s/gas O1s",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax.yaxis.offsetText.set_size(14)
    ax.xaxis.offsetText.set_size(14)
    xscale("log")
    yscale("log")
    legend(fontsize=14)
    tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)

    # savefig(string("../data/TK/","ratio_vs_model.png"))
    # savefig(string("../data/TK/","ratio_vs_model.pdf"))
end
