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
data_folderC1s = "../data/TK/C1s/"

# find the data files
folder_content = readdir(data_folderC1s);
list_match = match.(r"xlsx$",folder_content)
data_filesC1s = folder_content[list_match.!=nothing]
# fileName = data_filesC1s[end][1:end-5];
α_noiseC1s = Dict();
α_ratioC1s = Dict();
for idx_file in 1:length(data_filesC1s)
    fileName = data_filesC1s[idx_file][1:end-5];
    ρC1s_bulk = 1.0e-3parse(Cdouble,match(r"SDS_[0-9]*mM",fileName).match[5:end-2]); # extract concentration in mM

    # match.(r"SDS_[0-9]*mM",data_filesC1s)

    xf_data = XLSX.readxlsx(string(data_folderC1s,fileName,".xlsx"));
    xf_data_sheet_names = XLSX.sheetnames(xf_data);

    # TODO: find the data and the fitted cross sections!
    list_match = match.(r"Eph ?= ?",xf_data_sheet_names);
    list_match_fit = match.(r"[Ff]itt?ing(_| )?results?",xf_data_sheet_names);
    xf_raw_data_sheet_names = xf_data_sheet_names[list_match.!=nothing];
    xf_fit_data_sheet_names = xf_data_sheet_names[list_match_fit.!=nothing];

    Ndata = length(xf_raw_data_sheet_names);
    # # sort sheets by increasing photon energy
    # idx_hν = findfirst(XLSX.getdata(xf_data[xf_data_sheet_names[1]])[1,:].=="hν");
    # hν_list = zeros(Cdouble,Ndata);
    # for i in 1:Ndata
    #     hν_list[i] = XLSX.getdata(xf_data[xf_data_sheet_names[i]])[2,idx_hν]
    # end
    # idx_perm_data = sortperm(hν_list)
    # xf_data_sheet_names = xf_data_sheet_names[idx_perm_data];

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
        # plot(df.Wavelength,df.Curve3)
        # ylim(0.0)
    end


    x_fit = XLSX.getdata(xf_data[xf_fit_data_sheet_names[1]]);
    # remove missing columns
    x_fit = x_fit[:,broadcast(~,(ismissing.(x_fit[1, :])))]
    df_fit = DataFrame([col for col in eachcol(x_fit[2:end, :])], Symbol.(x_fit[1, :]))
    # remove missing rows
    df_fit = dropmissing(df_fit, disallowmissing=true)
    df_fit = df_fit |> @orderby(_[Symbol("Photon energy")]) |> @thenby(_[Symbol("Binding energy")]) |> DataFrame


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
    idx = 1

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
            local Be3 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Binding energy")][3])
        else
            local Be1 = dictPeak[plot_sym][!,Symbol("Binding energy")][1]
            local Be2 = dictPeak[plot_sym][!,Symbol("Binding energy")][2]
            local Be3 = dictPeak[plot_sym][!,Symbol("Binding energy")][3]
        end
        if (typeof(dictPeak[plot_sym][!,Symbol("Peak shift")][1])<:AbstractString)
            Be1 = Be1 + parse(Cdouble,dictPeak[plot_sym][!,Symbol("Peak shift")][1])
            Be2 = Be2 + parse(Cdouble,dictPeak[plot_sym][!,Symbol("Peak shift")][2])
            Be3 = Be3 + parse(Cdouble,dictPeak[plot_sym][!,Symbol("Peak shift")][3])
        else
            Be1 = Be1 + dictPeak[plot_sym][!,Symbol("Peak shift")][1]
            Be2 = Be2 + dictPeak[plot_sym][!,Symbol("Peak shift")][2]
            Be3 = Be3 + dictPeak[plot_sym][!,Symbol("Peak shift")][3]
        end
        if (typeof(dictPeak[plot_sym][!,Symbol("FWHM(G)")][1])<:AbstractString)
            local σe1 = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("FWHM(G)")][1])
            local σe2 = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("FWHM(G)")][2])
            local σe3 = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("FWHM(G)")][3])
        else
            local σe1 = 0.5*1.0e-3dictPeak[plot_sym][!,Symbol("FWHM(G)")][1];
            local σe2 = 0.5*1.0e-3dictPeak[plot_sym][!,Symbol("FWHM(G)")][2];
            local σe3 = 0.5*1.0e-3dictPeak[plot_sym][!,Symbol("FWHM(G)")][3];
        end
        if (typeof(dictPeak[plot_sym][!,Symbol("FWHM(L)")][1])<:AbstractString)
            local σe1_L = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("FWHM(L)")][1])
            local σe2_L = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("FWHM(L)")][2])
            local σe3_L = 0.5*1.0e-3parse(Cdouble,dictPeak[plot_sym][!,Symbol("FWHM(L)")][3])
        else
            local σe1_L = 0.5*1.0e-3dictPeak[plot_sym][!,Symbol("FWHM(L)")][1];
            local σe2_L = 0.5*1.0e-3dictPeak[plot_sym][!,Symbol("FWHM(L)")][2];
            local σe3_L = 0.5*1.0e-3dictPeak[plot_sym][!,Symbol("FWHM(L)")][3];
        end
        if (typeof(dictPeak[plot_sym][!,Symbol("Area")][1])<:AbstractString)
            local Ae1 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Area")][1])
            local Ae2 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Area")][2])
            local Ae3 = parse(Cdouble,dictPeak[plot_sym][!,Symbol("Area")][3])
        else
            local Ae1 = dictPeak[plot_sym][!,Symbol("Area")][1]
            local Ae2 = dictPeak[plot_sym][!,Symbol("Area")][2]
            local Ae3 = dictPeak[plot_sym][!,Symbol("Area")][3]
        end

        local σ_peak_1 = (1.0/sqrt(2π*σe1^2))*exp.(-0.5*((Be.-Be1)/σe1).^2)
        # σ_peak_1 = σ_peak_1 + (1.0/(π*σe1_L))./(1.0 .+ ((Be.-Be1)/σe1_L).^2)
        local σ_peak_2 = (1.0/sqrt(2π*σe2^2))*exp.(-0.5*((Be.-Be2)/σe2).^2)
        # σ_peak_2 = σ_peak_2 + (1.0/(π*σe2_L))./(1.0 .+ ((Be.-Be2)/σe2_L).^2)
        local σ_peak_3 = (1.0/sqrt(2π*σe3^2))*exp.(-0.5*((Be.-Be3)/σe3).^2)
        # σ_peak_3 = σ_peak_3 + (1.0/(π*σe3_L))./(1.0 .+ ((Be.-Be3)/σe3_L).^2)

        
        if false
            ax = subplot(2,2,i)
            # title(string("Eph = ",df_Eph[!,Symbol("Photon energy")][i]," [eV]"),fontsize=14)
            plot(Be,dictAllData[plot_sym].Raw_spectrum,label="Data")
            plot(Be,dictAllData[plot_sym].Background.+dKe*(Ae1*σ_peak_1+Ae2*σ_peak_2+Ae3*σ_peak_3),label="fitted spectra peaks")
            plot(Be,dictAllData[plot_sym].Background.+dictAllData[plot_sym].Curve1.+dictAllData[plot_sym].Curve2.+dictAllData[plot_sym].Curve3,label="fitted spectra curves")
            # plot(Be,dictAllData[plot_sym].Background.+dictAllData[plot_sym].Curve1,label="fitted spectra curve1")
            # plot(Be,dictAllData[plot_sym].Background.+dictAllData[plot_sym].Curve2,label="fitted spectra curve2")
            # plot(Be,dictAllData[plot_sym].Background.+dictAllData[plot_sym].Curve3,label="fitted spectra curve3")
            # plot(Be,dictAllData[plot_sym].Background,)
            # plot(Be,Ae1*σ_peak_1+Ae2*σ_peak_2+Ae3*σ_peak_3)
            ylim(0.0)
            xlabel("binding energy [eV]",fontsize=14); 
            ylabel("spectrum [count]",fontsize=14) 
            xticks(fontsize=14); yticks(fontsize=14); 
            legend(fontsize=14)
            ax.invert_xaxis()
            ax.text(0.1, 0.5, string("Eph = ",df_Eph[!,Symbol("Photon energy")][i]," [eV]"), transform=ax.transAxes,fontsize=14)
        end
        # TODO: compute the measurement model
        local λe = 1.0e-3dictPeak[plot_sym][!,:IMFP][2]
        σ_all = dictAllData[plot_sym].Curve1.+dictAllData[plot_sym].Curve2.+dictAllData[plot_sym].Curve3;
        σ_all = σ_all/(dKe*sum(σ_all));
        # compute the geomtry factor
        H_geom,H_rθy,Arn,Aθj,Ayk = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λe);
        if false
            plot(r,H_geom)
        end
        # estimate the alignment parameter
        Sbg = dictAllData[plot_sym].Background;
        S_noisy = dictAllData[plot_sym].Raw_spectrum;
        ρ = ρC1s_bulk*ones(Cdouble,Nr)
        α_al_noise[idx],_    = noiseAndParameterEstimation(σ_all,H_geom,S_noisy,Sbg,ρ)
        α_al_noise[idx]      = α_al_noise[idx]/(dictPeak[plot_sym][!,Symbol("Sigma")][1]*dictPeak[plot_sym][!,Symbol("Photon flux")][1])
        α_al_ratio[idx]      = dictPeak[plot_sym][!,Symbol("Alignment")][1]
        idx = idx + 1;
    end
    tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
    # savefig(string(fileName,"_plot.pdf"))
    # savefig(string(fileName,"_plot.png"))
    α_noiseC1s[Symbol(fileName)] = α_al_noise
    α_ratioC1s[Symbol(fileName)] = α_al_ratio
end

figure()
for i in 1:length(α_noiseC1s)
    local plot_sym = Symbol(data_filesC1s[i][1:end-5]);
    scatter(1.0e9mean(α_noiseC1s[plot_sym]),mean(α_ratioC1s[plot_sym])^2.5)
end

figure()
for i in 1:length(α_noiseC1s)
    local plot_sym = Symbol(data_filesC1s[i][1:end-5]);
    scatter(1.0e9α_noiseC1s[plot_sym],α_ratioC1s[plot_sym].^2.5)
end

if false

    key_symbol = Symbol.(keys(dictAllData));
    key_symbol = key_symbol[idx_perm_data]



    xf_model = XLSX.readxlsx(string(data_folderC1s,"model.xlsx"));
    xf_model_sheet_names = XLSX.sheetnames(xf_model);
    idx_hν = findfirst(XLSX.getdata(xf_model[xf_model_sheet_names[1]])[1,:].=="hν");
    for i in 1:Ndata
        hν_list[i] = XLSX.getdata(xf_model[xf_model_sheet_names[i]])[2,idx_hν]
    end
    idx_perm_model = sortperm(hν_list)
    xf_model_sheet_names = xf_model_sheet_names[idx_perm_model];

    dictAllGeom = Dict();
    for xf_name in xf_model_sheet_names
        local df = DataFrame(XLSX.getdata(xf_model[xf_name])[2:end,:],:auto)
        rename!(df,Symbol.(XLSX.getdata(xf_model[xf_name])[1,:]));
        dictAllGeom[Symbol(xf_name)] = (eachcol(df),names(df))
    end
    key_symbol_geom = Symbol.(keys(dictAllGeom));
    key_symbol_geom = key_symbol_geom[idx_perm_model]

    # TODO: write project for simultaneously compute the peak fitting and the background removal from the model:
    #       I_i = Poisson(aσ(p,K_i)+Sbg_i) ≃ aσ(p,K_i)+Sbg_i +ε_i, ε_i ∼ N(0,I_i) and optimize for p (parameter of the function to be fitted)

    ##
    ## background removal
    ##
    Sbg_est = Dict()
    plot_sym = :Be # :Ke
    for j_slect in 1:Ndata
        # j_slect = 2
        Ny = length(dictAllData[key_symbol[j_slect]][1][:Ke]);
        D_2nd = diagm(Ny-2,Ny,1 => 2ones(Cdouble,Ny-2), 0 => -ones(Cdouble,Ny-2) ,2 => -ones(Cdouble,Ny-2)); # D2nd(Ny);
        # D_2nd = D_2nd./(dictAllData[key_symbol[j_slect]][1][:Stot][2:end-1]/maximum(dictAllData[key_symbol[j_slect]][1][:Stot]))
        κb = 0.0; #0.01;
        λb = 1.0e6;
        z_baseline_1 = baseline_removal(dictAllData[key_symbol[j_slect]][1][:Snoisy],λb,κb,D_2nd);
        Sbg_est[key_symbol[j_slect]] = z_baseline_1
    end 
    plot_sym = :Ke # :Be # 
    include("plotBckRemoval.jl")


    ##
    ## cross section density estimation
    ##


    τm = [1.0/3.0; 1.0/3.0; 1.0/3.0]; #[0.85; 0.125; 1.0-0.85-0.125];  # 
    μm = [290.2; 292.0; 293.0]; # [290.3; 291.9; 293.5];
    σm = sqrt(2.0)*[0.45; 0.25; 0.6]; # [290.3; 291.9; 293.5]/500.0;
    τt = zeros(Cdouble,Ndata,3);
    μt = zeros(Cdouble,Ndata,3);
    σt = zeros(Cdouble,Ndata,3);

    dKe = dictAllData[key_symbol[Ndata]][1][:Ke][2]-dictAllData[key_symbol[Ndata]][1][:Ke][1];
    S_cs_dens = Dict();
    for j in 1:Ndata
        local symbol_h = key_symbol[j];
        local be = dictAllData[symbol_h][1][:Be];
        # local spectrum = dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg]; 
        # local spectrum = S_oi[symbol_h][1];
        local spectrum = dictAllData[symbol_h][1][:Snoisy]-Sbg_est[symbol_h]
        spectrum[spectrum.<0.0] .= 0.0
        # estimate the peaks centers and spreads
        # 
        if j<=3
            # τt[j,:],μt[j,:],σt[j,:] = EM_peaks(be,spectrum,τm,μm,σm,200)
            τt[j,:],μt[j,:],σt[j,:] = EM_peaks(be[(be.>287.5).*(be.<295.5)],spectrum[(be.>287.5).*(be.<295.5)],τm,μm,σm,200) # crop spectrum and be, but it's not much better
        else
            # τt[j,:],μt[j,:],σt[j,:] = EM_peaks(be,spectrum,τm,μm,σm,500)
            τt[j,:],μt[j,:],σt[j,:] = EM_peaks(be[(be.>287.5).*(be.<295.5)],spectrum[(be.>287.5).*(be.<295.5)],τm,μm,σm,500) # crop spectrum and be, but it's not much better
        end

        # estimation
        σ_est_1 = τt[j,1]*(1.0/sqrt(2.0π*σt[j,1]^2))*exp.(-(be.-μt[j,1]).^2/(2.0σt[j,1]^2));
        σ_est_2 = τt[j,2]*(1.0/sqrt(2.0π*σt[j,2]^2))*exp.(-(be.-μt[j,2]).^2/(2.0σt[j,2]^2));
        σ_est_3 = τt[j,3]*(1.0/sqrt(2.0π*σt[j,3]^2))*exp.(-(be.-μt[j,3]).^2/(2.0σt[j,3]^2));

        # push the density to dictionary
        S_cs_dens[symbol_h] = σ_est_1+σ_est_2+σ_est_3
    end


    μBe = [290.2; 292.0; 293.0]
    σ_be = sqrt(2.0)*[0.45; 0.25; 0.6];


    plot_sym = :Be #  :Ke # 
    include("plotCrossSection.jl")


    ##
    ## noise removal
    ##

    S_oi = Dict();
    for (ν_sym,λ_sym) in zip(key_symbol,key_symbol_geom)
        # SVD of the measurement model (geometry factor and cross section density)
        F_λ  = svd(S_cs_dens[ν_sym]*dictAllGeom[λ_sym][1][:H]');
        s_th = 2
        # noise projection
        UF_ν  = F_λ.U[:,s_th:end]'*(dictAllData[ν_sym][1][:Snoisy]-Sbg_est[ν_sym]);
        # noise estimation
        noise_ν = F_λ.U[:,s_th:end]*UF_ν;
        # substract the noise from the signal
        σ_ν  = dictAllData[ν_sym][1][:Snoisy]-noise_ν;
        σ_ν[σ_ν.<=0.0] .= 0.0 # just make sure the signal is positive
        # keep track of the operation (spectrum of interest, non-noisy spectrum, noise)
        S_oi[ν_sym] = (σ_ν-Sbg_est[ν_sym],σ_ν,noise_ν)
    end


    plot_sym = :Ke #   :Be #  
    include("plotNoiseRemoval.jl") 

    ##
    ## compute the peak area form 
    ##

    A_peak = zeros(Cdouble,Ndata)
    Aj_bg = zeros(Cdouble,Ndata)
    for i in 1:Ndata
        A_peak[i] = dKe*sum(S_oi[key_symbol[i]][1])
        Aj_bg[i]  = dKe*sum(S_oi[key_symbol[i]][2])
    end

    A_peak_1 = A_peak.*τt[:,1];
    A_peak_2 = A_peak.*τt[:,2];
    A_peak_3 = A_peak.*τt[:,3];

    R_1j_1  = A_peak_1[2:end]/A_peak_1[1];
    σR_ij_1 = (Aj_bg[2:end]+A_peak[2:end])./(A_peak_1[1]^2) 

    R_1j_2  = A_peak_2[2:end]/A_peak_2[1];
    σR_ij_2 = (Aj_bg[2:end]+A_peak[2:end])./(A_peak_2[1]^2) 

    R_1j_3  = A_peak_3[2:end]/A_peak_3[1];
    σR_ij_3 = (Aj_bg[2:end]+A_peak[2:end])./(A_peak_3[1]^2) 

    SNR1 = R_1j_1./sqrt.(σR_ij_1);
    SNR2 = R_1j_2./sqrt.(σR_ij_2);
    SNR3 = R_1j_3./sqrt.(σR_ij_3);

    ##
    ## alignement factor estimation
    ##

    τ_al_noise    = zeros(Cdouble,Ndata);
    α_al_noise    = zeros(Cdouble,Ndata);
    for i in 1:Ndata
        local symbol_h = key_symbol[i]; 
        local simbol_λ = key_symbol_geom[i];
        local Nr = length(dictAllGeom[simbol_λ][1][:H])
        τ_al_noise[i],_    = noiseAndParameterEstimation(S_cs_dens[symbol_h],Array{Cdouble,1}(dictAllGeom[simbol_λ][1][:H]),Array{Cdouble,1}(dictAllData[symbol_h][1][:Snoisy]),Sbg_est[symbol_h],ones(Cdouble,Nr))
        α_al_noise[i] = τ_al_noise[i]/(dictAllData[symbol_h][1][:T][1]*dictAllData[symbol_h][1][:F][1]*XPSpack.σ_C1s_interp[dictAllData[symbol_h][1][:hν][1]])
    end


    ##
    ## sensitivity matrix
    ##

    function modelSensitivity_j(τ1::Cdouble,τj::Cdouble,H1::Array{Cdouble,1},Hj::Array{Cdouble,1},ρ::Array{Cdouble,1})
        ((τj/(τ1*H1'*ρ))^2)*(Hj - (Hj'*ρ/(H1'*ρ))*H1)*((Hj - (Hj'*ρ/(H1'*ρ))*H1)')
    end

    function modelSensitivity(τ::Array{Cdouble},H::Array{Cdouble,2},ρ::Array{Cdouble,1},σR_ij::Array{Cdouble,1})
        Nr = length(ρ);
        J_all = zeros(Cdouble,Nr,Nr);
        for j in 2:Ndata
            J_all[:,:] = J_all[:,:] + σR_ij[j-1]*modelSensitivity_j(τ[1],τ[j],H[1,:],H[j,:],ρ);
        end
        J_all
    end

    ρ_gt = Array{Cdouble,1}(dictAllGeom[key_symbol_geom[1]][1][:ρ]);
    r = Array{Cdouble,1}(dictAllGeom[key_symbol_geom[1]][1][:r]);
    μ0 = Array{Cdouble,1}(dictAllGeom[key_symbol_geom[1]][1][:radius])[1];
    Nr = length(ρ_gt);
    H = zeros(Cdouble,Ndata,Nr);
    for i in 1:Ndata
        H[i,:] = dictAllGeom[key_symbol_geom[i]][1][:H]
    end
    J_all = modelSensitivity(τ_al_noise,H,ρ_gt,σR_ij_1);

    include("plotSensitivityMatrix.jl")
end
