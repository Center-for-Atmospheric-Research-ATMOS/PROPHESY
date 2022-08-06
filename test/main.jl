# load simulated data, and fit the peaks, estimate the alignment factor, remove the noise and estimate the sensitivity matrices
## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)
rc("figure",max_open_warning=50)
using myPlot

# data manipulation (loading, writing, etc)
using Printf
using XLSX # CSV does not deal with multiple sheets
using DataFrames

# scientific package from the official Julia repositories
using LinearAlgebra
using StatsBase
using Interpolations

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)


data_folder = "../data/cylinder_radius_10.0/eal_5_restricted_range/"

FLAG_0001 = false;
FLAG_0002 = false;
FLAG_0003 = false;
FLAG_0004 = true;

if FLAG_0001
    PROFILE_TAG = "0001/"
end
if FLAG_0002
    PROFILE_TAG = "0002/"
end
if FLAG_0003
    PROFILE_TAG = "0003/"
end
if FLAG_0004
    PROFILE_TAG = "0004/"
end

data_folder = string(data_folder,PROFILE_TAG);


xf_data = XLSX.readxlsx(string(data_folder,"data.xlsx"));
xf_data_sheet_names = XLSX.sheetnames(xf_data);
Ndata = length(xf_data_sheet_names);
# sort sheets by increasing photon energy
idx_hν = findfirst(XLSX.getdata(xf_data[xf_data_sheet_names[1]])[1,:].=="hν");
hν_list = zeros(Cdouble,Ndata);
for i in 1:Ndata
    hν_list[i] = XLSX.getdata(xf_data[xf_data_sheet_names[i]])[2,idx_hν]
end
idx_perm_data = sortperm(hν_list)
xf_data_sheet_names = xf_data_sheet_names[idx_perm_data];

dictAllData = Dict();
for xf_name in xf_data_sheet_names
    local df = DataFrame(Matrix{Cdouble}(XLSX.getdata(xf_data[xf_name])[2:end,:]),:auto)
    rename!(df,Symbol.(XLSX.getdata(xf_data[xf_name])[1,:]));
    dictAllData[Symbol(xf_name)] = (eachcol(df),names(df))
end
key_symbol = Symbol.(keys(dictAllData));
key_symbol = key_symbol[idx_perm_data]



xf_model = XLSX.readxlsx(string(data_folder,"model.xlsx"));
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

#TODO: plot the different stages (select 4 spectra) 

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
    λb = 1.0e8 # 1.0e10; #1.0e5;
    z_baseline_1 = baseline_removal(dictAllData[key_symbol[j_slect]][1][:Snoisy],λb,κb,D_2nd);
    Sbg_est[key_symbol[j_slect]] = z_baseline_1
    figure();
    scatter(dictAllData[key_symbol[j_slect]][1][plot_sym], dictAllData[key_symbol[j_slect]][1][:Snoisy])
    plot(dictAllData[key_symbol[j_slect]][1][plot_sym], dictAllData[key_symbol[j_slect]][1][:Snoisy]-z_baseline_1)
    plot(dictAllData[key_symbol[j_slect]][1][plot_sym], z_baseline_1)
    plot(dictAllData[key_symbol[j_slect]][1][plot_sym], dictAllData[key_symbol[j_slect]][1][:Sbg])
    plot(dictAllData[key_symbol[j_slect]][1][plot_sym], dictAllData[key_symbol[j_slect]][1][:SpectrumA_1])
    plot(dictAllData[key_symbol[j_slect]][1][plot_sym], dictAllData[key_symbol[j_slect]][1][:Stot])
    if (plot_sym==:Be)
        ax = gca()
        ax.invert_xaxis();
    end
end


##
## noise removal
##

S_oi = Dict();
for (ν_sym,λ_sym) in zip(key_symbol,key_symbol_geom)
    # SVD of the measurement model (geometry factor and cross section density)
    F_λ  = svd(dictAllData[ν_sym][1][:σ_cs_dens]*dictAllGeom[λ_sym][1][:H]');
    s_th = 2

    # noise projection
    # UF_ν  = F_λ.U[:,s_th:end]'*(dictAllData[ν_sym][1][:Snoisy]-dictAllData[ν_sym][1][:Sbg]); # with the true background
    UF_ν  = F_λ.U[:,s_th:end]'*(dictAllData[ν_sym][1][:Snoisy]-Sbg_est[ν_sym]);
    # noise estimation
    noise_ν = F_λ.U[:,s_th:end]*UF_ν;

    # substract the noise from the signal
    σ_ν  = dictAllData[ν_sym][1][:Snoisy]-noise_ν;
    σ_ν[σ_ν.<=0.0] .= 0.0 # just make sure the signal is positive

    # keep track of the operation (spectrum of interest, non-noisy spectrum, noise)
    S_oi[ν_sym] = (σ_ν-Sbg_est[ν_sym],σ_ν,noise_ν)

    x_symbol = :Ke;
    figure(); 
    plot(dictAllData[ν_sym][1][x_symbol],dictAllData[ν_sym][1][:Stot]) 
    scatter(dictAllData[ν_sym][1][x_symbol],σ_ν)

    figure(); 
    plot(dictAllData[ν_sym][1][x_symbol],dictAllData[ν_sym][1][:SpectrumA_1]) 
    scatter(dictAllData[ν_sym][1][x_symbol],σ_ν-Sbg_est[ν_sym])
    scatter(dictAllData[ν_sym][1][x_symbol],σ_ν-dictAllData[ν_sym][1][:Sbg])


    figure()
    scatter(collect(1:241),noise_ν,color="tab:blue")
    fill_between(collect(1:241),-sqrt.(σ_ν),sqrt.(σ_ν),alpha=0.5,color="tab:blue")
end


##
## cross section density estimation
##


τm = [0.85; 0.125; 1.0-0.85-0.125];  # [1.0/3.0; 1.0/3.0; 1.0/3.0];
μm = [290.2; 292.0; 293.0]; # [290.3; 291.9; 293.5];
σm = sqrt(2.0)*[0.45; 0.25; 0.6]; # [290.3; 291.9; 293.5]/500.0;
τt = zeros(Cdouble,Ndata,3);
μt = zeros(Cdouble,Ndata,3);
σt = zeros(Cdouble,Ndata,3);

dKe = dictAllData[key_symbol[Ndata]][1][:Ke][2]-dictAllData[key_symbol[Ndata]][1][:Ke][1];
for j in 1:Ndata
    local be = dictAllData[key_symbol[j]][1][:Be];
    # local spectrum = dictAllData[key_symbol[j]][1][:Snoisy]-dictAllData[key_symbol[j]][1][:Sbg]; 
    local spectrum = S_oi[key_symbol[j]][1];
    # estimate the peaks centers and spreads
    # τt[j,:],μt[j,:],σt[j,:] = EM_peaks(be[(be.>287.5).*(be.<295.5)],spectrum[(be.>287.5).*(be.<295.5)],τm,μm,σm,200) # crop spectrum and be, but it's not much better
    τt[j,:],μt[j,:],σt[j,:] = EM_peaks(be,spectrum,τm,μm,σm,200)
end


μBe = [290.2; 292.0; 293.0]
σ_be = sqrt(2.0)*[0.45; 0.25; 0.6];


figure(); 
scatter(collect(1:Ndata),abs.(μt[:,1].-μBe[1]))
scatter(collect(1:Ndata),abs.(μt[:,2].-μBe[2]))
scatter(collect(1:Ndata),abs.(μt[:,3].-μBe[3]))

figure(); 
scatter(collect(1:Ndata),σt[:,1])
scatter(collect(1:Ndata),σt[:,2])
scatter(collect(1:Ndata),σt[:,3])

figure(); 
scatter(collect(1:Ndata),τt[:,1])
scatter(collect(1:Ndata),τt[:,2])
scatter(collect(1:Ndata),τt[:,3])


color_array = ["tab:blue"; "tab:orange"; "tab:green"; "tab:red"; "tab:purple"; "tab:brown"; "tab:pink"; "tab:gray"; "tab:olive"; "tab:cyan"; "magenta"; "yellow"; "hotpink"; "darkmagenta"; "chartreuse"; "deepskyblue"; "navy"; "darkcyan"; "crimson"; "firebrick"]; 
for j in  1:Ndata # 1:5:Ndata
    local μKe0=50.0
    local μKe1=1200.0
    local symbol_h = key_symbol[j] ; # Symbol(string("hν_",Int64(round(hν[j]))));
    local be = dictAllData[symbol_h][1][:Be];
    local μKe = dictAllData[symbol_h][1][:μKe];
    # partial cross section (one for each chemical state)
    σ_peak_1 = (1.0/sqrt(2.0π*σ_be[1]^2))*exp.(-(be.-μBe[1]).^2/(2.0σ_be[1]^2));
    σ_peak_2 = (1.0/sqrt(2.0π*σ_be[2]^2))*exp.(-(be.-μBe[2]).^2/(2.0σ_be[2]^2));
    σ_peak_3 = (1.0/sqrt(2.0π*σ_be[3]^2))*exp.(-(be.-μBe[3]).^2/(2.0σ_be[3]^2));
    # quantity of chemical states
    p1 = 0.85 .+ (0.77-0.85)*(μKe[1].-μKe0)./(μKe1-μKe0);
    p2 = 0.125 .+ (0.12-0.125)*(μKe[1].-μKe0)./(μKe1-μKe0);
    p3 = 1.0-(p1+p2);

    # estimation
    σ_est_1 = τt[j,1]*(1.0/sqrt(2.0π*σt[j,1]^2))*exp.(-(be.-μt[j,1]).^2/(2.0σt[j,1]^2));
    σ_est_2 = τt[j,2]*(1.0/sqrt(2.0π*σt[j,2]^2))*exp.(-(be.-μt[j,2]).^2/(2.0σt[j,2]^2));
    σ_est_3 = τt[j,3]*(1.0/sqrt(2.0π*σt[j,3]^2))*exp.(-(be.-μt[j,3]).^2/(2.0σt[j,3]^2));

    println(dKe*sum(p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3))
    println(dKe*sum(σ_est_1+σ_est_2+σ_est_3),"\n")

    figure()
    plot(be,p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3,color=color_array[j])
    scatter(be,σ_est_1+σ_est_2+σ_est_3,color=color_array[j])

    plot(be,(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])/(dKe*sum(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])),color=color_array[j])
end




# TODO: alignment factor and  sensitivity matrix computation as well and SA log cost function 