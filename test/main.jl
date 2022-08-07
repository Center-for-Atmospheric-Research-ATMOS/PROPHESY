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

# TODO: plot the different stages (select 4 spectra) 
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

figure(figsize=[12, 10]); 
plot_sym = :Ke # :Be # 
for (i_plot,i_data) in zip(1:4,Int64.(round.(collect(LinRange(1,Ndata,4)))))
    ax1 = subplot(2,2,i_plot)
    symbol_h = key_symbol[i_data]
    title(string("Eph =", Int64(round(dictAllData[symbol_h][1][:hν][1]))," [eV]"),fontsize=14)
    plot(dictAllData[symbol_h][1][plot_sym],Sbg_est[symbol_h],label="estimated background"); 
    plot(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1].Sbg,label="background gt"); 
    plot(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1].SpectrumA_1,label="SOI gt"); 
    scatter(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1].Snoisy-Sbg_est[symbol_h],label="estimated SOI"); 
    xlabel("kinetic energy [eV]",fontsize=14); ylabel("spectrum [count]",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    legend(fontsize=14)
    if (plot_sym==:Be)
        xlabel("binding energy [eV]",fontsize=14); 
        ax = gca()
        ax.invert_xaxis();
    end
end
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)


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



color_array = ["tab:blue"; "tab:orange"; "tab:green"; "tab:red"; "tab:purple"; "tab:brown"; "tab:pink"; "tab:gray"; "tab:olive"; "tab:cyan"; "magenta"; "yellow"; "hotpink"; "darkmagenta"; "chartreuse"; "deepskyblue"; "navy"; "darkcyan"; "crimson"; "firebrick"]; 
plot_sym = :Be #  :Ke # 
figure(figsize=[12, 10]); 
# for j in  1:Ndata # 1:5:Ndata
for (i_plot,i_data) in zip(1:4,Int64.(round.(collect(LinRange(1,Ndata,4)))))
    local μKe0=50.0
    local μKe1=1200.0
    local symbol_h = key_symbol[i_data]; 
    local be = dictAllData[symbol_h][1][plot_sym] #dictAllData[symbol_h][1][:Be];
    
    local μKe = dictAllData[symbol_h][1][:μKe];
    # partial cross section (one for each chemical state)
    σ_peak_1 = (1.0/sqrt(2.0π*σ_be[1]^2))*exp.(-(be.-μBe[1]).^2/(2.0σ_be[1]^2));
    σ_peak_2 = (1.0/sqrt(2.0π*σ_be[2]^2))*exp.(-(be.-μBe[2]).^2/(2.0σ_be[2]^2));
    σ_peak_3 = (1.0/sqrt(2.0π*σ_be[3]^2))*exp.(-(be.-μBe[3]).^2/(2.0σ_be[3]^2));
    # quantity of chemical states
    p1 = 0.85 .+ (0.77-0.85)*(μKe[1].-μKe0)./(μKe1-μKe0);
    p2 = 0.125 .+ (0.12-0.125)*(μKe[1].-μKe0)./(μKe1-μKe0);
    p3 = 1.0-(p1+p2);

    # plotting
    ax1 = subplot(2,2,i_plot)
    title(string("Eph =", Int64(round(dictAllData[symbol_h][1][:hν][1]))," [eV]"),fontsize=14)
    plot(be,p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3,color=color_array[3],label="GT") # i_plot
    plot(be,S_cs_dens[symbol_h],color=color_array[1],label="estimation") # i_plot
    scatter(be,(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])/(dKe*sum(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])),color=color_array[2],label="data") # i_plot

    xlabel("kinetic energy [eV]",fontsize=14); 
    ylabel("spectrum [count]",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    legend(fontsize=14)

    if (plot_sym==:Be)
        xlabel("binding energy [eV]",fontsize=14); 
        ax1.invert_xaxis();
    end
end
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)



##
## noise removal
##

S_oi = Dict();
for (ν_sym,λ_sym) in zip(key_symbol,key_symbol_geom)
    # SVD of the measurement model (geometry factor and cross section density)
    F_λ  = svd(S_cs_dens[ν_sym]*dictAllGeom[λ_sym][1][:H]');
    
    # figure(); imshow(F_λ.U); colorbar()
    # figure(); plot(F_λ.U[:,1])
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

    figure()
    
end


plot_sym = :Ke #   :Be #   
figure(figsize=[12, 10]); 
# for j in  1:Ndata # 1:5:Ndata
for (i_plot,i_data) in zip(1:4,Int64.(round.(collect(LinRange(1,Ndata,4)))))
    local symbol_h = key_symbol[i_data]; 
    local be = dictAllData[symbol_h][1][plot_sym] 
    local noise_ν = S_oi[symbol_h][3]

    # plotting
    ax1 = subplot(2,2,i_plot)
    title(string("Eph =", Int64(round(dictAllData[symbol_h][1][:hν][1]))," [eV]"),fontsize=14)
    scatter(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1][:Snoisy]-Sbg_est[symbol_h],color=color_array[2],label="spectrum-background")
    plot(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1][:Snoisy]-Sbg_est[symbol_h]-noise_ν,color=color_array[1],label="spectrum-background-noise")
    plot(dictAllData[symbol_h][1][plot_sym],noise_ν,color=color_array[4],label="estimated noise")
    plot(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1][:SpectrumA_1],color=color_array[3],label="GT")
    xlabel("kinetic energy [eV]",fontsize=14); 
    ylabel("spectrum [count]",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    legend(fontsize=14)

    if (plot_sym==:Be)
        xlabel("binding energy [eV]",fontsize=14); 
        ax1.invert_xaxis();
    end
end
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)




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



# TODO: SA log cost function 

function modelSensitivity_j(τ1::Cdouble,τj::Cdouble,H1::Array{Cdouble,1},Hj::Array{Cdouble,1},ρ::Array{Cdouble,1})
    ((τj/(τ1*H1'*ρ))^2)*(Hj - (Hj'*ρ/(H1'*ρ))*H1)*((Hj - (Hj'*ρ/(H1'*ρ))*H1)')
end



H1 = Array{Cdouble,1}(dictAllGeom[key_symbol_geom[1]][1][:H])
ρ_gt = H1 = Array{Cdouble,1}(dictAllGeom[key_symbol_geom[1]][1][:ρ]); # ones(Cdouble,Nr)
Nr = length(H1);
J_all = zeros(Cdouble,Nr,Nr);
for j in 2:Ndata
    Hj = Array{Cdouble,1}(dictAllGeom[key_symbol_geom[j]][1][:H])
    J_j = modelSensitivity_j(τ_al_noise[1],τ_al_noise[j],H1,Hj,ρ_gt)
    global J_all = J_all + J_j;
end

F_j_all  = svd(J_all);

figure();
imshow(J_all)
colorbar()


figure();
imshow(F_j_all.U)
colorbar()

# left and right singular vectors
figure()
semilogy(abs.(F_j_all.U[:,1:Ndata-1].*F_j_all.S[1:Ndata-1]'))
figure()
semilogy(abs.(F_j_all.Vt[1:Ndata-1,:].*F_j_all.S[1:Ndata-1])')


figure(); semilogy(F_j_all.S)

S_approx = zeros(Cdouble,Nr,Nr);
S_approx[1,1] = F_j_all.S[1];
J_approx = F_j_all.U*S_approx*F_j_all.Vt

figure();
imshow(abs.(J_approx-J_all))
colorbar()

figure();
imshow(abs.(F_j_all.U[:,1:Ndata-1]*F_j_all.Vt[1:Ndata-1,:]))
colorbar()



F_j_all.U[:,1]