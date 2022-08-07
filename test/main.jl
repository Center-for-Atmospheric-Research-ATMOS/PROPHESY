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

R_1j = A_peak[2:end]/A_peak[1];
σR_ij = Aj_bg[2:end]./(A_peak[1]^2) # needa check again the definition of the variance... it requires splitting the different peak, i.e. A_peak_1 = A_peak.*τt[:,1], A_peak_2 = A_peak.*τt[:,2] and A_peak_3 = A_peak.*τt[:,3]

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

function modelSensitivity(τ::Array{Cdouble},H::Array{Cdouble,2},ρ::Array{Cdouble,1})
    Nr = length(ρ);
    J_all = zeros(Cdouble,Nr,Nr);
    for j in 2:Ndata
        J_all[:,:] = J_all[:,:] + modelSensitivity_j(τ[1],τ[j],H[1,:],H[j,:],ρ);
    end
    J_all
end

ρ_gt = Array{Cdouble,1}(dictAllGeom[key_symbol_geom[1]][1][:ρ]);
Nr = length(ρ_gt);
H = zeros(Cdouble,Ndata,Nr);
for i in 1:Ndata
    H[i,:] = dictAllGeom[key_symbol_geom[i]][1][:H]
end
J_all = modelSensitivity(τ_al_noise,H,ρ_gt);
F_j_all  = svd(J_all);
S_approx = zeros(Cdouble,Nr,Nr);
S_approx[1,1] = F_j_all.S[1];
J_approx = F_j_all.U*S_approx*F_j_all.Vt

figure();
imshow(J_all)
colorbar()

figure();
imshow(J_approx)
colorbar()


figure();
imshow(abs.(J_approx-J_all))
colorbar()

figure();
imshow(abs.(F_j_all.U[:,1:Ndata-1]*F_j_all.Vt[1:Ndata-1,:]))
colorbar()


# TODO: SA log cost function 


function acceptSampleRatio(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},
    y::Array{Cdouble,1},ΓI::Array{Cdouble,1},τH::Array{Cdouble,2},
    Dprior::Array{Cdouble,2})
    # if the posterior probability is larger for the proposed state ρ_prop than the current state ρ_cur, then accept the state, otherwise, reject it with probability p

    r_cp = 0.0
    # likelihood
    [global r_cp = r_cp + 0.5*((y[j]-((τH[j+1,:]*ρ_cur) /(τH[1,:]*ρ_cur)))^2) /ΓI[j] for j in 1:Ndata-1]
    [global r_cp = r_cp - 0.5*((y[j]-((τH[j+1,:]*ρ_prop)/(τH[1,:]*ρ_prop)))^2)/ΓI[j] for j in 1:Ndata-1]
    # smoothness
    r_cp = r_cp + 0.5ρ_cur'Dprior*ρ_cur - 0.5ρ_prop'Dprior*ρ_prop;

    ρ_new = ρ_cur;
    if (r_cp>=0.0)
        # unconditionally accept the new state
        ρ_new = ρ_prop
    else
        # accept the state with probability e^{r_cp}
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        else
            r_cp = 0.0
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
end

function samplePosterior(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},ΓI::Array{Cdouble,1},τH::Array{Cdouble,2},Dprior::Array{Cdouble,2};Ns::Int64=10000)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanismSmooth(ρ_all[i,:],Γsqrt)
        
        # accept or reject the sample
        ρ_all[i+1,:],_ = acceptSample(ρ_all[i,:],ρ_all[i+1,:],y,ΓI,τH,Dprior)
    end
    ρ_all
end


# TODO: use R_1j as input for sampling the posterior
# 