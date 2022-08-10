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
using XPSsampling


data_folder = "../data/cylinder_radius_10.0/eal_5_restricted_range/"

FLAG_0001 = false;
FLAG_0002 = false;
FLAG_0003 = true;
FLAG_0004 = false;

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

function modelSensitivity(τ::Array{Cdouble},H::Array{Cdouble,2},ρ::Array{Cdouble,1})
    Nr = length(ρ);
    J_all = zeros(Cdouble,Nr,Nr);
    for j in 2:Ndata
        J_all[:,:] = J_all[:,:] + modelSensitivity_j(τ[1],τ[j],H[1,:],H[j,:],ρ);
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
J_all = modelSensitivity(τ_al_noise,H,ρ_gt);

include("plotSensitivityMatrix.jl")

# TODO: SA log cost function 

idx_surf = findfirst(r.-10.0 .>-0.005);

function acceptSampleRatio(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},
    y::Array{Cdouble,1},ΓI::Array{Cdouble,1},τH::Array{Cdouble,2},
    Dprior::Array{Cdouble,2})
    # if the posterior probability is larger for the proposed state ρ_prop than the current state ρ_cur, then accept the state, otherwise, reject it with probability p

    r_cp = 0.0
    # likelihood
    for j in 1:length(y)
        r_cp = r_cp + 0.5*((y[j]-((τH[j+1,:]'*ρ_cur) /(τH[1,:]'*ρ_cur)))^2) /ΓI[j]
        r_cp = r_cp - 0.5*((y[j]-((τH[j+1,:]'*ρ_prop)/(τH[1,:]'*ρ_prop)))^2)/ΓI[j]
    end
    # smoothness
    r_cp = r_cp + 0.5ρ_cur'Dprior*ρ_cur - 0.5ρ_prop'Dprior*ρ_prop;

    ρ_new = ρ_cur;
    keepSample = false;
    if (r_cp>=0.0)
        # unconditionally accept the new state
        ρ_new = ρ_prop
        keepSample = true;
    else
        # accept the state with probability e^{r_cp}
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
            keepSample = true;
        else
            r_cp = 0.0
            keepSample = false;
        end
    end
    ρ_new,r_cp,keepSample # maybe it could return the computed values
end

Nρ = 20;
ρ_max = 2.0;
ρ_val = collect(LinRange(0.0,ρ_max,Nρ));
P_j = zeros(Cdouble,Nρ,Nρ);
δidx = 2.0;
for k in 1:Nρ
    P_j[:,k] = exp.(-abs.(collect(1:Nρ).-k)/δidx)
    P_j[k,k] = 0.0 # stating that the kernel makes sure the new value is different from the current one
    P_j[:,k] = P_j[:,k]./sum(P_j[:,k]);
end
figure(); imshow(P_j); colorbar()
cumsumP_j = cumsum(P_j,dims=1)
cumsumP_j[end,:] .= 1.0 


function newIdxValue(p_jk::Array{Cdouble,1})
    # normally, length(p_jk) = Nl+1
    Nρ = length(p_jk);
    v = rand();
    idx_new = -1;
    if (v<=p_jk[1])
        idx_new = 1;
    end
    if (idx_new==-1)
        for j in 1:Nρ-1
            if ((v>p_jk[j]) & (v<=p_jk[j+1]))
                idx_new = j+1;
            end
        end
    end
    idx_new
end

function newValue(p_jk::Array{Cdouble,1})
    ρ_val[newIdxValue(p_jk,Nρ)]
end


function kernelSingle(idx_curr::Array{Int64,1},cum_prod_j::Array{Cdouble,2})
    # draw index
    Nl = length(idx_curr);
    idx_elt = rand(1:Nl);
    # draw new value
    idx_next = idx_curr;
    idx_next[idx_elt] = newIdxValue(cum_prod_j[:,idx_curr[idx_elt]])
    # return 
    idx_next
end

# TODO: test kernel
testKer = ones(Int64,25,100);
for i in 2:100
    testKer[:,i] = kernelSingle(testKer[:,i-1],cumsumP_j)
end
figure(); imshow(testKer); colorbar()


function kernelSingle(ρ_curr::Array{Cdouble,1},σ::Cdouble;ρ_max::Cdouble=2.0)
    # draw index
    Nl = length(ρ_curr);
    idx_elt = rand(1:Nl)
    # draw new value
    ρ_next = ρ_curr;
    ρ_next[idx_elt] = ρ_curr[idx_elt] + σ*randn()
    # apply the constraint
    ρ_next[idx_elt] = max(ρ_next[idx_elt],-ρ_next[idx_elt])
    ρ_next[idx_elt] = min(ρ_next[idx_elt],ρ_max)
    ρ_next
end

function kernelTriple(ρ_curr::Array{Cdouble,1},σ::Cdouble;ρ_max::Cdouble=2.0)
    # draw index (the central index of the trhee)
    Nl = length(ρ_curr);
    idx_elt = rand(1:Nl)
    # modify the state
    ρ_next = ρ_curr;
    if ((idx_elt>1) && (idx_elt<Nl))
        # draw new values
        ρ_next[idx_elt] = ρ_curr[idx_elt] + σ*randn()
        ρ_next[idx_elt+1] = ρ_next[idx_elt+1] + 0.5σ*randn()
        ρ_next[idx_elt-1] = ρ_next[idx_elt-1] + 0.5σ*randn()
        # apply constraint
        ρ_next[idx_elt] = max(ρ_next[idx_elt],-ρ_next[idx_elt])
        ρ_next[idx_elt] = min(ρ_next[idx_elt],ρ_max)
        ρ_next[idx_elt+1] = max(ρ_next[idx_elt+1],-ρ_next[idx_elt+1])
        ρ_next[idx_elt+1] = min(ρ_next[idx_elt+1],ρ_max)
        ρ_next[idx_elt-1] = max(ρ_next[idx_elt-1],-ρ_next[idx_elt-1])
        ρ_next[idx_elt-1] = min(ρ_next[idx_elt-1],ρ_max)
    elseif (idx_elt==1)
        # draw new values
        ρ_next[idx_elt] = ρ_curr[idx_elt] + σ*randn()
        ρ_next[idx_elt+1] = ρ_next[idx_elt+1] + 0.5σ*randn()
        # apply constraint
        ρ_next[idx_elt] = max(ρ_next[idx_elt],-ρ_next[idx_elt])
        ρ_next[idx_elt] = min(ρ_next[idx_elt],ρ_max)
        ρ_next[idx_elt+1] = max(ρ_next[idx_elt+1],-ρ_next[idx_elt+1])
        ρ_next[idx_elt+1] = min(ρ_next[idx_elt+1],ρ_max)
    else
        # draw new values
        ρ_next[idx_elt] = ρ_curr[idx_elt] + σ*randn()
        ρ_next[idx_elt-1] = ρ_next[idx_elt-1] + 0.5σ*randn()
        # apply constraint
        ρ_next[idx_elt] = max(ρ_next[idx_elt],-ρ_next[idx_elt])
        ρ_next[idx_elt] = min(ρ_next[idx_elt],ρ_max)
        ρ_next[idx_elt-1] = max(ρ_next[idx_elt-1],-ρ_next[idx_elt-1])
        ρ_next[idx_elt-1] = min(ρ_next[idx_elt-1],ρ_max)
    end
    
    # return new state
    ρ_next
end

KERNEL_SMOOTH  = false
KERNEL_SINGLE  = false
KERNEL_SINGLE_IDX = true
KERNEL_TRIPLE  = false
function samplePosterior(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},ΓI::Array{Cdouble,1},τH::Array{Cdouble,2},Dprior::Array{Cdouble,2};Ns::Int64=10000)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start));
    r_cp = zeros(Cdouble,Ns);
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        if KERNEL_SMOOTH
            ρ_all[i+1,:] = XPSsampling.transmissionMechanismSmooth(ρ_all[i,:],Γsqrt)
            ρ_all[i+1,1:idx_surf-1] .= 1.0
            ρ_all[i+1,end] = 0.0
            ρ_all[ρ_all.<0.0] .= 0.0
            ρ_all[ρ_all.>3.0] .= 3.0
        end
        if KERNEL_SINGLE
            ρ_all[i+1,:] = [ρ_all[i,1:137]; kernelSingle(ρ_all[i,idx_surf:end-1],0.2); 0.0]
        end
        if KERNEL_TRIPLE
            ρ_all[i+1,:] = [ρ_all[i,1:137]; kernelTriple(ρ_all[i,idx_surf:end-1],0.2); 0.0]
        end

        # accept or reject the sample
        ρ_all[i+1,:],r_cp[i] = acceptSampleRatio(ρ_all[i,:],ρ_all[i+1,:],y,ΓI,τH,Dprior)
    end
    ρ_all,r_cp
end


function samplePosteriorIdx(idx_start::Array{Int64,1},ρ_val::Array{Cdouble,1},y::Array{Cdouble,1},ΓI::Array{Cdouble,1},τH::Array{Cdouble,2},Dprior::Array{Cdouble,2};Ns::Int64=10000,δidx::Cdouble=2.0)
    # init transistion probability 
    Nρ = length(ρ_val);
    P_j = zeros(Cdouble,Nρ,Nρ);
    for k in 1:Nρ
        P_j[:,k] = exp.(-abs.(collect(1:Nρ).-k)/δidx)
        P_j[k,k] = 0.0 # stating that the kernel makes sure the new value is different from the current one
        P_j[:,k] = P_j[:,k]./sum(P_j[:,k]);
    end
    cumsumP_j = cumsum(P_j,dims=1)
    cumsumP_j[end,:] .= 1.0

    # all samples
    ρ_all        = zeros(Cdouble,Ns+1,length(idx_start));
    idx_all      = zeros(Int64,Ns+1,length(idx_start));
    r_cp         = zeros(Cdouble,Ns);
    ρ_all[1,:]   = ρ_val[idx_start];
    idx_all[1,:] = idx_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        idx_all[i+1,:] = [idx_start[1:idx_surf-1]; kernelSingle(idx_all[i,idx_surf:end-1],cumsumP_j); 1];
        ρ_all[i+1,:] = ρ_val[idx_all[i+1,:]];

        # accept or reject the sample
        ρ_all[i+1,:],r_cp[i],keepSample = acceptSampleRatio(ρ_all[i,:],ρ_all[i+1,:],y,ΓI,τH,Dprior)
        if (!keepSample)
            idx_all[i+1,:] = idx_all[i,:];
        end
    end
    ρ_all,r_cp
end


# TODO: use R_1j as input for sampling the posterior
# data: R_1j_1, varaince: σR_ij_1

σw = 0.1*5.0e-4 # small compared with the amplitude of the state 
w = σw*ones(Cdouble,Nr); # not optimal because we know that the concentration varies more in the region near the surface rather than deep in the sample
Γsqrt = real(sqrt(corrCovariance(w;cor_len=10.0)));
p0 = 0.099 # shameful artifact
Ns      = 100000#0;
Ns_burn = 10000#0;
D2nd = diagm(Nr-2,Nr,1 => 2ones(Cdouble,Nr-2), 0 => -ones(Cdouble,Nr-2) ,2 => -ones(Cdouble,Nr-2));
Dprior = D2nd'*D2nd;

if KERNEL_SINGLE_IDX
    idx_start = [10*ones(Int64,137); round.(Int64,collect(LinRange(10,1,Nr-138+1)))]
    ρ_all,r_cp = samplePosteriorIdx(idx_start,ρ_val,R_1j_1,σR_ij_1,τ_al_noise.*H,1.0e3Dprior;Ns=Ns,δidx=2.0)
else
    ρ_start = ones(Cdouble,Nr); #  1.0ρ_gt
    ρ_start[138:end] = collect(LinRange(1,0,Nr-138+1));
    # ρ_start = ρ_gt
    ρ_all,r_cp = samplePosterior(ρ_start,Γsqrt,R_1j_1,σR_ij_1,τ_al_noise.*H,1.0e3Dprior;Ns=Ns); # 10000
end

figure(); plot(cumsum(r_cp));
minVal,idx_min = findmin(cumsum(r_cp))


μρ = dropdims(mean(ρ_all,dims=1),dims=1)

Γρ = cov(ρ_all)
figure(); ## plot(ρ_all[1:1000:end,:]')
plot(r,μρ,color=color_array[1],label = "mean")
plot(r,ρ_gt,color=color_array[3],label = "GT")
plot(r,ρ_all[idx_min+1,:],color=color_array[5],label = "last")
fill_between(r,μρ-sqrt.(diag(Γρ)),μρ+sqrt.(diag(Γρ)),alpha=0.5,color=color_array[1])

