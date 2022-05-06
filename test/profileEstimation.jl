## load the packages used in the estimation
# plotting
using PyPlot
fm = PyPlot.matplotlib.font_manager.json_load("/home/matthew/.cache/matplotlib/fontlist-v310.json") # TODO: look for the path automatically
rc("text", usetex=true)
rc("figure",max_open_warning=50)
using myPlot

# data manipulation (loading, writing, etc)
using Printf
using CSV
using DataFrames

# scientific package from the official Julia repositories
using LinearAlgebra
using StatsBase

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack
using XPSinv

SAVE_FIG = false;

# tags
FLAG_0001 = true             # selection of the profile (one must be true and the others false)
FLAG_0002 = false
FLAG_0003 = false
FLAG_0004 = false

MARG_UN   = true             # set to true to load the mean measurement operator as well as the covarainces

MODEL_ERROR_1 = true         # selection of the error level in the measurement model (1->0.5%, 2->1%, 3->2.5%)
MODEL_ERROR_2 = false
MODEL_ERROR_3 = false


SHORT_RANGE = false          # select either wide range of attenuation lengths (false) or a restricted range more similar to experimental setup (true)

MODEL_5   = false            # select the number of attenuation lengths probed
MODEL_10  = false
MODEL_20  = true

FLAG_NOISE_1 = false          # selection of the noise level (one must be true and the others false)
FLAG_NOISE_2 = false
FLAG_NOISE_3 = false
FLAG_NOISE_4 = false
FLAG_NOISE_5 = true
FLAG_NOISE_6 = false

STRONG_PRIOR = false

filename_save = "truncated_model_hard_values"; 

# WARNING: σw depends on the strength of the prior (if the prior is weak, the posterior distribution will tend to the likelyhood, and so, the communication kernel based on the prior distribution will not be as efficient as in the strong a priori case)
σw = 5.0e-4 # small compared with the amplitude of the state 

if STRONG_PRIOR
    σd = 0.01;
    prior_strength = "strong_prior"
else
    σd = 0.1
    prior_strength = "weak_prior"
end
cor_len_lowres = 2.5 # this quantity should be set in terms of length instead of number of discretization points

# std of the known values (bulk and boundary)
σB = 0.01; 

# where to load the data from
model_folder = "../data/";

if (MODEL_5 & !SHORT_RANGE)
    model_folder = string(model_folder,"eal_5/")
    model_type = "5_datapoints_wide_range"
end
if (MODEL_10 & !SHORT_RANGE)
    model_folder = string(model_folder,"eal_10/")
    model_type = "10_datapoints_wide_range"
end
if (MODEL_20 & !SHORT_RANGE)
    model_folder = string(model_folder,"eal_20/")
    model_type = "20_datapoints_wide_range"
end
if (MODEL_5 & SHORT_RANGE)
    model_folder = string(model_folder,"eal_5_restricted_range/")
    model_type = "5_datapoints_short_range"
end
if (MODEL_10 & SHORT_RANGE)
    model_folder = string(model_folder,"eal_10_restricted_range/")
    model_type = "10_datapoints_short_range"
end



if FLAG_0001
    data_folder  = string(model_folder,"0001/")
    profile_flag = "00001"
end

if FLAG_0002
    data_folder  = string(model_folder,"0002/")
    profile_flag = "00002"
end

if FLAG_0003
    data_folder  = string(model_folder,"0003/")
    profile_flag = "00003"
end

if FLAG_0004
    data_folder  = string(model_folder,"0004/")
    profile_flag = "00004"
end



if FLAG_NOISE_1
    data_folder  = string(data_folder,"noise_level_0.001/")
    σnoise = 0.001;
    noise_level = "noise_0.001"
end
if FLAG_NOISE_2
    data_folder  = string(data_folder,"noise_level_0.005/")
    σnoise = 0.005;
    noise_level = "noise_0.005"
end
if FLAG_NOISE_3
    data_folder  = string(data_folder,"noise_level_0.01/")
    σnoise = 0.01;
    noise_level = "noise_0.01"
end
if FLAG_NOISE_4
    data_folder  = string(data_folder,"noise_level_0.05/")
    σnoise = 0.05;
    noise_level = "noise_0.05"
end
if FLAG_NOISE_5
    data_folder  = string(data_folder,"noise_level_0.1/")
    σnoise = 0.1;
    noise_level = "noise_0.1"
end
if FLAG_NOISE_6
    data_folder  = string(data_folder,"noise_level_0.5/")
    σnoise = 0.5;
    noise_level = "noise_0.5"
end

model_folder_lowres = string(model_folder,"lowres/");

if MODEL_ERROR_1
    model_folder_lowres_un = string(model_folder_lowres,"error_model_0.005_percent/")
    model_error = "model_error_0.5"
end
if MODEL_ERROR_2
    model_folder_lowres_un = string(model_folder_lowres,"error_model_0.01_percent/")
    model_error = "model_error_1.0"
end
if MODEL_ERROR_3
    model_folder_lowres_un = string(model_folder_lowres,"error_model_0.025_percent/")
    model_error = "model_error_2.5"
end

filename_save = string(profile_flag,"_",filename_save,"_",model_type,"_",model_error,"_",noise_level,"_",prior_strength);

# load some data
dfRepData = CSV.File(string(data_folder,"repeated_data.csv");header=true) |> DataFrame;
repData   = Matrix{Cdouble}(dfRepData)[2:end,:];
λe        = Matrix{Cdouble}(dfRepData)[1,:];
Ndata     = length(λe);
Nrep      = size(repData,1);

# load some measurement model
dfr_lowres = CSV.File(string(model_folder_lowres,"radial_discretization_lowres.csv");header=true) |> DataFrame
r_lowres   = dropdims(Matrix{Cdouble}(dfr_lowres),dims=1)
Nr_lowres = length(r_lowres);
dfH_lowres = CSV.File(string(model_folder_lowres_un,"H_lowres.csv");header=true) |> DataFrame
H_lowres   = Matrix{Cdouble}(dfH_lowres);

if MARG_UN
    global dfμH = CSV.File(string(model_folder_lowres_un,"mean_H_lowres.csv");header=true) |> DataFrame
    global μH = Matrix{Cdouble}(Matrix{Cdouble}(dfμH)');
    global ΓH = Array{Cdouble}(undef,Nr_lowres,Nr_lowres,Ndata);
    for i in 1:Ndata
        local df = CSV.File(string(model_folder_lowres_un,"cov/cov_H_lowres_",i,".csv");header=false) |> DataFrame
        ΓH[:,:,i] = Matrix{Cdouble}(df)
    end
end

# load the GT
dfr   = CSV.File(string(model_folder,"radial_discretization.csv");header=true) |> DataFrame
if FLAG_0001
    dfRho = CSV.File(string(model_folder,"/0001/concentration_profile.csv");header=true) |> DataFrame
end
if FLAG_0002
    dfRho = CSV.File(string(model_folder,"/0002/concentration_profile.csv");header=true) |> DataFrame
end
if FLAG_0003
    dfRho = CSV.File(string(model_folder,"/0003/concentration_profile.csv");header=true) |> DataFrame
end
if FLAG_0004
    dfRho = CSV.File(string(model_folder,"/0004/concentration_profile.csv");header=true) |> DataFrame
end

r = dropdims(Matrix{Cdouble}(dfr),dims=1);
μ0 = r[1]; # for now, the first radial discretization distance is exactly on the boundary of the cylinder
ρA_1 = reverse(dropdims(Matrix{Cdouble}(dfRho),dims=1));



# deeper than some distance, the signal is not likely to be disantangled
d0 = 5.0e-3 # 15.0e-3 # NOTE: this value should depend on the penetration depth
N0 = findfirst(r.-μ0.<=-d0);
N0_lowres = findfirst(r_lowres.-μ0.<=-d0);
N_lowres = N0_lowres-1;

# slice the model (3 terms: boundary, surface and bulk)
H0 = H_lowres[:,1];
H_tilde = H_lowres[:,2:N0_lowres];
Hb = H_lowres[:,N0_lowres+1:end];
if MARG_UN
    global μH0 = μH[:,1];
    global μH_tilde = μH[:,2:N0_lowres];
    global μHb = μH[:,N0_lowres+1:end];
end

# data correction
Δy = dropdims(sum(Hb,dims=2),dims=2)*ρA_1[end];
δy = H0*ρA_1[1];
y_tilde  = repData.-(Δy+δy)';

if MARG_UN
    global Δyμ = dropdims(sum(μHb,dims=2),dims=2)*ρA_1[end];
    global δyμ = μH0*ρA_1[1];
    global y_tildeμ = repData.-(Δyμ+δyμ)';
end


# regularization (smoothness: applied as sparsity in the second order difference)
DN = D2nd(N_lowres+3);
D0 = DN[:,1];
D_tilde = DN[:,2:N_lowres+1];
Db = DN[:,N_lowres+2:N_lowres+3];

# correction of the regularization "data"
Δyd = -dropdims(sum(Db,dims=2),dims=2)*ρA_1[end];
δyd = -D0*ρA_1[1];
yd = Δyd+δyd;


# smoosh together the several part of the model into augmented operators
Htrunc = [H_tilde; D_tilde];                                     # conditional to data and measurement model
Htrunc_un = [μH_tilde; D_tilde; diagm(ones(Cdouble,N_lowres))];  # conditional to data

#
# covariances: the crafed Bayesian models assumes that some covariance are known, i.e. measurement noise, smoothness, known values and measurement operator (the last only in the marginalized case)
#

# measurement noise covariance
ΓI = σnoise^2*diagm(ones(Cdouble,Ndata));

# covariance of the measurement model not in the surface layers
if MARG_UN
    global ΓH0b = Array{Cdouble,3}(undef,Nr_lowres-N_lowres,Nr_lowres-N_lowres,Ndata);
    for i in 1:Ndata
        ΓH0b[1,1,i] = ΓH[1,1,i]
        ΓH0b[1,2:end,i] = ΓH[1,N0_lowres+1:end,i]
        ΓH0b[2:end,1,i] = ΓH[N0_lowres+1:end,1,i]
        ΓH0b[2:end,2:end,i] = ΓH[N0_lowres+1:end,N0_lowres+1:end,i]
    end
    # accumulated measurement model error (not in the surface layers)
    global σεH = sqrt.(dropdims(sum(ΓH0b,dims=(1,2)),dims=(1,2)));
end


# covariance matrix for the a priori distribution (second order difference)
Γprior_lowres = zeros(Cdouble,Nr_lowres,Nr_lowres);
for i in 1:Nr_lowres
    Γprior_lowres[i,i] =  1.0;
    for j in i+1:Nr_lowres
        Γprior_lowres[i,j] = Γprior_lowres[i,i]*exp(-(i-j)^2/(0.5*cor_len_lowres^2));
        Γprior_lowres[j,i] = Γprior_lowres[i,j];
    end
end
Γd_lowres = (N_lowres/Ndata)*(σd^2)*Γprior_lowres[2:N_lowres+2,2:N_lowres+2];  # scale the a priori strength with the quantity of data, so that it is possible to compare the results
Γd_lowres_inv = inv(Γd_lowres);

# modifed model to noise ratio: the extra term comes from the truncation of the model, it accounts for the model uncertainty and the known values uncertainty
ΓHΓyinv           = zeros(Cdouble,N_lowres,N_lowres);
for k in 1:Ndata
    # global ΓHΓyinv = ΓHΓyinv + ΓH[2:N0_lowres,2:N0_lowres,k]/(ΓI[k,k]+(σB*σεH[k])^2+1.0e-2) # test: does not seem to stabilize the case very low noise
    global ΓHΓyinv = ΓHΓyinv + ΓH[2:N0_lowres,2:N0_lowres,k]/(ΓI[k,k]+(σB*σεH[k])^2)
end

#
# data inversion: estimation of the concentration profile using CP algorithm
#

W_stop_lowres = ones(Cdouble,N_lowres);
τ0 = 1.0e1 # 
x00 = 0.5ones(Cdouble,N_lowres); # since the concentration is normalized by the bulk concentration, the initial state is taken as uniform with value 1/2
N_max_iter = 20000#00; #NOTE: for very low noise levels, the convergence can be fairly slow, so, one might consider increasing the maximum number of iterations
r_n_tol=0.001;
r_y_tol=0.001;
r_y_tol_un=r_y_tol; 


# for each noise sample
Nsample = min(20,Nrep);
ρ_est_block       = zeros(Cdouble,Nsample,Nr_lowres);
ρ_est_cp_block    = zeros(Cdouble,Nsample,Nr_lowres);
ρ_est_cp_block_un = zeros(Cdouble,Nsample,Nr_lowres);
for i in 1:Nsample
    println(i,"/",Nsample)
    # augmented data
    local Y = [y_tilde[i,:]; yd];

    # naive inversion (no constraints)
    local ρ_trunc = inv(Htrunc'*Htrunc)*Htrunc'*Y;
    ρ_est_block[i,:] = [ρA_1[1]; ρ_trunc; ρA_1[end]*ones(Cdouble,Nr_lowres-N0_lowres)]

    # 
    local ρ_est,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_quad(x00,y_tilde[i,:],yd,Htrunc,ΓI,Γd_lowres,W_stop_lowres;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)
    ρ_est_cp_block[i,:] = [ρA_1[1]; ρ_est; ρA_1[end]*ones(Cdouble,Nr_lowres-N0_lowres)]
    println(N_last,"/",N_max_iter)
    
    # marginalization of uncertainty
    local ρ_est_un,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_quad_un(x00,y_tildeμ[i,:],yd,Htrunc_un,ΓI,Γd_lowres,ΓHΓyinv,W_stop_lowres;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol_un)
    ρ_est_cp_block_un[i,:] = [ρA_1[1]; ρ_est_un; ρA_1[end]*ones(Cdouble,Nr_lowres-N0_lowres)]
    println(N_last,"/",N_max_iter)
end


μρ = dropdims(mean(ρ_est_cp_block[1:Nsample,:],dims=1),dims=1);
Γρ = cov(ρ_est_cp_block[1:Nsample,:]);

μρ_un = dropdims(mean(ρ_est_cp_block_un[1:Nsample,:],dims=1),dims=1);
Γρ_un = cov(ρ_est_cp_block_un[1:Nsample,:]);


#
# posterior covariance estimation (get an idea of how good the estimation can be)
#

Nsample = 1
# σw = 1.0e-5 # 0.5*1.0e-3 # 0.001; # small compared with the amplitude of the state 
w = σw*ones(Cdouble,N_lowres); # not optimal because we know that the concentration varies more in the region near the surface rather than deep in the sample
Γsqrt = real(sqrt(corrCovariance(w;cor_len=10.0)));
p0 = 0.099 # shameful artifact
Ns      = 1000000;
Ns_burn = 100000;

ΓIinv = zeros(Cdouble,Ndata,Ndata);
[ΓIinv[k,k] = 1.0/(ΓI[k,k]+(σB*σεH[k])^2) for k in 1:Ndata]

μρ_I = zeros(Cdouble,N_lowres,Nsample);
Γρ_I = zeros(Cdouble,N_lowres,N_lowres,Nsample);
μρ_I_un = zeros(Cdouble,N_lowres,Nsample);
Γρ_I_un = zeros(Cdouble,N_lowres,N_lowres,Nsample);
deltaU = zeros(Cdouble,Ns);
deltaUun = zeros(Cdouble,Ns);
ρ_all = zeros(Cdouble,Ns+1,N_lowres);
ρ_all_un = zeros(Cdouble,Ns+1,N_lowres);
for i in 1:Nsample
    println(i,"/",Nsample)
    # conditional to data and model
    ρ_all[:], deltaU[:] = samplePosterior(ρ_est_cp_block[i,2:N0_lowres],Γsqrt,p0*ones(Cdouble,Ns),y_tilde[i,:],yd,ΓIinv,Γd_lowres_inv,H_tilde,D_tilde;Ns=Ns);
    
    # error marginalization
    ρ_all_un[:], deltaUun[:] = samplePosteriorMargin(ρ_est_cp_block_un[i,2:N0_lowres],Γsqrt,p0*ones(Cdouble,Ns),y_tildeμ[i,:],yd,ΓIinv,Γd_lowres_inv,μH_tilde,D_tilde,ΓHΓyinv;Ns=Ns);

    # compute a covariance matrix from the samples 
    μρ_I[:,i] = dropdims(mean(ρ_all[Ns_burn:Ns,:],dims=1),dims=1);
    Γρ_I[:,:,i] = cov(ρ_all[Ns_burn:Ns,:]);
    μρ_I_un[:,i] = dropdims(mean(ρ_all_un[Ns_burn:Ns,:],dims=1),dims=1);
    Γρ_I_un[:,:,i] = cov(ρ_all_un[Ns_burn:Ns,:]);
end


μμρ_I = dropdims(mean(μρ_I,dims=2),dims=2);
μΓρ_I = dropdims(mean(Γρ_I,dims=3),dims=3);

μμρ_I_un = dropdims(mean(μρ_I_un,dims=2),dims=2);
μΓρ_I_un = dropdims(mean(Γρ_I_un,dims=3),dims=3);


figure(); plot(cumsum(-deltaU)); plot(cumsum(-deltaUun))
if SAVE_FIG
    savefig(string(filename_save,"_energy_values.png"))
    savefig(string(filename_save,"_energy_values.pdf"))
end

# plot the estimation for both version (conditional to data and model, and conditional to data only) showing the covariance of the posterior 
# and the variability due to the noise in the data. For each profile, plot in one figure different level of noise and different level of model uncertainty (2 of each)

figure(figsize=[10,6])
ax1 = subplot(121)
l_cp_post_est,   = plot(1000.0(μ0.-r_lowres[2:N0_lowres]),μμρ_I,color="tab:red") # μρ
l_cp_post_cov    = fill_between(1000.0(μ0.-r_lowres[2:N0_lowres]),μμρ_I-sqrt.(diag(μΓρ_I)),μμρ_I+sqrt.(diag(μΓρ_I)),alpha=0.5,color="tab:red")
l_cp_noise_mean, = plot(1000.0(μ0.-r_lowres),μρ,color="tab:blue")
l_cp_noise_cov   = fill_between(1000.0(μ0.-r_lowres),μρ-sqrt.(diag(Γρ)),μρ+sqrt.(diag(Γρ)),alpha=0.5,color="tab:blue")
l_gt,            = plot(1000.0(μ0.-r),ρA_1,color="tab:green")
legend([(l_cp_post_est,l_cp_post_cov),(l_cp_noise_mean,l_cp_noise_cov),l_gt],["sampled posterior","est.+noise variability","GT"],fontsize=14,loc="lower right")
xlim(0.0,10.0)
ylim(0.0,1.5maximum(ρA_1))
xlabel("depth [nm]",fontsize=14)
xticks(fontsize=14)
ylabel("concentration [a.u.]",fontsize=14)
yticks(fontsize=14)

ax2 = subplot(122)
l_cp_post_est_un,   = plot(1000.0(μ0.-r_lowres[2:N0_lowres]),μμρ_I_un,color="tab:red")
l_cp_post_cov_un    = fill_between(1000.0(μ0.-r_lowres[2:N0_lowres]),μμρ_I_un-sqrt.(diag(μΓρ_I_un)),μμρ_I_un+sqrt.(diag(μΓρ_I_un)),alpha=0.5,color="tab:red")
l_cp_noise_mean_un, = plot(1000.0(μ0.-r_lowres),μρ_un,color="tab:blue",label="mean value marginal")
l_cp_noise_cov_un   = fill_between(1000.0(μ0.-r_lowres),μρ_un-sqrt.(diag(Γρ_un)),μρ_un+sqrt.(diag(Γρ_un)),alpha=0.5,color="tab:blue",label="uncertainty marginal")
l_gt_un,            = plot(1000.0(μ0.-r),ρA_1,color="tab:green",label="GT")
legend([(l_cp_post_est_un,l_cp_post_cov_un),(l_cp_noise_mean_un,l_cp_noise_cov_un),l_gt_un],["sampled posterior","est.+noise variability","GT"],fontsize=14,loc="lower right")
xlim(0.0,10.0)
ylim(0.0,1.5maximum(ρA_1))
xlabel("depth [nm]",fontsize=14)
xticks(fontsize=14)
ylabel("concentration [a.u.]",fontsize=14)
yticks(fontsize=14)

tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)
ax1.annotate("P\$(\\rho|H,I)\$", xy=(3, 1),  xycoords="data", xytext=(0.6, 0.85), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("P\$(\\rho|I)\$",   xy=(3, 1),  xycoords="data", xytext=(0.6, 0.85), textcoords="axes fraction", color="black",fontsize=14)

if SAVE_FIG
    println("saving: ",filename_save)
    savefig(string(filename_save,".png"))
    savefig(string(filename_save,".pdf"))
end