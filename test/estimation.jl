## load the packages used in the estimation
# plotting
using PyPlot
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
# using XPSpack
using XPSinv
using XPSsampling

# checking number of threads: launch Julia with the argument --threads 16 to run julia using 16 threads
println("you're running this script with ",Threads.nthreads()," threads") #WARNING: set ntasks to 1 when loading CSV files (it seems that multithreading is not safe with CSV.File)

SAVE_FIG = false;

# tags
FLAG_0001 = false             # selection of the profile (one must be true and the others false)
FLAG_0002 = false
FLAG_0003 = false
FLAG_0004 = true

MODEL_ERROR_1 = false         # selection of the error level in the measurement model (1->0.5%, 2->1%, 3->2.5%)
MODEL_ERROR_2 = false
MODEL_ERROR_3 = true


SHORT_RANGE = false          # select either wide range of attenuation lengths (false) or a restricted range more similar to experimental setup (true)

MODEL_5   = true            # select the number of attenuation lengths probed
MODEL_10  = false
MODEL_20  = false

FLAG_NOISE_1 = false         # selection of the noise level (one must be true and the others false)
FLAG_NOISE_2 = false
FLAG_NOISE_3 = false
FLAG_NOISE_4 = false
FLAG_NOISE_5 = false
FLAG_NOISE_6 = true

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
dfRepData = CSV.File(string(data_folder,"repeated_data.csv");header=true,ntasks=1) |> DataFrame;
repData   = Matrix{Cdouble}(dfRepData)[2:end,:];
λe        = Matrix{Cdouble}(dfRepData)[1,:];
Ndata     = length(λe);
Nrep      = size(repData,1);

# load some measurement model
dfr_lowres = CSV.File(string(model_folder_lowres,"radial_discretization_lowres.csv");header=true,ntasks=1) |> DataFrame
r_lowres   = dropdims(Matrix{Cdouble}(dfr_lowres),dims=1)
Nr_lowres = length(r_lowres);
dfH_lowres = CSV.File(string(model_folder_lowres_un,"H_lowres.csv");header=true,ntasks=1) |> DataFrame
H_lowres   = Matrix{Cdouble}(dfH_lowres);

# load the GT
dfr   = CSV.File(string(model_folder,"radial_discretization.csv");header=true,ntasks=1) |> DataFrame
if FLAG_0001
    dfRho = CSV.File(string(model_folder,"/0001/concentration_profile.csv");header=true,ntasks=1) |> DataFrame
end
if FLAG_0002
    dfRho = CSV.File(string(model_folder,"/0002/concentration_profile.csv");header=true,ntasks=1) |> DataFrame
end
if FLAG_0003
    dfRho = CSV.File(string(model_folder,"/0003/concentration_profile.csv");header=true,ntasks=1) |> DataFrame
end
if FLAG_0004
    dfRho = CSV.File(string(model_folder,"/0004/concentration_profile.csv");header=true,ntasks=1) |> DataFrame
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
Hnot = [H0 Hb];
Hnot1 = sum(Hnot;dims=2);


# data correction
Δy = dropdims(sum(Hb,dims=2),dims=2)*ρA_1[end];
δy = H0*ρA_1[1];
y_tilde  = repData.-(Δy+δy)';



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

#
# covariances: the crafed Bayesian models assumes that some covariance are known, i.e. measurement noise, smoothness, known values and measurement operator (the last only in the marginalized case)
#

# measurement noise covariance
ΓI = σnoise^2*diagm(ones(Cdouble,Ndata));
ΓItrunc = ΓI + σB^2*Hnot1*Hnot1';
ΓIinv = inv(ΓItrunc);

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

w = σw*ones(Cdouble,N_lowres); # not optimal because we know that the concentration varies more in the region near the surface rather than deep in the sample
Γsqrt = real(sqrt(corrCovariance(w;cor_len=10.0)));
Ns      = 1000000;
Ns_burn =  100000;


# for each noise sample
Nsample = min(30,Nrep); #10
ρ_cp    = zeros(Cdouble,Nsample,Nr_lowres);
μρ_HI = zeros(Cdouble,N_lowres,Nsample);
Γρ_HI = zeros(Cdouble,N_lowres,N_lowres,Nsample);
deltaU = zeros(Cdouble,Ns,Nsample);

Threads.@threads  for i in 1:Nsample
    println(i,"/",Nsample)
    # argmax estimate
    local ρ_est,_,_,_,_,_,N_last = alg2_cp_quad(x00,y_tilde[i,:],yd,Htrunc,ΓItrunc,Γd_lowres,W_stop_lowres;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)
    global ρ_cp[i,:] = [ρA_1[1]; ρ_est; ρA_1[end]*ones(Cdouble,Nr_lowres-N0_lowres)];
    println(N_last,"/",N_max_iter)

    # posterior covariance estimation (conditional to data and model)
    println("posterior sampling: ",i)
    local ρ_all, deltaU[:,i] = samplePosterior(ρ_cp[i,2:N0_lowres],Γsqrt,y_tilde[i,:],yd,ΓIinv,Γd_lowres_inv,H_tilde,D_tilde;Ns=Ns);

    # compute a covariance matrix from the samples 
    global μρ_HI[:,i] = dropdims(mean(ρ_all[Ns_burn:Ns,:],dims=1),dims=1);
    global Γρ_HI[:,:,i] = cov(ρ_all[Ns_burn:Ns,:]);
end

# variability due to noise in data
# marginalization P(ρ|H) = ∫P(ρ|H,y)P(y)dy -> mean and covariance of ρ_HI = argmax P(ρ|H,y)
mean_ρ_H = dropdims(mean(ρ_cp[1:Nsample,:],dims=1),dims=1);
var_ρ_H  = cov(ρ_cp[1:Nsample,:]);

# marginalization P(ρ|H) = ∫P(ρ|H,y)P(y)dy -> mean and covariance of ρ∼P(ρ|H)
μρ_H = dropdims(mean(μρ_HI,dims=2),dims=2); # μ_{ρ|H}
Γρ_H = dropdims(mean(Γρ_HI,dims=3),dims=3); # Γ_{ρ|H}







##
## marginalization over the meas. op. space
##

i_sample       = min(2,Nsample);
N_model_sample = 30;
ρ_cp_HI        = zeros(Cdouble,N_model_sample,Nr_lowres);
deltaUh        = zeros(Cdouble,Ns,N_model_sample);
μρ_HI_sample   = zeros(Cdouble,N_lowres,N_model_sample);
Γρ_HI_sample   = zeros(Cdouble,N_lowres,N_lowres,N_model_sample);
# 
Threads.@threads for m in 1:min(N_model_sample,100) # this loop is just meant for sampling the model, so, each iteration is independent
    println(m,"/",min(N_model_sample,100))
    local num_sample;
    if (m<10)
        num_sample = string("00",m,"/")
    elseif ((m>=10) & (m<100))
        num_sample = string("0",m,"/")
    else
        num_sample = string(m,"/")
    end
    local dfH_lowres = CSV.File(string(model_folder_lowres_un,"/",num_sample,"H_lowres.csv");header=true,ntasks=1) |> DataFrame
    local H_lowres   = Matrix{Cdouble}(dfH_lowres);



    # slice the model (3 terms: boundary, surface and bulk)
    local H0 = H_lowres[:,1];
    local H_tilde = H_lowres[:,2:N0_lowres];
    local Hb = H_lowres[:,N0_lowres+1:end];

    # data and noise correction 
    local Δy      = dropdims(sum(Hb,dims=2),dims=2)*ρA_1[end];
    local δy      = H0*ρA_1[1];
    local y_tilde = repData.-(Δy+δy)';
    local Hnot    = [H0 Hb];
    local Hnot1   = sum(Hnot;dims=2);
    local ΓItrunc = ΓI + σB^2*Hnot1*Hnot1';
    local ΓIinv   = inv(ΓItrunc);


    # smoosh together the several part of the model into augmented operators
    local Htrunc = [H_tilde; D_tilde];                                     # conditional to data and measurement model

    # data inversion: estimation of the concentration profile using CP algorithm
    println(i_sample," data, model: ",m)
    local ρ_est,_,_,_,_,_,N_last = alg2_cp_quad(x00,y_tilde[i_sample,:],yd,Htrunc,ΓItrunc,Γd_lowres,W_stop_lowres;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)
    global ρ_cp_HI[m,:] = [ρA_1[1]; ρ_est; ρA_1[end]*ones(Cdouble,Nr_lowres-N0_lowres)]
    println(N_last,"/",N_max_iter," model ",m)

    # sample the posterior distribution for the fixed data y_tilde[i_sample,:]
    # posterior covariance estimation (conditional to data and model)
    println("posterior sampling: ",i_sample," model ",m)
    local ρ_all, deltaUh[:,m]  = samplePosterior(ρ_cp_HI[m,2:N0_lowres],Γsqrt,y_tilde[i_sample,:],yd,ΓIinv,Γd_lowres_inv,H_tilde,D_tilde;Ns=Ns);

    # compute a covariance matrix from the samples 
    global μρ_HI_sample[:,m]   = dropdims(mean(ρ_all[Ns_burn:Ns,:],dims=1),dims=1);
    global Γρ_HI_sample[:,:,m] = cov(ρ_all[Ns_burn:Ns,:]);
end

# compute the estimate variability due to the variation in the model
mean_ρ_y = dropdims(mean(ρ_cp_HI[1:min(N_model_sample,100),:],dims=1),dims=1);
var_ρ_y = cov(ρ_cp_HI[1:min(N_model_sample,100),:]);

# marginal mean and covariance
μρ_y = dropdims(mean(μρ_HI_sample,dims=2),dims=2); # μ_{ρ|y}


Γρ_y_mean = dropdims(mean(Γρ_HI_sample,dims=3),dims=3)
ΔΓρ_y1 = (1.0/N_model_sample^2)*μρ_y*μρ_y';
ΔΓρ_y2 = μρ_HI_sample[:,1]*μρ_HI_sample[:,1]';
[global ΔΓρ_y2 = ΔΓρ_y2 + μρ_HI_sample[:,i]*μρ_HI_sample[:,i]'; for i in 2:N_model_sample]
ΔΓρ_y2 = (1.0/N_model_sample)*ΔΓρ_y2;

Γρ_y = Γρ_y_mean - ΔΓρ_y1 + ΔΓρ_y2; # Γ_{ρ|y}

# figure(); imshow(Γρ_y); colorbar()
# figure(); imshow(Γρ_y_mean); colorbar()
# figure(); imshow(ΔΓρ_y1); colorbar()
# figure(); imshow(ΔΓρ_y2); colorbar()







figure(figsize=[10,6]); # plot(1000.0(μ0.-r_lowres),ρ_cp_HI');


ax1         = subplot(121) # variability of the estimates
l_ρ_cp,     = plot(1000.0(μ0.-r_lowres),ρ_cp[i_sample,:],color="tab:blue")
# conditional to y: variability
l_mean_ρ_y, = plot(1000.0(μ0.-r_lowres),mean_ρ_y,color="tab:red")
l_var_ρ_y   = fill_between(1000.0(μ0.-r_lowres),mean_ρ_y-sqrt.(diag(var_ρ_y)),mean_ρ_y+sqrt.(diag(var_ρ_y)),alpha=0.65,color="tab:red")
# conditional to H: variability
l_mean_ρ_H, = plot(1000.0(μ0.-r_lowres),mean_ρ_H,color="tab:orange")
l_var_ρ_H   = fill_between(1000.0(μ0.-r_lowres),mean_ρ_H-sqrt.(diag(var_ρ_H)),mean_ρ_H+sqrt.(diag(var_ρ_H)),alpha=0.25,color="tab:orange")
# ground truth
l_gt,       = plot(1000.0(μ0.-r),ρA_1,color="tab:green")
l_variability = [l_ρ_cp,(l_mean_ρ_y,l_var_ρ_y),(l_mean_ρ_H,l_var_ρ_H),l_gt];
l_names     = ["estimate \$\\hat{\\rho}|H,y\$","meas. op. var. \$\\hat{\\rho}|y\$","meas. noise var. \$\\hat{\\rho}|H\$","GT"]
legend(l_variability,l_names,fontsize=14,loc="upper right")
# legend([(l_cp_post_est,l_cp_post_cov),(l_cp_post_est_y,l_cp_post_cov_y),(l_cp_noise_mean,l_cp_noise_cov),(l_cp_var_mod,l_cp_var_mod_cov),l_gt],["sampled posterior P\$(\\rho|H,y)\$","sampled posterior P\$(\\rho|y)\$","\$E[\\hat{\\rho}|H]\\pm\\Gamma[\\hat{\\rho}|H]\$","\$E[\\hat{\\rho}|y]\\pm\\Gamma[\\hat{\\rho}|y]\$","GT"],fontsize=14,loc="upper right")

xlim(0.0,8.5)
ylim(-0.1,1.5maximum(ρA_1))
xlabel("depth [nm]",fontsize=14)
xticks(fontsize=14)
ylabel("concentration [a.u.]",fontsize=14)
yticks(fontsize=14)

ax2 = subplot(122) # distributions and their marginalizations
# conditional to y and H
l_μρ_HI, = plot(1000.0(μ0.-r_lowres)[2:N0_lowres],μρ_HI[:,i_sample],color="tab:blue")
l_Γρ_HI  = fill_between(1000.0(μ0.-r_lowres)[2:N0_lowres],μρ_HI[:,i_sample]-sqrt.(diag(Γρ_HI[:,:,i_sample])),μρ_HI[:,i_sample]+sqrt.(diag(Γρ_HI[:,:,i_sample])),alpha=0.25,color="tab:blue")
# conditional to y: mean and covariance of the distribution
l_μρ_y, = plot(1000.0(μ0.-r_lowres)[2:N0_lowres],μρ_y,color="tab:red")
l_Γρ_y  = fill_between(1000.0(μ0.-r_lowres)[2:N0_lowres],μρ_y-sqrt.(diag(Γρ_y)),μρ_y+sqrt.(diag(Γρ_y)),alpha=0.25,color="tab:red")
# conditional to H: mean and covariance of the distribution
l_μρ_H,   = plot(1000.0(μ0.-r_lowres[2:N0_lowres]),μρ_H,color="tab:orange")
l_Γρ_H    = fill_between(1000.0(μ0.-r_lowres[2:N0_lowres]),μρ_H-sqrt.(diag(Γρ_H)),μρ_H+sqrt.(diag(Γρ_H)),alpha=0.5,color="tab:orange")
# ground truth
l_gt,       = plot(1000.0(μ0.-r),ρA_1,color="tab:green")

l_post    = [(l_μρ_HI,l_Γρ_HI),(l_μρ_y,l_Γρ_y),(l_μρ_H,l_Γρ_H),l_gt];
l_names   = ["P\$(\\rho|H,y)\$","P\$(\\rho|y)\$","P\$(\\rho|H)\$","GT"]
legend(l_post,l_names,fontsize=14,loc="upper right")

xlim(0.0,8.5)
ylim(-0.1,1.5maximum(ρA_1))
xlabel("depth [nm]",fontsize=14)
xticks(fontsize=14)
ylabel("concentration [a.u.]",fontsize=14)
yticks(fontsize=14)


tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)
# ax1.annotate("P\$(\\rho|H,I)\$", xy=(3, 1),  xycoords="data", xytext=(0.6, 0.85), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)
# ax2.annotate("P\$(\\rho|I)\$",   xy=(3, 1),  xycoords="data", xytext=(0.6, 0.85), textcoords="axes fraction", color="black",fontsize=14)


# CSV.write: r_lowres, ρ_cp, mean_ρ_y, var_ρ_y, ρ_cp_HI, mean_ρ_H, var_ρ_H, r, ρA_1  # estimates, mean and variability
# CSV.write: μρ_HI_sample, Γρ_HI_sample, μρ_y, Γρ_y                                  # sampling the meas. op. to compute the conditional to meas. noise
# CSV.write: μρ_HI, Γρ_HI, μρ_H, Γρ_H                                                # sampling the meas. noise to compute the conditional to meas. op.

function displayCov(figNum::Int64,r::Array{Cdouble,1},Γ::Array{Cdouble,2};_sub::Int64=111,fontsize_ticks::Int64=14)
    min_Γ,max_Γ = extrema(Γ);
    fig,ax,pcm = imshowData(figNum,r,r,Γ;_norm=:NoNorm,_vmin=min_Γ,_vmax=max_Γ,_edgecolors="face",_shading="None", _sub=_sub);
    xticks(fontsize=fontsize_ticks)
    yticks(fontsize=fontsize_ticks)
    fig,ax,pcm
end

function setVerticalColorbar(fig::Figure,pcm::PyPlot.PyObject,x::Cdouble,y::Cdouble,dx::Cdouble,dy::Cdouble,slabel::String;fontsize_label::Int64=10,fontsize_ticks::Int64=10,color::String="white",_power_lim::Bool=true)
    rc("ytick",color=color)
    cax = fig.add_axes([x, y, dx, dy])
    cb = fig.colorbar(pcm, orientation="vertical", cax=cax)
    cb.set_label(slabel, fontsize=fontsize_label, color=color) # 
    cb.ax.yaxis.set_tick_params(color=color)
    cb.ax.tick_params(labelsize=fontsize_ticks)
    cb.outline.set_edgecolor(color)
    if _power_lim
        cb.formatter.set_powerlimits((-1,2))
        cb.ax.yaxis.offsetText.set_size(fontsize_ticks)
    end
    cb.update_ticks()
    rc("ytick",color="black") # back to black
    cb
end

fig = figure(128,figsize=[5,10]);
_,ax1,pcm1 = displayCov(128,1000.0(μ0.-r_lowres[2:N0_lowres]),sqrt(Γρ_HI[:,:,1]);_sub=311)
ylabel("depth [nm]",fontsize=14)
_,ax2,pcm2 = displayCov(128,1000.0(μ0.-r_lowres[2:N0_lowres]),sqrt(Γρ_y);_sub=312)
ylabel("depth [nm]",fontsize=14)
_,ax3,pcm3 = displayCov(128,1000.0(μ0.-r_lowres[2:N0_lowres]),sqrt(Γρ_H);_sub=313)
xlabel("depth [nm]",fontsize=14)
ylabel("depth [nm]",fontsize=14)

cb1 = setVerticalColorbar(fig,pcm1,0.17,0.815,0.04,0.15,"concentration [a.u.]";fontsize_label=14,fontsize_ticks=12) # \$^2\$ # 
cb2 = setVerticalColorbar(fig,pcm2,0.17,0.495,0.04,0.15,"concentration [a.u.]";fontsize_label=14,fontsize_ticks=12)
cb3 = setVerticalColorbar(fig,pcm3,0.17,0.178,0.04,0.15,"concentration [a.u.]";fontsize_label=14,fontsize_ticks=12)

tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)

ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)
ax3.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)


fig = figure(129,figsize=[5,10]);
_,ax1,pcm1 = displayCov(129,1000.0(μ0.-r_lowres[2:N0_lowres]),Γρ_HI[:,:,1];_sub=311)
ylabel("depth [nm]",fontsize=14)
_,ax2,pcm2 = displayCov(129,1000.0(μ0.-r_lowres[2:N0_lowres]),Γρ_y;_sub=312)
ylabel("depth [nm]",fontsize=14)
_,ax3,pcm3 = displayCov(129,1000.0(μ0.-r_lowres[2:N0_lowres]),Γρ_H;_sub=313)
xlabel("depth [nm]",fontsize=14)
ylabel("depth [nm]",fontsize=14)

cb1 = setVerticalColorbar(fig,pcm1,0.17,0.815,0.04,0.15,"concentration\$^2\$ [a.u.]";fontsize_label=14,fontsize_ticks=12) 
cb2 = setVerticalColorbar(fig,pcm2,0.17,0.495,0.04,0.15,"concentration\$^2\$ [a.u.]";fontsize_label=14,fontsize_ticks=12)
cb3 = setVerticalColorbar(fig,pcm3,0.17,0.178,0.04,0.15,"concentration\$^2\$ [a.u.]";fontsize_label=14,fontsize_ticks=12)

tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)

ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)
ax3.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)

