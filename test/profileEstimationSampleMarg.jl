## can run this script with multithread: julia --threads NB_THREAD_YOU_WANT
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

# checking number of threads: launch Julia with the argument --threads 16 to run julia using 16 threads
println("you're running this script with ",Threads.nthreads()," threads") #WARNING: set ntasks to 1 when loading CSV files (it seems that multithreading is not safe with CSV.File)

SAVE_FIG = false;

# tags
FLAG_0001 = false             # selection of the profile (one must be true and the others false)
FLAG_0002 = false
FLAG_0003 = true
FLAG_0004 = false

MARG_UN   = true             # set to true to load the mean measurement operator as well as the covarainces

MODEL_ERROR_1 = false         # selection of the error level in the measurement model (1->0.5%, 2->1%, 3->2.5%)
MODEL_ERROR_2 = false
MODEL_ERROR_3 = true


SHORT_RANGE = false          # select either wide range of attenuation lengths (false) or a restricted range more similar to experimental setup (true)

MODEL_5   = false            # select the number of attenuation lengths probed
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
if (MODEL_20 & !SHORT_RANGE)
    model_folder = string(model_folder,"eal_20/")
    model_type = "20_datapoints_wide_range"
end
if (MODEL_5 & SHORT_RANGE)
    model_folder = string(model_folder,"eal_5_restricted_range/")
    model_type = "5_datapoints_short_range"
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

# truncation boundary
N0 = findfirst(r.-μ0.<=-d0);
N0_lowres = findfirst(r_lowres.-μ0.<=-d0);
N_lowres = N0_lowres-1;

# regularization (smoothness: applied as sparsity in the second order difference)
DN = D2nd(N_lowres+3);
D0 = DN[:,1];
D_tilde = DN[:,2:N_lowres+1];
Db = DN[:,N_lowres+2:N_lowres+3];

# correction of the regularization "data"
Δyd = -dropdims(sum(Db,dims=2),dims=2)*ρA_1[end];
δyd = -D0*ρA_1[1];
yd = Δyd+δyd;

#
# covariances: the crafed Bayesian models assumes that some covariance are known, i.e. measurement noise, smoothness and known values
#

# measurement noise covariance
ΓI = σnoise^2*diagm(ones(Cdouble,Ndata));

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

# parameters needed for the optimization
W_stop_lowres = ones(Cdouble,N_lowres);
N_max_iter = 20000#00; #NOTE: for very low noise levels, the convergence can be fairly slow, so, one might consider increasing the maximum number of iterations
r_n_tol=0.001;
r_y_tol=0.001; 
τ0 = 1.0e1 # 
x00 = 0.5ones(Cdouble,N_lowres); # since the concentration is normalized by the bulk concentration, the initial state is taken as uniform with value 1/2


##start loop model sampling here: start loading the model samples
N_model_sample = 100;
ρ_est_cp_block    = zeros(Cdouble,N_model_sample,Nr_lowres);
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

    # data correction
    local Δy = dropdims(sum(Hb,dims=2),dims=2)*ρA_1[end];
    local δy = H0*ρA_1[1];
    local y_tilde  = repData.-(Δy+δy)';

    local Γyb = σB^2*dropdims(sum([H0 Hb],dims=2),dims=2)*dropdims(sum([H0 Hb],dims=2),dims=2)';
    local Γy_tilde = ΓI + diagm(diag(Γyb))#  = Γyb+


    # smoosh together the several part of the model into augmented operators
    local Htrunc = [H_tilde; D_tilde];                                     # conditional to data and measurement model

    #
    # data inversion: estimation of the concentration profile using CP algorithm
    #


    # for each noise sample
    local Nsample = min(1,Nrep);
    # for i in 1:Nsample
        local i = 1 # just the first sample which is not corrupted by noise
        println(i,"/",Nsample)
        # augmented data
        # local Y = [y_tilde[i,:]; yd];

        # 
        local ρ_est,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_quad(x00,y_tilde[i,:],yd,Htrunc,Γy_tilde,Γd_lowres,W_stop_lowres;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)
        global ρ_est_cp_block[m,:] = [ρA_1[1]; ρ_est; ρA_1[end]*ones(Cdouble,Nr_lowres-N0_lowres)]
        println(N_last,"/",N_max_iter)
    # end

end

μρ_I_un = dropdims(mean(ρ_est_cp_block[1:min(N_model_sample,100),:],dims=1),dims=1);
Γρ_I_un = cov(ρ_est_cp_block[1:min(N_model_sample,100),:]);

figure(figsize=[10,6])
# ax1 = subplot(121)
l_cp_noise_mean, = plot(1000.0(μ0.-r_lowres),μρ_I_un,color="tab:blue")
l_cp_noise_cov   = fill_between(1000.0(μ0.-r_lowres),μρ_I_un-sqrt.(diag(Γρ_I_un)),μρ_I_un+sqrt.(diag(Γρ_I_un)),alpha=0.5,color="tab:blue")
l_gt,            = plot(1000.0(μ0.-r),ρA_1,color="tab:green")
# legend([(l_cp_post_est,l_cp_post_cov),(l_cp_noise_mean,l_cp_noise_cov),l_gt],["sampled posterior","est.+noise variability","GT"],fontsize=14,loc="lower right")
xlim(0.0,10.0)
ylim(0.0,1.5maximum(ρA_1))
xlabel("depth [nm]",fontsize=14)
xticks(fontsize=14)
ylabel("concentration [a.u.]",fontsize=14)
yticks(fontsize=14)
