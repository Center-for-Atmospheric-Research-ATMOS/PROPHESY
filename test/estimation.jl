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

# inversion and sampling package
using XPSinv
using XPSsampling

# checking number of threads: launch Julia with the argument --threads 16 to run julia using 16 threads
println("you're running this script with ",Threads.nthreads()," threads") #WARNING: set ntasks to 1 when loading CSV files (it seems that multithreading is not safe with CSV.File)

# SAVE_FIG = false;

# tags
FLAG_0001 = false             # selection of the profile (one must be true and the others false)
FLAG_0002 = true
FLAG_0003 = false
FLAG_0004 = false

STRONG_PRIOR = false

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

n_noise_sample = 50; # number of noise samples for the noise marginalization
n_model_sample = 50; # number of noise samples for the noise marginalization


if FLAG_0001
    profile_flag = "0001"
end

if FLAG_0002
    profile_flag = "0002"
end

if FLAG_0003
    profile_flag = "0003"
end

if FLAG_0004
    profile_flag = "0004"
end


## loop
MODEL_TYPE_LIST = ("5_datapoints_short_range","5_datapoints_wide_range","20_datapoints_wide_range")
MODEL_FOLD_LIST = ("eal_5_restricted_range/","eal_5/","eal_20/")

MODEL_ERROR_LIST = ("model_error_0.5","model_error_1.0","model_error_2.5")
MODEL_ERROR_FOLD = ("error_model_0.005_percent/","error_model_0.01_percent/","error_model_0.025_percent/")
NOISE_LEVEL_LIST = ("noise_level_0.01/","noise_level_0.1/","noise_level_0.5/")
noise_sigma      = (0.01,0.1,0.5)

for (mod_str,mod_fol) in zip(MODEL_TYPE_LIST,MODEL_FOLD_LIST)
    for (mod_err,mod_err_fold) in zip(MODEL_ERROR_LIST,MODEL_ERROR_FOLD)
        for (dat_noise,σnoise) in zip(NOISE_LEVEL_LIST,noise_sigma)
            # paths: load and save
            global model_folder = string("../data/",mod_fol)                         # the common depth discretization should be in this folder
            global model_folder_lowres = string(model_folder,"lowres/")              # load the low resolution models from this folder
            global model_folder_lowres_un = string(model_folder_lowres,mod_err_fold) # 
            global data_folder  = string(model_folder,profile_flag,"/",dat_noise)    # load the data from this folder
            global save_folder  = string("results/",profile_flag,"/",mod_str,"/",mod_err,"/",dat_noise,prior_strength,"/") # save the results in this folder
            mkpath(save_folder)

            #
            # load some data
            #

            include("load_data.jl")

            #
            # truncation of the model and the data
            #

            # deeper than some distance, the signal is not likely to be disantangled
            global d0 = 5.0e-3 # 15.0e-3 # NOTE: this value should depend on the penetration depth
            global N0 = findfirst(r.-μ0.<=-d0);
            global N0_lowres = findfirst(r_lowres.-μ0.<=-d0);
            global N_lowres = N0_lowres-1;

            # slice the model (3 terms: boundary, surface and bulk)
            global H0 = H_lowres[:,1];
            global H_tilde = H_lowres[:,2:N0_lowres];
            global Hb = H_lowres[:,N0_lowres+1:end];
            global Hnot = [H0 Hb];
            global Hnot1 = sum(Hnot;dims=2);


            # data correction
            global Δy = dropdims(sum(Hb,dims=2),dims=2)*ρA_1[end];
            global δy = H0*ρA_1[1];
            global y_tilde  = repData.-(Δy+δy)';



            # regularization (smoothness: applied as sparsity in the second order difference)
            global DN = D2nd(N_lowres+3);
            global D0 = DN[:,1];
            global D_tilde = DN[:,2:N_lowres+1];
            global Db = DN[:,N_lowres+2:N_lowres+3];

            # correction of the regularization "data"
            global Δyd = -dropdims(sum(Db,dims=2),dims=2)*ρA_1[end];
            global δyd = -D0*ρA_1[1];
            global yd = Δyd+δyd;


            # smoosh together the several part of the model into augmented operators
            global Htrunc = [H_tilde; D_tilde];                                     # conditional to data and measurement model

            #
            # covariances: the crafed Bayesian models assumes that some covariance are known, i.e. measurement noise, smoothness, known values and measurement operator (the last only in the marginalized case)
            #

            # measurement noise covariance
            global ΓI = σnoise^2*diagm(ones(Cdouble,Ndata));
            global ΓItrunc = ΓI + σB^2*Hnot1*Hnot1';
            global ΓIinv = inv(ΓItrunc);

            # covariance matrix for the a priori distribution (second order difference)
            global Γprior_lowres = zeros(Cdouble,Nr_lowres,Nr_lowres);
            for i in 1:Nr_lowres
                Γprior_lowres[i,i] =  1.0;
                for j in i+1:Nr_lowres
                    Γprior_lowres[i,j] = Γprior_lowres[i,i]*exp(-(i-j)^2/(0.5*cor_len_lowres^2));
                    Γprior_lowres[j,i] = Γprior_lowres[i,j];
                end
            end
            global Γd_lowres = (N_lowres/Ndata)*(σd^2)*Γprior_lowres[2:N_lowres+2,2:N_lowres+2];  # scale the a priori strength with the quantity of data, so that it is possible to compare the results
            global Γd_lowres_inv = inv(Γd_lowres);


            #
            # data inversion: estimation of the concentration profile using CP algorithm
            #
            include("data_inversion_noise_marginalization.jl")
            

            #
            # marginalization over the meas. op. space
            #
            include("data_inversion_model_marginalization.jl")     
            
            #
            # save results
            #
            include("save_results.jl")

            #
            # plot and save figures
            #
            include("plot_estimation.jl")
            close("all")
        end
    end
end

