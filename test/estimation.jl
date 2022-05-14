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


# TDTO: save files (see end of file), and loop over FLAG_NOISE, MODEL and MODEL_ERROR

# checking number of threads: launch Julia with the argument --threads 16 to run julia using 16 threads
println("you're running this script with ",Threads.nthreads()," threads") #WARNING: set ntasks to 1 when loading CSV files (it seems that multithreading is not safe with CSV.File)

SAVE_FIG = false;

# tags
FLAG_0001 = false             # selection of the profile (one must be true and the others false)
FLAG_0002 = false
FLAG_0003 = true
FLAG_0004 = false

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

n_noise_sample = 50; # number of noise samples for the noise marginalization
n_model_sample = 50; # number of noise samples for the noise marginalization

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
    profile_flag = "0001"
end

if FLAG_0002
    data_folder  = string(model_folder,"0002/")
    profile_flag = "0002"
end

if FLAG_0003
    data_folder  = string(model_folder,"0003/")
    profile_flag = "0003"
end

if FLAG_0004
    data_folder  = string(model_folder,"0004/")
    profile_flag = "0004"
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


## loop
MODEL_TYPE_LIST = ("5_datapoints_short_range","5_datapoints_wide_range","20_datapoints_wide_range")
MODEL_FOLD_LIST = ("eal_5_restricted_range/","eal_5/","eal_20/")

MODEL_ERROR_LIST = ("model_error_0.5","model_error_1.0","model_error_2.5")
MODEL_ERROR_FOLD = ("error_model_0.005_percent/","error_model_0.01_percent/","error_model_0.025_percent/")
NOISE_LEVEL_LIST = ("noise_level_0.01/","noise_level_0.1/","noise_level_0.5/")

for (mod_str,mod_fol) in zip(MODEL_TYPE_LIST,MODEL_FOLD_LIST)
    for (mod_err,mod_err_fold) in zip(MODEL_ERROR_LIST,MODEL_ERROR_FOLD)
        for dat_noise in NOISE_LEVEL_LIST
            # paths: load and save
            model_folder = string("../data/",mod_fol)                         # the common depth discretization should be in this folder
            model_folder_lowres = string(model_folder,"lowres/")              # load the low resolution models from this folder
            model_folder_lowres_un = string(model_folder_lowres,mod_err_fold) # 
            data_folder  = string(model_folder,profile_flag,"/",dat_noise)    # load the data from this folder
            save_folder  = string("results/",profile_flag,"/",mod_str,"/",mod_err,"/",dat_noise,prior_strength,"/") # save the results in this folder
            mkpath(save_folder)

            #
            # load some data
            #

            include("load_data.jl")

            #
            # truncation of the model and the data
            #

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

##

