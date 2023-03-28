# data manipulation (loading, writing, etc)
using Printf

# implemented scientific packages
using Statistics
using DataFrames
using CSV

# flags
PROFILE          = ("0001","0002","0003","0004");
MODEL_TYPE_LIST  = ("5_datapoints_short_range","5_datapoints_wide_range","20_datapoints_wide_range")
MODEL_ERROR_LIST = ("model_error_0.5","model_error_1.0","model_error_2.5")
NOISE_LEVEL_LIST = ("noise_level_0.01/","noise_level_0.1/","noise_level_0.5/")
noise_sigma      = (0.01,0.1,0.5)

# flag selection
i_profile  = 3;
i_mod_type = 3;
i_mod_err  = 1;
i_noise    = 1;

# result folder (to load data from)
data_folder = string("results/",PROFILE[i_profile],"/",MODEL_TYPE_LIST[i_mod_type],"/",MODEL_ERROR_LIST[i_mod_err],"/",NOISE_LEVEL_LIST[i_noise],"weak_prior/");

# load GT
df_r_hr = CSV.File(string(data_folder,"depth.csv");header=true,ntasks=1) |> DataFrame;
r_hr = dropdims(Array{Cdouble}(df_r_hr),dims=1);
df_ρ_GT = CSV.File(string(data_folder,"concentration_profile.csv");header=true,ntasks=1) |> DataFrame;
ρ_GT = dropdims(Array{Cdouble}(df_ρ_GT),dims=1);

# load estimates
df_r_lr = CSV.File(string(data_folder,"depth_lowres.csv");header=true,ntasks=1) |> DataFrame;
r_lr = dropdims(Array{Cdouble}(df_r_lr),dims=1);

df_ρ_cp = CSV.File(string(data_folder,"concentration_estimates.csv");header=true,ntasks=1) |> DataFrame;
ρ_cp = Matrix{Cdouble}(df_ρ_cp);

Nr = length(r_lr);
Ns = size(ρ_cp,1);


# compute amount of matter in the outter layer
idx_bb_lr = findfirst(r_lr.<(20.0-0.005)); # find first index where the bulk begins
idx_bb_hr = findfirst(r_hr.<(20.0-0.005))+1; # same for high resolution GT profile (to establish the GT)

dr_hr = r_hr[1]-r_hr[2]
q_GT  = sum(ρ_GT[1:idx_bb_hr])*dr_hr

dr_lr = r_lr[1]-r_lr[2]
q_cp = sum(ρ_cp[:,1:idx_bb_lr],dims=2)*dr_lr
μq,σq = mean(q_cp),std(q_cp);

println(data_folder)
println("The GT quantity of matter in the surface most layer:        ",q_GT)
println("The estimated quantity of matter in the surface most layer: ",μq," +/- ",σq)
println("the relative quantity and confidence interval:              ",μq/q_GT," +/-",σq/q_GT)

