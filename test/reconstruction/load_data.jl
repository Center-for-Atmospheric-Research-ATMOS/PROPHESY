dfRepData = CSV.File(string(data_folder,"repeated_data.csv");header=true,ntasks=1) |> DataFrame;
repData   = Matrix{Cdouble}(dfRepData)[2:end,:];
λe        = Matrix{Cdouble}(dfRepData)[1,:];
Ndata     = length(λe);
Nrep      = size(repData,1);

# number of samples used for the estimation of the variability and marginal distributions
Nsample        = min(n_noise_sample,Nrep);   # measurement repetition
i_sample       = min(2,Nsample);             # select a sample with noise (the first sample is not corrupted by noise, but it's not a good choice because not realistic)
N_model_sample = min(n_model_sample,100);    # I have computed 100 model sampels, so any natural smaller than 100 is good (none of the samples is the true model)


# load some measurement model
dfr_lowres = CSV.File(string(model_folder_lowres,"radial_discretization_lowres.csv");header=true,ntasks=1) |> DataFrame
r_lowres   = dropdims(Matrix{Cdouble}(dfr_lowres),dims=1)
Nr_lowres = length(r_lowres);
dfH_lowres = CSV.File(string(model_folder_lowres_un,"H_lowres.csv");header=true,ntasks=1) |> DataFrame
H_lowres   = Matrix{Cdouble}(dfH_lowres);

# load the GT
dfr        = CSV.File(string(model_folder,"radial_discretization.csv");header=true,ntasks=1) |> DataFrame
dfRho      = CSV.File(string(model_folder,profile_flag,"/concentration_profile.csv");header=true,ntasks=1) |> DataFrame
r = dropdims(Matrix{Cdouble}(dfr),dims=1);
μ0 = r[1]; # for now, the first radial discretization distance is exactly on the boundary of the cylinder
ρA_1 = reverse(dropdims(Matrix{Cdouble}(dfRho),dims=1));