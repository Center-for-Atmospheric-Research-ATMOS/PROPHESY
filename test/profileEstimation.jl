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

# tags
MARG_UN   = true   # set to true to load the mean measurement operator as well as the covarainces


FLAG_0001 = false              # selection of the profile (one must be true and the others false)
FLAG_0002 = false
FLAG_0003 = false
FLAG_0004 = true

data_folder  = "../data/eal_10/0001/noise_level_0.001/";
model_folder = "../data/eal_10/";

# eal: attenuation_length.csv (already in the data)
# radii: lowres/radial_discretization_lowres.csv
# low resolution model with random error in the eal: lowres/error_model_0.005_percent/H_lowres.csv
# low resolution average model: lowres/error_model_0.005_percent/mean_H_lowres.csv
# low resolution covariance model: lowres/error_model_0.005_percent/cov/cov_H_lowres_i.csv


# load some data
dfRepData = CSV.File(string(data_folder,"repeated_data.csv");header=true) |> DataFrame;
repData   = Matrix{Cdouble}(dfRepData)[2:end,:];
λe        = Matrix{Cdouble}(dfRepData)[1,:];
Ndata     = length(λe);

# load some measurement model
dfr_lowres = CSV.File(string(model_folder,"lowres/radial_discretization_lowres.csv");header=true) |> DataFrame
r_lowres   = dropdims(Matrix{Cdouble}(dfr_lowres),dims=1)
Nr_lowres = length(r_lowres);
dfH_lowres = CSV.File(string(model_folder,"lowres/error_model_0.005_percent/H_lowres.csv");header=true) |> DataFrame
H_lowres   = Matrix{Cdouble}(dfH_lowres);
dfμH = CSV.File(string(model_folder,"lowres/error_model_0.005_percent/mean_H_lowres.csv");header=true) |> DataFrame
μH = Matrix{Cdouble}(Matrix{Cdouble}(dfμH)');
ΓH = Array{Cdouble}(undef,Nr_lowres,Nr_lowres,Ndata);
for i in 1:Ndata
    local df = CSV.File(string(model_folder,"lowres/error_model_0.005_percent/cov/cov_H_lowres_",i,".csv");header=false) |> DataFrame
    ΓH[:,:,i] = Matrix{Cdouble}(df)
end

# load the GT
dfr   = CSV.File(string(model_folder,"radial_discretization.csv");header=true) |> DataFrame
dfRho = CSV.File(string(model_folder,"concentration_profile.csv");header=true) |> DataFrame


r = dropdims(Matrix{Cdouble}(dfr),dims=1);
μ0 = r[1] # for now, the first radial discretization distance is exactly on the boundary of the cylinder
ρA_1 = dropdims(Matrix{Cdouble}(dfRho),dims=1);



# deeper than some distance, the signal is not likely to be disantangled
d0 = 15.0e-3 # NOTE: this value should depend on the penetration depth
N0 = findfirst(r.-μ0.<=-d0) 
N0_lowres = findfirst(r_lowres.-μ0.<=-d0) 