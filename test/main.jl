# load simulated data, and fit the peaks, estimate the alignment factor, remove the noise and estimate the sensitivity matrices
## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)
rc("figure",max_open_warning=50)
using myPlot

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

FLAG_0001 = true;
FLAG_0002 = false;
FLAG_0003 = false;
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

dictAllData = Dict();
for xf_name in xf_data_sheet_names
    local df = DataFrame(Matrix{Cdouble}(XLSX.getdata(xf_data[xf_name])[2:end,:]),:auto)
    rename!(df,Symbol.(XLSX.getdata(xf_data[xf_name])[1,:]));
    dictAllData[Symbol(xf_name)] = (eachcol(df),names(df))
end

# xf_model = XLSX.readxlsx(string(data_folder,"model.xlsx"));
# xf_model_sheet_names = XLSX.sheetnames(xf_model);
# dictAllGeom = Dict();
# for xf_name in xf_model_sheet_names
#     local df = DataFrame(XLSX.getdata(xf_model[xf_name])[2:end,:],:auto)
#     rename!(df,Symbol.(XLSX.getdata(xf_model[xf_name])[1,:]));
#     dictAllGeom[Symbol(xf_name)] = (eachcol(df),names(df))
# end


Ndata = length(dictAllData);
τm = [0.85; 0.125; 1.0-0.85-0.125];  # [1.0/3.0; 1.0/3.0; 1.0/3.0];
μm = [290.2; 292.0; 293.0]; # [290.3; 291.9; 293.5];
σm = sqrt(2.0)*[0.45; 0.25; 0.6]; # [290.3; 291.9; 293.5]/500.0;
τt = zeros(Cdouble,Ndata,3);
μt = zeros(Cdouble,Ndata,3);
σt = zeros(Cdouble,Ndata,3);
key_symbol = Symbol.(keys(dictAllData))
for j in 1:Ndata
    local be = dictAllData[key_symbol[j]][1][:Be];
    local spectrum = dictAllData[key_symbol[j]][1][:Snoisy]-dictAllData[key_symbol[j]][1][:Sbg];
   # estimate the peaks centers and spreads
   τt[j,:],μt[j,:],σt[j,:] = EM_peaks(be,spectrum,τm,μm,σm,200)
end

