## load package
using PyPlot
rc("text", usetex=true)
rc("figure",max_open_warning=50)
color_array = ["tab:blue"; "tab:orange"; "tab:green"; "tab:red"; "tab:purple"; "tab:brown"; "tab:pink"; "tab:gray"; "tab:olive"; "tab:cyan"; "magenta"; "yellow"; "hotpink"; "darkmagenta"; "chartreuse"; "deepskyblue"; "navy"; "darkcyan"; "crimson"; "firebrick"]; 

# data manipulation (loading, writing, etc)
using Printf
using XPSfile
using DataFrames
using Query

# scientific package from the official Julia repositories
using LinearAlgebra
using StatsBase

# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)

# inversion package
using XPSinv

PLOT_FIG  = true
PLOT_DATA = false
BETTER_MODEL = true # makes a big difference (it could not be possible to compute a reconstruction if the model would not include the low density liquid around the sharp edge volume)
SAMPLING = false
if SAMPLING
    using XPSsampling
end

data_folder = "../../data/cylinder_radius_10.0/peak_shift/eal_5_restricted_range/"
FLAG_0001 = false
FLAG_0002 = false
FLAG_0003 = true
FLAG_0004 = false

##
## load data 
##
include("load_data.jl") # load both C1s and O1s, but only use the C1s data for profile reconstruction

##
## create measurement data from meta data
##
include("measurementModel.jl")


# normalize data # will be filled with the blood of my enemies... or more ethically, just the peak areas
y_data_1 = y_peak_1./(t_all.*α_all.*T_all.*F_all.*σ_all.*ρB_al*κ_units);
y_data_2 = y_peak_2./(t_all.*α_all.*T_all.*F_all.*σ_all.*ρB_al*κ_units);
y_data_3 = y_peak_3./(t_all.*α_all.*T_all.*F_all.*σ_all.*ρB_al*κ_units);
A_data   = A_tot./(t_all.*α_all.*T_all.*F_all.*σ_all.*ρB_al*κ_units);
σnoise   = σnoise./(t_all.*α_all.*T_all.*F_all.*σ_all.*ρB_al*κ_units);

##
## run the inversion
##

include("profileReconstruction.jl")

##
## plot the results
##
if PLOT_FIG
    figure()
    plot(r,ρ_cp,color=color_array[1])
    if SAMPLING
        plot(r[N0:end-1],μρ_HI,color=color_array[2])
        fill_between(r[N0:end-1],μρ_HI-stdρ_HI,μρ_HI+stdρ_HI,alpha=0.5,color=color_array[2])
    end
    plot(dictAllGeomC1s[symbolDictC1s[:hν_958]].r,dictAllGeomC1s[symbolDictC1s[:hν_958]].ρ/0.01,color=color_array[3])

    if SAMPLING
        figure()
        plot(deltaU)
    end
end

# println("percentage difference between data and predicted data ",100.0*(A_data-(Hgeom*dictAllGeomC1s[symbolDictC1s[:hν_958]].ρ[1:5:end]/0.01))./A_data)