# load simulated data, and fit the peaks, estimate the alignment factor, remove the noise and estimate the sensitivity matrices
## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)
rc("figure",max_open_warning=50)
# using myPlot
color_array = ["tab:blue"; "tab:orange"; "tab:green"; "tab:red"; "tab:purple"; "tab:brown"; "tab:pink"; "tab:gray"; "tab:olive"; "tab:cyan"; "magenta"; "yellow"; "hotpink"; "darkmagenta"; "chartreuse"; "deepskyblue"; "navy"; "darkcyan"; "crimson"; "firebrick"]; 

# data manipulation (loading, writing, etc)
using Printf
using XPSfile
using DataFrames
using Query

# scientific package from the official Julia repositories
using StatsBase

# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)


# import functions for importing and formatting data
include("loadFunctionXLSX2DataFrame.jl")


# folders where the data files are
data_folder = "../data/TK/";
data_folderC1s = string(data_folder,"C1s/")
data_folderO1s = string(data_folder,"O1s/")
data_folderS2p = string(data_folder,"S2p/")

# flags
FLAG_PLOT = true;
FLAG_SAVE_PLOT = true;

# geometry setup
λe0 = 2.0e-3;         # reference penetration depth in μm
δr = 2.0e-3           # transition to vacuum layer thickness (let's set about 1 nm)
k0 = 10;              # compute the measurement model over a distance of k0*λe0
Nr = 201;             # number of discretization points in the radial dimension
Nθ = 256;             # number of discretization points in the polar angle dimension
Ny = 256;             # number of discretization points in the cylinder axis dimension
L  = 50.0;            # height of the irradiated sample (the vertical extent of the beam is more like 20μm instead of 100μm)
μ0 = 10.0;            # radius of the cylinder
x0 = sqrt(2.0)*5000.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
y0 = 0.0;
z0 = 5000.0;

# spacial discretization 
r = collect(range(μ0-k0*λe0,μ0+δr,length=Nr));
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));



# column names fo the fitted results
photon_sym = Symbol("Photon energy");
ph_flu_sym = Symbol("Photon flux");
bind_sym   = Symbol("Binding energy");
shift_sym  = Symbol("Peak shift");
gauss_sym  = Symbol("FWHM(G)");
loren_sym  = Symbol("FWHM(L)");
area_sym   = Symbol("Area");
σtot_sym   = Symbol("Sigma");
α_al_sym   = Symbol("Alignment");

# regular expression to navigate the files
regData    = r"Eph ?= ?"                  # pattern found in the data sheet's name that does not appear in other sheets
regFit     = r"[Ff]itt?ing(_| )?results?" # pattern found in fit sheet's name
regEph     = r"=[0-9]*eV"                 # pattern for reading the photon energy
sortBy     = Symbol("Photon energy")      # column name for the photon energy  (sort the data by increasing order of photon energy)
thenSortBy = Symbol("Binding energy")     # colmun name for the binding energy (for the data sharing the same photon energy, sort them by binding energy)


# alignement from C1s data
include("alignmentFromDataC1s.jl")

# alignement from O1s data
include("alignmentFromDataO1s.jl")

# alignement from S2p data
include("alignmentFromDataS2p.jl")

# plotting
if FLAG_PLOT
    figure(figsize=[10, 10])
    ax = subplot(221)
    for i in 1:length(α_noiseO1s)
        local plot_sym = Symbol(data_filesO1s[i][1:end-5]);
        # scatter(1.0e9mean(α_noise[plot_sym]),mean(α_ratio[plot_sym]))
        scatter(1.0e9α_noiseO1s[plot_sym],α_ratioO1s[plot_sym],label=replace(data_filesO1s[i][1:end-5],"_"=>" ")) #.^2.5
    end
    xlim(0.1,10.0)
    ylim(0.1,10.0)
    xlabel("model estimation [x\$10^{9}\$]",fontsize=14); 
    ylabel("liq O1s/gas O1s",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax.yaxis.offsetText.set_size(14)
    ax.xaxis.offsetText.set_size(14)
    xscale("log")
    yscale("log")
    legend(fontsize=14)

    ax = subplot(222)
    for i in 1:length(α_noiseC1s)
        local plot_sym = Symbol(data_filesC1s[i][1:end-5]);
        scatter(1.0e9mean(α_noiseC1s[plot_sym]),mean(α_ratioC1s[plot_sym]),label=replace(data_filesC1s[i][1:end-5],"_"=>" ")) # .^2.5
    end
    xlim(0.1,10.0)
    ylim(0.1,10.0)
    xlabel("model estimation [x\$10^{9}\$]",fontsize=14); 
    ylabel("liq O1s/gas O1s",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax.yaxis.offsetText.set_size(14)
    ax.xaxis.offsetText.set_size(14)
    xscale("log")
    yscale("log")
    legend(fontsize=14)

    ax = subplot(224)
    for i in 1:length(α_noiseS2p)
        local plot_sym = Symbol(data_filesS2p[i][1:end-5]);
        scatter(1.0e9mean(α_noiseS2p[plot_sym]),mean(α_ratioS2p[plot_sym]),label=replace(data_filesS2p[i][1:end-5],"_"=>" ")) # .^2.5
    end
    xlim(3.0e-4,1.5e-2)
    ylim(0.1,10.0)
    xlabel("model estimation [x\$10^{9}\$]",fontsize=14); 
    ylabel("liq O1s/gas O1s",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax.yaxis.offsetText.set_size(14)
    ax.xaxis.offsetText.set_size(14)
    xscale("log")
    yscale("log")
    legend(fontsize=14)
    tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)

    if FLAG_SAVE_PLOT
        savefig(string("../data/TK/","ratio_vs_model_O1s_C1s_S2p_mean.png"))
        savefig(string("../data/TK/","ratio_vs_model_O1s_C1s_S2p_mean.pdf"))
    end
end
