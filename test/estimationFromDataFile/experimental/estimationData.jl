# load simulated data, and fit the peaks, estimate the alignment factor, remove the noise and estimate the sensitivity matrices
## load the packages used in the estimation
# plotting
using PyPlot
using PyCall
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
using LinearAlgebra
using StatsBase


# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)

# inversion package
using XPSinv

SAMPLING = true
if SAMPLING
    using XPSsampling
end

PLOT_FIG  = true
PLOT_DATA = false

# folders where the data files are
data_folder = "../../../data/TK/";
data_folderC1s = string(data_folder,"C1s/")
data_folderC1s = string(data_folder,"S2p/")

# flags
FLAG_PLOT = true;
FLAG_SAVE_PLOT = false;

# geometry setup
λe0 = 2.0e-3;         # reference penetration depth in μm
δr = 2.0e-3           # transition to vacuum layer thickness (let's set about 1 nm)
k0 = 10;              # compute the measurement model over a distance of k0*λe0
Nr = 101; # 201;      # number of discretization points in the radial dimension
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

# unit conversion constant (some of the quantities are in μm, some in L and some in Mbarn)
NA = 6.022e23;
κ_simple_units = 1.0e-37*NA; # simplified model
κ_units        = 1.0e-25*NA; # original model



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
include("loadDataAndAlignment.jl")

# now, you can run the inversion with the data

y_data_1 # SNR OK
y_data_2 # SNR OK ish
y_data_3 # SNR... well, there's no need to try this one
σ_noise  # standard deviation of the normalized noise


##
## run the inversion using only y_data_1
##
include("profileReconstruction.jl")

##
## plot the results
##
if PLOT_FIG
    figure(figsize=[10, 5])
    ax1 = subplot(121)
    l_ρ, = plot(1.0e3*(r.-μ0),ρC1s_bulk*ρ_cp,color=color_array[1])
    if SAMPLING
        l_ρ_std =  fill_between(1.0e3*(r[N0:end-1].-μ0),ρC1s_bulk*(ρ_est-stdρ_HI),ρC1s_bulk*(ρ_est+stdρ_HI),alpha=0.5,color=color_array[1])
    end
    l_λ = Array{PyObject,1}(undef,Ndata)
    for i in 1:Ndata
        l_λ[i], = plot(-[λ_all[i]; λ_all[i]],[0.0; 1.0], color=color_array[i+1])
    end
    # l_s, = plot(-[0.0; 0.0],[0.0; 1.0], color=color_array[Ndata+2]) # ;"sharp edge surface"
    xlim(1.0e3*(r[1].-μ0),2.5e3*(r[end].-μ0))
    ylim(-0.01ρC1s_bulk,max(1.6ρC1s_bulk,1.1maximum(ρC1s_bulk*(ρ_est+stdρ_HI))))
    xlabel("depth [nm]",fontsize=14); 
    ylabel("molar concentraion [M]",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    # ax1 = gca()
    ax1.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax1.yaxis.offsetText.set_size(14)
    ax1.xaxis.offsetText.set_size(14)
    # legend(fontsize=12)
    legend([(l_ρ,l_ρ_std);l_λ],["estimates\$\\pm\\sigma\$ [M]"; string.("\$\\lambda_e\$ = ",floor.(100λ_all)/100," [nm]")],fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)

    ax2 = subplot(122)
    scatter(λ_all,ρC1s_bulk*y_data_1)
    ylim(0.9ρC1s_bulk*minimum(y_data_1),1.1ρC1s_bulk*maximum(y_data_1))
    xlabel("IMFP [nm]",fontsize=14); 
    ylabel("peak area [mol]",fontsize=14)
    xticks(fontsize=14); yticks(fontsize=14); 
    ax2.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax2.yaxis.offsetText.set_size(14)
    ax2.xaxis.offsetText.set_size(14)
    legend(["data"],fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
    

    tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
    ax1.text(-0.14, 0.97, "a)", transform=ax1.transAxes,fontsize=16)
    ax2.text(0.05, 0.1+0.75, replace(data_filesC1s[idx_file][1:end-5],"_"=>" "), transform=ax2.transAxes,fontsize=16)
    ax2.text(-0.14, 0.97, "b)", transform=ax2.transAxes,fontsize=16)

    savefig(string(data_filesC1s[idx_file][1:end-5],"_reconstruction_and_data_new.png"))
    savefig(string(data_filesC1s[idx_file][1:end-5],"_reconstruction_and_data_new.pdf"))

    if SAMPLING
        figure()
        plot(deltaU)
    end
end

