## load package
using PyPlot
using PyCall
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
PLOT_DATA = true
BETTER_MODEL = true # makes a big difference (it could not be possible to compute a reconstruction if the model would not include the low density liquid around the sharp edge volume)
SAMPLING = true
if SAMPLING
    using XPSsampling
end

data_folder = "../../../../data/cylinder_radius_10.0/peak_shift/eal_5_restricted_range/"
FLAG_0001 = false
FLAG_0002 = false
FLAG_0003 = false
FLAG_0004 = false
FLAG_0005 = true

##
## load data 
##
include("load_data.jl") # load both S2p and O1s, but only use the S2p data for profile reconstruction

##
## create measurement data from meta data
##
include("measurementModel.jl")


# normalize data # will be filled with the blood of my enemies... or more ethically, just the peak areas
y_data_1 = y_peak_1./(EssS2p[1,3]*t_all.*α_all.*T_all.*F_all.*σ_all.*ρB_al*κ_units);
y_data_2 = y_peak_2./(EssS2p[2,3]*t_all.*α_all.*T_all.*F_all.*σ_all.*ρB_al*κ_units);
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
    plot(r.-μ0,ρ_cp,color=color_array[1])
    if SAMPLING
        # plot(r[N0:end-1],μρ_HI,color=color_array[2])
        # fill_between(r[N0:end-1],μρ_HI-stdρ_HI,μρ_HI+stdρ_HI,alpha=0.5,color=color_array[2])
        fill_between(r[N0:end-1].-μ0,ρ_cp[N0:end-1]-stdρ_HI,ρ_cp[N0:end-1]+stdρ_HI,alpha=0.5,color=color_array[2])
    end
    plot(dictAllGeomS2p[symbolDictS2p[Symbol(string("hν_",round(Int64,EssS2p[1,1])))]].r.-μ0,dictAllGeomS2p[symbolDictS2p[Symbol(string("hν_",round(Int64,EssS2p[1,1])))]].ρ/0.01,color=color_array[3])

    if SAMPLING
        figure()
        plot(deltaU)
    end
end

# println("percentage difference between data and predicted data ",100.0*(A_data-(Hgeom*dictAllGeomS2p[symbolDictS2p[:hν_958]].ρ[1:5:end]/0.01))./A_data)

# df_α =  DataFrame( "α_gt"=> α_gt_mean, "α_ratio" =>  α_ratio_mean, "x_offset" => xc_off[1:5:end]);
# if SAVE_DATA
#     # XLSX.writetable(string(data_folder,"/alignment_list.xlsx"), collect(eachcol(df_α)), names(df_α));
#     if USING_GT_PROFILE
#         XLSX.writetable(string(data_folder,"/alignment_list_new_eal_gt.xlsx"), collect(eachcol(df_α)), names(df_α));
#     else
#         XLSX.writetable(string(data_folder,"/alignment_list_new_eal.xlsx"), collect(eachcol(df_α)), names(df_α));
#     end
# end

if false
    figure(figsize=[10, 5])
    ax1 = subplot(121)
    l_ρ, = plot(1.0e3*(r.-μ0),ρB_al[1]*ρ_cp,color=color_array[1])
    if SAMPLING
        l_ρ_std =  fill_between(1.0e3*(r[N0:end-1].-μ0),ρB_al[1]*(ρ_est-stdρ_HI),ρB_al[1]*(ρ_est+stdρ_HI),alpha=0.5,color=color_array[1])
    end
    plot(1000.0*(dictAllGeomS2p[symbolDictS2p[:hν_532]].r.-μ0),dictAllGeomS2p[symbolDictS2p[:hν_532]].ρ,color=color_array[3])
    # plot(r.-μ0,ρB_al[1]*ρ_cp,color="darkgreen")
    l_λ = Array{PyObject,1}(undef,NdataS2p)
    for i in 1:NdataS2p
        l_λ[i], = plot(-[EssS2p[1:2:end,2][i]; EssS2p[1:2:end,2][i]],[0.0; 1.0], color=color_array[i+3])
    end
    xlim(1.0e3*(r[1].-μ0),2.5e3*(r[end].-μ0))
    ylim(-0.01ρB_al[1],max(1.6ρB_al[1],1.1maximum(ρB_al[1]*(ρ_est+stdρ_HI))))
    xlabel("depth [nm]",fontsize=14); 
    ylabel("molar concentraion [M]",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    # ax1 = gca()
    ax1.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax1.yaxis.offsetText.set_size(14)
    ax1.xaxis.offsetText.set_size(14)
    # legend(fontsize=12)
    legend([(l_ρ,l_ρ_std);l_λ],["estimates\$\\pm\\sigma\$ [M]"; string.("\$\\lambda_e\$ = ",floor.(100EssS2p[1:2:end,2])/100," [nm]")],fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)

    ax2 = subplot(122)
    scatter(EssS2p[1:2:end,2],ρB_al[1]*y_data_1)
    ylim(0.9ρB_al[1]*minimum(y_data_1),1.1ρB_al[1]*maximum(y_data_1))
    xlabel("IMFP [nm]",fontsize=14); 
    ylabel("peak area [mol]",fontsize=14)
    xticks(fontsize=14); yticks(fontsize=14); 
    ax2.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax2.yaxis.offsetText.set_size(14)
    ax2.xaxis.offsetText.set_size(14)
    legend(["data"],fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)


    tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
    ax1.text(-0.14, 0.97, "a)", transform=ax1.transAxes,fontsize=16)
    ax2.text(0.05, 0.1+0.75, "S 2p 10mM", transform=ax2.transAxes,fontsize=16)
    ax2.text(-0.14, 0.97, "b)", transform=ax2.transAxes,fontsize=16)

    if FLAG_0001
        save_path = string(data_folder,folder_tag,"/new_eal/0001/S2p/")
    end
    if FLAG_0002
        save_path = string(data_folder,folder_tag,"/new_eal/0002/S2p/")
    end
    if FLAG_0003
        save_path = string(data_folder,folder_tag,"/new_eal/0003/S2p/")
    end
    if FLAG_0004
        save_path = string(data_folder,folder_tag,"/new_eal/0004/S2p/")
    end
    if FLAG_0005
        save_path = string(data_folder,folder_tag,"/new_eal/0005/S2p/")
    end

    # savefig(string(save_path,"reconstruction_and_data.png"))
    # savefig(string(save_path,"reconstruction_and_data.pdf"))
end

