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
# using StatsBase


COMPUTE_SENSITIVITY = false;


## load the results 
FLAG_0001 = true             # selection of the profile (one must be true and the others false)
FLAG_0002 = false
FLAG_0003 = false
FLAG_0004 = false
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

FLAG_LIST = ("0001","0002","0003","0004");
i_flag = 1

MODEL_TYPE_LIST  = ("5_datapoints_short_range","5_datapoints_wide_range","20_datapoints_wide_range")
MODEL_FOLD_LIST  = ("eal_5_restricted_range/","eal_5/","eal_20/")
i_eal = 3

MODEL_ERROR_LIST = ("model_error_0.5","model_error_1.0","model_error_2.5")
MODEL_ERROR_FOLD = ("error_model_0.005_percent/","error_model_0.01_percent/","error_model_0.025_percent/")
i_err = 1

NOISE_LEVEL_LIST = ("noise_level_0.01/","noise_level_0.1/","noise_level_0.5/")
noise_sigma      = (0.01,0.1,0.5)
i_noise = 2




prior_strength   = "weak_prior"


LEGEND_FONTSIZE = 14
LABEL_FONTSIZE  = 16
TICKS_FONTSIZE  = 16
ANOTA_FONTSIZE  = 16
fig0 = figure(1,figsize=[10,15]); #[10,6]

for i_flag in 2:4
    mod_str      = MODEL_TYPE_LIST[i_eal];
    mod_fol      = MODEL_FOLD_LIST[i_eal];
    mod_err      = MODEL_ERROR_LIST[i_err];
    mod_err_fold = MODEL_ERROR_FOLD[i_err];
    dat_noise    = NOISE_LEVEL_LIST[i_noise];
    σnoise       = noise_sigma[i_noise]
    σd           = 0.1
    cor_len_lowres = 2.5
    save_folder  = string("results/",FLAG_LIST[i_flag],"/",mod_str,"/",mod_err,"/",dat_noise,prior_strength,"/")
    model_folder           = string("../data/",mod_fol)
    model_folder_lowres    = string(model_folder,"lowres/")              # load the low resolution models from this folder
    model_folder_lowres_un = string(model_folder_lowres,mod_err_fold)

    ## model 
    dfH_lowres = CSV.File(string(model_folder_lowres_un,"H_lowres.csv");header=true,ntasks=1) |> DataFrame
    H_lowres   = Matrix{Cdouble}(dfH_lowres);

    ## data (for computing the SNR)
    data_folder  = string(model_folder,FLAG_LIST[i_flag],"/",dat_noise)
    dfRepData = CSV.File(string(data_folder,"repeated_data.csv");header=true,ntasks=1) |> DataFrame;
    repData   = Matrix{Cdouble}(dfRepData)[2:end,:];
    λe        = Matrix{Cdouble}(dfRepData)[1,:];
    Ndata     = length(λe);
    Nrep      = size(repData,1);
    Nsample   = min(50,Nrep);
    i_sample  = min(2,Nsample);

    ## results
    dfr_lowres = CSV.File(string(save_folder,"depth_lowres.csv");header=true,ntasks=1) |> DataFrame ;
    r_lowres   = dropdims(Matrix{Cdouble}(dfr_lowres),dims=1);
    Nr_lowres  = length(r_lowres);
    df_ρ       = CSV.File(string(save_folder,"concentration_profile.csv");header=true,ntasks=1) |> DataFrame ;
    ρA_1       = dropdims(Matrix{Cdouble}(df_ρ),dims=1);
    df_r       = CSV.File(string(save_folder,"depth.csv");header=true,ntasks=1) |> DataFrame ;
    r          = dropdims(Matrix{Cdouble}(df_r),dims=1);
    μ0         = r[1]; 
    d0         = 5.0e-3 # 15.0e-3 # NOTE: this value should depend on the penetration depth
    N0         = findfirst(r.-μ0.<=-d0);
    N0_lowres  = findfirst(r_lowres.-μ0.<=-d0);
    N_lowres   = N0_lowres-1;

    df_ρ_cp = CSV.File(string(save_folder,"concentration_estimates.csv");header=true,ntasks=1) |> DataFrame ;
    ρ_cp = Matrix{Cdouble}(df_ρ_cp);

    df_mean_ρ_H = CSV.File(string(save_folder,"concentration_noise_variability_mean.csv");header=true,ntasks=1) |> DataFrame ;
    mean_ρ_H    = dropdims(Matrix{Cdouble}(df_mean_ρ_H),dims=1);
    df_var_ρ_H = CSV.File(string(save_folder,"concentration_noise_variability_cov.csv");header=true,ntasks=1) |> DataFrame ;
    var_ρ_H    = Matrix{Cdouble}(df_var_ρ_H);

    df_μρ_H = CSV.File(string(save_folder,"concentration_distribution_mean_conditional_to_model.csv");header=true,ntasks=1) |> DataFrame ;
    μρ_H    = dropdims(Matrix{Cdouble}(df_μρ_H),dims=1);
    df_Γρ_H = CSV.File(string(save_folder,"concentration_distribution_cov_conditional_to_model.csv");header=true,ntasks=1) |> DataFrame ;
    Γρ_H    = Matrix{Cdouble}(df_Γρ_H);

    df_μρ_HI = CSV.File(string(save_folder,"concentration_distribution_mean_one_model.csv");header=true,ntasks=1) |> DataFrame ;
    μρ_HI = Matrix{Cdouble}(df_μρ_HI)';
    Γρ_HI = zeros(Cdouble,N_lowres,N_lowres,Nsample);
    for sample_flag in 1:Nsample
        local df_Γρ_HI = CSV.File(string(save_folder,"concentration_distribution_mean_one_model_sample_",sample_flag,".csv");header=true,ntasks=1) |> DataFrame ;
        Γρ_HI[:,:,sample_flag] = Matrix{Cdouble}(df_Γρ_HI);
    end


    df_ρ_cp_HI = CSV.File(string(save_folder,"concentration_model_variability_estimates.csv");header=true,ntasks=1) |> DataFrame;
    ρ_cp_HI = Matrix{Cdouble}(df_ρ_cp_HI);

    df_mean_ρ_y = CSV.File(string(save_folder,"concentration_model_variability_mean.csv");header=true,ntasks=1) |> DataFrame;
    mean_ρ_y    = dropdims(Matrix{Cdouble}(df_mean_ρ_y),dims=1);
    df_var_ρ_y = CSV.File(string(save_folder,"concentration_model_variability_cov.csv");header=true,ntasks=1) |> DataFrame;
    var_ρ_y = Matrix{Cdouble}(df_var_ρ_y);

    df_μρ_y = CSV.File(string(save_folder,"concentration_distribution_mean_conditional_to_noise.csv");header=true,ntasks=1) |> DataFrame;
    μρ_y    = dropdims(Matrix{Cdouble}(df_μρ_y),dims=1);
    df_Γρ_y = CSV.File(string(save_folder,"concentration_distribution_cov_conditional_to_noise.csv");header=true,ntasks=1) |> DataFrame;
    Γρ_y = Matrix{Cdouble}(df_Γρ_y);

    df_μρ_HI_sample = CSV.File(string(save_folder,"concentration_distribution_mean_one_data.csv");header=true,ntasks=1) |> DataFrame;
    μρ_HI_sample = Matrix{Cdouble}(df_μρ_HI_sample);
    Γρ_HI_sample = zeros(Cdouble,N_lowres,N_lowres,Nsample);
    for sample_flag in 1:Nsample
        local df_Γρ_HI_sample = CSV.File(string(save_folder,"concentration_distribution_mean_one_data_sample_",sample_flag,".csv");header=true,ntasks=1) |> DataFrame ;
        Γρ_HI_sample[:,:,sample_flag] = Matrix{Cdouble}(df_Γρ_HI_sample);
    end



    ## sensitivity matrix
    if COMPUTE_SENSITIVITY
        Ndata = size(H_lowres,1);
        H_tilde = H_lowres[:,2:N0_lowres];
        ΓI = σnoise^2*diagm(ones(Cdouble,Ndata));
        # ΓI = σnoise^2*diagm(collect(Ndata:-1.0:1.0)); # that would be a good experimental case (if the noise level decreases with the attenuation length, we an probe deeper)
        # ΓI = σnoise^2*diagm(collect(1.0:Ndata));
        DN = D2nd(N_lowres+3)
        D_tilde = DN[:,2:N_lowres+1]
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
        J_theory = H_tilde'inv(ΓI)*H_tilde + D_tilde'*Γd_lowres_inv*D_tilde
        Γ_theory = 2.0inv(J_theory)
        figure(); imshow(Γ_theory); colorbar()
        # F_theory = svd(J_theory);
        F_theory = svd(H_tilde'inv(ΓI)*H_tilde);
        S_th = F_theory.S;
        S_th[2:end] .= 0.0;
        J_trunc = F_theory.U*diagm(S_th)*F_theory.Vt;
        figure(); imshow(J_trunc); colorbar()
    end


    ## let's get plotting
    #
    # estimates and uncertainties
    #
    ##
    ## ax1         = subplot(121) # variability of the estimates
    ##
    if i_flag==2
        ax1         = subplot(321)
    elseif i_flag==3
        ax1         = subplot(323)
    else
        ax1         = subplot(325)
    end

    l_ρ_cp,     = plot(1000.0(μ0.-r_lowres),ρ_cp[i_sample,:],color="tab:blue")
    # conditional to y: variability
    l_mean_ρ_y, = plot(1000.0(μ0.-r_lowres),mean_ρ_y,color="tab:red")
    l_var_ρ_y   = fill_between(1000.0(μ0.-r_lowres),mean_ρ_y-sqrt.(diag(var_ρ_y)),mean_ρ_y+sqrt.(diag(var_ρ_y)),alpha=0.65,color="tab:red")
    # conditional to H: variability
    l_mean_ρ_H, = plot(1000.0(μ0.-r_lowres),mean_ρ_H,color="tab:orange")
    l_var_ρ_H   = fill_between(1000.0(μ0.-r_lowres),mean_ρ_H-sqrt.(diag(var_ρ_H)),mean_ρ_H+sqrt.(diag(var_ρ_H)),alpha=0.25,color="tab:orange")
    # ground truth
    l_gt,       = plot(1000.0(μ0.-r),ρA_1,color="tab:green")
    l_variability = [l_ρ_cp,(l_mean_ρ_y,l_var_ρ_y),(l_mean_ρ_H,l_var_ρ_H),l_gt];
    l_names     = ["estimate \$\\hat{\\rho}|H,\\mathbf{y}\$","meas. op. var. \$\\hat{\\rho}|\\mathbf{y}\$","meas. noise var. \$\\hat{\\rho}|H\$","GT"]
    legend(l_variability,l_names,fontsize=LEGEND_FONTSIZE,loc="upper right",borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
    # legend([(l_cp_post_est,l_cp_post_cov),(l_cp_post_est_y,l_cp_post_cov_y),(l_cp_noise_mean,l_cp_noise_cov),(l_cp_var_mod,l_cp_var_mod_cov),l_gt],["sampled posterior P\$(\\rho|H,y)\$","sampled posterior P\$(\\rho|y)\$","\$E[\\hat{\\rho}|H]\\pm\\Gamma[\\hat{\\rho}|H]\$","\$E[\\hat{\\rho}|y]\\pm\\Gamma[\\hat{\\rho}|y]\$","GT"],fontsize=14,loc="upper right")

    # xlim(0.0,8.5)
    xlim(0.0,6.0)
    # ylim(-0.1,1.5maximum(ρA_1))
    # ylim(-0.1,3.75maximum(ρA_1))
    if i_flag==1 #FLAG_0001
        ylim(-0.1,3.75maximum(ρA_1))
    elseif i_flag==2 #FLAG_0002
        ylim(-0.1,1.5maximum(ρA_1))
    elseif i_flag==3 #FLAG_0003
        ylim(-0.1,2.25maximum(ρA_1))
    else
        ylim(-0.1,1.49maximum(ρA_1))
    end
    if i_flag==4 #FLAG_0004
        xlabel("depth [nm]",fontsize=LABEL_FONTSIZE)
    end
    xticks(fontsize=TICKS_FONTSIZE)
    ylabel("concentration [a.u.]",fontsize=LABEL_FONTSIZE)
    yticks(fontsize=TICKS_FONTSIZE)

    ##
    ## ax2 = subplot(122) # distributions and their marginalizations
    ##
    if i_flag==2 # FLAG_0002
        ax2         = subplot(322)
    elseif i_flag==3 # FLAG_0003
        ax2         = subplot(324)
    else
        ax2         = subplot(326)
    end
    # conditional to y and H
    l_μρ_HI, = plot(1000.0(μ0.-r_lowres)[2:N0_lowres],μρ_HI[:,i_sample],color="tab:blue")
    l_Γρ_HI  = fill_between(1000.0(μ0.-r_lowres)[2:N0_lowres],μρ_HI[:,i_sample]-sqrt.(diag(Γρ_HI[:,:,i_sample])),μρ_HI[:,i_sample]+sqrt.(diag(Γρ_HI[:,:,i_sample])),alpha=0.25,color="tab:blue")
    # conditional to y: mean and covariance of the distribution
    l_μρ_y, = plot(1000.0(μ0.-r_lowres)[2:N0_lowres],μρ_y,color="tab:red")
    l_Γρ_y  = fill_between(1000.0(μ0.-r_lowres)[2:N0_lowres],μρ_y-sqrt.(diag(Γρ_y)),μρ_y+sqrt.(diag(Γρ_y)),alpha=0.25,color="tab:red")
    # conditional to H: mean and covariance of the distribution
    l_μρ_H,   = plot(1000.0(μ0.-r_lowres[2:N0_lowres]),μρ_H,color="tab:orange")
    # (Ns-Ns_burn)/(Ns-Ns_burn-1)*(Γρ_H+μρ_H*μρ_H') - μρ_H*μρ_H'
    l_Γρ_H    = fill_between(1000.0(μ0.-r_lowres[2:N0_lowres]),μρ_H-sqrt.(diag(Γρ_H)),μρ_H+sqrt.(diag(Γρ_H)),alpha=0.5,color="tab:orange")
    # ground truth
    l_gt,       = plot(1000.0(μ0.-r),ρA_1,color="tab:green")

    l_post    = [(l_μρ_HI,l_Γρ_HI),(l_μρ_y,l_Γρ_y),(l_μρ_H,l_Γρ_H),l_gt];
    l_names   = ["P\$(\\rho|H,\\mathbf{y})\$","P\$(\\rho|\\mathbf{y})\$","P\$(\\rho|H)\$","GT"]
    legend(l_post,l_names,fontsize=LEGEND_FONTSIZE,loc="upper right",borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)

    # xlim(0.0,8.5)
    xlim(0.0,6.0)
    # ylim(-0.1,1.5maximum(ρA_1))
    # ylim(-0.1,3.75maximum(ρA_1))
    if i_flag==1 # FLAG_0001
        ylim(-0.1,3.75maximum(ρA_1))
    elseif i_flag==2 # FLAG_0002
        ylim(-0.1,1.5maximum(ρA_1))
    elseif i_flag==3 # FLAG_0003
        ylim(-0.1,2.25maximum(ρA_1))
    else
        ylim(-0.1,1.49maximum(ρA_1))
    end
    if i_flag==4 #FLAG_0004
        xlabel("depth [nm]",fontsize=LABEL_FONTSIZE)
    end
    xticks(fontsize=TICKS_FONTSIZE)
    # ylabel("concentration [a.u.]",fontsize=LABEL_FONTSIZE)
    yticks(fontsize=TICKS_FONTSIZE)


    if i_flag==2 #FLAG_0002
        ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.1-0.05, 0.975+0.01), textcoords="axes fraction", color="black",fontsize=ANOTA_FONTSIZE)
        ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.1-0.015, 0.975+0.01), textcoords="axes fraction", color="black",fontsize=ANOTA_FONTSIZE)
    elseif i_flag==3 # FLAG_0003
        ax1.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(-0.1-0.05, 0.975+0.01), textcoords="axes fraction", color="black",fontsize=ANOTA_FONTSIZE)
        ax2.annotate("d)", xy=(3, 1),  xycoords="data", xytext=(-0.1-0.015, 0.975+0.01), textcoords="axes fraction", color="black",fontsize=ANOTA_FONTSIZE)
    else
        ax1.annotate("e)", xy=(3, 1),  xycoords="data", xytext=(-0.1-0.05, 0.975+0.01), textcoords="axes fraction", color="black",fontsize=ANOTA_FONTSIZE)
        ax2.annotate("f)", xy=(3, 1),  xycoords="data", xytext=(-0.1-0.015, 0.975+0.01), textcoords="axes fraction", color="black",fontsize=ANOTA_FONTSIZE)
    end
end
tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
# fig0.savefig(string(save_folder,"estimates_and_uncertainties.png"))
# fig0.savefig(string(save_folder,"estimates_and_uncertainties.pdf"))
# fig0.savefig(string("results/0001/estimates_and_uncertainties_N5_W5_W20.png"))
# fig0.savefig(string("results/0001/estimates_and_uncertainties_N5_W5_W20.pdf"))
# fig0.savefig(string("results/0001/estimates_and_uncertainties_W20_model_err.png"))
# fig0.savefig(string("results/0001/estimates_and_uncertainties_W20_model_err.pdf"))
# fig0.savefig(string("results/0001/estimates_and_uncertainties_W20_noise_err.png"))
# fig0.savefig(string("results/0001/estimates_and_uncertainties_W20_noise_err.pdf"))
# fig0.savefig(string("results/estimates_and_uncertainties_W20_FLAG234.png"))
# fig0.savefig(string("results/estimates_and_uncertainties_W20_FLAG234.pdf"))



#
# distribution covariances
#
if false

    
    fig = figure(2,figsize=[5,10]);
    _,ax1,pcm1 = displayCov(2,1000.0(μ0.-r_lowres[2:N0_lowres]),Γρ_HI[:,:,1];_sub=311)
    ylabel("depth [nm]",fontsize=14)
    _,ax2,pcm2 = displayCov(2,1000.0(μ0.-r_lowres[2:N0_lowres]),Γρ_y;_sub=312)
    ylabel("depth [nm]",fontsize=14)
    _,ax3,pcm3 = displayCov(2,1000.0(μ0.-r_lowres[2:N0_lowres]),Γρ_H;_sub=313)
    xlabel("depth [nm]",fontsize=14)
    ylabel("depth [nm]",fontsize=14)

    cb1 = setVerticalColorbar(fig,pcm1,0.17,0.815,0.04,0.15,"concentration\$^2\$ [a.u.]";fontsize_label=14,fontsize_ticks=12) 
    cb2 = setVerticalColorbar(fig,pcm2,0.17,0.495,0.04,0.15,"concentration\$^2\$ [a.u.]";fontsize_label=14,fontsize_ticks=12)
    cb3 = setVerticalColorbar(fig,pcm3,0.17,0.178,0.04,0.15,"concentration\$^2\$ [a.u.]";fontsize_label=14,fontsize_ticks=12)

    tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)

    ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)
    ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)
    ax3.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)

    # fig.savefig(string(save_folder,"distribution_covariance.png"))
    # fig.savefig(string(save_folder,"distribution_covariance.pdf"))
end

# mean(repData,dims=1)
# std(repData,dims=1)

# dropdims((mean(repData,dims=1).^2)./var(repData,dims=1),dims=1)


