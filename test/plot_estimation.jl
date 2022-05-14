#
# estimates and uncertainties
#
fig0 = figure(1,figsize=[10,6]); 
ax1         = subplot(121) # variability of the estimates
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
l_names     = ["estimate \$\\hat{\\rho}|H,y\$","meas. op. var. \$\\hat{\\rho}|y\$","meas. noise var. \$\\hat{\\rho}|H\$","GT"]
legend(l_variability,l_names,fontsize=14,loc="upper right")
# legend([(l_cp_post_est,l_cp_post_cov),(l_cp_post_est_y,l_cp_post_cov_y),(l_cp_noise_mean,l_cp_noise_cov),(l_cp_var_mod,l_cp_var_mod_cov),l_gt],["sampled posterior P\$(\\rho|H,y)\$","sampled posterior P\$(\\rho|y)\$","\$E[\\hat{\\rho}|H]\\pm\\Gamma[\\hat{\\rho}|H]\$","\$E[\\hat{\\rho}|y]\\pm\\Gamma[\\hat{\\rho}|y]\$","GT"],fontsize=14,loc="upper right")

xlim(0.0,8.5)
ylim(-0.1,1.5maximum(ρA_1))
xlabel("depth [nm]",fontsize=14)
xticks(fontsize=14)
ylabel("concentration [a.u.]",fontsize=14)
yticks(fontsize=14)

ax2 = subplot(122) # distributions and their marginalizations
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
l_names   = ["P\$(\\rho|H,y)\$","P\$(\\rho|y)\$","P\$(\\rho|H)\$","GT"]
legend(l_post,l_names,fontsize=14,loc="upper right")

xlim(0.0,8.5)
ylim(-0.1,1.5maximum(ρA_1))
xlabel("depth [nm]",fontsize=14)
xticks(fontsize=14)
ylabel("concentration [a.u.]",fontsize=14)
yticks(fontsize=14)


tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.975), textcoords="axes fraction", color="black",fontsize=14)

fig0.savefig(string(save_folder,"estimates_and_uncertainties.png"))
fig0.savefig(string(save_folder,"estimates_and_uncertainties.pdf"))








#
# distribution covariances
#
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

fig.savefig(string(save_folder,"distribution_covariance.png"))
fig.savefig(string(save_folder,"distribution_covariance.pdf"))
