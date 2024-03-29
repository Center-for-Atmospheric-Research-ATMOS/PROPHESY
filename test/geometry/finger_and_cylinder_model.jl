## compute the complete model
H_r_far,H_rθy_far,Arn,Aθj_far,Ayk = cylinder_gain_H(r_surf,θ_far,Y,x0_far,y0_far,z0_far,μ0,λe);
H_r_near,H_rθy_near,Arn,Aθj_near,Ayk = cylinder_gain_H(r_surf,θ_near,Y,x0_near,y0_near,z0_near,μ0,λe);
H_r_finger_far,_          = finger_gain_H(x0_far,y0_far,r_surf.-μ0,x0_far,y0_far,z0_far,λe);
H_r_finger_near,_         = finger_gain_H(x0_near,y0_near,r_surf.-μ0,x0_near,y0_near,z0_near,λe);

relative_eal_far  = log.(H_r_finger_far[2:end-1]/H_r_finger_far[end-1])./log.(H_r_far[2:end-1]/H_r_far[end-1]);
relative_eal_near = log.(H_r_finger_near[2:end-1]/H_r_finger_near[end-1])./log.(H_r_near[2:end-1]/H_r_near[end-1]);
r_eal_min_far,r_eal_max_far   = extrema(relative_eal_far[1:end-1])
r_eal_min_near,r_eal_max_near = extrema(relative_eal_near[1:end-1])

H_z,_                  = finger_gain_H(x0_far,y0_far,r_surf.-μ0,x0_far,y0_far,z0_far,λe);
H_z_far_min,_          = finger_gain_H(x0_far,y0_far,r_surf.-μ0,x0_far,y0_far,z0_far,r_eal_min_far*λe);
H_z_far_max,_          = finger_gain_H(x0_far,y0_far,r_surf.-μ0,x0_far,y0_far,z0_far,r_eal_max_far*λe);
H_z_near_min,_         = finger_gain_H(x0_near,y0_near,r_surf.-μ0,x0_near,y0_near,z0_near,r_eal_min_near*λe);
H_z_near_max,_         = finger_gain_H(x0_near,y0_near,r_surf.-μ0,x0_near,y0_near,z0_near,r_eal_max_near*λe);

## plot
figure();
plot(r_surf[2:end-1],relative_eal_far, label="far")
plot(r_surf[2:end-1],relative_eal_near, label="near")
legend()
xlabel("radius [\$\\mu m\$]")
ylabel("log ratio: \$\\frac{\\log(\\frac{I^{\\mathrm{planar}}}{I^{\\mathrm{planar}}_0} )}{\\log(\\frac{I^{\\mathrm{cylindrical}}}{I^{\\mathrm{cylindrical}}_0})}\$")

## plot
## compare the cylindrical acquisition gain density in the near and far cases
fig = figure(figsize=[9,5])
# near analyzer
ax1 = subplot(121,polar=true)
ax1.set_rlabel_position(-22.5)
ax1.set_rticks([μ0-0.03, μ0-0.02, μ0-0.01, μ0])
ax1.plot(atan(x0_far,z0_far)*ones(Cdouble,2),[μ0-0.03; μ0], color="red",label="\$\\theta\\simeq54.7\$")
pcm1 = ax1.pcolormesh(θ_near,r_surf,H_rθy_near[:,:,128],edgecolors="face")
ylim(μ0-0.03,μ0)
ax1.legend(loc="lower left", bbox_to_anchor=(.5 + cos(atan(x0_far,z0_far))/2, .5 + sin(atan(x0_far,z0_far))/2))
cax1 = fig.add_axes([0.15, .37, 0.02, 0.25])
cb1 = fig.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
cb1.set_label("gain [a.u.]", fontsize=10) # , color="white"

#far analyzer
ax2 = subplot(122,polar=true)
ax2.set_rlabel_position(-22.5)
ax2.set_rticks([μ0-0.03, μ0-0.02, μ0-0.01, μ0])
ax2.plot(atan(x0_far,z0_far)*ones(Cdouble,2),[μ0-0.03; μ0], color="red",label="\$\\theta\\simeq54.7\$")
pcm2 = ax2.pcolormesh(θ_far,r_surf,H_rθy_far[:,:,128],edgecolors="face")
ylim(μ0-0.03,μ0)
ax2.legend(loc="lower left", bbox_to_anchor=(.5 + cos(atan(x0_far,z0_far))/2, .5 + sin(atan(x0_far,z0_far))/2))
cax2 = fig.add_axes([0.65, .37, 0.02, 0.25])
cb2 = fig.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
cb2.set_label("gain [a.u.]", fontsize=10) # , color="white"

# annotation and squeezing
tight_layout(pad=1.0, w_pad=2.0, h_pad=0.2)
ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)

# savefig("gain_cylinder_near.png")
# savefig("gain_cylinder_near.pdf")

## plot
figure();
plot(r_surf,H_r_far/maximum(H_r_far), label="cylinder: far \$\\lambda_e\$", color="tab:blue");
plot(r_surf,H_r_near/maximum(H_r_near), label="cylinder: near \$\\lambda_e\$", color="tab:green");
plot(r_surf,H_z/maximum(H_z), label="planar approximation \$\\lambda_e\$", color="tab:orange");
s = @sprintf "planar approximation limits (far)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_far r_eal_max_far
fill_between(r_surf,H_z_far_min/maximum(H_z_far_min),H_z_far_max/maximum(H_z_far_max),alpha=0.5,color="tab:red",label=s)
s = @sprintf "planar approximation limits (near)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_near r_eal_max_near
fill_between(r_surf,H_z_near_min/maximum(H_z_near_min),H_z_near_max/maximum(H_z_near_max),alpha=0.5,color="tab:pink",label=s)
legend()
xlabel("radius [\$\\mu m\$]")
ylabel("normalized gain")
