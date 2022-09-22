# for each cases (cylindrical near and far, and finger), simulate the acquisition of one datum for each profile
M_far  = [H_r_far'*reverse(ρA_1); H_r_far'*reverse(ρA_2); H_r_far'*reverse(ρA_3); H_r_far'*reverse(ρA_4)];
M_near = [H_r_near'*reverse(ρA_1); H_r_near'*reverse(ρA_2); H_r_near'*reverse(ρA_3); H_r_near'*reverse(ρA_4)];
M_far_sphere  = [H_r_sphere_far'*reverse(ρA_1); H_r_sphere_far'*reverse(ρA_2); H_r_sphere_far'*reverse(ρA_3); H_r_sphere_far'*reverse(ρA_4)];
M_near_sphere = [H_r_sphere_near'*reverse(ρA_1); H_r_sphere_near'*reverse(ρA_2); H_r_sphere_near'*reverse(ρA_3); H_r_sphere_near'*reverse(ρA_4)];
M_z    = [H_z'*reverse(ρA_1); H_z'*reverse(ρA_2); H_z'*reverse(ρA_3); H_z'*reverse(ρA_4)];

## plot
figure(figsize=[10,5])
ax1 = subplot(121)
l_cylinder_far,    = plot(r_surf,H_r_far/maximum(H_r_far), label="cylinder", color="tab:blue"); # : \$\\lambda_e\$
# l_cylinder_near,   = plot(r_surf,H_r_near/maximum(H_r_near), label="cylinder: near \$\\lambda_e\$", color="tab:orange");
# l_sphere_near,     = plot(r_surf,H_r_sphere_near/maximum(H_r_sphere_near), label="sphere: near \$\\lambda_e\$", color="tab:red");
l_sphere_far,      = plot(r_surf,H_r_sphere_far/maximum(H_r_sphere_far), label="sphere", color="tab:orange"); # purple # : far \$\\lambda_e\$
l_finger,          = plot(r_surf,H_z/maximum(H_z), label="pointwise", color="tab:green"); #  \$\\lambda_e\$
# s = @sprintf "planar approximation limits (far)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_far r_eal_max_far
# l_cylinder_far_un  = fill_between(r_surf,H_z_far_min/maximum(H_z_far_min),H_z_far_max/maximum(H_z_far_max),alpha=0.5,color="tab:blue",label=s)
# s = @sprintf "planar approximation limits (near)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_near r_eal_max_near
# l_cylinder_near_un = fill_between(r_surf,H_z_near_min/maximum(H_z_near_min),H_z_near_max/maximum(H_z_near_max),alpha=0.5,color="tab:orange",label=s)
# legend([(l_cylinder_far,l_cylinder_far_un),(l_cylinder_near,l_cylinder_near_un),l_sphere_near,l_sphere_far,l_finger],["cylinder: far \$\\lambda_e\$","cylinder: near \$\\lambda_e\$","sphere: far \$\\lambda_e\$","sphere: near \$\\lambda_e\$","pointwise"],fontsize=14)
# legend([l_cylinder_far,l_cylinder_near,l_sphere_near,l_sphere_far,l_finger],["cylinder: far \$\\lambda_e\$","cylinder: near \$\\lambda_e\$","sphere: far \$\\lambda_e\$","sphere: near \$\\lambda_e\$","pointwise"],fontsize=14)
legend(fontsize=14,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
xlabel("radial distance [\$\\mu\$m]",fontsize=16)
ylabel("normalized gain [a.u.]",fontsize=16)
xticks(fontsize=16)
yticks(fontsize=16)
ax1.xaxis.offsetText.set_size(12)

ax2 = subplot(122)
scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_far./M_z, label="cylinder vs pointwise", color="tab:blue")
# scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_near./M_z, label="cylinder near vs pointwise", color="tab:orange")
scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_far_sphere./M_z, label="sphere vs pointwise", color="tab:orange") # purple
# scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_near_sphere./M_z, label="sphere near vs pointwise", color="tab:red")
xlabel("concentration profiles",fontsize=16)
ylabel("relative acquisition",fontsize=16)
xticks(fontsize=16)
yticks(fontsize=16)
ylim(0.0)
legend(fontsize=14,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
ax2.ticklabel_format(axis="y",style="sci", scilimits=(-2, 2))
ax2.yaxis.offsetText.set_size(12)

tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
ax2.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-1.27-0.08-0.05, 0.99), textcoords="axes fraction", color="black",fontsize=16)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.07-0.01-0.10, 0.99), textcoords="axes fraction", color="black",fontsize=16)
ax1.text(0.25, 0.4, "\$\\lambda_e\$ = 2[nm]", transform=ax1.transAxes,fontsize=20)
ax2.text(0.15, 0.6, "\$\\frac{H_{\\mathrm{cylinder}}(\\rho_i,\\lambda_e)}{H_{\\mathrm{pointwise}}(\\rho_i,\\lambda_e)}\$ and \$\\frac{H_{\\mathrm{sphere}}(\\rho_i,\\lambda_e)}{H_{\\mathrm{pointwise}}(\\rho_i,\\lambda_e)}\$", transform=ax2.transAxes,fontsize=20)


# savefig("plans_vs_cylinder.png")
# savefig("plans_vs_cylinder.pdf")

# savefig("pointwise_vs_cylinder_vs_sphere_radius_20_microns.png")
# savefig("pointwise_vs_cylinder_vs_sphere_radius_20_microns.pdf")
