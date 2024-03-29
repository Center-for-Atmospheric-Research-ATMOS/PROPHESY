# for each cases (cylindrical near and far, and finger), simulate the acquisition of one datum for each profile
M_far  = [H_r_far'*reverse(ρA_1); H_r_far'*reverse(ρA_2); H_r_far'*reverse(ρA_3); H_r_far'*reverse(ρA_4)];
M_near = [H_r_near'*reverse(ρA_1); H_r_near'*reverse(ρA_2); H_r_near'*reverse(ρA_3); H_r_near'*reverse(ρA_4)];
M_z    = [H_z'*reverse(ρA_1); H_z'*reverse(ρA_2); H_z'*reverse(ρA_3); H_z'*reverse(ρA_4)];

## plot
figure(figsize=[10,5])
ax1 = subplot(121)
l_cylinder_far,    = plot(r_surf,H_r_far/maximum(H_r_far), label="cylinder: far \$\\lambda_e\$", color="tab:blue");
l_cylinder_near,   = plot(r_surf,H_r_near/maximum(H_r_near), label="cylinder: near \$\\lambda_e\$", color="tab:orange");
l_finger,          = plot(r_surf,H_z/maximum(H_z), label="finger approximation \$\\lambda_e\$", color="tab:green");
s = @sprintf "planar approximation limits (far)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_far r_eal_max_far
l_cylinder_far_un  = fill_between(r_surf,H_z_far_min/maximum(H_z_far_min),H_z_far_max/maximum(H_z_far_max),alpha=0.5,color="tab:blue",label=s)
s = @sprintf "planar approximation limits (near)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_near r_eal_max_near
l_cylinder_near_un = fill_between(r_surf,H_z_near_min/maximum(H_z_near_min),H_z_near_max/maximum(H_z_near_max),alpha=0.5,color="tab:orange",label=s)
legend([(l_cylinder_far,l_cylinder_far_un),(l_cylinder_near,l_cylinder_near_un),l_finger],["cylinder: far \$\\lambda_e\$","cylinder: near \$\\lambda_e\$","finger approx."],fontsize=14)
xlabel("distance to center [\$\\mu m\$]",fontsize=14)
ylabel("normalized gain",fontsize=14)
xticks(fontsize=14)
yticks(fontsize=14)

ax2 = subplot(122)
scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_far./M_z, label="cylinder far vs planar", color="tab:blue")
scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_near./M_z, label="cylinder near vs planar", color="tab:orange")
xlabel("concentration profiles",fontsize=14)
ylabel("relative acquisition: \$\\frac{I_{\\mathrm{cylinder}}}{I_{\\mathrm{plane}}}\$",fontsize=14)
xticks(fontsize=14)
yticks(fontsize=14)
ylim(0.0)
legend(fontsize=14)

tight_layout(pad=1.0, w_pad=2.0, h_pad=0.2)
ax2.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-1.27-0.08, 0.99), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.07-0.01, 0.99), textcoords="axes fraction", color="black",fontsize=14)

# savefig("plans_vs_cylinder.png")
# savefig("plans_vs_cylinder.pdf")
