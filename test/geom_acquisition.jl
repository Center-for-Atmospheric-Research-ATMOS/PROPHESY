# for each cases (cylindrical near and far, and finger), simulate the acquisition of one datum for each profile
M_far  = [H_r_far'*reverse(ρA_1); H_r_far'*reverse(ρA_2); H_r_far'*reverse(ρA_3); H_r_far'*reverse(ρA_4)];
M_near = [H_r_near'*reverse(ρA_1); H_r_near'*reverse(ρA_2); H_r_near'*reverse(ρA_3); H_r_near'*reverse(ρA_4)];
M_z    = [H_z'*reverse(ρA_1); H_z'*reverse(ρA_2); H_z'*reverse(ρA_3); H_z'*reverse(ρA_4)];

## plot
figure(figsize=[10,5])
ax1 = subplot(121)
plot(r_surf,H_r_far/maximum(H_r_far), label="cylinder: far \$\\lambda_e\$", color="tab:blue");
plot(r_surf,H_r_near/maximum(H_r_near), label="cylinder: near \$\\lambda_e\$", color="tab:green");
plot(r_surf,H_z/maximum(H_z), label="finger approximation \$\\lambda_e\$", color="tab:orange");
s = @sprintf "planar approximation limits (far)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_far r_eal_max_far
fill_between(r_surf,H_z_far_min/maximum(H_z_far_min),H_z_far_max/maximum(H_z_far_max),alpha=0.5,color="tab:red",label=s)
s = @sprintf "planar approximation limits (near)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_near r_eal_max_near
fill_between(r_surf,H_z_near_min/maximum(H_z_near_min),H_z_near_max/maximum(H_z_near_max),alpha=0.5,color="tab:pink",label=s)
legend()
xlabel("radius [\$\\mu m\$]")
ylabel("normalized gain")

ax2 = subplot(122)
scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_far./M_z, label="cylinder far vs planar", color="tab:blue")
scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_near./M_z, label="cylinder near vs planar", color="tab:green")
xlabel("concentration profiles")
ylabel("relative acquisition: \$\\frac{I_{\\mathrm{cylinder}}}{I_{\\mathrm{plane}}}\$")
ylim(0.0)
legend()

tight_layout(pad=1.0, w_pad=2.0, h_pad=0.2)
ax2.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-1.27, 0.99), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.07, 0.99), textcoords="axes fraction", color="black",fontsize=14)

# savefig("plans_vs_cylinder.png")
# savefig("plans_vs_cylinder.pdf")
