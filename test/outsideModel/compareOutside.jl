## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)

# modeling XPS
using XPSpack


##
## distances
##

μ0 = 10.0; # radius of the microjet
δr = 2.0e-3 # transition to vacuum layer thickness (let's set about 1 nm)
λe = 2.0e-3; # μ0; # EAL
L = 3μ0;

N = 51;
K = 256;
J = 257
r = collect(range(0.0,μ0+δr,length=N));
θ = collect(range(0.0,2π,length=J)) ; #collect(0:0.1:2π);
y = 0.0
Y = 0.0μ0 .+ collect(range(-0.5L,0.5L,length=K));


# far away from the analyzer
x0_far = 200.0*sqrt(2.0);
y0_far = 0.0;
z0_far = 200.0;

r_surf     = collect(range(μ0-5λe,μ0,length=N));
r_surf_out = collect(range(μ0-5λe,μ0+δr,length=N));
θ0_far     = atan(x0_far,z0_far);
θ_far      = collect(range(θ0_far-π/2.0,θ0_far+π/2.0,length=J));

H_r_far,H_rθy_far,Arn,Aθj_far,Ayk = cylinder_gain_H(r_surf,θ_far,Y,x0_far,y0_far,z0_far,μ0,λe);
H_r_far_out,H_rθy_far_out,Arn_out,Aθj_far_out,Ayk_out = cylinder_gain_H(r_surf_out,θ_far,Y,x0_far,y0_far,z0_far,μ0,λe);

# some profile
# simulate some data (one point in the kinetic energy spectrum for four different concentration profiles)
ρ0 = 1.0
ρ_vac = 0.0
r_th_out  = 2.346; # 2.0 #
r_th      = 0.346;

ρA_1 = logistic.(1000.0reverse(μ0.-r_surf).-r_th,ρ_vac,ρ0,2.0);
ρA_2 = logistic.(1000.0reverse(μ0.-r_surf).-r_th,ρ_vac,ρ0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r_surf).-1.0).^2. /(2.0*0.25^2));
ρA_3 = logistic.(1000.0reverse(μ0.-r_surf).-r_th,ρ_vac,ρ0,2.0) .+ exp.(-(1000.0reverse(μ0.-r_surf).-1.5).^2. /(2.0*0.5^2));
ρA_4 = exp.(-(1000.0reverse(μ0.-r_surf).-0.0).^2. /(2.0*0.5^2));


ρA_1_out = logistic.(1000.0reverse(δr.+μ0.-r_surf_out).-r_th_out,ρ_vac,ρ0,2.0);
ρA_2_out = logistic.(1000.0reverse(δr.+μ0.-r_surf_out).-r_th_out,ρ_vac,ρ0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r_surf_out).-0.0).^2. /(2.0*0.25^2));
ρA_3_out = logistic.(1000.0reverse(δr.+μ0.-r_surf_out).-r_th_out,ρ_vac,ρ0,2.0) .+ exp.(-(1000.0reverse(μ0.-r_surf_out).-0.5).^2. /(2.0*0.5^2));
ρA_4_out = exp.(-(1000.0(reverse(δr.+μ0.-r_surf_out).-δr)).^2. /(2.0*0.5^2));



figure(); 
plot(reverse(μ0.-r_surf),ρA_1)
plot(reverse(μ0.-r_surf),ρA_2)
plot(reverse(μ0.-r_surf),ρA_3)
plot(reverse(μ0.-r_surf),ρA_4)

figure(); 
plot(reverse(μ0.-r_surf_out),ρA_1_out)
plot(reverse(μ0.-r_surf_out),ρA_2_out)
plot(reverse(μ0.-r_surf_out),ρA_3_out)
plot(reverse(μ0.-r_surf_out),ρA_4_out)



# simulate one acquisition per profile
M_far      = [H_r_far'*reverse(ρA_1); H_r_far'*reverse(ρA_2); H_r_far'*reverse(ρA_3); H_r_far'*reverse(ρA_4)];
M_far_out  = [H_r_far_out'*reverse(ρA_1_out); H_r_far_out'*reverse(ρA_2_out); H_r_far_out'*reverse(ρA_3_out); H_r_far_out'*reverse(ρA_4_out)];



figure(figsize=[10,5])
ax1 = subplot(121)
l_cylinder_in,    = plot(r_surf,    H_r_far,     label="cylinder inside \$\\lambda_e\$", color="tab:blue");
l_cylinder_inout, = plot(r_surf_out,H_r_far_out, label="cylinder inside out outside \$\\lambda_e\$", color="tab:orange");
legend(fontsize=14)
xlabel("distance to center [\$\\mu m\$]",fontsize=14)
ylabel("normalized gain",fontsize=14)
xticks(fontsize=14)
yticks(fontsize=14)

ax2 = subplot(122)
scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_far,     label="sharp edge: inside", color="tab:blue")
scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_far_out, label="sharp edge: inside+outside", color="tab:orange")
xlabel("concentration profiles",fontsize=14)
ylabel("acquisition values",fontsize=14)
xticks(fontsize=14)
yticks(fontsize=14)
ylim(0.0)
legend(fontsize=14)

tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
ax2.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-1.27-0.08, 0.99), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.07-0.01, 0.99), textcoords="axes fraction", color="black",fontsize=14)


# savefig("compare_acquisition_in_out_sharp.pdf")
# savefig("compare_acquisition_in_out_sharp.png")

