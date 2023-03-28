## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)

# data manipulation (loading, writing, etc)
using CSV
using DataFrames # for loading and saving the results, e.g. CSV.write("measured_psd.csv",DataFrame(PSD_meas'); writeheader=false)


# scientific package from the official Julia repositories
using LinearAlgebra

# modeling XPS
using XPSpack

# plotting background removal 
PLOT_FIG   = true;

# load data
d_c1s_Ek = CSV.File(string("../data/spectra_example.csv"); delim=",", header=true) |> DataFrame
Ny = length(d_c1s_Ek[:,1]);


##
## background removal
##

# regularization operator and parameters
D_2nd = diagm(Ny-2,Ny,1 => 2ones(Cdouble,Ny-2), 0 => -ones(Cdouble,Ny-2) ,2 => -ones(Cdouble,Ny-2)); 
κb = 0.01;
λb = 1.0e5;

# background estimation
z_baseline_1 = baseline_removal(d_c1s_Ek[:,2],λb,κb,D_2nd);
z_baseline_2 = baseline_removal(d_c1s_Ek[:,4],λb,κb,D_2nd);
z_baseline_3 = baseline_removal(d_c1s_Ek[:,6],λb,κb,D_2nd);
z_baseline_4 = baseline_removal(d_c1s_Ek[:,8],λb,κb,D_2nd);
z_baseline_5 = baseline_removal(d_c1s_Ek[:,10],λb,κb,D_2nd);


##
## cross section probability density
##

# standard deviation of the noise in the PE signal (can be estimated using SVD once the entire model is put together, but a rough approximation is enough)
σ_I1 = 100.0;
σ_I2 = 400.0;
σ_I3 = 500.0;
σ_I4 = 10.0;
σ_I5 = 20.0;

# number of iteration (main loop and inner loop of BFGS implementation)
Nbfgs = 10;
Nsearch = 10;
Nloop  = 10;

# estimate the cross section spread function and the overall integral for each spectrum)
Xend1,μ_XR1,Γ_XR1,R1,σ_R1,R_samples1 = cross_section_spread_function_sample(reverse(d_c1s_Ek[:,2])-reverse(z_baseline_1),reverse(d_c1s_Ek[:,1]),σ_I1;Nbfgs=Nbfgs,Nsearch=Nsearch,N_sample=Nloop,σ_2nd=0.001);
Xend2,μ_XR2,Γ_XR2,R2,σ_R2,R_samples2 = cross_section_spread_function_sample(reverse(d_c1s_Ek[:,4])-reverse(z_baseline_2),reverse(d_c1s_Ek[:,3]),σ_I2;Nbfgs=Nbfgs,Nsearch=Nsearch,N_sample=Nloop,σ_2nd=0.001);
Xend3,μ_XR3,Γ_XR3,R3,σ_R3,R_samples3 = cross_section_spread_function_sample(reverse(d_c1s_Ek[:,6])-reverse(z_baseline_3),reverse(d_c1s_Ek[:,5]),σ_I3;Nbfgs=Nbfgs,Nsearch=Nsearch,N_sample=Nloop,σ_2nd=0.002);
Xend4,μ_XR4,Γ_XR4,R4,σ_R4,R_samples4 = cross_section_spread_function_sample(reverse(d_c1s_Ek[:,8])-reverse(z_baseline_4),reverse(d_c1s_Ek[:,7]),σ_I4;Nbfgs=Nbfgs,Nsearch=Nsearch,N_sample=Nloop,σ_2nd=0.002);
Xend5,μ_XR5,Γ_XR5,R5,σ_R5,R_samples5 = cross_section_spread_function_sample(reverse(d_c1s_Ek[:,10])-reverse(z_baseline_5),reverse(d_c1s_Ek[:,9]),σ_I5;Nbfgs=Nbfgs,Nsearch=Nsearch,N_sample=Nloop,σ_2nd=0.005);


# plot the cross-section probability density
figure();
l_plot_1,       = plot(reverse(d_c1s_Ek[:,1]),Xend1[:,1])
l_fill_1        = fill_between(reverse(d_c1s_Ek[:,1]),Xend1[:,1]-sqrt.(diag(Γ_XR1)),Xend1[:,1]+sqrt.(diag(Γ_XR1)),alpha=0.5,color="tab:blue")
l_plot_2,       = plot(reverse(d_c1s_Ek[:,3]),Xend2[:,1])
l_fill_2        = fill_between(reverse(d_c1s_Ek[:,3]),Xend2[:,1]-sqrt.(diag(Γ_XR2)),Xend2[:,1]+sqrt.(diag(Γ_XR2)),alpha=0.5,color="tab:orange")
l_plot_3,       = plot(reverse(d_c1s_Ek[:,5]),Xend3[:,1])
l_fill_3        = fill_between(reverse(d_c1s_Ek[:,5]),Xend3[:,1]-sqrt.(diag(Γ_XR3)),Xend3[:,1]+sqrt.(diag(Γ_XR3)),alpha=0.5,color="tab:green")
l_plot_4,       = plot(reverse(d_c1s_Ek[:,7]),Xend4[:,1])
l_fill_4        = fill_between(reverse(d_c1s_Ek[:,7]),Xend4[:,1]-sqrt.(diag(Γ_XR4)),Xend4[:,1]+sqrt.(diag(Γ_XR4)),alpha=0.5,color="tab:red")
l_plot_5,       = plot(reverse(d_c1s_Ek[:,9]),Xend5[:,1])
l_fill_5        = fill_between(reverse(d_c1s_Ek[:,9]),Xend5[:,1]-sqrt.(diag(Γ_XR5)),Xend5[:,1]+sqrt.(diag(Γ_XR5)),alpha=0.5,color="tab:purple")
xlabel("binding energy [eV]",fontsize=16)
ylabel("cross section density [eV\$^{-1}\$]",fontsize=16)
xticks(fontsize=16)
yticks(fontsize=16)
xlim(288.0,295.0)
ax = gca()
ax.invert_xaxis()
legend([(l_plot_1,l_fill_1),(l_plot_2,l_fill_2),(l_plot_3,l_fill_3),(l_plot_4,l_fill_4),(l_plot_5,l_fill_5)],["\$h\\nu\$ = 350 [eV]"; "\$h\\nu\$ = 500 [eV]"; "\$h\\nu\$ = 700 [eV]"; "\$h\\nu\$ = 900 [eV]"; "\$h\\nu\$ = 1200 [eV]";],fontsize=14,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
ax.text(0.5, 0.12, "\$\\tilde{\\sigma}_{\\mathrm{C1s}}^k(K_e)\$", transform=ax.transAxes,fontsize=20)
# savefig("results/cross_section_spread_function.png")
# savefig("results/cross_section_spread_function.pdf")


# plotting the output of the background removal algorithm
if PLOT_FIG
   #  350 eV
   figure();
   scatter(d_c1s_Ek[:,1], d_c1s_Ek[:,2])
   plot(d_c1s_Ek[:,1], d_c1s_Ek[:,2]-z_baseline_1)
   plot(d_c1s_Ek[:,1], z_baseline_1)

   #  500 eV
   figure();
   scatter(d_c1s_Ek[:,3], d_c1s_Ek[:,4])
   plot(d_c1s_Ek[:,3], d_c1s_Ek[:,4]-z_baseline_2)
   plot(d_c1s_Ek[:,3], z_baseline_2)

   #  700 eV
   figure();
   scatter(d_c1s_Ek[:,5], d_c1s_Ek[:,6])
   plot(d_c1s_Ek[:,5], d_c1s_Ek[:,6]-z_baseline_3)
   plot(d_c1s_Ek[:,5], z_baseline_3)

   #  900 eV
   figure();
   scatter(d_c1s_Ek[:,7], d_c1s_Ek[:,8])
   plot(d_c1s_Ek[:,7], d_c1s_Ek[:,8]-z_baseline_4)
   plot(d_c1s_Ek[:,7], z_baseline_4)

   # 1200 eV
   figure();
   scatter(d_c1s_Ek[:,9], d_c1s_Ek[:,10])
   plot(d_c1s_Ek[:,9], d_c1s_Ek[:,10]-z_baseline_5)
   plot(d_c1s_Ek[:,9], z_baseline_5)

   # all in one figure
   x_all = [reverse(d_c1s_Ek[:,1]).+60.0; reverse(d_c1s_Ek[:,3]).+210.0; reverse(d_c1s_Ek[:,5]).+410.0; reverse(d_c1s_Ek[:,7]).+610.0; reverse(d_c1s_Ek[:,9]).+910.0];
   y_all_with_baseline = [reverse(d_c1s_Ek[:,2]); reverse(d_c1s_Ek[:,4]); reverse(d_c1s_Ek[:,6]); reverse(d_c1s_Ek[:,8]); reverse(d_c1s_Ek[:,10])];
   z_baseline = [reverse(z_baseline_1); reverse(z_baseline_2); reverse(z_baseline_3); reverse(z_baseline_4); reverse(z_baseline_5)];
   figure();
   scatter(x_all, y_all_with_baseline)
   plot(x_all,z_baseline)
   plot(x_all, y_all_with_baseline-z_baseline)
   plot(x_all, z_baseline)
   plot([x_all[1]; x_all[end]],zeros(Cdouble,2))
end