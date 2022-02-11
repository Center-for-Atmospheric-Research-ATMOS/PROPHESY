## load the packages used in the estimation
# plotting
using PyPlot
# set the font in the plots
# PyPlot.matplotlib.rcParams["font.family"]="serif"
# rc("font", family="serif")
# PyPlot.matplotlib.rcParams["mathtext.fontset"]="Computer Modern Roman"
# rc("font", serif="Computer Modern Roman")
# PyPlot.matplotlib.rcParams["text.usetex"]
# fm = matplotlib.font_manager.json_load("/path/to/fontlist-v300.json")
# fm.findfont("serif", rebuild_if_missing=False)
# fm.findfont("serif", fontext="afm", rebuild_if_missing=False)
fm = PyPlot.matplotlib.font_manager.json_load("/home/mattoz/.cache/matplotlib/fontlist-v310.json")
fm.findfont("serif", rebuild_if_missing=false)
fm.findfont("serif", fontext="afm", rebuild_if_missing=false)
rc("font",family="serif",serif="Computer Modern Roman")
rc("text", usetex=true)
using myPlot

# data manipulation (loading, writing, etc)
using CSV
using Dates
using DataFrames # for loading and saving the results, e.g. CSV.write("measured_psd.csv",DataFrame(PSD_meas'); writeheader=false)
using Printf

# scientific package from the official Julia repositories
using LinearAlgebra
using Statistics
using DSP
using SpecialMatrices
using Polynomials
using StatsBase

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack

# what to save
SAVE_MODEL = false;
SOME_INV   = false; #WARNING: not the peak area version, requires XPSinv (not yet committed)
PLOT_FIG   = false


d_c1s_Ek = CSV.File(string("data/","C1s on Ek, Na-prop and Prop.csv"); delim=",", header=true) |> DataFrame
Ny = length(d_c1s_Ek[:,1]);
D_2nd = diagm(Ny-2,Ny,1 => 2ones(Cdouble,Ny-2), 0 => -ones(Cdouble,Ny-2) ,2 => -ones(Cdouble,Ny-2)); # D2nd(Ny);
κb = 0.01;
λb = 1.0e5;
z_baseline_1 = baseline_removal(d_c1s_Ek[:,2],λb,κb,D_2nd);
if PLOT_FIG
   figure();
   scatter(d_c1s_Ek[:,1], d_c1s_Ek[:,2])
   plot(d_c1s_Ek[:,1], d_c1s_Ek[:,2]-z_baseline_1)
   plot(d_c1s_Ek[:,1], z_baseline_1)
end

# κb = 0.01;
# λb = 1.0e5; # 1.0e3
z_baseline_2 = baseline_removal(d_c1s_Ek[:,4],λb,κb,D_2nd);
if PLOT_FIG
   figure();
   scatter(d_c1s_Ek[:,3], d_c1s_Ek[:,4])
   plot(d_c1s_Ek[:,3], d_c1s_Ek[:,4]-z_baseline_2)
   plot(d_c1s_Ek[:,3], z_baseline_2)
end

# κb = 0.01;
# λb = 1.0e5;
z_baseline_3 = baseline_removal(d_c1s_Ek[:,6],λb,κb,D_2nd);
if PLOT_FIG
   figure();
   scatter(d_c1s_Ek[:,5], d_c1s_Ek[:,6])
   plot(d_c1s_Ek[:,5], d_c1s_Ek[:,6]-z_baseline_3)
   plot(d_c1s_Ek[:,5], z_baseline_3)
end

# κb = 0.01;
# λb = 1.0e5; # 1.0e4
z_baseline_4 = baseline_removal(d_c1s_Ek[:,8],λb,κb,D_2nd);
if PLOT_FIG
   figure();
   scatter(d_c1s_Ek[:,7], d_c1s_Ek[:,8])
   plot(d_c1s_Ek[:,7], d_c1s_Ek[:,8]-z_baseline_4)
   plot(d_c1s_Ek[:,7], z_baseline_4)
end

# κb = 0.01;
# λb = 1.0e5;
z_baseline_5 = baseline_removal(d_c1s_Ek[:,10],λb,κb,D_2nd);
if PLOT_FIG
   figure();
   scatter(d_c1s_Ek[:,9], d_c1s_Ek[:,10])
   plot(d_c1s_Ek[:,9], d_c1s_Ek[:,10]-z_baseline_5)
   plot(d_c1s_Ek[:,9], z_baseline_5)
end

x_all = [reverse(d_c1s_Ek[:,1]).+60.0; reverse(d_c1s_Ek[:,3]).+210.0; reverse(d_c1s_Ek[:,5]).+410.0; reverse(d_c1s_Ek[:,7]).+610.0; reverse(d_c1s_Ek[:,9]).+910.0];
y_all_with_baseline = [reverse(d_c1s_Ek[:,2]); reverse(d_c1s_Ek[:,4]); reverse(d_c1s_Ek[:,6]); reverse(d_c1s_Ek[:,8]); reverse(d_c1s_Ek[:,10])];
z_baseline = [reverse(z_baseline_1); reverse(z_baseline_2); reverse(z_baseline_3); reverse(z_baseline_4); reverse(z_baseline_5)];
y_all = y_all_with_baseline - z_baseline;
if PLOT_FIG
   figure();
   scatter(x_all, y_all_with_baseline)
   plot(x_all,z_baseline)
   plot(x_all, y_all_with_baseline-z_baseline)
   plot(x_all, z_baseline)
   plot([x_all[1]; x_all[end]],zeros(Cdouble,2))
end

function Ψ_from_spectra(Be::Array{Cdouble,2},ħν_exp::Array{Cdouble,1},μKe_exp::Array{Cdouble,1},Zi::Array{Cdouble,1},I_C1s::Array{Cdouble,2},σ_I::Array{Cdouble,1},τm::Array{Cdouble,1},μm::Array{Cdouble,1},σm::Array{Cdouble,1};Ke_min::Cdouble=288.5,Ke_max::Cdouble=294.5,Nbfgs::Int64=10,Nsearch::Int64=10,Nloop::Int64=10)
   #check dimensions
   Nke,Nbe = size(Be)
   if (Nke!=length(ħν_exp))
      throw("Model from spectra: inputs' dimensions.")
   end
   if (Nke!=length(μKe_exp))
      throw("Model from spectra: inputs' dimensions.")
   end
   if (size(Be)!=size(I_C1s))
      throw("Model from spectra: inputs' dimensions.")
   end

   # baseline removal
   D_2nd = diagm(Nbe-2,Nbe,1 => 2ones(Cdouble,Nbe-2), 0 => -ones(Cdouble,Nbe-2) ,2 => -ones(Cdouble,Nbe-2));
   κb = 0.01;
   λb = 1.0e5;
   z_baseline = zeros(Cdouble,Nke,Nbe);
   for i in 1:Nke
      z_baseline[i,:] = baseline_removal(I_C1s[i,:],λb,κb,D_2nd);
   end

   # estimate the cross section spread function and the overall integral for each spectrum)
   μ_XR = zeros(Cdouble,Nke,Nbe);
   Γ_XR = zeros(Cdouble,Nke,Nbe,Nbe);
   for i in 1:Nke
      _,μ_XR[i,:],Γ_XR[i,:,:],R[i],σ_R[i],_ = cross_section_spread_function_sample(I_C1s[i,:]-z_baseline[i,:],Ke[i,:],σ_I[i];Nbfgs=Nbfgs,Nsearch=Nsearch,N_sample=Nloop,σ_2nd=0.001);
   end

   # estimate the Np Gaussian fits for each spectrum
   Np = length(τm);
   τt = zeros(Cdouble,Nke,Np);
   μt = zeros(Cdouble,Nke,Np);
   σt = zeros(Cdouble,Nke,Np);

   for i in 1:Nke
      # estimate the peaks centers and spreads
      idx_up  = findlast(Ke[i,:].>Ke_min);
      idx_low = findfirst(Ke[i,:].<Ke_max);
      τt[i,:],μt[i,:],σt[i,:] = EM_peaks(Ke[i,idx_low:idx_up],μ_XR[i,idx_low:idx_up],τm,μm,σm,100)
   end

   # depth profile: normalized depth operator
   Fν_exp = ones(Cdouble,Nke);
   T_exp  = ones(Cdouble,Nke);
   σν_exp = ones(Cdouble,Nke,Nbe);
   wsXPS = XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σν_exp_1;α_exp=1.0);
   ## number of integration points for the Simpsons rule
   Nz0 = 50;
   ## ρ_tot_int: total concentration. It should be mainly just water concentration
   σ_z0 = 0.5; # [5 Å] width of the transition region
   z00 = 0.5;  # half height depth
   H     = Ψ_lin_peaks(Zi,wsXPS;Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=0.0,κ_eal=0.0);
   H_mean,H_std = Ψ_lin_peaks_mean_and_std(Zi,wsXPS;Nz=Nz0,κ_cs=0.0,κ_eal=0.05,σ_z=σ_z0,z0=z00);
   # TODO: sample EAL in a nice way: for each spectrum, the eal error is very much correlated, be it's not correlated between spectra
   σ_cs_orb.(ħν_exp,"C1s");
end

##
## cross section spread
##

# discretization step length for the different kinetic energy levels
dKe1 = d_c1s_Ek[1,1]-d_c1s_Ek[2,1];
dKe2 = d_c1s_Ek[1,3]-d_c1s_Ek[2,3];
dKe3 = d_c1s_Ek[1,5]-d_c1s_Ek[2,5];
dKe4 = d_c1s_Ek[1,7]-d_c1s_Ek[2,7];
dKe5 = d_c1s_Ek[1,9]-d_c1s_Ek[2,9];

# standard deviation of the noise in the PE signal (can be estimated using SVD once the entire model is put together, but a rough approximation should be enough (make sure it's not too unstable))
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
xlabel("binding energy [eV]",fontsize=12)
ylabel("cross section spread function [eV\$^{-1}\$]",fontsize=12)
legend([(l_plot_1,l_fill_1),(l_plot_2,l_fill_2),(l_plot_3,l_fill_3),(l_plot_4,l_fill_4),(l_plot_5,l_fill_5)],["\$\\hbar\\nu\$ = 60 eV"; "\$\\hbar\\nu\$ = 210 eV"; "\$\\hbar\\nu\$ = 410 eV"; "\$\\hbar\\nu\$ = 610 eV"; "\$\\hbar\\nu\$ = 910 eV";])
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
# savefig("cross_section_spread_function.png")
# savefig("cross_section_spread_function.pdf")



figure(); plot(reverse(d_c1s_Ek[:,1]),Xend1[:,1]); plot(reverse(d_c1s_Ek[:,1]),(reverse(d_c1s_Ek[:,2])-reverse(z_baseline_1))/R1)
figure(); plot(reverse(d_c1s_Ek[:,3]),Xend2[:,1]); plot(reverse(d_c1s_Ek[:,3]),(reverse(d_c1s_Ek[:,4])-reverse(z_baseline_2))/R2)
figure(); plot(reverse(d_c1s_Ek[:,5]),Xend3[:,1]); plot(reverse(d_c1s_Ek[:,5]),(reverse(d_c1s_Ek[:,6])-reverse(z_baseline_3))/R3)
figure(); plot(reverse(d_c1s_Ek[:,7]),Xend4[:,1]); plot(reverse(d_c1s_Ek[:,7]),(reverse(d_c1s_Ek[:,8])-reverse(z_baseline_4))/R4)
figure(); plot(reverse(d_c1s_Ek[:,9]),Xend5[:,1]); plot(reverse(d_c1s_Ek[:,9]),(reverse(d_c1s_Ek[:,10])-reverse(z_baseline_5))/R5)


##
## load the computed integrals
##
d_Na = CSV.File(string("data/","na_prop_prop_2_output(1).csv"); delim=",", header=true) |> DataFrame

# check the computed areas
sum(d_Na[1:3,3]) # /(R1/dKe1)
sum(d_Na[4:6,3]) # /(R2/dKe2)
sum(d_Na[7:9,3]) # /(R3/dKe3)
sum(d_Na[10:12,3]) # /(R4/dKe4)
sum(d_Na[13:15,3]) # /(R5/dKe5)

# normalized area
sum(d_Na[1:3,10])
sum(d_Na[4:6,10])
sum(d_Na[7:9,10])
sum(d_Na[10:12,10])
sum(d_Na[13:15,10])

# normalized computed values of the model
# d_Na[1:3:end,4].*d_Na[1:3:end,7].*d_Na[1:3:end,9]
(R1/dKe1)/(d_Na[1,4].*d_Na[1,7].*d_Na[1,9])
(R2/dKe2)/(d_Na[4,4].*d_Na[4,7].*d_Na[4,9])
(R3/dKe3)/(d_Na[7,4].*d_Na[7,7].*d_Na[7,9])
(R4/dKe4)/(d_Na[10,4].*d_Na[10,7].*d_Na[10,9])
(R5/dKe5)/(d_Na[13,4].*d_Na[13,7].*d_Na[13,9])

(σ_R1/dKe1)/(d_Na[1,4].*d_Na[1,7].*d_Na[1,9])
(σ_R2/dKe2)/(d_Na[4,4].*d_Na[4,7].*d_Na[4,9])
(σ_R3/dKe3)/(d_Na[7,4].*d_Na[7,7].*d_Na[7,9])
(σ_R4/dKe4)/(d_Na[10,4].*d_Na[10,7].*d_Na[10,9])
(σ_R5/dKe5)/(d_Na[13,4].*d_Na[13,7].*d_Na[13,9])

##
## fit the spectra with Gaussians
##
NN = length(d_c1s_Ek[:,1]);
μ_XR = [(R1/dKe1)*reverse(μ_XR1) (R2/dKe2)*reverse(μ_XR2) (R3/dKe3)*reverse(μ_XR3) (R4/dKe4)*reverse(μ_XR4) (R5/dKe5)*reverse(μ_XR5)];

if false
   # one Gaussian per peak: does not seem to fit well enough... or maybe it does... the residus are very similar in both cases
   τm = [1.0/3.0; 1.0/3.0; 1.0/3.0];
   μm = [290.3; 291.9; 293.5];
   σm = [290.3; 291.9; 293.5]/500.0;

   τt = zeros(Cdouble,5,3);
   μt = zeros(Cdouble,5,3);
   σt = zeros(Cdouble,5,3);

   for i in 1:5
      # estimate the peaks centers and spreads
      idx_up  = findlast(d_c1s_Ek[:,(2i)-1].>288.5);
      idx_low = findfirst(d_c1s_Ek[:,(2i)-1].<294.5);
      τt[i,:],μt[i,:],σt[i,:] = EM_peaks(d_c1s_Ek[idx_low:idx_up,(2i)-1],μ_XR[idx_low:idx_up,i],τm,μm,σm,100) # too many near zero values are skewing the estimation here: TODO: change the support of μ_XRi
      # τt[i,:],μt[i,:],σt[i,:] = EM_peaks(d_c1s_Ek[μ_XR[:,i].>0.05maximum(μ_XR[:,i]),(2i)-1],μ_XR[μ_XR[:,i].>0.05maximum(μ_XR[:,i]),i],τm,μm,σm,100)
      # τt[i,:],μt[i,:],σt[i,:] = EM_peaks(d_c1s_Ek[:,(2i)-1],d_c1s_Ek[:,2i]-reverse(z_baseline[(i-1)*NN+1:i*NN]),τm,μm,σm,100)
      # plt
      x_dist_1 = (τt[i,1]/σt[i,1])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,1])/σt[i,1]).^2);
      x_dist_2 = (τt[i,2]/σt[i,2])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,2])/σt[i,2]).^2);
      x_dist_3 = (τt[i,3]/σt[i,3])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,3])/σt[i,3]).^2);
      if PLOT_FIG
         figure();
         plot(d_c1s_Ek[:,(2i)-1],x_dist_1)
         plot(d_c1s_Ek[:,(2i)-1],x_dist_2)
         plot(d_c1s_Ek[:,(2i)-1],x_dist_3)
         plot(d_c1s_Ek[:,(2i)-1],x_dist_1+x_dist_2+x_dist_3)
         max_xx = maximum(x_dist_1+x_dist_2+x_dist_3);
         # max_dd = maximum(d_c1s_Ek[:,2i]-reverse(z_baseline[(i-1)*NN+1:i*NN]));
         # plot(d_c1s_Ek[:,(2i)-1],(max_xx/max_dd)*(d_c1s_Ek[:,2i]-reverse(z_baseline[(i-1)*NN+1:i*NN])))
         max_dd = maximum(μ_XR[:,i]);
         # plot(d_c1s_Ek[:,(2i)-1],(max_xx/max_dd)*μ_XR[:,i]-(x_dist_1+x_dist_2+x_dist_3))
         plot(d_c1s_Ek[:,(2i)-1],(max_xx/max_dd)*μ_XR[:,i])
      end
   end
else
   # two gaussian per peak
   τm = [1.0/6.0; 1.0/6.0; 1.0/6.0; 1.0/6.0; 1.0/6.0; 1.0/6.0];
   μm = [290.3-0.1; 290.3+0.1; 291.9-0.2; 291.9+0.2; 293.3; 294.1];
   σm = [290.3/1000.0; 290.3/1000.0; 291.9/1000.0; 291.9/1000.0; 293.5/1000.0; 293.5/1000.0];

   τt = zeros(Cdouble,5,6);
   μt = zeros(Cdouble,5,6);
   σt = zeros(Cdouble,5,6);


   for i in 1:5
      # estimate the peaks centers and spreads
      idx_up  = findlast(d_c1s_Ek[:,(2i)-1].>288.5);
      idx_low = findfirst(d_c1s_Ek[:,(2i)-1].<294.5);
      τt[i,:],μt[i,:],σt[i,:] = EM_peaks(d_c1s_Ek[idx_low:idx_up,(2i)-1],μ_XR[idx_low:idx_up,i],τm,μm,σm,100) # too many near zero values are skewing the estimation here: TODO: change the support of μ_XRi
      # τt[i,:],μt[i,:],σt[i,:] = EM_peaks(d_c1s_Ek[μ_XR[:,i].>0.05maximum(μ_XR[:,i]),(2i)-1],μ_XR[μ_XR[:,i].>0.05maximum(μ_XR[:,i]),i],τm,μm,σm,100)
      # τt[i,:],μt[i,:],σt[i,:] = EM_peaks(d_c1s_Ek[:,(2i)-1],d_c1s_Ek[:,2i]-reverse(z_baseline[(i-1)*NN+1:i*NN]),τm,μm,σm,100)
      # plt
      x_dist_1 = (τt[i,1]/σt[i,1])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,1])/σt[i,1]).^2);
      x_dist_2 = (τt[i,2]/σt[i,2])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,2])/σt[i,2]).^2);
      x_dist_3 = (τt[i,3]/σt[i,3])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,3])/σt[i,3]).^2);
      x_dist_4 = (τt[i,4]/σt[i,4])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,4])/σt[i,4]).^2);
      x_dist_5 = (τt[i,5]/σt[i,5])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,5])/σt[i,5]).^2);
      x_dist_6 = (τt[i,6]/σt[i,6])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,6])/σt[i,6]).^2);
      if PLOT_FIG
         figure();
         plot(d_c1s_Ek[:,(2i)-1],x_dist_1)
         plot(d_c1s_Ek[:,(2i)-1],x_dist_2)
         plot(d_c1s_Ek[:,(2i)-1],x_dist_3)
         plot(d_c1s_Ek[:,(2i)-1],x_dist_4)
         plot(d_c1s_Ek[:,(2i)-1],x_dist_5)
         plot(d_c1s_Ek[:,(2i)-1],x_dist_6)
         plot(d_c1s_Ek[:,(2i)-1],x_dist_1+x_dist_2+x_dist_3+x_dist_4+x_dist_5+x_dist_6)
         max_xx = maximum(x_dist_1+x_dist_2+x_dist_3+x_dist_4+x_dist_5+x_dist_6);
         # max_dd = maximum(d_c1s_Ek[:,2i]-reverse(z_baseline[(i-1)*NN+1:i*NN]));
         # plot(d_c1s_Ek[:,(2i)-1],(max_xx/max_dd)*(d_c1s_Ek[:,2i]-reverse(z_baseline[(i-1)*NN+1:i*NN])))
         max_dd = maximum(μ_XR[:,i]);
         plot(d_c1s_Ek[:,(2i)-1],(max_xx/max_dd)*μ_XR[:,i])
         # plot(d_c1s_Ek[:,(2i)-1],(max_xx/max_dd)*μ_XR[:,i]-(x_dist_1+x_dist_2+x_dist_3+x_dist_4+x_dist_5+x_dist_6))
      end
   end
end

## number of integration points for the Simpsons rule
Nz0 = 50;
## ρ_tot_int: total concentration. It should be mainly just water concentration
σ_z0 = 0.5; # [5 Å] width of the transition region
z00 = 0.5;  # half height depth

ħν_exp = 1.0*d_Na[1:3:end,5];
Fν_exp = 1.0*d_Na[1:3:end,4]; # 1.0e-9
#WARNING: rectifying mirror current, maybe the gain to transform the current to photon flux is not the same for all photon energy!!!!!
Fν_exp[4] = Fν_exp[4]/20.0;
Fν_exp[5] = Fν_exp[5]/45.0;
T_exp  = d_Na[1:3:end,9];
μKe_exp = d_Na[1:3:end,8];
Be_exp = [d_c1s_Ek[:,1]'; d_c1s_Ek[:,3]'; d_c1s_Ek[:,5]'; d_c1s_Ek[:,7]'; d_c1s_Ek[:,9]'];
σν_exp_1 = zeros(Cdouble,5,length(d_c1s_Ek[:,1]));
σν_exp_2 = zeros(Cdouble,5,length(d_c1s_Ek[:,1]));
σν_exp_3 = zeros(Cdouble,5,length(d_c1s_Ek[:,1]));


##
##TODO: use one peak by one peak to create models and then data. Later, create bad models in the last two channels
##
for i in 1:5
   x_dist_1 = (τt[i,1]/σt[i,1])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,1])/σt[i,1]).^2);
   x_dist_2 = (τt[i,2]/σt[i,2])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,2])/σt[i,2]).^2);
   x_dist_3 = (τt[i,3]/σt[i,3])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,3])/σt[i,3]).^2);
   max_xx = maximum(x_dist_1+x_dist_2+x_dist_3);
   σν_exp_1[i,:] = (σ_cs_orb(ħν_exp[i],"C1s")/max_xx)*x_dist_1
   σν_exp_2[i,:] = (σ_cs_orb(ħν_exp[i],"C1s")/max_xx)*x_dist_2
   σν_exp_3[i,:] = (σ_cs_orb(ħν_exp[i],"C1s")/max_xx)*x_dist_3
end

Fν_exp   .= 1.0;
T_exp    .= 1.0;
σν_exp_1 = σν_exp_1./σ_cs_orb.(ħν_exp,"C1s");
σν_exp_2 = σν_exp_2./σ_cs_orb.(ħν_exp,"C1s");
σν_exp_3 = σν_exp_3./σ_cs_orb.(ħν_exp,"C1s");

wsXPS_1 = XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σν_exp_1;α_exp=1.0);
wsXPS_2 = XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σν_exp_2;α_exp=1.0);
wsXPS_3 = XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σν_exp_3;α_exp=1.0);


# depth discretization
N = 50;
Z_max = 10.0;
Zi = collect(range(0.0,Z_max,length=N));
H_1 = Ψ_lin_peaks(Zi,wsXPS_1;Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=0.0,κ_eal=0.0);
H_2 = Ψ_lin_peaks(Zi,wsXPS_2;Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=0.0,κ_eal=0.0);
H_3 = Ψ_lin_peaks(Zi,wsXPS_3;Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=0.0,κ_eal=0.0);
# figure(); imshow(H_1); colorbar()
# figure(); imshow(H_2); colorbar()
# figure(); imshow(H_3); colorbar()
H_mean_1,H_std_1 = Ψ_lin_peaks_mean_and_std(Zi,wsXPS_1;Nz=Nz0,κ_cs=0.05,κ_eal=0.05,σ_z=σ_z0,z0=z00);
H_mean_2,H_std_2 = Ψ_lin_peaks_mean_and_std(Zi,wsXPS_2;Nz=Nz0,κ_cs=0.05,κ_eal=0.05,σ_z=σ_z0,z0=z00);
H_mean_3,H_std_3 = Ψ_lin_peaks_mean_and_std(Zi,wsXPS_3;Nz=Nz0,κ_cs=0.05,κ_eal=0.05,σ_z=σ_z0,z0=z00);

##
## peak area model
##

A1 = Matrix{Cdouble}([dKe1*dropdims(sum(H_1[1:wsXPS_1.Nbe,:],dims=1),dims=1) dKe2*dropdims(sum(H_1[wsXPS_1.Nbe+1:2wsXPS_1.Nbe,:],dims=1),dims=1) dKe3*dropdims(sum(H_1[2wsXPS_1.Nbe+1:3wsXPS_1.Nbe,:],dims=1),dims=1) dKe4*dropdims(sum(H_1[3wsXPS_1.Nbe+1:4wsXPS_1.Nbe,:],dims=1),dims=1) dKe5*dropdims(sum(H_1[4wsXPS_1.Nbe+1:5wsXPS_1.Nbe,:],dims=1),dims=1)]');
A2 = Matrix{Cdouble}([dKe1*dropdims(sum(H_2[1:wsXPS_1.Nbe,:],dims=1),dims=1) dKe2*dropdims(sum(H_2[wsXPS_1.Nbe+1:2wsXPS_1.Nbe,:],dims=1),dims=1) dKe3*dropdims(sum(H_2[2wsXPS_1.Nbe+1:3wsXPS_1.Nbe,:],dims=1),dims=1) dKe4*dropdims(sum(H_2[3wsXPS_1.Nbe+1:4wsXPS_1.Nbe,:],dims=1),dims=1) dKe5*dropdims(sum(H_2[4wsXPS_1.Nbe+1:5wsXPS_1.Nbe,:],dims=1),dims=1)]');
A3 = Matrix{Cdouble}([dKe1*dropdims(sum(H_3[1:wsXPS_1.Nbe,:],dims=1),dims=1) dKe2*dropdims(sum(H_3[wsXPS_1.Nbe+1:2wsXPS_1.Nbe,:],dims=1),dims=1) dKe3*dropdims(sum(H_3[2wsXPS_1.Nbe+1:3wsXPS_1.Nbe,:],dims=1),dims=1) dKe4*dropdims(sum(H_3[3wsXPS_1.Nbe+1:4wsXPS_1.Nbe,:],dims=1),dims=1) dKe5*dropdims(sum(H_3[4wsXPS_1.Nbe+1:5wsXPS_1.Nbe,:],dims=1),dims=1)]');

A_std_1 = Matrix{Cdouble}([dKe1*sqrt.(dropdims(sum(H_std_1[1:wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe2*sqrt.(dropdims(sum(H_std_1[wsXPS_1.Nbe+1:2wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe3*sqrt.(dropdims(sum(H_std_1[2wsXPS_1.Nbe+1:3wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe4*sqrt.(dropdims(sum(H_std_1[3wsXPS_1.Nbe+1:4wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe5*sqrt.(dropdims(sum(H_std_1[4wsXPS_1.Nbe+1:5wsXPS_1.Nbe,:].^2,dims=1),dims=1))]');
A_std_2 = Matrix{Cdouble}([dKe1*sqrt.(dropdims(sum(H_std_2[1:wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe2*sqrt.(dropdims(sum(H_std_2[wsXPS_1.Nbe+1:2wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe3*sqrt.(dropdims(sum(H_std_2[2wsXPS_1.Nbe+1:3wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe4*sqrt.(dropdims(sum(H_std_2[3wsXPS_1.Nbe+1:4wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe5*sqrt.(dropdims(sum(H_std_2[4wsXPS_1.Nbe+1:5wsXPS_1.Nbe,:].^2,dims=1),dims=1))]');
A_std_3 = Matrix{Cdouble}([dKe1*sqrt.(dropdims(sum(H_std_3[1:wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe2*sqrt.(dropdims(sum(H_std_3[wsXPS_1.Nbe+1:2wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe3*sqrt.(dropdims(sum(H_std_3[2wsXPS_1.Nbe+1:3wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe4*sqrt.(dropdims(sum(H_std_3[3wsXPS_1.Nbe+1:4wsXPS_1.Nbe,:].^2,dims=1),dims=1)) dKe5*sqrt.(dropdims(sum(H_std_3[4wsXPS_1.Nbe+1:5wsXPS_1.Nbe,:].^2,dims=1),dims=1))]');

A = A1+A2+A3;
A_std = sqrt.(A_std_1.^2 + A_std_2.^2 + A_std_3.^2);

figure()
for i in 1:5
   plot(Zi,A1[i,:])
   fill_between(Zi,A1[i,:]-A_std_1[i,:],A1[i,:]+A_std_1[i,:],alpha=0.5)
end

figure()
for i in 1:5
   plot(Zi,A2[i,:])
   fill_between(Zi,A2[i,:]-A_std_2[i,:],A2[i,:]+A_std_2[i,:],alpha=0.5)
end

figure()
for i in 1:5
   plot(Zi,A3[i,:])
   fill_between(Zi,A3[i,:]-A_std_3[i,:],A3[i,:]+A_std_3[i,:],alpha=0.5)
end

# normalized area operator
figure()
for i in 1:5
   plot(Zi,A[i,:])
   fill_between(Zi,A[i,:]-A_std[i,:],A[i,:]+A_std[i,:],alpha=0.5)
end

# non-normalized area operator
figure()
for i in 1:5
   plot(Zi,d_Na[3*(i-1)+1,4]*d_Na[3*(i-1)+1,9]*σ_cs_orb(ħν_exp[i],"C1s")*A[i,:])
   fill_between(Zi,d_Na[3*(i-1)+1,4]*d_Na[3*(i-1)+1,9]*σ_cs_orb(ħν_exp[i],"C1s")*(A[i,:]-A_std[i,:]),d_Na[3*(i-1)+1,4]*d_Na[3*(i-1)+1,9]*σ_cs_orb(ħν_exp[i],"C1s")*(A[i,:]+A_std[i,:]),alpha=0.5)
end


##
## profiles
##

ρA_1 = logistic.(Zi.-2.0,0.0,1.0,2.0);
ρA_2 = logistic.(Zi.-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(Zi.-1.0).^2. /(2.0*0.25^2));
ρA_3 = logistic.(Zi.-2.0,0.0,1.0,2.0) .+ exp.(-(Zi.-1.5).^2. /(2.0*0.5^2));
ρA_4 = exp.(-(Zi.-2.5).^2. /(2.0*0.5^2));


##
## generate spectra
##
IA_1_clean_1 = H_1*ρA_1;
IA_2_clean_1 = H_1*ρA_2;
IA_3_clean_1 = H_1*ρA_3;
IA_4_clean_1 = H_1*ρA_4;

R1s = A1*ρA_1
r1s = R1s[2:end]./R1s[1:end-1];
RR1s = r1s.*A1[2:end,:]-A1[1:end-1,:];

σ_noise = 0.0
IA_1_1 = IA_1_clean_1 + σ_noise*randn(wsXPS_1.Nke*wsXPS_1.Nbe);
IA_2_1 = IA_2_clean_1 + σ_noise*randn(wsXPS_1.Nke*wsXPS_1.Nbe);
IA_3_1 = IA_3_clean_1 + σ_noise*randn(wsXPS_1.Nke*wsXPS_1.Nbe);
IA_4_1 = IA_4_clean_1 + σ_noise*randn(wsXPS_1.Nke*wsXPS_1.Nbe);

figure();
for i in 1:wsXPS_1.Nν
   scatter(wsXPS_1.Be[i,:],IA_1_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
   plot(wsXPS_1.Be[i,:],IA_1_clean_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
end
xlabel("B\$_e\$ [a.u.]")
ylabel("PE signal [a.u.]")
xlim(minimum(wsXPS_1.Be),maximum(wsXPS_1.Be))
title("simulated spectra")
ylim(0.0,1.1maximum(IA_1_1))


figure();
for i in 1:wsXPS_1.Nν
   scatter(wsXPS_1.Be[i,:],IA_2_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
   plot(wsXPS_1.Be[i,:],IA_2_clean_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
end
xlabel("B\$_e\$ [a.u.]")
ylabel("PE signal [a.u.]")
xlim(minimum(wsXPS_1.Be),maximum(wsXPS_1.Be))
title("simulated spectra")
ylim(0.0,1.1maximum(IA_2_1))


figure();
for i in 1:wsXPS_1.Nν
   scatter(wsXPS_1.Be[i,:],IA_3_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
   plot(wsXPS_1.Be[i,:],IA_3_clean_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
end
xlabel("B\$_e\$ [a.u.]")
ylabel("PE signal [a.u.]")
xlim(minimum(wsXPS_1.Be),maximum(wsXPS_1.Be))
title("simulated spectra")
ylim(0.0,1.1maximum(IA_3_1))


figure();
for i in 1:wsXPS_1.Nν
   scatter(wsXPS_1.Be[i,:],IA_4_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
   plot(wsXPS_1.Be[i,:],IA_4_clean_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
end
xlabel("B\$_e\$ [a.u.]")
ylabel("PE signal [a.u.]")
xlim(minimum(wsXPS_1.Be),maximum(wsXPS_1.Be))
title("simulated spectra")
ylim(0.0,1.1maximum(IA_4_1))
# end

if SAVE_MODEL
   # sample the relative error in cross section and eal to generate measurement model with modelling error
   model_folder = "./data/lin/real/low_res_50/more/";
   κ_css  = collect(range(-0.05,0.05,length=11));
   κ_eals = collect(range(-0.05,0.05,length=11));
   i = 0
   t_elapsed = @elapsed for (κ1,κ2) in Iterators.product(κ_css,κ_eals)
      global i = i + 1;
      # compute the measurment model with some error
      # H,λeal,σ_cs = Ψ_lin_peaks(Zi,Kes,μ_σ_ke;Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=κ1,κ_eal=κ2);
      Hκ = Ψ_lin_peaks(Zi,wsXPS;Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=κ1,κ_eal=κ2);
      # save the model, the intended and true crosss section and eal
      s = @sprintf "%s/%i/" model_folder i
      mkpath(s)
      CSV.write(string(s,"kinetic_energy.csv"),DataFrame(wsXPS.Ke',:auto); header=false);
      CSV.write(string(s,"depth.csv"),DataFrame(Zi',:auto); header=false);
      CSV.write(string(s,"H.csv"),DataFrame(Hκ,:auto); header=false);
   end
end

if false
   # rescale H
   # y_all = [d_c1s_Ek[:,2]; d_c1s_Ek[:,4]; d_c1s_Ek[:,6]; d_c1s_Ek[:,8];d_c1s_Ek[:,10]]
   y_all = [d_c1s_Ek[:,2].-z_baseline_1; d_c1s_Ek[:,4].-z_baseline_2; d_c1s_Ek[:,6].-z_baseline_3; d_c1s_Ek[:,8].-z_baseline_4;d_c1s_Ek[:,10].-z_baseline_5]
   H = (maximum(y_all[4wsXPS.Nbe+1:5wsXPS.Nbe])/maximum(H[4wsXPS.Nbe+1:5wsXPS.Nbe,2]))*H

   F = svd(H, full=true);
   F1 = svd(H[1:wsXPS.Nbe,:], full=true);
   F2 = svd(H[wsXPS.Nbe+1:2wsXPS.Nbe,:], full=true);
   F3 = svd(H[2wsXPS.Nbe+1:3wsXPS.Nbe,:], full=true);
   F4 = svd(H[3wsXPS.Nbe+1:4wsXPS.Nbe,:], full=true);
   F5 = svd(H[4wsXPS.Nbe+1:5wsXPS.Nbe,:], full=true);

   abs_ub_1            = abs.(F.U'*IA_1);
   abs_ub_2            = abs.(F.U'*IA_2);
   abs_ub_3            = abs.(F.U'*IA_3);
   min_ub,max_ub = extrema([abs_ub_1; abs_ub_2; abs_ub_3])
   l_scat_4_ub_1        = scatter(collect(1:wsXPS.Nbe*wsXPS.Nν),abs.(abs_ub_1),color="tab:blue")
   l_scat_4_ub_2        = scatter(collect(1:wsXPS.Nbe*wsXPS.Nν),abs.(abs_ub_2),color="tab:orange")
   l_scat_4_ub_3        = scatter(collect(1:wsXPS.Nbe*wsXPS.Nν),abs.(abs_ub_3),color="tab:pink")
   ylim(0.5min_ub,5.0max_ub)
   yscale("log")
   xlabel("SVD left index",fontsize=12)
   ylabel("SVD left coefficient",fontsize=12)
   legend([l_scat_4_ub_1,l_scat_4_ub_2,l_scat_4_ub_3],["data 1","data 2","darta 3"])


   figure()
   abs_ub            = abs.(F.U'*y_all);
   min_ub,max_ub = extrema(abs_ub)
   l_scat_4_ub        = scatter(collect(1:wsXPS.Nbe*wsXPS.Nν),abs.(abs_ub),color="tab:blue")
   ylim(0.5min_ub,5.0max_ub)
   yscale("log")
   xlabel("SVD left index",fontsize=12)
   ylabel("SVD left coefficient",fontsize=12)
   legend([l_scat_4_ub],["data"])

   figure()
   abs_ub_1 = abs.(F1.U'*(d_c1s_Ek[:,2]-z_baseline_1))
   abs_ub_2 = abs.(F2.U'*(d_c1s_Ek[:,4]-z_baseline_2))
   abs_ub_3 = abs.(F3.U'*(d_c1s_Ek[:,6]-z_baseline_3))
   abs_ub_4 = abs.(F4.U'*(d_c1s_Ek[:,8]-z_baseline_4))
   abs_ub_5 = abs.(F5.U'*(d_c1s_Ek[:,10]-z_baseline_5))
   σ_1 = mean(abs_ub_1[50:end])
   σ_2 = mean(abs_ub_2[50:end])
   σ_3 = mean(abs_ub_3[50:end])
   σ_4 = mean(abs_ub_4[50:end])
   σ_5 = mean(abs_ub_5[50:end])
   min_ub,max_ub = extrema([abs_ub_1; abs_ub_2; abs_ub_3; abs_ub_4; abs_ub_5;])
   l_scat_ub_1        = scatter(collect(1:wsXPS.Nbe),abs.(abs_ub_1),color="tab:blue")
   plot([1;wsXPS.Nbe], [σ_1;σ_1],color="tab:blue")
   l_scat_ub_2        = scatter(collect(wsXPS.Nbe+1:2wsXPS.Nbe),abs.(abs_ub_2),color="tab:orange")
   plot([wsXPS.Nbe+1;2wsXPS.Nbe], [σ_2;σ_2],color="tab:orange")
   l_scat_ub_3        = scatter(collect(2wsXPS.Nbe+1:3wsXPS.Nbe),abs.(abs_ub_3),color="tab:green")
   plot([2wsXPS.Nbe+1;3wsXPS.Nbe], [σ_3;σ_3],color="tab:green")
   l_scat_ub_4        = scatter(collect(3wsXPS.Nbe+1:4wsXPS.Nbe),abs.(abs_ub_4),color="tab:red")
   plot([3wsXPS.Nbe+1;4wsXPS.Nbe], [σ_4;σ_4],color="tab:red")
   l_scat_ub_5        = scatter(collect(4wsXPS.Nbe+1:5wsXPS.Nbe),abs.(abs_ub_5),color="tab:cyan")
   plot([4wsXPS.Nbe+1;5wsXPS.Nbe], [σ_5;σ_5],color="tab:cyan")
   ylim(0.5min_ub,5.0max_ub)
   yscale("log")
   xlabel("SVD left index",fontsize=12)
   ylabel("SVD left coefficient",fontsize=12)
   legend([l_scat_ub_1,l_scat_ub_2,l_scat_ub_3,l_scat_ub_4,l_scat_ub_5],["data 1"; "data 2"; "data 3"; "data 4"; "data 5"])






   figure()
   l_scat_1        = scatter(collect(1:wsXPS.Nbe),d_c1s_Ek[:,2]-z_baseline_1,color="tab:blue")
   l_scat_2        = scatter(collect(wsXPS.Nbe+1:2wsXPS.Nbe),d_c1s_Ek[:,4]-z_baseline_2,color="tab:orange")
   l_scat_3        = scatter(collect(2wsXPS.Nbe+1:3wsXPS.Nbe),d_c1s_Ek[:,6]-z_baseline_3,color="tab:green")
   l_scat_4        = scatter(collect(3wsXPS.Nbe+1:4wsXPS.Nbe),d_c1s_Ek[:,8]-z_baseline_4,color="tab:red")
   l_scat_5        = scatter(collect(4wsXPS.Nbe+1:5wsXPS.Nbe),d_c1s_Ek[:,10]-z_baseline_5,color="tab:cyan")
   xlabel("SVD left index",fontsize=12)
   ylabel("PE signal [a.u.]",fontsize=12)
   legend([l_scat_1,l_scat_2,l_scat_3,l_scat_4,l_scat_5],["data 1"; "data 2"; "data 3"; "data 4"; "data 5"])



   figure()
   plot(wsXPS.Ke[4wsXPS.Nbe+1:5wsXPS.Nbe],y_all[4wsXPS.Nbe+1:5wsXPS.Nbe])
   plot(wsXPS.Ke[4wsXPS.Nbe+1:5wsXPS.Nbe],(maximum(y_all[4wsXPS.Nbe+1:5wsXPS.Nbe])/maximum(H[4wsXPS.Nbe+1:5wsXPS.Nbe,2]))*H[4wsXPS.Nbe+1:5wsXPS.Nbe,2])


   figure()
   plot(wsXPS.Ke[wsXPS.Nbe+1:2wsXPS.Nbe],y_all[wsXPS.Nbe+1:2wsXPS.Nbe])
   plot(wsXPS.Ke[wsXPS.Nbe+1:2wsXPS.Nbe],(maximum(y_all[wsXPS.Nbe+1:2wsXPS.Nbe])/maximum(H[wsXPS.Nbe+1:2wsXPS.Nbe,2]))*H[wsXPS.Nbe+1:2wsXPS.Nbe,2])


   minH,maxH = extrema(H)
   fig1, ax1, pcm1 = imshowData(1254,Zi,wsXPS.Ke[wsXPS.Nbe+1:2wsXPS.Nbe],H[wsXPS.Nbe+1:2wsXPS.Nbe,:],_norm=:Normalize,_vmin=minH,_vmax=maxH,_edgecolors="face",_sub=221)
   xlabel("depth \$z\$",fontsize=14) # [nm]
   ylabel("Kinetic energy \$K_e\$",fontsize=14)
   xlim(Zi[1],Zi[end])
   ylim(wsXPS.Ke[wsXPS.Nbe+1],wsXPS.Ke[2wsXPS.Nbe])
   rc("ytick",color="white")
   cax1 = fig1.add_axes([0.08+0.33, .75, 0.02, 0.2])
   cb1 = fig1.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
   cb1.set_label("measurement operator [a.u.]", color="white", fontsize=10)
   cb1.ax.yaxis.set_tick_params(color="white")
   cb1.outline.set_edgecolor("white")
   rc("ytick",color="black")


   minU,maxU = extrema(F.U)
   minU,maxU = extrema(F.U[wsXPS.Nbe+1:2wsXPS.Nbe,:])
   fig1, ax2, pcm2 = imshowData(1254,collect(1.0:wsXPS.Nbe*wsXPS.Nν),wsXPS.Ke[wsXPS.Nbe+1:2wsXPS.Nbe],F.U[wsXPS.Nbe+1:2wsXPS.Nbe,:],_norm=:Normalize,_vmin=minU,_vmax=maxU,_edgecolors="face",_sub=222)
   xlabel("basis index [data space]",fontsize=14) # [nm]
   ylabel("Kinetic energy \$K_e\$",fontsize=14)
   xlim(1,length(wsXPS.Ke))
   ylim(wsXPS.Ke[wsXPS.Nbe+1],wsXPS.Ke[2wsXPS.Nbe])
   rc("ytick",color="white")
   cax2 = fig1.add_axes([0.57+0.31, .77, 0.02, 0.2])
   cb2 = fig1.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
   cb2.set_label("left vector entries [a.u.]", color="white", fontsize=10)
   cb2.ax.yaxis.set_tick_params(color="white")
   cb2.outline.set_edgecolor("white")
   rc("ytick",color="black")


   minV,maxV = extrema(F.Vt)
   fig1, ax3, pcm3 = imshowData(1254,Zi,collect(1.0:N),F.Vt,_norm=:Normalize,_vmin=minV,_vmax=maxV,_edgecolors="face",_sub=223)
   xlabel("depth \$z\$",fontsize=14) # [nm]
   ylabel("basis index [depth space]",fontsize=14)
   xlim(Zi[1],Zi[end])
   ylim(1,N)
   rc("ytick",color="white")
   cax3 = fig1.add_axes([0.08, .27, 0.02, 0.2])
   cb3 = fig1.colorbar(pcm3, orientation="vertical", cax=cax3, shrink=0.6)
   cb3.set_label("right vector entries [a.u.]", color="white", fontsize=10)
   cb3.ax.yaxis.set_tick_params(color="white")
   cb3.outline.set_edgecolor("white")
   rc("ytick",color="black")


   ax4 = subplot(224)
   semilogy(F.S)
   xlabel("basis index [depth space]",fontsize=14)
   ylabel("singular values [a.u.]",fontsize=14)


   fig1.set_figwidth(10.0) # (5.6)
   fig1.set_figheight(9.0)
   tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)

   ax3.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.12, 2.12), textcoords="axes fraction", color="black",fontsize=14)
   ax3.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(1.07, 2.12), textcoords="axes fraction", color="black",fontsize=14)
   ax3.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(-0.12, 1.0), textcoords="axes fraction", color="black",fontsize=14)
   ax3.annotate("d)", xy=(3, 1),  xycoords="data", xytext=(1.07, 1.0), textcoords="axes fraction", color="black",fontsize=14)
end







##
## try some inversion
##
if SOME_INV
   using XPSinv


   D_2nd = D2nd(N);
   d_2nd = 2;
   B = zeros(Cdouble,2,N); B[1,1] = 1.0; B[2,end] = 1.0
   γ = σ_1^2*ones(Cdouble,wsXPS.Nbe*wsXPS.Nke);
   γ[wsXPS.Nbe+1:2wsXPS.Nbe] .= σ_2^2;
   γ[2wsXPS.Nbe+1:3wsXPS.Nbe] .= σ_3^2;
   γ[3wsXPS.Nbe+1:4wsXPS.Nbe] .= σ_4^2;
   γ[4wsXPS.Nbe+1:5wsXPS.Nbe] .= σ_5^2;
   Γ_inv = diagm(1.0./γ);
   # z00  = 2.0;
   σ_z0 = 1.0;

   γ_D = logistic.(Zi.-3.0,0.2,0.0,2.0σ_z0)[2:end-1];
   # γ_D = 0.1logistic.(Zi.-3.0,0.2,0.0,2.0σ_z0)[2:end-1].+10.0;
   γ_D = (γ_D .+ 0.01maximum(γ_D)).^2;
   Γ_D_inv = diagm(1.0./γ_D);

   ρ0 = ρA_1[1] # 0.0;
   γ0 = 0.01^2 # 1.0^2 # 0.01^2 # 0.1^2
   ρB = ρA_1[end] # 1.0;
   γB = 0.01^2 # 0.01^2 # 0.1^2
   ρ_bar = [ρ0;ρB];
   γ_ρ = [γ0;γB];
   Γ_ρ_inv = diagm(1.0./γ_ρ);
   rankK = rank(H);



   ##
   ## draw samples from data distribution
   ##
   # μ_ρ,σ_ρ,ρ_samples,F,W_inv = data_sample_nso_un(1000,rankK,H,0.0H,σ_noise*ones(Cdouble,Nke),Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,IA_1_clean,ρ_bar)

   # ρ_clean,_,_ = iterative_nso(rankK,H,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,y_all,ρ_bar);
   # ρ_clean,_,_ = iterative_nso(rankK,H,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,H*ρA_1,ρ_bar);


   # ρ_clean,σ_ρ,ρ_samples,F,W_inv = data_sample_nso_un(100,rankK,H,0.01H,sqrt.(γ),Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,H*ρA_1,ρ_bar)
   Nmin = 50; # 20
   ρ_clean,σ_ρ,ρ_samples,F,W_inv = shuffle_data_sample_nso_un(100,wsXPS.Nke,rankK,H,1.0e-2H,γ,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,H*ρA_1,ρ_bar;Nmin=Nmin);


   figure(); plot(Zi,ρ_clean)
   figure(); plot(Zi,ρ_samples')

   # ρ_clean .= 1.0
   # ρ_clean[:] = ρA_3/1.5
   #
   # ρ_clean[:] = 1.2ρA_1
   #
   #
   # ρ_clean[:] = ρA_3-0.9ρA_1

   figure()
   l_scat_1        = scatter(collect(1:wsXPS.Nbe),d_c1s_Ek[:,2]-z_baseline_1,color="tab:blue")
   l_plot_1,       = plot(collect(1:wsXPS.Nbe),H[1:wsXPS.Nbe,:]*ρ_clean,color="tab:blue")
   l_fill_1        = fill_between(collect(1:wsXPS.Nbe),H[1:wsXPS.Nbe,:]*ρ_clean-σ_1*ones(Cdouble,wsXPS.Nbe),H[1:wsXPS.Nbe,:]*ρ_clean+σ_1*ones(Cdouble,wsXPS.Nbe),alpha=0.5,color="tab:blue")
   l_scat_2        = scatter(collect(wsXPS.Nbe+1:2wsXPS.Nbe),d_c1s_Ek[:,4]-z_baseline_2,color="tab:orange")
   l_plot_2,       = plot(collect(wsXPS.Nbe+1:2wsXPS.Nbe),H[wsXPS.Nbe+1:2wsXPS.Nbe,:]*ρ_clean,color="tab:orange")
   l_fill_2        = fill_between(collect(wsXPS.Nbe+1:2wsXPS.Nbe),H[wsXPS.Nbe+1:2wsXPS.Nbe,:]*ρ_clean-σ_2*ones(Cdouble,wsXPS.Nbe),H[wsXPS.Nbe+1:2wsXPS.Nbe,:]*ρ_clean+σ_2*ones(Cdouble,wsXPS.Nbe),alpha=0.5,color="tab:orange")
   l_scat_3        = scatter(collect(2wsXPS.Nbe+1:3wsXPS.Nbe),d_c1s_Ek[:,6]-z_baseline_3,color="tab:green")
   l_plot_3,       = plot(collect(2wsXPS.Nbe+1:3wsXPS.Nbe),H[2wsXPS.Nbe+1:3wsXPS.Nbe,:]*ρ_clean,color="tab:green")
   l_fill_3        = fill_between(collect(2wsXPS.Nbe+1:3wsXPS.Nbe),H[2wsXPS.Nbe+1:3wsXPS.Nbe,:]*ρ_clean-σ_3*ones(Cdouble,wsXPS.Nbe),H[2wsXPS.Nbe+1:3wsXPS.Nbe,:]*ρ_clean+σ_3*ones(Cdouble,wsXPS.Nbe),alpha=0.5,color="tab:green")
   l_scat_4        = scatter(collect(3wsXPS.Nbe+1:4wsXPS.Nbe),d_c1s_Ek[:,8]-z_baseline_4,color="tab:red")
   l_plot_4,       = plot(collect(3wsXPS.Nbe+1:4wsXPS.Nbe),H[3wsXPS.Nbe+1:4wsXPS.Nbe,:]*ρ_clean,color="tab:red")
   l_fill_4        = fill_between(collect(3wsXPS.Nbe+1:4wsXPS.Nbe),H[3wsXPS.Nbe+1:4wsXPS.Nbe,:]*ρ_clean-σ_4*ones(Cdouble,wsXPS.Nbe),H[3wsXPS.Nbe+1:4wsXPS.Nbe,:]*ρ_clean+σ_4*ones(Cdouble,wsXPS.Nbe),alpha=0.5,color="tab:red")
   l_scat_5        = scatter(collect(4wsXPS.Nbe+1:5wsXPS.Nbe),d_c1s_Ek[:,10]-z_baseline_5,color="tab:cyan")
   l_plot_5,       = plot(collect(4wsXPS.Nbe+1:5wsXPS.Nbe),H[4wsXPS.Nbe+1:5wsXPS.Nbe,:]*ρ_clean,color="tab:cyan")
   l_fill_5        = fill_between(collect(4wsXPS.Nbe+1:5wsXPS.Nbe),H[4wsXPS.Nbe+1:5wsXPS.Nbe,:]*ρ_clean-σ_5*ones(Cdouble,wsXPS.Nbe),H[4wsXPS.Nbe+1:5wsXPS.Nbe,:]*ρ_clean+σ_5*ones(Cdouble,wsXPS.Nbe),alpha=0.5,color="tab:cyan")

   xlabel("SVD left index",fontsize=12)
   ylabel("PE signal [a.u.]",fontsize=12)
   legend([l_scat_1,(l_plot_1,l_fill_1),l_scat_2,(l_plot_2,l_fill_2),l_scat_3,(l_plot_3,l_fill_3),l_scat_4,(l_plot_4,l_fill_4),l_scat_5,(l_plot_5,l_fill_5)],["data 1"; "reconstruction 1"; "data 2"; "reconstruction 2"; "data 3"; "reconstruction 3"; "data 4"; "reconstruction 4"; "data 5"; "reconstruction 5"])






   ##
   ## cp
   ##



   A = [H; Matrix{Cdouble}(I,N,N); D_2nd; B];
   W_stop = ones(Cdouble,N);
   τ0 = 1.0e1 # 4
   x0 = zeros(Cdouble,N);
   s0 = A*x0;
   N_max_iter = 2*1000#0 # 00;
   r_n_tol=1.0
   r_y_tol=0.005;
   γ_H = 1.0e-2H;

   # xn,sn,taun,N_last = alg2_cp_gaussian_un_no_mem(x0,s0,IA_1,A,γ,γ_H,γ_D,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
   # xn,sn,taun,F_ALL,G_ALL,F_STAR_ALL,G_STAR_ALL,INNER_PROD,T_ALL,N_last= alg2_cp_gaussian_un_no_mem(x0,s0,IA_1,A,γ,γ_H,γ_D,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);

   Ns = 5 # 00;
   ρ_samples = zeros(Cdouble,Ns,N);

   ρ_samples[1,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,H*ρA_1,ρ0,ρB,A,γ,γ_H,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
   t_elapsed = @elapsed for i in 2:Ns
      idx = shuffle_data(wsXPS.Nke*wsXPS.Nbe,wsXPS.Nke;Nmin=Nmin);
      # A = [H[idx,:]; Matrix{Cdouble}(I,N,N); D_2nd; B];
      x0 = zeros(Cdouble,N);
      s0 = A*x0
      # ρ_samples[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,IA_1[idx],ρ0,ρB,A,γ[idx],γ_H[idx,:],γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
      ρ_samples[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,H*ρA_1+sqrt.(γ).*randn(wsXPS.Nke*wsXPS.Nbe),ρ0,ρB,A,γ,γ_H,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
   end

   μ_ρ = dropdims(mean(ρ_samples,dims=1),dims=1);
   σ_ρ = dropdims(std(ρ_samples,dims=1),dims=1);






   fig1 = figure()
   # state reconstruction
   ax1 = subplot(121)
   l_scat_1               = scatter(Zi,ρA_1)
   # l_plot_1_gt,           = plot(Zi,F.V[:,1:rankK]*F.Vt[1:rankK,:]*ρA_1[idx_res],color="tab:blue");
   l_plot_1_nso_un,       = plot(Zi,ρ_samples[1,:],color="tab:orange")
   l_plot_1_nso_un_sam,   = plot(Zi,μ_ρ,color="tab:pink")
   l_plot_1_nso_un_sam_σ  = fill_between(Zi,μ_ρ-σ_ρ,μ_ρ+σ_ρ,alpha=0.5,color="tab:pink")
   # l_plot_1_nso_un_sam_σ  = fill_between(Zi,μ_ρ-2.0σ_ρ,μ_ρ+2.0σ_ρ,alpha=0.5,color="tab:pink")
   xlabel("depth [nm]",fontsize=12)
   ylabel("concentration [a.u.]",fontsize=12)
   # legend([l_scat_1,l_plot_1_gt,l_plot_1_nso_un,(l_plot_1_nso_un_sam,l_plot_1_nso_un_sam_σ)],["GT","GT low rank","NSO model uncertainty","NSO sampling"])
   legend([l_scat_1,l_plot_1_nso_un,(l_plot_1_nso_un_sam,l_plot_1_nso_un_sam_σ)],["GT","CP model uncertainty","CP sampling"])


   # reconstructed data
   ax2 = subplot(122)
   # l_scat_2               = scatter(Ke,IA_1_clean,color="tab:green")
   l_scat_2_data          = scatter(wsXPS.Ke,H*ρA_1,color="tab:blue")
   l_plot_2_nso_un,       = plot(wsXPS.Ke,H*ρ_samples[1,:],color="tab:orange")
   l_plot_2_nso_un_sam,   = plot(wsXPS.Ke,H*μ_ρ,color="tab:pink")
   l_plot_2_nso_un_sam_σ  = fill_between(wsXPS.Ke,H*μ_ρ-sqrt.(γ),H*μ_ρ+sqrt.(γ),alpha=0.5,color="tab:pink")
   xlabel("kinetic energy [eV]",fontsize=12)
   ylabel("PE signal [a.u.]",fontsize=12)
   # legend([l_scat_2_data,l_scat_2,l_plot_2_nso_un,(l_plot_2_nso_un_sam,l_plot_2_nso_un_sam_σ)],["data","GT","NSO model uncertainty","NSO sampling"])
   legend([l_scat_2_data,l_plot_2_nso_un,(l_plot_2_nso_un_sam,l_plot_2_nso_un_sam_σ)],["data","CP model uncertainty","CP sampling"])


   fig1.set_figwidth(10.0) # (5.6)
   fig1.set_figheight(4.5)
   # tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
   tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)

   ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.97), textcoords="axes fraction", color="black",fontsize=14)
   ax1.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(1.02, 0.97), textcoords="axes fraction", color="black",fontsize=14)
end
