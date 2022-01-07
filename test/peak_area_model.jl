## load the packages used in the estimation
# plotting
using PyPlot
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
using XPSinv


##
## load some data
##
σ_noise = 0.002;
data_folder = "./data/lin/4peaks/0.0020/";
folder_res = "./results/lin/4peaks/0.0020/"

# σ_noise = 0.02;
# data_folder = "./data/lin/4peaks/0.0200/";
# folder_res = "./results/lin/4peaks/0.0200/"

DATA1 = false;
DATA2 = false;
DATA3 = true;

PLOT_MEAS_OP = false;

mkpath(folder_res)
if DATA1
   d_IA_1 = CSV.File(string(data_folder,"IA_1.csv"); delim=",", header=false) |> DataFrame;
   IA_1 = Matrix{Cdouble}(d_IA_1)[1,:];
   d_IA_1_clean = CSV.File(string(data_folder,"IA_1_clean.csv"); delim=",", header=false) |> DataFrame;
   IA_1_clean = Matrix{Cdouble}(d_IA_1_clean)[1,:];
   d_ρA_1 = CSV.File(string(data_folder,"rhoA_1.csv"); delim=",", header=false) |> DataFrame;
   ρA_1 = Matrix{Cdouble}(d_ρA_1)[1,:];
end

if DATA2
   d_IA_1 = CSV.File(string(data_folder,"IA_2.csv"); delim=",", header=false) |> DataFrame;
   IA_1 = Matrix{Cdouble}(d_IA_1)[1,:];
   d_IA_1_clean = CSV.File(string(data_folder,"IA_2_clean.csv"); delim=",", header=false) |> DataFrame;
   IA_1_clean = Matrix{Cdouble}(d_IA_1_clean)[1,:];
   d_ρA_1 = CSV.File(string(data_folder,"rhoA_2.csv"); delim=",", header=false) |> DataFrame;
   ρA_1 = Matrix{Cdouble}(d_ρA_1)[1,:];
end

if DATA3
   d_IA_1 = CSV.File(string(data_folder,"IA_3.csv"); delim=",", header=false) |> DataFrame;
   IA_1 = Matrix{Cdouble}(d_IA_1)[1,:];
   d_IA_1_clean = CSV.File(string(data_folder,"IA_3_clean.csv"); delim=",", header=false) |> DataFrame;
   IA_1_clean = Matrix{Cdouble}(d_IA_1_clean)[1,:];
   d_ρA_1 = CSV.File(string(data_folder,"rhoA_3.csv"); delim=",", header=false) |> DataFrame;
   ρA_1 = Matrix{Cdouble}(d_ρA_1)[1,:];
end
d_Zi_high_res = CSV.File(string(data_folder,"depth.csv"); delim=",", header=false) |> DataFrame;
Zi_high = Matrix{Cdouble}(d_Zi_high_res)[1,:];


##
## load the model
##
model_folder = "./data/lin/4peaks/low_res_50/";
d_H = CSV.File(string(model_folder,"H.csv"); delim=",", header=false) |> DataFrame;
H = Matrix{Cdouble}(d_H);
d_H_std = CSV.File(string(model_folder,"H_standard_deviation.csv"); delim=",", header=false) |> DataFrame;
H_std = Matrix{Cdouble}(d_H_std);
d_Zi = CSV.File(string(model_folder,"depth.csv"); delim=",", header=false) |> DataFrame;
Zi = Matrix{Cdouble}(d_Zi)[1,:];
d_Ke = CSV.File(string(model_folder,"kinetic_energy.csv"); delim=",", header=false) |> DataFrame;
Ke = Matrix{Cdouble}(d_Ke)[1,:];




##
## peak area model
##
dKe1 = Ke[2]-Ke[1];
dKe2 = Ke[102]-Ke[101];
dKe3 = Ke[202]-Ke[201];
dKe4 = Ke[302]-Ke[301];

A1 = 1.0dKe1*dropdims(sum(H[1:100,:],dims=1),dims=1);
A2 = 1.0dKe2*dropdims(sum(H[101:200,:],dims=1),dims=1);
A3 = 2.0dKe3*dropdims(sum(H[201:300,:],dims=1),dims=1);
A4 = 2.0dKe4*dropdims(sum(H[301:400,:],dims=1),dims=1);

A = Matrix{Cdouble}([A1 A2 A3 A4]');

figure(); plot(Zi,A')
figure(); plot(Zi,A[2:4,:]'-A[1:3,:]')

N = length(Zi);
Nke = length(Ke);
idx_res = zeros(Int64,N);
for i in 1:N-1
   idx_res[i] = findfirst(Zi_high.>=Zi[i])
end
idx_res[N] = length(Zi_high)

##
## compute the peaks area from the spectra
##
y_all = [sum(IA_1[1:100])*dKe1; sum(IA_1[101:200])*dKe2; sum(IA_1[201:300])*dKe3; sum(IA_1[301:400])*dKe4];
y_all_clean = A*ρA_1[idx_res];
σ_all = σ_noise*sqrt(100.0)*[dKe1; dKe2; dKe3; dKe4];
figure(); scatter(A*ρA_1[idx_res],y_all-A*ρA_1[idx_res]); fill_between(A*ρA_1[idx_res], -σ_all, σ_all)
figure(); scatter(collect(1:Nke),H*ρA_1[idx_res]); scatter(collect(1:Nke),IA_1); fill_between(collect(1:Nke),H*ρA_1[idx_res].-σ_noise,H*ρA_1[idx_res].+σ_noise,alpha=0.5,color="tab:pink")


##
## elements of inversion
##
D_1st = D1st(N);
D_2nd = D2nd(N);
d_2nd = 2;
B = zeros(Cdouble,2,N); B[1,1] = 1.0; B[2,end] = 1.0
γ = σ_all.^2 # σ_noise^2*ones(Cdouble,Nke);
Γ_inv = diagm(1.0./γ);
z00  = 2.0;
σ_z0 = 1.0;
# amplitude of the first and second order derivatives (δ_D and ι_D)
if DATA1
   δ_D = 10.0*(abs.(D_1st*logistic.(Zi.-z00,0.0,1.0,σ_z0)) + abs.(D_1st*logistic.(Zi.-1.0,0.0,1.0,σ_z0)) + abs.(D_1st*logistic.(Zi.-3.0,0.0,1.0,σ_z0)));
   ι_D = 10.0*(abs.(D_2nd*logistic.(Zi.-z00,0.0,1.0,σ_z0)) + abs.(D_2nd*logistic.(Zi.-1.0,0.0,1.0,σ_z0)) + abs.(D_2nd*logistic.(Zi.-3.0,0.0,1.0,σ_z0)));
else
   δ_D = 5.0*logistic.(Zi.-3.0,0.2,0.0,2.0σ_z0)[1:end-1];
   ι_D = 10.0*logistic.(Zi.-3.0,0.2,0.0,2.0σ_z0)[2:end-1];
   if DATA2
      δ_D = 2.0δ_D
      # ι_D = 20.0ι_D
      ι_D = 2.0ι_D
   end
end
γ_D = (ι_D .+ 0.01maximum(ι_D)).^2;
Γ_D_inv = diagm(1.0./γ_D);
ρ0 = ρA_1[idx_res[1]];
γ0 = 0.01^2
ρB = ρA_1[idx_res[end]];
γB = 0.01^2
ρ_bar = [ρ0;ρB];
γ_ρ = [γ0;γB];
Γ_ρ_inv = diagm(1.0./γ_ρ);
rankK = rank(A);

A1_std = 1.0*dKe1*sqrt.(dropdims(sum(H_std[1:100,:].^2,dims=1),dims=1));
A2_std = 1.0*dKe2*sqrt.(dropdims(sum(H_std[101:200,:].^2,dims=1),dims=1));
A3_std = 2.0*dKe3*sqrt.(dropdims(sum(H_std[201:300,:].^2,dims=1),dims=1));
A4_std = 2.0*dKe4*sqrt.(dropdims(sum(H_std[301:400,:].^2,dims=1),dims=1));
A_std = Matrix{Cdouble}([A1_std A2_std A3_std A4_std]');
γ_H = A_std # 1.0e-2A_std;
γ = 1.0*σ_all # σ_noise*ones(Cdouble,Nke);


##
## NSO (conditional to data and model)
##
# ρ_clean_nso,_,F = iterative_nso(rankK,A,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,y_all_clean,ρ_bar);
ρ_clean_nso,_,F = iterative_nso(rankK,A,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,y_all,ρ_bar);


figure(); plot(Zi,ρ_clean_nso); plot(Zi,ρA_1[idx_res])

##
## NSO (sampling model and data)
##
Nh = 121
Ns = 25
ρ_samples_nso_h = zeros(Cdouble,Nh*Ns,N);
t_elapsed = @elapsed for i in 1:Nh
   d_H = CSV.File(string(model_folder,i,"/","H.csv"); delim=",", header=false) |> DataFrame;
   H = 10.0Matrix{Cdouble}(d_H);

   A1 = 1.0*dKe1*dropdims(sum(H[1:100,:],dims=1),dims=1);
   A2 = 1.0*dKe2*dropdims(sum(H[101:200,:],dims=1),dims=1);
   A3 = 2.0*dKe3*dropdims(sum(H[201:300,:],dims=1),dims=1);
   A4 = 2.0*dKe4*dropdims(sum(H[301:400,:],dims=1),dims=1);

   A = Matrix{Cdouble}([A1 A2 A3 A4]');

   # compute reconstruction with the given model
   for j in 1:Ns
      ρ_samples_nso_h[(i-1)*Ns+j,:],_,_ = iterative_nso(rankK,A,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,y_all_clean+γ.*randn(4),ρ_bar);
   end
end

μ_ρ_nso_h = dropdims(mean(ρ_samples_nso_h,dims=1),dims=1)
σ_ρ_nso_h = dropdims(std(ρ_samples_nso_h,dims=1),dims=1)

figure(); plot(Zi,μ_ρ_nso_h); plot(Zi,ρA_1[idx_res]); fill_between(Zi,μ_ρ_nso_h-σ_ρ_nso_h,μ_ρ_nso_h+σ_ρ_nso_h,alpha=0.5,color="tab:pink")

##
## CP
##

AA = [A; Matrix{Cdouble}(I,N,N); D_2nd; B];
W_stop = ones(Cdouble,N);
τ0 = 1.0e1 # 4
x0 = 0.0*ρA_1[idx_res]
s0 = AA*x0
N_max_iter = 20000 # 00;
r_n_tol=0.1*1.0
r_y_tol=0.1*0.005;


Ns = 5 # 00;
ρ_samples = zeros(Cdouble,Ns,N);

ρ_samples[1,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all_clean,ρ0,ρB,AA,γ,γ_H,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
t_elapsed = @elapsed for i in 2:Ns
   # idx = shuffle_data(Nke,Npeak;Nmin=Nmin);
   # A = [H[idx,:]; Matrix{Cdouble}(I,N,N); D_2nd; B];
   x0 = 0.0*ρA_1[idx_res]
   s0 = AA*x0
   # ρ_samples[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,IA_1[idx],ρ0,ρB,A,γ[idx],γ_H[idx,:],γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
   ρ_samples[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all_clean+γ.*randn(4),ρ0,ρB,AA,γ,γ_H,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
end

μ_ρ = dropdims(mean(ρ_samples,dims=1),dims=1);
σ_ρ = dropdims(std(ρ_samples,dims=1),dims=1);

figure(); plot(Zi,μ_ρ); plot(Zi,ρA_1[idx_res]); fill_between(Zi,μ_ρ-σ_ρ,μ_ρ+σ_ρ,alpha=0.5,color="tab:pink")




##
## more plotting
##
if PLOT_MEAS_OP
   minH,maxH = extrema(A)
   fig1, ax1, pcm1 = imshowData(1254,Zi,collect(1.0:4.0),A,_norm=:Normalize,_vmin=minH,_vmax=maxH,_edgecolors="face",_sub=221)
   xlabel("depth \$z\$",fontsize=14) # [nm]
   ylabel("Kinetic energy \$K_e\$",fontsize=14)
   xlim(Zi[1],Zi[end])
   ylim(1.0,5.0)
   rc("ytick",color="white")
   cax1 = fig1.add_axes([0.08+0.33, .75, 0.02, 0.2])
   cb1 = fig1.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
   cb1.set_label("measurement operator [a.u.]", color="white", fontsize=10)
   cb1.ax.yaxis.set_tick_params(color="white")
   cb1.outline.set_edgecolor("white")
   rc("ytick",color="black")


   minU,maxU = extrema(F.U)
   fig1, ax2, pcm2 = imshowData(1254,collect(1.0:4.0),collect(1.0:4.0),F.U,_norm=:Normalize,_vmin=minU,_vmax=maxU,_edgecolors="face",_sub=222)
   xlabel("basis index [data space]",fontsize=14) # [nm]
   ylabel("Kinetic energy \$K_e\$",fontsize=14)
   xlim(1.0,5.0)
   ylim(1.0,5.0)
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
