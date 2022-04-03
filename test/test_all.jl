## load the packages used in the estimation
# plotting
using PyPlot
fm = PyPlot.matplotlib.font_manager.json_load("/home/matthew/.cache/matplotlib/fontlist-v310.json")
# fm.findfont("serif", rebuild_if_missing=false)
# fm.findfont("serif", fontext="afm", rebuild_if_missing=false)
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
data_folder = "/home/matthew/Data/XPS/lin/4peaks/0.0020/"; # "./data/lin/4peaks/0.0020/";
folder_res = "./results/lin/4peaks/0.0020/"

# σ_noise = 0.02;
# data_folder = "./data/lin/4peaks/0.0200/";
# folder_res = "./results/lin/4peaks/0.0200/"

DATA1 = false;
DATA2 = false;
DATA3 = true;

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
model_folder = "/home/matthew/Data/XPS/lin/4peaks/low_res_50/";  #"./data/lin/4peaks/low_res_50/";
d_H = CSV.File(string(model_folder,"H.csv"); delim=",", header=false) |> DataFrame;
H = Matrix{Cdouble}(d_H);
d_H_std = CSV.File(string(model_folder,"H_standard_deviation.csv"); delim=",", header=false) |> DataFrame;
H_std = Matrix{Cdouble}(d_H_std);
d_Zi = CSV.File(string(model_folder,"depth.csv"); delim=",", header=false) |> DataFrame;
Zi = Matrix{Cdouble}(d_Zi)[1,:];
d_Ke = CSV.File(string(model_folder,"kinetic_energy.csv"); delim=",", header=false) |> DataFrame;
Ke = Matrix{Cdouble}(d_Ke)[1,:];


##
## define the matrices and vector of the problem
##
N = length(Zi);
Nke = length(Ke);
idx_res = zeros(Int64,N);
for i in 1:N-1
   idx_res[i] = findfirst(Zi_high.>=Zi[i])
end
idx_res[N] = length(Zi_high) # findfirst(Zi_high.>=Zi[N])
D_1st = D1st(N);
D_2nd = D2nd(N);
d_2nd = 2;
B = zeros(Cdouble,2,N); B[1,1] = 1.0; B[2,end] = 1.0
γ = σ_noise^2*ones(Cdouble,Nke);
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
rankK = rank(H);
γ_H = 1.0e-2H_std;
γ = σ_noise*ones(Cdouble,Nke);


##
## NSO
##


##
## draw samples from data distribution
##
Npeak = 4
Nmin = 50;

# the deal
μ_ρ_nso,σ_ρ_nso,ρ_samples_nso,F,W_inv = shuffle_data_sample_nso_un(100,Npeak,rankK,H,γ_H,σ_noise*ones(Cdouble,Nke),Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,IA_1,ρ_bar;Nmin=Nmin);

##
## draw samples from data and model distribution
##
ρ_clean_nso,_,F = iterative_nso(rankK,H,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,IA_1_clean,ρ_bar);

Nh = 121
Ns = 25
ρ_samples_nso_h = zeros(Cdouble,Nh*Ns,N);
t_elapsed = @elapsed for i in 1:Nh
   # model_folder = "./data/lin/4peaks/low_res_50/";
   local d_H = CSV.File(string(model_folder,i,"/","H.csv"); delim=",", header=false) |> DataFrame;
   local H = Matrix{Cdouble}(d_H);

   # compute reconstruction with the given model
   for j in 1:Ns
      ρ_samples_nso_h[(i-1)*Ns+j,:],_,_ = iterative_nso(rankK,H,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,IA_1_clean+σ_noise*randn(Nke),ρ_bar);
   end
end

μ_ρ_nso_h = dropdims(mean(ρ_samples_nso_h,dims=1),dims=1)
σ_ρ_nso_h = dropdims(std(ρ_samples_nso_h,dims=1),dims=1)


##
## CP
##

A = [H; Matrix{Cdouble}(I,N,N); D_2nd; B];
W_stop = ones(Cdouble,N);
τ0 = 1.0e1 # 4
x0 = 0.0*ρA_1[idx_res]
s0 = A*x0
N_max_iter = 2000#0 # 00;
r_n_tol=1.0
r_y_tol=0.005;


# xn,sn,taun,N_last = alg2_cp_gaussian_un_no_mem(x0,s0,IA_1,A,γ,γ_H,γ_D,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
# xn,sn,taun,F_ALL,G_ALL,F_STAR_ALL,G_STAR_ALL,INNER_PROD,T_ALL,N_last= alg2_cp_gaussian_un_no_mem(x0,s0,IA_1,A,γ,γ_H,γ_D,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);

Ns = 5 # 00;
ρ_samples = zeros(Cdouble,Ns,N);

ρ_samples[1,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,IA_1_clean,ρ0,ρB,A,γ,γ_H,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
t_elapsed = @elapsed for i in 2:Ns
   idx = shuffle_data(Nke,Npeak;Nmin=Nmin);
   # A = [H[idx,:]; Matrix{Cdouble}(I,N,N); D_2nd; B];
   local x0 = 0.0*ρA_1[idx_res]
   local s0 = A*x0
   # ρ_samples[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,IA_1[idx],ρ0,ρB,A,γ[idx],γ_H[idx,:],γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
   ρ_samples[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,IA_1_clean+σ_noise*randn(Nke),ρ0,ρB,A,γ,γ_H,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
end

μ_ρ = dropdims(mean(ρ_samples,dims=1),dims=1);
σ_ρ = dropdims(std(ρ_samples,dims=1),dims=1);


##
## draw samples from data and model distribution
##
ρ_clean_cp,_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,IA_1_clean,ρ0,ρB,A,γ,γ_H,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);

Nh = 121
Ns = 25
ρ_samples_cp_h = zeros(Cdouble,Nh*Ns,N);
t_elapsed = @elapsed for i in 1:Nh
   # model_folder = "./data/lin/4peaks/low_res_50/";
   local d_H = CSV.File(string(model_folder,i,"/","H.csv"); delim=",", header=false) |> DataFrame;
   local H = Matrix{Cdouble}(d_H);
   local A = [H; Matrix{Cdouble}(I,N,N); D_2nd; B];

   # compute reconstruction with the given model
   for j in 1:Ns
      ρ_samples_cp_h[(i-1)*Ns+j,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,IA_1_clean+σ_noise*randn(Nke),ρ0,ρB,A,γ,γ_H,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
   end
end

μ_ρ_cp_h = dropdims(mean(ρ_samples_cp_h,dims=1),dims=1)
σ_ρ_cp_h = dropdims(std(ρ_samples_cp_h,dims=1),dims=1)


##
## plot some results
##
fig1 = figure()
# state reconstruction
ax1 = subplot(121)
l_scat_1                 = scatter(Zi,ρA_1[idx_res])
l_plot_1_cp_un_sam,      = plot(Zi,μ_ρ,color="tab:pink")
l_plot_1_cp_h_un_sam,    = plot(Zi,μ_ρ_cp_h,color="tab:olive")
l_plot_1_nso_un_sam,     = plot(Zi,μ_ρ_nso,color="tab:cyan")
l_plot_1_nso_h_un_sam,   = plot(Zi,μ_ρ_nso_h,color="magenta")
l_plot_1_cp_un_sam_σ     = fill_between(Zi,μ_ρ-σ_ρ,μ_ρ+σ_ρ,alpha=0.5,color="tab:pink")
l_plot_1_cp_h_un_sam_σ   = fill_between(Zi,μ_ρ_cp_h-σ_ρ_cp_h,μ_ρ_cp_h+σ_ρ_cp_h,alpha=0.5,color="tab:olive")
l_plot_1_nso_un_sam_σ    = fill_between(Zi,μ_ρ_nso-σ_ρ_nso,μ_ρ_nso+σ_ρ_nso,alpha=0.5,color="tab:cyan")
l_plot_1_nso_h_un_sam_σ  = fill_between(Zi,μ_ρ_nso_h-σ_ρ_nso_h,μ_ρ_nso_h+σ_ρ_nso_h,alpha=0.5,color="magenta")
xlabel("depth [nm]",fontsize=12)
ylabel("concentration [a.u.]",fontsize=12)
legend([l_scat_1,(l_plot_1_cp_un_sam,l_plot_1_cp_un_sam_σ),(l_plot_1_cp_h_un_sam,l_plot_1_cp_h_un_sam_σ),(l_plot_1_nso_un_sam,l_plot_1_nso_un_sam_σ),(l_plot_1_nso_h_un_sam,l_plot_1_nso_h_un_sam_σ)],["GT","CP sampling","CP sampling (data and model)","NSO sampling","NSO sampling (data and model)"])


# reconstructed data
ax2 = subplot(122)
l_scat_2_data            = scatter(Ke,IA_1,color="tab:blue")
l_plot_2_cp_un_sam,      = plot(Ke,H*μ_ρ,color="tab:pink")
l_plot_2_cp_h_un_sam,    = plot(Ke,H*μ_ρ_cp_h,color="tab:olive")
l_plot_2_nso_un_sam,     = plot(Ke,H*μ_ρ_nso,color="tab:cyan")
l_plot_2_nso_h_un_sam,   = plot(Ke,H*μ_ρ_nso_h,color="magenta")
l_plot_2_cp_un_sam_σ     = fill_between(Ke,H*μ_ρ-σ_noise*ones(Cdouble,Nke),H*μ_ρ+σ_noise*ones(Cdouble,Nke),alpha=0.5,color="tab:pink")
l_plot_2_cp_h_un_sam_σ   = fill_between(Ke,H*μ_ρ_cp_h-σ_noise*ones(Cdouble,Nke),H*μ_ρ_cp_h+σ_noise*ones(Cdouble,Nke),alpha=0.5,color="tab:olive")
l_plot_2_nso_un_sam_σ    = fill_between(Ke,H*μ_ρ_nso-σ_noise*ones(Cdouble,Nke),H*μ_ρ_nso+σ_noise*ones(Cdouble,Nke),alpha=0.5,color="tab:cyan")
l_plot_2_nso_h_un_sam_σ  = fill_between(Ke,H*μ_ρ_nso_h-σ_noise*ones(Cdouble,Nke),H*μ_ρ_nso_h+σ_noise*ones(Cdouble,Nke),alpha=0.5,color="magenta")
xlabel("kinetic energy [eV]",fontsize=12)
ylabel("PE signal [a.u.]",fontsize=12)
legend([l_scat_2_data,(l_plot_2_cp_un_sam,l_plot_2_cp_un_sam_σ),(l_plot_2_cp_h_un_sam,l_plot_2_cp_h_un_sam_σ),(l_plot_2_nso_un_sam,l_plot_2_nso_un_sam_σ),(l_plot_2_nso_h_un_sam,l_plot_2_nso_h_un_sam_σ)],["data","CP sampling","CP sampling (data and model)","NSO sampling","NSO sampling (data and model)"])


fig1.set_figwidth(10.0)
fig1.set_figheight(4.5)
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)

ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.97), textcoords="axes fraction", color="black",fontsize=14)
ax1.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(1.02, 0.97), textcoords="axes fraction", color="black",fontsize=14)

# savefig("cp_un_reconstruction_resample_y_rho1.png")
# savefig("cp_un_reconstruction_resample_y_rho1.pdf")
# savefig("cp_un_reconstruction_generate_y_rho1.png")
# savefig("cp_un_reconstruction_generate_y_rho1.pdf")
# savefig("cp_un_reconstruction_resample_y_rho2.png")
# savefig("cp_un_reconstruction_resample_y_rho2.pdf")
# savefig("cp_un_reconstruction_generate_y_rho2.png")
# savefig("cp_un_reconstruction_generate_y_rho2.pdf")
# savefig("cp_un_reconstruction_resample_y_rho3.png")
# savefig("cp_un_reconstruction_resample_y_rho3.pdf")
# savefig("cp_un_reconstruction_generate_y_rho3.png")
# savefig("cp_un_reconstruction_generate_y_rho3.pdf")
