## load the packages used in the estimation
# plotting
using PyPlot
fm = PyPlot.matplotlib.font_manager.json_load("/home/matthew/.cache/matplotlib/fontlist-v310.json") # TODO: look for the path automatically
# fm = PyPlot.matplotlib.font_manager.json_load("/home/mattoz/.cache/matplotlib/fontlist-v310.json")
# fm.findfont("serif", rebuild_if_missing=false)
# fm.findfont("serif", fontext="afm", rebuild_if_missing=false)
# rc("font",family="serif",serif="Computer Modern Roman")
rc("text", usetex=true)
rc("figure",max_open_warning=50)
using myPlot

# data manipulation (loading, writing, etc)
using Printf

# scientific package from the official Julia repositories
using LinearAlgebra
# using Statistics
# using DSP
# using SpecialMatrices
# using Polynomials
using StatsBase

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack
using XPSinv


# create some data with a fairly fine discretization
include("dataGenCylinder.jl")


##TODO: compute the correction with the other measurement model (different discretization and slightly different penetration depth (maybe add a random element in lambda?))
# data correction due to the signal coming from the deep strates of the surface for which we will assume that the bulk value has been reach for the concentration
Δy = H_better[:,N0+1:end]*ρA_1[N0+1:end];
# NOTE: maybe this part should be discussed (what should be the value at the "surface"?)
δy = H_better[:,1]*ρA_1[1];

figure(); plot(y_data'); plot((y_data.-Δy')'); plot((y_data.-(Δy+δy)')')
plot(Δy); plot(δy)

y_tilde = y_data.-(Δy+δy)';
H_tilde = H_better[:,2:N0]



N = size(H_tilde,2);
N = 51;
N_lowres = 13;

D_2nd = D2nd(N);
D_2nd_lowres = D2nd(N_lowres);

# A = [H_tilde; Matrix{Cdouble}(I,N,N); 10.0D_2nd];
# A_lowres = [H_tilde; Matrix{Cdouble}(I,N,N); 10.0D_2nd];
# w = (1.0/σnoise[1])*ones(Cdouble,Ndata);
# W_stop = ones(Cdouble,N);
# τ0 = 1.0e1 # 4
# x00 = 0.5ones(Cdouble,N); # 0.0*ρA_1[idx_res]
# # s0 = A*x00;
# N_max_iter = 2000#0 # 00;
# r_n_tol=0.1*1.0
# r_y_tol=0.005;


##
## other truncation
##

if false
    H0 = H_better[:,1];
    H_tilde = H_better[:,2:N0];
    Hb = H_better[:,N0+1:end];
    Δy = Hb*ρA_1[N0+1:end];
    δy = H0*ρA_1[1];
    # dropdims(sum(H_better[:,N0+1:end],dims=2),dims=2)
else
    H0 = H_lowres[:,1];
    H_tilde = H_lowres[:,2:N0_lowres];
    Hb = H_lowres[:,N0_lowres+1:end];
    Δy = dropdims(sum(Hb,dims=2),dims=2) # Hb*ρA_1[N0+1:end];
    δy = H0*ρA_1[1];
end


y_tilde = y_data.-(Δy+δy)';


if false
    DN = D2nd(N+3);
    D0 = DN[:,1];
    D_tilde = DN[:,2:N+1];
    Db = DN[:,N+2:N+3];

    Δyd = -Db*ρA_1[N+2:N+3]
    δyd = -D0*ρA_1[1];
else
    DN = D2nd(N_lowres+3); #+3
    D0 = DN[:,1];
    D_tilde = DN[:,2:N_lowres+1];
    Db = DN[:,N_lowres+2:N_lowres+3];

    Δyd = -dropdims(sum(Db,dims=2),dims=2) # Db*ρA_1[N+2:N+3]
    δyd = -D0*ρA_1[1];
end


yd = Δyd+δyd;



Htrunc = [H_tilde; D_tilde];

# covariance matrix for the a priori distribution

Γprior = zeros(Cdouble,Nr,Nr);
Γprior_lowres = zeros(Cdouble,Nr_lowres,Nr_lowres);
cor_len = 5.0;
cor_len_lowres = 5.0/4.0;
for i in 1:Nr
    # Γprior[i,i] = (1.0-ρA_1[i]+0.1)^2;
    Γprior[i,i] =  1.0
    for j in i+1:Nr
        Γprior[i,j] = Γprior[i,i]*exp(-(i-j)^2/(0.5*cor_len^2))
        Γprior[j,i] = Γprior[i,j]
    end
end
for i in 1:Nr_lowres
    # Γprior[i,i] = (1.0-ρA_1[i]+0.1)^2;
    Γprior_lowres[i,i] =  1.0
    for j in i+1:Nr_lowres
        Γprior_lowres[i,j] = Γprior_lowres[i,i]*exp(-(i-j)^2/(0.5*cor_len_lowres^2))
        Γprior_lowres[j,i] = Γprior_lowres[i,j]
    end
end
# Γytilde = ΓI[:,:,1] + 1.0e-6*(Γprior[1,1]*H0*H0' + Hb*Γprior[N+2:end,N+2:end]*Hb');
Γd = 0.0001*Γprior[2:N+2,2:N+2]; # inv(D_tilde'*inv(Γprior[2:N+1,2:N+1])*D_tilde)
Γd_lowres = 0.0001*Γprior_lowres[2:N_lowres+2,2:N_lowres+2];


# 
if false
    W_stop = ones(Cdouble,N);
else 
    W_stop_lowres = ones(Cdouble,N_lowres);
end
τ0 = 1.0e1
if false
    x00 = 0.5ones(Cdouble,N);
else
    x00 = 0.5ones(Cdouble,N_lowres);
end
N_max_iter = 200000 # 0;
r_n_tol=0.0001;
r_y_tol=0.005;

# for each noise sample
if false
    ρ_est_block    = zeros(Cdouble,Nnoise,Nr);
    ρ_est_cp_block = zeros(Cdouble,Nnoise,Nr);
else
    ρ_est_block    = zeros(Cdouble,Nnoise,Nr_lowres);
    ρ_est_cp_block = zeros(Cdouble,Nnoise,Nr_lowres);
end


for i in 1:Nnoise
    local Y = [y_tilde[i,:]; yd];

    local ρ_trunc = inv(Htrunc'*Htrunc)*Htrunc'*Y;

    # ρ_est_block[i,:] = [ρA_1[1]; ρ_trunc; ρA_1[N_lowres+2:end]]
    ρ_est_block[i,:] = [ρA_1[1]; ρ_trunc; ones(Cdouble,Nr_lowres-N0_lowres)]

    # local w = [(1.0/σnoise[i])*ones(Cdouble,Ndata); (1.0/0.1)*ones(Cdouble,N+1)];
    
    # local ρ_est,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_gaussian(0.5ones(Cdouble,N),Y,Htrunc,w,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
    # local ρ_est,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_quad(0.5ones(Cdouble,N),y_tilde[i,:],yd,Htrunc,ΓI[:,:,i],Γd::Array{Cdouble,2},W_stop;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)
    local ρ_est,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_quad(0.5ones(Cdouble,N_lowres),y_tilde[i,:],yd,Htrunc,ΓI[:,:,i],Γd_lowres,W_stop_lowres;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)

    # ρ_est_cp_block[i,:] = [ρA_1[1]; ρ_est; ρA_1[N_lowres+2:end]]
    ρ_est_cp_block[i,:] = [ρA_1[1]; ρ_est; ones(Cdouble,Nr_lowres-N0_lowres)]
    
end


μρ = dropdims(mean(ρ_est_cp_block,dims=1),dims=1);
Γρ = cov(ρ_est_cp_block);

figure()
plot(ρ_est_block')
plot(ρA_1)

figure()
# plot(ρ_est_block')
plot(ρ_est_cp_block')
plot(ρA_1)

if false
figure()
plot(μ0.-r,μρ)
fill_between(μ0.-r,μρ-sqrt.(diag(Γρ)),μρ+sqrt.(diag(Γρ)),alpha=0.5,color="tab:blue",label="uncertainty")
plot(μ0.-r,ρA_1)
end