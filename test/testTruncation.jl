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

D_2nd = D2nd(N);

A = [H_tilde; Matrix{Cdouble}(I,N,N); 10.0D_2nd]; #  B];
w = (1.0/σnoise[1])*ones(Cdouble,Ndata);
W_stop = ones(Cdouble,N);
τ0 = 1.0e1 # 4
x00 = 0.5ones(Cdouble,N); # 0.0*ρA_1[idx_res]
s0 = A*x00;
N_max_iter = 2000#0 # 00;
r_n_tol=0.1*1.0
r_y_tol=0.005;


##
## other truncation
##

H0 = H_better[:,1];
H_tilde = H_better[:,2:N+1];
Hb = H_better[:,N+2:end];

Δy = Hb*ρA_1[N+2:end];
δy = H0*ρA_1[1];
y_tilde = y_data.-(Δy+δy)';


DN = 1.0D2nd(N+3);
D0 = DN[:,1];
D_tilde = DN[:,2:N+1];
Db = DN[:,N+2:N+3];

Δyd = -Db*ρA_1[N+2:N+3]
δyd = -D0*ρA_1[1];
yd = Δyd+δyd;

Y = [y_tilde[1,:]; yd];
Htrunc = [H_tilde; D_tilde];


# covariance matrix for the a priori distribution

# Γprior = zeros(Cdouble,Nr,Nr);
# cor_len = 5.0;
# for i in 1:Nr
#     # Γprior[i,i] = (1.0-ρA_1[i]+0.1)^2;
#     Γprior[i,i] =  1.0
#     for j in i+1:Nr
#         Γprior[i,j] = Γprior[i,i]*exp(-(i-j)^2/(0.5*cor_len^2))
#         Γprior[j,i] = Γprior[i,j]
#     end
# end
# Γytilde = ΓI[:,:,1] + 1.0e-6*(Γprior[1,1]*H0*H0' + Hb*Γprior[N+2:end,N+2:end]*Hb');
# Γd = inv(D_tilde'*inv(Γprior)*D_tilde)

ρ_trunc = inv(Htrunc'*Htrunc)*Htrunc'*Y;
ρ_est_block = [ρA_1[1]; ρ_trunc; ρA_1[N+2:end]]

figure()
plot(ρ_est_block)
plot(ρA_1)

w = [(1.0/σnoise[1])*ones(Cdouble,Ndata); (1.0/0.1)*ones(Cdouble,N+1)];
W_stop = ones(Cdouble,N);
τ0 = 1.0e1 # 4
x00 = 0.5ones(Cdouble,N); # 0.0*ρA_1[idx_res]
N_max_iter = 2000000 # 00;
r_n_tol=0.0001*1.0
r_y_tol=0.005;

ρ_est,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_gaussian(0.5ones(Cdouble,N),Y,Htrunc,w,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);

ρ_est_cp_block = [ρA_1[1]; ρ_est; ρA_1[N+2:end]]

figure()
plot(ρ_est_block)
plot(ρ_est_cp_block)
plot(ρA_1)

